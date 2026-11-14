#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Utils.h>
#include <R_ext/Visibility.h>
#include <zlib.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace {

using std::string_view;
namespace fs = std::filesystem;

void check_interrupt(void*) { R_CheckUserInterrupt(); }

void checkpoint(std::size_t row) {
  if (row % 65536 == 0 && !R_ToplevelExec(check_interrupt, nullptr))
    throw std::runtime_error("Z matrix preparation interrupted.");
}

// zlib also reads uncompressed files. Buffer input so each call handles many
// records without crossing the R interface or materializing a study table.
class Input {
 public:
  explicit Input(const std::string& path) : path_(path) {
#ifdef _WIN32
    file_ = gzopen_w(fs::u8path(path).c_str(), "rb");
#else
    file_ = gzopen(path.c_str(), "rb");
#endif
    if (!file_) fail("Cannot open file");
  }
  ~Input() { if (file_) gzclose(file_); }
  Input(const Input&) = delete;
  Input& operator=(const Input&) = delete;

  bool line(std::string& result) {
    result.clear();
    while (true) {
      if (start_ == end_) {
        const int n = gzread(file_, buffer_.data(), buffer_.size());
        if (n <= 0) {
          int status = Z_OK;
          const char* message = gzerror(file_, &status);
          if (n < 0 || (status != Z_OK && status != Z_STREAM_END))
            fail(std::string("Cannot read input: ") + message);
          if (result.empty()) return false;
          break;
        }
        start_ = 0;
        end_ = static_cast<std::size_t>(n);
      }
      const char* begin = buffer_.data() + start_;
      const auto newline = static_cast<const char*>(std::memchr(begin, '\n', end_ - start_));
      if (newline) {
        result.append(begin, static_cast<std::size_t>(newline - begin));
        start_ += static_cast<std::size_t>(newline - begin) + 1;
        break;
      }
      result.append(begin, end_ - start_);
      start_ = end_;
    }
    if (!result.empty() && result.back() == '\r') result.pop_back();
    return true;
  }

  [[noreturn]] void fail(const std::string& why) const {
    throw std::runtime_error(why + " in " + path_);
  }

  void close() {
    const auto status = gzclose(file_);
    file_ = nullptr;
    if (status != Z_OK) fail("Cannot close input");
  }

 private:
  std::string path_;
  gzFile file_;
  std::array<char, 65536> buffer_;
  std::size_t start_ = 0, end_ = 0;
};

// gact writes one record per line. Quoted fields may contain delimiters,
// including doubled quotes in fields the writer does not use.
void split(string_view line, char delimiter, std::vector<string_view>& fields, const Input& input) {
  fields.clear();
  std::size_t start = 0;
  while (true) {
    std::size_t end;
    if (start < line.size() && line[start] == '"') {
      end = start + 1;
      while (true) {
        end = line.find('"', end);
        if (end == string_view::npos) input.fail("Unclosed quoted field");
        ++end;
        if (end == line.size() || line[end] != '"') break;
        ++end;
      }
      if (end != line.size() && line[end] != delimiter) input.fail("Invalid quoted field");
      fields.push_back(line.substr(start + 1, end - start - 2));
    } else {
      end = line.find(delimiter, start);
      if (end == string_view::npos) end = line.size();
      fields.push_back(line.substr(start, end - start));
    }
    if (end == line.size()) return;
    start = end + 1;
  }
}

class Table {
 public:
  explicit Table(const std::string& path) : input(path) {
    if (!input.line(line)) input.fail("Missing header");
    if (line.compare(0, 3, "\xef\xbb\xbf") == 0) line.erase(0, 3);
    delimiter = line.find('\t') != std::string::npos ? '\t' : ',';
    split(line, delimiter, fields, input);
    width = fields.size();
    for (std::size_t i = 0; i < width; ++i) {
      if (!header.emplace(fields[i], i).second) input.fail("Duplicate column name");
    }
  }

  std::size_t column(const std::string& name) const {
    const auto at = header.find(name);
    return at == header.end() ? missing : at->second;
  }

  std::size_t required(const std::string& name) const {
    const auto index = column(name);
    if (index == missing) input.fail("Missing column " + name);
    return index;
  }

  bool next() {
    while (input.line(line)) {
      if (line.empty()) continue;
      split(line, delimiter, fields, input);
      if (fields.size() != width) input.fail("Unexpected column count");
      checkpoint(++row);
      return true;
    }
    return false;
  }

  string_view operator[](std::size_t index) const {
    return index == missing ? string_view{} : fields[index];
  }

  static constexpr std::size_t missing = std::numeric_limits<std::size_t>::max();
  string_view raw_line() const { return line; }
  Input input;

 private:
  std::string line;
  char delimiter;
  std::size_t width, row = 0;
  std::unordered_map<std::string, std::size_t> header;
  std::vector<string_view> fields;
};

std::ofstream output_file(const fs::path& path) {
  std::ofstream out;
  out.exceptions(std::ios::failbit | std::ios::badbit);
  out.open(path, std::ios::binary);
  return out;
}

class TextOutput {
 public:
  TextOutput(const std::string& path, bool gzip) : path_(path) {
    if (gzip) {
#ifdef _WIN32
      file_ = gzopen_w(fs::u8path(path).c_str(), "wb");
#else
      file_ = gzopen(path.c_str(), "wb");
#endif
      if (!file_) throw std::runtime_error("Cannot create " + path);
    } else {
      out_ = output_file(fs::u8path(path));
    }
  }
  ~TextOutput() { if (file_) gzclose(file_); }

  void line(string_view line) {
    if (file_) {
      if (line.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
          gzwrite(file_, line.data(), static_cast<unsigned>(line.size())) != static_cast<int>(line.size()) ||
          gzputc(file_, '\n') < 0)
        throw std::runtime_error("Cannot write " + path_);
    } else {
      out_ << line << '\n';
    }
  }

  void close() {
    if (file_) {
      const auto status = gzclose(file_);
      file_ = nullptr;
      if (status != Z_OK) throw std::runtime_error("Cannot finish " + path_);
    } else {
      out_.close();
    }
  }

 private:
  std::string path_;
  std::ofstream out_;
  gzFile file_ = nullptr;
};

bool valid_field(string_view value) {
  return !value.empty() && value.find_first_of("\t\r\n") == string_view::npos;
}

using Index = std::unordered_map<string_view, std::size_t>;

double number(string_view token) {
  if (token.empty() || token == "NA" || token == ".")
    return std::numeric_limits<double>::quiet_NaN();
  // R's native parser uses a decimal point independently of the locale.
  // The field ends at a delimiter, closing quote, or the string's NUL.
  char* end = nullptr;
  const double value = R_strtod(token.data(), &end);
  if (end == token.data()) return std::numeric_limits<double>::quiet_NaN();
  const auto limit = token.data() + token.size();
  while (end < limit && *end == ' ') ++end;
  return end == limit ? value : std::numeric_limits<double>::quiet_NaN();
}

char uppercase(char letter) {
  return letter >= 'a' && letter <= 'z' ? letter - 'a' + 'A' : letter;
}

struct CaseInsensitive {
  using is_transparent = void;

  bool operator()(string_view left, string_view right) const {
    return std::lexicographical_compare(left.begin(), left.end(), right.begin(), right.end(),
        [](char a, char b) { return uppercase(a) < uppercase(b); });
  }
};

// Intern repeated chromosome and allele strings instead of copying them per row.
class Labels {
 public:
  std::uint32_t find(string_view text) const {
    const auto at = values_.find(text);
    return at == values_.end() ? 0 : at->second;
  }

  std::uint32_t add(string_view text) {
    const auto key = find(text);
    if (key != 0) return key;
    return values_.emplace(text, values_.size() + 1).first->second;
  }

 private:
  std::map<std::string, std::uint32_t, CaseInsensitive> values_;
};

string_view chromosome(string_view name) {
  if (name.size() > 3 && uppercase(name[0]) == 'C' &&
      uppercase(name[1]) == 'H' && uppercase(name[2]) == 'R')
    name.remove_prefix(3);
  if (name == "X" || name == "x") return "23";
  if (name == "Y" || name == "y") return "24";
  if (name == "MT" || name == "mt" || name == "M" || name == "m") return "25";
  return name;
}

struct Marker {
  std::uint64_t pos = 0;
  std::uint32_t chr = 0, ea = 0, nea = 0;
};

struct MarkerFields {
  MarkerFields(const Table& table, bool ld)
      : chr(table.required(ld ? "chromosome" : "chr")),
        pos(table.required(ld ? "position_bp" : "pos")),
        ea(table.required(ld ? "allele1" : "ea")),
        nea(table.required(ld ? "allele2" : "nea")) {}

  std::uint64_t position(const Table& table) const {
    const double value = number(table[pos]);
    if (!(value > 0 && value <= 9007199254740992.0 && std::floor(value) == value)) return 0;
    return static_cast<std::uint64_t>(value);
  }

  std::size_t chr, pos, ea, nea;
};

int allele_sign(const Marker& reference, const MarkerFields& fields, const Table& table,
                const Labels& labels, string_view id) {
  if (reference.pos != fields.position(table) ||
      reference.chr != labels.find(chromosome(table[fields.chr])))
    table.input.fail("Coordinate mismatch for " + std::string(id));
  const auto ea = labels.find(table[fields.ea]);
  const auto nea = labels.find(table[fields.nea]);
  if (ea == reference.ea && nea == reference.nea) return 1;
  if (ea == reference.nea && nea == reference.ea) return -1;
  table.input.fail("Allele mismatch for " + std::string(id));
}

void write_markers(const std::string& path, const fs::path& output, const Index& index,
                   std::vector<Marker>& markers, Labels& labels) {
  auto out = output_file(output);
  out << "rsids\tchr\tpos\tea\tnea\n";
  // Most GAlists retain catalog order. Stream those rows directly and only
  // collect the catalog in memory if the caller supplies a different order.
  for (bool reorder : {false, true}) {
    Table table(path);
    const MarkerFields identity(table, false);
    const std::array<std::size_t, 5> fields = {
        table.required("rsids"), identity.chr, identity.pos, identity.ea, identity.nea};
    std::vector<std::string> records;
    if (reorder) {
      out.seekp(0, std::ios::beg);
      out << "rsids\tchr\tpos\tea\tnea\n";
      records.resize(index.size());
    }
    std::size_t count = 0;
    bool ordered = true;
    std::string record;
    while (table.next()) {
      const auto id = table[fields[0]];
      const auto match = index.find(id);
      if (match == index.end()) continue;
      if (!reorder && match->second != count) {
        ordered = false;
        break;
      }
      record.clear();
      for (const auto field : fields) {
        const auto value = table[field];
        if (!valid_field(value)) table.input.fail("Invalid marker metadata for " + std::string(id));
        if (!record.empty()) record += '\t';
        record += value;
      }
      const Marker marker{identity.position(table), labels.add(chromosome(table[identity.chr])),
                          labels.add(table[identity.ea]), labels.add(table[identity.nea])};
      if (marker.pos == 0 || marker.ea == marker.nea)
        table.input.fail("Invalid marker position or alleles for " + std::string(id));
      markers[match->second] = marker;
      record += '\n';
      if (reorder) {
        auto& target = records[match->second];
        if (!target.empty()) table.input.fail("Repeated marker " + std::string(id));
        target = record;
      } else {
        out << record;
      }
      ++count;
    }
    if (!ordered) continue;
    if (count != index.size()) table.input.fail("Marker catalog does not cover GAlist$rsids");
    for (const auto& record : records) out << record;
    out.close();
    return;
  }
}

std::size_t study_id_column(const Table& table) {
  return table.column("rsids") == Table::missing ? table.required("marker") : table.column("rsids");
}

std::unordered_map<std::string, Marker> study_markers(const std::string& study,
    const std::string& catalog, Labels& labels) {
  std::unordered_map<std::string, Marker> markers;
  {
    Table table(study);
    const auto id = study_id_column(table);
    while (table.next()) {
      const auto rsid = table[id];
      if (!valid_field(rsid) || !markers.try_emplace(std::string(rsid)).second)
        table.input.fail("Invalid or repeated study marker " + std::string(rsid));
    }
  }
  Table table(catalog);
  const MarkerFields fields(table, false);
  const auto id = table.required("rsids");
  while (table.next()) {
    const auto at = markers.find(std::string(table[id]));
    if (at == markers.end()) continue;
    if (at->second.pos) table.input.fail("Repeated catalog marker " + at->first);
    const auto pos = fields.position(table);
    if (!pos || !valid_field(table[fields.chr]) ||
        !valid_field(table[fields.ea]) || !valid_field(table[fields.nea]))
      table.input.fail("Invalid catalog metadata for " + at->first);
    at->second = {pos, labels.add(chromosome(table[fields.chr])),
                  labels.add(table[fields.ea]), labels.add(table[fields.nea])};
    if (at->second.ea == at->second.nea)
      table.input.fail("Invalid catalog alleles for " + at->first);
  }
  return markers;
}

// Copy unchanged fields verbatim, including quoted CSV fields.
void replace_field(std::string& line, const char* origin, string_view old,
                   string_view replacement) {
  line.replace(static_cast<std::size_t>(old.data() - origin), old.size(), replacement);
}

std::string as_text(double value) {
  std::ostringstream out;
  out.imbue(std::locale::classic());
  out << std::setprecision(std::numeric_limits<double>::max_digits10) << value;
  return out.str();
}

SEXP recode_sumstat(SEXP source, SEXP catalog, SEXP stage, SEXP destination, SEXP gzip) {
  try {
    const std::string study_path = CHAR(STRING_ELT(source, 0));
    const std::string temporary = CHAR(STRING_ELT(stage, 0));
    Labels labels;
    const auto markers = study_markers(study_path, CHAR(STRING_ELT(catalog, 0)), labels);
    Table table(study_path);
    const MarkerFields fields(table, false);
    const auto id = study_id_column(table);
    const auto b = table.column("b"), se = table.column("seb");
    const auto eaf = table.column("eaf");
    const std::array<std::size_t, 2> z_columns = {table.column("z"), table.column("z_score")};
    TextOutput out(temporary, Rf_asLogical(gzip) == TRUE);
    out.line(table.raw_line());
    std::size_t reversed = 0, repaired = 0;
    std::string line;
    std::vector<std::pair<string_view, std::string>> edits;
    while (table.next()) {
      const auto key = table[id];
      const auto at = markers.find(std::string(key));
      if (at == markers.end() || !at->second.pos)
        table.input.fail("Study marker missing from allele catalog: " + std::string(key));
      const auto sign = allele_sign(at->second, fields, table, labels, key);
      edits.clear();
      const double beta = number(table[b]), error = number(table[se]);
      if (sign < 0) {
        ++reversed;
        edits.emplace_back(table[fields.ea], table[fields.nea]);
        edits.emplace_back(table[fields.nea], table[fields.ea]);
        if (std::isfinite(beta)) edits.emplace_back(table[b], as_text(-beta));
        const double frequency = number(table[eaf]);
        if (std::isfinite(frequency)) edits.emplace_back(table[eaf], as_text(1 - frequency));
      }
      for (const auto col : z_columns) {
        const double old_z = number(table[col]);
        if (!std::isfinite(old_z)) continue;
        double new_z = sign * old_z;
        const double expected = std::isfinite(beta) && error > 0 && std::isfinite(error) ?
            beta / error : 0;
        if (std::isfinite(expected) && expected != 0 &&
            std::signbit(expected) != std::signbit(old_z) &&
            std::abs(expected + old_z) <= 1e-6 * std::abs(expected)) {
          // Older ingestion could change beta and leave a stored Z unchanged.
          new_z = -new_z;
          ++repaired;
        }
        if (new_z != old_z) edits.emplace_back(table[col], as_text(new_z));
      }
      const auto raw = table.raw_line();
      if (edits.empty()) {
        out.line(raw);
        continue;
      }
      line = raw;
      std::sort(edits.begin(), edits.end(), [](const auto& a, const auto& b) {
        return a.first.data() > b.first.data();
      });
      for (const auto& edit : edits) replace_field(line, raw.data(), edit.first, edit.second);
      out.line(line);
    }
    out.close();
    table.input.close();
    const std::string target = CHAR(STRING_ELT(destination, 0));
    if (study_path == target) {
      if (reversed || repaired)
        fs::rename(fs::u8path(temporary), fs::u8path(target));
    } else {
      fs::create_hard_link(fs::u8path(temporary), fs::u8path(target));
    }
    SEXP counts = PROTECT(Rf_allocVector(REALSXP, 2));
    REAL(counts)[0] = static_cast<double>(reversed);
    REAL(counts)[1] = static_cast<double>(repaired);
    UNPROTECT(1);
    return counts;
  } catch (const std::exception& e) {
    Rf_error("%s", e.what());
  }
  return R_NilValue;
}

std::string ld_map_name(const std::string& panel) {
  return "ld-" + panel + ".u32le.bin";
}

// Zero denotes an uncovered LD row. Positive entries are one-based Z rows.
void write_ld_row(std::ofstream& out, std::uint32_t row) {
  const unsigned char bytes[] = {static_cast<unsigned char>(row),
      static_cast<unsigned char>(row >> 8), static_cast<unsigned char>(row >> 16),
      static_cast<unsigned char>(row >> 24)};
  out.write(reinterpret_cast<const char*>(bytes), sizeof(bytes));
}

struct LdCoverage {
  std::size_t rows, matched;
};

LdCoverage write_ld_map(const std::string& path, const fs::path& output, const Index& index,
                       const std::vector<Marker>& markers, const Labels& labels) {
  Table table(path);
  const MarkerFields fields(table, true);
  const auto idx = table.required("index"), id = table.required("marker_id");
  auto out = output_file(output);
  std::size_t row = 0, matched = 0;
  std::vector<bool> seen(markers.size());
  while (table.next()) {
    if (number(table[idx]) != static_cast<double>(row))
      table.input.fail("Expected consecutive zero-based sparse LD rows");
    const auto key = table[id];
    const auto at = index.find(key);
    std::uint32_t z_index = 0;
    if (at != index.end()) {
      const auto z_row = at->second;
      if (seen[z_row]) table.input.fail("Repeated sparse LD marker " + std::string(key));
      seen[z_row] = true;
      if (allele_sign(markers[z_row], fields, table, labels, key) < 0)
        table.input.fail("Reversed LD dosage allele for " + std::string(key) +
                         "; align the BED panel and its LD to the gact marker catalog first");
      z_index = static_cast<std::uint32_t>(z_row + 1);
      ++matched;
    }
    write_ld_row(out, z_index);
    ++row;
  }
  if (matched == 0) table.input.fail("No sparse LD markers overlap GAlist$rsids");
  out.close();
  return {row, matched};
}

std::size_t read_study(const std::string& path, const Index& index,
                       const std::vector<Marker>& markers, const Labels& labels,
                       std::vector<double>& column, std::vector<bool>& seen) {
  const auto nan = std::numeric_limits<double>::quiet_NaN();
  std::fill(column.begin(), column.end(), nan);
  std::fill(seen.begin(), seen.end(), false);
  Table table(path);
  const MarkerFields fields(table, false);
  const auto id = study_id_column(table);
  const auto z = table.column("z"), b = table.column("b"), se = table.column("seb");
  if (z == Table::missing && (b == Table::missing || se == Table::missing))
    table.input.fail("Expected z or both b and seb");
  std::size_t observed = 0;
  while (table.next()) {
    const auto key = table[id];
    const auto match = index.find(key);
    if (match == index.end()) continue;
    const auto row = match->second;
    if (seen[row]) table.input.fail("Repeated marker " + std::string(key));
    seen[row] = true;
    const auto sign = allele_sign(markers[row], fields, table, labels, key);
    const double effect = number(table[b]), error = number(table[se]);
    double value = std::isfinite(effect) && std::isfinite(error) && error > 0 ?
        effect / error : nan;
    if (!std::isfinite(value)) value = number(table[z]);
    if (std::isfinite(value)) {
      column[row] = sign * value;
      ++observed;
    }
  }
  return observed;
}

void write_column(std::ofstream& out, std::vector<double>& column) {
  static_assert(sizeof(double) == 8 && std::numeric_limits<double>::is_iec559,
                "The Z matrix requires IEEE FP64.");
  const std::uint16_t one = 1;
  auto bytes = reinterpret_cast<char*>(column.data());
  if (*reinterpret_cast<const unsigned char*>(&one) != 1) {
    for (std::size_t i = 0; i < column.size(); ++i)
      std::reverse(bytes + i * 8, bytes + (i + 1) * 8);
  }
  out.write(bytes, static_cast<std::streamsize>(column.size() * 8));
}

SEXP write_z_matrix(SEXP rsids, SEXP marker_file, SEXP paths, SEXP ids, SEXP ancestry,
                    SEXP directory, SEXP panels, SEXP panel_paths) {
  try {
    const auto size = static_cast<std::size_t>(XLENGTH(rsids));
    if (XLENGTH(panels) && size > std::numeric_limits<std::uint32_t>::max())
      throw std::runtime_error("The LD row map supports up to 2^32 - 1 Z rows.");
    Index index;
    index.reserve(size);
    for (std::size_t i = 0; i < size; ++i) {
      checkpoint(i);
      const auto id = STRING_ELT(rsids, i);
      if (id == NA_STRING || !valid_field(CHAR(id)) || !index.emplace(CHAR(id), i).second)
        throw std::runtime_error("GAlist$rsids must contain unique, nonempty marker IDs.");
    }
    const auto output = fs::u8path(CHAR(STRING_ELT(directory, 0)));
    Labels labels;
    std::vector<Marker> markers(size);
    write_markers(CHAR(STRING_ELT(marker_file, 0)), output / "markers.tsv", index, markers, labels);
    if (XLENGTH(panels) != 0) {
      auto catalog = output_file(output / "ld.tsv");
      catalog << "ancestry\tmap\tld_markers\tld_rows\tmatched\n";
      for (R_xlen_t i = 0; i < XLENGTH(panels); ++i) {
        const std::string panel = CHAR(STRING_ELT(panels, i));
        const std::string source = CHAR(STRING_ELT(panel_paths, i));
        const auto map = ld_map_name(panel);
        const auto count = write_ld_map(source, output / map, index, markers, labels);
        catalog << panel << '\t' << map << '\t' << source << '\t'
                << count.rows << '\t' << count.matched << '\n';
      }
      catalog.close();
    }
    auto z = output_file(output / "z.f64.bin");
    auto studies = output_file(output / "studies.tsv");
    studies << "study_id\tancestry\tobserved\tmissing\tld_map\n";
    std::vector<double> column(size);
    std::vector<bool> seen(size);
    for (R_xlen_t j = 0; j < XLENGTH(ids); ++j) {
      checkpoint(0);
      const std::string panel = CHAR(STRING_ELT(ancestry, j));
      const auto observed = read_study(CHAR(STRING_ELT(paths, j)), index, markers, labels,
                                       column, seen);
      write_column(z, column);
      studies << CHAR(STRING_ELT(ids, j)) << '\t' << panel
              << '\t' << observed << '\t' << size - observed << '\t';
      for (R_xlen_t i = 0; i < XLENGTH(panels); ++i) {
        if (panel == CHAR(STRING_ELT(panels, i))) {
          studies << ld_map_name(panel);
          break;
        }
      }
      studies << '\n';
    }
    z.close();
    studies.close();
    return R_NilValue;
  } catch (const std::exception& e) {
    Rf_error("%s", e.what());
  }
  return R_NilValue;
}

const R_CallMethodDef calls[] = {
    {"C_gact_write_z", reinterpret_cast<DL_FUNC>(&write_z_matrix), 8},
    {"C_gact_recode_sumstat", reinterpret_cast<DL_FUNC>(&recode_sumstat), 5},
    {nullptr, nullptr, 0}
};
}  // namespace

extern "C" void attribute_visible R_init_gact(DllInfo* dll) {
  R_registerRoutines(dll, nullptr, calls, nullptr, nullptr);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
