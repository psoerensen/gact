// Streaming harmonization of verified GWAS summary statistics.
#include <zlib.h>
#ifndef _WIN32
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/file.h>
#include <unistd.h>
#else
#define NOMINMAX
#include <windows.h>
#endif

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <limits>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#define R_NO_REMAP
#include <Rinternals.h>

namespace {
using std::string;
using std::string_view;
using Key = uint64_t;
constexpr uint64_t magic = 0x4741435450414e31ULL;
constexpr double missing = std::numeric_limits<double>::quiet_NaN();
constexpr int missing_integer = std::numeric_limits<int>::min();

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

class Gzip {
    gzFile file;
public:
    Gzip(const char* path, const char* mode) : file(gzopen(path, mode)) {
        require(file != nullptr, "Cannot open gzip file");
        require(gzbuffer(file, 1 << 20) == 0, "Cannot allocate gzip buffer");
    }
    ~Gzip() { if (file) gzclose(file); }
    bool read(string& line, bool prefix = false) {
        line.clear();
        char buffer[4096];
        while (gzgets(file, buffer, sizeof buffer)) {
            require(buffer[0] != '\0', "Unexpected NUL in gzip text");
            line.append(buffer);
            if (line.back() == '\n') {
                line.pop_back();
                if (!line.empty() && line.back() == '\r') line.pop_back();
                return true;
            }
        }
        int error;
        gzerror(file, &error);
        if (prefix && (error == Z_BUF_ERROR || error == Z_OK || error == Z_STREAM_END)) return false;
        require(error == Z_OK || error == Z_STREAM_END, "Invalid or truncated gzip input");
        return !line.empty();
    }
    void write(string_view value) {
        require(gzwrite(file, value.data(), value.size()) == static_cast<int>(value.size()),
                "Cannot write gzip output");
    }
    void finish() {
        int result = gzclose(file);
        file = nullptr;
        require(result == Z_OK, "Cannot close gzip output");
    }
    Gzip(const Gzip&) = delete;
};

template <size_t N>
std::array<string_view, N> fields(string_view line, char delim) {
    std::array<string_view, N> result;
    for (size_t i = 0; i + 1 < N; ++i) {
        auto end = line.find(delim);
        require(end != string_view::npos, "Input row has too few fields");
        result[i] = line.substr(0, end);
        line.remove_prefix(end + 1);
    }
    require(line.find(delim) == string_view::npos, "Input row has too many fields");
    result.back() = line;
    return result;
}

int integer(string_view field) {
    int value;
    auto result = std::from_chars(field.data(), field.data() + field.size(), value);
    if (result.ec == std::errc() && result.ptr == field.data() + field.size()) return value;
    // fwrite can serialize integer-valued numeric coordinates as 1e+06.
    if (field.empty() || field.find_first_not_of("0123456789+-.eE") != string_view::npos)
        return missing_integer;
    char* end;
    double numeric = std::strtod(field.data(), &end);
    return end == field.data() + field.size() && std::isfinite(numeric) &&
        numeric >= std::numeric_limits<int>::min() && numeric <= std::numeric_limits<int>::max() &&
        std::trunc(numeric) == numeric ? static_cast<int>(numeric) : missing_integer;
}

template <size_t N>
std::array<string_view, N> csv_fields(string_view line, char delimiter = ',') {
    std::array<string_view, N> result;
    for (size_t i = 0; i < N; ++i) {
        size_t end;
        if (!line.empty() && line.front() == '"') {
            end = line.find('"', 1);
            require(end != string_view::npos, "Unclosed quoted catalog field");
            result[i] = line.substr(1, end - 1);
            ++end;
            require(end == line.size() || line[end] == delimiter, "Invalid quoted catalog field");
        } else {
            end = line.find(delimiter);
            if (end == string_view::npos) end = line.size();
            result[i] = line.substr(0, end);
            require(result[i].find('"') == string_view::npos, "Invalid catalog quote");
        }
        if (i + 1 == N) require(end == line.size(), "Catalog row has too many fields");
        else {
            require(end < line.size(), "Catalog row has too few fields");
            line.remove_prefix(end + 1);
        }
    }
    return result;
}

int chromosome(string_view field) {
    if (field.size() >= 3 && (field[0] == 'c' || field[0] == 'C') &&
        (field[1] == 'h' || field[1] == 'H') && (field[2] == 'r' || field[2] == 'R'))
        field.remove_prefix(3);
    if (field == "X" || field == "x") return 23;
    if (field == "Y" || field == "y") return 24;
    if (field == "M" || field == "m" || field == "MT" || field == "mt") return 25;
    return integer(field);
}

bool absent(string_view value) { return value.empty() || value == "NA"; }

double real(string_view value) {
    if (absent(value)) return missing;
    double result;
    const auto parsed = std::from_chars(value.data(), value.data() + value.size(), result);
    if (parsed.ec == std::errc() && parsed.ptr == value.data() + value.size() &&
        std::isfinite(result)) return result;
    // Preserve the original parser's handling of leading '+', underflow and errors.
    char* end;
    result = std::strtod(value.data(), &end);
    require(end == value.data() + value.size() && std::isfinite(result),
            "Non-finite or invalid summary statistic");
    return result;
}

int allele(string_view value) {
    if (value.size() != 1) return -1;
    switch (value[0]) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return -1;
    }
}

Key key(int chr, int pos, int effect, int other) {
    return (Key(chr) << 35) | (Key(pos) << 4) | (effect << 2) | other;
}

Key swap_alleles(Key variant) {
    return (variant & ~Key{15}) | ((variant & 3) << 2) | ((variant >> 2) & 3);
}

struct Header { uint64_t signature, size, pool_size; };
constexpr uint64_t lift_magic = 0x474143544C465431ULL;
struct LiftHeader { uint64_t signature, size; };
struct Lift { Key source, target; };
struct Marker {
    Key key;
    double eaf;
    uint32_t id;
    uint32_t shared;
};
struct Lookup { Key key; uint32_t row; uint32_t unused; };
static_assert(sizeof(Header) == 24 && sizeof(Marker) == 24 && sizeof(Lookup) == 16 &&
              sizeof(LiftHeader) == 16 && sizeof(Lift) == 16);

std::vector<Lift> read_liftover(SEXP paths) {
    std::vector<Lift> entries;
    for (R_xlen_t i = 0; i < XLENGTH(paths); ++i) {
        Gzip input(CHAR(STRING_ELT(paths, i)), "rb");
        string line;
        require(input.read(line) && line == "chr\tpos\tref\talt\ttarget_pos\tnea\tea\treversed\trsids",
                "Invalid liftover header");
        while (input.read(line)) {
            auto col = fields<9>(line, '\t');
            int chr = chromosome(col[0]), pos = integer(col[1]), target = integer(col[4]);
            int ea = allele(col[3]), nea = allele(col[2]);
            int target_ea = allele(col[6]), target_nea = allele(col[5]);
            require(chr >= 1 && chr <= 22 && pos > 0 && target > 0 &&
                    ea >= 0 && nea >= 0 && ea != nea && target_ea >= 0 && target_nea >= 0 &&
                    ((col[7] == "TRUE" && target_ea == 3 - ea && target_nea == 3 - nea) ||
                     (col[7] == "FALSE" && target_ea == ea && target_nea == nea)),
                    "Invalid liftover mapping");
            entries.push_back({key(chr, pos, ea, nea), key(chr, target, target_ea, target_nea)});
        }
    }
    std::sort(entries.begin(), entries.end(), [](const auto& a, const auto& b) { return a.source < b.source; });
    auto duplicate = std::adjacent_find(entries.begin(), entries.end(),
        [](const auto& a, const auto& b) { return a.source == b.source; });
    require(duplicate == entries.end(), "Duplicate source key in liftover map");
    return entries;
}

void write_liftover(SEXP paths, const char* output) {
    auto entries = read_liftover(paths);
    require(!entries.empty(), "Empty liftover map");
    std::ofstream out(output, std::ios::binary);
    LiftHeader header{lift_magic, entries.size()};
    out.write(reinterpret_cast<const char*>(&header), sizeof header);
    out.write(reinterpret_cast<const char*>(entries.data()), sizeof(Lift) * entries.size());
    out.close();
    require(out.good(), "Cannot write liftover index");
}

void index(const char* markers, const char* cpra, const char* output) {
    Gzip input(markers, "rb");
    Gzip keys(cpra, "rb");
    string line, cpra_line, pool;
    std::vector<Marker> catalog;
    std::vector<Lookup> lookup;
    require(input.read(line), "Missing marker header");
    if (line.compare(0, 3, "\xef\xbb\xbf") == 0) line.erase(0, 3);
    const char delimiter = line.find('\t') != string::npos ? '\t' : ',';
    const std::array<string_view, 6> expected_header = {"rsids", "chr", "pos", "ea", "nea", "eaf"};
    require(csv_fields<6>(line, delimiter) == expected_header, "Unsupported marker header");
    bool has_key = keys.read(cpra_line);
    if (has_key && cpra_line.compare(0, 3, "\xef\xbb\xbf") == 0) cpra_line.erase(0, 3);
    if (has_key && csv_fields<1>(cpra_line)[0] == "cpra") has_key = keys.read(cpra_line);
    uint64_t row = 0;
    while (input.read(line)) {
      ++row;
      try {
        auto col = csv_fields<6>(line, delimiter);
        int chr = chromosome(col[1]), pos = integer(col[2]);
        int effect = allele(col[3]), other = allele(col[4]);
        if (!(chr > 0 && chr <= 25 && pos > 0 && effect >= 0 && other >= 0 && effect != other))
            throw std::runtime_error("Invalid catalog coordinate or alleles: chr=" + string(col[1]) +
                ", pos=" + string(col[2]) + ", ea=" + string(col[3]) + ", nea=" + string(col[4]));
        require(has_key, "CPRA has fewer rows than the marker catalog");
        auto cpra_fields = fields<4>(csv_fields<1>(cpra_line)[0], '_');
        require(chromosome(cpra_fields[0]) == chr && integer(cpra_fields[1]) == pos &&
                allele(cpra_fields[2]) == effect && allele(cpra_fields[3]) == other,
                "CPRA and marker catalog disagree");
        require(!absent(col[0]) && pool.size() + col[0].size() < UINT32_MAX &&
                catalog.size() < UINT32_MAX, "Invalid or oversized marker ID table");
        Key variant = key(chr, pos, effect, other);
        lookup.push_back({variant, static_cast<uint32_t>(catalog.size()), 0});
        catalog.push_back({variant, real(col[5]), static_cast<uint32_t>(pool.size()), 0});
        pool.append(col[0]);
        pool.push_back('\0');
        has_key = keys.read(cpra_line);
      } catch (const std::exception& error) {
        throw std::runtime_error("Marker catalog row " + std::to_string(row) + ": " + error.what());
      }
    }
    require(!has_key && !catalog.empty(), "CPRA and catalog row counts disagree");
    std::vector<uint32_t> ids(catalog.size());
    std::iota(ids.begin(), ids.end(), 0);
    auto id = [&](uint32_t i) { return string_view(pool.c_str() + catalog[i].id); };
    std::sort(ids.begin(), ids.end(), [&](uint32_t a, uint32_t b) { return id(a) < id(b); });
    for (size_t i = 1; i < ids.size(); ++i) {
        if (id(ids[i - 1]) == id(ids[i])) catalog[ids[i - 1]].shared = catalog[ids[i]].shared = 1;
    }
    std::sort(lookup.begin(), lookup.end(), [](const auto& a, const auto& b) { return a.key < b.key; });
    std::ofstream out(output, std::ios::binary);
    Header header{magic, catalog.size(), pool.size()};
    out.write(reinterpret_cast<const char*>(&header), sizeof header);
    out.write(reinterpret_cast<const char*>(catalog.data()), sizeof(Marker) * catalog.size());
    out.write(reinterpret_cast<const char*>(lookup.data()), sizeof(Lookup) * lookup.size());
    out.write(pool.data(), pool.size());
    out.close();
    require(out.good(), "Cannot write marker index");
}

class MappedFile {
#ifdef _WIN32
    std::vector<uint64_t> storage;
#endif
public:
    void* mapping;
    size_t bytes;
    MappedFile(const char* path, size_t min_size) : mapping(nullptr), bytes(0) {
        if (!path) return;
#ifdef _WIN32
        std::ifstream in(path, std::ios::binary | std::ios::ate);
        require(in.good() && in.tellg() >= static_cast<std::streamoff>(min_size),
                "Invalid alignment index size");
        bytes = static_cast<size_t>(in.tellg());
        storage.resize((bytes + sizeof(uint64_t) - 1) / sizeof(uint64_t));
        in.seekg(0);
        in.read(reinterpret_cast<char*>(storage.data()), bytes);
        require(in.good(), "Cannot read alignment index");
        mapping = storage.data();
#else
        int fd = open(path, O_RDONLY);
        require(fd >= 0, "Cannot open alignment index");
        struct stat info;
        if (fstat(fd, &info) != 0 || info.st_size < static_cast<off_t>(min_size)) {
            close(fd);
            throw std::runtime_error("Invalid alignment index size");
        }
        bytes = info.st_size;
        mapping = mmap(nullptr, bytes, PROT_READ, MAP_SHARED, fd, 0);
        close(fd);
        require(mapping != MAP_FAILED, "Cannot memory-map alignment index");
#endif
    }
    ~MappedFile() {
#ifndef _WIN32
        if (mapping) munmap(mapping, bytes);
#endif
    }
    MappedFile(const MappedFile&) = delete;
};

class Liftover {
    MappedFile file;
    const Lift* entries = nullptr;
    size_t size = 0;
public:
    explicit Liftover(const char* path) : file(path, sizeof(LiftHeader)) {
        if (!path) return;
        auto h = static_cast<const LiftHeader*>(file.mapping);
        require(h->signature == lift_magic && h->size && h->size < UINT32_MAX &&
                file.bytes == sizeof(LiftHeader) + h->size * sizeof(Lift),
                "Invalid liftover index");
        entries = reinterpret_cast<const Lift*>(h + 1);
        size = h->size;
    }
    Key find(Key source) const {
        auto end = entries + size;
        auto hit = std::lower_bound(entries, end, source,
            [](const Lift& row, Key target) { return row.source < target; });
        return hit != end && hit->source == source ? hit->target : 0;
    }
};

class Catalog {
    MappedFile file;
    const Lookup* lookup;
    const char* pool;
    size_t cursor = 0;
public:
    const Marker* marker;
    size_t size;
    explicit Catalog(const char* path) : file(path, sizeof(Header)) {
        auto h = static_cast<const Header*>(file.mapping);
        require(h->signature == magic && h->size < UINT32_MAX && h->pool_size < UINT32_MAX &&
                file.bytes == sizeof(Header) + h->size * (sizeof(Marker) + sizeof(Lookup)) + h->pool_size,
                "Invalid marker index");
        size = h->size;
        marker = reinterpret_cast<const Marker*>(h + 1);
        lookup = reinterpret_cast<const Lookup*>(marker + size);
        pool = reinterpret_cast<const char*>(lookup + size);
    }
    std::pair<const Lookup*, const Lookup*> orientations(Key variant) {
        // All allele pairs for a position are adjacent. Common sorted input
        // usually advances zero or one entries. Seek on jumps and unsorted rows.
        const Key start = variant & ~Key{15};
        if ((cursor && lookup[cursor - 1].key >= start) ||
            (size - cursor > 8 && lookup[cursor + 8].key < start)) {
            cursor = std::lower_bound(lookup, lookup + size, start,
                [](const Lookup& row, Key target) { return row.key < target; }) - lookup;
        } else {
            while (cursor < size && lookup[cursor].key < start) ++cursor;
        }
        const Key swapped = swap_alleles(variant);
        const Lookup* direct = nullptr;
        const Lookup* reverse = nullptr;
        for (size_t j = cursor; j < size && lookup[j].key <= (start | 15); ++j) {
            if (!direct && lookup[j].key == variant) direct = lookup + j;
            if (!reverse && lookup[j].key == swapped) reverse = lookup + j;
        }
        return {direct, reverse};
    }
    bool duplicate(const Lookup* hit) const {
        return hit && hit + 1 < lookup + size && hit[1].key == hit->key;
    }
    const char* id(uint32_t row) const { return pool + marker[row].id; }
};

enum Drop { Missing, Unsplit, Strand, Unmatched, Ambiguous, Shared, Conflict, Redundant, Maf,
            LiftUnmapped, LowConfidence, Info, Drops };
enum State : uint32_t { Empty, Direct, Flipped, Conflicting };

bool equal(double a, double b) { return a == b || (std::isnan(a) && std::isnan(b)); }

enum Column { Chr, Pos, EA, NEA, Beta, SE, P, EAF, N, Columns };
const std::array<string_view, Columns> canonical = {"chr","pos","ea","nea","b","seb","p","eaf","n"};

void split(string_view line, std::vector<string_view>& columns) {
    columns.clear();
    for (;;) {
        auto end = line.find('\t');
        columns.push_back(line.substr(0, end));
        if (end == string_view::npos) return;
        line.remove_prefix(end + 1);
    }
}

void csv(string& row, string_view value) {
    if (absent(value)) return;
    if (value.find_first_of(",\"\r\n") == string_view::npos) {
        row.append(value);
        return;
    }
    row.push_back('"');
    for (char ch : value) {
        if (ch == '"') row.push_back('"');
        row.push_back(ch);
    }
    row.push_back('"');
}

void number(string& row, double value) {
    if (std::isnan(value)) return;
    char buffer[64];
    auto formatted = std::to_chars(buffer, buffer + sizeof buffer, value);
    require(formatted.ec == std::errc(), "Cannot format estimate");
    row.append(buffer, formatted.ptr);
}

struct Schema {
    string header;
    std::vector<string_view> names;
    std::array<size_t, Columns> field;
    size_t low_confidence, info;
    double min_info;
    bool forward;
    Schema(SEXP columns, SEXP mapping, bool strand) : forward(strand) {
        for (R_xlen_t i = 0; i < XLENGTH(columns); ++i) {
            if (i) header.push_back('\t');
            header += CHAR(STRING_ELT(columns, i));
        }
        split(header, names);
        require(XLENGTH(mapping) == Columns || XLENGTH(mapping) == Columns - 1, "Invalid native mapping");
        for (size_t i = 0; i < field.size(); ++i) {
            auto name = i < static_cast<size_t>(XLENGTH(mapping)) ? STRING_ELT(mapping, i) : NA_STRING;
            auto at = name == NA_STRING ? names.end() :
                std::find(names.begin(), names.end(), CHAR(name));
            require(at != names.end() || i >= EAF, "Required column absent from native mapping");
            field[i] = at - names.begin();
        }
        auto quality = [&](const char* key) {
            SEXP column = Rf_getAttrib(mapping, Rf_install(key));
            if (column == R_NilValue) return names.size();
            auto at = std::find(names.begin(), names.end(), CHAR(STRING_ELT(column, 0)));
            require(at != names.end(), "Quality filter column absent from verified header");
            return static_cast<size_t>(at - names.begin());
        };
        low_confidence = quality("low_confidence");
        info = quality("info");
        min_info = Rf_asReal(Rf_getAttrib(mapping, Rf_install("min_info")));
    }
};

bool precedes(const Schema& schema, const std::vector<string_view>& row, const string& saved,
              std::vector<string_view>& old) {
    split(saved, old);
    for (size_t j = 0; j < row.size(); ++j) {
        if (row[j] == old[j] || j == schema.field[Chr] || j == schema.field[Pos]) continue;
        if (absent(row[j]) || absent(old[j])) {
            if (absent(row[j]) != absent(old[j])) return !absent(row[j]);
            continue;
        }
        if (j == schema.field[Beta] || j == schema.field[SE] ||
            j == schema.field[P] || j == schema.field[EAF] || j == schema.field[N]) {
            auto a = real(row[j]), b = real(old[j]);
            if (a != b) return a < b;
        } else {
            return row[j] < old[j];
        }
    }
    return false;
}

enum Output { Raw, MarkerID, Chromosome, Position, Effect, Other, Coefficient, Error, Probability, Frequency, SampleSize };

void process(const char* index_path, const char* source, const char* output, const char* qc_path,
             const Schema& schema, const char* lift_path, bool neglog10p, const char* z_path) {
    Catalog catalog(index_path);
    Liftover liftover(lift_path);
    Gzip input(source, "rb");
    string line;
    require(input.read(line) && line == schema.header, "Input header differs from verified recipe");
    struct Estimate { double beta, se, p, eaf, n; State state; uint32_t copies; };
    std::vector<Estimate> values(catalog.size);
    std::vector<string> original(output ? catalog.size : 0);
    std::vector<string_view> columns, previous;
    std::array<uint64_t, Drops> drops{};
    uint64_t rows = 0, retained = 0, flipped = 0, mapped = 0, lifted = 0, with_n = 0;
    bool supplied_eaf = schema.field[EAF] < schema.names.size();
    bool supplied_n = schema.field[N] < schema.names.size();
    while (input.read(line)) {
        split(line, columns);
        require(columns.size() == schema.names.size(), "Input row width differs from verified header");
        ++rows;
        if (schema.low_confidence != schema.names.size()) {
            auto flag = columns[schema.low_confidence];
            if (absent(flag) || flag == "true" || flag == "TRUE" || flag == "True" || flag == "1") {
                ++drops[LowConfidence];
                continue;
            }
            require(flag == "false" || flag == "FALSE" || flag == "False" || flag == "0",
                    "Invalid low-confidence flag");
        }
        if (schema.info != schema.names.size()) {
            double info = real(columns[schema.info]);
            require(std::isnan(info) || (info >= 0 && info <= 1), "INFO outside [0,1]");
            if (!std::isfinite(info) || info < schema.min_info) { ++drops[Info]; continue; }
        }
        auto value = [&](Column j) { return columns[schema.field[j]]; };
        int chr = chromosome(value(Chr)), pos = integer(value(Pos));
        double beta = real(value(Beta)), se = real(value(SE)), p = real(value(P));
        double eaf = supplied_eaf ? real(value(EAF)) : missing;
        double n = missing;
        if (supplied_n && !absent(value(N))) {
            char* end;
            n = std::strtod(value(N).data(), &end);
            require(end == value(N).data() + value(N).size() && std::isfinite(n) && n > 0,
                    "Per-variant n must be finite and strictly positive, or missing");
        }
        require((std::isnan(se) || se > 0) && (std::isnan(p) || (p >= 0 && (neglog10p || p <= 1))) &&
                (std::isnan(eaf) || (eaf >= 0 && eaf <= 1)), "SE, p-value or EAF outside the allowed range");
        if (chr == missing_integer || pos == missing_integer ||
            absent(value(EA)) || absent(value(NEA)) || std::isnan(p)) {
            ++drops[Missing];
            continue;
        }
        if (value(EA).find(',') != string_view::npos || value(NEA).find(',') != string_view::npos) {
            ++drops[Unsplit];
            continue;
        }
        int ea = allele(value(EA)), nea = allele(value(NEA));
        if (!schema.forward && ea >= 0 && nea >= 0 && ea + nea == 3) {
            ++drops[Strand];
            continue;
        }
        if (ea < 0 || nea < 0 || chr < 1 || chr > 25 || pos < 1) {
            ++drops[Unmatched];
            continue;
        }
        Key variant = key(chr, pos, ea, nea);
        if (lift_path) {
            auto lifted_direct = liftover.find(variant);
            auto lifted_swap = liftover.find(swap_alleles(variant));
            if (lifted_direct && lifted_swap) { ++drops[Ambiguous]; continue; }
            variant = lifted_direct ? lifted_direct : swap_alleles(lifted_swap);
            if (!variant) { ++drops[LiftUnmapped]; continue; }
            ++lifted;
        }
        auto [direct, swap] = catalog.orientations(variant);
        if ((direct && swap) || catalog.duplicate(direct) || catalog.duplicate(swap)) {
            ++drops[Ambiguous];
            continue;
        }
        auto match = direct ? direct : swap;
        if (!match) { ++drops[Unmatched]; continue; }
        auto& marker = catalog.marker[match->row];
        if (marker.shared) { ++drops[Shared]; continue; }
        ++mapped;
        if (neglog10p) p = std::pow(10.0, -p);
        if (p == 0) p = std::numeric_limits<double>::min();
        if (supplied_eaf && !direct) eaf = 1 - eaf;
        Estimate candidate{direct ? beta : -beta, se, p, eaf, n, direct ? Direct : Flipped, 0};
        auto& best = values[match->row];
        auto raw = output ? &original[match->row] : nullptr;
        auto copies = best.copies + 1;
        if (best.state == Empty) {
            best = candidate;
            if (raw) *raw = line;
        } else if (best.state != Conflicting) {
            if (!equal(best.beta, candidate.beta) || !equal(best.se, se) || best.p != p ||
                !equal(std::nearbyint(best.eaf * 1e12), std::nearbyint(eaf * 1e12)) || !equal(best.n, n)) {
                best.state = Conflicting;
                if (raw) raw->clear();
            } else if ((best.state == Flipped && candidate.state == Direct) ||
                       (raw && best.state == candidate.state && precedes(schema, columns, *raw, previous))) {
                best = candidate;
                if (raw) *raw = line;
            }
        }
        best.copies = copies;
    }

    std::vector<string_view> header = schema.names;
    std::vector<Output> output_fields(header.size(), Raw);
    auto add = [&](string_view name, Output value) {
        auto found = std::find(header.begin(), header.end(), name);
        if (found == header.end()) {
            output_fields.push_back(value);
            header.push_back(name);
        } else {
            output_fields[found - header.begin()] = value;
        }
    };
    add("marker", MarkerID);
    const std::array<Output, Columns> outputs = {
        Chromosome, Position, Effect, Other, Coefficient, Error, Probability, Frequency, SampleSize
    };
    for (size_t j = 0; j < EAF; ++j) {
        if (j == Beta && supplied_eaf) add("eaf", Frequency);
        add(canonical[j], outputs[j]);
    }
    add("rsids", MarkerID);
    if (!supplied_eaf) add("eaf", Frequency);
    add("n", SampleSize); // Keep missing per-variant n explicit even when study totals are known.

    std::unique_ptr<Gzip> out;
    if (output) out = std::make_unique<Gzip>(output, "wb1");
    std::ofstream z;
    if (z_path) {
        z.exceptions(std::ios::failbit | std::ios::badbit);
        z.open(z_path, std::ios::binary);
    }
    static_assert(sizeof(double) == 8 && std::numeric_limits<double>::is_iec559);
    const uint16_t one = 1;
    uint64_t observed_z = 0;
    auto write_z = [&](double value) {
        if (!z_path) return;
        auto bytes = reinterpret_cast<char*>(&value);
        if (*reinterpret_cast<const char*>(&one) == 0) std::reverse(bytes, bytes + sizeof(value));
        z.write(bytes, sizeof(value));
    };
    string buffer;
    for (auto field : header) {
        if (!buffer.empty()) buffer.push_back(',');
        csv(buffer, field);
    }
    buffer.push_back('\n');
    constexpr char bases[] = "ACGT";
    for (uint32_t i = 0; i < catalog.size; ++i) {
        const auto& v = values[i];
        if (v.state == Empty) { write_z(missing); continue; }
        if (v.state == Conflicting) { drops[Conflict] += v.copies; write_z(missing); continue; }
        drops[Redundant] += v.copies - 1;
        const auto& marker = catalog.marker[i];
        ++retained;
        with_n += !std::isnan(v.n);
        flipped += v.state == Flipped;
        const double score = v.beta == 0 && v.se > 0 ? 0.0 : v.beta / v.se;
        observed_z += std::isfinite(score);
        write_z(std::isfinite(score) ? score : missing);
        if (!out) continue;
        split(original[i], columns);
        for (size_t j = 0; j < output_fields.size(); ++j) {
            if (j) buffer.push_back(',');
            switch (output_fields[j]) {
                case Raw: csv(buffer, columns[j]); break;
                case MarkerID: csv(buffer, catalog.id(i)); break;
                case Chromosome: buffer += std::to_string(marker.key >> 35); break;
                case Position: buffer += std::to_string((marker.key >> 4) & 0x7fffffff); break;
                case Effect: buffer.push_back(bases[(marker.key >> 2) & 3]); break;
                case Other: buffer.push_back(bases[marker.key & 3]); break;
                case Coefficient: number(buffer, v.beta); break;
                case Error: number(buffer, v.se); break;
                case Probability: number(buffer, v.p); break;
                case Frequency: number(buffer, v.eaf); break;
                case SampleSize: number(buffer, v.n); break;
            }
        }
        buffer.push_back('\n');
        if (buffer.size() >= 1 << 20) { out->write(buffer); buffer.clear(); }
    }
    if (out) {
        if (!buffer.empty()) out->write(buffer);
        out->finish();
    }
    if (z_path) z.close();
    auto dropped = std::accumulate(drops.begin(), drops.end(), uint64_t{0});
    require(rows == retained + dropped, "Inconsistent QC counters");
    constexpr const char* reasons[] = {"missing_required", "unsplit_multiallelic",
        "unresolved_strand", "unmatched", "ambiguous_catalog_match", "shared_catalog_rsid",
        "conflicting_duplicate", "redundant_duplicate", "eaf_mismatch", "liftover_unmapped",
        "provider_low_confidence", "provider_info"};
    std::ofstream qc(qc_path);
    qc << "{\"strand\":\"" << (schema.forward ? "forward" : "unknown")
       << "\",\"n_input\":" << rows << ",\"n_mapped\":" << mapped << ",\"n_final\":" << retained
       << ",\"n_dropped\":" << dropped << ",\"n_flipped\":" << flipped
       << ",\"n_removed_maf\":" << drops[Maf]
       << ",\"n_with_sample_size\":" << with_n << ",\"n_missing_sample_size\":" << retained - with_n;
    if (lift_path) qc << ",\"n_lifted\":" << lifted;
    qc << ",\"n_observed_z\":" << observed_z
       << ",\"n_unusable_effect\":" << retained - observed_z
       << ",\"n_missing_backbone\":" << catalog.size - retained
       << ",\"n_backbone\":" << catalog.size;
    qc << ",\"dropped_by_reason\":{";
    for (size_t i = 0; i < Drops; ++i) {
        if (i == LiftUnmapped && !lift_path) continue;
        if (i) qc << ',';
        qc << '"' << reasons[i] << "\":" << drops[i];
    }
    qc << "}}\n";
    qc.close();
    require(qc.good(), "Cannot write QC");
}

const char* path(SEXP arg) { return CHAR(STRING_ELT(arg, 0)); }
}  // namespace

extern "C" SEXP gact_gwas_index(SEXP markers, SEXP cpra, SEXP output) {
    try { index(path(markers), path(cpra), path(output)); }
    catch (const std::exception& e) { Rf_error("%s", e.what()); }
    return R_NilValue;
}

extern "C" SEXP gact_gwas_check_ids(SEXP index_path, SEXP ids) {
    try {
        Catalog catalog(path(index_path));
        require(XLENGTH(ids) == static_cast<R_xlen_t>(catalog.size), "Marker catalog and GAlist disagree");
        for (size_t i = 0; i < catalog.size; ++i) {
            SEXP id = STRING_ELT(ids, i);
            require(id != NA_STRING && std::strcmp(CHAR(id), catalog.id(i)) == 0,
                    "Marker catalog and GAlist disagree");
        }
    } catch (const std::exception& e) { Rf_error("%s", e.what()); }
    return R_NilValue;
}

extern "C" SEXP gact_gwas_liftover_index(SEXP paths, SEXP output) {
    try { write_liftover(paths, path(output)); }
    catch (const std::exception& e) { Rf_error("%s", e.what()); }
    return R_NilValue;
}

extern "C" SEXP gact_gwas_stage(SEXP index_path, SEXP input, SEXP output, SEXP qc,
                                SEXP headers, SEXP mapping, SEXP forward, SEXP liftover,
                                SEXP neglog10p) {
    try {
        const Schema schema(headers, mapping, Rf_asLogical(forward) == TRUE);
        const char* full = STRING_ELT(output, 0) == NA_STRING ? nullptr : path(output);
        const char* z = XLENGTH(output) == 2 ? CHAR(STRING_ELT(output, 1)) : nullptr;
        require(full || z, "No GWAS output requested");
        process(path(index_path), path(input), full, path(qc), schema,
                liftover == R_NilValue ? nullptr : path(liftover), Rf_asLogical(neglog10p) == TRUE, z);
    } catch (const std::exception& e) { Rf_error("%s", e.what()); }
    return R_NilValue;
}

extern "C" SEXP gact_gwas_preview(SEXP file, SEXP limit, SEXP prefix) {
    try {
        Gzip input(path(file), "rb");
        string line, text;
        int count = Rf_asInteger(limit);
        const bool partial = Rf_asLogical(prefix) == TRUE;
        for (int i = 0; i <= count && input.read(line, partial); ++i) {
            text += line;
            text.push_back('\n');
        }
        return Rf_mkString(text.c_str());
    } catch (const std::exception& e) { Rf_error("%s", e.what()); }
    return R_NilValue;
}

namespace {
#ifdef _WIN32
using Guard = HANDLE;
#else
using Guard = int;
#endif

void close_guard(SEXP token) {
    auto guard = static_cast<Guard*>(R_ExternalPtrAddr(token));
    if (!guard) return;
#ifdef _WIN32
    CloseHandle(*guard);
#else
    close(*guard);
#endif
    delete guard;
    R_ClearExternalPtr(token);
}
}

extern "C" SEXP gact_gwas_lock(SEXP name) {
    SEXP token = PROTECT(R_MakeExternalPtr(nullptr, R_NilValue, R_NilValue));
    R_RegisterCFinalizerEx(token, close_guard, TRUE);
#ifdef _WIN32
    auto handle = CreateFileA(path(name), GENERIC_READ | GENERIC_WRITE,
        FILE_SHARE_READ | FILE_SHARE_WRITE, nullptr, OPEN_ALWAYS, FILE_ATTRIBUTE_NORMAL, nullptr);
    if (handle == INVALID_HANDLE_VALUE) Rf_error("Cannot open the GWAS import guard");
    OVERLAPPED region{};
    if (!LockFileEx(handle, LOCKFILE_EXCLUSIVE_LOCK | LOCKFILE_FAIL_IMMEDIATELY, 0, 1, 0, &region)) {
        CloseHandle(handle);
        Rf_error("Another GWAS import holds the database lock, or file locking is unavailable");
    }
#else
    int handle = open(path(name), O_RDWR | O_CREAT | O_CLOEXEC, 0600);
    if (handle < 0) Rf_error("Cannot open the GWAS import guard");
    if (flock(handle, LOCK_EX | LOCK_NB) != 0) {
        close(handle);
        Rf_error("Another GWAS import holds the database lock, or file locking is unavailable");
    }
#endif
    R_SetExternalPtrAddr(token, new Guard(handle));
    UNPROTECT(1);
    return token;
}

extern "C" SEXP gact_gwas_unlock(SEXP token) {
    close_guard(token);
    return R_NilValue;
}

extern "C" SEXP gact_gwas_sync(SEXP paths, SEXP directory) {
    const bool is_dir = Rf_asLogical(directory) == TRUE;
    for (R_xlen_t i = 0; i < XLENGTH(paths); ++i) {
        const char* name = CHAR(STRING_ELT(paths, i));
#ifdef _WIN32
        if (is_dir) continue;
        auto file = CreateFileA(name, GENERIC_WRITE, FILE_SHARE_READ | FILE_SHARE_WRITE, nullptr,
            OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, nullptr);
        if (file == INVALID_HANDLE_VALUE) Rf_error("Cannot sync GWAS transaction file: %s", name);
        int ok = FlushFileBuffers(file);
        CloseHandle(file);
#else
        int file = open(name, O_RDONLY | O_CLOEXEC);
        if (file < 0) Rf_error("Cannot sync GWAS transaction path: %s", name);
        bool ok = fsync(file) == 0;
        close(file);
#endif
        if (!ok) Rf_error("Cannot sync GWAS transaction path: %s", name);
    }
    return R_NilValue;
}
