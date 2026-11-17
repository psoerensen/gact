// Streaming harmonization of verified GWAS summary statistics.
#include <zlib.h>
#ifndef _WIN32
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#endif

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <limits>
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
    return result.ec == std::errc() && result.ptr == field.data() + field.size()
        ? value : missing_integer;
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
    char* end;
    double result = std::strtod(value.data(), &end);
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

struct Header { uint64_t signature, size, pool_size; };
struct Marker {
    Key key;
    double eaf;
    uint32_t id;
    uint32_t shared;
};
struct Lookup { Key key; uint32_t row; uint32_t unused; };
static_assert(sizeof(Header) == 24 && sizeof(Marker) == 24 && sizeof(Lookup) == 16);

void index(const char* markers, const char* cpra, const char* output) {
    Gzip input(markers, "rb");
    Gzip keys(cpra, "rb");
    string line, cpra_line, pool;
    std::vector<Marker> catalog;
    std::vector<Lookup> lookup;
    require(input.read(line) && line == "rsids,chr,pos,ea,nea,eaf", "Unsupported marker header");
    bool has_key = keys.read(cpra_line);
    if (has_key && cpra_line == "cpra") has_key = keys.read(cpra_line);
    while (input.read(line)) {
        auto col = fields<6>(line, ',');
        int chr = chromosome(col[1]), pos = integer(col[2]);
        int effect = allele(col[3]), other = allele(col[4]);
        require(chr > 0 && chr <= 25 && pos > 0 && effect >= 0 && other >= 0,
                "Catalog must contain single-base alleles on chromosomes 1–25");
        string expected = string(col[1]) + "_" + string(col[2]) + "_" + string(col[3]) + "_" + string(col[4]);
        require(has_key && cpra_line == expected, "CPRA and marker catalog disagree");
        require(!absent(col[0]) && pool.size() + col[0].size() < UINT32_MAX &&
                catalog.size() < UINT32_MAX, "Invalid or oversized marker ID table");
        Key variant = key(chr, pos, effect, other);
        lookup.push_back({variant, static_cast<uint32_t>(catalog.size()), 0});
        catalog.push_back({variant, real(col[5]), static_cast<uint32_t>(pool.size()), 0});
        pool.append(col[0]);
        pool.push_back('\0');
        has_key = keys.read(cpra_line);
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

class Catalog {
#ifdef _WIN32
    std::vector<uint64_t> storage;
#endif
    void* mapping;
    size_t bytes;
    const Lookup* lookup;
    const char* pool;
public:
    const Marker* marker;
    size_t size;
    explicit Catalog(const char* path) {
#ifdef _WIN32
        std::ifstream in(path, std::ios::binary | std::ios::ate);
        require(in.good() && in.tellg() >= static_cast<std::streamoff>(sizeof(Header)),
                "Invalid marker index size");
        bytes = static_cast<size_t>(in.tellg());
        storage.resize((bytes + sizeof(uint64_t) - 1) / sizeof(uint64_t));
        in.seekg(0);
        in.read(reinterpret_cast<char*>(storage.data()), bytes);
        require(in.good(), "Cannot read marker index");
        mapping = storage.data();
#else
        int fd = open(path, O_RDONLY);
        require(fd >= 0, "Cannot open marker index");
        struct stat info;
        if (fstat(fd, &info) != 0 || info.st_size < static_cast<off_t>(sizeof(Header))) {
            close(fd);
            throw std::runtime_error("Invalid marker index size");
        }
        bytes = info.st_size;
        mapping = mmap(nullptr, bytes, PROT_READ, MAP_SHARED, fd, 0);
        close(fd);
        require(mapping != MAP_FAILED, "Cannot memory-map marker index");
#endif
        auto h = static_cast<const Header*>(mapping);
        require(h->signature == magic && h->size < UINT32_MAX && h->pool_size < UINT32_MAX &&
                bytes == sizeof(Header) + h->size * (sizeof(Marker) + sizeof(Lookup)) + h->pool_size,
                "Invalid marker index");
        size = h->size;
        marker = reinterpret_cast<const Marker*>(h + 1);
        lookup = reinterpret_cast<const Lookup*>(marker + size);
        pool = reinterpret_cast<const char*>(lookup + size);
    }
    ~Catalog() {
#ifndef _WIN32
        munmap(mapping, bytes);
#endif
    }
    Catalog(const Catalog&) = delete;
    const Lookup* find(Key variant) const {
        auto end = lookup + size;
        auto hit = std::lower_bound(lookup, end, variant,
            [](const Lookup& row, Key target) { return row.key < target; });
        return hit != end && hit->key == variant ? hit : nullptr;
    }
    bool duplicate(const Lookup* hit) const {
        return hit && hit + 1 < lookup + size && hit[1].key == hit->key;
    }
    const char* id(uint32_t row) const { return pool + marker[row].id; }
};

enum Drop { Missing, Unsplit, Strand, Unmatched, Ambiguous, Shared, Conflict, Redundant, Maf, Drops };
enum State : uint32_t { Empty, Direct, Flipped, Conflicting };

bool equal(double a, double b) { return a == b || (std::isnan(a) && std::isnan(b)); }

enum Column { Chr, Pos, EA, NEA, Beta, SE, P, EAF, Columns };
const std::array<string_view, Columns> canonical = {"chr","pos","ea","nea","b","seb","p","eaf"};

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
    bool forward;
    Schema(SEXP columns, SEXP mapping, bool strand) : forward(strand) {
        for (R_xlen_t i = 0; i < XLENGTH(columns); ++i) {
            if (i) header.push_back('\t');
            header += CHAR(STRING_ELT(columns, i));
        }
        split(header, names);
        require(XLENGTH(mapping) == Columns, "Invalid native mapping");
        for (size_t i = 0; i < field.size(); ++i) {
            auto name = STRING_ELT(mapping, i);
            auto at = name == NA_STRING ? names.end() :
                std::find(names.begin(), names.end(), CHAR(name));
            require(at != names.end() || i == EAF, "Required column absent from native mapping");
            field[i] = at - names.begin();
        }
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
            j == schema.field[P] || j == schema.field[EAF] || schema.names[j] == "N_sample_size") {
            auto a = real(row[j]), b = real(old[j]);
            if (a != b) return a < b;
        } else {
            return row[j] < old[j];
        }
    }
    return false;
}

enum Output { Raw, MarkerID, Chromosome, Position, Effect, Other, Coefficient, Error, Probability, Frequency };

void process(const char* index_path, const char* source, const char* output, const char* qc_path,
             const Schema& schema) {
    Catalog catalog(index_path);
    Gzip input(source, "rb");
    string line;
    require(input.read(line) && line == schema.header, "Input header differs from verified recipe");
    struct Estimate { double beta, se, p, eaf; State state; uint32_t copies; };
    std::vector<Estimate> values(catalog.size);
    std::vector<string> original(catalog.size);
    std::vector<string_view> columns, previous;
    std::array<uint64_t, Drops> drops{};
    uint64_t rows = 0, retained = 0, flipped = 0, mapped = 0;
    std::array<bool, EAF> present{};
    bool supplied_eaf = schema.field[EAF] < schema.names.size();
    while (input.read(line)) {
        split(line, columns);
        require(columns.size() == schema.names.size(), "Input row width differs from verified header");
        auto value = [&](Column j) { return columns[schema.field[j]]; };
        int chr = chromosome(value(Chr)), pos = integer(value(Pos));
        double beta = real(value(Beta)), se = real(value(SE)), p = real(value(P));
        double eaf = supplied_eaf ? real(value(EAF)) : missing;
        require((std::isnan(se) || se > 0) && (std::isnan(p) || (p >= 0 && p <= 1)) &&
                (std::isnan(eaf) || (eaf >= 0 && eaf <= 1)), "SE, p-value or EAF outside the allowed range");
        bool observed[] = {chr != missing_integer, pos != missing_integer,
                           !absent(value(EA)), !absent(value(NEA)), !std::isnan(beta),
                           !std::isnan(se), !std::isnan(p)};
        for (size_t j = 0; j < present.size(); ++j) present[j] = present[j] || observed[j];
        ++rows;
        if (!(observed[Chr] && observed[Pos] && observed[EA] && observed[NEA] && observed[P])) {
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
        auto direct = catalog.find(key(chr, pos, ea, nea));
        auto swap = catalog.find(key(chr, pos, nea, ea));
        if ((direct && swap) || catalog.duplicate(direct) || catalog.duplicate(swap)) {
            ++drops[Ambiguous];
            continue;
        }
        auto match = direct ? direct : swap;
        if (!match) { ++drops[Unmatched]; continue; }
        auto& marker = catalog.marker[match->row];
        if (marker.shared) { ++drops[Shared]; continue; }
        ++mapped;
        if (p == 0) p = std::numeric_limits<double>::min();
        if (!supplied_eaf) eaf = marker.eaf;
        else if (!direct) eaf = 1 - eaf;
        Estimate candidate{direct ? beta : -beta, se, p, eaf, direct ? Direct : Flipped, 0};
        auto& best = values[match->row];
        auto& raw = original[match->row];
        auto copies = best.copies + 1;
        if (best.state == Empty) {
            best = candidate;
            raw = line;
        } else if (best.state != Conflicting) {
            if (!equal(best.beta, candidate.beta) || !equal(best.se, se) || best.p != p ||
                !equal(std::nearbyint(best.eaf * 1e12), std::nearbyint(eaf * 1e12))) {
                best.state = Conflicting;
                raw.clear();
            } else if ((best.state == Flipped && candidate.state == Direct) ||
                       (best.state == candidate.state && precedes(schema, columns, raw, previous))) {
                best = candidate;
                raw = line;
            }
        }
        best.copies = copies;
    }
    require(std::all_of(present.begin(), present.end(), [](bool seen) { return seen; }),
            "Required summary statistic is wholly missing");

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
        Chromosome, Position, Effect, Other, Coefficient, Error, Probability, Frequency
    };
    for (size_t j = 0; j < EAF; ++j) {
        if (j == Beta && supplied_eaf) add("eaf", Frequency);
        add(canonical[j], outputs[j]);
    }
    add("rsids", MarkerID);
    if (!supplied_eaf) add("eaf", Frequency);

    Gzip out(output, "wb6");
    string buffer;
    for (auto field : header) {
        if (!buffer.empty()) buffer.push_back(',');
        csv(buffer, field);
    }
    buffer.push_back('\n');
    constexpr char bases[] = "ACGT";
    for (uint32_t i = 0; i < catalog.size; ++i) {
        const auto& v = values[i];
        if (v.state == Empty) continue;
        if (v.state == Conflicting) { drops[Conflict] += v.copies; continue; }
        drops[Redundant] += v.copies - 1;
        const auto& marker = catalog.marker[i];
        if (std::abs(marker.eaf - v.eaf) > 0.5) { ++drops[Maf]; continue; }
        ++retained;
        flipped += v.state == Flipped;
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
            }
        }
        buffer.push_back('\n');
        if (buffer.size() >= 1 << 20) { out.write(buffer); buffer.clear(); }
    }
    if (!buffer.empty()) out.write(buffer);
    out.finish();
    auto dropped = std::accumulate(drops.begin(), drops.end(), uint64_t{0});
    require(retained && rows == retained + dropped, "No retained variants or inconsistent QC counters");
    constexpr const char* reasons[] = {"missing_required", "unsplit_multiallelic",
        "unresolved_strand", "unmatched", "ambiguous_catalog_match", "shared_catalog_rsid",
        "conflicting_duplicate", "redundant_duplicate", "eaf_mismatch"};
    std::ofstream qc(qc_path);
    qc << "{\"strand\":\"" << (schema.forward ? "forward" : "unknown")
       << "\",\"n_input\":" << rows << ",\"n_mapped\":" << mapped << ",\"n_final\":" << retained
       << ",\"n_dropped\":" << dropped << ",\"n_flipped\":" << flipped
       << ",\"n_removed_maf\":" << drops[Maf] << ",\"dropped_by_reason\":{";
    for (size_t i = 0; i < Drops; ++i) {
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

extern "C" SEXP gact_gwas_stage(SEXP index_path, SEXP input, SEXP output, SEXP qc,
                                SEXP headers, SEXP mapping, SEXP forward) {
    try {
        const Schema schema(headers, mapping, Rf_asLogical(forward) == TRUE);
        process(path(index_path), path(input), path(output), path(qc), schema);
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
