#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include <R_ext/Error.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <climits>
#include <cstdint>
#include <cstring>
#include <functional>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
#include "protected-total.h"

namespace {
// Every allocating R API call runs below an unwind context. Its cleanup
// converts an R long jump into a C++ exception; the .Call boundary resumes the
// original R condition after all native containers and PROTECT guards unwind.
struct RJump { SEXP token; };
static SEXP active_unwind_token = R_NilValue;
void r_unwind_cleanup(void *data, Rboolean jump) {
    if (jump) throw RJump{static_cast<SEXP>(data)};
}
template <typename Function>
auto r_api(Function function) -> decltype(function()) {
    typedef decltype(function()) Result;
    Result result;
    struct Context { Function *function; Result *result; } context = {&function, &result};
    R_UnwindProtect([](void *data) -> SEXP {
        Context *context = static_cast<Context *>(data);
        *context->result = (*context->function)();
        return R_NilValue;
    }, &context, r_unwind_cleanup, active_unwind_token, active_unwind_token);
    return result;
}
#define R_API(expression) r_api([&]() { return (expression); })
struct UnwindToken {
    SEXP previous;
    explicit UnwindToken(SEXP token) : previous(active_unwind_token) { active_unwind_token = token; }
    ~UnwindToken() { active_unwind_token = previous; }
};
struct Protect {
    int count;
    Protect() : count(0) {}
    SEXP add(SEXP x) { PROTECT(x); ++count; return x; }
    ~Protect() { UNPROTECT(count); }
};

void fail(const char *message) { throw std::runtime_error(message); }
void check_interrupt() {
    r_api([]() { R_CheckUserInterrupt(); return R_NilValue; });
}

SEXP field(SEXP object, const char *name) {
    if (TYPEOF(object) != VECSXP) fail("Native grid specification must be a named list.");
    SEXP names = R_API(Rf_getAttrib(object, R_NamesSymbol));
    if (TYPEOF(names) != STRSXP) fail("Native grid specification must be named.");
    for (R_xlen_t i = 0; i < XLENGTH(object); ++i)
        if (!std::strcmp(CHAR(STRING_ELT(names, i)), name)) return VECTOR_ELT(object, i);
    throw std::runtime_error(std::string("Missing native grid field: ") + name);
}

void finite_vector(SEXP x, const char *name) {
    if (TYPEOF(x) != REALSXP) throw std::runtime_error(std::string(name) + " must be numeric.");
    const double *values = R_API(REAL(x));
    for (R_xlen_t i = 0; i < XLENGTH(x); ++i)
        if (!R_FINITE(values[i])) throw std::runtime_error(std::string(name) + " must be finite.");
}

void matrix_shape(SEXP x, int rows, int columns, const char *name) {
    finite_vector(x, name);
    SEXP dimensions = R_API(Rf_getAttrib(x, R_DimSymbol));
    if (TYPEOF(dimensions) != INTSXP || XLENGTH(dimensions) != 2 ||
        INTEGER(dimensions)[0] != rows || INTEGER(dimensions)[1] != columns)
        throw std::runtime_error(std::string(name) + " has incompatible dimensions.");
}

double scalar(SEXP x, const char *name) {
    finite_vector(x, name);
    if (XLENGTH(x) != 1) throw std::runtime_error(std::string(name) + " must be a scalar.");
    return R_API(REAL(x))[0];
}

// Build the same row-scaled jump kernels as make_markov_kernel(), followed by
// the tiny row normalization in prepare_markov_recursion(). The factored
// representation avoids allocating/validating a separate R matrix at every
// jump of every baseline profile. All supplied values remain checked here.
void factored_kernels(SEXP factors, int jumps, int k,
                      std::vector<double>& storage, std::vector<const double *>& kernels,
                      double& maximum_row_error, size_t& adjusted_rows) {
    SEXP base_r = field(factors, "base_increments");
    SEXP risk_r = field(factors, "relative_risk");
    SEXP from_r = field(factors, "from"), to_r = field(factors, "to");
    SEXP absorbing_r = field(factors, "absorbing");
    finite_vector(risk_r, "relative_risk");
    finite_vector(from_r, "from"); finite_vector(to_r, "to");
    finite_vector(absorbing_r, "absorbing");
    if (XLENGTH(risk_r) > INT_MAX ||
        R_API(Rf_getAttrib(risk_r, R_DimSymbol)) != R_NilValue ||
        R_API(Rf_getAttrib(from_r, R_DimSymbol)) != R_NilValue ||
        R_API(Rf_getAttrib(to_r, R_DimSymbol)) != R_NilValue ||
        R_API(Rf_getAttrib(absorbing_r, R_DimSymbol)) != R_NilValue)
        fail("Factored kernel risks and state indices must be numeric vectors.");
    const int edges = static_cast<int>(XLENGTH(risk_r));
    if (XLENGTH(from_r) != edges || XLENGTH(to_r) != edges)
        fail("Factored kernel edges and relative risks have incompatible lengths.");
    matrix_shape(base_r, jumps, edges, "base_increments");
    const double *base = REAL(base_r), *risk = REAL(risk_r);
    const double *from = REAL(from_r), *to = REAL(to_r), *absorbing = REAL(absorbing_r);
    std::vector<unsigned char> terminal(k);
    for (R_xlen_t a = 0; a < XLENGTH(absorbing_r); ++a) {
        const double state = absorbing[a];
        if (state < 1 || state > k || state != std::floor(state))
            fail("Factored kernel absorbing states must index fitted states.");
        const int s = static_cast<int>(state) - 1;
        if (terminal[s]) fail("Factored kernel absorbing states must be distinct.");
        terminal[s] = 1;
    }
    if (static_cast<size_t>(k) > std::numeric_limits<size_t>::max() / static_cast<size_t>(k))
        fail("Factored kernel dimensions overflow the memory address range.");
    const size_t layer = static_cast<size_t>(k) * k;
    if (jumps && layer > std::numeric_limits<size_t>::max() / sizeof(double) / static_cast<size_t>(jumps))
        fail("Factored kernel dimensions overflow the memory address range.");
    std::vector<size_t> positions(edges);
    std::unordered_map<size_t, unsigned char> unique_edges;
    unique_edges.reserve(edges);
    for (int e = 0; e < edges; ++e) {
        if (from[e] < 1 || from[e] > k || from[e] != std::floor(from[e]) ||
            to[e] < 1 || to[e] > k || to[e] != std::floor(to[e]) || from[e] == to[e])
            fail("Factored kernel edges must join distinct fitted states.");
        if (!(risk[e] > 0)) fail("Factored kernel relative risks must be positive.");
        positions[e] = static_cast<size_t>(from[e] - 1) + static_cast<size_t>(k) * static_cast<size_t>(to[e] - 1);
        if (!unique_edges.emplace(positions[e], 1).second)
            fail("Factored kernel edges must be distinct.");
    }
    storage.assign(layer * static_cast<size_t>(jumps), 0);
    std::vector<double> row_sums(k);
    for (int j = 0; j < jumps; ++j) {
        if ((j & 255) == 0) check_interrupt();
        double *kernel = storage.data() + layer * j;
        for (int e = 0; e < edges; ++e) {
            const double increment = base[j + static_cast<R_xlen_t>(jumps) * e];
            if (increment < 0) fail("Factored kernel baseline increments must be nonnegative.");
            const double value = increment * risk[e];
            if (!R_FINITE(value)) fail("Factored Markov intensity increments overflow.");
            if (increment > 0 && value == 0)
                fail("A positive factored Markov intensity increment underflowed to zero.");
            if (terminal[static_cast<int>(from[e]) - 1] && value > 0)
                fail("Declared absorbing states have outgoing intensity increments.");
            kernel[positions[e]] = value;
        }
        for (int s = 0; s < k; ++s) {
            long double total = 0;
            for (int t = 0; t < k; ++t) total += kernel[s + static_cast<size_t>(k)*t];
            const double outgoing = static_cast<double>(total);
            if (!R_FINITE(outgoing)) fail("Factored Markov intensity row sums overflow.");
            if (outgoing > 1) {
                const double row_scale = 1 / outgoing;
                for (int t = 0; t < k; ++t) {
                    double& probability = kernel[s + static_cast<size_t>(k)*t];
                    const double positive = probability;
                    probability *= row_scale;
                    if (positive > 0 && probability == 0)
                        fail("A positive factored Markov kernel probability underflowed to zero.");
                }
                ++adjusted_rows;
            }
            total = 0;
            for (int t = 0; t < k; ++t) total += kernel[s + static_cast<size_t>(k)*t];
            kernel[s + static_cast<size_t>(k)*s] = std::max(0.0, 1 - static_cast<double>(total));
        }
        double layer_error = 0;
        for (int s = 0; s < k; ++s) {
            long double total = 0;
            for (int t = 0; t < k; ++t) total += kernel[s + static_cast<size_t>(k)*t];
            row_sums[s] = static_cast<double>(total);
            layer_error = std::max(layer_error, std::fabs(row_sums[s] - 1));
        }
        if (layer_error > 1e-12) fail("Factored Markov kernel rows lost unit mass.");
        maximum_row_error = std::max(maximum_row_error, layer_error);
        if (layer_error > 0)
            for (int t = 0; t < k; ++t) for (int s = 0; s < k; ++s) {
                double& probability = kernel[s + static_cast<size_t>(k)*t];
                const double positive = probability;
                probability /= row_sums[s];
                if (positive > 0 && probability == 0)
                    fail("A positive factored Markov kernel probability underflowed to zero.");
            }
        kernels[j] = kernel;
    }
}

struct Distribution {
    size_t n;
    int d;
    std::vector<double> mass, value, lower, upper, variance;
    std::vector<unsigned char> protected_atom, variance_reliable, constant;
    Distribution(size_t count, int dimensions) : n(count), d(dimensions),
        mass(count), value(count * dimensions), lower(count * dimensions),
        upper(count * dimensions), variance(count * dimensions), protected_atom(count),
        variance_reliable(count * dimensions, 1), constant(count * dimensions, 1) {}
};

struct CoordinateHash {
    size_t operator()(const std::vector<double>& x) const {
        size_t h = 0;
        for (size_t j = 0; j < x.size(); ++j)
            h ^= std::hash<double>()(x[j]) + static_cast<size_t>(0x9e3779b9U) + (h << 6) + (h >> 2);
        return h;
    }
};

// Work arrays belong to one native call. Reuse capacity across state layers
// and queries, while keeping recursive R callbacks and later calls isolated.
struct GridWorkspace {
    std::vector<double> cells, minima, maxima, strides, coordinate;
    std::vector<size_t> group, dense_groups;
    std::vector<std::uint64_t> flat_keys;
    std::vector<size_t> flat_groups;
    std::vector<std::uint32_t> flat_stamps;
    std::uint32_t flat_epoch;
    std::unordered_map<std::vector<double>, size_t, CoordinateHash> coordinate_groups;
    std::unordered_map<std::vector<double>, size_t, CoordinateHash> point_groups;
    std::vector<unsigned char> heavy_point;
    std::vector<double> standardized, minimum, maximum, deviation;
    std::vector<unsigned char> ambiguous;
    size_t dense_layers, hashed_layers, maximum_point_groups;
    double largest_unprotected_point_mass;
    const double point_threshold;
    // Per-call partition control. Intrinsic protected atoms are unaffected;
    // a threshold above one disables only additional exact-point retention.
    explicit GridWorkspace(double threshold) : flat_epoch(0), dense_layers(0), hashed_layers(0), maximum_point_groups(0), largest_unprotected_point_mass(0), point_threshold(threshold) {}
    void prepare_flat(size_t n) {
        size_t capacity = 1;
        while (capacity < 2 * n) capacity *= 2;
        if (flat_keys.size() < capacity) {
            flat_keys.resize(capacity); flat_groups.resize(capacity);
            flat_stamps.assign(capacity, 0); flat_epoch = 1;
        } else if (++flat_epoch == 0) {
            std::fill(flat_stamps.begin(), flat_stamps.end(), 0); flat_epoch = 1;
        }
    }
    size_t flat_group(double numeric_key, size_t& groups) {
        const std::uint64_t key = static_cast<std::uint64_t>(numeric_key);
        // Integer cell keys are exact below 2^52. Mixing before masking avoids
        // collisions from regular grid strides while preserving insertion order.
        std::uint64_t hash = key;
        hash ^= hash >> 33; hash *= UINT64_C(0xff51afd7ed558ccd);
        hash ^= hash >> 33; hash *= UINT64_C(0xc4ceb9fe1a85ec53);
        hash ^= hash >> 33;
        const size_t mask = flat_keys.size() - 1;
        size_t position = static_cast<size_t>(hash) & mask;
        while (flat_stamps[position] == flat_epoch && flat_keys[position] != key)
            position = (position + 1) & mask;
        if (flat_stamps[position] != flat_epoch) {
            flat_stamps[position] = flat_epoch; flat_keys[position] = key;
            flat_groups[position] = groups++;
        }
        return flat_groups[position];
    }
};

std::shared_ptr<Distribution> coalesce(const std::shared_ptr<Distribution>& input,
                                      const double *step, const std::vector<double>& origin,
                                      GridWorkspace& work) {
    const size_t n = input->n;
    const int d = input->d;
    std::vector<double>& cells = work.cells;
    std::vector<double>& minima = work.minima;
    std::vector<double>& maxima = work.maxima;
    std::vector<double>& strides = work.strides;
    cells.resize(n * d); minima.assign(d, R_PosInf);
    maxima.assign(d, R_NegInf); strides.resize(d);
    work.heavy_point.resize(n);
    size_t heavy_points = 0;
    for (size_t a = 0; a < n; ++a) {
        bool exact_point = !input->protected_atom[a];
        for (int j = 0; j < d; ++j) {
            double cell = std::isinf(step[j]) ? 0 : std::floor((input->value[a*d+j] - origin[j]) / step[j]);
            if (!R_FINITE(cell)) fail("Native Markov grid cell coordinates overflow.");
            cells[a*d+j] = cell;
            minima[j] = std::min(minima[j], cell);
            maxima[j] = std::max(maxima[j], cell);
            exact_point = exact_point && input->lower[a*d+j] == input->upper[a*d+j];
        }
        if (exact_point) work.largest_unprotected_point_mass =
            std::max(work.largest_unprotected_point_mass, input->mass[a]);
        // Preserve important discrete atoms exactly instead of mixing them
        // with a neighboring continuous cell. This status is recomputed from
        // the physical point envelope and current mass at every coalescence;
        // low-mass descendants are not permanently protected.
        work.heavy_point[a] = exact_point && input->mass[a] >= work.point_threshold;
        heavy_points += work.heavy_point[a];
    }
    double product = 1;
    for (int j = 0; j < d; ++j) {
        strides[j] = product;
        // Different component proofs cannot share a cell: a constant
        // envelope endpoint requires the same canonical query arithmetic.
        product *= 2 * (maxima[j] - minima[j] + 1);
    }
    const bool numeric = R_FINITE(product) && product < 4503599627370496.0;
    // A compact mixed-radix domain needs one lookup per node instead of one
    // separately allocated hash entry. Bound scratch space by both occupied
    // node count and a fixed one-million-entry ceiling; sparse/high-dimensional
    // domains keep the existing hash representation.
    const bool dense = numeric && product <= 4.0 * static_cast<double>(n) &&
        product <= 1048576.0;
    std::unordered_map<std::vector<double>, size_t, CoordinateHash>& coordinate_groups = work.coordinate_groups;
    if (dense) {
        work.dense_groups.assign(static_cast<size_t>(product), n);
        ++work.dense_layers;
    } else if (numeric) {
        work.prepare_flat(n); ++work.hashed_layers;
    } else {
        coordinate_groups.clear();
        if (coordinate_groups.bucket_count() < n) coordinate_groups.reserve(n);
        ++work.hashed_layers;
    }
    std::vector<size_t>& group = work.group;
    std::vector<double>& coordinate = work.coordinate;
    group.resize(n); coordinate.resize(2*static_cast<size_t>(d));
    auto& point_groups = work.point_groups;
    point_groups.clear();
    if (point_groups.bucket_count() < heavy_points) point_groups.reserve(heavy_points);
    size_t groups = 0, protected_group = n;
    for (size_t a = 0; a < n; ++a) {
        if (input->protected_atom[a]) {
            if (protected_group == n) protected_group = groups++;
            group[a] = protected_group;
        } else if (work.heavy_point[a]) {
            for (int j = 0; j < d; ++j) {
                coordinate[j] = input->value[a*d+j];
                coordinate[static_cast<size_t>(d)+j] = input->constant[a*d+j];
            }
            const auto inserted = point_groups.emplace(coordinate, groups);
            if (inserted.second) ++groups;
            group[a] = inserted.first->second;
        } else if (numeric) {
            double key = 0;
            for (int j = 0; j < d; ++j)
                key += (2 * (cells[a*d+j] - minima[j]) + input->constant[a*d+j]) * strides[j];
            if (dense) {
                size_t& destination = work.dense_groups[static_cast<size_t>(key)];
                if (destination == n) destination = groups++;
                group[a] = destination;
            } else {
                group[a] = work.flat_group(key, groups);
            }
        } else {
            for (int j = 0; j < d; ++j) {
                coordinate[j] = cells[a*d+j];
                coordinate[static_cast<size_t>(d)+j] = input->constant[a*d+j];
            }
            std::pair<std::unordered_map<std::vector<double>, size_t, CoordinateHash>::iterator, bool> inserted =
                coordinate_groups.emplace(coordinate, groups);
            if (inserted.second) ++groups;
            group[a] = inserted.first->second;
        }
    }
    work.maximum_point_groups = std::max(work.maximum_point_groups, point_groups.size());
    if (groups == n) return input;
    std::shared_ptr<Distribution> output(new Distribution(groups, d));
    std::fill(output->lower.begin(), output->lower.end(), R_PosInf);
    std::fill(output->upper.begin(), output->upper.end(), R_NegInf);
    for (size_t a = 0; a < n; ++a) {
        const size_t g = group[a];
        const double p = input->mass[a];
        output->mass[g] += p;
        output->protected_atom[g] |= input->protected_atom[a];
        for (int j = 0; j < d; ++j) {
            output->constant[g*d+j] &= input->constant[a*d+j];
            output->value[g*d+j] += input->value[a*d+j] * p;
            output->lower[g*d+j] = std::min(output->lower[g*d+j], input->lower[a*d+j]);
            output->upper[g*d+j] = std::max(output->upper[g*d+j], input->upper[a*d+j]);
        }
    }
    for (size_t g = 0; g < groups; ++g) {
        if (!(output->mass[g] > 0) || !R_FINITE(output->mass[g])) fail("Invalid coalesced Markov grid mass.");
        for (int j = 0; j < d; ++j) {
            double value = output->value[g*d+j] / output->mass[g];
            output->value[g*d+j] = std::max(output->lower[g*d+j], std::min(output->upper[g*d+j], value));
        }
    }
    for (size_t a = 0; a < n; ++a) {
        const size_t g = group[a];
        for (int j = 0; j < d; ++j) {
            const double delta = input->value[a*d+j] - output->value[g*d+j];
            const double squared = delta*delta;
            const double summand = input->variance[a*d+j] + squared;
            const double weighted = summand * input->mass[a];
            if (!input->variance_reliable[a*d+j] ||
                (delta != 0 && squared < std::numeric_limits<double>::min()) ||
                (summand > 0 && weighted < std::numeric_limits<double>::min()))
                output->variance_reliable[g*d+j] = 0;
            output->variance[g*d+j] += weighted;
        }
    }
    for (size_t g = 0; g < groups; ++g)
        for (int j = 0; j < d; ++j) {
            output->variance[g*d+j] /= output->mass[g];
            if (!R_FINITE(output->variance[g*d+j]) || output->variance[g*d+j] < 0)
                fail("Within-node reward variances overflow; check reward units.");
            if (output->variance[g*d+j] > 0 &&
                output->variance[g*d+j] < std::numeric_limits<double>::min())
                output->variance_reliable[g*d+j] = 0;
            if (output->constant[g*d+j]) {
                output->variance[g*d+j] = 0;
                output->variance_reliable[g*d+j] = 1;
            }
        }
    return output;
}

struct Coverage {
    std::vector<double> probability, lower, upper;
    Coverage(size_t length) : probability(length), lower(length), upper(length) {}
};

void compensated_add(double value, double& sum, double& correction) {
    const double updated = sum + value;
    if (std::fabs(sum) >= std::fabs(value)) correction += (sum - updated) + value;
    else correction += (value - updated) + sum;
    sum = updated;
}

void prefix_probability(std::vector<double>& x) {
    long double sum = 0;
    for (size_t a = 0; a < x.size(); ++a) {
        sum += x[a];
        x[a] = std::max(0.0, std::min(1.0, static_cast<double>(sum)));
    }
}

Coverage coverage(const Distribution& nodes, const std::vector<double>& accrued,
                  const std::vector<double>& center, const std::vector<double>& scale,
                  const std::vector<double>& shift, const std::vector<double>& protected_total,
                  const double *q, size_t length, GridWorkspace& work) {
    const size_t n = nodes.n;
    const int d = nodes.d;
    Coverage output(length);
    std::vector<double>& standardized = work.standardized;
    std::vector<double>& minimum = work.minimum;
    std::vector<double>& maximum = work.maximum;
    standardized.resize(n*d); minimum.resize(n); maximum.resize(n);
    for (size_t a = 0; a < n; ++a) {
        double score = 0, low_score = 0, high_score = 0;
        for (int j = 0; j < d; ++j) {
            // Form the complete target before evaluating its frozen score.
            // Reassociating accrued-center can change inclusive atom ties.
            // Query times translate a view of the outcome-jump-anchored
            // distribution. They must not split a sojourn in the stored law.
            // Keep (future + shift) + accrued, including at exact atoms.
            const double total = nodes.constant[a*d+j] ? protected_total[j] : (nodes.value[a*d+j] + shift[j]) + accrued[j];
            const double total_low = nodes.constant[a*d+j] ? protected_total[j] : (nodes.lower[a*d+j] + shift[j]) + accrued[j];
            const double total_high = nodes.constant[a*d+j] ? protected_total[j] : (nodes.upper[a*d+j] + shift[j]) + accrued[j];
            const double value = (total - center[j]) / scale[j];
            const double low = (total_low - center[j]) / scale[j];
            const double high = (total_high - center[j]) / scale[j];
            if (!R_FINITE(value) || !R_FINITE(low) || !R_FINITE(high))
                fail("Standardized reward scores overflow; check centers and scales.");
            standardized[a*d+j] = value;
            score = std::max(score, std::fabs(value));
            low_score = std::max(low_score, std::max(low, -high));
            high_score = std::max(high_score, std::max(std::fabs(low), std::fabs(high)));
        }
        minimum[a] = low_score;
        maximum[a] = high_score;
        size_t index = std::lower_bound(q, q + length, score) - q;
        if (index < length) output.probability[index] += nodes.mass[a];
        index = std::lower_bound(q, q + length, high_score) - q;
        if (index < length) output.lower[index] += nodes.mass[a];
        index = std::lower_bound(q, q + length, low_score) - q;
        if (index < length) output.upper[index] += nodes.mass[a];
    }
    prefix_probability(output.probability);
    prefix_probability(output.lower);
    prefix_probability(output.upper);
    std::vector<unsigned char>& ambiguous = work.ambiguous;
    ambiguous.resize(length);
    bool any_ambiguous = false;
    for (size_t a = 0; a < length; ++a) {
        ambiguous[a] = output.upper[a] - output.lower[a] > 1e-14;
        any_ambiguous = any_ambiguous || ambiguous[a];
    }
    if (any_ambiguous) {
        std::vector<double>& deviation = work.deviation;
        deviation.resize(d);
        for (size_t a = 0; a < n; ++a) {
            const size_t first = std::lower_bound(q, q + length, minimum[a]) - q;
            const size_t end = std::lower_bound(q, q + length, maximum[a]) - q;
            // Variances are checked when a node is created/coalesced. A node
            // whose complete envelope settles every requested threshold needs
            // no square roots or moment-bound calculations at this query.
            if (first == end) continue;
            for (int j = 0; j < d; ++j) {
                const double variance = nodes.variance[a*d+j];
                deviation[j] = nodes.variance_reliable[a*d+j] ?
                    std::sqrt(variance) / scale[j] : R_PosInf;
                if (variance > 0 && deviation[j] < std::numeric_limits<double>::min())
                    deviation[j] = R_PosInf;
            }
            for (size_t b = first; b < end; ++b) {
                if (!ambiguous[b]) continue;
                double failure = 0, inside_upper = 1;
                for (int j = 0; j < d; ++j) {
                    const double mean = standardized[a*d+j];
                    const double distances[2] = {mean + q[b], q[b] - mean};
                    for (int face = 0; face < 2; ++face) {
                        const double distance = distances[face];
                        double ratio = 0;
                        if (deviation[j] > 0) {
                            const double z = distance / deviation[j];
                            ratio = 1 / (1 + z*z);
                        }
                        failure += distance >= 0 ? ratio : 1;
                        inside_upper = std::min(inside_upper, distance >= 0 ? 1 : ratio);
                    }
                }
                output.lower[b] += nodes.mass[a] * std::max(0.0, 1 - failure);
                output.upper[b] -= nodes.mass[a] * (1 - inside_upper);
            }
            if ((a & 4095U) == 0) check_interrupt();
        }
    }
    for (size_t a = 0; a < length; ++a)
        output.probability[a] = std::max(output.lower[a], std::min(output.upper[a], output.probability[a]));
    return output;
}

SEXP named_list(Protect& protection, const std::vector<const char *>& names) {
    SEXP result = protection.add(R_API(Rf_allocVector(VECSXP, names.size())));
    SEXP labels = protection.add(R_API(Rf_allocVector(STRSXP, names.size())));
    for (size_t i = 0; i < names.size(); ++i) SET_STRING_ELT(labels, i, R_API(Rf_mkChar(names[i])));
    R_API(Rf_setAttrib(result, R_NamesSymbol, labels));
    return result;
}

SEXP numeric_result(Protect& protection, const std::vector<double>& x) {
    SEXP result = protection.add(R_API(Rf_allocVector(REALSXP, x.size())));
    std::copy(x.begin(), x.end(), REAL(result));
    return result;
}

SEXP evaluated(SEXP call) {
    return R_API(Rf_eval(call, R_BaseEnv));
}

void transformed_coverage(SEXP transform, int index, Coverage& value) {
    if (transform == R_NilValue) return;
    Protect protection;
    SEXP values = named_list(protection, {"probability", "lower", "upper"});
    SET_VECTOR_ELT(values, 0, numeric_result(protection, value.probability));
    SET_VECTOR_ELT(values, 1, numeric_result(protection, value.lower));
    SET_VECTOR_ELT(values, 2, numeric_result(protection, value.upper));
    SEXP row = protection.add(R_API(Rf_ScalarInteger(index + 1)));
    SEXP call = protection.add(R_API(Rf_lang3(transform, row, values)));
    SEXP result = protection.add(evaluated(call));
    const char *names[3] = {"probability", "lower", "upper"};
    std::vector<double> *destinations[3] = {&value.probability, &value.lower, &value.upper};
    for (int j = 0; j < 3; ++j) {
        SEXP x = R_NilValue;
        try { x = field(result, names[j]); }
        catch (const std::exception&) { fail("query_transform returned invalid probability brackets."); }
        if ((TYPEOF(x) != REALSXP && TYPEOF(x) != INTSXP) ||
            XLENGTH(x) != static_cast<R_xlen_t>(value.probability.size()) ||
            R_API(Rf_getAttrib(x, R_DimSymbol)) != R_NilValue)
            fail("query_transform returned invalid probability brackets.");
        const double *real = TYPEOF(x) == REALSXP ? R_API(REAL(x)) : NULL;
        const int *integer = TYPEOF(x) == INTSXP ? R_API(INTEGER(x)) : NULL;
        for (size_t a = 0; a < value.probability.size(); ++a) {
            double z = real ? real[a] : integer[a];
            if (!R_FINITE(z) || (integer && integer[a] == NA_INTEGER))
                fail("query_transform returned invalid probability brackets.");
            (*destinations[j])[a] = z;
        }
    }
    for (size_t a = 0; a < value.probability.size(); ++a)
        if (value.lower[a] < -1e-12 || value.upper[a] > 1 + 1e-12 ||
            value.probability[a] < value.lower[a] - 1e-12 ||
            value.probability[a] > value.upper[a] + 1e-12)
            fail("query_transform returned invalid probability brackets.");
}

struct Event {
    double time;
    int jump;
    std::vector<int> right, left;
    Event(double t) : time(t), jump(-1) {}
};

SEXP run_grid(SEXP specification) {
    Protect protection;
    SEXP times_r = field(specification, "jump_times");
    SEXP kernels_r = field(specification, "kernels");
    SEXP rewards_r = field(specification, "rewards");
    SEXP query_times_r = field(specification, "query_times");
    SEXP states_r = field(specification, "query_state");
    SEXP accrued_r = field(specification, "accrued");
    SEXP center_r = field(specification, "center");
    SEXP scale_r = field(specification, "scale");
    SEXP q_r = field(specification, "q");
    SEXP side_r = field(specification, "side");
    SEXP step_r = field(specification, "step");
    SEXP origin_r = field(specification, "origin");
    SEXP transform = field(specification, "query_transform");
    SEXP reducer = field(specification, "reducer");
    SEXP known_r = field(specification, "query_known");
    SEXP weights_r = field(specification, "query_weights");
    const double tau = scalar(field(specification, "tau"), "tau");
    const double budget = scalar(field(specification, "max_nodes"), "max_nodes");
    const double exact = scalar(field(specification, "exact_limit"), "exact_limit");
    const double point_threshold = scalar(field(specification, "point_mass_threshold"), "point_mass_threshold");
    if (point_threshold < 0.001 || point_threshold > 2)
        fail("point_mass_threshold must lie between 0.001 and 2.");
    finite_vector(times_r, "jump_times"); finite_vector(query_times_r, "query_times");
    finite_vector(states_r, "query_state"); finite_vector(q_r, "q");
    finite_vector(rewards_r, "rewards"); finite_vector(origin_r, "origin");
    SEXP reward_dimensions = R_API(Rf_getAttrib(rewards_r, R_DimSymbol));
    if (TYPEOF(reward_dimensions) != INTSXP || XLENGTH(reward_dimensions) != 2)
        fail("rewards must be a state-by-target matrix.");
    const int k = INTEGER(reward_dimensions)[0], d = INTEGER(reward_dimensions)[1];
    if (k < 1 || d < 1 || !(tau > 0)) fail("Invalid reward dimensions or prediction horizon.");
    if (XLENGTH(query_times_r) < 1 || XLENGTH(query_times_r) > INT_MAX || XLENGTH(q_r) < 1 || XLENGTH(q_r) > INT_MAX ||
        XLENGTH(times_r) > INT_MAX) fail("Native grid query dimensions exceed supported limits.");
    const int h = static_cast<int>(XLENGTH(query_times_r));
    const int length = static_cast<int>(XLENGTH(q_r));
    const int jumps = static_cast<int>(XLENGTH(times_r));
    const bool factored = kernels_r == R_NilValue;
    if (!factored && (TYPEOF(kernels_r) != VECSXP || XLENGTH(kernels_r) != jumps))
        fail("kernels must supply one matrix per jump.");
    matrix_shape(accrued_r, h, d, "accrued"); matrix_shape(center_r, h, d, "center"); matrix_shape(scale_r, h, d, "scale");
    if (XLENGTH(states_r) != h || TYPEOF(side_r) != INTSXP || XLENGTH(side_r) != h ||
        TYPEOF(step_r) != REALSXP || XLENGTH(step_r) != d || XLENGTH(origin_r) != d)
        fail("Native grid query vectors have incompatible lengths.");
    if (budget < k || budget != std::floor(budget) || budget > INT_MAX ||
        exact < 1 || exact != std::floor(exact) || exact > INT_MAX)
        fail("Native Markov node budgets must be valid integers.");
    const size_t max_nodes = static_cast<size_t>(budget), exact_limit = static_cast<size_t>(exact);
    if (max_nodes > std::numeric_limits<size_t>::max() / static_cast<size_t>(d) / sizeof(double) / 4)
        fail("Native Markov node dimensions overflow the memory address range.");
    const double *times = R_API(REAL(times_r)), *rewards = R_API(REAL(rewards_r)), *query_times = R_API(REAL(query_times_r));
    const double *states = R_API(REAL(states_r)), *q = R_API(REAL(q_r)), *step = R_API(REAL(step_r)), *origin = R_API(REAL(origin_r));
    const double *accrued = R_API(REAL(accrued_r)), *centers = R_API(REAL(center_r)), *scales = R_API(REAL(scale_r));
    const int *side = R_API(INTEGER(side_r));
    const double *known = NULL, *weights = NULL;
    if (known_r != R_NilValue) {
        SEXP dimensions = R_API(Rf_getAttrib(known_r, R_DimSymbol));
        if (TYPEOF(known_r) != REALSXP || TYPEOF(dimensions) != INTSXP || XLENGTH(dimensions) != 2 ||
            INTEGER(dimensions)[0] != h || INTEGER(dimensions)[1] != length)
            fail("query_known must be a query-by-threshold matrix.");
        known = R_API(REAL(known_r));
        for (R_xlen_t a = 0; a < XLENGTH(known_r); ++a)
            if (!ISNAN(known[a]) && known[a] != 0 && known[a] != 1)
                fail("query_known entries must be NA, zero, or one.");
    }
    if (weights_r != R_NilValue) {
        finite_vector(weights_r, "query_weights");
        if (XLENGTH(weights_r) != h || reducer != R_NilValue)
            fail("query_weights must have one value per query and cannot accompany a reducer.");
        weights = R_API(REAL(weights_r));
    }
    for (int j = 0; j < d; ++j) if (!(step[j] > 0) || ISNAN(step[j])) fail("Grid cell widths must be positive, allowing infinity for constant rewards.");
    for (int a = 0; a < length; ++a) if (q[a] < 0 || (a && q[a] <= q[a-1])) fail("q must be nonnegative and strictly increasing.");
    for (int a = 0; a < h; ++a) {
        if (query_times[a] < 0 || query_times[a] > tau || states[a] < 1 || states[a] > k || states[a] != std::floor(states[a]) ||
            (side[a] != 0 && side[a] != 1)) fail("Invalid native Markov grid query.");
        for (int j = 0; j < d; ++j) if (!(scales[a+h*j] > 0)) fail("Query scales must be positive.");
    }
    std::vector<const double *> kernels(jumps);
    std::vector<double> kernel_storage;
    double native_kernel_row_error = 0;
    size_t adjusted_kernel_rows = 0;
    if (factored)
        factored_kernels(field(specification, "kernel_factors"), jumps, k,
                          kernel_storage, kernels, native_kernel_row_error, adjusted_kernel_rows);
    for (int j = 0; j < jumps; ++j) {
        if (!(times[j] > 0) || times[j] > tau || (j && times[j] <= times[j-1])) fail("Jump times must be positive, ordered, and at most tau.");
        if (factored) continue;
        SEXP kernel = VECTOR_ELT(kernels_r, j);
        matrix_shape(kernel, k, k, "kernel");
        kernels[j] = R_API(REAL(kernel));
        for (int s = 0; s < k; ++s) {
            long double sum = 0;
            for (int t = 0; t < k; ++t) {
                const double p = kernels[j][s+k*t];
                if (p < 0) fail("Kernel probabilities must be nonnegative.");
                sum += p;
            }
            if (std::fabs(static_cast<double>(sum) - 1) > 1e-12) fail("Kernel rows must sum to one.");
        }
    }
    if (transform != R_NilValue && !Rf_isFunction(transform)) fail("query_transform must be NULL or a function.");
    SEXP initialize = R_NilValue, add = R_NilValue, finalize = R_NilValue;
    SEXP state_box = protection.add(R_API(Rf_allocVector(VECSXP, 1)));
    if (reducer != R_NilValue) {
        initialize = field(reducer, "initialize"); add = field(reducer, "add"); finalize = field(reducer, "finalize");
        if (!Rf_isFunction(initialize) || !Rf_isFunction(add) || !Rf_isFunction(finalize)) fail("reducer must supply initialize, add, and finalize functions.");
        Protect temporary;
        SEXP call = temporary.add(R_API(Rf_lang1(initialize)));
        SEXP initialized = temporary.add(evaluated(call));
        SET_VECTOR_ELT(state_box, 0, initialized);
    }
    SEXP result = named_list(protection, {"probability", "lower", "upper", "reduced", "diagnostics"});
    double *probability_output = NULL, *lower_output = NULL, *upper_output = NULL;
    if (reducer == R_NilValue && weights == NULL) {
        SEXP p = protection.add(R_API(Rf_allocMatrix(REALSXP, h, length)));
        SEXP lo = protection.add(R_API(Rf_allocMatrix(REALSXP, h, length)));
        SEXP hi = protection.add(R_API(Rf_allocMatrix(REALSXP, h, length)));
        SET_VECTOR_ELT(result, 0, p); SET_VECTOR_ELT(result, 1, lo); SET_VECTOR_ELT(result, 2, hi);
        probability_output = REAL(p); lower_output = REAL(lo); upper_output = REAL(hi);
    }
    std::vector<double> event_times(query_times, query_times+h);
    const double minimum_query = *std::min_element(query_times, query_times+h);
    for (int j = 0; j < jumps; ++j) if (times[j] >= minimum_query) event_times.push_back(times[j]);
    std::sort(event_times.begin(), event_times.end(), std::greater<double>());
    event_times.erase(std::unique(event_times.begin(), event_times.end()), event_times.end());
    std::vector<Event> events;
    events.reserve(event_times.size());
    for (size_t a = 0; a < event_times.size(); ++a) {
        events.emplace_back(event_times[a]);
        const double *found = std::lower_bound(times, times+jumps, event_times[a]);
        if (found != times+jumps && *found == event_times[a]) events.back().jump = static_cast<int>(found-times);
    }
    for (int a = 0; a < h; ++a) {
        size_t index = std::lower_bound(event_times.begin(), event_times.end(), query_times[a], std::greater<double>()) - event_times.begin();
        if (side[a]) events[index].left.push_back(a); else events[index].right.push_back(a);
    }
    std::vector<std::shared_ptr<Distribution> > distributions;
    for (int s = 0; s < k; ++s) {
        std::shared_ptr<Distribution> initial(new Distribution(1, d));
        initial->mass[0] = 1; initial->protected_atom[0] = 1;
        distributions.push_back(initial);
    }
    std::vector<double> query_accrued(d), query_center(d), query_scale(d), query_shift(d), protected_total(d), shifted_origin(d);
    GridWorkspace workspace(point_threshold);
    std::vector<unsigned char> changed_state(k);
    double anchor = tau, maximum_width = 0, mass_error = 0, traversed = 0;
    int maximum_width_query = -1, maximum_width_threshold = -1;
    size_t peak_nodes = k;
    int coalesced_layers = 0;
    std::vector<double> weighted_probability(weights ? length : 0), weighted_lower(weights ? length : 0), weighted_upper(weights ? length : 0);
    std::vector<double> probability_correction(weights ? length : 0), lower_correction(weights ? length : 0), upper_correction(weights ? length : 0);
    std::vector<double> active_q;
    std::vector<int> active_indices;
    active_q.reserve(length); active_indices.reserve(length);
    const auto record = [&](int a) {
        const int state = static_cast<int>(states[a])-1;
        for (int j = 0; j < d; ++j) {
            query_accrued[j] = accrued[a+h*j]; query_center[j] = centers[a+h*j]; query_scale[j] = scales[a+h*j];
            query_shift[j] = (anchor-query_times[a]) * rewards[state+k*j];
            protected_total[j] = life_protected_total_reward(tau, query_times[a], rewards[state+k*j], query_accrued[j]);
        }
        Coverage value(length);
        if (known == NULL) {
            value = coverage(*distributions[state], query_accrued, query_center, query_scale, query_shift, protected_total, q, length, workspace);
        } else {
            active_q.clear(); active_indices.clear();
            for (int b = 0; b < length; ++b) {
                const double membership = known[a + static_cast<R_xlen_t>(h)*b];
                if (ISNAN(membership)) { active_q.push_back(q[b]); active_indices.push_back(b); }
                else value.probability[b] = value.lower[b] = value.upper[b] = membership;
            }
            if (!active_q.empty()) {
                Coverage unresolved = coverage(*distributions[state], query_accrued, query_center, query_scale, query_shift, protected_total,
                                               active_q.data(), active_q.size(), workspace);
                for (size_t b = 0; b < active_indices.size(); ++b) {
                    const int index = active_indices[b];
                    value.probability[index] = unresolved.probability[b];
                    value.lower[index] = unresolved.lower[b]; value.upper[index] = unresolved.upper[b];
                }
            }
        }
        transformed_coverage(transform, a, value);
        for (int b = 0; b < length; ++b) {
            const double width = value.upper[b] - value.lower[b];
            if (width > maximum_width) {
                maximum_width = width;
                maximum_width_query = a;
                maximum_width_threshold = b;
            }
        }
        if (weights != NULL) {
            const double w = weights[a];
            for (int b = 0; b < length; ++b) {
                compensated_add(w * value.probability[b], weighted_probability[b], probability_correction[b]);
                compensated_add(w * (w >= 0 ? value.lower[b] : value.upper[b]), weighted_lower[b], lower_correction[b]);
                compensated_add(w * (w >= 0 ? value.upper[b] : value.lower[b]), weighted_upper[b], upper_correction[b]);
            }
        } else if (reducer == R_NilValue) {
            for (int b = 0; b < length; ++b) {
                const R_xlen_t index = a + static_cast<R_xlen_t>(h)*b;
                probability_output[index] = value.probability[b]; lower_output[index] = value.lower[b]; upper_output[index] = value.upper[b];
            }
        } else {
            Protect temporary;
            SEXP row = temporary.add(R_API(Rf_ScalarInteger(a+1)));
            SEXP p = numeric_result(temporary, value.probability);
            SEXP lo = numeric_result(temporary, value.lower);
            SEXP hi = numeric_result(temporary, value.upper);
            SEXP call = temporary.add(R_API(Rf_lang6(add, VECTOR_ELT(state_box, 0), row, p, lo, hi)));
            SEXP updated = temporary.add(evaluated(call));
            SET_VECTOR_ELT(state_box, 0, updated);
        }
    };
    for (size_t event = 0; event < events.size(); ++event) {
        check_interrupt();
        const double time = events[event].time, duration = anchor - time;
        const int jump = events[event].jump;
        // Only fitted outcome jumps advance the stored distribution anchor.
        // Adding, removing or reordering history queries leaves that law fixed.
        if (jump >= 0 && duration > 0) {
            for (int s = 0; s < k; ++s) {
                Distribution& z = *distributions[s];
                for (int j = 0; j < d; ++j) {
                    if (rewards[s+k*j] == 0) continue;
                    const double shift = duration * rewards[s+k*j];
                    const double canonical = (tau-time) * rewards[s+k*j];
                    for (size_t a = 0; a < z.n; ++a) {
                        if (z.constant[a*d+j]) {
                            z.value[a*d+j] = z.lower[a*d+j] = z.upper[a*d+j] = canonical;
                            z.variance[a*d+j] = 0; z.variance_reliable[a*d+j] = 1;
                        } else if (shift != 0) {
                            z.value[a*d+j] += shift; z.lower[a*d+j] += shift; z.upper[a*d+j] += shift;
                        }
                    }
                }
            }
        }
        if (jump >= 0) anchor = time;
        for (size_t a = 0; a < events[event].right.size(); ++a) record(events[event].right[a]);
        if (jump >= 0) {
            const double *kernel = kernels[jump];
            std::vector<std::shared_ptr<Distribution> > next = distributions;
            std::fill(changed_state.begin(), changed_state.end(), 0);
            for (int s = 0; s < k; ++s) {
                size_t count = 0;
                int destinations = 0, sole = -1;
                for (int t = 0; t < k; ++t) if (kernel[s+k*t] > 0) {
                    count += distributions[t]->n; ++destinations; sole = t;
                    if (count > max_nodes) fail("Native Markov grid node budget exceeded.");
                }
                if (destinations == 1 && sole == s) continue;
                changed_state[s] = 1;
                std::shared_ptr<Distribution> mixed(new Distribution(count, d));
                size_t position = 0;
                for (int t = 0; t < k; ++t) {
                    const double weight = kernel[s+k*t];
                    if (!(weight > 0)) continue;
                    const Distribution& z = *distributions[t];
                    bool equivalent = true;
                    // A destination change at tau has zero occupation time.
                    // This flag proves a numerical atom, not endpoint-state
                    // equality or the private block departure category.
                    if (time != tau)
                        for (int j = 0; j < d; ++j) equivalent = equivalent && rewards[s+k*j] == rewards[t+k*j];
                    for (size_t a = 0; a < z.n; ++a, ++position) {
                        const double p = z.mass[a] * weight;
                        if (p == 0 && z.mass[a] > 0) fail("A positive Markov grid mass underflowed to zero.");
                        mixed->mass[position] = p;
                        mixed->protected_atom[position] = z.protected_atom[a] && equivalent;
                        for (int j = 0; j < d; ++j) {
                            mixed->value[position*d+j] = z.value[a*d+j];
                            mixed->lower[position*d+j] = z.lower[a*d+j];
                            mixed->upper[position*d+j] = z.upper[a*d+j];
                            mixed->variance[position*d+j] = z.variance[a*d+j];
                            mixed->variance_reliable[position*d+j] = z.variance_reliable[a*d+j];
                            mixed->constant[position*d+j] = z.constant[a*d+j] &&
                                (time == tau || rewards[s+k*j] == rewards[t+k*j]);
                        }
                    }
                }
                traversed += count;
                if (count > exact_limit) {
                    for (int j = 0; j < d; ++j) shifted_origin[j] = origin[j] - time*rewards[s+k*j];
                    mixed = coalesce(mixed, step, shifted_origin, workspace);
                    ++coalesced_layers;
                }
                next[s] = mixed;
            }
            distributions.swap(next);
            size_t total_nodes = 0;
            for (int s = 0; s < k; ++s) {
                total_nodes += distributions[s]->n;
                if (!changed_state[s]) continue;
                long double mass = 0;
                for (size_t a = 0; a < distributions[s]->n; ++a) mass += distributions[s]->mass[a];
                const double error = std::fabs(static_cast<double>(mass) - 1);
                if (!R_FINITE(error) || error > 1e-10) fail("Markov grid probabilities lost unit mass.");
                mass_error = std::max(mass_error, error);
            }
            peak_nodes = std::max(peak_nodes, total_nodes);
            if (total_nodes > max_nodes) fail("Native Markov grid node budget exceeded.");
        }
        for (size_t a = 0; a < events[event].left.size(); ++a) record(events[event].left[a]);
    }
    if (reducer != R_NilValue) {
        Protect temporary;
        SEXP call = temporary.add(R_API(Rf_lang2(finalize, VECTOR_ELT(state_box, 0))));
        SEXP reduced = temporary.add(evaluated(call));
        SET_VECTOR_ELT(result, 3, reduced);
    } else if (weights != NULL) {
        for (int b = 0; b < length; ++b) {
            weighted_probability[b] += probability_correction[b];
            weighted_lower[b] += lower_correction[b];
            weighted_upper[b] += upper_correction[b];
            if (!R_FINITE(weighted_probability[b]) || !R_FINITE(weighted_lower[b]) || !R_FINITE(weighted_upper[b]))
                fail("Signed Markov grid aggregation overflowed; check observation weights.");
        }
        SEXP reduced = named_list(protection, {"probability", "lower", "upper"});
        SET_VECTOR_ELT(reduced, 0, numeric_result(protection, weighted_probability));
        SET_VECTOR_ELT(reduced, 1, numeric_result(protection, weighted_lower));
        SET_VECTOR_ELT(reduced, 2, numeric_result(protection, weighted_upper));
        SET_VECTOR_ELT(result, 3, reduced);
    }
    SEXP diagnostics = named_list(protection, {"method", "peak_nodes", "coalesced_state_layers", "traversed_nodes",
        "mass_error", "maximum_probability_error", "query_count", "jump_count", "engine",
        "dense_grouping_layers", "hashed_grouping_layers", "native_kernel_row_error", "adjusted_kernel_rows",
        "maximum_error_query", "maximum_error_threshold", "maximum_error_time", "maximum_error_q",
        "largest_unprotected_point_mass", "point_mass_threshold", "maximum_retained_point_groups"});
    SET_VECTOR_ELT(diagnostics, 0, protection.add(R_API(Rf_mkString("deterministic_joint_grid_envelopes"))));
    SET_VECTOR_ELT(diagnostics, 1, protection.add(R_API(Rf_ScalarReal(peak_nodes))));
    SET_VECTOR_ELT(diagnostics, 2, protection.add(R_API(Rf_ScalarInteger(coalesced_layers))));
    SET_VECTOR_ELT(diagnostics, 3, protection.add(R_API(Rf_ScalarReal(traversed))));
    SET_VECTOR_ELT(diagnostics, 4, protection.add(R_API(Rf_ScalarReal(mass_error))));
    SET_VECTOR_ELT(diagnostics, 5, protection.add(R_API(Rf_ScalarReal(maximum_width))));
    SET_VECTOR_ELT(diagnostics, 6, protection.add(R_API(Rf_ScalarInteger(h))));
    SET_VECTOR_ELT(diagnostics, 7, protection.add(R_API(Rf_ScalarInteger(jumps))));
    SET_VECTOR_ELT(diagnostics, 8, protection.add(R_API(Rf_mkString("native_cpp"))));
    SET_VECTOR_ELT(diagnostics, 9, protection.add(R_API(Rf_ScalarReal(workspace.dense_layers))));
    SET_VECTOR_ELT(diagnostics, 10, protection.add(R_API(Rf_ScalarReal(workspace.hashed_layers))));
    SET_VECTOR_ELT(diagnostics, 11, protection.add(R_API(Rf_ScalarReal(native_kernel_row_error))));
    SET_VECTOR_ELT(diagnostics, 12, protection.add(R_API(Rf_ScalarReal(adjusted_kernel_rows))));
    SET_VECTOR_ELT(diagnostics, 13, protection.add(R_API(Rf_ScalarInteger(maximum_width_query + 1))));
    SET_VECTOR_ELT(diagnostics, 14, protection.add(R_API(Rf_ScalarInteger(maximum_width_threshold + 1))));
    SET_VECTOR_ELT(diagnostics, 15, protection.add(R_API(Rf_ScalarReal(maximum_width_query < 0 ? NA_REAL : query_times[maximum_width_query]))));
    SET_VECTOR_ELT(diagnostics, 16, protection.add(R_API(Rf_ScalarReal(maximum_width_threshold < 0 ? NA_REAL : q[maximum_width_threshold]))));
    SET_VECTOR_ELT(diagnostics, 17, protection.add(R_API(Rf_ScalarReal(workspace.largest_unprotected_point_mass))));
    SET_VECTOR_ELT(diagnostics, 18, protection.add(R_API(Rf_ScalarReal(point_threshold))));
    SET_VECTOR_ELT(diagnostics, 19, protection.add(R_API(Rf_ScalarReal(workspace.maximum_point_groups))));
    SET_VECTOR_ELT(result, 4, diagnostics);
    return result;
}
// Bounded scalar first-exit query recursion for marginal threshold initialization.
// Reuses this translation unit's factored kernels and unwind-safe R helpers.
struct MarginalQuery { int index, bin; };

int marginal_first_ge(const double *x, int n, double value) {
    int lo=0, hi=n;
    while(lo<hi) { int mid=lo+(hi-lo)/2;
        if(x[mid]<value) lo=mid+1; else hi=mid; }
    return lo;
}
int marginal_count_le(const double *x, int n, double value) {
    int lo=0, hi=n;
    while(lo<hi) { int mid=lo+(hi-lo)/2;
        if(x[mid]<=value) lo=mid+1; else hi=mid; }
    return lo;
}
int marginal_length(SEXP x, const char *name) {
    finite_vector(x,name);
    if(XLENGTH(x)>INT_MAX || R_API(Rf_getAttrib(x,R_DimSymbol))!=R_NilValue)
        fail("Scalar monotone vectors must have supported lengths and no dimensions.");
    return static_cast<int>(XLENGTH(x));
}
// Check between bounded row/edge/query blocks. Every block is at most L
// operations (L <= one million); small public profiles amortize the checks.
void marginal_work_check(size_t work, size_t& pending) {
    pending+=work;
    if(pending>=65536) { check_interrupt(); pending=0; }
}

SEXP run_monotone_queries(SEXP spec) {
    SEXP jt=field(spec,"jump_times"), rw=field(spec,"rewards");
    SEXP ti=field(spec,"query_times"), st=field(spec,"query_state");
    SEXP ac=field(spec,"accrued"), qs=field(spec,"q");
    SEXP factors=field(spec,"kernel_factors");
    const int jnum=marginal_length(jt,"jump_times");
    const int h=marginal_length(ti,"query_times");
    const int lnum=marginal_length(qs,"q");
    const int k=marginal_length(rw,"rewards");
    if(h<1 || lnum<1 || k<1 || marginal_length(st,"query_state")!=h ||
       marginal_length(ac,"accrued")!=h) fail("Invalid scalar monotone query dimensions.");
    const double tau=scalar(field(spec,"tau"),"tau");
    const double center=scalar(field(spec,"center"),"center");
    const double scale=scalar(field(spec,"scale"),"scale");
    const double cap=scalar(field(spec,"max_cells"),"max_cells");
    const double kernel_cap=scalar(field(spec,"max_kernel_cells"),"max_kernel_cells");
    if(tau<=0 || scale<=0 || cap<1 || cap>1000000 || cap!=std::floor(cap) ||
       kernel_cap<1 || kernel_cap>1000000 || kernel_cap!=std::floor(kernel_cap))
        fail("Invalid scalar monotone scale, horizon or memory cap.");
    const double lower_endpoint_score=std::fabs(center/scale);
    const double terminal_score=std::fabs((tau-center)/scale);
    if(!R_FINITE(lower_endpoint_score) || !R_FINITE(terminal_score))
        fail("Scalar monotone standardized scores overflow; use the ordinary grid.");
    if(h>cap/lnum || k>cap/lnum || k>kernel_cap/k ||
       (jnum>0 && static_cast<double>(k)*k>kernel_cap/jnum))
        fail("Scalar monotone query, state or kernel cell cap exceeded.");
    // Bound every factor index ledger before REAL/materialization or the
    // copied constructor's positions/hash allocations, including J=0.
    const R_xlen_t maximum_edges=static_cast<R_xlen_t>(k)*(k-1);
    if(XLENGTH(field(factors,"relative_risk"))>maximum_edges ||
       XLENGTH(field(factors,"from"))>maximum_edges ||
       XLENGTH(field(factors,"to"))>maximum_edges ||
       XLENGTH(field(factors,"absorbing"))>k)
        fail("Scalar monotone factor edge or absorbing count exceeds the state graph.");
    SEXP debug=field(spec,"inspect_kernels");
    if(TYPEOF(debug)!=LGLSXP || XLENGTH(debug)!=1 || LOGICAL(debug)[0]==NA_LOGICAL)
        fail("Invalid scalar kernel inspection flag.");
    const double *times=REAL(jt), *reward=REAL(rw), *query_time=REAL(ti);
    const double *state=REAL(st), *accrued=REAL(ac), *q=REAL(qs);
    for(int j=0;j<jnum;++j) if(times[j]<=0 || times[j]>tau || (j && times[j]<=times[j-1]))
        fail("Invalid scalar monotone jump times.");
    for(int a=0;a<lnum;++a) {
        if((a&4095)==0) check_interrupt();
        if(q[a]<0 || (a && q[a]<=q[a-1])) fail("Invalid scalar monotone thresholds.");
    }
    for(int s=0;s<k;++s) if(reward[s]!=0 && reward[s]!=1)
        fail("Scalar monotone rewards must be binary.");
    std::vector<MarginalQuery> queries(h);
    for(int a=0;a<h;++a) {
        if((a&4095)==0) check_interrupt();
        if(state[a]<1 || state[a]>k || state[a]!=std::floor(state[a]) ||
           query_time[a]<0 || query_time[a]>tau)
            fail("Invalid scalar monotone query.");
        if(reward[static_cast<int>(state[a])-1]==1 && accrued[a]!=query_time[a])
            fail("Active scalar monotone target must equal the query time.");
        if(reward[static_cast<int>(state[a])-1]==0 &&
           !R_FINITE((accrued[a]-center)/scale))
            fail("Scalar monotone standardized scores overflow; use the ordinary grid.");
        queries[a]={a,marginal_count_le(times,jnum,query_time[a])};
    }
    std::sort(queries.begin(),queries.end(),[](const MarginalQuery& a,const MarginalQuery& b) {
        return a.bin!=b.bin ? a.bin>b.bin : a.index<b.index;
    });
    double row_error=0; size_t adjusted=0;
    std::vector<double> storage;
    std::vector<const double *> kernels(jnum);
    // Use the same constructor and operation order as the general grid.
    factored_kernels(factors,jnum,k,storage,kernels,row_error,adjusted);
    SEXP fr=field(factors,"from"), tr=field(factors,"to");
    const int e=static_cast<int>(XLENGTH(fr));
    const double *from=REAL(fr), *to=REAL(tr);
    for(int a=0;a<e;++a) {
        if((a&4095)==0) check_interrupt();
        if(reward[static_cast<int>(to[a])-1]>reward[static_cast<int>(from[a])-1])
            fail("Scalar monotone edge increases a target reward.");
    }
    // Match the independent R clock oracle's destination-major edge order.
    std::vector<int> edge_order(e);
    for(int a=0;a<e;++a) edge_order[a]=a;
    std::sort(edge_order.begin(),edge_order.end(),[&](int a,int b) {
        return to[a]!=to[b] ? to[a]<to[b] : from[a]<from[b];
    });
    std::vector<double> buffer_a(static_cast<size_t>(k)*lnum),buffer_b(buffer_a.size(),1);
    std::vector<double *> value(k), next(k);
    const double full_score=std::max(lower_endpoint_score,terminal_score);
    const int active_lnum=marginal_first_ge(q,lnum,full_score);
    size_t work_since_interrupt=0;
    for(int s=0;s<k;++s) {
        marginal_work_check(lnum,work_since_interrupt);
        value[s]=buffer_a.data()+static_cast<size_t>(s)*lnum;
        next[s]=buffer_b.data()+static_cast<size_t>(s)*lnum;
        const int cutoff=reward[s] ? marginal_first_ge(q,lnum,terminal_score) : 0;
        for(int a=0;a<lnum;++a) value[s][a]=(a>=cutoff);
    }
    Protect protection;
    SEXP out=protection.add(R_API(Rf_allocMatrix(REALSXP,h,lnum)));
    double *answer=REAL(out);
    int cursor=0,bin=jnum;
    double state_updates=0;
    while(cursor<h) {
        while(cursor<h && queries[cursor].bin==bin) {
            marginal_work_check(lnum,work_since_interrupt);
            const int z=queries[cursor].index,s=static_cast<int>(state[z])-1;
            const double fixed=reward[s] ? 0 : std::fabs((accrued[z]-center)/scale);
            const int cutoff=marginal_first_ge(q,lnum,fixed);
            for(int a=0;a<lnum;++a) {
                double prob=a>=cutoff ? value[s][a] : 0;
                if(!R_FINITE(prob) || prob < -1e-10 || prob > 1+1e-10)
                    fail("Scalar monotone probability is outside [0,1].");
                if(prob<0) prob=0; else if(prob>1) prob=1;
                answer[z+static_cast<R_xlen_t>(h)*a]=prob;
            }
            ++cursor;
        }
        if(cursor>=h) break;
        const int j=--bin;
        if(j<0) fail("Scalar monotone query ordering failed.");
        if((j&255)==0) check_interrupt();
        const double *kernel=kernels[j];
        // Every row uses exactly its prepared diagonal, including rows with
        // no outgoing edge. Simultaneous updates retain rolling state buffers.
        for(int s=0;s<k;++s) {
            marginal_work_check(active_lnum,work_since_interrupt);
            const double stay=kernel[s+static_cast<size_t>(k)*s];
            for(int a=0;a<active_lnum;++a) next[s][a]=stay*value[s][a];
        }
        for(int index=0;index<e;++index) {
            marginal_work_check(1,work_since_interrupt);
            const int a=edge_order[index],f=static_cast<int>(from[a])-1,t=static_cast<int>(to[a])-1;
            const double probability=kernel[f+static_cast<size_t>(k)*t];
            if(probability==0) continue;
            const double threshold=reward[f]>reward[t] ? std::fabs((times[j]-center)/scale) : 0;
            const int cutoff=marginal_first_ge(q,lnum,threshold);
            marginal_work_check(std::max(0,active_lnum-cutoff),work_since_interrupt);
            for(int b=cutoff;b<active_lnum;++b) next[f][b]+=probability*value[t][b];
        }
        for(int s=0;s<k;++s) std::swap(value[s],next[s]);
        state_updates+=static_cast<double>(k)*active_lnum;
    }
    SEXP result=protection.add(R_API(Rf_allocVector(VECSXP,5)));
    SEXP names=protection.add(R_API(Rf_allocVector(STRSXP,5)));
    const char *labels[]={"probability","kernel_row_error","adjusted_rows","state_updates","kernel_storage"};
    for(int a=0;a<5;++a) SET_STRING_ELT(names,a,R_API(Rf_mkChar(labels[a])));
    SET_VECTOR_ELT(result,0,out);
    SET_VECTOR_ELT(result,1,protection.add(R_API(Rf_ScalarReal(row_error))));
    SET_VECTOR_ELT(result,2,protection.add(R_API(Rf_ScalarReal(adjusted))));
    SET_VECTOR_ELT(result,3,protection.add(R_API(Rf_ScalarReal(state_updates))));
    if(LOGICAL(debug)[0]) {
        SEXP inspect=protection.add(R_API(Rf_allocVector(REALSXP,storage.size())));
        std::copy(storage.begin(),storage.end(),REAL(inspect));
        SET_VECTOR_ELT(result,4,inspect);
    }
    R_API(Rf_setAttrib(result,R_NamesSymbol,names));
    return result;
}

#include "markov-alive-conditional.h"

} // namespace

extern "C" SEXP life_markov_grid_native(SEXP specification) {
    char message[1024];
    SEXP token = PROTECT(R_MakeUnwindCont());
    try {
        UnwindToken context(token);
        SEXP result = run_grid(specification);
        UNPROTECT(1);
        return result;
    }
    catch (const RJump& jump) { R_ContinueUnwind(jump.token); }
    catch (const std::bad_alloc&) { std::snprintf(message, sizeof(message), "%s", "Native Markov grid allocation failed; lower the explicit node budget or query output size."); }
    catch (const std::exception& error) { std::snprintf(message, sizeof(message), "%s", error.what()); }
    catch (...) { std::snprintf(message, sizeof(message), "%s", "Unknown native Markov grid error."); }
    UNPROTECT(1);
    Rf_error("%s", message);
    return R_NilValue;
}

extern "C" SEXP life_markov_monotone_queries(SEXP specification) {
    char message[1024];
    SEXP token=PROTECT(R_MakeUnwindCont());
    try {
        UnwindToken context(token);
        SEXP result=run_monotone_queries(specification);
        UNPROTECT(1); return result;
    }
    catch(const RJump& jump) { R_ContinueUnwind(jump.token); }
    catch(const std::bad_alloc&) { std::snprintf(message,sizeof(message),"%s","Scalar monotone allocation failed."); }
    catch(const std::exception& error) { std::snprintf(message,sizeof(message),"%s",error.what()); }
    catch(...) { std::snprintf(message,sizeof(message),"%s","Unknown scalar monotone error."); }
    UNPROTECT(1); Rf_error("%s",message); return R_NilValue;
}

extern "C" SEXP life_markov_alive_conditional(SEXP specification, SEXP gate) {
    char message[1024];
    SEXP token=PROTECT(R_MakeUnwindCont());
    try {
        UnwindToken context(token);
        if(TYPEOF(gate)!=LGLSXP || XLENGTH(gate)!=1 || LOGICAL(gate)[0]==NA_LOGICAL)
            fail("Invalid alive-time evaluator selection.");
        AliveQuery query(specification);
        SEXP answer=LOGICAL(gate)[0] ? alive_conditioned_kernels(specification,query) :
            alive_failure_window(specification,query);
        UNPROTECT(1); return answer;
    }
    catch(const RJump& jump) { R_ContinueUnwind(jump.token); }
    catch(const std::bad_alloc&) { std::snprintf(message,sizeof(message),"Alive-time kernel allocation failed."); }
    catch(const std::exception& error) { std::snprintf(message,sizeof(message),"%s",error.what()); }
    catch(...) { std::snprintf(message,sizeof(message),"Unknown alive-time calculation error."); }
    UNPROTECT(1); Rf_error("%s",message); return R_NilValue;
}
