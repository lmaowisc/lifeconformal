#include "flexible-forest-view.h"

/* Each helper makes the binary64 store boundary explicit. R evaluates these
 * vector operations separately; contraction/reassociation would change ties. */
static double plus(double a, double b) { volatile double z = a + b; return z; }
static double minus(double a, double b) { volatile double z = a - b; return z; }
static double times(double a, double b) { volatile double z = a * b; return z; }
static double divide(double a, double b) { volatile double z = a / b; return z; }

SEXP life_forest_block_info(void) {
    SEXP out = PROTECT(Rf_allocVector(INTSXP, 4));
    INTEGER(out)[0] = DBL_MANT_DIG; INTEGER(out)[1] = LDBL_MANT_DIG;
    INTEGER(out)[2] = DBL_MAX_EXP; INTEGER(out)[3] = LDBL_MAX_EXP;
    UNPROTECT(1); return out;
}

static int scalar_int(SEXP x, const char *name) {
    if (TYPEOF(x) != INTSXP || XLENGTH(x) != 1 || INTEGER(x)[0] == NA_INTEGER)
        Rf_error("Invalid forest continuation %s.", name);
    return INTEGER(x)[0];
}

typedef struct {
    int n, p, d, states, edges, k, B, precision;
    long double workspace;
    const int *from, *to;
    const double *reward;
    forest_context *forest;
    double *constant;
    int *terminal;
} block_model;

static int prepare_model(SEXP adapter, int n, int B, int precision,
                         block_model *model, int materialize) {
    if (n < 0 || n > (64 * 1024 * 1024) / (64 * sizeof(double)) || B < 1 || n % B ||
        (precision != DBL_MANT_DIG && precision != LDBL_MANT_DIG)) return 0;
    SEXP fits = field(adapter, "fits", VECSXP), edges = field(adapter, "edges", VECSXP);
    SEXP problem = field(adapter, "problem", VECSXP);
    SEXP rewards = field(problem, "rewards", REALSXP);
    SEXP rd = Rf_getAttrib(rewards, R_DimSymbol);
    if (TYPEOF(rd) != INTSXP || XLENGTH(rd) != 2 || INTEGER(rd)[0] < 1 || INTEGER(rd)[1] < 1)
        Rf_error("Invalid forest continuation reward dimensions.");
    int states = INTEGER(rd)[0], d = INTEGER(rd)[1];
    SEXP from = field(edges, "from", INTSXP), to = field(edges, "to", INTSXP);
    R_xlen_t E = XLENGTH(from);
    if (E > INT_MAX || XLENGTH(to) != E || XLENGTH(fits) != E)
        Rf_error("Invalid forest continuation edge dimensions.");
    SEXP predictors = field(adapter, "n_predictors", INTSXP);
    int p = scalar_int(predictors, "predictors");
    if (p < 0 || p > INT_MAX - 3 - d) return 0;
    if ((long double)E * (sizeof(forest_context) + sizeof(double)) +
        (long double)states * sizeof(int) > 64 * 1024 * 1024) return 0;
    int k = 0;
    forest_context *views = materialize ? (forest_context *)R_alloc(E, sizeof(forest_context)) : NULL;
    double *constants = materialize ? (double *)R_alloc(E, sizeof(double)) : NULL;
    for (R_xlen_t e = 0; e < E; ++e) {
        if ((e & 65535) == 0) R_CheckUserInterrupt();
        if (INTEGER(from)[e] < 1 || INTEGER(from)[e] > states ||
            INTEGER(to)[e] < 1 || INTEGER(to)[e] > states)
            Rf_error("Invalid forest continuation edge.");
        SEXP fit = VECTOR_ELT(fits, e), constant = optional(fit, "constant");
        if (constant != R_NilValue) {
            if ((TYPEOF(constant) != REALSXP && TYPEOF(constant) != INTSXP) || XLENGTH(constant) != 1)
                Rf_error("Invalid forest continuation constant rate.");
            double c = Rf_asReal(constant);
            if (!R_FINITE(c) || c < 0) Rf_error("Nonfinite flexible transition rate.");
            if (materialize) constants[e] = c;
        } else {
            SEXP forest = optional(fit, "compiled");
            if (forest == R_NilValue) return 0;
            forest_context local;
            if (!prepare_forest(forest, &local)) return 0;
            if (local.p != p + 3 + d) Rf_error("Invalid forest continuation feature count.");
            if (local.variables > k) k = local.variables;
            if (materialize) { views[e] = local; constants[e] = -1; }
        }
    }
    /* Conservative total: existing bank + original/normalized/output state + reusable numeric
     * and index buffers. This bound is checked before any path is advanced. */
    long double bytes = (long double)n *
        (64 * sizeof(double) + 3 * ((5 + (long double)d) * sizeof(double) + 2 * sizeof(int)) +
         ((long double)p + 3 + d + E + 4) * sizeof(double) +
         ((long double)k + 8) * sizeof(int)) +
         (long double)E * (sizeof(forest_context) + sizeof(double)) +
         (long double)states * sizeof(int);
    if (bytes > 64 * 1024 * 1024) return 0;
    if (model) model->workspace = bytes;
    if (materialize) {
        int *terminal = (int *)R_alloc(states, sizeof(int));
        for (int s = 0; s < states; ++s) terminal[s] = 1;
        for (R_xlen_t e = 0; e < E; ++e) terminal[INTEGER(from)[e] - 1] = 0;
        SEXP absorbing = field(problem, "absorbing", INTSXP);
        for (R_xlen_t i = 0; i < XLENGTH(absorbing); ++i) {
            int s = INTEGER(absorbing)[i];
            if (s < 1 || s > states) Rf_error("Invalid forest continuation absorbing state.");
            terminal[s - 1] = 1;
        }
        model->n = n; model->B = B; model->p = p; model->d = d;
        model->states = states; model->edges = (int)E; model->k = k;
        model->precision = precision; model->from = INTEGER(from); model->to = INTEGER(to);
        model->reward = REAL(rewards); model->forest = views;
        model->constant = constants; model->terminal = terminal;
    }
    return 1;
}

SEXP life_forest_block_supported(SEXP adapter, SEXP n, SEXP B, SEXP precision) {
    return Rf_ScalarLogical(prepare_model(adapter, scalar_int(n, "rows"),
        scalar_int(B, "draws"), scalar_int(precision, "precision"), NULL, 0));
}

static SEXP copied_vector(SEXP path, const char *name, SEXPTYPE type, R_xlen_t length) {
    SEXP z = field(path, name, type);
    if (XLENGTH(z) != length) Rf_error("Invalid forest continuation state dimensions.");
    return z;
}

static SEXP trace_step(SEXP path, int iteration, const int *active, int na,
        const int *selected, int ns, const int *representative, const double *rates,
        int edges, const double *total, const double *wait, const int *event) {
    const char *names[] = {"path", "iteration", "active", "selected", "rates", "total", "wait", "event"};
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 8));
    SEXP nms = PROTECT(Rf_allocVector(STRSXP, 8));
    for (int j = 0; j < 8; ++j) SET_STRING_ELT(nms, j, Rf_mkChar(names[j]));
    Rf_setAttrib(out, R_NamesSymbol, nms);
    SET_VECTOR_ELT(out, 0, Rf_duplicate(path));
    SET_VECTOR_ELT(out, 1, Rf_ScalarInteger(iteration));
    SET_VECTOR_ELT(out, 2, Rf_allocVector(INTSXP, na));
    SET_VECTOR_ELT(out, 3, Rf_allocVector(INTSXP, ns));
    SET_VECTOR_ELT(out, 4, Rf_allocMatrix(REALSXP, na, edges));
    SET_VECTOR_ELT(out, 5, Rf_allocVector(REALSXP, na));
    SET_VECTOR_ELT(out, 6, Rf_allocVector(REALSXP, na));
    SET_VECTOR_ELT(out, 7, Rf_allocVector(LGLSXP, na));
    for (int i = 0; i < na; ++i) {
        int r = active[i];
        INTEGER(VECTOR_ELT(out, 2))[i] = r + 1;
        REAL(VECTOR_ELT(out, 5))[i] = total[r]; REAL(VECTOR_ELT(out, 6))[i] = wait[r];
        LOGICAL(VECTOR_ELT(out, 7))[i] = event[r];
        for (int e = 0; e < edges; ++e)
            REAL(VECTOR_ELT(out, 4))[i + (R_xlen_t)na * e] = rates[representative[r] + (R_xlen_t)ns * e];
    }
    for (int i = 0; i < ns; ++i) INTEGER(VECTOR_ELT(out, 3))[i] = selected[i] + 1;
    UNPROTECT(2); return out;
}

static SEXP forest_advance_impl(SEXP adapter, SEXP x, SEXP group, SEXP path,
        SEXP bank, SEXP controls, int trusted_bank) {
    int n = (int)XLENGTH(group);
    if (TYPEOF(group) != INTSXP || XLENGTH(group) > INT_MAX)
        Rf_error("Invalid forest continuation profile indices.");
    int B = scalar_int(field(controls, "draws", INTSXP), "draws");
    int precision = scalar_int(field(controls, "precision", INTSXP), "precision");
    int steps = scalar_int(field(controls, "steps", INTSXP), "steps");
    int maximum = scalar_int(field(controls, "max_iter", INTSXP), "iteration budget");
    int iteration = scalar_int(field(path, "iteration", INTSXP), "iteration");
    int tracing = scalar_int(field(controls, "trace", INTSXP), "trace");
    if (steps < 1 || steps > 32 || iteration < 0 || maximum < 0 ||
        steps > 32 - iteration % 32 || iteration > INT_MAX - steps || (tracing != 0 && tracing != 1))
        Rf_error("Invalid forest continuation iteration range.");
    block_model model;
    if (!prepare_model(adapter, n, B, precision, &model, 1)) return R_NilValue;
    const int p = model.p, d = model.d, E = model.edges, states = model.states;
    SEXP xd = Rf_getAttrib(x, R_DimSymbol), bd = Rf_getAttrib(bank, R_DimSymbol);
    if ((TYPEOF(x) != REALSXP && TYPEOF(x) != INTSXP) || TYPEOF(xd) != INTSXP || XLENGTH(xd) != 2 ||
        INTEGER(xd)[0] < 0 || INTEGER(xd)[1] != p ||
        XLENGTH(x) != (R_xlen_t)INTEGER(xd)[0] * p ||
        TYPEOF(bank) != REALSXP || TYPEOF(bd) != INTSXP || XLENGTH(bd) != 2 ||
        INTEGER(bd)[0] != n || INTEGER(bd)[1] != 64 || XLENGTH(bank) != (R_xlen_t)n * 64)
        Rf_error("Invalid forest continuation input dimensions.");
    int nx = INTEGER(xd)[0], real_x = TYPEOF(x) == REALSXP;
    for (int i = 0; i < n; ++i)
        if (INTEGER(group)[i] < 1 || INTEGER(group)[i] > nx ||
            INTEGER(group)[i] != INTEGER(group)[(i / B) * B] ||
            (i > 0 && i % B == 0 && INTEGER(group)[i] <= INTEGER(group)[i - B]))
            Rf_error("Invalid forest continuation profile layout.");
    if (!trusted_bank) {
    for (R_xlen_t i = 0; i < XLENGTH(bank); ++i)
        if (!R_FINITE(REAL(bank)[i]) || REAL(bank)[i] < 0 || REAL(bank)[i] > 1)
            Rf_error("Invalid forest continuation bank.");
    }
    SEXP problem = field(adapter, "problem", VECSXP), options = field(adapter, "options", VECSXP);
    double tau = Rf_asReal(optional(problem, "tau"));
    double step = Rf_asReal(optional(options, "step_size"));
    SEXP budget = optional(options, "max_events");
    double max_events = Rf_asReal(budget);
    if (!R_FINITE(tau) || tau <= 0 || !R_FINITE(step) || step <= 0 ||
        !R_FINITE(max_events) || max_events < 0 || max_events != floor(max_events))
        Rf_error("Invalid forest continuation controls.");

    SEXP result = PROTECT(Rf_duplicate(path));
    double *tt = REAL(copied_vector(result, "tt", REALSXP, n));
    int *state = INTEGER(copied_vector(result, "state", INTSXP, n));
    double *entry = REAL(copied_vector(result, "entry", REALSXP, n));
    double *cell = REAL(copied_vector(result, "cell", REALSXP, n));
    double *cell_end = REAL(copied_vector(result, "cell_end", REALSXP, n));
    double *previous = REAL(copied_vector(result, "previous", REALSXP, n));
    double *accrued = REAL(copied_vector(result, "accrued", REALSXP, (R_xlen_t)n * d));
    int *count = INTEGER(copied_vector(result, "count", INTSXP, n));
    int *iter = INTEGER(copied_vector(result, "iteration", INTSXP, 1));
    for (int i = 0; i < n; ++i)
        if (!R_FINITE(tt[i]) || tt[i] < 0 || tt[i] > tau || state[i] < 1 || state[i] > states ||
            !R_FINITE(entry[i]) || !R_FINITE(cell[i]) || cell[i] < 0 || cell[i] != floor(cell[i]) ||
            ISNAN(cell_end[i]) || !R_FINITE(previous[i]) || previous[i] < 0 ||
            previous[i] != floor(previous[i]) || count[i] < 0 || count[i] > iteration)
            Rf_error("Invalid forest continuation path state.");
    if (tracing && (long double)n * (d + (long double)E + 20) * 8 * steps > 4 * 1024 * 1024)
        Rf_error("Forest continuation trace budget exceeded.");
    SEXP traces = PROTECT(Rf_allocVector(VECSXP, tracing ? steps : 0));
    int trace_count = 0;
    int *active = (int *)R_alloc(n, sizeof(int)), *selected = (int *)R_alloc(n, sizeof(int));
    int *representative = (int *)R_alloc(n, sizeof(int)), *rows = (int *)R_alloc(n, sizeof(int));
    int *first = (int *)R_alloc(n / B, sizeof(int)), *event = (int *)R_alloc(n, sizeof(int));
    int *destination = (int *)R_alloc(n, sizeof(int));
    int *bins = (int *)R_alloc((R_xlen_t)n * model.k, sizeof(int));
    double *features = (double *)R_alloc((R_xlen_t)n * (p + 3 + d), sizeof(double));
    double *rates = (double *)R_alloc((R_xlen_t)n * E, sizeof(double));
    double *h = (double *)R_alloc(n, sizeof(double)), *total = (double *)R_alloc(n, sizeof(double));
    double *wait = (double *)R_alloc(n, sizeof(double)), *advanced = (double *)R_alloc(n, sizeof(double));

    for (int tick = 0; tick < steps; ++tick) {
        R_CheckUserInterrupt();
        int na = 0;
        for (int r = 0; r < n; ++r) {
            if (model.terminal[state[r] - 1] && tt[r] < tau) {
                for (int j = 0; j < d; ++j)
                    accrued[r + (R_xlen_t)n * j] = plus(accrued[r + (R_xlen_t)n * j],
                        times(model.reward[state[r] - 1 + (R_xlen_t)states * j], minus(tau, tt[r])));
                tt[r] = tau;
            }
            if (tt[r] < tau) active[na++] = r;
        }
        if (!na) break;
        if (++iteration > maximum) Rf_error("Flexible continuation iteration budget exceeded.");
        *iter = iteration;
        int column = (iteration - 1) % 32;
        for (int i = 0; i < na; ++i) {
            int r = active[i]; h[r] = minus(fmin(cell_end[r], tau), tt[r]);
            if (h[r] <= 0) Rf_error("The flexible intensity-cell clock failed to advance.");
        }
        for (int i = 0; i < n / B; ++i) first[i] = -1;
        int ns = 0;
        for (int i = 0; i < na; ++i) {
            int r = active[i], profile = r / B;
            if (!count[r] && first[profile] >= 0) { representative[r] = first[profile]; continue; }
            representative[r] = ns;
            if (!count[r]) first[profile] = ns;
            selected[ns++] = r;
        }
        for (int s = 0; s < ns; ++s) {
            int r = selected[s];
            for (int j = 0; j < p; ++j) {
                R_xlen_t at = INTEGER(group)[r] - 1 + (R_xlen_t)nx * j;
                /* Element access avoids materializing a whole ALTREP baseline
                 * matrix, and integer inputs need no whole-matrix conversion. */
                double value;
                if (real_x) value = REAL_ELT(x, at);
                else {
                    int z = INTEGER_ELT(x, at);
                    value = z == NA_INTEGER ? NA_REAL : (double)z;
                }
                features[s + (R_xlen_t)ns * j] = value;
            }
            double dur = times(plus(cell[r], .5), step);
            double midpoint = plus(entry[r], dur), shift = minus(midpoint, tt[r]);
            features[s + (R_xlen_t)ns * p] = midpoint;
            features[s + (R_xlen_t)ns * (p + 1)] = dur;
            features[s + (R_xlen_t)ns * (p + 2)] = previous[r];
            for (int j = 0; j < d; ++j)
                features[s + (R_xlen_t)ns * (p + 3 + j)] = plus(accrued[r + (R_xlen_t)n * j],
                    times(model.reward[state[r] - 1 + (R_xlen_t)states * j], shift));
        }
        for (R_xlen_t i = 0; i < (R_xlen_t)ns * (p + 3 + d); ++i)
            if (!R_FINITE(features[i])) Rf_error("Invalid flexible-model history features.");
        for (R_xlen_t i = 0; i < (R_xlen_t)ns * E; ++i) rates[i] = 0;
        for (int e = 0; e < E; ++e) {
            R_CheckUserInterrupt();
            int nr = 0;
            for (int s = 0; s < ns; ++s) if (state[selected[s]] == model.from[e]) rows[nr++] = s;
            if (!nr) continue;
            if (model.constant[e] >= 0) {
                for (int i = 0; i < nr; ++i) rates[rows[i] + (R_xlen_t)ns * e] = model.constant[e];
                continue;
            }
            forest_context *view = model.forest + e;
            for (int j = 0; j < view->variables; ++j)
                for (int i = 0; i < nr; ++i)
                    bins[i + (R_xlen_t)nr * j] = bin(features[rows[i] + (R_xlen_t)ns * (view->columns[j] - 1)],
                        view->cuts + view->starts[j], view->counts[j]);
            for (R_xlen_t g = 0; g < view->groups; ++g) {
                if ((g & 63) == 0) R_CheckUserInterrupt();
                for (int i = 0; i < nr; ++i) {
                    int one = view->first[g] ? bins[i + (R_xlen_t)nr * (view->first[g] - 1)] : 0;
                    int two = view->second[g] ? bins[i + (R_xlen_t)nr * (view->second[g] - 1)] : 0;
                    R_xlen_t at = rows[i] + (R_xlen_t)ns * e;
                    rates[at] = plus(rates[at], view->values[(R_xlen_t)view->offsets[g] + one + (R_xlen_t)view->sizes[g] * two]);
                }
            }
            for (int i = 0; i < nr; ++i) {
                R_xlen_t at = rows[i] + (R_xlen_t)ns * e;
                rates[at] = divide(rates[at], (double)view->trees);
                if (!R_FINITE(rates[at])) Rf_error("Nonfinite flexible transition rate.");
            }
        }
        for (int i = 0; i < na; ++i) {
            int r = active[i], s = representative[r];
            if (precision != DBL_MANT_DIG) {
                long double sum = 0;
                for (int e = 0; e < E; ++e) sum += rates[s + (R_xlen_t)ns * e];
                volatile double rounded = (double)sum; total[r] = rounded;
            } else {
                double sum = 0;
                for (int e = 0; e < E; ++e) sum = plus(sum, rates[s + (R_xlen_t)ns * e]);
                total[r] = sum;
            }
            if (!R_FINITE(total[r])) Rf_error("Total outgoing transition intensity exceeds the numerical range.");
            wait[r] = total[r] > 0 ? divide(-log(REAL(bank)[r + (R_xlen_t)n * (2 * column)]), total[r]) : R_PosInf;
            event[r] = wait[r] < h[r];
            double elapsed = fmin(wait[r], h[r]);
            for (int j = 0; j < d; ++j)
                accrued[r + (R_xlen_t)n * j] = plus(accrued[r + (R_xlen_t)n * j],
                    times(model.reward[state[r] - 1 + (R_xlen_t)states * j], elapsed));
            advanced[r] = plus(tt[r], elapsed);
        }
        for (int i = 0; i < na; ++i) {
            int r = active[i];
            if (event[r] && advanced[r] <= tt[r])
                Rf_error("A simulated waiting time is below the numerical resolution of the time scale.");
        }
        for (int i = 0; i < na; ++i) {
            int r = active[i]; tt[r] = advanced[r];
            if (!event[r]) {
                tt[r] = fmin(cell_end[r], tau); cell[r] = plus(cell[r], 1);
                cell_end[r] = plus(entry[r], times(plus(cell[r], 1), step));
            } else {
                if (count[r] == INT_MAX)
                    Rf_error("Flexible continuation exceeded max_events; inspect fitted intensities or increase the explicit budget.");
                ++count[r];
                if (count[r] > max_events)
                    Rf_error("Flexible continuation exceeded max_events; inspect fitted intensities or increase the explicit budget.");
            }
        }
        for (int i = 0; i < na; ++i) if (event[active[i]]) {
            int r = active[i], s = representative[r];
            double threshold = times(REAL(bank)[r + (R_xlen_t)n * (2 * column + 1)], total[r]);
            double cumulative = 0; destination[r] = 0;
            for (int e = 0; e < E; ++e) {
                cumulative = plus(cumulative, rates[s + (R_xlen_t)ns * e]);
                if (!destination[r] && threshold < cumulative) destination[r] = model.to[e];
            }
            if (!destination[r]) Rf_error("Failed competing-transition probability lookup.");
        }
        for (int i = 0; i < na; ++i) if (event[active[i]]) {
            int r = active[i]; state[r] = destination[r]; entry[r] = tt[r]; cell[r] = 0;
            cell_end[r] = plus(entry[r], step); previous[r] = plus(previous[r], 1);
        }
        if (tracing) SET_VECTOR_ELT(traces, trace_count++, trace_step(result, iteration,
            active, na, selected, ns, representative, rates, E, total, wait, event));
    }
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
    SEXP names = PROTECT(Rf_allocVector(STRSXP, 2));
    SET_STRING_ELT(names, 0, Rf_mkChar("path")); SET_STRING_ELT(names, 1, Rf_mkChar("trace"));
    Rf_setAttrib(out, R_NamesSymbol, names);
    SET_VECTOR_ELT(out, 0, result);
    SEXP trace_out = PROTECT(Rf_allocVector(VECSXP, trace_count));
    for (int i = 0; i < trace_count; ++i) SET_VECTOR_ELT(trace_out, i, VECTOR_ELT(traces, i));
    SET_VECTOR_ELT(out, 1, trace_out);
    UNPROTECT(5); return out;
}

/* Ordinary matrices retain the accepted checks, including the complete scan. */
SEXP life_forest_block_advance(SEXP adapter, SEXP x, SEXP group, SEXP path,
        SEXP bank, SEXP controls) {
    return forest_advance_impl(adapter, x, group, path, bank, controls, 0);
}

#include <Rmath.h>
#include <R_ext/Random.h>
#include <R_ext/Utils.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#define BANK_MAGIC UINT64_C(0x4c46424b50525631)
#define BANK_COLUMNS 64
#define BANK_META_RESERVE 4096
#define BANK_MAX_LIVE 2
enum { BANK_BUILDING = 1, BANK_READY = 2, BANK_INVALID = 3 };
typedef struct {
    uint64_t magic;
    SEXP bank, group, bundle, adapter, weak;
    int B, profiles, block, rows, completed, state, busy;
    int audit, copied, precision;
    long double reserved;
    size_t payload;
} bank_owner;
static bank_owner *owners[BANK_MAX_LIVE];
static int bank_active_calls = 0;

static SEXP bank_tag(void) { return Rf_install("life_generated_forest_bank_v1"); }
static int slot_of(const void *address) {
    for (int i = 0; i < BANK_MAX_LIVE; ++i)
        if (address && owners[i] == address) return i;
    return -1;
}
static int free_slot(void) {
    for (int i = 0; i < BANK_MAX_LIVE; ++i) if (!owners[i]) return i;
    return -1;
}
static int integer_scalar(SEXP value, const char *label) {
    if (TYPEOF(value) != INTSXP || XLENGTH(value) != 1 || INTEGER(value)[0] == NA_INTEGER)
        Rf_error("Invalid forest uniform bank %s.", label);
    return INTEGER(value)[0];
}
static bank_owner *owner_of(SEXP token) {
    if (TYPEOF(token) != EXTPTRSXP || R_ExternalPtrTag(token) != bank_tag())
        Rf_error("Foreign forest uniform bank handle.");
    void *address = R_ExternalPtrAddr(token);
    int index = slot_of(address);
    if (index < 0) Rf_error("Released, restored or foreign forest uniform bank handle.");
    bank_owner *owner = owners[index]; /* Registry membership precedes dereference. */
    SEXP bundle = R_ExternalPtrProtected(token);
    if (owner->magic != BANK_MAGIC || bundle != owner->bundle || TYPEOF(bundle) != VECSXP ||
        XLENGTH(bundle) != 4 || VECTOR_ELT(bundle, 0) != owner->bank ||
        VECTOR_ELT(bundle, 1) != owner->group || VECTOR_ELT(bundle, 2) != owner->adapter ||
        VECTOR_ELT(bundle, 3) != owner->weak || R_WeakRefKey(owner->weak) != token)
        Rf_error("Invalid forest uniform bank ownership.");
    SEXP bank = owner->bank;
    if (TYPEOF(bank) != REALSXP || TYPEOF(owner->group) != INTSXP ||
        XLENGTH(owner->group) != owner->rows)
        Rf_error("Invalid forest uniform bank ownership.");
    SEXP dim = Rf_getAttrib(bank, R_DimSymbol);
    if (TYPEOF(dim) != INTSXP || XLENGTH(dim) != 2 || INTEGER(dim)[0] != owner->rows ||
        INTEGER(dim)[1] != BANK_COLUMNS || XLENGTH(bank) != (R_xlen_t)owner->rows * BANK_COLUMNS)
        Rf_error("Invalid forest uniform bank storage.");
    return owner;
}
static void finalize_bank(SEXP token) {
    if (TYPEOF(token) != EXTPTRSXP) return;
    int index = slot_of(R_ExternalPtrAddr(token));
    bank_owner *owner = index < 0 ? NULL : owners[index];
    if (index >= 0) owners[index] = NULL;
    R_ClearExternalPtr(token);
    R_SetExternalPtrProtected(token, R_NilValue);
    if (owner) {
        owner->state = BANK_INVALID;
        owner->magic = 0;
        free(owner);
    }
}

/* The accepted W includes the new bank. Reserve another bank/group payload
 * for replacement overlap (or the actual larger live payload), stored group
 * copies, metadata, and optionally one audit copy plus the accepted 4 MiB trace
 * cap. This is bounded call-local live payload, not process RSS or arbitrary
 * retained audit outputs. Every check runs before allocating/generating draws. */
/* Pure C ledger: safe to repeat after weakref allocation runs finalizers. */
static long double bank_required_bytes(int n, int audit, long double workspace) {
    long double payload = (long double)n * BANK_COLUMNS * sizeof(double);
    long double groups = (long double)n * sizeof(int), live = 0;
    for (int i = 0; i < BANK_MAX_LIVE; ++i)
        if (owners[i]) live += owners[i]->payload + (long double)owners[i]->rows * sizeof(int) +
            owners[i]->audit * (owners[i]->payload + 4 * 1024 * 1024);
    long double overlap = payload + groups;
    if (live > overlap) overlap = live;
    return workspace + overlap + groups + BANK_META_RESERVE +
        audit * (payload + 4 * 1024 * 1024);
}
static int bank_reservation(SEXP adapter, SEXP group, int B, int precision,
        int audit, double budget, long double *required, long double *workspace) {
    if (TYPEOF(group) != INTSXP || XLENGTH(group) < 1 || XLENGTH(group) > 16000 ||
        B < 1 || B > 1000 || XLENGTH(group) % B || XLENGTH(group) / B > 16 ||
        (audit != 0 && audit != 1) || !R_FINITE(budget) || budget < 0 ||
        budget > 64 * 1024 * 1024) return 0;
    int n = (int)XLENGTH(group);
    for (int i = 0; i < n; ++i)
        if (INTEGER(group)[i] < 1 || INTEGER(group)[i] != INTEGER(group)[(i / B) * B] ||
            (i > 0 && i % B == 0 && INTEGER(group)[i] <= INTEGER(group)[i - B])) return 0;
    block_model model;
    if (!prepare_model(adapter, n, B, precision, &model, 0)) return 0;
    *required = bank_required_bytes(n, audit, model.workspace);
    if (workspace) *workspace = model.workspace;
    return *required <= budget && free_slot() >= 0 &&
        sizeof(bank_owner) * BANK_MAX_LIVE <= BANK_META_RESERVE;
}
static SEXP bank_supported_impl(SEXP adapter, SEXP group, SEXP draws,
        SEXP precision, SEXP budget, SEXP audit) {
    long double required;
    if (TYPEOF(budget) != REALSXP || XLENGTH(budget) != 1) return Rf_ScalarLogical(0);
    return Rf_ScalarLogical(bank_reservation(adapter, group,
        integer_scalar(draws, "draw count"), integer_scalar(precision, "precision"),
        integer_scalar(audit, "audit allowance"), REAL(budget)[0], &required, NULL));
}
static SEXP bank_create_impl(SEXP adapter, SEXP group, SEXP draws,
        SEXP precision, SEXP block, SEXP budget, SEXP audit) {
    int B = integer_scalar(draws, "draw count");
    int prec = integer_scalar(precision, "precision");
    int block_id = integer_scalar(block, "block index");
    int allow_audit = integer_scalar(audit, "audit allowance");
    long double required, workspace;
    if (block_id < 1 || block_id > INT_MAX / 32 + 1 || TYPEOF(budget) != REALSXP ||
        XLENGTH(budget) != 1 || !bank_reservation(adapter, group, B, prec,
            allow_audit, REAL(budget)[0], &required, &workspace))
        Rf_error("Forest uniform-bank consumer dimensions, model or payload budget unsupported.");
    double max_bytes = REAL(budget)[0];
    int n = (int)XLENGTH(group), G = n / B;
    size_t payload = (size_t)n * BANK_COLUMNS * sizeof(double);
    SEXP token = PROTECT(R_MakeExternalPtr(NULL, bank_tag(), R_NilValue));
    SEXP bundle = PROTECT(Rf_allocVector(VECSXP, 4));
    SEXP bank = PROTECT(Rf_allocMatrix(REALSXP, n, BANK_COLUMNS));
    double *initial = REAL(bank);
    for (R_xlen_t i = 0; i < (R_xlen_t)n * BANK_COLUMNS; ++i) initial[i] = 0.0;
    SET_VECTOR_ELT(bundle, 0, bank);
    SET_VECTOR_ELT(bundle, 1, Rf_duplicate(group));
    SET_VECTOR_ELT(bundle, 2, adapter);
    R_SetExternalPtrProtected(token, bundle);
    /* Recheck live payload after allocation/finalizers, before publication. */
    if (!bank_reservation(adapter, group, B, prec, allow_audit, REAL(budget)[0], &required, &workspace))
        Rf_error("Forest uniform-bank consumer payload eligibility changed before generation.");
    /* Reserve the protection slot before creating the callback. After the
     * weakref exists, every failure explicitly runs it before leaving C. */
    PROTECT_INDEX weak_index;
    PROTECT_WITH_INDEX(R_NilValue, &weak_index);
    SEXP weak = R_MakeWeakRefC(token, R_NilValue, finalize_bank, TRUE);
    REPROTECT(weak, weak_index);
    SET_VECTOR_ELT(bundle, 3, weak);
    required = bank_required_bytes(n, allow_audit, workspace);
    int index = free_slot();
    if (index < 0 || required > max_bytes) {
        R_RunWeakRefFinalizer(weak);
        Rf_error("Forest uniform-bank payload eligibility changed before generation.");
    }
    bank_owner *owner = (bank_owner *)calloc(1, sizeof(bank_owner));
    if (!owner) {
        R_RunWeakRefFinalizer(weak);
        Rf_error("Unable to allocate forest uniform-bank descriptor.");
    }
    owner->magic = BANK_MAGIC; owner->bank = bank;
    owner->group = VECTOR_ELT(bundle, 1); owner->bundle = bundle; owner->weak = weak;
    owner->B = B; owner->profiles = G; owner->block = block_id;
    owner->rows = n; owner->state = BANK_BUILDING; owner->payload = payload;
    owner->reserved = required; owner->audit = allow_audit;
    owner->adapter = adapter; owner->precision = prec;
    owners[index] = owner;
    R_SetExternalPtrAddr(token, owner);
    UNPROTECT(4);
    return token;
}

typedef struct {
    bank_owner *owner;
    double *values;
    SEXP rng_unwind;
    int profile, acquired, committed, written;
} fill_context;

/* Each bounded profile uses exactly the running R runif(0,1) sequence. The
 * caller retains its original profile/block flex_with_seed scope. */
static SEXP generate_profile(void *data) {
    fill_context *context = (fill_context *)data;
    bank_owner *owner = context->owner;
    R_CheckUserInterrupt();
    GetRNGstate();
    context->acquired = 1;
    int count = owner->B * BANK_COLUMNS;
    for (int j = 0; j < count; ++j) {
        double value = runif(0.0, 1.0);
        if (!R_FINITE(value) || value < 0 || value > 1)
            Rf_error("Invalid generated forest uniform-bank value.");
        int row = context->profile * owner->B + j / BANK_COLUMNS;
        int column = j % BANK_COLUMNS;
        context->values[row + (R_xlen_t)owner->rows * column] = value;
        context->written += 1;
    }
    return R_NilValue;
}
static void finish_rng(void *data, Rboolean jump) {
    fill_context *context = (fill_context *)data;
    bank_owner *owner = context->owner;
    owner->state = BANK_INVALID; /* Before any cleanup operation can jump. */
    if (context->acquired) {
        context->acquired = 0; /* Clear BEFORE Put: no re-entry/double commit. */
        PutRNGstate();
        context->committed = 1;
    }
    (void)jump;
}
static SEXP complete_fill(void *data) {
    fill_context *context = (fill_context *)data;
    bank_owner *owner = context->owner;
    R_UnwindProtect(generate_profile, context, finish_rng, context, context->rng_unwind);
    if (!context->committed || context->written != owner->B * BANK_COLUMNS)
        Rf_error("Forest uniform bank fill did not commit completely.");
    owner->completed += 1;
    owner->state = BANK_BUILDING;
    return R_NilValue;
}
static void unlock_fill(void *data, Rboolean jump) {
    fill_context *context = (fill_context *)data;
    if (jump) context->owner->state = BANK_INVALID;
    context->owner->busy = 0;
}
static SEXP bank_fill_impl(SEXP token, SEXP profile) {
    (void)owner_of(token);
    int index = integer_scalar(profile, "profile index");
    SEXP rng_unwind = PROTECT(R_MakeUnwindCont());
    SEXP owner_unwind = PROTECT(R_MakeUnwindCont());
    /* Allocation can run finalizers. Reacquire before retaining the owner. */
    bank_owner *owner = owner_of(token);
    if (owner->busy || owner->state != BANK_BUILDING || index != owner->completed + 1 || index > owner->profiles)
        Rf_error("Forest uniform bank fill is incomplete, duplicated, out of order or sealed.");
    fill_context context = {owner, REAL(owner->bank), rng_unwind, index - 1, 0, 0, 0};
    owner->busy = 1;
    owner->state = BANK_INVALID; /* Non-consumable before GetRNGstate. */
    R_UnwindProtect(complete_fill, &context, unlock_fill, &context, owner_unwind);
    UNPROTECT(2);
    return R_NilValue;
}
static SEXP bank_seal_impl(SEXP token) {
    bank_owner *owner = owner_of(token);
    if (owner->state != BANK_BUILDING || owner->completed != owner->profiles)
        Rf_error("Forest uniform bank is not complete and cannot seal.");
    owner->state = BANK_READY;
    return R_NilValue;
}
static SEXP bank_copy_impl(SEXP token) {
    bank_owner *owner = owner_of(token);
    if (owner->state != BANK_READY || !owner->audit || owner->copied)
        Rf_error("Forest uniform bank audit copy is unavailable or already consumed.");
    owner->copied = 1; /* Reserve the one permitted independent copy before allocation. */
    SEXP bank = PROTECT(owner->bank);
    SEXP out = PROTECT(Rf_duplicate(bank));
    UNPROTECT(2);
    return out;
}
static SEXP bank_release_impl(SEXP token) {
    if (TYPEOF(token) != EXTPTRSXP || R_ExternalPtrTag(token) != bank_tag())
        Rf_error("Foreign forest uniform bank handle.");
    if (R_ExternalPtrAddr(token)) {
        bank_owner *owner = owner_of(token);
        if (bank_active_calls > 1)
            Rf_error("Cannot release a forest uniform bank during another bank operation.");
        if (owner->busy) Rf_error("Forest uniform bank is active and cannot release.");
        /* Clears the callback before destroying the owner. A later GC, including
         * after DLL unload, can no longer invoke this library's finalizer. */
        R_RunWeakRefFinalizer(owner->weak);
    }
    return R_NilValue;
}
static SEXP bank_inspect_impl(SEXP token) {
    bank_owner *owner = owner_of(token);
    const char *names[] = {"state", "completed", "draws", "profiles", "block", "payload", "reserved_payload"};
    double values[] = {owner->state, owner->completed, owner->B, owner->profiles, owner->block,
        (double)owner->payload, (double)owner->reserved};
    SEXP out = PROTECT(Rf_allocVector(REALSXP, 7));
    SEXP labels = PROTECT(Rf_allocVector(STRSXP, 7));
    for (int i = 0; i < 7; ++i) {
        REAL(out)[i] = values[i]; SET_STRING_ELT(labels, i, Rf_mkChar(names[i]));
    }
    Rf_setAttrib(out, R_NamesSymbol, labels); UNPROTECT(2); return out;
}
static SEXP bank_live_impl(void) {
    int count = 0;
    for (int i = 0; i < BANK_MAX_LIVE; ++i) if (owners[i]) count += 1;
    return Rf_ScalarInteger(count);
}

/* This is the only route to trusted_bank=1. The native registry and hidden
 * protected storage prove generation/completion; R attributes cannot grant it. */
typedef struct {
    SEXP token, adapter, x, group, path, controls;
    bank_owner *owner;
    int was_ready;
} consume_context;
static SEXP consume_ready(void *data) {
    consume_context *context = (consume_context *)data;
    bank_owner *owner = context->owner;
    if (!context->was_ready) Rf_error("Forest uniform-bank consumer bank is not READY.");
    int B = scalar_int(field(context->controls, "draws", INTSXP), "draws");
    int precision = scalar_int(field(context->controls, "precision", INTSXP), "precision");
    int iteration = scalar_int(field(context->path, "iteration", INTSXP), "iteration");
    int tracing = scalar_int(field(context->controls, "trace", INTSXP), "trace");
    if (context->adapter != owner->adapter || precision != owner->precision ||
        B != owner->B || TYPEOF(context->group) != INTSXP ||
        XLENGTH(context->group) != owner->rows || iteration < 0 || iteration % 32 ||
        iteration / 32 + 1 != owner->block || (tracing && !owner->audit))
        Rf_error("Forest uniform-bank consumer bank identity or block mismatch.");
    for (int i = 0; i < owner->rows; ++i)
        if (INTEGER(context->group)[i] != INTEGER(owner->group)[i])
            Rf_error("Forest uniform-bank consumer bank profile identity mismatch.");
    return forest_advance_impl(context->adapter, context->x, context->group,
        context->path, owner->bank, context->controls, 1);
}
static void consume_release(void *data, Rboolean jump) {
    consume_context *context = (consume_context *)data;
    /* Clear the weakref callback, then registry/address, before freeing. */
    R_RunWeakRefFinalizer(context->owner->weak);
    (void)jump;
}
static SEXP bank_advance_impl(SEXP adapter, SEXP x, SEXP group, SEXP path,
        SEXP token, SEXP controls) {
    (void)owner_of(token);
    SEXP unwind = PROTECT(R_MakeUnwindCont());
    bank_owner *owner = owner_of(token);
    if (owner->busy) Rf_error("Forest uniform-bank consumer bank is already active.");
    consume_context context = {token, adapter, x, group, path, controls, owner,
        owner->state == BANK_READY};
    owner->busy = 1; owner->state = BANK_INVALID;
    SEXP out = R_UnwindProtect(consume_ready, &context, consume_release, &context, unwind);
    UNPROTECT(1);
    return out;
}

/* The cleanup context is established without an allocating continuation.
 * Mark every bank entry active before its first allocation, including creation
 * and read-only diagnostics, so reentrant DLL unload cannot strand a callback. */
typedef struct { int operation, entered; SEXP args[7]; } bank_call;
static SEXP bank_call_body(void *data) {
    bank_call *call = (bank_call *)data;
    if (bank_active_calls == INT_MAX) Rf_error("Too many active forest uniform-bank operations.");
    ++bank_active_calls;
    call->entered = 1;
    switch (call->operation) {
    case 0: return bank_supported_impl(call->args[0], call->args[1], call->args[2], call->args[3], call->args[4], call->args[5]);
    case 1: return bank_create_impl(call->args[0], call->args[1], call->args[2], call->args[3], call->args[4], call->args[5], call->args[6]);
    case 2: return bank_fill_impl(call->args[0], call->args[1]);
    case 3: return bank_seal_impl(call->args[0]);
    case 4: return bank_copy_impl(call->args[0]);
    case 5: return bank_release_impl(call->args[0]);
    case 6: return bank_inspect_impl(call->args[0]);
    case 7: return bank_live_impl();
    case 8: return bank_advance_impl(call->args[0], call->args[1], call->args[2], call->args[3], call->args[4], call->args[5]);
    default: Rf_error("Invalid forest uniform-bank operation.");
    }
    return R_NilValue;
}
static void bank_call_cleanup(void *data) {
    bank_call *call = (bank_call *)data;
    if (call->entered) {
        call->entered = 0;
        --bank_active_calls;
    }
}
SEXP life_forest_bank_supported(SEXP a0, SEXP a1, SEXP a2, SEXP a3, SEXP a4, SEXP a5) {
    bank_call call = {0, 0, {a0, a1, a2, a3, a4, a5}};
    return R_ExecWithCleanup(bank_call_body, &call, bank_call_cleanup, &call);
}
SEXP life_forest_bank_create(SEXP a0, SEXP a1, SEXP a2, SEXP a3, SEXP a4, SEXP a5, SEXP a6) {
    bank_call call = {1, 0, {a0, a1, a2, a3, a4, a5, a6}};
    return R_ExecWithCleanup(bank_call_body, &call, bank_call_cleanup, &call);
}
SEXP life_forest_bank_fill(SEXP a0, SEXP a1) {
    bank_call call = {2, 0, {a0, a1}};
    return R_ExecWithCleanup(bank_call_body, &call, bank_call_cleanup, &call);
}
SEXP life_forest_bank_seal(SEXP a0) {
    bank_call call = {3, 0, {a0}};
    return R_ExecWithCleanup(bank_call_body, &call, bank_call_cleanup, &call);
}
SEXP life_forest_bank_copy(SEXP a0) {
    bank_call call = {4, 0, {a0}};
    return R_ExecWithCleanup(bank_call_body, &call, bank_call_cleanup, &call);
}
SEXP life_forest_bank_release(SEXP a0) {
    bank_call call = {5, 0, {a0}};
    return R_ExecWithCleanup(bank_call_body, &call, bank_call_cleanup, &call);
}
SEXP life_forest_bank_inspect(SEXP a0) {
    bank_call call = {6, 0, {a0}};
    return R_ExecWithCleanup(bank_call_body, &call, bank_call_cleanup, &call);
}
SEXP life_forest_bank_live(void) {
    bank_call call = {7, 0, {R_NilValue}};
    return R_ExecWithCleanup(bank_call_body, &call, bank_call_cleanup, &call);
}
SEXP life_forest_bank_advance(SEXP a0, SEXP a1, SEXP a2, SEXP a3, SEXP a4, SEXP a5) {
    bank_call call = {8, 0, {a0, a1, a2, a3, a4, a5}};
    return R_ExecWithCleanup(bank_call_body, &call, bank_call_cleanup, &call);
}

/* The loader checks every native registry before draining either one. */
void life_forest_bank_unload_check(void) {
    if (bank_active_calls) Rf_error("Cannot unload an active forest uniform-bank operation.");
    for (int i = 0; i < BANK_MAX_LIVE; ++i)
        if (owners[i] && owners[i]->busy)
            Rf_error("Cannot unload an active forest uniform-bank operation.");
}
void life_forest_bank_finalize_all(void) {
    for (int i = 0; i < BANK_MAX_LIVE; ++i)
        if (owners[i]) R_RunWeakRefFinalizer(owners[i]->weak);
}
