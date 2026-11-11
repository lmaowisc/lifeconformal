/* Validated read-only additive forest views, local to one continuation call. */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <math.h>

#define CONTEXT_MAGIC UINT64_C(0x4c46435458563141)
typedef struct {
    uint64_t magic;
    int p, variables;
    R_xlen_t trees, groups;
    const int *columns, *starts, *counts;
    const double *cuts;
    const int *first, *second, *sizes, *offsets;
    const double *values;
} forest_context;
static SEXP field(SEXP x, const char *name, SEXPTYPE type) {
    if (TYPEOF(x) != VECSXP) Rf_error("Invalid forest continuation view source.");
    SEXP names = Rf_getAttrib(x, R_NamesSymbol);
    if (TYPEOF(names) != STRSXP || XLENGTH(names) != XLENGTH(x))
        Rf_error("Forest continuation view fields must be named.");
    for (R_xlen_t i = 0; i < XLENGTH(x); ++i)
        if (STRING_ELT(names, i) != NA_STRING && strcmp(CHAR(STRING_ELT(names, i)), name) == 0) {
            SEXP out = VECTOR_ELT(x, i);
            if (TYPEOF(out) != type) Rf_error("Invalid forest continuation view field '%s'.", name);
            return out;
        }
    Rf_error("Missing forest continuation view field '%s'.", name);
    return R_NilValue;
}

static SEXP optional(SEXP x, const char *name) {
    if (TYPEOF(x) != VECSXP) Rf_error("Invalid forest continuation view source.");
    SEXP names = Rf_getAttrib(x, R_NamesSymbol);
    if (TYPEOF(names) != STRSXP || XLENGTH(names) != XLENGTH(x))
        Rf_error("Forest continuation view fields must be named.");
    for (R_xlen_t i = 0; i < XLENGTH(x); ++i)
        if (STRING_ELT(names, i) != NA_STRING && strcmp(CHAR(STRING_ELT(names, i)), name) == 0)
            return VECTOR_ELT(x, i);
    return R_NilValue;
}

static int prepare_forest(SEXP forest, forest_context *context) {
    SEXP lookup = optional(forest, "lookup"), additive = optional(forest, "additive");
    if (lookup == R_NilValue || additive == R_NilValue) return 0;
    SEXP roots = field(forest, "roots", INTSXP), feature = field(forest, "feature", INTSXP);
    SEXP left = field(forest, "left", INTSXP), right = field(forest, "right", INTSXP);
    SEXP direction = field(forest, "direction", INTSXP), cut = field(forest, "cut", REALSXP);
    SEXP value = field(forest, "value", REALSXP), n_features = field(forest, "n_features", INTSXP);
    R_xlen_t nodes = XLENGTH(feature), trees = XLENGTH(roots);
    if (!trees || !nodes || nodes > INT_MAX || XLENGTH(n_features) != 1 ||
        INTEGER(n_features)[0] < 0 || XLENGTH(left) != nodes || XLENGTH(right) != nodes ||
        XLENGTH(direction) != nodes || XLENGTH(cut) != nodes || XLENGTH(value) != nodes)
        Rf_error("Invalid forest continuation view dimensions.");
    int p = INTEGER(n_features)[0], forward = 1;
    const int *root = INTEGER(roots), *column = INTEGER(feature);
    const int *lc = INTEGER(left), *rc = INTEGER(right), *parity = INTEGER(direction);
    const double *split = REAL(cut), *node_value = REAL(value);
    for (R_xlen_t j = 0; j < nodes; ++j) {
        if (!R_FINITE(node_value[j]) || node_value[j] < 0)
            Rf_error("Invalid forest continuation view node rate.");
        if (column[j] != 0) {
            if (column[j] < 1 || column[j] > p ||
                (parity[j] != -1 && parity[j] != 1) || !R_FINITE(split[j]))
                Rf_error("Invalid forest continuation view split.");
            if (lc[j] < 1 || lc[j] > nodes || rc[j] < 1 || rc[j] > nodes)
                Rf_error("Invalid forest continuation view child index.");
            if (lc[j] <= j + 1 || rc[j] <= j + 1) forward = 0;
        }
        if ((j & 65535) == 0) R_CheckUserInterrupt();
    }
    for (R_xlen_t b = 0; b < trees; ++b)
        if (root[b] < 1 || root[b] > nodes) Rf_error("Invalid forest continuation view root.");
    /* Unusual old/private metadata stays on the original checked traversal. */
    if (!forward) return 0;

    SEXP variables = field(lookup, "variables", INTSXP), starts = field(lookup, "cut_start", INTSXP);
    SEXP counts = field(lookup, "cut_count", INTSXP), cuts = field(lookup, "cuts", REALSXP);
    SEXP f1 = field(lookup, "feature1", INTSXP), f2 = field(lookup, "feature2", INTSXP);
    SEXP offsets = field(lookup, "offset", INTSXP), sizes1 = field(lookup, "size1", INTSXP);
    SEXP sizes2 = field(lookup, "size2", INTSXP), maps1 = field(lookup, "map_start1", INTSXP);
    SEXP maps2 = field(lookup, "map_start2", INTSXP), maps = field(lookup, "map", INTSXP);
    SEXP leaves = field(lookup, "leaf", INTSXP);
    R_xlen_t k = XLENGTH(variables), nc = XLENGTH(cuts), nm = XLENGTH(maps), nl = XLENGTH(leaves);
    if (!k || k > p || XLENGTH(starts) != k || XLENGTH(counts) != k ||
        XLENGTH(f1) != trees || XLENGTH(f2) != trees || XLENGTH(offsets) != trees ||
        XLENGTH(sizes1) != trees || XLENGTH(sizes2) != trees ||
        XLENGTH(maps1) != trees || XLENGTH(maps2) != trees)
        Rf_error("Invalid forest continuation view lookup dimensions.");
    const int *vars = INTEGER(variables), *cs = INTEGER(starts), *cc = INTEGER(counts);
    const int *lf1 = INTEGER(f1), *lf2 = INTEGER(f2), *lo = INTEGER(offsets);
    const int *ls1 = INTEGER(sizes1), *ls2 = INTEGER(sizes2), *lm1 = INTEGER(maps1), *lm2 = INTEGER(maps2);
    const double *cut_values = REAL(cuts);
    for (R_xlen_t j = 0; j < k; ++j) {
        if (vars[j] < 1 || vars[j] > p || (j && vars[j] <= vars[j - 1]) ||
            cs[j] < 0 || cc[j] < 1 || (R_xlen_t)cs[j] + cc[j] > nc)
            Rf_error("Invalid forest continuation view cuts.");
        for (int h = 0; h < cc[j]; ++h) {
            double z = cut_values[(R_xlen_t)cs[j] + h];
            if (!R_FINITE(z) || (h && z <= cut_values[(R_xlen_t)cs[j] + h - 1]))
                Rf_error("Invalid forest continuation view cuts.");
            if ((h & 65535) == 0) R_CheckUserInterrupt();
        }
    }
    SEXP tree_start = field(forest, "tree_start", INTSXP), tree_end = field(forest, "tree_end", INTSXP);
    if (XLENGTH(tree_start) != trees || XLENGTH(tree_end) != trees)
        Rf_error("Invalid forest continuation view tree dimensions.");
    const int *ts = INTEGER(tree_start), *te = INTEGER(tree_end);
    for (R_xlen_t b = 0; b < trees; ++b) {
        if (ts[b] < 1 || te[b] > nodes || ts[b] > te[b] || root[b] < ts[b] || root[b] > te[b])
            Rf_error("Invalid forest continuation view tree bounds.");
        if (lo[b] == -1) {
            if (column[root[b] - 1] != 0) return 0;
            continue;
        }
        if (lo[b] < 0 || lf1[b] < 1 || lf1[b] > k || lf2[b] < 0 || lf2[b] > k ||
            ls1[b] < 1 || ls2[b] < 1 || (!lf2[b] && ls2[b] != 1) ||
            (R_xlen_t)ls1[b] * ls2[b] > nl - (R_xlen_t)lo[b] ||
            lm1[b] < 0 || (R_xlen_t)lm1[b] + cc[lf1[b] - 1] + 1 > nm ||
            (lf2[b] && (lm2[b] < 0 || (R_xlen_t)lm2[b] + cc[lf2[b] - 1] + 1 > nm)))
            Rf_error("Invalid forest continuation view lookup offsets.");
    }
    SEXP af1 = field(additive, "feature1", INTSXP), af2 = field(additive, "feature2", INTSXP);
    SEXP as1 = field(additive, "size1", INTSXP), as2 = field(additive, "size2", INTSXP);
    SEXP ao = field(additive, "offset", INTSXP), av = field(additive, "value", REALSXP);
    SEXP atrees = field(additive, "trees", INTSXP);
    R_xlen_t groups = XLENGTH(af1), values = XLENGTH(av);
    if (!groups || XLENGTH(af2) != groups || XLENGTH(as1) != groups ||
        XLENGTH(as2) != groups || XLENGTH(ao) != groups ||
        XLENGTH(atrees) != 1 || INTEGER(atrees)[0] != trees)
        Rf_error("Invalid forest continuation view additive dimensions.");
    const int *a1 = INTEGER(af1), *a2 = INTEGER(af2), *s1 = INTEGER(as1), *s2 = INTEGER(as2), *off = INTEGER(ao);
    const double *rates = REAL(av);
    for (R_xlen_t g = 0; g < groups; ++g) {
        if (a1[g] < 0 || a1[g] > k || a2[g] < 0 || a2[g] > k || (!a1[g] && a2[g]) ||
            s1[g] < 1 || s2[g] < 1 ||
            (R_xlen_t)s1[g] != (a1[g] ? (R_xlen_t)cc[a1[g] - 1] + 1 : 1) ||
            (R_xlen_t)s2[g] != (a2[g] ? (R_xlen_t)cc[a2[g] - 1] + 1 : 1) ||
            off[g] < 0 || (R_xlen_t)s1[g] * s2[g] > values - (R_xlen_t)off[g])
            Rf_error("Invalid forest continuation view additive table.");
    }
    for (R_xlen_t i = 0; i < values; ++i) {
        if (!R_FINITE(rates[i]) || rates[i] < 0) Rf_error("Invalid forest continuation view additive rate.");
        if ((i & 65535) == 0) R_CheckUserInterrupt();
    }

    context->magic = CONTEXT_MAGIC; context->p = p; context->variables = (int)k;
    context->trees = trees; context->groups = groups;
    context->columns = vars; context->starts = cs; context->counts = cc; context->cuts = cut_values;
    context->first = a1; context->second = a2; context->sizes = s1;
    context->offsets = off; context->values = rates;
    return 1;
}

static int bin(double x, const double *cuts, int count) {
    int low = 0, high = count;
    while (low < high) {
        int middle = low + (high - low) / 2;
        if (x < cuts[middle]) high = middle;
        else low = middle + 1;
    }
    return low;
}

