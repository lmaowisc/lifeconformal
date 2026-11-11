#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include <limits.h>
#include <string.h>

/* Direct traversal for the complete numeric trees fitted by this package.
 * No backend model frames, node-label searches, or training rows are needed.
 * Traversal preserves tree summation order. Optional additive tables preserve
 * the finite ensemble while reassociating nonnegative floating-point sums. */
static SEXP forest_field(SEXP object, const char *name, SEXPTYPE type) {
    if (TYPEOF(object) != VECSXP)
        Rf_error("Invalid compiled transition forest.");
    SEXP names = Rf_getAttrib(object, R_NamesSymbol);
    if (TYPEOF(names) != STRSXP || XLENGTH(names) != XLENGTH(object))
        Rf_error("Compiled transition forest fields must be named.");
    for (R_xlen_t i = 0; i < XLENGTH(object); ++i) {
        if (STRING_ELT(names, i) != NA_STRING &&
            strcmp(CHAR(STRING_ELT(names, i)), name) == 0) {
            SEXP value = VECTOR_ELT(object, i);
            if (TYPEOF(value) != type)
                Rf_error("Invalid type for compiled transition forest field '%s'.", name);
            return value;
        }
    }
    Rf_error("Missing compiled transition forest field '%s'.", name);
    return R_NilValue;
}

/* Optional exact lookup; old serialized forests retain traversal. */
static SEXP forest_optional(SEXP object, const char *name) {
    SEXP names = Rf_getAttrib(object, R_NamesSymbol);
    for (R_xlen_t i = 0; i < XLENGTH(object); ++i)
        if (STRING_ELT(names, i) != NA_STRING &&
            strcmp(CHAR(STRING_ELT(names, i)), name) == 0)
            return VECTOR_ELT(object, i);
    return R_NilValue;
}

static int forest_bin(double x, const double *cuts, int count) {
    int low = 0, high = count;
    while (low < high) {
        int middle = low + (high - low) / 2;
        if (x < cuts[middle]) high = middle;
        else low = middle + 1;
    }
    return low;
}

static SEXP forest_rates(SEXP forest, SEXP features, SEXP rows,
                         SEXP cache, SEXP cache_rows) {
    SEXP roots = forest_field(forest, "roots", INTSXP);
    SEXP feature = forest_field(forest, "feature", INTSXP);
    SEXP left = forest_field(forest, "left", INTSXP);
    SEXP right = forest_field(forest, "right", INTSXP);
    SEXP direction = forest_field(forest, "direction", INTSXP);
    SEXP cut = forest_field(forest, "cut", REALSXP);
    SEXP value = forest_field(forest, "value", REALSXP);
    SEXP n_features = forest_field(forest, "n_features", INTSXP);
    R_xlen_t nodes = XLENGTH(feature), trees = XLENGTH(roots), nr = XLENGTH(rows);
    if (!trees || !nodes || nodes > INT_MAX || XLENGTH(n_features) != 1 ||
        XLENGTH(left) != nodes || XLENGTH(right) != nodes ||
        XLENGTH(direction) != nodes || XLENGTH(cut) != nodes ||
        XLENGTH(value) != nodes || TYPEOF(rows) != INTSXP)
        Rf_error("Invalid compiled transition forest dimensions.");
    SEXP dim = Rf_getAttrib(features, R_DimSymbol);
    if ((TYPEOF(features) != REALSXP && TYPEOF(features) != INTSXP) ||
        TYPEOF(dim) != INTSXP || XLENGTH(dim) != 2 ||
        INTEGER(dim)[1] != INTEGER(n_features)[0])
        Rf_error("Numeric features do not match the compiled transition forest.");
    const int n = INTEGER(dim)[0], p = INTEGER(dim)[1];
    if (n < 0 || p < 0 || XLENGTH(features) != (R_xlen_t)n * p)
        Rf_error("Invalid transition forest feature dimensions.");

    const int *root = INTEGER(roots), *column = INTEGER(feature);
    const int *lc = INTEGER(left), *rc = INTEGER(right), *parity = INTEGER(direction);
    const double *split = REAL(cut), *rate = REAL(value);
    int forward = 1;
    /* Immutable metadata is checked once per call, not at every visited node
     * for every draw. Root/child validity makes later pointer access safe;
     * the traversal counter below still detects malformed cycles. */
    for (R_xlen_t j = 0; j < nodes; ++j) {
        if (!R_FINITE(rate[j]) || rate[j] < 0)
            Rf_error("Invalid terminal transition forest rate.");
        if (column[j] != 0) {
            if (column[j] < 1 || column[j] > p ||
                (parity[j] != -1 && parity[j] != 1) || !R_FINITE(split[j]))
                Rf_error("Invalid numeric split in compiled transition forest.");
            if (lc[j] < 1 || lc[j] > nodes || rc[j] < 1 || rc[j] > nodes)
                Rf_error("Invalid child index or cycle in compiled transition forest.");
            if (lc[j] <= j + 1 || rc[j] <= j + 1) forward = 0;
        }
        if ((j & 65535) == 0) R_CheckUserInterrupt();
    }
    for (R_xlen_t b = 0; b < trees; ++b) {
        if (root[b] < 1 || root[b] > nodes)
            Rf_error("Invalid child index or cycle in compiled transition forest.");
    }
    SEXP numeric_features = PROTECT(Rf_coerceVector(features, REALSXP));
    int nprotect = 1;
    const double *x = REAL(numeric_features);
    for (R_xlen_t j = 0; j < XLENGTH(numeric_features); ++j) {
        if (!R_FINITE(x[j])) Rf_error("Transition forest features must be finite.");
        if ((j & 65535) == 0) R_CheckUserInterrupt();
    }
    int *selected = (int *) R_alloc(nr, sizeof(int));
    for (R_xlen_t i = 0; i < nr; ++i) {
        int row = INTEGER(rows)[i];
        if (row == NA_INTEGER || row < 1 || row > n)
            Rf_error("Invalid transition forest row index.");
        selected[i] = row - 1;
    }

    int *cached = NULL, *selected_cache = NULL;
    const int *tree_start = NULL, *tree_end = NULL, *box_start = NULL;
    const int *box_count = NULL, *box_feature = NULL;
    const double *box_lower = NULL, *box_upper = NULL;
    int cache_n = 0;
    if (cache != R_NilValue) {
        if (TYPEOF(cache) != ENVSXP || TYPEOF(cache_rows) != INTSXP ||
            XLENGTH(cache_rows) != n)
            Rf_error("Invalid private transition forest cache.");
        /* base::get0 with inherits=FALSE preserves local-frame semantics
         * using public APIs already present in R 4.1. Protect the call and a
         * promise/active-binding result through subsequent native allocations.
         * A private sentinel distinguishes missing nodes from a malformed
         * local NULL value without consulting any parent environment. */
        SEXP absent_nodes = PROTECT(Rf_allocVector(VECSXP, 1));
        SEXP nodes_name = PROTECT(Rf_mkString("nodes"));
        SEXP inherit = PROTECT(Rf_ScalarLogical(FALSE));
        SEXP get_call = PROTECT(Rf_lang5(Rf_findFun(Rf_install("get0"), R_BaseEnv),
                                        nodes_name, cache, inherit, absent_nodes));
        SET_TAG(CDDDR(get_call), Rf_install("inherits"));
        SET_TAG(CDR(CDDDR(get_call)), Rf_install("ifnotfound"));
        SEXP cached_nodes = PROTECT(Rf_eval(get_call, R_BaseEnv));
        nprotect += 5;
        if (cached_nodes == absent_nodes)
            Rf_error("Invalid private transition forest cache: missing local nodes.");
        SEXP cache_dim = Rf_getAttrib(cached_nodes, R_DimSymbol);
        if (TYPEOF(cached_nodes) != INTSXP || TYPEOF(cache_dim) != INTSXP ||
            XLENGTH(cache_dim) != 2 || INTEGER(cache_dim)[0] < 0 ||
            INTEGER(cache_dim)[1] != trees)
            Rf_error("Invalid private transition forest cache dimensions.");
        cache_n = INTEGER(cache_dim)[0];
        if (XLENGTH(cached_nodes) != (R_xlen_t)cache_n * trees)
            Rf_error("Invalid private transition forest cache length.");
        cached = INTEGER(cached_nodes);
        SEXP ts = forest_field(forest, "tree_start", INTSXP);
        SEXP te = forest_field(forest, "tree_end", INTSXP);
        SEXP bs = forest_field(forest, "box_start", INTSXP);
        SEXP bc = forest_field(forest, "box_count", INTSXP);
        SEXP bf = forest_field(forest, "box_feature", INTSXP);
        SEXP bl = forest_field(forest, "box_lower", REALSXP);
        SEXP bu = forest_field(forest, "box_upper", REALSXP);
        R_xlen_t bounds = XLENGTH(bf);
        if (XLENGTH(ts) != trees || XLENGTH(te) != trees ||
            XLENGTH(bs) != nodes || XLENGTH(bc) != nodes ||
            XLENGTH(bl) != bounds || XLENGTH(bu) != bounds)
            Rf_error("Invalid compiled transition forest cache bounds.");
        tree_start = INTEGER(ts); tree_end = INTEGER(te);
        box_start = INTEGER(bs); box_count = INTEGER(bc); box_feature = INTEGER(bf);
        box_lower = REAL(bl); box_upper = REAL(bu);
        for (R_xlen_t b = 0; b < trees; ++b) {
            if (tree_start[b] < 1 || tree_end[b] > nodes || tree_start[b] > tree_end[b] ||
                root[b] < tree_start[b] || root[b] > tree_end[b])
                Rf_error("Invalid compiled transition forest tree bounds.");
        }
        for (R_xlen_t j = 0; j < nodes; ++j) {
            if (box_start[j] < 0 || box_count[j] < 0 ||
                (R_xlen_t)box_start[j] + box_count[j] > bounds)
                Rf_error("Invalid compiled transition forest leaf bounds.");
        }
        for (R_xlen_t k = 0; k < bounds; ++k) {
            if (box_feature[k] < 1 || box_feature[k] > p ||
                ISNAN(box_lower[k]) || ISNAN(box_upper[k]) || box_lower[k] >= box_upper[k])
                Rf_error("Invalid compiled transition forest leaf box.");
            if ((k & 65535) == 0) R_CheckUserInterrupt();
        }
        selected_cache = (int *) R_alloc(nr, sizeof(int));
        for (R_xlen_t i = 0; i < nr; ++i) {
            int row = INTEGER(cache_rows)[selected[i]];
            if (row == NA_INTEGER || row < 1 || row > cache_n)
                Rf_error("Invalid private transition forest cache row.");
            selected_cache[i] = row - 1;
        }
    }
    /* Acyclic preorder is guaranteed by our compiler. Unusual old/private
     * arrays instead use checked traversal, including its cycle detection. */
    SEXP lookup = forward ? forest_optional(forest, "lookup") : R_NilValue;
    int *bins = NULL;
    const int *lf1 = NULL, *lf2 = NULL, *lo = NULL, *ls1 = NULL, *ls2 = NULL;
    const int *lm1 = NULL, *lm2 = NULL, *lmap = NULL, *lleaf = NULL;
    R_xlen_t leaf_length = 0, bin_variables = 0;
    const int *bin_counts = NULL, *bin_columns = NULL, *bin_starts = NULL;
    const double *bin_cuts = NULL;
    int lookup_complete = 1;
    if (lookup != R_NilValue) {
        SEXP variables = forest_field(lookup, "variables", INTSXP);
        SEXP starts = forest_field(lookup, "cut_start", INTSXP);
        SEXP counts = forest_field(lookup, "cut_count", INTSXP);
        SEXP cuts = forest_field(lookup, "cuts", REALSXP);
        SEXP f1 = forest_field(lookup, "feature1", INTSXP);
        SEXP f2 = forest_field(lookup, "feature2", INTSXP);
        SEXP offsets = forest_field(lookup, "offset", INTSXP);
        SEXP sizes1 = forest_field(lookup, "size1", INTSXP);
        SEXP sizes2 = forest_field(lookup, "size2", INTSXP);
        SEXP maps1 = forest_field(lookup, "map_start1", INTSXP);
        SEXP maps2 = forest_field(lookup, "map_start2", INTSXP);
        SEXP maps = forest_field(lookup, "map", INTSXP);
        SEXP leaves = forest_field(lookup, "leaf", INTSXP);
        R_xlen_t k = XLENGTH(variables), nc = XLENGTH(cuts), nm = XLENGTH(maps);
        bin_variables = k;
        leaf_length = XLENGTH(leaves);
        if (!k || k > p || XLENGTH(starts) != k || XLENGTH(counts) != k ||
            XLENGTH(f1) != trees || XLENGTH(f2) != trees || XLENGTH(offsets) != trees ||
            XLENGTH(sizes1) != trees || XLENGTH(sizes2) != trees ||
            XLENGTH(maps1) != trees || XLENGTH(maps2) != trees)
            Rf_error("Invalid compiled transition forest lookup dimensions.");
        const int *vars = INTEGER(variables), *cs = INTEGER(starts), *cc = INTEGER(counts);
        bin_counts = cc; bin_columns = vars; bin_starts = cs;
        const double *cuts_ptr = REAL(cuts);
        bin_cuts = cuts_ptr;
        lf1 = INTEGER(f1); lf2 = INTEGER(f2); lo = INTEGER(offsets);
        ls1 = INTEGER(sizes1); ls2 = INTEGER(sizes2);
        lm1 = INTEGER(maps1); lm2 = INTEGER(maps2);
        lmap = INTEGER(maps); lleaf = INTEGER(leaves);
        for (R_xlen_t j = 0; j < k; ++j) {
            if (vars[j] < 1 || vars[j] > p || (j && vars[j] <= vars[j - 1]) ||
                cs[j] < 0 || cc[j] < 1 || (R_xlen_t)cs[j] + cc[j] > nc)
                Rf_error("Invalid compiled transition forest lookup cuts.");
            for (int h = 0; h < cc[j]; ++h) {
                double z = cuts_ptr[(R_xlen_t)cs[j] + h];
                if (!R_FINITE(z) || (h && z <= cuts_ptr[(R_xlen_t)cs[j] + h - 1]))
                    Rf_error("Invalid compiled transition forest lookup cuts.");
            }
        }
        if (tree_start == NULL) {
            SEXP ts = forest_field(forest, "tree_start", INTSXP);
            SEXP te = forest_field(forest, "tree_end", INTSXP);
            if (XLENGTH(ts) != trees || XLENGTH(te) != trees)
                Rf_error("Invalid compiled transition forest tree bounds.");
            tree_start = INTEGER(ts); tree_end = INTEGER(te);
        }
        for (R_xlen_t b = 0; b < trees; ++b) {
            if (tree_start[b] < 1 || tree_end[b] > nodes || tree_start[b] > tree_end[b] ||
                root[b] < tree_start[b] || root[b] > tree_end[b])
                Rf_error("Invalid compiled transition forest tree bounds.");
            if (lo[b] == -1) {
                if (column[root[b] - 1] != 0) lookup_complete = 0;
                continue;
            }
            if (lo[b] < 0 || lf1[b] < 1 || lf1[b] > k || lf2[b] < 0 || lf2[b] > k ||
                ls1[b] < 1 || ls2[b] < 1 || (!lf2[b] && ls2[b] != 1) ||
                (R_xlen_t)ls1[b] * ls2[b] > leaf_length - (R_xlen_t)lo[b] ||
                lm1[b] < 0 || (R_xlen_t)lm1[b] + cc[lf1[b] - 1] + 1 > nm ||
                (lf2[b] && (lm2[b] < 0 || (R_xlen_t)lm2[b] + cc[lf2[b] - 1] + 1 > nm)))
                Rf_error("Invalid compiled transition forest lookup offsets.");
        }
        /* Query bins are temporary and capped at 16 MiB; large queries use
         * on-the-fly additive lookup when present, otherwise checked traversal.
         * Neither route allocates a row-by-tree prediction matrix. */
        if (nr <= (16 * 1024 * 1024 / sizeof(int)) / k) {
            bins = (int *) R_alloc(nr * k, sizeof(int));
            for (R_xlen_t j = 0; j < k; ++j) {
                for (R_xlen_t i = 0; i < nr; ++i)
                    bins[i + nr * j] = forest_bin(x[selected[i] + (R_xlen_t)n * (vars[j] - 1)],
                                                  cuts_ptr + cs[j], cc[j]);
                R_CheckUserInterrupt();
            }
        }
    }
    SEXP answer = PROTECT(Rf_allocVector(REALSXP, nr));
    ++nprotect;
    double *out = REAL(answer);
    for (R_xlen_t i = 0; i < nr; ++i) out[i] = 0;
    /* Trees using the same feature pair form an additive table on the exact
     * union of their split bins. This preserves the finite ensemble function,
     * with only nonnegative floating-point addition reassociated. Cached and
     * old serialized fits keep the checked original-order tree route. */
    SEXP additive = lookup != R_NilValue && lookup_complete && cached == NULL ?
        forest_optional(forest, "additive") : R_NilValue;
    if (additive != R_NilValue) {
        SEXP af1 = forest_field(additive, "feature1", INTSXP);
        SEXP af2 = forest_field(additive, "feature2", INTSXP);
        SEXP as1 = forest_field(additive, "size1", INTSXP);
        SEXP as2 = forest_field(additive, "size2", INTSXP);
        SEXP ao = forest_field(additive, "offset", INTSXP);
        SEXP av = forest_field(additive, "value", REALSXP);
        SEXP atrees = forest_field(additive, "trees", INTSXP);
        R_xlen_t groups = XLENGTH(af1), values = XLENGTH(av);
        if (!groups || XLENGTH(af2) != groups || XLENGTH(as1) != groups ||
            XLENGTH(as2) != groups || XLENGTH(ao) != groups ||
            XLENGTH(atrees) != 1 || INTEGER(atrees)[0] != trees)
            Rf_error("Invalid additive transition forest dimensions.");
        const int *a1 = INTEGER(af1), *a2 = INTEGER(af2), *s1 = INTEGER(as1);
        const int *s2 = INTEGER(as2), *offset = INTEGER(ao);
        const double *rates = REAL(av);
        for (R_xlen_t g = 0; g < groups; ++g) {
            if (a1[g] < 0 || a1[g] > bin_variables || a2[g] < 0 ||
                a2[g] > bin_variables || (!a1[g] && a2[g]) ||
                s1[g] < 1 || s2[g] < 1 ||
                (R_xlen_t)s1[g] != (a1[g] ? (R_xlen_t)bin_counts[a1[g] - 1] + 1 : 1) ||
                (R_xlen_t)s2[g] != (a2[g] ? (R_xlen_t)bin_counts[a2[g] - 1] + 1 : 1) ||
                offset[g] < 0 || (R_xlen_t)s1[g] * s2[g] > values - (R_xlen_t)offset[g])
                Rf_error("Invalid additive transition forest table.");
        }
        /* Validate only queried table entries: all compiled entries are
         * nonnegative finite sums, and every read is range-checked above. */
        for (R_xlen_t g = 0; g < groups; ++g) {
            R_CheckUserInterrupt();
            for (R_xlen_t i = 0; i < nr; ++i) {
                int one = 0, two = 0;
                if (a1[g]) {
                    int j = a1[g] - 1;
                    one = bins != NULL ? bins[i + nr * j] :
                        forest_bin(x[selected[i] + (R_xlen_t)n * (bin_columns[j] - 1)],
                                   bin_cuts + bin_starts[j], bin_counts[j]);
                }
                if (a2[g]) {
                    int j = a2[g] - 1;
                    two = bins != NULL ? bins[i + nr * j] :
                        forest_bin(x[selected[i] + (R_xlen_t)n * (bin_columns[j] - 1)],
                                   bin_cuts + bin_starts[j], bin_counts[j]);
                }
                double current = rates[(R_xlen_t)offset[g] + one + (R_xlen_t)s1[g] * two];
                if (!R_FINITE(current) || current < 0)
                    Rf_error("Invalid additive transition forest rate.");
                out[i] += current;
                if ((i & 65535) == 0) R_CheckUserInterrupt();
            }
        }
        for (R_xlen_t i = 0; i < nr; ++i) {
            out[i] /= (double)trees;
            if (!R_FINITE(out[i])) Rf_error("Nonfinite flexible transition rate.");
        }
        UNPROTECT(nprotect);
        return answer;
    }
    /* Each row still adds trees in their original order. Tree-major traversal
     * keeps node metadata and the column-major leaf indices local in memory. */
    for (R_xlen_t b = 0; b < trees; ++b) {
        R_CheckUserInterrupt();
        for (R_xlen_t i = 0; i < nr; ++i) {
            int row = selected[i], node = root[b];
            R_xlen_t cache_at = 0;
            if (cached != NULL) {
                cache_at = selected_cache[i] + (R_xlen_t)cache_n * b;
                int old = cached[cache_at];
                if (old != 0) {
                    if (old < 1 || old > nodes || old < tree_start[b] ||
                        old > tree_end[b] || column[old - 1] != 0)
                        Rf_error("Invalid leaf in private transition forest cache.");
                    int at = old - 1, inside = 1;
                    R_xlen_t end = (R_xlen_t)box_start[at] + box_count[at];
                    for (R_xlen_t k = box_start[at]; k < end; ++k) {
                        double current = x[row + (R_xlen_t)n * (box_feature[k] - 1)];
                        if (current < box_lower[k] || current >= box_upper[k]) {
                            inside = 0;
                            break;
                        }
                    }
                    if (inside) node = old;
                }
            }
            if (bins != NULL && lo[b] >= 0) {
                int a = lmap[(R_xlen_t)lm1[b] + bins[i + nr * (lf1[b] - 1)]];
                int d = lf2[b] ? lmap[(R_xlen_t)lm2[b] + bins[i + nr * (lf2[b] - 1)]] : 0;
                if (a < 0 || a >= ls1[b] || d < 0 || d >= ls2[b])
                    Rf_error("Invalid compiled transition forest lookup map.");
                node = lleaf[(R_xlen_t)lo[b] + a + (R_xlen_t)ls1[b] * d];
                if (node < tree_start[b] || node > tree_end[b] || column[node - 1] != 0)
                    Rf_error("Invalid compiled transition forest lookup leaf.");
            }
            R_xlen_t visited = 0;
            while (column[node - 1] != 0) {
                if (++visited > nodes)
                    Rf_error("Invalid child index or cycle in compiled transition forest.");
                int at = node - 1;
                double current = x[row + (R_xlen_t)n * (column[at] - 1)];
                int below = current < split[at];
                node = (parity[at] == -1 ? below : !below) ? lc[at] : rc[at];
                if ((visited & 65535) == 0) R_CheckUserInterrupt();
            }
            out[i] += rate[node - 1];
            if (cached != NULL) cached[cache_at] = node;
            if ((i & 4095) == 0) R_CheckUserInterrupt();
        }
    }
    for (R_xlen_t i = 0; i < nr; ++i) {
        out[i] /= (double)trees;
        if (!R_FINITE(out[i])) Rf_error("Nonfinite flexible transition rate.");
    }
    UNPROTECT(nprotect);
    return answer;
}

SEXP life_forest_rates(SEXP forest, SEXP features, SEXP rows) {
    return forest_rates(forest, features, rows, R_NilValue, R_NilValue);
}

SEXP life_forest_rates_cached(SEXP forest, SEXP features, SEXP rows,
                              SEXP cache, SEXP cache_rows) {
    return forest_rates(forest, features, rows, cache, cache_rows);
}

/* Construct the exact standardized design without two sweep/aperm arrays.
 * BLAS products and the nonlinear link remain in R, in their original order. */
SEXP life_flexible_design(SEXP features, SEXP rows, SEXP center, SEXP scale) {
    SEXP dim = Rf_getAttrib(features, R_DimSymbol);
    if ((TYPEOF(features) != REALSXP && TYPEOF(features) != INTSXP) || TYPEOF(dim) != INTSXP || XLENGTH(dim) != 2 ||
        TYPEOF(rows) != INTSXP || TYPEOF(center) != REALSXP || TYPEOF(scale) != REALSXP)
        Rf_error("Invalid flexible design inputs.");
    int n = INTEGER(dim)[0], p = INTEGER(dim)[1];
    R_xlen_t nr = XLENGTH(rows);
    if (n < 0 || p < 0 || p == INT_MAX || nr > INT_MAX ||
        XLENGTH(features) != (R_xlen_t)n * p ||
        XLENGTH(center) != p || XLENGTH(scale) != p)
        Rf_error("Invalid flexible design dimensions.");
    SEXP answer = PROTECT(Rf_allocMatrix(REALSXP, (int)nr, p + 1));
    for (R_xlen_t i = 0; i < nr; ++i) {
        if (INTEGER(rows)[i] == NA_INTEGER || INTEGER(rows)[i] < 1 || INTEGER(rows)[i] > n)
            Rf_error("Invalid flexible design row index.");
        REAL(answer)[i] = 1;
    }
    for (int j = 0; j < p; ++j) {
        if (!R_FINITE(REAL(center)[j]) || !R_FINITE(REAL(scale)[j]) || REAL(scale)[j] <= 0)
            Rf_error("Invalid flexible design standardization.");
        R_CheckUserInterrupt();
        for (R_xlen_t i = 0; i < nr; ++i) {
            R_xlen_t at = INTEGER(rows)[i] - 1 + (R_xlen_t)n * j;
            double x = TYPEOF(features) == REALSXP ? REAL(features)[at] :
                (INTEGER(features)[at] == NA_INTEGER ? NA_REAL : INTEGER(features)[at]);
            if (!R_FINITE(x)) Rf_error("Flexible design features must be finite.");
            double centered = x - REAL(center)[j];
            REAL(answer)[i + nr * (j + 1)] = centered / REAL(scale)[j];
            if ((i & 65535) == 0) R_CheckUserInterrupt();
        }
    }
    UNPROTECT(1);
    return answer;
}
