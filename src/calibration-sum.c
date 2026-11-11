/* Call-local calibration accumulator; never included in a returned object.
 * Match R colSums/colMeans: binary64 products, input-row-order LDOUBLE sums,
 * and extended division before the final binary64 conversion. R's LDOUBLE
 * macro is private; the caller supplies the runtime mantissa precision.
 */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Memory.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

typedef union { double narrow; long double wide; } sum_cell;
typedef struct sum_state {
    struct sum_state *next;
    SEXP finalizer;
    int ncols;
    int extended;
    R_xlen_t nrows;
    sum_cell sum[];
} sum_state;

static sum_state *states = NULL;
static int active_sum_calls = 0;

static SEXP state_tag(void) {
    return Rf_install("lifeconformal.calibration_sum.v1");
}

/* Check registry membership before dereferencing an untrusted pointer. */
static sum_state *validated_state(SEXP state) {
    if (TYPEOF(state) != EXTPTRSXP || R_ExternalPtrTag(state) != state_tag())
        Rf_error("state must be a live extended-sum external pointer");
    void *address = R_ExternalPtrAddr(state);
    for (sum_state *s = states; s != NULL; s = s->next)
        if ((void *) s == address) return s;
    Rf_error("state must be a live extended-sum external pointer");
    return NULL;
}

static void finalize_state(SEXP state) {
    if (TYPEOF(state) != EXTPTRSXP || R_ExternalPtrTag(state) != state_tag())
        return;
    void *address = R_ExternalPtrAddr(state);
    sum_state **link = &states;
    while (*link != NULL) {
        if ((void *) *link == address) {
            sum_state *old = *link;
            *link = old->next;
            free(old);
            break;
        }
        link = &(*link)->next;
    }
    R_ClearExternalPtr(state);
    R_SetExternalPtrProtected(state, R_NilValue);
}

static int integer_scalar(SEXP x, const char *name) {
    if (TYPEOF(x) != INTSXP || XLENGTH(x) != 1 ||
        Rf_getAttrib(x, R_DimSymbol) != R_NilValue || INTEGER(x)[0] == NA_INTEGER)
        Rf_error("%s must be one nonmissing integer", name);
    return INTEGER(x)[0];
}

static SEXP calibration_sum_create_impl(SEXP ncols, SEXP precision_digits) {
    int nc = integer_scalar(ncols, "ncols");
    int digits = integer_scalar(precision_digits, "precision_digits");
    if (nc < 0) Rf_error("ncols must be nonnegative");
    if (digits != DBL_MANT_DIG && digits != LDBL_MANT_DIG)
        Rf_error("runtime precision does not match compiled double/long double");
    if ((size_t) nc > (SIZE_MAX - sizeof(sum_state)) / sizeof(sum_cell))
        Rf_error("ncols exceeds allocation limit");
    SEXP state = PROTECT(R_MakeExternalPtr(NULL, state_tag(), R_NilValue));
    /* Retain the registration so close/unload can clear its callback before
     * the DLL disappears. Allocate all R objects before the descriptor. */
    PROTECT_INDEX finalizer_index;
    PROTECT_WITH_INDEX(R_NilValue, &finalizer_index);
    SEXP finalizer = R_MakeWeakRefC(state, R_NilValue, finalize_state, TRUE);
    REPROTECT(finalizer, finalizer_index);
    R_SetExternalPtrProtected(state, finalizer);
    sum_state *s = (sum_state *) calloc(1, sizeof(sum_state) +
                                      (size_t) nc * sizeof(sum_cell));
    if (s == NULL) {
        R_RunWeakRefFinalizer(finalizer);
        Rf_error("cannot allocate extended-sum state");
    }
    s->finalizer = finalizer;
    s->ncols = nc;
    s->extended = digits != DBL_MANT_DIG;
    s->nrows = 0;
    s->next = states;
    states = s;
    R_SetExternalPtrAddr(state, s);
    /* Explicit initializations also cover non-all-bits-zero float platforms. */
    for (int j = 0; j < nc; ++j) {
        if ((j & 65535) == 0) R_CheckUserInterrupt();
        if (s->extended) s->sum[j].wide = 0.;
        else s->sum[j].narrow = 0.;
    }
    UNPROTECT(2);
    return state;
}

static SEXP calibration_sum_add_impl(SEXP state, SEXP values, SEXP weights) {
    sum_state *s = validated_state(state);
    SEXP dim = Rf_getAttrib(values, R_DimSymbol);
    if (TYPEOF(values) != REALSXP || TYPEOF(dim) != INTSXP || XLENGTH(dim) != 2)
        Rf_error("values must be a double matrix");
    int nr = INTEGER(dim)[0], nc = INTEGER(dim)[1];
    if (nr < 0 || nc < 0 || nc != s->ncols ||
        (nr > 0 && (R_xlen_t) nc > R_XLEN_T_MAX / nr))
        Rf_error("values dimensions do not match state");
    if ((R_xlen_t) nr * nc != XLENGTH(values))
        Rf_error("values dimensions do not match state");
    if (TYPEOF(weights) != REALSXP || XLENGTH(weights) != nr ||
        Rf_getAttrib(weights, R_DimSymbol) != R_NilValue)
        Rf_error("row_weights must be a double vector with one weight per row");
    if (s->nrows > R_XLEN_T_MAX - nr) Rf_error("row count exceeds R limit");
    R_CheckUserInterrupt();
    const double *w = REAL(weights), *x = REAL(values);
    for (int i = 0; i < nr; ++i) {
        if (i && (i & 65535) == 0) R_CheckUserInterrupt();
        if (!R_FINITE(w[i])) Rf_error("row_weights must be finite");
    }
    /* Atomic on validation errors. Temporary memory depends on columns only. */
    long double *next = R_allocLD((size_t) nc);
    unsigned int products_since_interrupt = 0;
    for (int j = 0; j < nc; ++j) {
        if (j && (j & 65535) == 0) R_CheckUserInterrupt();
        const double *column = x + (R_xlen_t) nr * j;
        if (s->extended) {
            long double total = s->sum[j].wide;
            for (int i = 0; i < nr; ++i) {
                if (++products_since_interrupt == 65536) {
                    R_CheckUserInterrupt();
                    products_since_interrupt = 0;
                }
                if (!R_FINITE(column[i])) Rf_error("values must be finite");
                /* Force the same rounding as R's materialized x * weights. */
                volatile double product = column[i] * w[i];
                if (!R_FINITE(product)) Rf_error("weighted products must be finite");
                total += product;
            }
            next[j] = total;
        } else {
            double total = s->sum[j].narrow;
            for (int i = 0; i < nr; ++i) {
                if (++products_since_interrupt == 65536) {
                    R_CheckUserInterrupt();
                    products_since_interrupt = 0;
                }
                if (!R_FINITE(column[i])) Rf_error("values must be finite");
                volatile double product = column[i] * w[i];
                if (!R_FINITE(product)) Rf_error("weighted products must be finite");
                volatile double rounded = total + product;
                total = rounded;
            }
            next[j] = total;
        }
    }
    for (int j = 0; j < nc; ++j) {
        if (s->extended) s->sum[j].wide = next[j];
        else s->sum[j].narrow = (double) next[j];
    }
    s->nrows += nr;
    return R_NilValue;
}

static SEXP calibration_sum_value_impl(SEXP state, SEXP divisor) {
    sum_state *s = validated_state(state);
    if (TYPEOF(divisor) != REALSXP || XLENGTH(divisor) != 1 ||
        Rf_getAttrib(divisor, R_DimSymbol) != R_NilValue ||
        !R_FINITE(REAL(divisor)[0]) || REAL(divisor)[0] < 0. ||
        (REAL(divisor)[0] == 0. && s->nrows != 0))
        Rf_error("divisor must be finite positive double, or zero for empty state");
    double denominator = REAL(divisor)[0];
    SEXP result = PROTECT(Rf_allocVector(REALSXP, s->ncols));
    for (int j = 0; j < s->ncols; ++j) {
        if ((j & 65535) == 0) R_CheckUserInterrupt();
        if (s->extended) REAL(result)[j] = (double) (s->sum[j].wide / denominator);
        else REAL(result)[j] = s->sum[j].narrow / denominator;
    }
    UNPROTECT(1);
    return result;
}

static SEXP calibration_sum_reset_impl(SEXP state) {
    sum_state *s = validated_state(state);
    /* Reset is just an atomic commit, like add's final store loop. */
    R_CheckUserInterrupt();
    for (int j = 0; j < s->ncols; ++j) {
        if (s->extended) s->sum[j].wide = 0.;
        else s->sum[j].narrow = 0.;
    }
    s->nrows = 0;
    return R_NilValue;
}

/* Explicit release keeps native memory bounded across repeated calls. */
static SEXP calibration_sum_close_impl(SEXP state) {
    sum_state *s = validated_state(state);
    /* Another finalizer can hold an alias even while this key is protected.
     * Do not let a nested explicit close free storage used by an outer call. */
    if (active_sum_calls > 1)
        Rf_error("cannot close an extended-sum state during another sum operation");
    R_RunWeakRefFinalizer(s->finalizer);
    return R_NilValue;
}

static SEXP calibration_sum_info_impl(void) {
    int count = 0;
    for (sum_state *s = states; s != NULL; s = s->next) ++count;
    SEXP out = PROTECT(Rf_allocVector(INTSXP, 7));
    int info[] = {DBL_MANT_DIG, LDBL_MANT_DIG, DBL_MAX_EXP, LDBL_MAX_EXP,
                  (int) sizeof(double), (int) sizeof(long double), count};
    for (int i = 0; i < 7; ++i) INTEGER(out)[i] = info[i];
    UNPROTECT(1);
    return out;
}

/* R_ExecWithCleanup establishes a non-allocating cleanup context before the
 * body increments this counter. This also covers constructor allocations,
 * validation failures and interrupts without changing accumulation order. */
enum sum_operation { SUM_CREATE, SUM_ADD, SUM_VALUE, SUM_RESET, SUM_CLOSE, SUM_INFO };
typedef struct {
    enum sum_operation operation;
    SEXP first, second, third;
    int entered;
} sum_call;

static SEXP run_sum_body(void *data) {
    sum_call *call = (sum_call *) data;
    if (active_sum_calls == INT_MAX)
        Rf_error("too many active extended-sum operations");
    ++active_sum_calls;
    call->entered = 1;
    switch (call->operation) {
    case SUM_CREATE: return calibration_sum_create_impl(call->first, call->second);
    case SUM_ADD: return calibration_sum_add_impl(call->first, call->second, call->third);
    case SUM_VALUE: return calibration_sum_value_impl(call->first, call->second);
    case SUM_RESET: return calibration_sum_reset_impl(call->first);
    case SUM_CLOSE: return calibration_sum_close_impl(call->first);
    case SUM_INFO: return calibration_sum_info_impl();
    }
    Rf_error("invalid extended-sum operation");
    return R_NilValue;
}

static void finish_sum_call(void *data) {
    sum_call *call = (sum_call *) data;
    if (call->entered) {
        call->entered = 0;
        --active_sum_calls;
    }
}

static SEXP run_sum_call(enum sum_operation operation, SEXP first, SEXP second, SEXP third) {
    sum_call call = {operation, first, second, third, 0};
    return R_ExecWithCleanup(run_sum_body, &call, finish_sum_call, &call);
}

SEXP life_calibration_sum_create(SEXP ncols, SEXP precision_digits) {
    return run_sum_call(SUM_CREATE, ncols, precision_digits, R_NilValue);
}
SEXP life_calibration_sum_add(SEXP state, SEXP values, SEXP weights) {
    return run_sum_call(SUM_ADD, state, values, weights);
}
SEXP life_calibration_sum_value(SEXP state, SEXP divisor) {
    return run_sum_call(SUM_VALUE, state, divisor, R_NilValue);
}
SEXP life_calibration_sum_reset(SEXP state) {
    return run_sum_call(SUM_RESET, state, R_NilValue, R_NilValue);
}
SEXP life_calibration_sum_close(SEXP state) {
    return run_sum_call(SUM_CLOSE, state, R_NilValue, R_NilValue);
}
SEXP life_calibration_sum_info(void) {
    return run_sum_call(SUM_INFO, R_NilValue, R_NilValue, R_NilValue);
}

/* The shared DLL hook checks every registry before it drains any of them. */
void life_calibration_sum_unload_check(void) {
    if (active_sum_calls)
        Rf_error("cannot unload lifeconformal during an extended-sum operation");
}
void life_calibration_sum_finalize_all(void) {
    while (states != NULL)
        R_RunWeakRefFinalizer(states->finalizer);
}
