#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include "protected-total.h"

/* Complete aligned query occupations without first rounding end-start.
 * Both matrices are query by target, in R column-major order. */
SEXP life_occupation_complete(SEXP end, SEXP start, SEXP reward, SEXP accrued) {
    if (TYPEOF(end) != REALSXP || TYPEOF(start) != REALSXP ||
        TYPEOF(reward) != REALSXP || TYPEOF(accrued) != REALSXP)
        Rf_error("Occupation completion inputs must be double arrays.");
    SEXP dim = Rf_getAttrib(accrued, R_DimSymbol);
    SEXP rd = Rf_getAttrib(reward, R_DimSymbol);
    if (TYPEOF(dim) != INTSXP || XLENGTH(dim) != 2 ||
        TYPEOF(rd) != INTSXP || XLENGTH(rd) != 2 ||
        INTEGER(dim)[0] < 0 || INTEGER(dim)[1] < 0 ||
        INTEGER(dim)[0] != INTEGER(rd)[0] || INTEGER(dim)[1] != INTEGER(rd)[1])
        Rf_error("Occupation reward and accrued matrices must have matching dimensions.");
    R_xlen_t n = INTEGER(dim)[0], d = INTEGER(dim)[1];
    if (XLENGTH(accrued) != n*d || XLENGTH(reward) != n*d ||
        XLENGTH(end) != n || XLENGTH(start) != n ||
        Rf_getAttrib(end, R_DimSymbol) != R_NilValue ||
        Rf_getAttrib(start, R_DimSymbol) != R_NilValue)
        Rf_error("Occupation times must have one entry per matrix row.");
    for (R_xlen_t i = 0; i < n; ++i)
        if (!R_FINITE(REAL(start)[i]) || !R_FINITE(REAL(end)[i]) ||
            REAL(start)[i] < 0 || REAL(end)[i] < REAL(start)[i])
            Rf_error("Occupation times must be finite with 0 <= start <= end.");
    SEXP result = PROTECT(Rf_duplicate(accrued));
    for (R_xlen_t j = 0; j < d; ++j) for (R_xlen_t i = 0; i < n; ++i) {
        R_xlen_t k = i + n*j;
        double r = REAL(reward)[k], a = REAL(accrued)[k];
        if (!R_FINITE(r) || (r != 0 && r != 1) || !R_FINITE(a) || a < 0)
            Rf_error("Occupation rewards must be binary and accrued values finite and nonnegative.");
        double value = life_protected_total_reward(REAL(end)[i], REAL(start)[i], r, a);
        if (!R_FINITE(value)) Rf_error("Occupation completion overflowed.");
        REAL(result)[k] = value;
        if ((k & 65535) == 0) R_CheckUserInterrupt();
    }
    UNPROTECT(1);
    return result;
}
