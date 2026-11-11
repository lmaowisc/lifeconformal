#include <R.h>
#include <Rinternals.h>
#include <math.h>

/* The immediate successor representable double, not a user-scale epsilon.
   The R caller always verifies the structural certificate again. */
SEXP life_nextafter(SEXP x) {
    if (TYPEOF(x) != REALSXP)
        Rf_error("life_nextafter requires a numeric vector");
    R_xlen_t n = XLENGTH(x);
    SEXP answer = PROTECT(Rf_allocVector(REALSXP, n));
    const double *input = REAL(x);
    double *output = REAL(answer);
    for (R_xlen_t i = 0; i < n; ++i)
        output[i] = nextafter(input[i], R_PosInf);
    UNPROTECT(1);
    return answer;
}
