#ifndef LIFECONFORMAL_PROTECTED_TOTAL_H
#define LIFECONFORMAL_PROTECTED_TOTAL_H
#include <float.h>
#include <math.h>

/* Nonoverlapping expansion for tau - time + accrued. FastTwoSum is exact
 * after the magnitude swap. The final residual-sign correction handles a
 * halfway rounding decision influenced by a still smaller partial.
 * The input order avoids spurious overflow for valid occupation histories.
 * Callers must reject nonfinite results, as for their other reward arithmetic.
 */
static inline double life_occupation_sum3(double tau, double time, double accrued) {
    double partials[3], inputs[3] = {tau, -time, accrued};
    int n = 0;
    for (int index = 0; index < 3; ++index) {
        double x = inputs[index];
        int i = 0;
        for (int j = 0; j < n; ++j) {
            double y = partials[j];
            if (fabs(x) < fabs(y)) {
                const double temporary = x; x = y; y = temporary;
            }
            const double hi = x + y;
            if (hi > DBL_MAX || hi < -DBL_MAX) return hi;
            const double lo = y - (hi - x);
            if (lo != 0) partials[i++] = lo;
            x = hi;
        }
        n = i;
        if (x != 0) partials[n++] = x;
    }
    if (!n) return 0;
    double hi = partials[--n], lo = 0;
    while (n) {
        const double x = hi, y = partials[--n];
        hi = x + y;
        lo = y - (hi - x);
        if (lo != 0) break;
    }
    if (n && ((lo < 0 && partials[n-1] < 0) ||
              (lo > 0 && partials[n-1] > 0))) {
        const double y = lo * 2, x = hi + y;
        if (y == x - hi) hi = x;
    }
    return hi;
}

/* Binary occupation rewards are the public contract. Arbitrary internal
 * numeric rewards retain the original formula without an exact-product claim.
 */
static inline double life_protected_total_reward(double tau, double time,
                                                double reward, double accrued) {
    if (reward == 0) return accrued;
    if (reward == 1) return life_occupation_sum3(tau, time, accrued);
    return accrued + (tau-time)*reward;
}
#endif
