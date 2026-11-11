#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include "protected-total.h"

/* Shared-threshold sweep. State and target matrices follow R column-major
   order; subject row-index lists preserve interleaved input without sorting.
   Bounds depend only on accrued history and declared structural reachability.
   The R implementation in cdt.R is retained as an independent execution path. */

typedef struct {
    int n, d, K, L;
    const double *reward, *minimum, *maximum, *center, *scale, *q;
    double tau;
    double *time, *certificate;
    int *member, *later;
    double segments, evaluations, roots, adjustments, adjacent, refinements;
    double fallbacks, delay, complete;
} Cdt;

typedef struct {
    double *lower, *upper;
    double low, high;
    int open;
} Bounds;

static SEXP named(SEXP x, const char *key) {
    if (TYPEOF(x) != VECSXP) Rf_error("Invalid CDT input list");
    SEXP names = Rf_getAttrib(x, R_NamesSymbol);
    for (R_xlen_t j = 0; j < XLENGTH(x); ++j)
        if (names != R_NilValue && !strcmp(CHAR(STRING_ELT(names, j)), key))
            return VECTOR_ELT(x, j);
    Rf_error("Missing CDT input field: %s", key);
    return R_NilValue;
}

static int threshold_count(const double *q, int L, double value, int inclusive) {
    int left = 0, right = L;
    while (left < right) {
        int middle = left + (right - left) / 2;
        if (q[middle] < value || (inclusive && q[middle] == value)) left = middle + 1;
        else right = middle;
    }
    return left;
}

static void bounds(Cdt *c, int person, int state, const double *accrued,
                   double t, int use_open, Bounds *b) {
    double remain = fmax(0, c->tau - t);
    b->low = 0;
    b->high = 0;
    b->open = 0;
    for (int j = 0; j < c->d; ++j) {
        R_xlen_t sk = state + (R_xlen_t)c->K * j;
        R_xlen_t pj = person + (R_xlen_t)c->n * j;
        b->lower[j] = life_protected_total_reward(c->tau,t,c->minimum[sk],accrued[j]);
        b->upper[j] = life_protected_total_reward(c->tau,t,c->maximum[sk],accrued[j]);
        if (!R_FINITE(b->lower[j]) || !R_FINITE(b->upper[j]))
            Rf_error("CDT occupation bounds overflowed.");
        double dl = (b->lower[j] - c->center[pj]) / c->scale[pj];
        double du = (c->center[pj] - b->upper[j]) / c->scale[pj];
        int lop = use_open && remain > 0 && c->reward[sk] != c->minimum[sk];
        int uop = use_open && remain > 0 && c->reward[sk] != c->maximum[sk];
        if (dl > b->low) { b->low = dl; b->open = lop; }
        else if (dl == b->low && lop) b->open = 1;
        if (du > b->low) { b->low = du; b->open = uop; }
        else if (du == b->low && uop) b->open = 1;
        b->high = fmax(b->high, fmax(fabs(dl), fabs(du)));
    }
    c->evaluations += 1;
}

static void record(Cdt *c, int person, int k, double t, int value, int later) {
    R_xlen_t at = person + (R_xlen_t)c->n * k;
    c->time[at] = t;
    c->certificate[at] = t;
    c->member[at] = value;
    c->later[at] = later;
}

static void classify(Cdt *c, int person, Bounds *b, double t, int *lo, int *hi) {
    int last = threshold_count(c->q, c->L, b->low, b->open) - 1;
    if (last > *hi) last = *hi;
    while (*lo <= last) { record(c, person, *lo, t, 0, 0); *lo += 1; }
    int first = threshold_count(c->q, c->L, b->high, 0);
    if (first < *lo) first = *lo;
    while (*hi >= first) { record(c, person, *hi, t, 1, 0); *hi -= 1; }
}

static int verify(Cdt *c, int person, int state, Bounds *b, const double *accrued, double start,
                  double end, int event, int k, double t, int inside) {
    c->evaluations += 1;
    if (!R_FINITE(t) || t < start || t > end) return 0;
    int can_open = t < c->tau && !(event && t == end);
    double q = c->q[k];
    for (int j = 0; j < c->d; ++j) {
        R_xlen_t sk = state + (R_xlen_t)c->K * j;
        R_xlen_t pj = person + (R_xlen_t)c->n * j;
        double rise = c->reward[sk] - c->minimum[sk];
        double fall = c->maximum[sk] - c->reward[sk];
        double lower = rise == 0 ? b->lower[j] : life_occupation_sum3(t,start,accrued[j]);
        double upper = fall == 0 ? b->upper[j] : life_occupation_sum3(c->tau,t,accrued[j]);
        if (!R_FINITE(lower) || !R_FINITE(upper))
            Rf_error("CDT occupation bounds overflowed.");
        double dl = (lower - c->center[pj]) / c->scale[pj];
        double du = (c->center[pj] - upper) / c->scale[pj];
        if (inside) {
            if (fmax(fabs(dl), fabs(du)) > q) return 0;
        } else if (dl > q || du > q ||
                   (dl == q && can_open && rise > 0) ||
                   (du == q && can_open && fall > 0)) return 1;
    }
    return inside;
}

static void solve(Cdt *c, int person, int state, Bounds *b, const double *accrued, double start,
                  double end, int event, int k, int inside) {
    double elapsed = inside ? 0 : R_PosInf;
    for (int j = 0; j < c->d; ++j) {
        R_xlen_t sk = state + (R_xlen_t)c->K * j;
        R_xlen_t pj = person + (R_xlen_t)c->n * j;
        double rise = c->reward[sk] - c->minimum[sk];
        double fall = c->maximum[sk] - c->reward[sk];
        double width = c->q[k] * c->scale[pj];
        if (inside) {
            double required_low = c->center[pj] - width - b->lower[j];
            double required_high = b->upper[j] - c->center[pj] - width;
            /* A stationary endpoint has no root to solve. Test its signed
               standardized distance just as in bounds() and verify(). At
               an inclusive boundary, multiplying q by scale can round back
               to a smaller distance and incorrectly return infinity. */
            double tl = rise > 0 ? required_low / rise :
                ((c->center[pj] - b->lower[j]) / c->scale[pj] <= c->q[k] ? 0 : R_PosInf);
            double tu = fall > 0 ? required_high / fall :
                ((b->upper[j] - c->center[pj]) / c->scale[pj] <= c->q[k] ? 0 : R_PosInf);
            elapsed = fmax(elapsed, fmax(tl, tu));
        } else {
            double tl = rise > 0 ? (c->center[pj] + width - b->lower[j]) / rise : R_PosInf;
            double tu = fall > 0 ? (b->upper[j] - c->center[pj] + width) / fall : R_PosInf;
            elapsed = fmin(elapsed, fmin(tl, tu));
        }
    }
    c->roots += 1;
    if (((R_xlen_t)c->roots & 65535) == 0) R_CheckUserInterrupt();
    double initial = start + fmax(0, elapsed);
    double root = initial;
    int ok = verify(c, person, state, b, accrued, start, end, event, k, root, inside);
    for (int step = 0; step < 8 && !ok && R_FINITE(root) && root >= start && root < end; ++step) {
        root = nextafter(root, R_PosInf);
        ok = verify(c, person, state, b, accrued, start, end, event, k, root, inside);
    }
    if (ok && root != initial) c->adjacent += 1;
    if (!ok && R_FINITE(root) && root >= start && root < end &&
        verify(c, person, state, b, accrued, start, end, event, k, end, inside)) {
        double lower = root, upper = end;
        for (int iteration = 0; iteration < 64; ++iteration) {
            double middle = lower + (upper - lower) / 2;
            if (!(middle > lower && middle < upper)) break;
            if (verify(c, person, state, b, accrued, start, end, event, k, middle, inside)) upper = middle;
            else lower = middle;
        }
        root = upper;
        ok = 1;
        c->refinements += 1;
    }
    if (!ok) { root = end; c->fallbacks += 1; }
    int adjusted = root != initial || !ok;
    if (adjusted) {
        c->adjustments += 1;
        if (R_FINITE(initial)) c->delay = fmax(c->delay, root - initial);
    }
    record(c, person, k, root, inside, adjusted);
}

SEXP life_cdt_sweep(SEXP prepared, SEXP thresholds) {
    int np = 0;
    SEXP history = named(prepared, "history");
    SEXP problem = named(prepared, "problem");
    SEXP subjects = named(history, "subjects");
    SEXP intervals = named(history, "intervals");
    SEXP support = named(problem, "structural_support");
    SEXP rows = named(history, "row_index");
    SEXP rewards = named(problem, "rewards");
    SEXP dim = Rf_getAttrib(rewards, R_DimSymbol);
    if (TYPEOF(dim) != INTSXP || XLENGTH(dim) != 2) Rf_error("Invalid CDT reward matrix");
    int K = INTEGER(dim)[0], d = INTEGER(dim)[1];
    if (K <= 0 || d <= 0 || TYPEOF(rows) != VECSXP || XLENGTH(rows) > INT_MAX ||
        XLENGTH(thresholds) <= 0 || XLENGTH(thresholds) > INT_MAX)
        Rf_error("Invalid CDT dimensions");
    int n = (int)XLENGTH(rows), L = (int)XLENGTH(thresholds);
#define REAL_INPUT(v, x) SEXP v = PROTECT(Rf_coerceVector((x), REALSXP)); ++np
#define INT_INPUT(v, x) SEXP v = PROTECT(Rf_coerceVector((x), INTSXP)); ++np
    REAL_INPUT(rr, rewards);
    REAL_INPUT(rmin, named(support, "reward_min"));
    REAL_INPUT(rmax, named(support, "reward_max"));
    REAL_INPUT(mu, named(prepared, "center"));
    REAL_INPUT(sc, named(prepared, "scale"));
    REAL_INPUT(q, thresholds);
    REAL_INPUT(observed, named(subjects, "observation_end"));
    REAL_INPUT(stops, named(intervals, "stop"));
    INT_INPUT(initial, named(subjects, "initial_state"));
    INT_INPUT(events, named(intervals, "event"));
    INT_INPUT(to, named(intervals, "to"));
    INT_INPUT(censored, named(intervals, "terminal_censor"));
    INT_INPUT(absorbing, named(problem, "absorbing"));
#undef REAL_INPUT
#undef INT_INPUT
    R_xlen_t m = XLENGTH(stops);
    if (XLENGTH(rmin) != (R_xlen_t)K*d || XLENGTH(rmax) != (R_xlen_t)K*d ||
        XLENGTH(mu) != (R_xlen_t)n*d || XLENGTH(sc) != (R_xlen_t)n*d ||
        XLENGTH(initial) != n || XLENGTH(observed) != n || XLENGTH(events) != m ||
        XLENGTH(to) != m || XLENGTH(censored) != m)
        Rf_error("Inconsistent CDT array lengths");
    int *is_absorbing = (int *)R_alloc(K, sizeof(int));
    for (int k = 0; k < K; ++k) is_absorbing[k] = 0;
    for (R_xlen_t j = 0; j < XLENGTH(absorbing); ++j) {
        int k = INTEGER(absorbing)[j];
        if (k < 1 || k > K) Rf_error("Invalid CDT absorbing index");
        is_absorbing[k-1] = 1;
    }
    SEXP times = PROTECT(Rf_allocMatrix(REALSXP, n, L)); ++np;
    SEXP members = PROTECT(Rf_allocMatrix(LGLSXP, n, L)); ++np;
    SEXP certificates = PROTECT(Rf_allocMatrix(REALSXP, n, L)); ++np;
    SEXP later = PROTECT(Rf_allocMatrix(LGLSXP, n, L)); ++np;
    SEXP resolved = PROTECT(Rf_allocMatrix(LGLSXP, n, L)); ++np;
    SEXP endpoint = PROTECT(Rf_allocVector(REALSXP, n)); ++np;
    Cdt c = {0};
    c.n=n; c.d=d; c.K=K; c.L=L; c.reward=REAL(rr); c.minimum=REAL(rmin);
    c.maximum=REAL(rmax); c.center=REAL(mu); c.scale=REAL(sc); c.q=REAL(q);
    c.tau=Rf_asReal(named(problem,"tau")); c.time=REAL(times);
    c.certificate=REAL(certificates); c.member=LOGICAL(members); c.later=LOGICAL(later);
    if (!R_FINITE(c.tau) || c.tau <= 0) Rf_error("Invalid CDT horizon");
    for (int k=0; k<L; ++k)
        if (ISNAN(c.q[k]) || c.q[k]<0 || (k>0 && c.q[k]<=c.q[k-1]))
            Rf_error("Invalid CDT threshold order");
    for (R_xlen_t j=0; j<(R_xlen_t)n*d; ++j)
        if (!R_FINITE(c.center[j]) || !R_FINITE(c.scale[j]) || c.scale[j]<=0)
            Rf_error("Invalid CDT score parameters");
    for (int i=0; i<n; ++i) {
        if (!R_FINITE(REAL(observed)[i]) || REAL(observed)[i] < 0)
            Rf_error("Invalid CDT observation endpoint");
        REAL(endpoint)[i]=fmin(REAL(observed)[i], c.tau);
    }
    for (int k=0; k<L; ++k) for (int i=0; i<n; ++i) {
        R_xlen_t at=i+(R_xlen_t)n*k;
        c.time[at]=REAL(endpoint)[i]; c.member[at]=NA_LOGICAL;
        c.certificate[at]=NA_REAL; c.later[at]=0;
    }
    double *accrued=(double*)R_alloc(d,sizeof(double));
    double *end_accrued=(double*)R_alloc(d,sizeof(double));
    int *constant_reward=(int*)R_alloc(d,sizeof(int));
    Bounds start_b={0}, end_b={0};
    start_b.lower=(double*)R_alloc(d,sizeof(double)); start_b.upper=(double*)R_alloc(d,sizeof(double));
    end_b.lower=(double*)R_alloc(d,sizeof(double)); end_b.upper=(double*)R_alloc(d,sizeof(double));
    for (int i=0; i<n; ++i) {
        if (i % 128 == 0) R_CheckUserInterrupt();
        int initial_index=INTEGER(initial)[i];
        if (initial_index<1 || initial_index>K) Rf_error("Invalid CDT initial-state index");
        int state=initial_index-1;
        for (int j=0; j<d; ++j) { accrued[j]=0; constant_reward[j]=1; }
        double t=0;
        int lo=0, hi=L-1;
        bounds(&c,i,state,accrued,t,1,&start_b);
        classify(&c,i,&start_b,t,&lo,&hi);
        SEXP indices=VECTOR_ELT(rows,i);
        if (TYPEOF(indices)!=INTSXP) Rf_error("Invalid CDT history row indices");
        for (R_xlen_t row=0; row<XLENGTH(indices) && lo<=hi && t<c.tau; ++row) {
            int row_index=INTEGER(indices)[row];
            if (row_index<1 || (R_xlen_t)row_index>m)
                Rf_error("CDT history row index out of bounds");
            int r=row_index-1;
            double stop=REAL(stops)[r], end=fmin(stop,c.tau);
            if (!R_FINITE(stop) || stop<t) Rf_error("Invalid CDT chronological endpoint");
            if ((INTEGER(events)[r]!=0 && INTEGER(events)[r]!=1) ||
                (INTEGER(censored)[r]!=0 && INTEGER(censored)[r]!=1))
                Rf_error("Invalid CDT event indicator");
            int event=INTEGER(events)[r] && stop<=c.tau;
            double duration=end-t;
            for (int j=0; j<d; ++j) {
                double reward=c.reward[state+(R_xlen_t)K*j];
                double initial_reward=c.reward[initial_index-1+(R_xlen_t)K*j];
                if (duration>0 && reward!=initial_reward) constant_reward[j]=0;
                end_accrued[j]=constant_reward[j] ? end*initial_reward : accrued[j]+duration*reward;
                if (!R_FINITE(end_accrued[j])) Rf_error("CDT accrued occupation overflowed.");
            }
            bounds(&c,i,state,end_accrued,end,!event,&end_b);
            c.segments+=1;
            int last=threshold_count(c.q,L,end_b.low,end_b.open)-1;
            if (last>hi) last=hi;
            int first=threshold_count(c.q,L,end_b.high,0);
            if (first<lo) first=lo;
            while (lo<=last) { solve(&c,i,state,&start_b,accrued,t,end,event,lo,0); ++lo; }
            if (first<lo) first=lo;
            while (hi>=first) { solve(&c,i,state,&start_b,accrued,t,end,event,hi,1); --hi; }
            for (int j=0; j<d; ++j) accrued[j]=end_accrued[j];
            t=end;
            if (event) {
                int destination=INTEGER(to)[r];
                if (destination<1 || destination>K) Rf_error("Invalid CDT destination index");
                state=destination-1;
            }
            bounds(&c,i,state,accrued,t,1,&start_b);
            if (lo<=hi) classify(&c,i,&start_b,t,&lo,&hi);
            if (t>=c.tau || is_absorbing[state] || INTEGER(censored)[r]) break;
        }
        if (lo>hi) c.complete+=1;
    }
    for (R_xlen_t j=0; j<(R_xlen_t)n*L; ++j) LOGICAL(resolved)[j]=c.member[j]!=NA_LOGICAL;
    const char *dn[]={"segments","bound_evaluations","analytic_roots","numerical_later_certificates",
                      "adjacent_float_adjustments","bracket_refinements","endpoint_fallbacks",
                      "max_root_delay","complete_subjects"};
    double dv[]={c.segments,c.evaluations,c.roots,c.adjustments,c.adjacent,c.refinements,
                 c.fallbacks,c.delay,c.complete};
    SEXP diagnostics=PROTECT(Rf_allocVector(VECSXP,9)); ++np;
    SEXP diagnostic_names=PROTECT(Rf_allocVector(STRSXP,9)); ++np;
    for (int j=0; j<9; ++j) {
        SET_VECTOR_ELT(diagnostics,j,Rf_ScalarReal(dv[j]));
        SET_STRING_ELT(diagnostic_names,j,Rf_mkChar(dn[j]));
    }
    Rf_setAttrib(diagnostics,R_NamesSymbol,diagnostic_names);
    SEXP result=PROTECT(Rf_allocVector(VECSXP,7)); ++np;
    SEXP result_names=PROTECT(Rf_allocVector(STRSXP,7)); ++np;
    const char *rn[]={"time","resolved","membership","certificate_time",
                      "numerical_later_certificate","endpoint","diagnostics"};
    SEXP rv[]={times,resolved,members,certificates,later,endpoint,diagnostics};
    for (int j=0; j<7; ++j) { SET_VECTOR_ELT(result,j,rv[j]); SET_STRING_ELT(result_names,j,Rf_mkChar(rn[j])); }
    Rf_setAttrib(result,R_NamesSymbol,result_names);
    UNPROTECT(np);
    return result;
}
