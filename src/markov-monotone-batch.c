#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include <math.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include "protected-total.h"

/* One record per requested through-time history. No H by L matrix exists. */
typedef struct { int index, subject, bin; } LifeQuery;

static SEXP field(SEXP x, const char *name) {
    SEXP names = Rf_getAttrib(x, R_NamesSymbol);
    if (TYPEOF(x) != VECSXP || TYPEOF(names) != STRSXP || XLENGTH(names)!=XLENGTH(x))
        Rf_error("Native monotone input must be a named list.");
    for (R_xlen_t i = 0; i < XLENGTH(x); ++i)
        if (strcmp(CHAR(STRING_ELT(names, i)), name) == 0) return VECTOR_ELT(x, i);
    Rf_error("Missing native monotone field '%s'.", name);
    return R_NilValue;
}

static void matrix_shape(SEXP x, int type, const char *name, int *rows, int *cols) {
    SEXP dim = Rf_getAttrib(x, R_DimSymbol);
    if (TYPEOF(x) != type || TYPEOF(dim) != INTSXP || XLENGTH(dim) != 2 ||
        INTEGER(dim)[0] < 0 || INTEGER(dim)[1] < 0)
        Rf_error("Invalid native monotone matrix '%s'.", name);
    *rows = INTEGER(dim)[0]; *cols = INTEGER(dim)[1];
    if (XLENGTH(x) != (R_xlen_t) *rows * *cols)
        Rf_error("Inconsistent native monotone dimensions for '%s'.", name);
}

static void require_shape(SEXP x, int type, const char *name, int nr, int nc) {
    int r, c; matrix_shape(x, type, name, &r, &c);
    if (r != nr || c != nc) Rf_error("Incorrect dimensions for '%s'.", name);
}

static int first_ge(const double *x, int n, double value) {
    int lo = 0, hi = n;
    while (lo < hi) { int mid = lo + (hi - lo) / 2;
        if (x[mid] < value) lo = mid + 1; else hi = mid; }
    return lo;
}

static int count_le(const double *x, int n, double value) {
    int lo = 0, hi = n;
    while (lo < hi) { int mid = lo + (hi - lo) / 2;
        if (x[mid] <= value) lo = mid + 1; else hi = mid; }
    return lo;
}

static int compare_queries(const void *aa, const void *bb) {
    const LifeQuery *a = (const LifeQuery *) aa, *b = (const LifeQuery *) bb;
    if (a->subject != b->subject) return a->subject < b->subject ? -1 : 1;
    if (a->bin != b->bin) return a->bin > b->bin ? -1 : 1;
    return a->index < b->index ? -1 : a->index > b->index;
}

SEXP life_markov_monotone_batch(SEXP spec) {
    SEXP jt = field(spec,"jump_times"), bi = field(spec,"base_increments");
    SEXP ed = field(spec,"edges"), rr = field(spec,"relative_risk");
    SEXP rw = field(spec,"rewards"), tt = field(spec,"tau");
    SEXP ce = field(spec,"center"), sc = field(spec,"scale"), qs = field(spec,"q");
    SEXP su = field(spec,"subject_index"), ti = field(spec,"query_times");
    SEXP st = field(spec,"query_state"), ac = field(spec,"accrued"), we = field(spec,"weight");
    SEXP dt = field(spec,"determination_time"), re = field(spec,"resolved");
    SEXP me = field(spec,"membership");
    SEXP fl = field(spec,"flip"), oc = field(spec,"original_center");
    SEXP oa = field(spec,"original_accrued");
    int n,e,k,d,jnum,cols,er,ec;
    matrix_shape(rr,REALSXP,"relative_risk",&n,&e);
    matrix_shape(rw,REALSXP,"rewards",&k,&d);
    matrix_shape(bi,REALSXP,"base_increments",&jnum,&cols);
    matrix_shape(ed,INTSXP,"edges",&er,&ec);
    if (n < 1 || k < 1 || d < 1 || cols != e || er != e || ec != 2)
        Rf_error("Incompatible native monotone model dimensions.");
    if (TYPEOF(jt)!=REALSXP || XLENGTH(jt)!=jnum || TYPEOF(qs)!=REALSXP ||
        XLENGTH(qs)<1 || XLENGTH(qs)>INT_MAX || TYPEOF(tt)!=REALSXP || XLENGTH(tt)!=1 ||
        !R_FINITE(REAL(tt)[0]) || REAL(tt)[0]<=0)
        Rf_error("Invalid native monotone time or threshold vectors.");
    int lnum=(int)XLENGTH(qs);
    double tau=REAL(tt)[0];
    require_shape(ce,REALSXP,"center",n,d); require_shape(sc,REALSXP,"scale",n,d);
    if (TYPEOF(ti)!=REALSXP || XLENGTH(ti)>INT_MAX) Rf_error("Invalid query times.");
    int h=(int)XLENGTH(ti);
    if (TYPEOF(su)!=INTSXP || TYPEOF(st)!=INTSXP || TYPEOF(we)!=REALSXP ||
        XLENGTH(su)!=h || XLENGTH(st)!=h || XLENGTH(we)!=h)
        Rf_error("Invalid native query vector lengths or types.");
    require_shape(ac,REALSXP,"accrued",h,d);
    require_shape(oc,REALSXP,"original_center",n,d);
    require_shape(oa,REALSXP,"original_accrued",h,d);
    if(TYPEOF(fl)!=LGLSXP || XLENGTH(fl)!=d) Rf_error("Invalid native orientation vector.");
    require_shape(dt,REALSXP,"determination_time",n,lnum);
    require_shape(re,LGLSXP,"resolved",n,lnum); require_shape(me,LGLSXP,"membership",n,lnum);
    const double *times=REAL(jt),*base=REAL(bi),*risk=REAL(rr),*reward=REAL(rw);
    const double *center=REAL(ce),*scale=REAL(sc),*q=REAL(qs),*query_time=REAL(ti);
    const double *accrued=REAL(ac),*weight=REAL(we),*determination=REAL(dt);
    const int *subject=INTEGER(su),*state=INTEGER(st),*edges=INTEGER(ed);
    const int *resolved=LOGICAL(re),*membership=LOGICAL(me);
    const int *flip=LOGICAL(fl);
    const double *original_center=REAL(oc),*original_accrued=REAL(oa);
    for(int j=0;j<jnum;++j) if(!R_FINITE(times[j]) || times[j]<=0 || times[j]>tau ||
        (j>0 && times[j]<=times[j-1])) Rf_error("Invalid native jump times.");
    for(int a=0;a<lnum;++a) if(!R_FINITE(q[a]) || q[a]<0 || (a>0 && q[a]<=q[a-1]))
        Rf_error("Invalid native thresholds.");
    for(R_xlen_t a=0;a<(R_xlen_t)n*d;++a) if(!R_FINITE(center[a]) ||
        !R_FINITE(original_center[a]) || !R_FINITE(scale[a]) || scale[a]<=0)
        Rf_error("Invalid native score parameters.");
    for(int b=0;b<d;++b) if(flip[b]==NA_LOGICAL) Rf_error("Missing native orientation.");
    for(R_xlen_t a=0;a<(R_xlen_t)k*d;++a) if(reward[a]!=0 && reward[a]!=1)
        Rf_error("Native monotone rewards must be binary.");
    for(int a=0;a<e;++a) {
        int f=edges[a]-1,t=edges[a+e]-1;
        if(f<0 || f>=k || t<0 || t>=k || f==t) Rf_error("Invalid native edge index.");
        for(int b=0;b<a;++b) if(edges[b]==edges[a] && edges[b+e]==edges[a+e])
            Rf_error("Duplicate native edge.");
        for(int b=0;b<d;++b) if(reward[t+(R_xlen_t)k*b]>reward[f+(R_xlen_t)k*b])
            Rf_error("Native monotone edge increases a target reward.");
    }
    for(R_xlen_t a=0;a<(R_xlen_t)jnum*e;++a) if(!R_FINITE(base[a]) || base[a]<0)
        Rf_error("Invalid native baseline increment.");
    for(R_xlen_t a=0;a<(R_xlen_t)n*e;++a) if(!R_FINITE(risk[a]) || risk[a]<=0)
        Rf_error("Invalid native relative risk.");
    for(int a=0;a<h;++a) {
        if(subject[a]<1 || subject[a]>n || state[a]<1 || state[a]>k ||
            !R_FINITE(query_time[a]) || query_time[a]<0 || query_time[a]>tau ||
            !R_FINITE(weight[a])) Rf_error("Invalid native query.");
        for(int b=0;b<d;++b) if(!R_FINITE(accrued[a+(R_xlen_t)h*b]) ||
            !R_FINITE(original_accrued[a+(R_xlen_t)h*b]) ||
            (flip[b] && reward[state[a]-1+(R_xlen_t)k*b]==1 &&
             original_accrued[a+(R_xlen_t)h*b]!=0) ||
            (reward[state[a]-1+(R_xlen_t)k*b]==1 && accrued[a+(R_xlen_t)h*b]!=query_time[a]))
            Rf_error("Active native monotone target must equal the query time.");
    }
    for(R_xlen_t a=0;a<(R_xlen_t)n*lnum;++a) {
        if(resolved[a]==NA_LOGICAL || (resolved[a] && (membership[a]==NA_LOGICAL ||
            !R_FINITE(determination[a]) || determination[a]<0 || determination[a]>tau)))
            Rf_error("Invalid native structural certificate.");
    }
    LifeQuery *queries=(LifeQuery *)R_alloc(h ? h : 1,sizeof(LifeQuery));
    int query_count=0;
    for(int a=0;a<h;++a) if(weight[a]!=0) {
        queries[query_count].index=a; queries[query_count].subject=subject[a]-1;
        queries[query_count].bin=count_le(times,jnum,query_time[a]); ++query_count;
    }
    qsort(queries,query_count,sizeof(LifeQuery),compare_queries);
    double *buffer_a=(double *)R_alloc((size_t)k*lnum,sizeof(double));
    double *buffer_b=(double *)R_alloc((size_t)k*lnum,sizeof(double));
    double **value=(double **)R_alloc(k,sizeof(double *));
    double **next=(double **)R_alloc(k,sizeof(double *));
    for(int s=0;s<k;++s) {value[s]=buffer_a+(size_t)s*lnum; next[s]=buffer_b+(size_t)s*lnum;}
    double *p=(double *)R_alloc(e ? e : 1,sizeof(double));
    double *total=(double *)R_alloc(k,sizeof(double));
    /* R stores the returned subject-by-threshold matrix by columns. Keep
       the repeatedly updated row and its certificate metadata contiguous
       while traversing one subject, then scatter that row only once. */
    double *subject_sum=(double *)R_alloc(lnum,sizeof(double));
    double *compensation=(double *)R_alloc(lnum,sizeof(double));
    double *subject_determination=(double *)R_alloc(lnum,sizeof(double));
    int *subject_resolved=(int *)R_alloc(lnum,sizeof(int));
    int *subject_membership=(int *)R_alloc(lnum,sizeof(int));
    SEXP out=PROTECT(Rf_allocMatrix(REALSXP,n,lnum));
    memset(REAL(out),0,(size_t)n*lnum*sizeof(double));
    double *answer=REAL(out);
    R_xlen_t adjusted=0; int cursor=0;
    for(int i=0;i<n;++i) {
        if((i&15)==0) R_CheckUserInterrupt();
        if(cursor>=query_count || queries[cursor].subject!=i) continue;
        int end=cursor; while(end<query_count && queries[end].subject==i) ++end;
        memset(subject_sum,0,(size_t)lnum*sizeof(double));
        memset(compensation,0,(size_t)lnum*sizeof(double));
        for(int a=0;a<lnum;++a) {
            R_xlen_t pos=i+(R_xlen_t)n*a;
            subject_determination[a]=determination[pos];
            subject_resolved[a]=resolved[pos];
            subject_membership[a]=membership[pos];
        }
        double full_score=0;
        for(int b=0;b<d;++b) {
            double lower=fabs(original_center[i+(R_xlen_t)n*b]/scale[i+(R_xlen_t)n*b]);
            double upper=fabs((tau-original_center[i+(R_xlen_t)n*b])/scale[i+(R_xlen_t)n*b]);
            if(lower>full_score) full_score=lower;
            if(upper>full_score) full_score=upper;
        }
        /* Above this threshold every exit clock in [0,tau] is covered.
           Those state values are identically one and need no recursion. */
        int active_lnum=first_ge(q,lnum,full_score);
        for(int s=0;s<k;++s) {
            double threshold=0;
            for(int b=0;b<d;++b) if(reward[s+(R_xlen_t)k*b]==1) {
                double terminal=flip[b] ? 0 : tau;
                double v=fabs((terminal-original_center[i+(R_xlen_t)n*b])/scale[i+(R_xlen_t)n*b]);
                if(v>threshold) threshold=v; }
            int cutoff=first_ge(q,lnum,threshold);
            for(int a=0;a<lnum;++a) {value[s][a]=(a>=cutoff); next[s][a]=1;}
        }
        int bin=jnum;
        while(cursor<end) {
            while(cursor<end && queries[cursor].bin==bin) {
                int z=queries[cursor].index,s=state[z]-1;
                double fixed=0;
                for(int b=0;b<d;++b) if(reward[s+(R_xlen_t)k*b]==0) {
                    double target=original_accrued[z+(R_xlen_t)h*b];
                    if(flip[b]) target=life_protected_total_reward(tau,query_time[z],1,target);
                    if(!R_FINITE(target)) Rf_error("Protected monotone target overflowed; check reward and time units.");
                    double v=fabs((target-original_center[i+(R_xlen_t)n*b])/scale[i+(R_xlen_t)n*b]);
                    if(v>fixed) fixed=v;
                }
                int cutoff=first_ge(q,lnum,fixed);
                for(int a=0;a<lnum;++a) {
                    double prob=a>=cutoff ? value[s][a] : 0;
                    if(subject_resolved[a] && query_time[z]>=subject_determination[a])
                        prob=subject_membership[a];
                    if(!R_FINITE(prob) || prob < -1e-10 || prob > 1+1e-10)
                        Rf_error("Native monotone probability is outside [0,1].");
                    if(prob<0) prob=0; else if(prob>1) prob=1;
                    double term=weight[z]*prob, previous=subject_sum[a], sum=previous+term;
                    if(fabs(previous)>=fabs(term)) compensation[a]+=(previous-sum)+term;
                    else compensation[a]+=(term-sum)+previous;
                    subject_sum[a]=sum;
                }
                ++cursor;
            }
            if(cursor>=end) break;
            int j=--bin;
            if(j<0) Rf_error("Internal native monotone query ordering failed.");
            if((j&255)==0) R_CheckUserInterrupt();
            memset(total,0,(size_t)k*sizeof(double));
            for(int a=0;a<e;++a) {
                p[a]=base[j+(R_xlen_t)jnum*a]*risk[i+(R_xlen_t)n*a];
                if(p[a]==0 && base[j+(R_xlen_t)jnum*a]>0)
                    Rf_error("Positive Markov intensity increments underflow; check the fitted hazards and prediction profiles.");
                total[edges[a]-1]+=p[a];
            }
            for(int s=0;s<k;++s) {
                if(!R_FINITE(total[s])) Rf_error("Native Markov outgoing increments overflow.");
                if(total[s]>1) ++adjusted;
                if(total[s]>0) {
                    double stay=total[s]>1 ? 0 : 1-total[s];
                    for(int a=0;a<active_lnum;++a) next[s][a]=stay*value[s][a];
                }
            }
            for(int a=0;a<e;++a) if(p[a]>0) {
                int f=edges[a]-1,t=edges[a+e]-1;
                double probability=p[a]/(total[f]>1 ? total[f] : 1), threshold=0;
                if(probability==0)
                    Rf_error("Positive Markov transition probabilities underflow during row normalization; check the fitted intensity increments.");
                for(int b=0;b<d;++b) if(reward[f+(R_xlen_t)k*b]>reward[t+(R_xlen_t)k*b]) {
                    double target=flip[b] ? tau-times[j] : times[j];
                    double v=fabs((target-original_center[i+(R_xlen_t)n*b])/scale[i+(R_xlen_t)n*b]);
                    if(v>threshold) threshold=v;
                }
                int cutoff=first_ge(q,lnum,threshold);
                for(int b=cutoff;b<active_lnum;++b)
                    next[f][b]+=probability*value[t][b];
            }
            /* Swap only changed origins after every simultaneous update. */
            for(int s=0;s<k;++s) if(total[s]>0) {
                double *swap=value[s]; value[s]=next[s]; next[s]=swap;
            }
        }
        for(int a=0;a<lnum;++a)
            answer[i+(R_xlen_t)n*a]=subject_sum[a]+compensation[a];
    }
    for(R_xlen_t a=0;a<(R_xlen_t)n*lnum;++a) if(!R_FINITE(answer[a]))
        Rf_error("Native signed coverage contributions overflow.");
    SEXP result=PROTECT(Rf_allocVector(VECSXP,2));
    SEXP names=PROTECT(Rf_allocVector(STRSXP,2));
    SET_STRING_ELT(names,0,Rf_mkChar("contributions")); SET_STRING_ELT(names,1,Rf_mkChar("adjusted_rows"));
    SET_VECTOR_ELT(result,0,out); SET_VECTOR_ELT(result,1,Rf_ScalarReal((double)adjusted));
    Rf_setAttrib(result,R_NamesSymbol,names);
    UNPROTECT(3); return result;
}
