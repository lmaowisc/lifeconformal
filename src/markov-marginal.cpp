// Same binary64 products, Neumaier updates and query order as the grid stream.
// Call-local R-managed storage; no persistent ownership or callback state.
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Arith.h>
#include <cmath>
#include <cstring>
#include <initializer_list>

namespace {
SEXP item(SEXP x, const char* name) {
    SEXP names = Rf_getAttrib(x, R_NamesSymbol);
    if (TYPEOF(names)!=STRSXP || XLENGTH(names)!=XLENGTH(x))
        Rf_error("Marginal reduction requires complete field names.");
    for (R_xlen_t i=0; i<XLENGTH(x); ++i)
        if (TYPEOF(names)==STRSXP && !std::strcmp(CHAR(STRING_ELT(names,i)),name))
            return VECTOR_ELT(x,i);
    Rf_error("Missing marginal reducer field: %s",name);
    return R_NilValue;
}
void compensated_add(double value, double& sum, double& correction) {
    const double updated = sum + value;
    if (std::fabs(sum) >= std::fabs(value)) correction += (sum - updated) + value;
    else correction += (value - updated) + sum;
    sum = updated;
}
}
extern "C" SEXP life_markov_marginal_reduce(SEXP specification) {
    if (TYPEOF(specification)!=VECSXP) Rf_error("Marginal reduction requires a named list.");
    SEXP p=item(specification,"probability"), lo=item(specification,"lower"), hi=item(specification,"upper");
    SEXP weights=item(specification,"weights"), order=item(specification,"order");
    SEXP dims=Rf_getAttrib(p,R_DimSymbol);
    if (TYPEOF(p)!=REALSXP || TYPEOF(dims)!=INTSXP || XLENGTH(dims)!=2)
        Rf_error("Marginal probability must be a double matrix.");
    const int h=INTEGER(dims)[0], length=INTEGER(dims)[1];
    if (h<0 || length<0 || (length>0 && static_cast<R_xlen_t>(h)>R_XLEN_T_MAX/length) ||
        XLENGTH(p)!=static_cast<R_xlen_t>(h)*length)
        Rf_error("Invalid marginal probability dimensions.");
    // Match the marginal R route's cell cap before any new native allocation.
    // Individual-dimension caps also protect empty matrices with huge margins.
    if (h>1000000 || length>1000000 || XLENGTH(p)>1000000)
        Rf_error("Marginal reduction exceeds its one-million-cell budget.");
    for (SEXP x : {lo,hi}) {
        SEXP d=Rf_getAttrib(x,R_DimSymbol);
        if (TYPEOF(x)!=REALSXP || XLENGTH(x)!=XLENGTH(p) || TYPEOF(d)!=INTSXP ||
            XLENGTH(d)!=2 || INTEGER(d)[0]!=h || INTEGER(d)[1]!=length)
            Rf_error("Marginal bounds must match the probability matrix.");
    }
    if (TYPEOF(weights)!=REALSXP || XLENGTH(weights)!=h || TYPEOF(order)!=INTSXP || XLENGTH(order)!=h)
        Rf_error("Marginal weights/order must have one entry per query.");
    unsigned char* seen=reinterpret_cast<unsigned char*>(R_alloc(h,sizeof(unsigned char)));
    if(h) std::memset(seen,0,h);
    for (int a=0; a<h; ++a) {
        if ((a & 65535)==0) R_CheckUserInterrupt();
        const int index=INTEGER(order)[a];
        if (index==NA_INTEGER || index<1 || index>h || seen[index-1])
            Rf_error("Marginal query order must be a permutation.");
        seen[index-1]=1;
        if (!R_FINITE(REAL(weights)[a])) Rf_error("Marginal weights must be finite.");
    }
    for (R_xlen_t i=0; i<XLENGTH(p); ++i) {
        if ((i & 65535)==0) R_CheckUserInterrupt();
        if (!R_FINITE(REAL(p)[i]) || !R_FINITE(REAL(lo)[i]) || !R_FINITE(REAL(hi)[i]) ||
            REAL(lo)[i]<0 || REAL(hi)[i]>1 || REAL(lo)[i]>REAL(p)[i] || REAL(p)[i]>REAL(hi)[i])
            Rf_error("Marginal probability/bounds must be finite ordered values in [0,1].");
    }
    SEXP result=PROTECT(Rf_allocVector(VECSXP,3));
    SEXP names=PROTECT(Rf_allocVector(STRSXP,3));
    const char* labels[] = {"probability","lower","upper"};
    for(int field=0;field<3;++field) {
        SET_STRING_ELT(names,field,Rf_mkChar(labels[field]));
        SEXP values=PROTECT(Rf_allocVector(REALSXP,length));
        SET_VECTOR_ELT(result,field,values); UNPROTECT(1);
    }
    Rf_setAttrib(result,R_NamesSymbol,names);
    double* correction=reinterpret_cast<double*>(R_alloc(static_cast<size_t>(length)*3,sizeof(double)));
    for(int field=0;field<3;++field)
        for(int b=0;b<length;++b) { REAL(VECTOR_ELT(result,field))[b]=0; correction[field*static_cast<size_t>(length)+b]=0; }
    for(int position=0;position<h;++position) {
        if ((position & 1023)==0) R_CheckUserInterrupt();
        const int a=INTEGER(order)[position]-1;
        const double w=REAL(weights)[a];
        for(int b=0;b<length;++b) {
            const R_xlen_t index=a+static_cast<R_xlen_t>(h)*b;
            compensated_add(w*REAL(p)[index],REAL(VECTOR_ELT(result,0))[b],correction[b]);
            compensated_add(w*(w>=0 ? REAL(lo)[index] : REAL(hi)[index]),
                REAL(VECTOR_ELT(result,1))[b],correction[static_cast<size_t>(length)+b]);
            compensated_add(w*(w>=0 ? REAL(hi)[index] : REAL(lo)[index]),
                REAL(VECTOR_ELT(result,2))[b],correction[2*static_cast<size_t>(length)+b]);
        }
    }
    for(int field=0;field<3;++field) for(int b=0;b<length;++b) {
        REAL(VECTOR_ELT(result,field))[b]+=correction[field*static_cast<size_t>(length)+b];
        if(!R_FINITE(REAL(VECTOR_ELT(result,field))[b]))
            Rf_error("Signed Markov grid aggregation overflowed; check observation weights.");
    }
    UNPROTECT(2); return result;
}
