// Included inside the native Markov engine's private namespace so all three
// calculations use the same validated, row-scaled fitted transition kernels.
struct AliveQuery {
    const double *times, *query_times, *states;
    int jumps, queries;
    double tau;
    std::vector<double> storage;
    std::vector<const double*> kernels;
    explicit AliveQuery(SEXP spec) {
        SEXP jt = field(spec, "jump_times"), qt = field(spec, "query_times");
        SEXP st = field(spec, "query_state");
        finite_vector(jt, "jump_times"); finite_vector(qt, "query_times");
        finite_vector(st, "query_state");
        if (XLENGTH(jt) > 1000000 || XLENGTH(qt) < 1 || XLENGTH(qt) > 1000000 ||
            XLENGTH(qt) != XLENGTH(st)) fail("Invalid alive-time query dimensions.");
        jumps = static_cast<int>(XLENGTH(jt));
        queries = static_cast<int>(XLENGTH(qt));
        times = REAL(jt); query_times = REAL(qt); states = REAL(st);
        tau = scalar(field(spec, "tau"), "tau");
        if (tau <= 0) fail("The horizon must be positive.");
        for (int j = 0; j < jumps; ++j)
            if (times[j] < 0 || times[j] > tau || (j && times[j] <= times[j-1]))
                fail("Invalid alive-time jump grid.");
        for (int i = 0; i < queries; ++i)
            if (query_times[i] < 0 || query_times[i] > tau || states[i] < 1 ||
                states[i] > 3 || states[i] != std::floor(states[i]))
                fail("Invalid alive-time query time or state.");
        SEXP factors = field(spec, "kernel_factors");
        if (scalar(field(factors,"absorbing"), "absorbing") != 3)
            fail("Expected the third state to be absorbing.");
        kernels.resize(jumps);
        double row_error = 0; size_t adjusted = 0;
        factored_kernels(factors, jumps, 3, storage, kernels, row_error, adjusted);
    }
};

// Upper bound on P(home target satisfied, alive target missed). On any such
// trajectory the final alive time lies in the home interval plus accumulated
// hospital time and at most one remaining horizon of additional hospital time.
SEXP alive_failure_window(SEXP spec, const AliveQuery& x) {
    SEXP qr = field(spec,"q"), ar = field(spec,"accrued"), cr = field(spec,"center");
    SEXP sr = field(spec,"scale");
    finite_vector(qr,"q");
    const int h = x.queries;
    if (XLENGTH(qr) < 1 || XLENGTH(qr) > 1000000/h) fail("Alive-time query matrix is too large.");
    const int nq = static_cast<int>(XLENGTH(qr));
    matrix_shape(ar,h,2,"accrued"); matrix_shape(cr,h,2,"center"); matrix_shape(sr,h,2,"scale");
    const double *q = REAL(qr), *accrued = REAL(ar), *center = REAL(cr), *scale = REAL(sr);
    for(int j=0;j<nq;++j) if(q[j]<0) fail("Score thresholds must be nonnegative.");
    for(int i=0;i<2*h;++i) if(scale[i]<=0) fail("Score scales must be positive.");
    Protect protect;
    SEXP answer = protect.add(R_API(Rf_allocMatrix(REALSXP,h,nq)));
    double* out = REAL(answer);
    std::fill(out,out+h*nq,1.);
    const double padding = 128*std::numeric_limits<double>::epsilon()*(x.jumps+2);
    for (int i=0; i<h; ++i) {
        if ((i & 63)==0) check_interrupt();
        const int state = static_cast<int>(x.states[i])-1;
        if (state == 2 || accrued[i+h] != x.query_times[i]) continue;
        const int first = std::upper_bound(x.times,x.times+x.jumps,x.query_times[i])-x.times;
        std::vector<double> atoms(x.times+first,x.times+x.jumps);
        atoms.push_back(x.tau);
        std::vector<long double> cdf(atoms.size()+1,0);
        long double home = state==0, hospital = state==1;
        for(int j=first;j<x.jumps;++j) {
            const double *k = x.kernels[j];
            const long double death = home*k[6]+hospital*k[7];
            const long double next_home = home*k[0]+hospital*k[1];
            const long double next_hospital = home*k[3]+hospital*k[4];
            cdf[j-first+1]=cdf[j-first]+death;
            home=next_home; hospital=next_hospital;
        }
        cdf.back()=cdf[cdf.size()-2]+home+hospital;
        if(std::fabs(cdf.back()-1)>padding) return R_NilValue;
        std::vector<double> scores(atoms.size());
        for(size_t j=0;j<atoms.size();++j) scores[j]=(atoms[j]-center[i+h])/scale[i+h];
        const double gap=accrued[i+h]-accrued[i];
        for(int j=0;j<nq;++j) {
            const double d=scale[i]*q[j];
            const double rounding=padding*(1+x.tau+std::fabs(center[i])+d);
            const double low=center[i]-d+gap-rounding;
            const double high=center[i]+d+gap+(x.tau-x.query_times[i])+rounding;
            const size_t left=std::lower_bound(atoms.begin(),atoms.end(),low)-atoms.begin();
            const size_t right=std::upper_bound(atoms.begin(),atoms.end(),high)-atoms.begin();
            const size_t below=std::lower_bound(scores.begin(),scores.end(),-q[j])-scores.begin();
            const size_t above=std::upper_bound(scores.begin(),scores.end(),q[j])-scores.begin();
            const size_t last=std::min(right,below), start=std::max(left,above);
            const long double failure=(last>left ? cdf[last]-cdf[left] : 0)+
                                      (right>start ? cdf[right]-cdf[start] : 0);
            out[i+h*j]=std::min(1.,static_cast<double>(failure)+padding);
        }
    }
    return answer;
}

// h[j,s] is the chance, just before jump j in living state s, of an allowed
// death date or of an allowed terminal survival. Conditional transition rows
// are K(s,t)h_after(t) / h_before(s). A zero-denominator row is unreachable
// from every conditioning state with positive success probability.
SEXP alive_conditioned_kernels(SEXP spec, const AliveQuery& x) {
    SEXP accept=field(spec,"allow");
    if(TYPEOF(accept)!=LGLSXP || XLENGTH(accept)!=x.jumps+1) fail("Invalid terminal-atom mask.");
    const int *allowed=LOGICAL(accept);
    for(int j=0;j<=x.jumps;++j) if(allowed[j]==NA_LOGICAL) fail("Missing terminal-atom mask.");
    Protect protect;
    SEXP answer=named_list(protect,{"kernels","h"});
    SEXP matrices=protect.add(R_API(Rf_allocVector(VECSXP,x.jumps)));
    SEXP probability=protect.add(R_API(Rf_allocVector(REALSXP,x.queries)));
    SET_VECTOR_ELT(answer,0,matrices); SET_VECTOR_ELT(answer,1,probability);
    std::vector<long double> h((x.jumps+1)*2);
    h[2*x.jumps]=h[2*x.jumps+1]=allowed[x.jumps] ? 1 : 0;
    for(int j=x.jumps-1;j>=0;--j) {
        if((j & 255)==0) check_interrupt();
        Protect local;
        SEXP m=local.add(R_API(Rf_allocMatrix(REALSXP,3,3)));
        const double *k=x.kernels[j];
        double *cond=REAL(m);
        std::fill(cond,cond+9,0.); cond[8]=1.;
        for(int a=0;a<2;++a) {
            const long double w[]={k[a]*h[2*(j+1)],k[a+3]*h[2*(j+1)+1],
                                   allowed[j] ? static_cast<long double>(k[a+6]) : 0};
            const long double sum=w[0]+w[1]+w[2];
            h[2*j+a]=sum;
            if(sum>0) for(int b=0;b<3;++b) cond[a+3*b]=static_cast<double>(w[b]/sum);
            else cond[a+3*a]=1.;
        }
        SET_VECTOR_ELT(matrices,j,m);
    }
    for(int i=0;i<x.queries;++i) {
        const int a=static_cast<int>(x.states[i])-1;
        const int j=std::upper_bound(x.times,x.times+x.jumps,x.query_times[i])-x.times;
        REAL(probability)[i]=a==2 ? 0 : static_cast<double>(h[2*j+a]);
    }
    return answer;
}
