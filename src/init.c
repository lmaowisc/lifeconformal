#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

extern SEXP life_nextafter(SEXP);
extern SEXP life_occupation_complete(SEXP, SEXP, SEXP, SEXP);
extern SEXP life_calibration_sum_create(SEXP, SEXP);
extern SEXP life_calibration_sum_add(SEXP, SEXP, SEXP);
extern SEXP life_calibration_sum_value(SEXP, SEXP);
extern SEXP life_calibration_sum_reset(SEXP);
extern SEXP life_calibration_sum_close(SEXP);
extern SEXP life_calibration_sum_info(void);
extern SEXP life_cdt_sweep(SEXP, SEXP);
extern SEXP life_markov_monotone_batch(SEXP);
extern SEXP life_markov_monotone_queries(SEXP);
extern SEXP life_markov_grid_native(SEXP);
extern SEXP life_markov_alive_conditional(SEXP, SEXP);
extern SEXP life_markov_marginal_reduce(SEXP);
extern SEXP life_forest_rates(SEXP, SEXP, SEXP);
extern SEXP life_forest_rates_cached(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP life_flexible_design(SEXP, SEXP, SEXP, SEXP);
extern SEXP life_forest_block_info(void);
extern SEXP life_forest_block_supported(SEXP, SEXP, SEXP, SEXP);
extern SEXP life_forest_block_advance(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP life_forest_bank_supported(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP life_forest_bank_create(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP life_forest_bank_fill(SEXP, SEXP);
extern SEXP life_forest_bank_seal(SEXP);
extern SEXP life_forest_bank_copy(SEXP);
extern SEXP life_forest_bank_release(SEXP);
extern SEXP life_forest_bank_inspect(SEXP);
extern SEXP life_forest_bank_live(void);
extern SEXP life_forest_bank_advance(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern void life_forest_bank_unload_check(void);
extern void life_forest_bank_finalize_all(void);
extern void life_calibration_sum_unload_check(void);
extern void life_calibration_sum_finalize_all(void);

void attribute_visible R_unload_lifeconformal(DllInfo *dll);

/* The loader resolves unload hooks through registered symbols when dynamic
 * lookup is disabled. Its void(DllInfo*) call has one pointer argument. */
static const R_CMethodDef CEntries[] = {
    {"R_unload_lifeconformal", (DL_FUNC) &R_unload_lifeconformal, 1, NULL},
    {NULL, NULL, 0, NULL}
};

static const R_CallMethodDef CallEntries[] = {
    {"life_nextafter", (DL_FUNC) &life_nextafter, 1},
    {"life_occupation_complete", (DL_FUNC) &life_occupation_complete, 4},
    {"life_calibration_sum_create", (DL_FUNC) &life_calibration_sum_create, 2},
    {"life_calibration_sum_add", (DL_FUNC) &life_calibration_sum_add, 3},
    {"life_calibration_sum_value", (DL_FUNC) &life_calibration_sum_value, 2},
    {"life_calibration_sum_reset", (DL_FUNC) &life_calibration_sum_reset, 1},
    {"life_calibration_sum_close", (DL_FUNC) &life_calibration_sum_close, 1},
    {"life_calibration_sum_info", (DL_FUNC) &life_calibration_sum_info, 0},
    {"life_cdt_sweep", (DL_FUNC) &life_cdt_sweep, 2},
    {"life_markov_monotone_batch", (DL_FUNC) &life_markov_monotone_batch, 1},
    {"life_markov_monotone_queries", (DL_FUNC) &life_markov_monotone_queries, 1},
    {"life_markov_alive_conditional", (DL_FUNC) &life_markov_alive_conditional, 2},
    {"life_markov_grid_native", (DL_FUNC) &life_markov_grid_native, 1},
    {"life_markov_marginal_reduce", (DL_FUNC) &life_markov_marginal_reduce, 1},
    {"life_forest_rates", (DL_FUNC) &life_forest_rates, 3},
    {"life_forest_rates_cached", (DL_FUNC) &life_forest_rates_cached, 5},
    {"life_flexible_design", (DL_FUNC) &life_flexible_design, 4},
    {"life_forest_block_info", (DL_FUNC) &life_forest_block_info, 0},
    {"life_forest_block_supported", (DL_FUNC) &life_forest_block_supported, 4},
    {"life_forest_block_advance", (DL_FUNC) &life_forest_block_advance, 6},
    {"life_forest_bank_supported", (DL_FUNC) &life_forest_bank_supported, 6},
    {"life_forest_bank_create", (DL_FUNC) &life_forest_bank_create, 7},
    {"life_forest_bank_fill", (DL_FUNC) &life_forest_bank_fill, 2},
    {"life_forest_bank_seal", (DL_FUNC) &life_forest_bank_seal, 1},
    {"life_forest_bank_copy", (DL_FUNC) &life_forest_bank_copy, 1},
    {"life_forest_bank_release", (DL_FUNC) &life_forest_bank_release, 1},
    {"life_forest_bank_inspect", (DL_FUNC) &life_forest_bank_inspect, 1},
    {"life_forest_bank_live", (DL_FUNC) &life_forest_bank_live, 0},
    {"life_forest_bank_advance", (DL_FUNC) &life_forest_bank_advance, 6},
    {NULL, NULL, 0}
};

void attribute_visible R_init_lifeconformal(DllInfo *dll) {
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

void attribute_visible R_unload_lifeconformal(DllInfo *dll) {
    /* Reject a reentrant unload before destroying either registry. Explicit
     * weakref finalization clears callbacks before their DLL is removed. */
    life_forest_bank_unload_check();
    life_calibration_sum_unload_check();
    life_forest_bank_finalize_all();
    life_calibration_sum_finalize_all();
    (void)dll;
}
