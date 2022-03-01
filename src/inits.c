#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
 Check these declarations against the C/Fortran source code.
 */

/* .Call calls */
extern SEXP _LCP_id_low_search(SEXP);
extern SEXP _LCP_LCP_construction_distance_loop(SEXP, SEXP, SEXP);
extern SEXP _LCP_LCP_construction_path_distance(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _LCP_q_low_compute(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_LCP_id_low_search",                  (DL_FUNC) &_LCP_id_low_search,                  1},
    {"_LCP_LCP_construction_distance_loop", (DL_FUNC) &_LCP_LCP_construction_distance_loop, 3},
    {"_LCP_LCP_construction_path_distance", (DL_FUNC) &_LCP_LCP_construction_path_distance, 6},
    {"_LCP_q_low_compute",                  (DL_FUNC) &_LCP_q_low_compute,                  2},
    {NULL, NULL, 0}
};

void R_init_LCP(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}