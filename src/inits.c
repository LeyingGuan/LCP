#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
 Check these declarations against the C/Fortran source code.
 */

/* .Call calls */
extern SEXP _LCPcore_id_low_search(SEXP);
extern SEXP _LCPcore_LCP_construction_distance_loop(SEXP, SEXP, SEXP);
extern SEXP _LCPcore_LCP_construction_path_distance(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _LCPcore_q_low_compute(SEXP, SEXP);
extern SEXP _LCPcore_rcpparma_bothproducts(SEXP);
extern SEXP _LCPcore_rcpparma_hello_world();
extern SEXP _LCPcore_rcpparma_innerproduct(SEXP);
extern SEXP _LCPcore_rcpparma_outerproduct(SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_LCPcore_id_low_search",                  (DL_FUNC) &_LCPcore_id_low_search,                  1},
    {"_LCPcore_LCP_construction_distance_loop", (DL_FUNC) &_LCPcore_LCP_construction_distance_loop, 3},
    {"_LCPcore_LCP_construction_path_distance", (DL_FUNC) &_LCPcore_LCP_construction_path_distance, 6},
    {"_LCPcore_q_low_compute",                  (DL_FUNC) &_LCPcore_q_low_compute,                  2},
    {"_LCPcore_rcpparma_bothproducts",          (DL_FUNC) &_LCPcore_rcpparma_bothproducts,          1},
    {"_LCPcore_rcpparma_hello_world",           (DL_FUNC) &_LCPcore_rcpparma_hello_world,           0},
    {"_LCPcore_rcpparma_innerproduct",          (DL_FUNC) &_LCPcore_rcpparma_innerproduct,          1},
    {"_LCPcore_rcpparma_outerproduct",          (DL_FUNC) &_LCPcore_rcpparma_outerproduct,          1},
    {NULL, NULL, 0}
};

void R_init_LCPcore(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}