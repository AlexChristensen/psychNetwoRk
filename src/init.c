#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

// Declare the C functions you want to make available to R here
extern SEXP r_ziggurat(SEXP n, SEXP r_seed);
extern SEXP r_xoshiro_seeds(SEXP n, SEXP r_seed);
extern SEXP r_xoshiro_uniform(SEXP n, SEXP r_seed);
extern SEXP r_xoshiro_shuffle(SEXP r_vector, SEXP r_seed);
extern SEXP r_xoshiro_weighted_shuffle(SEXP r_vector, SEXP r_prob, SEXP r_seed);
extern SEXP r_xoshiro_shuffle_replace(SEXP r_vector, SEXP r_size, SEXP r_seed);
extern SEXP proximity_pass_c(SEXP nodes_r, SEXP ring_r, SEXP budget_r, SEXP pair_rows_r, SEXP pair_cols_r, SEXP pair_counts_r);
extern SEXP swapping_pass_c(SEXP nodes_r, SEXP ring_r, SEXP budget_r, SEXP total_budget_r);


// Register native routine
static const R_CallMethodDef CallEntries[] = {

    {
        "r_ziggurat", // Name of function call in R
        (DL_FUNC)&r_ziggurat, // Name of C function
         2 // Number of arguments
    },
    {
        "r_xoshiro_uniform",
        (DL_FUNC)&r_xoshiro_uniform,
         2
    },
    {
        "r_xoshiro_seeds",
        (DL_FUNC)&r_xoshiro_seeds,
         2
    },
    {
        "r_xoshiro_shuffle",
        (DL_FUNC)&r_xoshiro_shuffle,
         2
    },
    {
        "r_xoshiro_weighted_shuffle",
        (DL_FUNC)&r_xoshiro_weighted_shuffle,
         3
    },
    {
        "r_xoshiro_shuffle_replace",
        (DL_FUNC)&r_xoshiro_shuffle_replace,
         3
    },
    {
        "proximity_pass_c",
        (DL_FUNC)&proximity_pass_c,
         6
    },
    {
        "swapping_pass_c",
        (DL_FUNC)&swapping_pass_c,
         4
    },
    {NULL, NULL, 0}

};

// Set up call in package
void R_init_psychNetwoRk(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
