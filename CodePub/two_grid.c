#include <stdio.h>
#include "projections.h"
#include "umfpk.h"
#include "two_grid.h"
#include "residual.h"
#include "gs.h"

double two_grid_method(int n, int m, int L, int *ia, int *ja, double *a, double *b, double *x, double *gs_x) {
    double *r = malloc(n * sizeof(double)); // A * x pour pouvoir facilement faire A * x - b par la suite
    double *gs_r = malloc(n * sizeof(double)); // A * x pour pouvoir facilement faire A * x - b par la suite
    double res_gs = residual(&n, &ia, &ja, &a, &b, &gs_x, &gs_r);
    //printf("Initial residual is %.10e\n", res_gs);
    sym_gs(m, L, &n, &ia, &ja, &a, &b, &gs_x);
    sym_gs(m, L, &n, &ia, &ja, &a, &b, &gs_x);
    res_gs = residual(&n, &ia, &ja, &a, &b, &gs_x, &gs_r);
    //printf("Pre-smoothing residual is %.10e\n", res_gs);

    if ((m - 1) % 2 != 0) {
        printf("old_m: %d\n", m);
        printf("Error: old_m must be odd\n");
        exit(1);
    }

    // restriction to the coarse grid
    int m_coarse = (m-1)/2 + 1;
    int q_coarse = (m_coarse-1) / 11; // nombre de fois que m-1 est multiple de 11
    int p_coarse = 4 * m_coarse - 4; // nombre de points sur le périmètre
    int n_coarse = m_coarse * m_coarse // nombre total de points dans le carré
            - (5 * q_coarse) * (8 * q_coarse) // nombre de points dans le rectangle supérieur droit
            - p_coarse; // number of points on the walls 
    //printf("n_coarse: %d\n", n_coarse);
    double *restr_r = malloc(n_coarse * sizeof(double));

    //for (int i = 0; i < n; i++) {
    //    gs_r[i] = 1.0; // test vector for debugging
    //}

    restriction(m_coarse, q_coarse, &n_coarse, &ia, &ja, &a, &b, &gs_x, &gs_r, &restr_r);
    
    //for (int i = 0; i < n_coarse; i++) { // debug
    //    printf("restr_r[%d] = %f\n", i, restr_r[i]);
    //}

    // solve the coarse problem
    // generate the coarse matrix

    int *ia_coarse, *ja_coarse; 
    double *a_coarse, *b_coarse, *r_coarse;
    if (prob(m_coarse, &n_coarse, &ia_coarse, &ja_coarse, &a_coarse, &b_coarse, 0))
        return 1;

    r_coarse = malloc(n_coarse * sizeof(double));

    if (solve_umfpack(n_coarse, ia_coarse, ja_coarse, a_coarse, restr_r, r_coarse)) { // rh side is restr_r, we are solving A * x = b_coarse - A * x
                                                                                      // which we will use to compute the prolongation and go back to the fine grid
        free(ia_coarse); free(ja_coarse); free(a_coarse); free(b_coarse); free(r_coarse); // prevents memory leak
        return 1;
    }

    double *r_prol = malloc(n * sizeof(double));

    //for (int i = 0; i < n_coarse; i++) { // debug
    //    r_coarse[i] = 1.0; // test vector for debugging
    //}

    prolongation(m_coarse, q_coarse, &n_coarse, &r_coarse, &r_prol);

    //construct the x vector from this improvement to the residual

    for (int i = 0; i < n; i++) {
        //printf("gs_x[%d]: %f\n", i, gs_x[i]); // debug
        //printf("r_prol[%d]: %f\n", i, r_prol[i]); // debug
        gs_x[i] += r_prol[i];
    }

    double res_prolongation = residual(&n, &ia, &ja, &a, &b, &gs_x, &gs_r);

    //printf("Post-prolongation residual is %.10e\n", res_prolongation);

    sym_gs(m, L, &n, &ia, &ja, &a, &b, &gs_x);
    sym_gs(m, L, &n, &ia, &ja, &a, &b, &gs_x);    

    double two_grid_residual = residual(&n, &ia, &ja, &a, &b, &gs_x, &gs_r); 

    printf("Two-grid residual: %.10e\n", two_grid_residual);

    free(ia_coarse); free(ja_coarse); free(a_coarse); free(b_coarse); free(r); free(gs_r); free(restr_r); 

    return two_grid_residual; // success
}

int generate_coarse_problem(int m, int *ia_coarse_ptr, int *ja_coarse_ptr, double *a_coarse_ptr, 
                    double *b_coarse_ptr, void **Symbolic, void **Numeric, double *Info, double *Control) 
{
    // generates the coarse matrices separately from the two-grid method
    // uses factorize umfpack to store the LU factorization of the coarse matrix in the corresponding pointers in this function's arguments
    // this way umfpack doesn't have to do as much heavy lifting for no purpose
    int m_coarse = (m-1)/2 + 1;
    int q_coarse = (m_coarse-1) / 11; // nombre de fois que m-1 est multiple de 11
    int p_coarse = 4 * m_coarse - 4; // nombre de points sur le périmètre
    int n_coarse = m_coarse * m_coarse // nombre total de points dans le carré
            - (5 * q_coarse) * (8 * q_coarse) // nombre de points dans le rectangle supérieur droit
            - p_coarse; // number of points on the walls 
    //printf("n_coarse: %d\n", n_coarse);

    if (prob(m_coarse, &n_coarse, &ia_coarse_ptr, &ja_coarse_ptr, &a_coarse_ptr, &b_coarse_ptr, 0))
        return 1;

    if (factorize_umfpack(n_coarse, ia_coarse_ptr, ja_coarse_ptr, a_coarse_ptr, Symbolic, Numeric, Info, Control)) 
        return 1;
    
    return 0;
}