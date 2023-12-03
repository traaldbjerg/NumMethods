#include <stdio.h>
#include "projections.h"
#include "umfpk.h"
#include "two_grid.h"
#include "residual.h"
#include "gs.h"

int two_grid_method(int n, int m, double L, int **ia_ptr, int **ja_ptr, double **a_ptr, double *b, double *x, void *Numeric) {
    // uses the factorized coarse matrix to solve the coarse problem
    // then uses the solution to the coarse problem to improve the solution to the fine problem
    // this is the two-grid method with factorized coarse matrix
    double *r = malloc(n * sizeof(double)); // A * x pour pouvoir facilement faire A * x - b par la suite
    //printf("Initial residual is %.10e\n", res_gs);
    fwd_gs(m, L, &n, &ia_ptr[0], &ja_ptr[0], &a_ptr[0], &b, &x); 
    (void) residual(&n, &ia_ptr[0], &ja_ptr[0], &a_ptr[0], &b, &x, &r);
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
            
    double *restr_r = malloc(n_coarse * sizeof(double));
    restriction(m_coarse, q_coarse, &n_coarse, &r, &restr_r);

    double *r_coarse = malloc(n_coarse * sizeof(double));

    if (solve_umfpack_factorized(n_coarse, ia_ptr[1], ja_ptr[1], a_ptr[1], restr_r, r_coarse, Numeric)) { // rh side is restr_r, we are solving A * x = b_coarse - A * x
                                                                                      // which we will use to compute the prolongation and go back to the fine grid  
        free(r_coarse); // prevents memory leak
        sleep(1); // debug
        return 1;
    }

    double *r_prol = malloc(n * sizeof(double));
    prolongation(m_coarse, q_coarse, &n_coarse, &r_coarse, &r_prol);

    //construct the x vector from this improvement to the residual

    for (int i = 0; i < n; i++) {
        x[i] += r_prol[i];
    }

    bwd_gs(m, L, &n, &ia_ptr[0], &ja_ptr[0], &a_ptr[0], &b, &x);

    //double two_grid_residual = residual(&n, &ia_ptr[0], &ja_ptr[0], &a_ptr[0], &b, &x, &r);

    //printf("Two-grid residual: %.10e\n", two_grid_residual);

    free(r); free(restr_r); free(r_coarse); free(r_prol);

    return 0; // success

}


int generate_two_grid_problem(int m_fine, int **ia_ptr, int **ja_ptr, double **a_ptr, 
                    double **b_ptr, void *Numeric) 
{
    // generates the coarse matrices separately from the two-grid method
    // uses factorize umfpack to store the LU factorization of the coarse matrix in the corresponding pointers in this function's arguments
    // this way umfpack doesn't have to do as much heavy lifting for no purpose
    int m, n, p, q;

    for (int i = 0; i < 2; i++) {
        if (prob(m, &n, &ia_ptr[i], &ja_ptr[i], &a_ptr[i], &b_ptr[i], 0))
            return 1;
        m = (m_fine-1)/2 + 1;
        q = (m-1) / 11; // nombre de fois que m-1 est multiple de 11
        p = 4 * m - 4; // nombre de points sur le périmètre
        n = m * m // nombre total de points dans le carré
            - (5 * q) * (8 * q) // nombre de points dans le rectangle supérieur droit
            - p; // number of points on the walls
   }
    
    Numeric = factorize_umfpack(n, ia_ptr[1], ja_ptr[1], a_ptr[1]);
    
    return 0;
}
