#include <stdio.h>
#include "projections.h"
#include "umfpk.h"
#include "two_grid.h"
#include "residual.h"
#include "gs.h"

double v_cycle(int max_recursion, int n, int m, double L, int *ia, int *ja, double *a, double *b, double *gs_x) {
    double *r = malloc(n * sizeof(double)); // A * x pour pouvoir facilement faire A * x - b par la suite
    double *gs_r = malloc(n * sizeof(double)); // A * x pour pouvoir facilement faire A * x - b par la suite
    double res_gs = residual(&n, &ia, &ja, &a, &b, &gs_x, &gs_r);

    int i;
    int q = (m-1) / 11; // nombre de fois que m-1 est multiple de 11
    int m_coarse = (m-1)/2 + 1;
    int q_coarse = (m_coarse-1) / 11; // nombre de fois que m-1 est multiple de 11
    int p_coarse = 4 * m_coarse - 4; // nombre de points sur le périmètre
    int n_coarse = m_coarse * m_coarse // nombre total de points dans le carré
            - (5 * q_coarse) * (8 * q_coarse) // nombre de points dans le rectangle supérieur droit
            - p_coarse; // number of points on the walls 

    fwd_gs(m, L, &n, &ia, &ja, &a, &b, &gs_x);
    res_gs = residual(&n, &ia, &ja, &a, &b, &gs_x, &gs_r);

    // check that m-1 is a power of 2

    int m_copy = m-1;
    for (int i = 0; i < max_recursion; i++) {
        if (m_copy % 2 != 0) {
            printf("Error: m must be a power of 2\n");
            exit(1);
        } else {
            m_copy /= 2; // go down a level
        }
    }
    // restriction to the coarse grid

    double *restr_r = malloc(n_coarse * sizeof(double));
    restriction(m_coarse, q_coarse, &n_coarse, &gs_r, &restr_r);
    
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

    prolongation(m_coarse, q_coarse, &n_coarse, &r_coarse, &r_prol);

   //construct the x vector from this improvement to the residual

    for (int i = 0; i < n; i++) {
        gs_x[i] += r_prol[i]; // source of bug ????
    }

    double res_prolongation = residual(&n, &ia, &ja, &a, &b, &gs_x, &gs_r);

    bwd_gs(m, L, &n, &ia, &ja, &a, &b, &gs_x);

    double two_grid_residual = residual(&n, &ia, &ja, &a, &b, &gs_x, &gs_r); 

    printf("Two-grid residual: %.10e\n", two_grid_residual);

    free(ia_coarse); free(ja_coarse); free(a_coarse); free(b_coarse); free(r); free(gs_r); free(restr_r); free(r_prol); free(r_coarse);

    return two_grid_residual; // success
}

int generate_multigrid_problem(int max_recursion, int m_fine, int **ia_ptr, int **ja_ptr, double **a_ptr, 
                    double **b_ptr, void **Numeric_ptr) 
{
    // this time ia_ptr are an array to pointers to  the ia arrays for each coarse level in the multigrid problem
    // generates the coarse matrices separately from the multigrid method
    // uses factorize umfpack to store the LU factorization of the coarse matrix in the corresponding pointers in this function's arguments
    // this way umfpack doesn't have to do as much heavy lifting for no purpose

    // check that m-1 is a power of 2

    int m_copy = m_fine-1;

    for (int i = 0; i < max_recursion; i++) { // separate loop to avoid doing a lot of work to find out at the end that m is not acceptable
        if (m_copy % 2 != 0) {
            printf("m_copy: %d\n", m_copy);
            printf("Error: m must be a power of 2\n");
            exit(1);
        } else {
            m_copy /= 2; // go down a level
        }
    }

    for (int i = 0; i < max_recursion; i++) { // iterate over the levels of the multigrid, generate all matrices

        // prepare the arrays for the coarse problem
        int *ia_coarse = ia_ptr[i];
        int *ja_coarse = ja_ptr[i];
        double *a_coarse = a_ptr[i];
        double *b_coarse = b_ptr[i];
        void *Numeric = Numeric_ptr[i];
        // prepare the values for the coarse problem
        int m_coarse = (m_fine-1)/2 + 1;
        int q_coarse = (m_coarse-1) / 11; // nombre de fois que m-1 est multiple de 11
        int p_coarse = 4 * m_coarse - 4; // nombre de points sur le périmètre
        int n_coarse = m_coarse * m_coarse // nombre total de points dans le carré
                - (5 * q_coarse) * (8 * q_coarse) // nombre de points dans le rectangle supérieur droit
                - p_coarse; // number of points on the walls 

        if (prob(m_coarse, &n_coarse, &ia_coarse, &ja_coarse, &a_coarse, &b_coarse, 0))
            return 1;

        Numeric = factorize_umfpack(n_coarse, ia_coarse, ja_coarse, a_coarse);

        m_fine = m_coarse; // we go down to the next level
    }

    return 0;
}
