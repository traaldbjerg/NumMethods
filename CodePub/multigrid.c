#include <stdio.h>
#include "projections.h"
#include "umfpk.h"
#include "multigrid.h"
#include "residual.h"
#include "gs.h"
#include "prob.h"

int v_cycle(int max_recursion, int c, int n, int m, double L, int **ia_ptr, int **ja_ptr,
                             double **a_ptr, double *b, double *x, void *Numeric) {

    // this function is called recursively, c indicates the level of recursion
    // the first recursion starts at c = 0, and stops at current recursion = max_recursion
    // the last recursion is the coarsest grid, and the first recursion is the finest grid
    // uses the factorized coarsest problem matrix to solve the coarsest problem
    // then uses the solution to the coarse problem to improve the solution to the fine problem

    // x is the solution vector
    // b is the right hand side
    // ia, ja, a are the CSR matrix
    // n is the number of rows (or columns) of the matrix
    // Numeric is the pointer to the umfpack factorization of the coarse grid matrix
    // max_recursion is the level of recursions to do
    // c is the current recursion level

    double *r = malloc(n * sizeof(double)); // A * x pour pouvoir facilement faire A * x - b par la suite
    double res_gs; //= residual(&n, &ia, &ja, &a, &b, &x, &r);
    fwd_gs(m, L, &n, &ia_ptr[c], &ja_ptr[c], &a_ptr[c], &b, &x); 
    res_gs = residual(&n, &ia_ptr[c], &ja_ptr[c], &a_ptr[c], &b, &x, &r);

    //printf("Post-smoothing residual: %.10e\n", res_gs);

    // restriction to the coarse grid
    int m_coarse = (m-1)/2 + 1;
    int q_coarse = (m_coarse-1) / 11; // nombre de fois que m-1 est multiple de 11
    int p_coarse = 4 * m_coarse - 4; // nombre de points sur le périmètre
    int n_coarse = m_coarse * m_coarse // nombre total de points dans le carré
            - (5 * q_coarse) * (8 * q_coarse) // nombre de points dans le rectangle supérieur droit
            - p_coarse; // number of points on the walls
            
    double *restr_r = malloc(n_coarse * sizeof(double));

    restriction(m_coarse, q_coarse, &n_coarse, &r, &restr_r);

    double *correction = calloc(1, n_coarse * sizeof(double));

    if (c == max_recursion) { // we are at the coarsest grid, we can solve the problem directly
        //printf("I am solving now, c = %d\n", c); //debug
        if (solve_umfpack_factorized(n_coarse, ia_ptr[c+1], ja_ptr[c+1], a_ptr[c+1], restr_r, correction, Numeric)) { // rh side is restr_r, we are solving A * x = b_coarse - A * x
                                                                                          // which we will use to compute the prolongation and go back to the fine grid  
            free(correction); // prevents memory leak
            //sleep(1); // debug
            return 1;
        }
    } else { // go to the coarser level
        if (v_cycle(max_recursion, c + 1, n_coarse, m_coarse, L, ia_ptr, ja_ptr, a_ptr, restr_r, correction, Numeric)) {
            free(correction); // prevents memory leak
            //sleep(1); // debug
            return 1;
        }
    }

    double *correction_prol = malloc(n * sizeof(double));

    prolongation(m_coarse, q_coarse, &n_coarse, &correction, &correction_prol);

    //construct the x vector from this improvement to the residual

    for (int i = 0; i < n; i++) {
        x[i] += correction_prol[i];
    }

    bwd_gs(m, L, &n, &ia_ptr[c], &ja_ptr[c], &a_ptr[c], &b, &x);

    //double multigrid_residual = residual(&n, &ia_ptr[c], &ja_ptr[c], &a_ptr[c], &b, &x, &r);

    //printf("Multigrid residual: %.10e\n", multigrid_residual);

    free(r); free(restr_r); free(correction); free(correction_prol);

    return 0; // success

}



int w_cycle(int max_recursion, int c, int n, int m, double L, int **ia_ptr, int **ja_ptr,
                             double **a_ptr, double *b, double *x, void *Numeric) {

    // this function is called recursively, c indicates the level of recursion
    // the first recursion starts at c = 0, and stops at current recursion = max_recursion - 1
    // the last recursion is the coarsest grid, and the first recursion is the finest grid
    // uses the factorized coarsest problem matrix to solve the coarsest problem
    // then uses the solution to the coarse problem to improve the solution to the fine problem

    // x is the solution vector
    // b is the right hand side
    // ia, ja, a are the CSR matrix
    // n is the number of rows (or columns) of the matrix
    // Numeric is the pointer to the umfpack factorization of the coarse grid matrix
    // max_recursion is the level of recursions to do
    // c is the current recursion level

    double *r = malloc(n * sizeof(double)); // A * x pour pouvoir facilement faire A * x - b par la suite
    double res_gs; //= residual(&n, &ia, &ja, &a, &b, &x, &r);
    fwd_gs(m, L, &n, &ia_ptr[c], &ja_ptr[c], &a_ptr[c], &b, &x); 
    res_gs = residual(&n, &ia_ptr[c], &ja_ptr[c], &a_ptr[c], &b, &x, &r);

    //printf("This is starting c: %d\n", c); // check to see if it is actually a w-cycle or not
    //printf("Post-smoothing residual: %.10e\n", res_gs);

    // restriction to the coarse grid
    int m_coarse = (m-1)/2 + 1;
    int q_coarse = (m_coarse-1) / 11; // nombre de fois que m-1 est multiple de 11
    int p_coarse = 4 * m_coarse - 4; // nombre de points sur le périmètre
    int n_coarse = m_coarse * m_coarse // nombre total de points dans le carré
            - (5 * q_coarse) * (8 * q_coarse) // nombre de points dans le rectangle supérieur droit
            - p_coarse; // number of points on the walls
            
    double *restr_r = malloc(n_coarse * sizeof(double));

    restriction(m_coarse, q_coarse, &n_coarse, &r, &restr_r);

    double *correction = calloc(1, n_coarse * sizeof(double));

    //printf("This is restricted c: %d\n", c); // check to see if it is actually a w-cycle or not

    if (c == max_recursion) { // we are at the coarsest grid, we can solve the problem directly
        //printf("This is solved c = %d\n", c + 1); // debug, +1 so that the output is less confusing
        if (solve_umfpack_factorized(n_coarse, ia_ptr[c+1], ja_ptr[c+1], a_ptr[c+1], restr_r, correction, Numeric)) { // rh side is restr_r, we are solving A * x = b_coarse - A * x
                                                                                          // which we will use to compute the prolongation and go back to the fine grid  
            free(correction); // prevents memory leak
            //sleep(1); // debug
            return 1;
        }
    } else { // go to the coarser level BUT DO IT TWICE FOR THE W CYCLE
        if (w_cycle(max_recursion, c + 1, n_coarse, m_coarse, L, ia_ptr, ja_ptr, a_ptr, restr_r, correction, Numeric)) {
            free(correction); // prevents memory leak
            //sleep(1); // debug
            return 1;
        }
        //printf("This is intermediary c = %d\n", c); // check to see if it is actually a w-cycle or not
                                                      // the output of this one is a bit confusing however, because it goes up to +1
                                                      // but nothing happens at that level
        if (w_cycle(max_recursion, c + 1, n_coarse, m_coarse, L, ia_ptr, ja_ptr, a_ptr, restr_r, correction, Numeric)) {
            free(correction); // prevents memory leak
            //sleep(1); // debug
            return 1;
        }
    }

    double *correction_prol = malloc(n * sizeof(double));

    prolongation(m_coarse, q_coarse, &n_coarse, &correction, &correction_prol);

    //printf("This is prolonged c: %d\n", c); // check to see if it is actually a w-cycle or not

    //construct the x vector from this improvement to the residual

    for (int i = 0; i < n; i++) {
        x[i] += correction_prol[i];
    }

    bwd_gs(m, L, &n, &ia_ptr[c], &ja_ptr[c], &a_ptr[c], &b, &x);

    //double multigrid_residual = residual(&n, &ia_ptr[c], &ja_ptr[c], &a_ptr[c], &b, &x, &r);

    //printf("Multigrid residual: %.10e\n", multigrid_residual);

    free(r); free(restr_r); free(correction); free(correction_prol);

    //printf("This is final c: %d\n", c); // check to see if it is actually a w-cycle or not

    return 0; // success

}



void *generate_multigrid_problem(int max_recursion, int m, int **ia_ptr, int **ja_ptr, double **a_ptr, 
                    double **b_ptr) 
{
    // this time ia_ptr are an array to pointers to  the ia arrays for each coarse level in the multigrid problem
    // generates the coarse matrices separately from the multigrid method
    // uses factorize umfpack to store the LU factorization of the coarse matrix in the corresponding pointers in this function's arguments
    // this way umfpack doesn't have to do as much heavy lifting for no purpose

    // ia_ptr[0], ja_ptr[0], a_ptr[0] are the fine grid matrices
    // ia_ptr[1], ja_ptr[1], a_ptr[1] are the coarse grid matrices
    // ia_ptr[2], ja_ptr[2], a_ptr[2] are the next level coarse grid matrices
    // etc...
    // b_ptr[0] is the fine grid rhs
    // b_ptr[1] is the coarse grid rhs (not used)
    // b_ptr[2] is the next level coarse grid rhs (not used)
    // etc...
    // max_recursion is the level of recursions to do, indicating how many coarse levels there are

    // check that m-1 is a power of 2

    int m_copy = m-1;
    void *Numeric; // will be used to store the LU factorization of the coarse matrix

    for (int i = 0; i < max_recursion + 1; i++) { // separate loop to avoid doing a lot of work to find out at the end that m is not acceptable

    //printf("hello: %d\n", m_copy);

        if (m_copy % 2 != 0) {
            printf("m_copy: %d\n", m_copy);
            printf("Error: m must be a power of 2\n");
            exit(1);
        } else {
            m_copy /= 2; // go down a level
        }
    }

    for (int i = 0; i <= max_recursion + 1; i++) { // iterate over the levels of the multigrid, generate all matrices
        // prepare the values for the coarse problem
        int q = (m-1) / 11; // nombre de fois que m-1 est multiple de 11
        int p = 4 * m - 4; // nombre de points sur le périmètre
        int n = m * m // nombre total de points dans le carré
                - (5 * q) * (8 * q) // nombre de points dans le rectangle supérieur droit
                - p; // number of points on the walls 

        if (prob(m, &n, &ia_ptr[i], &ja_ptr[i], &a_ptr[i], &b_ptr[i], 0))
            return 1;
        
        if (i == max_recursion + 1) // we only need to build the last LU factorization
            Numeric = factorize_umfpack(n, ia_ptr[i], ja_ptr[i], a_ptr[i]); 

        m = (m-1)/2 + 1; // go down a level, might be a problem if the lowest level is m = 12
    }

    return Numeric;

}
