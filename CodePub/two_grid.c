#include <stdio.h>
#include "projections.h"
#include "umfpk.h"
#include "two_grid.h"
#include "residual.h"
#include "gs.h"

double two_grid_method(int n, int m, double L, int *ia, int *ja, double *a, double *b, double *gs_x) {
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
    //printf("n_coarse: %d\n", n_coarse);   printf("Initial residual is %.10e\n", res_gs);

    //FILE *f_out_initial_residual = fopen("mat/out_initial_residual.dat", "w");

    //i = 0;

    //for (int iy = 1; iy < m-1; iy++) { // vertical
    //    for (int ix = 1; ix < m-1; ix++) { // horizontal
    //        if ((iy < 6 * q || ix < 3 * q)) { // si on n'est pas dans le rectangle supérieur droit
    //            fprintf(f_out_initial_residual, "%f %f %f\n", iy * L / (q * 11), ix * L / (q * 11), (gs_r)[i]);
    //            i++; // cycler à travers les éléments de x dans le même ordre qu'ils y ont été placés dans prob.c
    //        }
    //    }
    //    fprintf(f_out_initial_residual, "\n"); // requis par la syntaxe de gnuplot, ligne supplémentaire entre chaque changement de valeur de la 1e colonne (iy dans ce cas-ci)
    //}

    //fclose(f_out_initial_residual);

    fwd_gs(m, L, &n, &ia, &ja, &a, &b, &gs_x);
    //fwd_gs(m, L, &n, &ia, &ja, &a, &b, &gs_x);
    res_gs = residual(&n, &ia, &ja, &a, &b, &gs_x, &gs_r);
    //printf("Pre-smoothing residual is %.10e\n", res_gs);

    //for (int i = 0; i < n; i++) {
    //    printf("gs_r[%d] = %f", i, (gs_r)[i]); // test vector for debugging
    //}

    if ((m - 1) % 2 != 0) {
        printf("old_m: %d\n", m);
        printf("Error: old_m must be odd\n");
        exit(1);
    }

    // restriction to the coarse grid

    double *restr_r = malloc(n_coarse * sizeof(double));

    //for (int i = 0; i < n; i++) {
    //    gs_r[i] = 1.0; // test vector for debugging
    //}

    //FILE *f_out_residual = fopen("mat/out_residual.dat", "w");

    //i = 0;

    //for (int iy = 1; iy < m-1; iy++) { // vertical
    //    for (int ix = 1; ix < m-1; ix++) { // horizontal
    //        if ((iy < 6 * q || ix < 3 * q)) { // si on n'est pas dans le rectangle supérieur droit
    //            fprintf(f_out_residual, "%f %f %f\n", iy * L / (q * 11), ix * L / (q * 11), (gs_r)[i]);
    //            i++; // cycler à travers les éléments de x dans le même ordre qu'ils y ont été placés dans prob.c
    //        }
    //    }
    //    fprintf(f_out_residual, "\n"); // requis par la syntaxe de gnuplot, ligne supplémentaire entre chaque changement de valeur de la 1e colonne (iy dans ce cas-ci)
    //}

    //fclose(f_out_residual);    

    restriction(m_coarse, q_coarse, &n_coarse, &gs_r, &restr_r);
    
    // plot the restricted residue for debugging

    //FILE *f_out_restriction = fopen("mat/out_restriction.dat", "w");

    //i = 0;

    //for (int iy = 1; iy < m_coarse-1; iy++) { // vertical
    //    for (int ix = 1; ix < m_coarse-1; ix++) { // horizontal
    //        if ((iy < 6 * q_coarse || ix < 3 * q_coarse)) { // si on n'est pas dans le rectangle supérieur droit
    //            fprintf(f_out_restriction, "%f %f %f\n", iy * L / (q_coarse * 11), ix * L / (q_coarse * 11), (restr_r)[i]);
    //            i++; // cycler à travers les éléments de x dans le même ordre qu'ils y ont été placés dans prob.c
    //        }
    //    }
    //    fprintf(f_out_restriction, "\n"); // requis par la syntaxe de gnuplot, ligne supplémentaire entre chaque changement de valeur de la 1e colonne (iy dans ce cas-ci)
    //}

    //fclose(f_out_restriction);
    //for (int i = 0; i < n_coarse; i++) { // debug
    //    printf("restr_r[%d] = %f\n", i, restr_r[i]);
    //}

    // solve the coarse problem
    // generate the coarse matrix

    int *ia_coarse, *ja_coarse; 
    double *a_coarse, *b_coarse, *r_coarse;
    if (prob(m_coarse, &n_coarse, &ia_coarse, &ja_coarse, &a_coarse, &b_coarse, 0))
        return 1;

    //FILE *f_out_coarse = fopen("mat/out_coarse.dat", "w");

    //for (int i = 0; i < n_coarse; i++) { // debug
    //    for (int j = ia_coarse[i]; j < ia_coarse[i+1]; j++) {
    //        for (int k = 0; k < n_coarse; k++) {
    //            if (ja_coarse[j] == k) {
    //                fprintf(f_out_coarse, "%f\n", a_coarse[j]);
    //            } else {
    //                fprintf(f_out_coarse, "%f\n", 0.0);
    //            }
    //        }
    //    }
    //}

    //fclose(f_out_coarse);

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

    double *r_prol_no_solve = malloc(n * sizeof(double)); // debug

    //prolongation(m_coarse, q_coarse, &n_coarse, &restr_r, &r_prol_no_solve);

    prolongation(m_coarse, q_coarse, &n_coarse, &r_coarse, &r_prol);

    //FILE *f_out_prolongation = fopen("mat/out_prolongation.dat", "w");

    //i = 0;

    //for (int iy = 1; iy < m-1; iy++) { // vertical
    //    for (int ix = 1; ix < m-1; ix++) { // horizontal
    //        if ((iy < 6 * q || ix < 3 * q)) { // si on n'est pas dans le rectangle supérieur droit
    //            fprintf(f_out_prolongation, "%f %f %f\n", iy * L / (q * 11), ix * L / (q * 11), (r_prol_no_solve)[i]);
    //            i++; // cycler à travers les éléments de x dans le même ordre qu'ils y ont été placés dans prob.c
    //        }
    //    }
    //    fprintf(f_out_prolongation, "\n"); // requis par la syntaxe de gnuplot, ligne supplémentaire entre chaque changement de valeur de la 1e colonne (iy dans ce cas-ci)
    //}

    //fclose(f_out_prolongation); 

    //construct the x vector from this improvement to the residual

    for (int i = 0; i < n; i++) {
        //printf("gs_x[%d]: %f\n", i, gs_x[i]); // debug
        //printf("r_prol[%d]: %f\n", i, r_prol[i]); // debug
        gs_x[i] += r_prol[i]; // source of bug ????
    }

    double res_prolongation = residual(&n, &ia, &ja, &a, &b, &gs_x, &gs_r);

    //printf("Post-prolongation residual is %.10e\n", res_prolongation);

    bwd_gs(m, L, &n, &ia, &ja, &a, &b, &gs_x);
    //bwd_gs(m, L, &n, &ia, &ja, &a, &b, &gs_x);    

    double two_grid_residual = residual(&n, &ia, &ja, &a, &b, &gs_x, &gs_r); 

    printf("Two-grid residual: %.10e\n", two_grid_residual);

    free(ia_coarse); free(ja_coarse); free(a_coarse); free(b_coarse); free(r); free(gs_r); free(restr_r); free(r_prol); free(r_prol_no_solve); free(r_coarse);

    return two_grid_residual; // success
}

int generate_coarse_problem(int m_fine, int *ia_coarse_ptr, int *ja_coarse_ptr, double *a_coarse_ptr, 
                    double *b_coarse_ptr, void *Numeric) 
{
    // generates the coarse matrices separately from the two-grid method
    // uses factorize umfpack to store the LU factorization of the coarse matrix in the corresponding pointers in this function's arguments
    // this way umfpack doesn't have to do as much heavy lifting for no purpose
    int m_coarse = (m_fine-1)/2 + 1;
    int q_coarse = (m_coarse-1) / 11; // nombre de fois que m-1 est multiple de 11
    int p_coarse = 4 * m_coarse - 4; // nombre de points sur le périmètre
    int n_coarse = m_coarse * m_coarse // nombre total de points dans le carré
            - (5 * q_coarse) * (8 * q_coarse) // nombre de points dans le rectangle supérieur droit
            - p_coarse; // number of points on the walls 
    //printf("n_coarse: %d\n", n_coarse);

    if (prob(m_coarse, &n_coarse, &ia_coarse_ptr, &ja_coarse_ptr, &a_coarse_ptr, &b_coarse_ptr, 0))
        return 1;

    Numeric = factorize_umfpack(n_coarse, ia_coarse_ptr, ja_coarse_ptr, a_coarse_ptr);
    
    return 0;
}

double factorized_two_grid_method(int n, int m, double L, int *ia, int *ja, double *a, double *b, double *gs_x, int *ia_coarse, int *ja_coarse, double *a_coarse, 
                    double *b_coarse, void *Numeric) {
    // uses the factorized coarse matrix to solve the coarse problem
    // then uses the solution to the coarse problem to improve the solution to the fine problem
    // this is the two-grid method with factorized coarse matrix
    double *r = malloc(n * sizeof(double)); // A * x pour pouvoir facilement faire A * x - b par la suite
    double *gs_r = malloc(n * sizeof(double)); // A * x pour pouvoir facilement faire A * x - b par la suite
    double res_gs; //= residual(&n, &ia, &ja, &a, &b, &gs_x, &gs_r);
    //printf("Initial residual is %.10e\n", res_gs);
    fwd_gs(m, L, &n, &ia, &ja, &a, &b, &gs_x); 
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
            
    double *restr_r = malloc(n_coarse * sizeof(double));
    restriction(m_coarse, q_coarse, &n_coarse, &gs_r, &restr_r);

    double *r_coarse = malloc(n_coarse * sizeof(double));

    if (solve_umfpack_factorized(n_coarse, ia_coarse, ja_coarse, a_coarse, restr_r, r_coarse, Numeric)) { // rh side is restr_r, we are solving A * x = b_coarse - A * x
                                                                                      // which we will use to compute the prolongation and go back to the fine grid  
        free(r_coarse); // prevents memory leak
        sleep(1); // debug
        return 1;
    }

    double *r_prol = malloc(n * sizeof(double));
    prolongation(m_coarse, q_coarse, &n_coarse, &r_coarse, &r_prol);

    //construct the x vector from this improvement to the residual

    for (int i = 0; i < n; i++) {
        gs_x[i] += r_prol[i];
    }

    bwd_gs(m, L, &n, &ia, &ja, &a, &b, &gs_x);

    double two_grid_residual = residual(&n, &ia, &ja, &a, &b, &gs_x, &gs_r);

    printf("Two-grid residual: %.10e\n", two_grid_residual);

    free(r); free(gs_r); free(restr_r); free(r_coarse); free(r_prol);

    return two_grid_residual; // success

}