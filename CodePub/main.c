#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "prob.h"
#include "umfpk.h"
#include "time.h"
#include "residual.h"
#include "rho.h"
#include "heatflux.h"
#include "gs.h"
#include "projections.h"
#include "two_grid.h"
#include "multigrid.h"
//#include "umfpack.h

// Fonction main

int main(int argc, char *argv[])
{
    // déclarer les variables 

    // possible m values for multigrid :
    // 0 : 23
    // 1 : 45
    // 2 : 89
    // 3 : 177
    // 4 : 353
    // 5 : 705
    // 6 : 1409
    // 7 : 2817
    // 8 : 5633
    // 9 : 11265
    // 10 : 22529
    int m =  5633;
    //int two_grid_iter = 17;
    int q = (m-1) / 11;
    int i;
    double L = 5.5;
    double tc1, tc2, tc3, tc4, tc5, tc6, tw1, tw2, tw3, tw4, tw5, tw6, tc7, tw7, tc8, tw8; // mis à jour le 13/10/22
    int n, *ia, *ja; 
    double *a, *b, *x, *gs_x;
    // multigrid declarations
    int max_recursion = 8; // 0 = 2-grid because no recursion happens, we solve directly at the first coarse level
                           // 1 = smooth, restrict, smooth, restrict, solve, prolong, smooth, prolong, smooth
                           // 2 = ...
    int counter;
    int **ia_ptr, **ja_ptr;
    double **a_ptr, **b_ptr;
    void *Numeric;
    double *res_vector;
    int p = 4 * m - 4; // nombre de points sur le périmètre
    int nb_dir = p; // nombre de points sur une porte/fenêtre
    n = m * m // nombre total de points dans le carré
        - (5 * q) * (8 * q) // nombre de points dans le rectangle supérieur droit
        - nb_dir; // number of points on the walls
    ia_ptr = malloc((max_recursion + 2) * sizeof(int *)); // + 2 because max_recursion = 0 is a 2-grid, and then store all of the next grids
    ja_ptr = malloc((max_recursion + 2) * sizeof(int *));
    a_ptr = malloc((max_recursion + 2) * sizeof(double *));
    b_ptr = malloc((max_recursion + 2) * sizeof(double *));

    tc1 = mytimer_cpu(); tw1 = mytimer_wall(); // mis à jour le 13/10/22

    Numeric = generate_multigrid_problem(max_recursion, m, ia_ptr, ja_ptr, a_ptr, b_ptr);

    printf("\nPROBLEM: ");
    printf("m = %5d   n = %8d\n", m, n);

    // around 0.3 seconds for matrix generation for m = 2817 so this is not a costly step

    // allocate memory for the solution vector 

    x = calloc(1, n * sizeof(double));
    gs_x = calloc(1, n * sizeof(double));
    if ( x == NULL ) {
    printf("\n ERREUR : not enough memory for solution vector\n\n");
        return 1;
    }

    double multigrid_residual = 1.0;

    double *gs_r = malloc(n * sizeof(double)); // A * x pour pouvoir facilement faire A * x - b par la suite
    int status = 0;
    res_vector = malloc(50 * sizeof(double));
    res_vector[0] = multigrid_residual;

    // solve the problem using the multigrid method

    counter = 0;
    tc5 = mytimer_cpu(); tw5 = mytimer_wall();
    while ((multigrid_residual > 1.35e-14) && (status == 0)) {
        status = w_cycle(max_recursion, 0, n, m, L, ia_ptr, ja_ptr, a_ptr,
                                                            b_ptr[0], gs_x, Numeric);
        multigrid_residual = residual(&n, &ia_ptr[0], &ja_ptr[0], &a_ptr[0], &b_ptr[0], &gs_x, &gs_r);
        counter++;
        res_vector[counter] = multigrid_residual;
        printf("Iteration %d, multigrid residual = %.10e\n", counter, multigrid_residual);
    }

    tc6 = mytimer_cpu(); tw6 = mytimer_wall();
    printf("\n");

    printf("\nSolution time, multigrid method (CPU): %5.2f sec",tc6-tc5);
    printf("\nSolution time, multigrid method (clock): %5.2f sec \n",tw6-tw5);
    printf("\nSolution time, multigrid method + factorization (CPU): %5.2f sec",tc6-tc1);
    printf("\nSolution time, multigrid method + factorization (clock): %5.2f sec \n",tw6-tw1);
    printf("Number of iterations : %d\n", counter);

    // plot the solution of the multigrid method

    FILE *f_out_multi = fopen("mat/out_mutligrid.dat", "w");

    i = 0;

    for (int iy = 1; iy < m-1; iy++) { // vertical
        for (int ix = 1; ix < m-1; ix++) { // horizontal
            if ((iy < 6 * q || ix < 3 * q)) { // si on n'est pas dans le rectangle supérieur droit
                fprintf(f_out_multi, "%f %f %f\n", iy * L / (q * 11), ix * L / (q * 11), (gs_x)[i]);
                i++; // cycler à travers les éléments de x dans le même ordre qu'ils y ont été placés dans prob.c
            }
        }
        fprintf(f_out_multi, "\n"); // requis par la syntaxe de gnuplot, ligne supplémentaire entre chaque changement de valeur de la 1e colonne (iy dans ce cas-ci)
    }

    fclose(f_out_multi);

    // plot the evolution of the residual norm

    FILE *f_out_multi_res = fopen("mat/out_mutligrid_res.dat", "w");

    for (int i = 0; i < counter; i++) {
        fprintf(f_out_multi_res, "%d %f\n", i, res_vector[i]);
    }

   // RECONSTRUCT PROLONGATION AND RESTRICTION MATRICES IN MATLAB

    //FILE *f_r_iaa = fopen("misc/r_iaa.txt", "w");
    //FILE *f_r_ja = fopen("misc/r_ja.txt", "w");
    //FILE *f_r_a = fopen("misc/r_a.txt", "w");

    //printf("RESTRICTION MATRIX\n");

    //for (int i = 0; i < n; i++) { // debugging restriction and projection matrices
    //    double *e = calloc(1, n * sizeof(double));
    //    e[i] = 1.0; // reset the solution to compare the 2 methods
    //    double *e_restr = malloc(n_coarse * sizeof(double));
    //    restriction(m_coarse, q_coarse, &n_coarse, &e, &e_restr);

    //    for (int j = 0; j < n_coarse; j++) {
    //        if (e_restr[j] == 1.0) {
    //            fprintf(f_r_iaa, "%d\n", j);
    //            fprintf(f_r_ja, "%d\n", i);
    //            fprintf(f_r_a, "%f\n", e_restr[j]);
    //        }
    //    }
    //}

    //FILE *f_p_iaa = fopen("misc/p_iaa.txt", "w");
    //FILE *f_p_ja = fopen("misc/p_ja.txt", "w");
    //FILE *f_p_a = fopen("misc/p_a.txt", "w");

    //printf("PROLONGATION MATRIX\n");

    //for (int i = 0; i < n_coarse; i++) { // debugging restriction and projection matrices
    //    double *e = calloc(1, n_coarse * sizeof(double));
    //    e[i] = 1.0; // reset the solution to compare the 2 methods
    //    double *e_prol = malloc(n * sizeof(double));
    //    prolongation(m_coarse, q_coarse, &n_coarse, &e, &e_prol);

    //    for (int j = 0; j < n; j++) {
    //        if (e_prol[j] != 0.0) {
    //            fprintf(f_p_iaa, "%d\n", j);
    //            fprintf(f_p_ja, "%d\n", i);
    //            fprintf(f_p_a, "%f\n", e_prol[j]);
    //        }
    //    }
    //}

    //fclose(f_r_iaa); fclose(f_r_ja); fclose(f_r_a); fclose(f_p_iaa); fclose(f_p_ja); fclose(f_p_a);

    free(ia); free(ja); free(a); free(b); free(x); free(gs_x);
    //free(ia_ptr); free(ja_ptr); free(a_ptr); free(b_ptr); // causes a double free :(
    //free(ia_coarse); free(ja_coarse); free(a_coarse); free(b_coarse); free(res_vector); //free(x_coarse); free(r_coarse); free(gs_r); // prevents memory leak
    //system("gnuplot -persist \"scripts/heatmap.gnu\""); // laisser gnuplot afficher la température de la pièce
    //system("gnuplot -persist \"scripts/heatmap_two_grid.gnu\"");
    //system("gnuplot -persist \"scripts/heatmap_initial_residual.gnu\"");
    //system("gnuplot -persist \"scripts/heatmap_residual.gnu\""); // laisser gnuplot afficher la température de la pièce
    //system("gnuplot -persist \"scripts/heatmap_restriction.gnu\"");
    //system("gnuplot -persist \"scripts/heatmap_prolongation.gnu\"");
    //system("gnuplot -persist \"scripts/two_grid_residual_plot.gnu\"");

    return 0;
}
