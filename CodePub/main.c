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
//#include "umfpack.h

// Fonction main

int main(int argc, char *argv[])
{
    
    // déclarer les variables 

    //double (*rho_ptr)(double, double, double) = &rho;
    int m = 155;
    //int u = 0;
    //int iter_max = 0; // régler le nombre d'itérations, mettre à 0 si on ne cherche pas à minimiser std_dev/avrg
    int two_grid_iter = 17;
    int q = (m-1) / 11;
    int i;
    //int use_petsc = 0; // activer ou déactiver l'utilisation de PETSc
    //double source_value = 500.0 ; // permet d'itérer sur les différentes valeurs de rho pour 
    //double flux_x, flux_y, rad_flux;
    //double source_save;
    //double save_dev;
    //double avrg, std_dev;
    //int dim;
    double L = 5.5;
    int n, *ia, *ja; 
    double *a, *b, *x, *gs_x;
    double tc1, tc2, tc3, tc4, tc5, tc6, tw1, tw2, tw3, tw4, tw5, tw6, tc7, tw7, tc8, tw8; // mis à jour le 13/10/22

    if (prob(m, &n, &ia, &ja, &a, &b, 1)) // source_value permet de lancer des simus de problèmes différents
        return 1;
    printf("\nPROBLEM: ");
    printf("m = %5d   n = %8d  nnz = %9d\n", m, n, ia[n]);

    // allouer la mémoire pour le vecteur de solution 

    x = calloc(1, n * sizeof(double));
    gs_x = calloc(1, n * sizeof(double));
    if ( x == NULL ) {
    printf("\n ERREUR : pas de mémoire pour vecteur des solutions\n\n");
        return 1;
    }

    // résoudre et mesurer le temps de solution 

    tc1 = mytimer_cpu(); tw1 = mytimer_wall(); // mis à jour le 13/10/22
    if ( solve_umfpack(n, ia, ja, a, b, x) ) {
        free(ia); free(ja); free(a); free(b); free(x); // empêche leak de mémoire en cas d'erreur
        return 1;
    }
    tc2 = mytimer_cpu(); tw2 = mytimer_wall(); // mis à jour le 13/10/22


    // calculer flux par la fenêtre, la porte et la puissance du radiateur

    // sauvegarder le vecteur solution pour faciliter la comparaison, principalement pour debug

    //FILE *f_out = fopen("mat/out_umfpack.dat", "w");

    //// créer le fichier de sortie pour gnuplot

    //int i = 0;
    ////avrg = 0.0; // construit pour donner la moyenne
    ////dim = 0; // taille de l'échantillon

    //for (int iy = 1; iy < m-1; iy++) { // vertical
    //    for (int ix = 1; ix < m-1; ix++) { // horizontal
    //        if ((iy < 6 * q || ix < 3 * q)) { // si on n'est pas dans le rectangle supérieur droit
    //            fprintf(f_out, "%f %f %f\n", iy * L / (q * 11), ix * L / (q * 11), (x)[i]);
    //            i++; // cycler à travers les éléments de x dans le même ordre qu'ils y ont été placés dans prob.c
    //        }
    //    }
    //    fprintf(f_out, "\n"); // requis par la syntaxe de gnuplot, ligne supplémentaire entre chaque changement de valeur de la 1e colonne (iy dans ce cas-ci)
    //}

    //fclose(f_out); // très important, sinon affichage incomplet de out.dat par gnuplot (optimisations compilateur n'attendaient pas l'écriture du fichier?)

    printf("\nTemps de solution, UMFPACK (CPU): %5.1f sec",tc2-tc1); // mis à jour le 13/10/22 
    printf("\nTemps de solution, UMFPACK (horloge): %5.1f sec \n",tw2-tw1); // mis à jour le 13/10/22 

/*    //tc3 = mytimer_cpu(); tw3 = mytimer_wall();

    double *r = malloc(n * sizeof(double)); // A * x pour pouvoir facilement faire A * x - b par la suite
    double *gs_r = malloc(n * sizeof(double)); // A * x pour pouvoir facilement faire A * x - b par la suite
    //double res_umfpack = residual(&n, &ia, &ja, &a, &b, &x, &r);
    double res_gs = residual(&n, &ia, &ja, &a, &b, &gs_x, &gs_r);
    printf("Initial residual is %.10e\n", res_gs);
    sym_gs(m, L, &n, &ia, &ja, &a, &b, &gs_x);
    sym_gs(m, L, &n, &ia, &ja, &a, &b, &gs_x);
    res_gs = residual(&n, &ia, &ja, &a, &b, &gs_x, &gs_r);
    printf("Pre-smoothing residual is %.10e\n", res_gs);
    int w = 0;

    //tc4 = mytimer_cpu(); tw4 = mytimer_wall();

    //printf("\nTemps de solution, Gauss-Seidel (CPU): %5.1f sec",tc4-tc3);
    //printf("\nTemps de solution, Gauss-Seidel (horloge): %5.1f sec \n",tw4-tw3);

    //printf("Initial residual\n\n\n");
    //for (int i = 0; i < n; i++) {
    //    printf("%f\n", gs_r[i]);
    //}
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
    if (prob(m_coarse, &n_coarse, &ia_coarse, &ja_coarse, &a_coarse, &b_coarse, rho_ptr, source_value, 0))
        return 1;
    //printf("\nPROBLEM: ");
    //printf("m = %5d   n = %8d  nnz = %9d\n", m_coarse, n_coarse, ia_coarse[n_coarse]);

    //r_coarse = malloc(n_coarse * sizeof(double));

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

    printf("Post-prolongation residual is %.10e\n", res_prolongation);

    //sym_gs(m, L, &n, &ia, &ja, &a, &b, &gs_x);
    //sym_gs(m, L, &n, &ia, &ja, &a, &b, &gs_x);    

    double two_grid_residual = residual(&n, &ia, &ja, &a, &b, &gs_x, &gs_r); 

    printf("Two-grid residual: %.10e\n", two_grid_residual);

*/

    
    /* initialisation des paramètres par défaut */
    //double Info [UMFPACK_INFO], Control [UMFPACK_CONTROL];
    //void *Symbolic, *Numeric ;
    //double *ia_coarse, *ja_coarse, *a_coarse, *b_coarse;
    //generate_coarse_problem(m, &ia_coarse, &ja_coarse, &a_coarse, &b_coarse, Symbolic, Numeric, Info, Control);

    double two_grid_residual = 1.0;

    double *gs_r = malloc(n * sizeof(double)); // A * x pour pouvoir facilement faire A * x - b par la suite

    int counter = 0;
    tc5 = mytimer_cpu(); tw5 = mytimer_wall();
    while (two_grid_residual > 2.5e-15) {
        //fwd_gs(m, L, &n, &ia, &ja, &a, &b, &gs_x);
        //two_grid_residual = residual(&n, &ia, &ja, &a, &b, &gs_x, &gs_r);
        two_grid_residual = two_grid_method(n, m, L, ia, ja, a, b, gs_x);
        //two_grid_residual = factorized_two_grid_method(n, m, L, ia, ja, a, b, gs_x, &ia_coarse, &ja_coarse, &a_coarse,
        //                                                    &b_coarse, Symbolic, Numeric, Info, Control);
        counter++;
        //if (counter == 3)
        //    sleep(50);
        //printf("GS residual: %.10e\n", two_grid_residual);
    }

    tc6 = mytimer_cpu(); tw6 = mytimer_wall();
    printf("\nSolution time, two-grid method (CPU): %5.1f sec",tc6-tc5);
    printf("\nSolution time, two-grid method (clock): %5.1f sec \n",tw6-tw5);
    printf("Number of iterations : %d\n", counter);

    //for (int i = 0; i < n; i++) {
    //    gs_x[i] = 0.0; // reset the solution to compare the 2 methods
    //}
    //double *gs_r = malloc(n * sizeof(double)); // A * x pour pouvoir facilement faire A * x - b par la suite

    //tc7 = mytimer_cpu(); tw7 = mytimer_wall();
    //while (residual(&n, &ia, &ja, &a, &b, &gs_x, &gs_r) > 5e-15) {
    //    sym_gs(m, L, &n, &ia, &ja, &a, &b, &gs_x);
    //}
    //tc8 = mytimer_cpu(); tw8 = mytimer_wall();
    //printf("\nSolution time, pure Gauss-Seidel (CPU): %5.1f sec",tc8-tc7);
    //printf("\nSolution time, pure Gauss-Seidel (clock): %5.1f sec \n",tw8-tw7);

    // créer le fichier de sortie pour gnuplot 

    FILE *f_out_two = fopen("mat/out_two_grid.dat", "w");

    i = 0;

    for (int iy = 1; iy < m-1; iy++) { // vertical
        for (int ix = 1; ix < m-1; ix++) { // horizontal
            if ((iy < 6 * q || ix < 3 * q)) { // si on n'est pas dans le rectangle supérieur droit
                fprintf(f_out_two, "%f %f %f\n", iy * L / (q * 11), ix * L / (q * 11), (gs_x)[i]);
                i++; // cycler à travers les éléments de x dans le même ordre qu'ils y ont été placés dans prob.c
            }
        }
        fprintf(f_out_two, "\n"); // requis par la syntaxe de gnuplot, ligne supplémentaire entre chaque changement de valeur de la 1e colonne (iy dans ce cas-ci)
    }

    fclose(f_out_two); // très important, sinon affichage incomplet de out.dat par gnuplot (optimisations compilateur n'attendaient pas l'écriture du fichier?)
    //sleep(1); // permet à gnuplot de lire le fichier avant de le supprimer

    int m_coarse = (m-1)/2 + 1;
    int q_coarse = (m_coarse-1) / 11; // nombre de fois que m-1 est multiple de 11
    int p_coarse = 4 * m_coarse - 4; // nombre de points sur le périmètre
    int n_coarse = m_coarse * m_coarse // nombre total de points dans le carré
            - (5 * q_coarse) * (8 * q_coarse) // nombre de points dans le rectangle supérieur droit
            - p_coarse; // number of points on the walls


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
    //system("gnuplot -persist \"heatmap.gnu\""); // laisser gnuplot afficher la température de la pièce
    //system("gnuplot -persist \"heatmap_two_grid.gnu\"");
    //system("gnuplot -persist \"heatmap_initial_residual.gnu\"");
    //system("gnuplot -persist \"heatmap_residual.gnu\""); // laisser gnuplot afficher la température de la pièce
    //system("gnuplot -persist \"heatmap_restriction.gnu\"");
    //system("gnuplot -persist \"heatmap_prolongation.gnu\"");

    //return 0;
}
