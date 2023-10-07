#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "prob.h"
#include "umfpk.h"
#include "time.h"
#include "residue.h"
#include "rho.h"
#include "heatflux.h"
#include "gs.h"


// Fonction main

int main(int argc, char *argv[])
{
    
    // déclarer les variables 

    double (*rho_ptr)(double, double, double) = &rho;
    int m = 166;
    int u = 0;
    int iter_max = 0; // régler le nombre d'itérations, mettre à 0 si on ne cherche pas à minimiser std_dev/avrg
    int q = (m-1) / 11;
    int use_petsc = 0; // activer ou déactiver l'utilisation de PETSc
    double source_value = 500.0 ; // permet d'itérer sur les différentes valeurs de rho pour 
    double flux_x, flux_y, rad_flux;
    double source_save;
    double save_dev;
    double avrg, std_dev;
    int dim;
    double L = 5.5;
    int n, *ia, *ja; 
    double *a, *b, *x, *petsc_x;
    double tc1, tc2, tc3, tc4, tc5, tc6, tw1, tw2, tw3, tw4, tw5, tw6; // mis à jour le 13/10/22 
    //double *vec_dev = malloc(1000 * sizeof(double)); // sauvegarder les 1000 résultats obtenus

    // itérer et trouver la meilleure valeur de rho pour résoudre le problème

    for (u=0; u < iter_max; u++) {

        // générer le problème

        if (prob(m, &n, &ia, &ja, &a, &b, rho_ptr, source_value, 0)) // source_value permet de lancer des simus de problèmes différents
            return 1;

        if (u == 0) { // pour ne pas imprimer le nombre d'inconnues 1000 fois
            printf("\nPROBLEM: ");
            printf("m = %5d   n = %8d  nnz = %9d\n", m, n, ia[n] );
        }

        // allouer la mémoire pour le vecteur de solution 

        x = malloc(n * sizeof(double));
        if ( x == NULL ) {
        free(ia); free(ja); free(a); free(b); free(x); // empêche leak de mémoire en cas d'erreur
        printf("\n ERREUR : pas de mémoire pour vecteur des solutions\n\n");
            return 1;
        }

        // résoudre et mesurer le temps de solution 

        if ( solve_umfpack(n, ia, ja, a, b, x) ) {
            free(ia); free(ja); free(a); free(b); free(x); // empêche leak de mémoire en cas d'erreur
            return 1;
        }
        
        // construire la moyenne

        int i = 0;
        avrg = 0.0; // construit pour donner la moyenne
        dim = 0; // taille de l'échantillon

        for (int iy = 0; iy < m; iy++) { // vertical
            for (int ix = 0; ix < m; ix++) { // horizontal
                if ((iy <= 6 * q || ix <= 3 * q) && !((iy == 0 && (q * 3 <= ix) && (ix <= q * 8)) || (ix == 0 && (q * 6 <= iy) && (iy <= q * 8)))) { // si on n'est pas dans le rectangle supérieur droit ou sur une porte / fenetre
                    avrg += (x)[i];
                    i++; // cycler à travers les éléments de x dans le même ordre qu'ils y ont été placés dans prob.c
                } else if ((iy == 0 && (q * 3 <= ix) && (ix <= q * 8))) { // il faut aussi représenter la fenêtre, or celle-ci ne fait pas partie des n inconnues -> rajouter à part
                    //avrg += 0.0; // perte de temps de calcul
                    dim += 1;
                } else if (ix == 0 && (q * 6 <= iy) && (iy <= q * 8)) { // il faut aussi représenter la porte, or celle-ci ne fait pas partie des n inconnues -> rajouter à part
                    avrg += 20.0;
                    dim += 1;
                }
            }
        }

        dim += i;
        avrg /= (double) dim; // on obtient la moyenne en la divisant par la taille de l'échantillon

        // calcul de l'écart type de l'échantillon

        std_dev = 0.0;

        i = 0;

        for (int iy = 0; iy < m; iy++) { // vertical
            for (int ix = 0; ix < m; ix++) { // horizontal
                if ((iy <= 6 * q || ix <= 3 * q) && !((iy == 0 && (q * 3 <= ix) && (ix <= q * 8)) || (ix == 0 && (q * 6 <= iy) && (iy <= q * 8)))) { // si on n'est pas dans le rectangle supérieur droit ou sur une porte / fenetre
                    std_dev += pow((x)[i], 2.0);
                    i++; // cycler à travers les éléments de x dans le même ordre qu'ils y ont été placés dans prob.c
                } else if ((iy == 0 && (q * 3 <= ix) && (ix <= q * 8))) { // il faut aussi représenter la fenêtre, or celle-ci ne fait pas partie des n inconnues -> rajouter à part
                    //std_dev += 0.0; // inutile
                } else if (ix == 0 && (q * 6 <= iy) && (iy <= q * 8)) { // il faut aussi représenter la porte, or celle-ci ne fait pas partie des n inconnues -> rajouter à part
                    std_dev += 400.0; //20^2, on écrit directement 400 pour épargner un calcul supplémentaire au processeur
                }
            }
        }
        
        std_dev /= (double) dim;
        std_dev -= pow(avrg, 2.0);
        std_dev = sqrt(std_dev); // l'écart-type est finalisé
        std_dev /= avrg; // on divise par la moyenne pour ne pas trop pénaliser les hautes températures de radiateur

        //vec_dev[u] = std_dev; // utile pour affichage de la courbe des écarts-types

        printf("\n%f    %f\n", source_value, std_dev);

        if (u == 0) { // première fois
            source_save = source_value;
            save_dev = std_dev;
        } else if (std_dev < save_dev) { // si on a trouvé un nouveau minimum
            source_save = source_value;
            save_dev = std_dev;
        }

        source_value += 1.0;

        free(ia); free(ja); free(a); free(b); free(x);

    }

    if (iter_max != 0) // si on a itéré
        source_value = source_save; // on assigne source_value à la meilleure solution
    // sinon on va résoudre avec source_value

    //printf("Meilleure valeur de rho: %f [K/m2]   Meilleur écart-type/moyenne: %f [-]\n", source_save, save_dev);

    // on a trouvé la meilleure temp du radiateur, on la re-résout pour pouvoir l'afficher ensuite (plus économe que de faire plein de fois de l'écriture de fichier)
    // alternativement, on aurait pu créer un array save_x dont on remplace les valeurs à chaque fois que l'écart-type/moyenne associé à x est meilleur

    if (prob(m, &n, &ia, &ja, &a, &b, rho_ptr, source_value, 1)) // source_value permet de lancer des simus de problèmes différents
        return 1;
    printf("\nPROBLEM: ");
    printf("m = %5d   n = %8d  nnz = %9d\n", m, n, ia[n]);

    // allouer la mémoire pour le vecteur de solution 

    x = calloc(1, n * sizeof(double));
    if ( x == NULL ) {
    printf("\n ERREUR : pas de mémoire pour vecteur des solutions\n\n");
        return 1;
    }

    // résoudre et mesurer le temps de solution 

    tc1 = mytimer_cpu(); tw1 = mytimer_wall(); // mis à jour le 13/10/22
    if ( solve_umfpack(n, ia, ja, a, b, x) ) {
        free(ia); free(ja); free(a); free(b);free(x); // empêche leak de mémoire en cas d'erreur
        return 1;
    }
    tc2 = mytimer_cpu(); tw2 = mytimer_wall(); // mis à jour le 13/10/22

    // calculer flux par la fenêtre, la porte et la puissance du radiateur

    //compute_heat_flux(m, q, &x, source_save, &flux_x, &flux_y, &rad_flux, rho_ptr); // void qui modifie directement les variables de flux

    //printf("\nValeur des 2 flux: (%f, %f) [W]\n", flux_x, flux_y);
    //printf("Puissance du radiateur et flux entrant/sortant: %f, %f\n", rad_flux, fabs(flux_x + flux_y));
    //printf("Proportion puissance/flux: %f\n", rad_flux/fabs(flux_x + flux_y));

    // sauvegarder le vecteur solution pour faciliter la comparaison, principalement pour debug

    
    FILE *f_out = fopen("mat/out_umfpack.dat", "w");

    // créer le fichier de sortie pour gnuplot

    int i = 0;
    avrg = 0.0; // construit pour donner la moyenne
    dim = 0; // taille de l'échantillon

    for (int iy = 1; iy < m-1; iy++) { // vertical
        for (int ix = 1; ix < m-1; ix++) { // horizontal
            if ((iy < 6 * q || ix < 3 * q)) { // si on n'est pas dans le rectangle supérieur droit ou sur une porte / fenetre
                fprintf(f_out, "%f %f %f\n", iy * L / (q * 11), ix * L / (q * 11), (x)[i]);
                i++; // cycler à travers les éléments de x dans le même ordre qu'ils y ont été placés dans prob.c
            }
        }
        fprintf(f_out, "\n"); // requis par la syntaxe de gnuplot, ligne supplémentaire entre chaque changement de valeur de la 1e colonne (iy dans ce cas-ci)
    }

    fclose(f_out); // très important, sinon affichage incomplet de out.dat par gnuplot (optimisations compilateur n'attendaient pas l'écriture du fichier?)


    printf("\nTemps de solution, UMFPACK (CPU): %5.1f sec",tc2-tc1); // mis à jour le 13/10/22 
    printf("\nTemps de solution, UMFPACK (horloge): %5.1f sec \n",tw2-tw1); // mis à jour le 13/10/22 

    tc3 = mytimer_cpu(); tw3 = mytimer_wall();

    double res_umfpack = residue(&n, &ia, &ja, &a, &b, &x);
    printf("\nRésidu de la solution: %.10e\n", res_umfpack);

    tc4 = mytimer_cpu(); tw4 = mytimer_wall();

    printf("\nTemps de calcul du résidu UMFPACK (CPU): %5.1f sec",tc4-tc3);
    printf("\nTemps de calcul du résidu UMFPACK (horloge): %5.1f sec \n",tw4-tw3);

    // comparer avec les routines PETSc

    /*if (use_petsc) { // si on veut utiliser PETSc

        petsc_x = malloc(n * sizeof(double));
        solve_petsc(n, ia, ja, a, b, &petsc_x); // construire les vecteurs/matrices + résoudre
        

        // créer le fichier de sortie pour gnuplot

        FILE *f_out_petsc = fopen("mat/out_petsc.dat", "w");
        i = 0;
        avrg = 0.0; // construit pour donner la moyenne
        dim = 0; // taille de l'échantillon

        for (int iy = 0; iy < m; iy++) { // vertical
            for (int ix = 0; ix < m; ix++) { // horizontal
                if ((iy <= 6 * q || ix <= 3 * q) && !((iy == 0 && (q * 3 <= ix) && (ix <= q * 8)) || (ix == 0 && (q * 6 <= iy) && (iy <= q * 8)))) { // si on n'est pas dans le rectangle supérieur droit ou sur une porte / fenetre
                    fprintf(f_out_petsc, "%f %f %f\n", iy * L / (q * 11), ix * L / (q * 11), (petsc_x)[i]);
                    i++; // cycler à travers les éléments de x dans le même ordre qu'ils y ont été placés dans prob.c
                } else if ((iy == 0 && (q * 3 <= ix) && (ix <= q * 8))) { // il faut aussi représenter la fenêtre, or celle-ci ne fait pas partie des n inconnues -> rajouter à part
                    fprintf(f_out_petsc, "%f %f %f\n", iy * L / (q * 11), ix * L / (q * 11), 0.0); 
                } else if (ix == 0 && (q * 6 <= iy) && (iy <= q * 8)) { // il faut aussi représenter la porte, or celle-ci ne fait pas partie des n inconnues -> rajouter à part
                    fprintf(f_out_petsc, "%f %f %f\n", iy * L / (q * 11), ix * L / (q * 11), 20.0);
                }
            }
            fprintf(f_out_petsc, "\n"); // requis par la syntaxe de gnuplot, ligne supplémentaire entre chaque changement de valeur de la 1e colonne (iy dans ce cas-ci)
        }

        fclose(f_out_petsc); // très important, sinon affichage incomplet de out.dat par gnuplot (optimisations compilateur n'attendaient pas l'écriture du fichier?)


        tc5 = mytimer_cpu(); tw5 = mytimer_wall();

        double r_petsc = residue(&n, &ia, &ja, &a, &b, &petsc_x);
        printf("\nRésidu de la solution PETSc: %.10e\n", r_petsc);

        tc6 = mytimer_cpu(); tw6 = mytimer_wall();

        printf("\nTemps de calcul du résidu PETSc (CPU): %5.1f sec",tc4-tc3);
        printf("\nTemps de calcul du résidu PETSc (horloge): %5.1f sec \n",tw4-tw3);

        system("gnuplot -persist \"heatmap_petsc_solution.gnu\""); // laisser gnuplot afficher la température de la pièce

        //free(petsc_x); // bug si déallouage de petsc_x :(

    }*/

    free(ia); free(ja); free(a); free(b); free(x); //free(vec_dev); 
    system("gnuplot -persist \"heatmap.gnu\""); // laisser gnuplot afficher la température de la pièce
    
    return 0;
}
