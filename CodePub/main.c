#include <stdio.h>
#include <stdlib.h>
#include "prob.h"
#include "umfpk.h"
#include "time.h"
#include "residue.h"


/* Fonction main */

int main(int argc, char *argv[])
{

    /* déclarer les variables */

    int m = 661;
    int q = (m-1) / 11;
    double L = 5.5;
    int n, *ia, *ja; 
    double *a, *b, *x;
    double tc1, tc2, tc3, tc4, tw1, tw2, tw3, tw4; /* mis à jour le 13/10/22 */

    /* générér le problème */

    if (prob(m, &n, &ia, &ja, &a, &b)) // pas oublier de rajouter rho
        return 1;
    printf("\nPROBLEM: ");
    printf("m = %5d   n = %8d  nnz = %9d\n", m, n, ia[n] );


    /* allouer la mémoire pour le vecteur de solution */

    x = malloc(n * sizeof(double));
    if ( x == NULL ) {
    printf("\n ERREUR : pas de mémoire pour vecteur des solutions\n\n");
        return 1;
    }

    /* résoudre et mesurer le temps de solution */

    tc1 = mytimer_cpu(); tw1 = mytimer_wall(); /* mis à jour le 13/10/22 */
    if ( solve_umfpack(n, ia, ja, a, b, x) ) {
        free(ia); free(ja); free(a); free(b);free(x); // empêche leak de mémoire en cas d'erreur
        return 1;
    }
    tc2 = mytimer_cpu(); tw2 = mytimer_wall(); /* mis à jour le 13/10/22 */

    // sauvegarder le vecteur solution pour faciliter la comparaison, principalement pour debug

    FILE *f_x = fopen("mat/x.txt", "w"); // debug
    FILE *f_out = fopen("mat/out.dat", "w");

    for (int i = 0; i < n + 1; i++) {
        fprintf(f_x, "%f\n", (x)[i]); // debug
    }

    // créer le fichier de sortie pour gnuplot

    int i = 0;

    for (int iy = 0; iy < m; iy++) { // vertical
        for (int ix = 0; ix < m; ix++) { // horizontal
            if ((iy <= 6 * q || ix <= 3 * q) && !((iy == 0 && (q * 3 <= ix) && (ix <= q * 8)) || (ix == 0 && (q * 6 <= iy) && (iy <= q * 8)))) { // si on n'est pas dans le rectangle supérieur droit ou sur une porte / fenetre
                fprintf(f_out, "%f %f %f\n", iy * L / (q * 11), ix * L / (q * 11), (x)[i]);
                i++; // cycler à travers les éléments de x dans le même ordre qu'ils y ont été placés dans prob.c
            } else if ((iy == 0 && (q * 3 <= ix) && (ix <= q * 8))) { // il faut aussi représenter la fenêtre, or celle-ci ne fait pas partie des n inconnues -> rajouter à part
                fprintf(f_out, "%f %f %f\n", iy * L / (q * 11), ix * L / (q * 11), 0.0); 
            } else if (ix == 0 && (q * 6 <= iy) && (iy <= q * 8)) { // il faut aussi représenter la porte, or celle-ci ne fait pas partie des n inconnues -> rajouter à part
                fprintf(f_out, "%f %f %f\n", iy * L / (q * 11), ix * L / (q * 11), 20.0);
            }
        }
        fprintf(f_out, "\n"); // requis par la syntaxe de gnuplot, ligne supplémentaire entre chaque changement de valeur de la 1e colonne (iy dans ce cas-ci)
    }

    tc3 = mytimer_cpu(); tw3 = mytimer_wall();

    double r = residue(&n, &ia, &ja, &a, &b, &x);
    printf("\nRésidu de la solution: %.10e\n", r);

    tc4 = mytimer_cpu(); tw4 = mytimer_wall();


    printf("\nTemps de solution (CPU): %5.1f sec",tc2-tc1); /* mis à jour le 13/10/22 */
    printf("\nTemps de solution (horloge): %5.1f sec \n",tw2-tw1); /* mis à jour le 13/10/22 */
    printf("\nTemps de calcul du résidu (CPU): %5.1f sec",tc4-tc3);
    printf("\nTemps de calcul du résidu (horloge): %5.1f sec \n",tw4-tw3);

    /* libérér la mémoire */

    free(ia); free(ja); free(a); free(b); free(x);
    system("gnuplot -persist \"heatmap.gnu\""); // laisser gnuplot afficher la température de la pièce
    return 0;
}

