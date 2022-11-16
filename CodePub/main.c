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

    int m = 1112;
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

    FILE *f_x = fopen("mat/x.txt", "w"); 

    for (int i = 0; i < n + 1; i++) {
        fprintf(f_x, "%f\n", (x)[i]);
    }

    tc3 = mytimer_cpu(); tw3 = mytimer_wall();

    double r = residue(&n, &ia, &ja, &a, &b, &x);
    printf("\nRésidu de la solution: %.10e\n", r);

    tc4 = mytimer_cpu(); tw4 = mytimer_wall();


    printf("\nTemps de solution (CPU): %5.1f sec",tc2-tc1); /* mis à jour le 13/10/22 */
    printf("\nTemps de solution (horloge): %5.1f sec \n",tw2-tw1); /* mis à jour le 13/10/22 */
    printf("\nTemps de calcul du résidu (CPU): %5.1f sec",tc4-tc3); /* mis à jour le 13/10/22 */
    printf("\nTemps de calcul du résidu (horloge): %5.1f sec \n",tw4-tw3); /* mis à jour le 13/10/22 */
    

    /* libérér la mémoire */

    free(ia); free(ja); free(a); free(b); free(x);
    return 0;
}

