#include <stdio.h>
#include <stdlib.h>
#include "prob.h"
#include "umfpk.h"
#include "time.h"


/* Fonction main */

int main(int argc, char *argv[])
{

  /* déclarer les variables */

  int m = 1211;
  int n, *ia, *ja; 
  double *a, *b, *x;
  double tc1, tc2, tw1, tw2; /* mis à jour le 13/10/22 */

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
  printf("\nTemps de solution (CPU): %5.1f sec",tc2-tc1); /* mis à jour le 13/10/22 */
  printf("\nTemps de solution (horloge): %5.1f sec \n",tw2-tw1); /* mis à jour le 13/10/22 */

  /* libérér la mémoire */

  free(ia); free(ja); free(a); free(b); free(x);
  return 0;
}

