#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int prob(int m, int *n, int **ia, int **ja, double **a, double **b)
/*
   But
   ===
   Générer le système linéaire n x n 
                          
                             Au = b                                   

   qui correspond à la discrétisation sur une grille cartésienne 
   régulière m x m de l'équation de Poisson à deux dimensions
              
            d    d        d    d
         - == ( == u ) - == ( == u )  = 0     sur [0,1] x [0,1]
           dx   dx       dy   dy

  avec les conditions aux limites de Dirichlet
         
         u = 20     sur (0,y)                 , avec 0 <=  y  <= 1
         du/dn = 1  sur (1,y), (x,0) et (x,1) , avec 0 <= x,y <= 1 .

  La numérotation des inconnues est lexicographique, la direction x étant 
  parcourue avant celle de y. La matrice A est retournée dans le format 
  CRS qui est défini via trois tableaux : 'ia', 'ja' et 'a'.

  Arguments
  =========
  m  (input)  - nombre de points par direction dans la grille
                (les valeurs m inférieures à 2 ne sont pas valides) 
  n  (output) - pointeur vers le nombre d'inconnus dans le système
  ia (output) - pointeur vers le tableau 'ia' de la matrice A
  ja (output) - pointeur vers le tableau 'ja' de la matrice A
  a  (output) - pointeur vers le tableau 'a' de la matrice A
  b  (output) - pointeur vers le tableau 'b'

*/
{
    int  nnz, ix, iy, ind;
    double invh2, Tp = 20.0;

    if(m <= 1) {
        printf("\n ERREUR : m = %d n'est pas une valeur valide\n\n",m);
        return 1;
    }
    invh2 = (m-1)*(m-1); /* h^-2 pour L=1 */
    *n  = (m-1) * m; /* nombre d'inconnues */
    nnz = 5 * (m-1) * m - 4 * m+2; /* nombre d'éléments non nuls */

    /* allocation des tableaux */

    *ia  = malloc((*n + 1) * sizeof(int));
    *ja  = malloc(nnz * sizeof(int));
    *a   = malloc(nnz * sizeof(double));
    *b   = malloc(*n * sizeof(double));

    /* allocation réussite? */

    if (*ia == NULL || *ja == NULL || *a == NULL || *b == NULL ) {
        printf("\n ERREUR : pas assez de mémoire pour générer le système\n\n");
        return 1;
    }

    /* partie principale : remplissage de la matrice */

    ind = 0; /* au cas où m<=1 */
    nnz = 0;
    for (iy = 0; iy < m; iy++) {
        for (ix = 0; ix < m-1; ix++) {
            /* numéro de l'équation */
            ind = ix + (m-1) * iy;
            /* marquer le début de la ligne suivante dans le tableau 'ia' */
            (*ia)[ind] = nnz;

            /* calculer le membre de droite */
            (*b)[ind] = 0.0; /* rho = 0 */
   
            /* remplissage de la ligne : voisin sud */
            if (iy == m-1) { /* condition de Neumann, bord nord */              
                (*a)[nnz] = -2*invh2;
                (*ja)[nnz] = ind - m+1;
                nnz++;
            } else if (iy > 0) { /* milieu du domaine */
                (*a)[nnz] = -invh2; 
                (*ja)[nnz] = ind - m+1;
                nnz++;
            } 

            /* remplissage de la ligne : voisin ouest */
            if (ix == m-2) { /* condition de Neumann, bord est */   
                (*a)[nnz] = -2*invh2; 
                (*ja)[nnz] = ind - 1;
                nnz++;
            } else if (ix > 0) { /* milieu du domaine */
                (*a)[nnz] = -invh2; 
                (*ja)[nnz] = ind - 1;
                nnz++;
            } else if (ix == 0) { /* condition de Dirichlet, bord est */
                (*b)[ind] += Tp*invh2; 
            }

            /* remplissage de la ligne : élément diagonal */
            (*a)[nnz] = 4.0*invh2;
            (*ja)[nnz] = ind;
            nnz++;

            /* remplissage de la ligne : voisin est */
            if (ix < m-2) { /* milieu du domaine */
                (*a)[nnz] = -invh2;
                (*ja)[nnz] = ind + 1;
                nnz++;
            } 

            /* remplissage de la ligne : voisin nord */
            if (iy == 0) { /* condition de Neumann, bord sud */         
                (*a)[nnz] = -2*invh2; 
                (*ja)[nnz] = ind + m-1;
                nnz++;
            } else if (iy < m-1) { /* milieu du domaine */
                (*a)[nnz] = -invh2;
                (*ja)[nnz] = ind + m-1;
                nnz++;
            }
        }
    }

    /* dernier élément du tableau 'ia' */
    (*ia)[ind + 1] = nnz;

    /* retour habituel de fonction */
    return 0;
}
