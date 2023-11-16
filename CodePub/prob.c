#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int prob(int m, int *n, int **ia, int **ja, double **a, double **b, double (*source_func)(double, double, double),
         double source_value, int write) // on rajoute la fonction source en argument
/*
   But
   ===
   Générer le système linéaire n x n 
                          
                             Au = b 

   qui correspond à la discrétisation sur une grille cartésienne 
   régulière m x m de l'équation de Poisson à deux dimensions
              
            d    d        d    d
         - == ( == u ) - == ( == u )  = rho     sur [0,5.5] x [0,5.5]
           dx   dx       dy   dy

  avec les conditions aux limites de Dirichlet
         
         u = 20     sur (0,y)                 , avec 3 <=  y  <= 4 [m]
         u = 0      sur (x, 0)                , avec 1.5 <= x <= 4 [m]
         du/dn = 1  sur (1,y), (x,0) et (x,1) , avec 0 <= x,y <= 1 .

  La numérotation des inconnues est lexicographique, la direction x étant 
  parcourue avant celle de y. La matrice A est retournée dans le format 
  CRS qui est défini via trois tableaux : 'ia', 'ja' et 'a'.

  Arguments
  =========
  m  (input)  - nombre de points par direction dans la grille
                (les valeurs m inférieures à 2 ne sont pas valides) 
  n  (output) - pointeur vers le nombre d'inconnues dans le système
  ia (output) - pointeur vers le tableau 'ia' de la matrice A
  ja (output) - pointeur vers le tableau 'ja' de la matrice A
  a  (output) - pointeur vers le tableau 'a' de la matrice A
  b  (output) - pointeur vers le tableau 'b'
  rho - fonction source

*/

{
    int  nnz, ix, iy, ind, q, p, nb_next_to_wall, nb_dir;
    double invh2, L = 5.5, Tp = 20.0;

    if(m <= 1 || (m - 1) % 11 != 0) { // multiple de 11 pour éviter des descriptions non idéales des CB
        printf("\n ERREUR : m = %d n'est pas une valeur valide\n\n",m);
        return 1;
    }
    q = (m-1) / 11; // nombre de fois que m-1 est multiple de 11
    p = 4 * m - 4; // nombre de points sur le périmètre
    nb_next_to_wall = 4 * (m-2); // nombre de points adjacents à un mur
    nb_dir = p; // nombre de points sur une porte/fenêtre
    invh2 = (m-1)*(m-1) / (L*L); /* h^-2 */
    //*n  = (m-1) * m; /* nombre d'inconnues */
    *n = m * m // nombre total de points dans le carré
        - (5 * q) * (8 * q) // nombre de points dans le rectangle supérieur droit
        - nb_dir; // number of points on the walls
    //*n = 1;
    //printf("Hello\n");
    //printf("Value of n: %d\n", *n); // \n is necessary, otherwise the buffer is not flushed
    //nnz = 5 * (m-1) * m - 4 * m+2; /* nombre d'éléments non nuls */
    nnz = 5 * (m * m - (5 * q) * (8 * q) - p) // points soumis à l'équation de la chaleur classique sans conditions particulières (on retire le rectangle supérieur droit)
        - 5 * nb_next_to_wall // we also need to remove the points directly in front of a wall (because of the dirichlet condition)
        + 4 * nb_next_to_wall // points directly in front of a wall
        //- 4 * 5 // points in the corners of the room
        //+ 3 * 5 // points in the corners of the room
        //+ 4 * (p - nb_dir) - 4 * 6  - 4 * 4 // nombre de points soumis à une condition de neumann simple (murs sauf coins et points adjacents porte/fenêtre)
        //+ 3 * 4 // points du mur adjacents à une porte/fenêtre
        //+ 3 * 5 // nombre de points soumis à une condition de neumann double (les coins de la pièce)
        //+ 5 * 1; // coin "du milieu", pas de normale donc on n'applique pas les conditions de neumann
        ;

    //nnz = 5 * *n - 1 * 4 * (m - 2); // nombre d'éléments non nuls
    printf("Value of nnz: %d\n", nnz);
    /* allocation des tableaux */

    *ia  = malloc((*n + 1) * sizeof(int));
    *ja  = malloc(nnz * sizeof(int));
    *a   = malloc(nnz * sizeof(double));
    *b   = malloc(*n * sizeof(double));

    /* allocation réussie? */

    if (*ia == NULL || *ja == NULL || *a == NULL || *b == NULL ) {
        printf("\n ERREUR : pas assez de mémoire pour générer le système\n\n");
        return 1;
    }

    /* partie principale : remplissage de la matrice */

    ind = 0; /* au cas où m<=1 */
    nnz = 0;

    // remplissage de la matrice
    // m-1 becquse indices start at 0
    for (iy = 1; iy < m-1; iy++) { // vertical 
        for (ix = 1; ix < m-1; ix++) { // horizontal
            if ((iy < 6 * q || ix < 3 * q)) { // si on n'est pas dans le rectangle supérieur droit
                /* marquer le début de la ligne suivante dans le tableau 'ia' */
                (*ia)[ind] = nnz;

                /* calculer le membre de droite */
                (*b)[ind] = 0; /* rho = temp locale */

                /* remplissage de la ligne : voisin sud */
                    //printf("Hello this is south, ind = %d, nnz = %d\n", ind, nnz);
                    if (iy == 1) { // RIGHT ABOVE THE BOTTOM WALL
                        (*b)[ind] += sin(sqrt(((ix * ix + (iy-1) * (iy-1))/invh2)))*invh2; // -1 because it is the coordinates of the wall point
                                                                                           // and not of our room point we need
                            // no need for a ja modification because this directly goes into the b vector, not the matrix
                    } else {
                        (*a)[nnz] = -invh2;
                        if ((q * 6 < iy)) { // en face de la porte
                            (*ja)[nnz] = ind - q * 3 + 1; // il y a moins d'inconnues sur cette bande => il faut moins diminuer les indices
                        } else {
                            (*ja)[nnz] = ind - m + 2; // +2 because there are less points to go back over now that all the walls are dirichlet
                        }
                        nnz++;
                    }

                /* remplissage de la ligne : voisin ouest */ 
                        //printf("Hello this is west, ind = %d, nnz = %d\n", ind, nnz);
                    if ((ix == 1)) { // si on est juste à droite de la porte
                        (*b)[ind] += sin(sqrt((((ix-1) * (ix-1) + iy* iy)/invh2)))*invh2;
                    } else {
                        (*a)[nnz] = -invh2; 
                        (*ja)[nnz] = ind - 1;
                        nnz++;
                    }
                    
                /* remplissage de la ligne : élément diagonal */
                    //printf("Hello this is diag, ind = %d, nnz = %d\n", ind, nnz);
                    (*a)[nnz] = 4.0*invh2;
                    (*ja)[nnz] = ind;
                    nnz++;

                /* remplissage de la ligne : voisin est */
                    //printf("Hello this is east, ind = %d, nnz = %d\n", ind, nnz);
                    if (ix == m - 2 || (ix == 3 * q - 1 && iy >= 6 * q)) { // if we are next to the right wall
                        (*b)[ind] += sin(sqrt((((ix+1) * (ix+1) + iy* iy)/invh2)))*invh2;
                    } else { // on saute le point juste à gauche de la fenêtre, puisque Tf = 0 Dirichlet ne fait rien
                        (*a)[nnz] = -invh2;
                        (*ja)[nnz] = ind + 1;
                        nnz++;
                    }

                /* remplissage de la ligne : voisin nord */
                    //printf("Hello this is north, ind = %d, nnz = %d\n", ind, nnz);
                    if (iy == m-2 || (iy == q * 6 - 1 && ix >= 3 * q)) { // juste en dessous de la porte
                        (*b)[ind] += sin(sqrt(((ix * ix + (iy-1) * (iy-1))/invh2)))*invh2; // -1 because it is the coordinates of the wall point
                    } else { /* milieu du domaine */
                        (*a)[nnz] = -invh2;
                        if ((q * 6 <= iy)) {
                            (*ja)[nnz] = ind + q * 3 - 1; // il y a moins d'inconnues sur cette bande => il faut moins augmenter les indices
                        } else {
                            (*ja)[nnz] = ind + m - 2; // sinon il y a m-2 points dans la liste avant d'arriver au point spatialement au dessus
                        }
                        nnz++;
                    }
                // prochaine equation
                ind++;
            }
        }
    }
    

    printf("%d\n", nnz);

    /* dernier élément du tableau 'ia' */
    (*ia)[ind] = nnz;


    if (write) { // ne pas créer de fichiers quand on veut lancer beaucoup de simulations en un coup -> l'io sur des fichiers est lent

        // création des fichiers

        FILE *f_ia = fopen("mat/ia.txt", "w");
        FILE *f_ja = fopen("mat/ja.txt", "w");
        FILE *f_a = fopen("mat/a.txt", "w");
        FILE *f_b = fopen("mat/b.txt", "w");

        // écriture dans les fichiers respectifs

        for (int i = 0; i < *n + 1; i++) {
            fprintf(f_ia, "%d\n", (*ia)[i]);
        }

        for (int i = 0; i < nnz; i++) {
            fprintf(f_ja, "%d\n", (*ja)[i]);
        }

        for (int i = 0; i < nnz; i++) {
            fprintf(f_a, "%f\n", (*a)[i]);
        }

        for (int i = 0; i < *n; i++) {
            fprintf(f_b, "%f\n", (*b)[i]);
        }

        fclose(f_ia); fclose(f_ja); fclose(f_a); fclose(f_b); // fermer les fichiers proprement

    }

    /* retour habituel de fonction */
    return 0;
}
