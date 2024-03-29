#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int prob(int m, int *n, int **ia, int **ja, double **a, double **b, double (*source_func)(double, double, double), double source_value, int write) // on rajoute la fonction source en argument
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
    int  nnz, ix, iy, ind, q, p, nb_dir;
    double invh2, L = 5.5, Tp = 20.0;

    if(m <= 1 || (m - 1) % 11 != 0) { // multiple de 11 pour éviter des descriptions non idéales des CB
        printf("\n ERREUR : m = %d n'est pas une valeur valide\n\n",m);
        return 1;
    }
    q = (m-1) / 11; // nombre de fois que m-1 est multiple de 11
    p = 4 * m - 4; // nombre de points sur le périmètre
    nb_dir = (2 * q + 1) + (5 * q + 1); // nombre de points sur une porte/fenêtre
    invh2 = (m-1)*(m-1) / (L*L); /* h^-2 */
    //*n  = (m-1) * m; /* nombre d'inconnues */
    *n = m * m // nombre total de points dans le carré
        - (5 * q) * (8 * q) + // nombre de points dans le rectangle supérieur droit
        - nb_dir; // nombre de points sur la porte et la fenêtre
    //nnz = 5 * (m-1) * m - 4 * m+2; /* nombre d'éléments non nuls */
    nnz = 5 * (m * m - (5 * q) * (8 * q) - p) // points soumis à l'équation de la chaleur classique sans conditions particulières (on retire le rectangle supérieur droit)
        - 5 * nb_dir // il faut aussi retirer les nb_dir points directement en face des portes/fenetres de la somme ci-dessus
        + 4 * nb_dir // points directement en face d'une porte/fenêtre
        + 4 * (p - nb_dir) - 4 * 6  - 4 * 4 // nombre de points soumis à une condition de neumann simple (murs sauf coins et points adjacents porte/fenêtre)
        + 3 * 4 // points du mur adjacents à une porte/fenêtre
        + 3 * 5 // nombre de points soumis à une condition de neumann double (les coins de la pièce)
        + 5 * 1; // coin "du milieu", pas de normale donc on n'applique pas les conditions de neumann
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

    // remplissage de la matrice

    for (iy = 0; iy < m; iy++) { // vertical
        for (ix = 0; ix < m; ix++) { // horizontal
            if ((iy <= 6 * q || ix <= 3 * q) && !((iy == 0 && (q * 3 <= ix) && (ix <= q * 8)) || (ix == 0 && (q * 6 <= iy) && (iy <= q * 8)))) { // si on n'est pas dans le rectangle supérieur droit ou sur une porte / fenetre
                /* marquer le début de la ligne suivante dans le tableau 'ia' */
                (*ia)[ind] = nnz;

                /* calculer le membre de droite */
                (*b)[ind] = source_func(ix * L / (q * 11), iy * L / (q * 11), source_value); /* rho = temp locale */

                /* remplissage de la ligne : voisin sud */
                    if (iy == m-1 || (iy == q * 6 && ix > q * 3)) { /* condition de Neumann, bord nord + attention au rectangle retiré */
                                                        // strictement supérieur implique que le coin 'intérieur' sans normale n'est pas traité avec une condition de Neumann
                        (*a)[nnz] = -2*invh2;
                        if (iy == m-1)
                            (*ja)[nnz] = ind - 3*q - 1; // sinon il y a m points dans la liste avant d'arriver au point spatialement en dessous
                        else
                            (*ja)[nnz] = ind - m + 1; // +1 à cause de la porte
                        nnz++;
                    } else if (iy > 0) { /* milieu du domaine */
                        if (iy == q * 8 + 1 && ix == 0) { // JUSTE AU DESSUS DE LA PORTE
                            (*b)[ind] += Tp*invh2;
                        } else if (!(iy == 1 && (q * 3 <= ix) && (ix <= q * 8))) { // IL NE FAUT RIEN FAIRE AU DESSUS DE LA FENETRE PUISQUE TF = 0
                            (*a)[nnz] = -invh2;
                            if (iy == q * 6 || iy == q * 6 + 1) {
                                (*ja)[nnz] = ind - m + 1; // dirichlet de la porte retire une inconnue
                            } else if ((q * 6 < iy) && (iy <= q * 8)) { // en face de la porte
                                (*ja)[nnz] = ind - q * 3; // il y a moins d'inconnues sur cette bande => il faut moins diminuer les indices
                            } else if (q * 8 < iy) {
                                (*ja)[nnz] = ind - q * 3 - 1; // diminuer d'un de plus car la condition de Dirichlet ne fait pas perdre une inconnue
                            } else if (iy == 1 && ix < q * 3) {
                                (*ja)[nnz] = ind - 6 * q; // la fenêtre diminue le nombre de points sur lesquels il faut retourner
                            } else {
                                (*ja)[nnz] = ind - m; // sinon il y a m points dans la liste avant d'arriver au point spatialement en dessous
                            }
                            nnz++;
                        }
                    }

                /* remplissage de la ligne : voisin ouest */
                    if (ix == m-1 || ((iy > q * 6 && ix == q * 3))) { /* condition de Neumann, bord est + attention au rectangle retiré*/   
                        (*a)[nnz] = -2*invh2; 
                        (*ja)[nnz] = ind - 1;
                        nnz++;
                    } else if (ix == q * 8 + 1 && iy == 0) { // juste à droite de la fenêtre mais il ne faut rien faire puisque Tf = 0
                        (*b)[ind] += 0*invh2;
                    } else if (ix > 0) { /* milieu du domaine */
                        if ((ix == 1) && (q * 6 <= iy) && (iy <= q * 8)) { // si on est juste à droite de la porte
                            (*b)[ind] += Tp*invh2;
                        } else {
                        (*a)[nnz] = -invh2; 
                        (*ja)[nnz] = ind - 1;
                        nnz++;
                        }
                    }

                /* remplissage de la ligne : élément diagonal */
                    (*a)[nnz] = 4.0*invh2;
                    (*ja)[nnz] = ind;
                    nnz++;

                /* remplissage de la ligne : voisin est */
                // IL NE FAUT RIEN FAIRE A COTE DE LA FENETRE PUISQUE TF = 0
                    if (ix == m - 1 || (ix == q * 3 && iy > q * 6)) { // inégalité stricte sur iy pour ne pas prendre en compte le coin spécial
                        // ne rien faire car pas de voisin est
                    } else if (ix == 0) { // Neumann sur bord ouest
                        (*a)[nnz] = -2*invh2;
                        (*ja)[nnz] = ind + 1;
                        nnz++;
                    } else if (ix < m-1) { /* milieu du domaine */
                        if (!(iy == 0 && ix == 3 * q - 1)) { // on saute le point juste à gauche de la fenêtre, puisque Tf = 0 Dirichlet ne fait rien
                            (*a)[nnz] = -invh2;
                            (*ja)[nnz] = ind + 1;
                            nnz++;
                        }
                    }

                /* remplissage de la ligne : voisin nord */
                    if (iy == 0) { /* condition de Neumann, bord sud */
                        (*a)[nnz] = -2*invh2; 
                        if (iy == 0 && ix < q * 3) {
                            (*ja)[nnz] = ind + q * 6;
                        } else {
                            (*ja)[nnz] = ind + m; // sinon il y a m points dans la liste avant d'arriver au point spatialement au dessus
                        }
                        nnz++;
                    } else if (ix == 0 && iy == q * 6 - 1) { // juste en dessous de la porte
                        (*b)[ind] += Tp*invh2;
                    } else if (iy < m - 1 && !((iy == q*6) && (q*3 < ix))) { /* milieu du domaine */
                        (*a)[nnz] = -invh2;

                        if ((q * 6 < iy) && (iy < q * 8)) {
                            (*ja)[nnz] = ind + q * 3; // il y a moins d'inconnues sur cette bande => il faut moins augmenter les indices
                        } else if (q * 8 <= iy) {
                            (*ja)[nnz] = ind + q * 3 + 1; // augmenter d'un de plus car la condition de Dirichlet ne fait pas perdre une inconnue
                        } else if (iy == q * 6 - 1 || iy == q * 6) {
                            (*ja)[nnz] = ind + m - 1;
                        } else {
                            (*ja)[nnz] = ind + m; // sinon il y a m points dans la liste avant d'arriver au point spatialement au dessus
                        }
                        nnz++;
                    }
                // prochaine equation
                ind++;
            }
        }
    }

    //printf("%d", nnz);

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
