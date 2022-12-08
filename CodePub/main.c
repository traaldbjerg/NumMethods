#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "prob.h"
#include "umfpk.h"
#include "time.h"
#include "residue.h"
#include "rho.h"
#include "heatflux.h"

// Fonction main 

int main(int argc, char *argv[])
{
    
    // déclarer les variables 

    double (*rho_ptr)(double, double, double) = &rho;
    int m = 166;
    int u = 0;
    int q = (m-1) / 11;
    //int nb_dir = (2 * q + 1) + (5 * q + 1); // nombre de points sur une porte/fenêtre
    double temp_rad = 0.00; // permet d'itérer sur les différentes valeurs de rho pour 
    double flux_x, flux_y, rad_flux;
    double save_temp;
    double save_dev;
    double avrg, std_dev;
    int dim;
    double L = 5.5;
    int n, *ia, *ja; 
    double *a, *b, *x;
    //double h = L / (m-1); // longueur d'un pas
    double tc1, tc2, tc3, tc4, tw1, tw2, tw3, tw4; // mis à jour le 13/10/22 
    double *vec_dev = malloc(1000 * sizeof(double)); // sauvegarder les 1000 résultats obtenus

    for (u=0; u < 1000; u++) {

        // générér le problème

        if (prob(m, &n, &ia, &ja, &a, &b, rho_ptr, temp_rad, 0)) // temp_rad permet de lancer des simus de problèmes différents
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

        //tc1 = mytimer_cpu(); tw1 = mytimer_wall(); // mis à jour le 13/10/22
        if ( solve_umfpack(n, ia, ja, a, b, x) ) {
            free(ia); free(ja); free(a); free(b); free(x); // empêche leak de mémoire en cas d'erreur
            return 1;
        }
        //tc2 = mytimer_cpu(); tw2 = mytimer_wall(); // mis à jour le 13/10/22
        //

        // sauvegarder le vecteur solution pour faciliter la comparaison, principalement pour debug

        /*FILE *f_x = fopen("mat/x.txt", "w"); // debug
        FILE *f_out = fopen("mat/out.dat", "w");

        for (int i = 0; i < n + 1; i++) {
            fprintf(f_x, "%f\n", (x)[i]); // debug
        }

        fclose(f_x);
        */

        // créer le fichier de sortie pour gnuplot

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

        // calculer gradients sur toute la zone comme autre critère d'uniformité + utile pour flux de chaleur
/*
        int i = 0;
        int j = 0;
        avrg = 0.0; // construit pour donner la moyenne
        dim = 0; // taille de l'échantillon

        double *grad_x = malloc((n + nb_dir) * sizeof(double)); // créer le vecteur grad_x
        double *grad_y = malloc((n + nb_dir) * sizeof(double)); // créer le vecteur grad_y

        for (int iy = 0; iy < m; iy++) { // vertical
            for (int ix = 0; ix < m; ix++) { // horizontal
                // différences centrées, avant ou arrière?
                if ((iy <= 6 * q || ix <= 3 * q) //&& !((iy == 0 && (q * 3 <= ix) && (ix <= q * 8)) || (ix == 0 && (q * 6 <= iy) && (iy <= q * 8)))
                                                ) { // si on n'est pas dans le rectangle supérieur droit*
                    
                    // grad selon x
                    if (ix == 0){ // différence avant
                        if ((q * 6 <= iy) && (iy <= q * 8)) { // porte
                            grad_x[j] = (x[i + 1] - 20) / h;
                        } else {
                            
                        }
                    } else if (ix == 1) {
                        if ((q * 6 <= iy) && (iy <= q * 8)) { // centrée mais attention dirichlet

                        }
                    } else if (ix == m-1) { // différence arrière

                    } else if (ix == 3 * q && iy < q * 6) { // différence arrière

                    } else { // différence centrée

                    }


                    // grad selon y
                    if (iy == 0) { // différence avant
                        if ((q * 3 <= ix) && (ix <= q * 8)) { // fenêtre
                            
                        } else { // pas fenêtre

                        }
                    } else if (iy == 1) { // centrée mais attention dirichlet
                        if ((q * 3 <= ix) && (ix <= q * 8)) { // fenêtre
                            
                        }
                    } else if (iy == m-1) { // différence arrière

                    } else if (ix < q*3 && iy == 6) { // différence arrière
                        
                    } else { // différence centrée

                    }
                }
                


*/



        printf("\n%f    %f\n", temp_rad, std_dev);

        if (u == 0) { // première fois
            save_temp = temp_rad;
            save_dev = std_dev;
        } else if (std_dev < save_dev) { // si on a trouvé un nouveau minimum
            save_temp = temp_rad;
            save_dev = std_dev;
        }

        //free(grad_x); free(grad_y);

        temp_rad += 1.0;

        free(ia); free(ja); free(a); free(b); free(x);

    }

    
    // on a trouvé la meilleure temp du radiateur, on la re-résout pour pouvoir l'afficher ensuite (plus économe que de faire plein de fois de l'écriture de fichier)

    temp_rad = save_temp; // on assigne temp_rad à la meilleure solution

    printf("Meilleure température du radiateur: %f      Meilleur écart-type: %f", save_temp, save_dev);

    if (prob(m, &n, &ia, &ja, &a, &b, rho_ptr, temp_rad, 0)) // temp_rad permet de lancer des simus de problèmes différents
        return 1;
    //printf("\nPROBLEM: ");
    //printf("m = %5d   n = %8d  nnz = %9d\n", m, n, ia[n] );


    // allouer la mémoire pour le vecteur de solution 

    x = malloc(n * sizeof(double));
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

/*
    // calculer flux de chaleur par porte et fenêtre
    // intégrer le gradient de température sur la porte et la fenêtre

    // flux par la fenêtre, iy == 1; q*3 <= ix; ix <= q*8
    // quel indice de x correspond au point directement au-dessus du coin gauche de la fenêtre?
    // m - q*5 - 1 + q*3 (on ne rajoute pas un +1 supplémentaire car les indices commencent à 0) 
    // gradient selon x pas à considérer pour la fenêtre car 1n et 1x orthogonaux
    // donc pas nécessaire de considérer grad_x aux extrémités de la fenêtre

    double flux_x = 0.0, flux_y = 0.0;

    for (int ind = m - q*5 - 1 + q*3; ind <= m + q*3; ind++) {
        // différence avant
        printf("\n%f", x[ind]);
        flux_y += -k * (x[ind] - 0.0); //pow((x[ind] - 0.0) / h, 2.0); // on devrait rajouter un signe moins car normale de la surface mais perte de temps de calcul car de toute façon mis au carré
                                         // on ne divise pas par h car on remultiplie par h (élément dl d'intégration)
    }
  
    //flux_y = sqrt(flux_y); // on prend la racine pour obtenir la norme

    // flux par la porte, ix == 1; q*6 <= iy; iy <= q*8
    // quel indice de x correspond au point directement à droite du bord inférieur de la porte?
    // 6 * q * m - 5 * q - 1
    // indice final?
    // 6 * q * m - 5 * q - 1 + (m-1) + (2 * q - 1) * 3 * q

    for (int ind = 6 * q * m - 5 * q - 1; ind <= 6 * q * m - 5 * q - 1 + (m-1) + (2 * q - 1) * 3 * q; ind += 3 * q) {
        // différence avant
        printf("\n%f", x[ind]);
        flux_x += -k * (x[ind] - 20.0); //pow((x[ind] - 20.0) / h, 2.0); // on devrait rajouter un signe moins mais perte de temps de calcul car de toute façon mis au carré
                                          // on ne divise pas par h car on remultiplie par h (élément dl d'intégration)
        if (ind == m - q*5 + (6 * q - 1) * m - 1 - 1 + 1)
            ind += - 3 * q + m - 1; // la première fois il faut longer tout le mur supérieur droit avant d'atteindre le prochain point en face de la porte
    }

    //flux_x = sqrt(flux_x); // on prend la racine pour obtenir la norme

    printf("\nValeur des 2 flux: (%f, %f) [W]\n", flux_x, flux_y);

    // estimer le flux sortant du radiateur
    // rho de valeur constante sur le radiateur -> il suffit de faire rho * S * k pour obtenir les bonnes valeurs

    save_temp = temp_rad; // au cas où on n'effectue pas la recherche de minimum

    double rad_flux = save_temp * 0.1 * 2.5 * k;
                    
    printf("Flux du radiateur et flux entrant/sortant: %f, %f\n", rad_flux, fabs(flux_x + flux_y));
    printf("Proportion entre les 2 flux: %f\n", rad_flux/fabs(flux_x + flux_y));

    */

    // calculer flux par fenêtre porte et puissance du radiateur

    compute_heat_flux(m, q, &x, save_temp, &flux_x, &flux_y, &rad_flux, rho_ptr); // void qui modifie directement les variables de flux

    printf("\nValeur des 2 flux: (%f, %f) [W]\n", flux_x, flux_y);
    printf("Flux du radiateur et flux entrant/sortant: %f, %f\n", rad_flux, fabs(flux_x + flux_y));
    printf("Proportion entre les 2 flux: %f\n", rad_flux/fabs(flux_x + flux_y));

    // sauvegarder le vecteur solution pour faciliter la comparaison, principalement pour debug

    FILE *f_x = fopen("mat/x.txt", "w"); // debug
    FILE *f_out = fopen("mat/out.dat", "w");

    for (int i = 0; i < n + 1; i++) {
        fprintf(f_x, "%f\n", (x)[i]); // debug
    }

    fclose(f_x);
    
    // créer le fichier de sortie pour gnuplot

    int i = 0;
    avrg = 0.0; // construit pour donner la moyenne
    dim = 0; // taille de l'échantillon

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

    fclose(f_out); // très important, sinon affichage incomplet de out.dat par gnuplot (optimisations compilateur n'attendaient pas l'écriture du fichier?)

    tc3 = mytimer_cpu(); tw3 = mytimer_wall();

    double r = residue(&n, &ia, &ja, &a, &b, &x);
    printf("\nRésidu de la solution: %.10e\n", r);

    tc4 = mytimer_cpu(); tw4 = mytimer_wall();


    printf("\nTemps de solution (CPU): %5.1f sec",tc2-tc1); // mis à jour le 13/10/22 
    printf("\nTemps de solution (horloge): %5.1f sec \n",tw2-tw1); // mis à jour le 13/10/22 
    printf("\nTemps de calcul du résidu (CPU): %5.1f sec",tc4-tc3);
    printf("\nTemps de calcul du résidu (horloge): %5.1f sec \n",tw4-tw3);

    // libérer la mémoire 

    free(ia); free(ja); free(a); free(b); free(x); free(vec_dev);
    system("gnuplot -persist \"heatmap.gnu\""); // laisser gnuplot afficher la température de la pièce
    return 0;
}
