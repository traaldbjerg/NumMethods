#include <math.h>
#include <stdlib.h>
#include <stdio.h>

/*double residue(int *n, int **ia, int **ja, double **a, double **b, double **x) {
    
    double res = 0;
    double norm_b = 0;

    // allocation de la mémoire pour le vecteur A * x
    //double *a_x = malloc(*n * sizeof(double)); // A * x pour pouvoir facilement faire A * x - b par la suite
    double current_component_a_x; // plus efficace d'allouer un seul double qu'on réutilise plutôt qu'un vecteur complet

    // construction de a_x
    for (int i = 0; i < *n; i++) { // chaque compo de a_x
        current_component_a_x = 0; // nouvelle composante donc on réinitialise
        for (int k = (*ia)[i]; k < (*ia)[i + 1]; k++) { // k permet d'accéder aux bons indices de ja
            current_component_a_x += (*a)[k] * (*x)[(*ja)[k]]; // le fait d'accéder successivement à a et ja permet de les garder en cache et augmenter la vitesse d'exécution
        }
        norm_b += pow((*b)[i], 2.0); // peut-être pas optimal d'aussi accéder à b ici
        res += pow(current_component_a_x - (*b)[i], 2.0); // le numérateur de res
    }

    norm_b = sqrt(norm_b);
    res = sqrt(res);
    //printf("\nNorme de b: %.10e\n", norm_b); // debug
    res /= norm_b;

    return res;
}*/

double residue(int *n, int **ia, int **ja, double **a, double **b, double **x, double **r) {
    
    double res = 0;
    double norm_b = 0;

    // allocation de la mémoire pour le vecteur A * x
    //double *a_x = malloc(*n * sizeof(double)); // A * x pour pouvoir facilement faire A * x - b par la suite
    double current_component_a_x; // plus efficace d'allouer un seul double qu'on réutilise plutôt qu'un vecteur complet

    // construction de a_x
    for (int i = 0; i < *n; i++) { // chaque compo de a_x
        current_component_a_x = 0; // nouvelle composante donc on réinitialise
        for (int k = (*ia)[i]; k < (*ia)[i + 1]; k++) { // k permet d'accéder aux bons indices de ja
            current_component_a_x += (*a)[k] * (*x)[(*ja)[k]]; // le fait d'accéder successivement à a et ja permet de les garder en cache et augmenter la vitesse d'exécution
        }
        norm_b += (*b)[i] * (*b)[i]; // peut-être pas optimal d'aussi accéder à b ici
        (*r)[i] = current_component_a_x - (*b)[i]; // we store the current component in r
        res += (*r)[i] * (*r)[i]; // le numérateur de res
    }

    norm_b = sqrt(norm_b);
    res = sqrt(res);
    //printf("\nNorme de b: %.10e\n", norm_b); // debug
    res /= norm_b;

    return res;
}
