#include <math.h>
#include <stdlib.h>
#include <stdio.h>

double residue(int *n, int **ia, int **ja, double **a, double **b, double **x) {
    
    double res = 0;
    //int nnz = sizeof(ja) / sizeof(int);
    double norm_b = 0;

    // allocation de la mémoire pour le vecteur A * x
    double *a_x = malloc(*n * sizeof(double)); // A * x pour pouvoir facilement faire A * x - b par la suite

    /*for (int i = 0; i < n; i++) {
        norm_b += pow((*b)[i], 2.0);
    }
    norm_b = sqrt(norm_b);
    */

    // construction de a_x
    for (int i = 0; i < *n; i++) { // chaque compo de a_x
        a_x[i] = 0;
        for (int k = (*ia)[i]; k < (*ia)[i + 1]; k++) { // k permet d'accéder aux bons indices de ja
            a_x[i] += (*a)[k] * (*x)[(*ja)[k]]; // le fait d'accéder successivement à a et ja permet de les garder en cache et augmenter la vitesse d'exécution
        }
        norm_b += pow((*b)[i], 2.0); // peut-être pas optimal d'aussi accéder à b ici
        res += pow(a_x[i] - (*b)[i], 2.0); // le numérateur de res
    }

    norm_b = sqrt(norm_b);
    res = sqrt(res);
    printf("\nNorme de b: %.10e\n", norm_b); // debug, actuellement du genre 1e07
    res /= norm_b;

    return res;
}
