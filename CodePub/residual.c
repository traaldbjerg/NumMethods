#include <math.h>
#include <stdlib.h>
#include <stdio.h>


double residual(int *n, int **ia, int **ja, double **a, double **b, double **x, double **r) {
    // computes the relative norm of the residual and modifies the residual vector
    // ia, ja and a are the CSR representation of the matrix A
    // b is the rhs vector
    // x is the solution vector
    // r is the residual vector
    
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
        (*r)[i] = (*b)[i] - current_component_a_x; // we store the current component in r
        //printf("This is r[%d] = %f\n", i, (*r)[i]); // debug
        //printf("This is b[%d] = %f\n", i, (*b)[i]); // debug
        //printf("This is current_component_a_x = %f\n", current_component_a_x); // debug
        res += (*r)[i] * (*r)[i]; // le numérateur de res
    }

    norm_b = sqrt(norm_b);
    res = sqrt(res);
    //printf("\nNorme de b: %.10e\n", norm_b); // debug
    res /= norm_b;

    return res;
}

double residual_norm(int *n, int **ia, int **ja, double **a, double **b, double **x) {
    // computes the relative norm of the residual without modifying the residual vector
    
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
        res += ((*b)[i] - current_component_a_x) * ((*b)[i] - current_component_a_x); // we store the current component in r
        //printf("This is r[%d] = %f\n", i, (*r)[i]); // debug
        //printf("This is b[%d] = %f\n", i, (*b)[i]); // debug
        //printf("This is current_component_a_x = %f\n", current_component_a_x); // debug
    }

    norm_b = sqrt(norm_b);
    res = sqrt(res);
    //printf("res norm = %.10e\n", res); // debug
    //printf("b norm = %.10e\n", norm_b); // debug
    res = res / norm_b;
    //printf("res norm / b norm = %.10e\n", res); // debug

    return res;
}
