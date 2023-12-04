#include "matrix.h"

void CSR_matrix_multiplication(int n, int *ia, int *ja, double *a, double *x, double *out) {
    // computes out = A * x
    // A is a sparse matrix, stored in CSR format in ia, ja and a
    // x and b are vectors of size n
    // ia and ja are vectors of size n+1 and n respectively
    // a is a vector of size ia[n]
    // out is a vector of size n
    for (int i = 0; i < n; i++) { // chaque compo de a_x
        out[i] = 0.0; // nouvelle composante donc on réinitialise
        for (int k = (ia)[i]; k < (ia)[i + 1]; k++) { // k permet d'accéder aux bons indices de ja
            (out)[i] += (a)[k] * (x)[(ja)[k]]; 
        }
    }
}

double dot_product(int n, double *x, double *y) {
    // computes out = x^T * y
    // x and y are vectors of size n
    // out is a scalar
    double out = 0;
    for (int i = 0; i < n; i++) {
        out += x[i] * y[i];
    }
    return out;
}

void scalar_product(int n, double s, double *x, double *out) {
    // computes out = s * x
    // x and out are vectors of size n
    // s is a scalar
    for (int i = 0; i < n; i++) {
        out[i] = s * x[i];
    }
}

void vector_addition(int n, double *x, double *y, double *out) {
    // computes out = x + y
    // x, y and out are vectors of size n
    for (int i = 0; i < n; i++) {
        out[i] = x[i] + y[i];
    }
}

double norm(int n, double *x) {
    // computes the norm of x
    // x is a vector of size n
    double out = 0;
    for (int i = 0; i < n; i++) {
        out += x[i] * x[i];
    }
    return sqrt(out);
}