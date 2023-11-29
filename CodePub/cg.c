#include <stdio.h>
#include "multigrid.h"
#include "preconditionner.h"
#include "cg.h"
#include "matrix.h"

int flexible_cg(int max_recursion, int n, int m, double L, int **ia_ptr, int **ja_ptr,
                             double **a_ptr, double *b, double *x, double *r, double *d, void *Numeric) {

    // conjugate gradient method, flexible version
    // construct B^-1 r_m, where B is the preconditionner
    // preconditionner is a multigrid method
    (void) residual(&n, &ia_ptr[0], &ja_ptr[0], &a_ptr[0], &b, &x, &r); // use the finest grid to calculate the residual
    double *b_r = calloc(1, n * sizeof(double));

    //for (int i = 0; i < n; i++) {
    //    b_r[i] = r[i];
    //}

    w_cycle(max_recursion, 0, n, m, L, ia_ptr, ja_ptr, a_ptr, r, b_r, Numeric); // makes b_r = B^-1 r_m
    double *a_b_r = malloc(n * sizeof(double));

    //for (int i = 0; i < n; i++) {
    //    printf("b_r[%d] = %f\n", i, b_r[i]);
    //}

    CSR_matrix_multiplication(n, ia_ptr[0], ja_ptr[0], a_ptr[0], b_r, a_b_r); // makes a_b_r = A * b_r

    //for (int i = 0; i < n; i++) {
    //    printf("a_b_r[%d] = %f\n", i, a_b_r[i]);
    //}

    double *a_d = malloc(n * sizeof(double));
    CSR_matrix_multiplication(n, ia_ptr[0], ja_ptr[0], a_ptr[0], d, a_d); // makes a_d = A * d

    //for (int i = 0; i < n; i++) {
    //    printf("a_d[%d] = %f\n", i, a_d[i]);
    //}

    double beta;

    if (dot_product(n, d, a_d) == 0) { // first iteration because d = 0
        beta = 0;
    } else {
        beta = - dot_product(n, d, a_b_r) / dot_product(n, d, a_d); // makes beta = d^T * a_b_r / d^T * a_d
    }
    //printf("beta = %f\n", beta);
    double *beta_d = malloc(n * sizeof(double));
    scalar_product(n, beta, d, beta_d); // makes beta_d = beta * d

    //for (int i = 0; i < n; i++) {
        //printf("beta_d[%d] = %f\n", i, beta_d[i]);
    //}

    vector_addition(n, b_r, beta_d, d); // makes d = B^-1 r + beta * d

    //for (int i = 0; i < n; i++) {
    //    printf("d[%d] = %f\n", i, d[i]);
    //}

    CSR_matrix_multiplication(n, ia_ptr[0], ja_ptr[0], a_ptr[0], d, a_d);
    
    double alpha = dot_product(n, d, r) / dot_product(n, d, a_d); // makes alpha = d^T * a_b_r / d^T * d
    double *alpha_d = malloc(n * sizeof(double));
    scalar_product(n, alpha, d, alpha_d); // makes alpha_d = alpha * d
    vector_addition(n, x, alpha_d, x); // makes x = x + alpha * d, it is ok to modify directly x
                                       //because we don't need a component again after it has been modified
    scalar_product(n, -alpha, a_d, a_d); // makes a_d = - alpha * a_d
    vector_addition(n, r, a_d, r); // makes r = r - alpha * a_d

    free(b_r); free(a_b_r); free(a_d); free(beta_d); free(alpha_d);

    return 0;

}


int standard_cg(int max_recursion, int n, int m, double L, int **ia_ptr, int **ja_ptr,
                             double **a_ptr, double *b, double *x, double *r, double *r_save, double *b_r_save, double *d, void *Numeric) {

    // conjugate gradient method, flexible version
    // construct B^-1 r_m, where B is the preconditionner

    // preconditionner is a multigrid method
    (void) residual(&n, &ia_ptr[0], &ja_ptr[0], &a_ptr[0], &b, &x, &r); // use the finest grid to calculate the residual
    double *b_r = calloc(1, n * sizeof(double));
    w_cycle(max_recursion, 0, n, m, L, ia_ptr, ja_ptr, a_ptr, r, b_r, Numeric); // makes b_r = B^-1 r_m

    double *a_d = malloc(n * sizeof(double));

    double beta;
    if (dot_product(n, r_save, b_r_save) == 0) { // first iteration, the save vectors are calloced so this condition will be true
        beta = 0;
    } else {
        beta = - dot_product(n, r, b_r) / dot_product(n, r_save, b_r_save); // makes beta = d^T * a_b_r / d^T * a_d
    }

    for (int i = 0; i < n; i++) { // update the save vectors for the next iteration
        r_save[i] = r[i];
        b_r_save[i] = b_r[i];
    }

    double *beta_d = malloc(n * sizeof(double));
    scalar_product(n, beta, d, beta_d); // makes beta_d = beta * d
    vector_addition(n, b_r, beta_d, d); // makes d = B^-1 r + beta * d

    CSR_matrix_multiplication(n, ia_ptr[0], ja_ptr[0], a_ptr[0], d, a_d);
    double alpha = dot_product(n, r, b_r) / dot_product(n, d, a_d); // makes alpha = d^T * a_b_r / d^T * d
    double *alpha_d = malloc(n * sizeof(double));
    scalar_product(n, alpha, d, alpha_d); // makes alpha_d = alpha * d
    vector_addition(n, x, alpha_d, x); // makes x = x + alpha * d, it is ok to modify directly x
                                       // because we don't need a component again after it has been modified

    scalar_product(n, -alpha, a_d, a_d); // makes a_d = - alpha * a_d
    vector_addition(n, r, a_d, r); // makes r = r - alpha * a_d

    free(b_r); free(a_d); free(beta_d); free(alpha_d);

    return 0;

}