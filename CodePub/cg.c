#include <stdio.h>
#include "multigrid.h"
#include "cg.h"
#include "matrix.h"
#include "residual.h"

int flexible_cg(int max_recursion, int n, int m, double L, int **ia_ptr, int **ja_ptr,
                             double **a_ptr, double *b, double *x, double *r, double *d, void *Numeric, int v_w) {

    // conjugate gradient method, flexible version
    // construct B^-1 r_m, where B is the preconditionner
    // preconditionner is a multigrid method
    // v_w = 0 for v-cycle, v_w = 1 for w-cycle

    (void) residual(&n, &ia_ptr[0], &ja_ptr[0], &a_ptr[0], &b, &x, &r); // use the finest grid to calculate the residual
    double *B_r = calloc(1, n * sizeof(double));

    //for (int i = 0; i < n; i++) {
    //    B_r[i] = r[i];
    //}

    if (v_w == 0) {
        v_cycle(max_recursion, 0, n, m, L, ia_ptr, ja_ptr, a_ptr, r, B_r, Numeric); // makes B_r = B^-1 r_m
    } else if (v_w == 1) {
        w_cycle(max_recursion, 0, n, m, L, ia_ptr, ja_ptr, a_ptr, r, B_r, Numeric); // makes B_r = B^-1 r_m
    } else {
        printf("Error: v_w must be 0 or 1\n");
        return 1;
    }
    double *A_B_r = malloc(n * sizeof(double));

    //for (int i = 0; i < n; i++) {
    //   printf("B_r[%d] = %f\n", i, B_r[i]);
    //}

    CSR_matrix_multiplication(n, ia_ptr[0], ja_ptr[0], a_ptr[0], B_r, A_B_r); // makes a_B_r = A * B_r

    //for (int i = 0; i < n; i++) {
    //    printf("a_B_r[%d] = %f\n", i, a_B_r[i]);
    //}

    double *A_d = malloc(n * sizeof(double));
    CSR_matrix_multiplication(n, ia_ptr[0], ja_ptr[0], a_ptr[0], d, A_d); // makes a_d = A * d

    //for (int i = 0; i < n; i++) {
    //    printf("a_d[%d] = %f\n", i, a_d[i]);
    //}

    double beta;

    if (dot_product(n, d, A_d) == 0) { // first iteration because d = 0
        beta = 0;
    } else {
        //printf("dot_product(n, d, a_B_r) = %f\n", dot_product(n, d, A_B_r));
        //printf("dot_product(n, d, a_d) = %f\n", dot_product(n, d, A_d));
        beta = - dot_product(n, d, A_B_r) / dot_product(n, d, A_d); // makes beta = d^T * a_B_r / d^T * a_d
    }
    //printf("beta = %f\n", beta);
    double *beta_d = malloc(n * sizeof(double));
    scalar_product(n, beta, d, beta_d); // makes beta_d = beta * d

    //for (int i = 0; i < n; i++) {
        //printf("beta_d[%d] = %f\n", i, beta_d[i]);
    //}

    vector_addition(n, B_r, beta_d, d); // makes d = B^-1 r + beta * d

    //for (int i = 0; i < n; i++) {
    //    printf("d[%d] = %f\n", i, d[i]);
    //}

    CSR_matrix_multiplication(n, ia_ptr[0], ja_ptr[0], a_ptr[0], d, A_d);
    
    double alpha = dot_product(n, d, r) / dot_product(n, d, A_d); // makes alpha = d^T * a_B_r / d^T * d

    //printf("alpha = %f\n", alpha);

    double *alpha_d = malloc(n * sizeof(double));
    scalar_product(n, alpha, d, alpha_d); // makes alpha_d = alpha * d
    vector_addition(n, x, alpha_d, x); // makes x = x + alpha * d, it is ok to modify directly x
                                       //because we don't need a component again after it has been modified
    scalar_product(n, -alpha, A_d, A_d); // makes a_d = - alpha * a_d
    vector_addition(n, r, A_d, r); // makes r = r - alpha * a_d

    free(B_r); free(A_B_r); free(A_d); free(beta_d); free(alpha_d);

    return 0;

}


int standard_cg(int max_recursion, int n, int m, double L, int **ia_ptr, int **ja_ptr,
                             double **a_ptr, double *b, double *x, double *r, double *r_B_r_save, double *d, void *Numeric) {

    // conjugate gradient method, flexible version
    // construct B^-1 r_m, where B is the preconditionner

    // preconditionner is a multigrid method

    //for (int i = 0; i < n; i++) {
        //printf("r_save[%d] = %f\n", i, r_save[i]);
        //printf("B_r_save[%d] = %f\n", i, B_r_save[i]);
    //}

    (void) residual(&n, &ia_ptr[0], &ja_ptr[0], &a_ptr[0], &b, &x, &r); // update the residual
    //printf("res = %f\n", res);
    double *B_r = calloc(1, n * sizeof(double));
    v_cycle(max_recursion, 0, n, m, L, ia_ptr, ja_ptr, a_ptr, r, B_r, Numeric); // makes B_r = B^-1 r_m

    //printf("r_B_r_save = %f\n", *r_B_r_save);

    //for (int i = 0; i < n; i++) {
        //printf("B_r[%d] = %f\n", i, B_r[i]);
    //}

    double *A_d = malloc(n * sizeof(double));

    double r_B_r = dot_product(n, r, B_r); // makes r_B_r = r * B^-1 r

    double beta;
    if (*r_B_r_save == 0) { // first iteration, the save vectors are calloced so this condition will be true
        beta = 0;
    } else {
        //printf("dot_product(n, r, B_r) = %f\n", dot_product(n, r, B_r));
        //printf("dot_product(n, r_save, B_r_save) = %f\n", dot_product(n, r_save, B_r_save));
        beta = r_B_r / *r_B_r_save; // makes beta = r * B^-1 r / r_m-1 * B^-1 r_m-1
    }

    //printf("r_B_r = %f\n", r_B_r);

    //printf("beta = %f\n", beta);

    double *beta_d = malloc(n * sizeof(double));
    scalar_product(n, beta, d, beta_d); // makes beta_d = beta * d
    vector_addition(n, B_r, beta_d, d); // makes d = B^-1 r + beta * d

    CSR_matrix_multiplication(n, ia_ptr[0], ja_ptr[0], a_ptr[0], d, A_d);
    double alpha = dot_product(n, r, B_r) / dot_product(n, d, A_d); // makes alpha = d^T * a_B_r / d^T * d

    //printf("alpha = %f\n", alpha);

    double *alpha_d = malloc(n * sizeof(double));
    scalar_product(n, alpha, d, alpha_d); // makes alpha_d = alpha * d
    vector_addition(n, x, alpha_d, x); // makes x = x + alpha * d, it is ok to modify directly x
                                       // because we don't need a component again after it has been modified

    scalar_product(n, -alpha, A_d, A_d); // makes a_d = - alpha * a_d
    vector_addition(n, r, A_d, r); // makes r = r - alpha * a_d 

    *r_B_r_save = r_B_r;

    free(B_r); free(A_d); free(beta_d); free(alpha_d);

    return 0;

}