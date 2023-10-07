void fwd_gs(int m, double L, int *n, int **ia, int **ja, double **a, double **b, double **x) {
    /*Gauss-seidel smoothing, forward implementation*/
    double invh2 = (m-1)*(m-1)/(L*L);
    double sum;
    for (int i = 0; i < *n; i++) {
        sum = 0;
        for (int j = (*ia)[i]; j < (*ia)[i+1]; j++) {
            sum += (*a)[j] * (*x)[(*ja)[j]];
        }
        (*x)[i] = ((*b)[i] - sum) / (4 * invh2); // not very clean to not use the matrix value but would be difficult to find
    }
}

void bwd_gs(int m, double L, int *n, int **ia, int **ja, double **a, double **b, double **x) {
    /*Gauss-seidel smoothing, backward implementation*/
    double invh2 = (m-1)*(m-1)/(L*L);
    double sum;
    for (int i = *n-1; i >= 0; i--) {
        sum = 0;
        for (int j = (*ia)[i]; j < (*ia)[i+1]; j++) {
            sum += (*a)[j] * (*x)[(*ja)[j]];
        }
        (*x)[i] = ((*b)[i] - sum) / (4 *invh2);
    }
}

void sym_gs(int m, double L, int *n, int **ia, int **ja, double **a, double **b, double **x) {
    /*Symmetric Gauss-seidel smoothing*/
    fwd_gs(m, L, n, ia, ja, a, b, x);
    bwd_gs(m, L, n, ia, ja, a, b, x);
}
