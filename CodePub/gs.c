void fwd_gs(int m, double L, int *n, int **ia, int **ja, double **a, double **b, double **x) {
    // Gauss-seidel smoothing, forward implementation
    // x is the solution vector
    // b is the right hand side
    // ia, ja, a are the CSR matrix
    // n is the number of rows (or columns) of the matrix
    
    double invh2 = (double) (m-1)*(m-1)/(L*L);
    double diag_elem = (double) 4*invh2;
    //printf("diag_elem = %f\n", diag_elem);`
    double sum;
    for (int i = 0; i < *n; i++) {
        sum = -diag_elem*(*x)[i]; // start at -diag elem, because it gets added in the sum afterwards
                                  // less total operations to do this than to have an if check at every iteration
        for (int j = (*ia)[i]; j < (*ia)[i+1]; j++) {
            sum += (*a)[j] * (*x)[(*ja)[j]];
            //printf("sum = %f\n", sum);
            //printf("a[%d] = %f\n", j, (*a)[j]);
        }
        (*x)[i] = ((*b)[i] - sum) / diag_elem; 
        //printf("sum = %f\n", sum);
        //printf("b[%d] = %f\n", i, (*b)[i]);
        //printf("x[%d] = %f\n", i, (*x)[i]);   
    }
}

void bwd_gs(int m, double L, int *n, int **ia, int **ja, double **a, double **b, double **x) {
    /*Gauss-seidel smoothing, backward implementation*/
    // x is the solution vector
    // b is the right hand side
    // ia, ja, a are the CSR matrix
    // n is the number of rows (or columns) of the matrix

    double invh2 = (m-1)*(m-1)/(L*L);
    double diag_elem = 4*invh2;
    double sum;
    for (int i = *n-1; i >= 0; i--) {
        sum = -diag_elem*(*x)[i]; // start with - the diagonal value so that the final sum = all values of the row except the diagonal one
        for (int j = (*ia)[i]; j < (*ia)[i+1]; j++) {
            sum += (*a)[j] * (*x)[(*ja)[j]]; // problem here because we also take the diagonal value
        }
        (*x)[i] = ((*b)[i] - sum) / (diag_elem);
    }
}

void sym_gs(int m, double L, int *n, int **ia, int **ja, double **a, double **b, double **x) {
    /*Symmetric Gauss-seidel smoothing*/
    // x is the solution vector
    // b is the right hand side
    // ia, ja, a are the CSR matrix
    // n is the number of rows (or columns) of the matrix
    
    fwd_gs(m, L, n, ia, ja, a, b, x);
    bwd_gs(m, L, n, ia, ja, a, b, x);
}
