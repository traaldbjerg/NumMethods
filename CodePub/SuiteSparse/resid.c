void residual(int n, int *col_ptr, int *row_ind, double *val,
           double *x, double *b)
{
    int i, j;
    for(i = 0; i < n; i++) {
        for (j = col_ptr[i]; j < col_ptr[i + 1]; j++) {
//            printf("b = %f --> ",b[row_ind[j]]);
            b[row_ind[j]] -= val[j] * x[i];
//            printf("b = %f, i = %d, j= %d aij = %f\n",b[row_ind[j]], i, row_ind[j],val[j]);
        }
    }
}
