void restriction(int m, double L, int *n, int **ia, int **ja, double **a, double **b, double **x) {
    /*function which restricts the fine vectors to a coarser grid*/

    //restriction matrix
    int *ia_r = malloc((m*m+1)*sizeof(int));
    int *ja_r = malloc((m*m+1)*sizeof(int));
    double *a_r = malloc((m*m+1)*sizeof(double));
    double *b_r = malloc((m*m+1)*sizeof(double));


}