double* restriction(int m, int q, int *n, int **ia, int **ja, double **a, double **b, double **x, double **r, double **r_restr);

void prolongation(int m, int q, int *n, double **r_coarse, double **r_prol);