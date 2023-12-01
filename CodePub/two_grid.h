double two_grid_method(int n, int m, double L, int *ia, int *ja, double *a, double *b, double *gs_x);

int generate_coarse_problem(int m, int *ia_coarse_ptr, int *ja_coarse_ptr, double *a_coarse_ptr, 
                    double *b_coarse_ptr, void *Numeric); 

int factorized_two_grid_method(int n, int m, double L, int **ia_ptr, int **ja_ptr, double **a_ptr, double *b, double *x, void *Numeric);