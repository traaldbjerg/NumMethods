int generate_two_grid_problem(int m, int **ia_ptr, int **ja_ptr, double **a_ptr, 
                    double **b_ptr, void *Numeric); 

int two_grid_method(int n, int m, double L, int **ia_ptr, int **ja_ptr, double **a_ptr, double *b, double *x, void *Numeric);