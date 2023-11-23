double two_grid_method(int n, int m, int L, int *ia, int *ja, double *a, double *b, double *x, double *gs_x);

int generate_coarse_problem(int m, int *ia_coarse_ptr, int *ja_coarse_ptr, double *a_coarse_ptr, 
                    double *b_coarse_ptr, void **Symbolic, void **Numeric, double *Info, double *Control); 
