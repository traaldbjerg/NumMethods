double v_cycle(int max_recursion, int n, int m, double L, int *ia, int *ja, double *a, double *b, double *gs_x);

int generate_multigrid_problem(int max_recursion, int m_fine, int **ia_coarse_ptr, int **ja_coarse_ptr, double **a_coarse_ptr, 
                    double **b_coarse_ptr, void **Numeric_ptr);
