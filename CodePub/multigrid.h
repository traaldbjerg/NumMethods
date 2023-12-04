int v_cycle(int max_recursion, int c, int n, int m, double L,
                                     int **ia_ptr, int **ja_ptr, double **a_ptr, double *b, double *x, void *Numeric);

int w_cycle(int max_recursion, int c, int n, int m, double L,
                                     int **ia_ptr, int **ja_ptr, double **a_ptr, double *b, double *x, void *Numeric);

void *generate_multigrid_problem(int max_recursion, int m, int **ia_ptr, int **ja_ptr, double **a_ptr, 
                    double **b_ptr);
