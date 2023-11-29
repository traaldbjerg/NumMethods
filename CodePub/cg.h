int flexible_cg(int max_recursion, int n, int m, double L, int **ia_ptr, int **ja_ptr,
                             double **a_ptr, double *b, double *x, double *r, double *d, void *Numeric);

int standard_cg(int max_recursion, int n, int m, double L, int **ia_ptr, int **ja_ptr,
                             double **a_ptr, double *b, double *x, double *r, double *r_save, double *b_r_save, double *d, void *Numeric);