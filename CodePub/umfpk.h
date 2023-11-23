/* prototype */
int solve_umfpack(int n, int *ia, int *ja, double *a, 
                  double *b, double *x);

int factorize_umfpack(int n, int *ia, int *ja, double *a, void **Symbolic, void **Numeric, double *Info, double *Control);

int solve_umfpack_factorized(int n, int *ia, int *ja, double *a, 
                             double *b, double *x, void **Symbolic, void **Numeric, double *Info, double *Control);