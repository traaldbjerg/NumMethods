/* prototype */
int solve_umfpack(int n, int *ia, int *ja, double *a, 
                  double *b, double *x);

void* factorize_umfpack(int n, int *ia, int *ja, double *a);

int solve_umfpack_factorized(int n, int *ia, int *ja, double *a, 
                             double *b, double *x, void *Numeric);

