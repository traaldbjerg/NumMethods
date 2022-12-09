#include <stdio.h>
#include "time.h"
#include "petscksp.h"


//extern PetscErrorCode ComputeMatrix(KSP, Mat, Mat, void *);
//extern PetscErrorCode ComputeRHS(KSP, Vec, void *);
//extern PetscErrorCode ComputeInitialGuess(KSP, Vec, void *);

void solve_petsc(int n, int *ia, int *ja, double *a, double *b, double **x) {

    // utiliser les solveurs de PETSc pour résoudre le système linéaire Ax = b
    // ceci requiert d'abord de convertir nos arrays en un format utilisable par PETSc

    PetscErrorCode ierr;

    ierr = PetscInitializeNoArguments(); // instancier toutes les variables nécessaires à PETSc

    KSP ksp;
    //DM da;
    PC pc;
    

    int nnz = ia[n];
    int *ixn = malloc(n * sizeof(int));
    int *ixnnz = malloc(nnz * sizeof(int));

    for (int i = 0; i < n; i++) {ixn[i] = i; }

    //printf("Hello\n"); // debug

    // déclarer les variables
    Mat A;
    Vec vec_b;
    Vec vec_x;



    MatCreateSeqAIJWithArrays(PETSC_COMM_WORLD, n, n, ia, ja, a, &A);

    // créer le vecteur b
    VecCreate(PETSC_COMM_WORLD, &vec_b);
    VecSetSizes(vec_b, PETSC_DECIDE, n);
    VecSetFromOptions(vec_b);
    VecSetValues(vec_b, n, ixn, b, INSERT_VALUES);
    VecAssemblyBegin(vec_b);
    VecAssemblyEnd(vec_b);

    // initialiser le vecteur x (sinon KSP considère vec_x comme un pointeur nul)
    VecCreate(PETSC_COMM_WORLD, &vec_x);
    VecSetSizes(vec_x, PETSC_DECIDE, n);
    VecSetFromOptions(vec_x);

    KSPCreate(PETSC_COMM_WORLD, &ksp); // créer le contexte pour le solveur
    KSPSetOperators(ksp, A, A); // A et A => pas de préconditionnement pour le moment
    KSPSetFromOptions(ksp);

    KSPGetPC(ksp, &pc); // lancer contexte pour le préconditionneur

    double tc5 = mytimer_cpu(); double tw5 = mytimer_wall();

    KSPSolve(ksp, vec_b, vec_x);

    double tc6 = mytimer_cpu(); double tw6 = mytimer_wall();

    printf("\nTemps de solution, PETSc (CPU): %5.1f sec",tc6-tc5);  
    printf("\nTemps de solution, PETSc (horloge): %5.1f sec \n",tw6-tw5); 

    KSPDestroy(&ksp);

    VecGetArray(vec_x, x);

   ierr = PetscFinalize();
}
