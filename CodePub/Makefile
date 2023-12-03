# librairies de SuiteSparse
# -Lpetsc/lib -lpetscsys

PETSC_DIR = petsc
PETSC_ARCH = arch-linux-c-debug

L1 = SuiteSparse/UMFPACK/Lib/libumfpack.a
L2 = SuiteSparse/CHOLMOD/Lib/libcholmod.a 
L3 = SuiteSparse/AMD/Lib/libamd.a 
L4 = SuiteSparse/CAMD/Lib/libcamd.a  
L5 = SuiteSparse/COLAMD/Lib/libcolamd.a 
L6 = SuiteSparse/CCOLAMD/Lib/libccolamd.a 
L7 = SuiteSparse/metis-4.0/libmetis.a
L8 = SuiteSparse/SuiteSparse_config/libsuitesparseconfig.a
LIB1 = $(L1) $(L2) $(L3) $(L4) $(L5) $(L6) $(L7) $(L8) -lm -lblas -llapack

#L9 = petsc/arch-linux-c-debug/lib/libpetsc.so
#LIB2 = $(L9)

COPT = -O3 -Wall

#-Xlinker --verbose

#include ${PETSC_DIR}/lib/petsc/conf/variables

#  To access the PETSc variables for the build, including compilers, compiler flags, libraries etc but
#  manage the build rules yourself (rarely needed) comment out the next lines
#
#include ${PETSC_DIR}/lib/petsc/conf/rules
#include ${PETSC_DIR}/lib/petsc/conf/test

default: main

clean: 
	rm *.o 
	rm main

main: main.c prob.o time.o  umfpack.o residual.o gs.o projections.o two_grid.o multigrid.o cg.o matrix.o
	cc -g $(COPT) $^ -o $@ $(LIB1) $(LIB2) -ISuiteSpare/UMFPACK/Include \
  -ISuiteSparse/SuiteSparse_config  -ISuiteSparse/AMD/Include

umfpack.o: umfpack.c
	cc -g $(COPT) -c $< -o $@ -ISuiteSparse/UMFPACK/Include \
  -ISuiteSparse/SuiteSparse_config  -ISuiteSparse/AMD/Include

#petsc.o: petsc.c
#	cc -g $(COPT) -c $< -o $@ $(LIB2) -Ipetsc/include -Ipetsc/arch-linux-c-debug/include

#-Ipetsc/include -Ipetsc/arch-linux-c-debug/include -Lpetsc/arch-linux-c-debug/include/petscconf.h -Lpetsc/lib -lpetsc -lpetscsys -lpetscksp

%.o: %.c
	cc -g $(COPT) -c $< -o $@

#include ${PETSC_DIR}/lib/petsc/conf/variables
#include ${PETSC_DIR}/lib/petsc/conf/rules
#include ${PETSC_DIR}/lib/petsc/conf/tests
