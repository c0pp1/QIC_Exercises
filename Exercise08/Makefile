# modify these variables according to your setup
compiler=gfortran

# if you are running in a conda environment automatically detect it
ifdef CONDA_PREFIX
libsdir=${CONDA_PREFIX}/lib
libs=-llapack -lblas
else
libsdir=/usr/lib
libs=-llapack -lblas -ltmglib
endif
fftlibdir=/usr/local/lib
includedir=/usr/local/include
compileflags=-fbounds-check -O3 -floop-parallelize-all

all: exercise08a exercise08b
	@echo "Done!"

exercise08a: exercise08a.f90 checkpoint_mod.f90 complex_matrix_ops.f90 lapack_wrapper.f90
	@${compiler} ${compileflags} -L${libsdir} checkpoint_mod.f90 complex_matrix_ops.f90 lapack_wrapper.f90 exercise08a.f90 -o exercise08a ${libs}

exercise08b: exercise08b.f90 checkpoint_mod.f90 complex_matrix_ops.f90 lapack_wrapper.f90
	@${compiler} ${compileflags} -L${libsdir} checkpoint_mod.f90 complex_matrix_ops.f90 lapack_wrapper.f90 exercise08b.f90 -o exercise08b ${libs}

clean:
	@echo "Cleaning up..."
	@rm -f exercise08a exercise08b

.PHONY: all clean
