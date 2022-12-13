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

all: exercise06
	@echo "Done!"

exercise06: exercise06.f90 checkpoint_mod.f90 complex_matrix_ops.f90 lapack_wrapper.f90
	@${compiler} -L${libsdir} checkpoint_mod.f90 complex_matrix_ops.f90 lapack_wrapper.f90 exercise06.f90 -o exercise06 ${libs}

clean:
	@echo "Cleaning up..."
	@rm -f exercise06

.PHONY: all clean
