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

all: exercise04
	@echo "Done!"

exercise04: exercise04.f90 lapack_wrapper.f90
	@${compiler} -L${libsdir} lapack_wrapper.f90 exercise04.f90 -o exercise04 ${libs}

clean:
	@echo "Cleaning up..."
	@rm -f exercise04

.PHONY: all clean
