# modify these variables according to your setup
compiler=gfortran

# if you are running in a conda environment automatically detect it
ifdef CONDA_PREFIX
libsdir=${CONDA_PREFIX}/lib
libs=-llapack -lblas -lfftw3
else
libsdir=/usr/lib
libs=-llapack -lblas -ltmglib -lfftw3
endif
fftlibdir=/usr/local/lib
includedir=/usr/local/include

all: exercise05
	@echo "Done!"

exercise05: exercise05.f90 checkpoint_mod.f90 lapack_wrapper.f90
	@${compiler} -cpp -L${libsdir} -L${fftlibdir} -I${includedir} checkpoint_mod.f90 lapack_wrapper.f90 exercise05.f90 -o exercise05 ${libs}

clean:
	@echo "Cleaning up..."
	@rm -f exercise05

.PHONY: all clean
