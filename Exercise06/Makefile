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

all: exercise06
	@echo "Done!"

exercise06: exercise06.f90 checkpoint_mod.f90 
	@${compiler} checkpoint_mod.f90 exercise06.f90 -o exercise06

clean:
	@echo "Cleaning up..."
	@rm -f exercise06

.PHONY: all clean
