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

all: exercise03a exercise03b
	@echo "Done!"

exercise03a: advanced_ops.f90 exercise03a.f90
	@${compiler} advanced_ops.f90 exercise03a.f90 -o exercise03a

exercise03b: exercise03b.f90 lapack_wrapper.f90
	@${compiler} -L${libsdir} lapack_wrapper.f90 exercise03b.f90 -o exercise03b ${libs}

clean:
	@echo "Cleaning up..."
	@rm -f exercise03a exercise03b

.PHONY: all clean
