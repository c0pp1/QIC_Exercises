# modify these variables according to your setup
compiler=gfortran
libs=$CONDA_PREFIX/lib

all: exercise03a exercise03b
	@echo "Done!"

exercise03a: advanced_ops.f90 exercise03a.f90
	@${compiler} advanced_ops.f90 exercise03a.f90 -o exercise03a

exercise03b: exercise03b.f90 lapack_wrapper.f90
	@${compiler} -L ${libs} lapack_wrapper.f90 exercise03b.f90 -o exercise03b -llapack -lblas

clean:
	@echo "Cleaning up..."
	@rm -f exercise03a exercise03b

.PHONY: all clean
