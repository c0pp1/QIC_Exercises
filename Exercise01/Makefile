#pass OPT=-Ox as argument to activate optimizations 
compiler=gfortran

all: exercise01a exercise01b exercise01c
	@echo "Done!"

exercise01a: advanced_ops.f90 exercise01a.f90
	@${compiler} advanced_ops.f90 exercise01a.f90 -o exercise01a

exercise01b: exercise01b.f90
	@${compiler} -fno-range-check exercise01b.f90 -o exercise01b

exercise01c: advanced_ops.f90 exercise01c.f90
	@${compiler} ${OPT} advanced_ops.f90 exercise01c.f90 -o exercise01c

clean:
	@echo "Cleaning up..."
	@rm -f exercise01a exercise01b exercise01c

.PHONY: all clean
