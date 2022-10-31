#pass OPT=-Ox as argument to activate optimizations 

all: exercise01a exercise01b exercise01c
	@echo "Done!"

exercise01a: advanced_ops.f90 exercise01a.f90
	@gfortran advanced_ops.f90 exercise01a.f90 -o exercise01a

exercise01b: exercise01b.f90
	@gfortran -fno-range-check exercise01b.f90 -o exercise01b

exercise01c: advanced_ops.f90 exercise01c.f90
	@gfortran ${OPT} advanced_ops.f90 exercise01c.f90 -o exercise01c

clean:
	@echo "Cleaning up..."
	@rm -f exercise01a exercise01b exercise01c

.PHONY: all clean