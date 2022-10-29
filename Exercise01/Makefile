all: exercise01a exercise01b
	@echo "Done!"

exercise01a: advanced_ops.f90 exercise01a.f90
	@gfortran advanced_ops.f90 exercise01a.f90 -o exercise01a

exercise01b: exercise01b.f90
	@gfortran -fno-range-check exercise01b.f90 -o exercise01b

clean:
	@echo "Cleaning up..."
	@rm -f exercise01a exercise01b

.PHONY: all clean