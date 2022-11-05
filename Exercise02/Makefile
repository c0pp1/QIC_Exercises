compiler=gfortran

all: exercise02ab, exercise02c
	@echo "Done!"

exercise02ab: advanced_ops.f90 checkpoint_mod.f90 exercise02ab.f90
	@${compiler} advanced_ops.f90 checkpoint_mod.f90 exercise02ab.f90 -o exercise02ab

exercise02c: exercise02c.f90
	@${compiler} exercise02c.f90 -o exercise02c

clean:
	@echo "Cleaning up..."
	@rm -f exercise01ab exercise01b exercise01c

.PHONY: all clean
