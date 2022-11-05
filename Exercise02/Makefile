compiler=gfortran

all: exercise02a
	@echo "Done!"

exercise02ab: advanced_ops.f90 checkpoint_mod.f90 exercise02ab.f90
	@${compiler} advanced_ops.f90 checkpoint_mod.f90 exercise02ab.f90 -o exercise02ab

clean:
	@echo "Cleaning up..."
	@rm -f exercise01ab exercise01b exercise01c

.PHONY: all clean
