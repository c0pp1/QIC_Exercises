compiler=gfortran

all: exercise03a
	@echo "Done!"

exercise03a: advanced_ops.f90 exercise03a.f90
	@${compiler} advanced_ops.f90 exercise03a.f90 -o exercise03a

clean:
	@echo "Cleaning up..."
	@rm -f exercise03a

.PHONY: all clean
