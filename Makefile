#pass OPT=-Ox as argument to activate optimizations 
export OPT

all: exercise01 exercise02

exercise01:
	@make -C ./Exercise01

exercise02:
	@make -C ./Exercise02

clean:
	@make -C ./Exercise01 clean
	@make -C ./Exercise02 clean

.PHONY: all clean