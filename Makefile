#pass OPT=-Ox as argument to activate optimizations 
export OPT

all: exercise01

exercise01:
	@make -C ./Exercise01

clean:
	@make -C ./Exercise01 clean

.PHONY: all clean