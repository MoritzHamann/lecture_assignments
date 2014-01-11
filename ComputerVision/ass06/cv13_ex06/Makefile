CC = g++

all: pde_morphology

hsTemplate: pde_morphology.c
	${CC} -O3 -o $@ $< -lm
clean:
	rm -f pde_morphology
