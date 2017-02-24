sample: sample.o bingham.o
	gcc -o sample sample.o bingham.o -lm
sample.o: sample.c bingham.h
	gcc -c sample.c
bingham.o: bingham.c bingham.h
	gcc -c bingham.c
