#include <stdlib.h>
#include <stdio.h>
#include <conio.h>
#include <math.h>
#include "bingham.h"

int main(int argc, char** argv)
{
	/* Input array, containing the three eigenvalues of matrix B */
	double b[3];
	/* Result array, whose elements are in the order of Z, Z(x1^2), Z(x2^2), <x1^2>, <x2^2>, <x1^2*x2^2>, <x1^4>, <x2^4> */
	double Q[8]={0};
	
	/* Initialize the data structure and read the interpolation data. */
	binghamDataSpace* binghamDataPointer= initiateBingham();
	
	/* Calculation */
	int i;
	for(i=0;i<3;i++){
		b[0]= (rand()/(double)(RAND_MAX));
		b[1]= (rand()/(double)(RAND_MAX));
		b[2]= (rand()/(double)(RAND_MAX));
		
		bingham(b, binghamDataPointer, Q);
	
		printf("Z:\t\t%.10e\nZ(x1^2):\t%.10e\nZ(x2^2):\t%.10e\n<x1^2>:\t\t%.10e\n<x2^2>:\t\t%.10e\n<x1^2*x2^2>:\t%.10e\n<x1^4>:\t\t%.10e\n<x2^4>:\t\t%.10e\n\n",Q[0], Q[1], Q[2], Q[3], Q[4], Q[5], Q[6], Q[7]);
	}
	
	/* Destroy the data structure and free the heap memory. */
	binghamDataFree(binghamDataPointer);

	return 0;
}
