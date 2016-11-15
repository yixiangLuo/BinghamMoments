#include <stdlib.h>
#include <stdio.h>
#include <conio.h>
#include <math.h>
#include "bingham.h"

int main(int argc, char** argv)
{
	/* Input array */
	double b[3];
	/* Result structure */
	Moments M;
	
	/* Initialize the data structure and read the interpolation data. */
	binghamDataSpace* binghamDataPointer= initiateBingham();
	
	/* Calculation */
	int i;
	for(i=0;i<3;i++){
		b[0]= (rand()/(double)(RAND_MAX));
		b[1]= (rand()/(double)(RAND_MAX));
		b[2]= (rand()/(double)(RAND_MAX));
		
		M = bingham(b, binghamDataPointer);
	
		printf("Z:\t\t%.10e\nZ(x1^2):\t%.10e\nZ(x2^2):\t%.10e\n<x1^2>:\t\t%.10e\n<x2^2>:\t\t%.10e\n<x1^2*x2^2>:\t%.10e\n<x1^4>:\t\t%.10e\n<x2^4>:\t\t%.10e\n\n", M.Z, M.Z_x12, M.Z_x22, M.Q_x12, M.Q_x22, M.Q_x12x22, M.Q_x14, M.Q_x24);
	}
	
	/* Destroy the data structure and free the heap memory. */
	binghamDataFree(binghamDataPointer);

	return 0;
}