#include <stdlib.h>
#include <stdio.h>
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
	
	/* Example of case 2.1 in the paper */
	b[0]= -50;
	b[1]= -60;
	b[2]= 0;	
	M = bingham(b, binghamDataPointer);
	printf("Z:\t\t%.10e\nZ(x1^2):\t%.10e\nZ(x2^2):\t%.10e\n<x1^2>:\t\t%.10e\n<x2^2>:\t\t%.10e\n<x1^2*x2^2>:\t%.10e\n<x1^4>:\t\t%.10e\n<x2^4>:\t\t%.10e\n\n", M.Z, M.Z_x12, M.Z_x22, M.Q_x12, M.Q_x22, M.Q_x12x22, M.Q_x14, M.Q_x24);
	
	/* Example of case 2.2 in the paper */
	b[0]= -50;
	b[1]= -10;
	b[2]= 0;	
	M = bingham(b, binghamDataPointer);
	printf("Z:\t\t%.10e\nZ(x1^2):\t%.10e\nZ(x2^2):\t%.10e\n<x1^2>:\t\t%.10e\n<x2^2>:\t\t%.10e\n<x1^2*x2^2>:\t%.10e\n<x1^4>:\t\t%.10e\n<x2^4>:\t\t%.10e\n\n", M.Z, M.Z_x12, M.Z_x22, M.Q_x12, M.Q_x22, M.Q_x12x22, M.Q_x14, M.Q_x24);
	
	/* Example of case 2.3 in the paper */
	b[0]= -10;
	b[1]= -20;
	b[2]= 0;	
	M = bingham(b, binghamDataPointer);
	printf("Z:\t\t%.10e\nZ(x1^2):\t%.10e\nZ(x2^2):\t%.10e\n<x1^2>:\t\t%.10e\n<x2^2>:\t\t%.10e\n<x1^2*x2^2>:\t%.10e\n<x1^4>:\t\t%.10e\n<x2^4>:\t\t%.10e\n\n", M.Z, M.Z_x12, M.Z_x22, M.Q_x12, M.Q_x22, M.Q_x12x22, M.Q_x14, M.Q_x24);
	
	/* Destroy the data structure and free the heap memory. */
	binghamDataFree(binghamDataPointer);

	return 0;
}