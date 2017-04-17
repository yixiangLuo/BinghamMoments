# [Introduction](#introduction)

This library is a C port of the implementation of the method for computing the moment functions of Bingham distribution described in the paper [_A Fast Algorithm for the Moments of Bingham Distribution_](https://arxiv.org/abs/1612.01015).

The density function of the Bingham distribution is given by

![the density function of the Bingham distribution](http://i1149.photobucket.com/albums/o595/yiwanjiangyou/densOrg_zps7g3u8lge.png?t=1479611442)

This library provides a fast way to calculate its moment function

![moment function](http://i1149.photobucket.com/albums/o595/yiwanjiangyou/moments_zpspkn9etgi.png?t=1479611376)

up to 4th order.
The calculation is done by piecewise rational approximation, where interpolation and Gaussian integrals are utilized. It is shown by numerical test that an accuracy of 5e−8 is attained remarkably faster than direct numerical quadrature.

The library provides several features:
* **Cross platform**. The source code can be compiled on Microsoft Visual Studio, GNU C Compiler (gcc), etc.

# [Reference](#reference)

(reference of the paper)

# [How to use](#how-to-use)

## [Compiling](#compiling)

Save the files in a fresh subdirectory on your system. Use

`make`

or type

`gcc -c bingham.c`

`gcc -c sample.c`

`gcc -o sample sample.o bingham.o`

## [Sample code](#sample-code)

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

# [Documentation](#documentation)

## [Data Structures](#data-structures)

### [Constants](#constants)

> Divide

the dividing number, which is denoted as d in the referenced paper

> InterpNodesNumAxis_2D_Z

total number of 2D interpolation nodes along x(and y) axis, in calculating Z, Z(x1^2) and Z(x2^2)

> InterpNodesNumAxis_2D_Q

total number of 2D interpolation nodes along x(and y) axis, in calculating \<x1^2\>, \<x2^2\>, \<x1^2*x2^2\>, \<x1^4\> and \<x2^4\>

> AreaSpacing_Z

grid size for 2D interpolation, in calculating Z, Z(x1^2) and Z(x2^2)

> AreaSpacing_Q

grid size for 2D interpolation, in calculating \<x1^2\>, \<x2^2\>, \<x1^2*x2^2\>, \<x1^4\> and \<x2^4\>

> InterpNodesNumAxis_1D

total number of 1D interpolation nodes, in calculating the derivatives of g in case 2.2 in the referenced paper.

> LinearSpacing

grid size for 1D interpolation, in calculating the derivatives of g in case 2.2 in the referenced paper.

### [Structures](#structures)

> binghamDataSpace

Data Fields:

* double\*\* binghamVal_Z

A pointer, pointing to a heap memory which stores the values of Z on interpolation nodes as Divide\<b1\<0, Divide\<b2\<0 and b3=0.

* double\*\* dB1Bingham_Z

A pointer, pointing to a heap memory which stores the values of dZ/db1 on interpolation nodes as Divide\<b1\<0, Divide\<b2\<0 and b3=0.

* double\*\* dB2Bingham_Z

A pointer, pointing to a heap memory which stores the values of dZ/db2 on interpolation nodes as Divide\<b1\<0, Divide\<b2\<0 and b3=0.

* double\*\* binghamVal_Z_x12

A pointer, pointing to a heap memory which stores the values of Z(x1^2) on interpolation nodes as Divide\<b1\<0, Divide\<b2\<0 and b3=0.

* double\*\* dB1Bingham_Z_x12

A pointer, pointing to a heap memory which stores the values of dZ(x1^2)/db1 on interpolation nodes as Divide\<b1\<0, Divide\<b2\<0 and b3=0.

* double\*\* dB2Bingham_Z_x12

A pointer, pointing to a heap memory which stores the values of dZ(x1^2)/db2 on interpolation nodes as Divide\<b1\<0, Divide\<b2\<0 and b3=0.

* double\*\* binghamVal_Q_x14

A pointer, pointing to a heap memory which stores the values of \<x1^4\> on interpolation nodes as Divide\<b1\<0, Divide\<b2\<0 and b3=0.

* double\*\* dB1Bingham_Q_x14

A pointer, pointing to a heap memory which stores the values of d\<x1^4\>/db1 on interpolation nodes as Divide\<b1\<0, Divide\<b2\<0 and b3=0.

* double\*\* dB2Bingham_Q_x14

A pointer, pointing to a heap memory which stores the values of d\<x1^4\>/db2 on interpolation nodes as Divide\<b1\<0, Divide\<b2\<0 and b3=0.

* double\*\* binghamVal_Q_x12x22

A pointer, pointing to a heap memory which stores the values of \<x1^2*x2^2\> on interpolation nodes as Divide\<b1\<0, Divide\<b2\<0 and b3=0.

* double\*\* dB1Bingham_Q_x12x22

A pointer, pointing to a heap memory which stores the values of d\<x1^2*x2^2\>/db1 on interpolation nodes as Divide\<b1\<0, Divide\<b2\<0 and b3=0.

* double\*\* dB2Bingham_Q_x12x22

A pointer, pointing to a heap memory which stores the values of d\<x1^2*x2^2\>/db2 on interpolation nodes as Divide\<b1\<0, Divide\<b2\<0 and b3=0.

* double\*\* finite1D_x0

A pointer, pointing to a heap memory which stores the values of the derivatives of g(x2^0|b2,x1) on interpolation nodes as Divide\<b\<0. finite1D_x0[n][m] = d^{2n} g(x2^0|b_m,x1) / d x1^{2n} at x1=0.

* double\*\* finite1D_x2

A pointer, pointing to a heap memory which stores the values of the derivatives of g(x2^2|b2,x1) on interpolation nodes as Divide\<b\<0. finite1D_x0[n][m] = d^{2n} g(x2^2|b_m,x1) / d x1^{2n} at x1=0.

* double\*\* finite1D_x4

A pointer, pointing to a heap memory which stores the values of the derivatives of g(x4^0|b2,x1) on interpolation nodes as Divide\<b\<0. finite1D_x0[n][m] = d^{2n} g(x4^0|b_m,x1) / d x1^{2n} at x1=0.

> Moments

Data Fields:

* double Z

value of Z

* double Z_x12

value of Z(x1^2)

* double Z_x22

value of Z(x2^2)

* double Q_x12

value of \<x1^2\>

* double Q_x22

value of \<x2^2\>

* double Q_x12x22

value of \<x1^2*x1^2\>

* double Q_x14

value of \<x1^4\>

* double Q_x24

value of \<x2^4\>

## [Function Documentation](#function-documentation)

> double hermiteInterpolate(

> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; double &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; areaSpacing,		
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; double &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; x1, 	
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; double &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; x2, 	
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; double &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; y1, 	
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; double &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; y2, 	
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; double &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; f11, 	
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; double &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; f12, 	
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; double &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; f21, 	
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; double &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; f22, 	
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; double &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; f11dx, 	
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; double &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; f12dx, 	
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; double &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; f21dx, 	
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; double &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; f22dx, 	
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; double &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; f11dy, 	
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; double &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; f12dy, 	
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; double &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; f21dy, 	
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; double &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; f22dy, 	
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; double &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; x, 	
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; double &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; y

>)

implement 2D interpolation. Receive the values of a function and its derivatives along both x and y axis at four vertexes of a grid [x1, x2]×[y1, y2], and return the interpolation value at point (x, y).

Parameters:

* areaSpacing

the side length of the grid, whchi is equal to x2-x1=y2-y1.

* x1

the lower bound of the grid along x axis.

* x2

the up bound of the grid along x axis.

* y1

the lower bound of the grid along y axis.

* y2

the up bound of the grid along y axis.

* f11

function value at (x1, y1).

* f12

function value at (x2, y1).

* f21

function value at (x1, y2).

* f22

function value at (x2, y2).

* f11dx

derivative of the function along x at (x1, y1).

* f12dx

derivative of the function along x at (x1, y2).

* f21dx

derivative of the function along x at (x2, y1).

* f22dx

derivative of the function along x at (x2, y2).

* f11dy

derivative of the function along y at (x1, y1).

* f12dy

derivative of the function along y at (x1, y2).

* f21dy

derivative of the function along y at (x2, y1).

* f22dy

derivative of the function along y at (x2, y2).

* x

x component of the point on which the function value is desired. x1 <= x <= x2.

* y

y component of the point on which the function value is desired. y1 <= y <= y2.

> double calcDerivative(

> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; double &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; x,	
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; int &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; order,	
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; double** &nbsp; finite1D

> )

1D interpolation function, for calculating the derivatives of g in case 2.2 in the referenced paper.

Parameters:

* x

the desired point, which is equal to b2 in case 2.2 in the referenced paper.

* order

the order of the derivative of g on x1.

* finite1D

a pointer, pointing to a heap memory which stores the values of each order derivatives at each interpolation nodes. 

> void sort(

> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; double &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; (*array)[3],		
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; int &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; rowNum

> )

sort the first row of a n × 3 array in an increasing order. Corresponding elements in other rows move with the elements in the first row, making column vectors remain unchanged.

Parameters:

* (\*array)[3]

a pointer pointing to the array to be ordered.

* rowNum

row number of the matrix.

> void binghamDataFree(

> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; double &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; binghamDataSpace* binghamDataPointer

> )

release the heap memory allocated for the data structure "binghamDataSpace".

Parameters:

* binghamDataPointer

a pointer pointing to the structure "binghamDataSpace".

> binghamDataSpace* initiateBingham()

allocate heap memory for the data structure "binghamDataSpace" and read in the data from binghamData.bi

> Moments bingham(

> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; double* &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; b,	
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; binghamDataSpace* binghamDataPointer

> )

calculate the moments of Bingham distribution.

Parameters:

* b

a array containing the three eigenvalues of matrix B, which is equal to vector b in the introduction.

* binghamDataPointer

a pointer pointing to the structure "binghamDataSpace".

