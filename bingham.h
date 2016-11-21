#ifndef BINGHAM_H
#define BINGHAM_H

#define PI 3.14159265359

/* the dividing number, which is denoted as d in the paper */
#define Divide (-30)
/* total number of 2D interpolation nodes along x(and y) axis, in calculating Z, Z(x1^2) and Z(x2^2) */
#define InterpNodesNumAxis_2D_Z 1201
/* total number of 2D interpolation nodes along x(and y) axis, in calculating <x1^2>, <x2^2>, <x1^2*x2^2>, <x1^4> and <x2^4> */
#define InterpNodesNumAxis_2D_Q 301
/* grid size for 2D interpolation, in calculating Z, Z(x1^2) and Z(x2^2) */
#define AreaSpacing_Z 0.025
/* grid size for 2D interpolation, in calculating <x1^2>, <x2^2>, <x1^2*x2^2>, <x1^4> and <x2^4>. */
#define AreaSpacing_Q 0.1
/* total number of 1D interpolation nodes, in calculating the derivatives of g in case 2.2 in the paper. */
#define InterpNodesNumAxis_1D 30001
/* grid size for 1D interpolation, in calculating the derivatives of g in case 2.2 in the paper. */
#define LinearSpacing 0.001

/* Difine the data structure for storing the interpolation data. */
typedef struct binghamDataSpace{
	double** binghamVal_Z;
	double** dB1Bingham_Z;
	double** dB2Bingham_Z;

	double** binghamVal_Z_x12;
	double** dB1Bingham_Z_x12;
	double** dB2Bingham_Z_x12;

	double** binghamVal_Q_x14;
	double** dB1Bingham_Q_x14;
	double** dB2Bingham_Q_x14;

	double** binghamVal_Q_x12x22;
	double** dB1Bingham_Q_x12x22;
	double** dB2Bingham_Q_x12x22;

	double** finite1D_x0;
	double** finite1D_x2;
	double** finite1D_x4;
} binghamDataSpace;

/* Difine the data structure of the result. */
typedef struct Moments{
	double Z;
	double Z_x12;
	double Z_x22;
	double Q_x12;
	double Q_x22;
	double Q_x12x22;
	double Q_x14;
	double Q_x24;
} Moments;

/* 2D interpolation function */
double hermiteInterpolate(double areaSpacing, double x1, double x2, double y1, double y2, double f11, double f12, double f21, double f22, double f11dx, double f12dx, double f21dx, double f22dx, double f11dy, double f12dy, double f21dy, double f22dy, double x, double y);

/* 1D interpolation function */
double calcDerivative(double x, int order, double** finite1D);

/* sort the input array b in a increasing order and record its original order */
void sort(double (*array)[3], int rowNum);


/* release the heap memory allocated for the data structure. */
void binghamDataFree(binghamDataSpace* binghamDataPointer);

/* allocate heap memory for the data structure and read in the data from binghamData.bi */
binghamDataSpace* initiateBingham();


/* function for calculating the moments of Bingham distribution */
Moments bingham(double* b, binghamDataSpace* binghamDataPointer);

#endif