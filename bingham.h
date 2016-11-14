#define PI 3.14159265359

/* the dividing number, which is denoted as d in the paper */
#define Divide (-30)
/* total number of 2D interpolation nodes along x(and y) axis, in calculating Z, Z(x1^2) and Z(x2^2) */
#define InterpNodesNumAxis_2D_Z 1201
/* total number of 2D interpolation nodes along x(and y) axis, in calculating <x1^2>, <x2^2>, <x1^2*x2^2>, <x1^4> and <x2^4> */
#define InterpNodesNumAxis_2D_Q 301
/* grid size for 2D interpolation, in calculating Z, Z(x1^2) and Z(x2^2). Its half and square is also difined to avoid repeated calculations. */
#define AreaSpacing_Z 0.025
#define AreaSpacingHalf_Z 0.0125
#define AreaSpacingSquare_Z 0.000625
/* grid size for 2D interpolation, in calculating <x1^2>, <x2^2>, <x1^2*x2^2>, <x1^4> and <x2^4>. Its half and square is also difined to avoid repeated calculations. */
#define AreaSpacing_Q 0.1
#define AreaSpacingHalf_Q 0.05
#define AreaSpacingSquare_Q 0.01
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

/* 2D interpolation function */
double hermiteInterpolate(double areaSpacing, double areaSpacingHalf, double areaSpacingSquare, double x1, double x2, double y1, double y2, double f11, double f12, double f21, double f22, double f11dB1, double f12dB1, double f21dB1, double f22dB1, double f11dB2, double f12dB2, double f21dB2, double f22dB2, double x, double y){
	double sepaX1=x-x1, sepaX2=x-x2;
	double squareX1=sepaX2*sepaX2/areaSpacingSquare, squareX2=sepaX1*sepaX1/areaSpacingSquare;
	double betaX1=sepaX1*squareX1, betaX2=sepaX2*squareX2;
	double alphaX1=betaX1/areaSpacingHalf+squareX1, alphaX2=-betaX2/areaSpacingHalf+squareX2;

	double sepaY1=y-y1, sepaY2=y-y2;
	double squareY1=sepaY2*sepaY2/areaSpacingSquare, squareY2=sepaY1*sepaY1/areaSpacingSquare;
	double betaY1=sepaY1*squareY1, betaY2=sepaY2*squareY2;
	double alphaY1=betaY1/areaSpacingHalf+squareY1, alphaY2=-betaY2/areaSpacingHalf+squareY2;

	double fd=f11*alphaX1+f21*alphaX2+f11dB1*betaX1+f21dB1*betaX2;
	double fu=f12*alphaX1+f22*alphaX2+f12dB1*betaX1+f22dB1*betaX2;
	double fl=f11*alphaY1+f12*alphaY2+f11dB2*betaY1+f12dB2*betaY2;
	double fr=f21*alphaY1+f22*alphaY2+f21dB2*betaY1+f22dB2*betaY2;

	double fddB2=(f21dB2-f11dB2)/areaSpacing*sepaX1+f11dB2;
	double fudB2=(f22dB2-f12dB2)/areaSpacing*sepaX1+f12dB2;
	double fldB1=(f12dB1-f11dB1)/areaSpacing*sepaY1+f11dB1;
	double frdB1=(f22dB1-f21dB1)/areaSpacing*sepaY1+f21dB1;

	double rowVal=fl*alphaX1+fr*alphaX2+fldB1*betaX1+frdB1*betaX2;
	double colVal=fd*alphaY1+fu*alphaY2+fddB2*betaY1+fudB2*betaY2;

	return (rowVal+colVal)/2;
}

/* 1D interpolation function */
double calcDerivative(double x, int order, double** finite1D){
	int index=(int) (-x/LinearSpacing);
	order=order/2;
	return finite1D[order][index]+(-x-index*LinearSpacing)*(finite1D[order][index+1]-finite1D[order][index])/LinearSpacing;
}

/* sort the input array b in a increasing order and record its original order */
void sort(double (*array)[3], int rowNum){
	int i,j;
	int min_idx;
	double temp;
    for(j=0; j<2; j++){
		min_idx = j;
		for(i=j+1; i<3; i++){
			if(array[0][i] < array[0][min_idx]) min_idx = i;
		}
		for(i=0; i<rowNum; i++){
			temp=array[i][j]; array[i][j]=array[i][min_idx]; array[i][min_idx]=temp;
		}
	}
}


/* release the heap memory allocated for the data structure. */
void binghamDataFree(binghamDataSpace* binghamDataPointer){
	int i;

	if(!binghamDataPointer) return;

	if(binghamDataPointer->binghamVal_Z){
		for(i=0;i<InterpNodesNumAxis_2D_Z;i++){
			free((void *)binghamDataPointer->binghamVal_Z[i]);
		}
	}
	free((void *)binghamDataPointer->binghamVal_Z);

	if(binghamDataPointer->dB1Bingham_Z){
		for(i=0;i<InterpNodesNumAxis_2D_Z;i++){
			free((void *)binghamDataPointer->dB1Bingham_Z[i]);
		}
	}
	free((void *)binghamDataPointer->dB1Bingham_Z);

	if(binghamDataPointer->dB2Bingham_Z){
		for(i=0;i<InterpNodesNumAxis_2D_Z;i++){
			free((void *)binghamDataPointer->dB2Bingham_Z[i]);
		}
	}
	free((void *)binghamDataPointer->dB2Bingham_Z);


	if(binghamDataPointer->binghamVal_Z_x12){
		for(i=0;i<InterpNodesNumAxis_2D_Z;i++){
			free((void *)binghamDataPointer->binghamVal_Z_x12[i]);
		}
	}
	free((void *)binghamDataPointer->binghamVal_Z_x12);

	if(binghamDataPointer->dB1Bingham_Z_x12){
		for(i=0;i<InterpNodesNumAxis_2D_Z;i++){
			free((void *)binghamDataPointer->dB1Bingham_Z_x12[i]);
		}
	}
	free((void *)binghamDataPointer->dB1Bingham_Z_x12);

	if(binghamDataPointer->dB2Bingham_Z_x12){
		for(i=0;i<InterpNodesNumAxis_2D_Z;i++){
			free((void *)binghamDataPointer->dB2Bingham_Z_x12[i]);
		}
	}
	free((void *)binghamDataPointer->dB2Bingham_Z_x12);


	if(binghamDataPointer->binghamVal_Q_x14){
		for(i=0;i<InterpNodesNumAxis_2D_Q;i++){
			free((void *)binghamDataPointer->binghamVal_Q_x14[i]);
		}
	}
	free((void *)binghamDataPointer->binghamVal_Q_x14);

	if(binghamDataPointer->dB1Bingham_Q_x14){
		for(i=0;i<InterpNodesNumAxis_2D_Q;i++){
			free((void *)binghamDataPointer->dB1Bingham_Q_x14[i]);
		}
	}
	free((void *)binghamDataPointer->dB1Bingham_Q_x14);

	if(binghamDataPointer->dB2Bingham_Q_x14){
		for(i=0;i<InterpNodesNumAxis_2D_Q;i++){
			free((void *)binghamDataPointer->dB2Bingham_Q_x14[i]);
		}
	}
	free((void *)binghamDataPointer->dB2Bingham_Q_x14);


	if(binghamDataPointer->binghamVal_Q_x12x22){
		for(i=0;i<InterpNodesNumAxis_2D_Q;i++){
			free((void *)binghamDataPointer->binghamVal_Q_x12x22[i]);
		}
	}
	free((void *)binghamDataPointer->binghamVal_Q_x12x22);

	if(binghamDataPointer->dB1Bingham_Q_x12x22){
		for(i=0;i<InterpNodesNumAxis_2D_Q;i++){
			free((void *)binghamDataPointer->dB1Bingham_Q_x12x22[i]);
		}
	}
	free((void *)binghamDataPointer->dB1Bingham_Q_x12x22);

	if(binghamDataPointer->dB2Bingham_Q_x12x22){
		for(i=0;i<InterpNodesNumAxis_2D_Q;i++){
			free((void *)binghamDataPointer->dB2Bingham_Q_x12x22[i]);
		}
	}
	free((void *)binghamDataPointer->dB2Bingham_Q_x12x22);


	if(binghamDataPointer->finite1D_x0){
		for(i=0;i<7;i++){
			free((void *)binghamDataPointer->finite1D_x0[i]);
		}
	}
	free((void *)binghamDataPointer->finite1D_x0);

	if(binghamDataPointer->finite1D_x2){
		for(i=0;i<7;i++){
			free((void *)binghamDataPointer->finite1D_x2[i]);
		}
	}
	free((void *)binghamDataPointer->finite1D_x2);

	if(binghamDataPointer->finite1D_x4){
		for(i=0;i<7;i++){
			free((void *)binghamDataPointer->finite1D_x4[i]);
		}
	}
	free((void *)binghamDataPointer->finite1D_x4);


	free((void *)binghamDataPointer);

	return ;
}

/* allocate heap memory for the data structure and read in the data from binghamData.bi */
binghamDataSpace* initiateBingham(){
	FILE* fp;
	int i;

	/* allocate heap memory for the data structure pointer */
	binghamDataSpace* binghamDataPointer = (binghamDataSpace*) malloc(sizeof(binghamDataSpace));
	if(binghamDataPointer==NULL){
		printf("fail to initiate bingham data space: menmory not enough!\n");
		exit(1);
	}

	/* initialize the pointers */
	binghamDataPointer->binghamVal_Z=NULL;
	binghamDataPointer->dB1Bingham_Z=NULL;
	binghamDataPointer->dB2Bingham_Z=NULL;
	binghamDataPointer->binghamVal_Z_x12=NULL;
	binghamDataPointer->dB1Bingham_Z_x12=NULL;
	binghamDataPointer->dB2Bingham_Z_x12=NULL;
	binghamDataPointer->binghamVal_Q_x14=NULL;
	binghamDataPointer->dB1Bingham_Q_x14=NULL;
	binghamDataPointer->dB2Bingham_Q_x14=NULL;
	binghamDataPointer->binghamVal_Q_x12x22=NULL;
	binghamDataPointer->dB1Bingham_Q_x12x22=NULL;
	binghamDataPointer->dB2Bingham_Q_x12x22=NULL;
	binghamDataPointer->finite1D_x0=NULL;
	binghamDataPointer->finite1D_x2=NULL;
	binghamDataPointer->finite1D_x4=NULL;

	/* open binghamData.bi */
	if((fp=fopen("binghamData.bi","rb+"))==NULL){
		binghamDataFree(binghamDataPointer);
		printf("Cannot open file strike any key exit!");
		exit(1);
	}

	/* allocate heap memory for each pointer and read in the data. */
	binghamDataPointer->binghamVal_Z = (double **) malloc(InterpNodesNumAxis_2D_Z * sizeof(double * ));
	if(binghamDataPointer->binghamVal_Z==NULL){
		binghamDataFree(binghamDataPointer);
		printf("fail to initiate bingham data space: menmory not enough!\n");
		exit(1);
	}
	for(i=0;i<InterpNodesNumAxis_2D_Z;i++) binghamDataPointer->binghamVal_Z[i] = NULL;
	for(i=0;i<InterpNodesNumAxis_2D_Z;i++){
		binghamDataPointer->binghamVal_Z[i] = (double*) malloc(InterpNodesNumAxis_2D_Z * sizeof(double));
		if(binghamDataPointer->binghamVal_Z[i]==NULL){
			binghamDataFree(binghamDataPointer);
			printf("fail to initiate bingham data space: menmory not enough!\n");
			exit(1);
		}
	}
	
	for(i=0; i<InterpNodesNumAxis_2D_Z; i++){
		fread(binghamDataPointer->binghamVal_Z[i], sizeof(double), InterpNodesNumAxis_2D_Z, fp);
	}

	binghamDataPointer->dB1Bingham_Z = (double **) malloc(InterpNodesNumAxis_2D_Z * sizeof(double * ));
	if(binghamDataPointer->dB1Bingham_Z==NULL){
		binghamDataFree(binghamDataPointer);
		printf("fail to initiate bingham data space: menmory not enough!\n");
		exit(1);
	}
	for(i=0;i<InterpNodesNumAxis_2D_Z;i++) binghamDataPointer->dB1Bingham_Z[i] = NULL;
	for(i=0;i<InterpNodesNumAxis_2D_Z;i++){
		binghamDataPointer->dB1Bingham_Z[i] = (double*) malloc(InterpNodesNumAxis_2D_Z * sizeof(double));
		if(binghamDataPointer->dB1Bingham_Z[i]==NULL){
			binghamDataFree(binghamDataPointer);
			printf("fail to initiate bingham data space: menmory not enough!\n");
			exit(1);
		}
	}
	
	for(i=0; i<InterpNodesNumAxis_2D_Z; i++){
		fread(binghamDataPointer->dB1Bingham_Z[i], sizeof(double), InterpNodesNumAxis_2D_Z, fp);
	}

	binghamDataPointer->dB2Bingham_Z = (double **) malloc(InterpNodesNumAxis_2D_Z * sizeof(double * ));
	if(binghamDataPointer->dB2Bingham_Z==NULL){
		binghamDataFree(binghamDataPointer);
		printf("fail to initiate bingham data space: menmory not enough!\n");
		exit(1);
	}
	for(i=0;i<InterpNodesNumAxis_2D_Z;i++) binghamDataPointer->dB2Bingham_Z[i] = NULL;
	for(i=0;i<InterpNodesNumAxis_2D_Z;i++){
		binghamDataPointer->dB2Bingham_Z[i] = (double*) malloc(InterpNodesNumAxis_2D_Z * sizeof(double));
		if(binghamDataPointer->dB2Bingham_Z[i]==NULL){
			binghamDataFree(binghamDataPointer);
			printf("fail to initiate bingham data space: menmory not enough!\n");
			exit(1);
		}
	}
	
	for(i=0; i<InterpNodesNumAxis_2D_Z; i++){
		fread(binghamDataPointer->dB2Bingham_Z[i], sizeof(double), InterpNodesNumAxis_2D_Z, fp);
	}


	binghamDataPointer->binghamVal_Z_x12 = (double **) malloc(InterpNodesNumAxis_2D_Z * sizeof(double * ));
	if(binghamDataPointer->binghamVal_Z_x12==NULL){
		binghamDataFree(binghamDataPointer);
		printf("fail to initiate bingham data space: menmory not enough!\n");
		exit(1);
	}
	for(i=0;i<InterpNodesNumAxis_2D_Z;i++) binghamDataPointer->binghamVal_Z_x12[i] = NULL;
	for(i=0;i<InterpNodesNumAxis_2D_Z;i++){
		binghamDataPointer->binghamVal_Z_x12[i] = (double*) malloc(InterpNodesNumAxis_2D_Z * sizeof(double));
		if(binghamDataPointer->binghamVal_Z_x12[i]==NULL){
			binghamDataFree(binghamDataPointer);
			printf("fail to initiate bingham data space: menmory not enough!\n");
			exit(1);
		}
	}
	
	for(i=0; i<InterpNodesNumAxis_2D_Z; i++){
		fread(binghamDataPointer->binghamVal_Z_x12[i], sizeof(double), InterpNodesNumAxis_2D_Z, fp);
	}

	binghamDataPointer->dB1Bingham_Z_x12 = (double **) malloc(InterpNodesNumAxis_2D_Z * sizeof(double * ));
	if(binghamDataPointer->dB1Bingham_Z_x12==NULL){
		binghamDataFree(binghamDataPointer);
		printf("fail to initiate bingham data space: menmory not enough!\n");
		exit(1);
	}
	for(i=0;i<InterpNodesNumAxis_2D_Z;i++) binghamDataPointer->dB1Bingham_Z_x12[i] = NULL;
	for(i=0;i<InterpNodesNumAxis_2D_Z;i++){
		binghamDataPointer->dB1Bingham_Z_x12[i] = (double*) malloc(InterpNodesNumAxis_2D_Z * sizeof(double));
		if(binghamDataPointer->dB1Bingham_Z_x12[i]==NULL){
			binghamDataFree(binghamDataPointer);
			printf("fail to initiate bingham data space: menmory not enough!\n");
			exit(1);
		}
	}
	
	for(i=0; i<InterpNodesNumAxis_2D_Z; i++){
		fread(binghamDataPointer->dB1Bingham_Z_x12[i], sizeof(double), InterpNodesNumAxis_2D_Z, fp);
	}

	binghamDataPointer->dB2Bingham_Z_x12 = (double **) malloc(InterpNodesNumAxis_2D_Z * sizeof(double * ));
	if(binghamDataPointer->dB2Bingham_Z_x12==NULL){
		binghamDataFree(binghamDataPointer);
		printf("fail to initiate bingham data space: menmory not enough!\n");
		exit(1);
	}
	for(i=0;i<InterpNodesNumAxis_2D_Z;i++) binghamDataPointer->dB2Bingham_Z_x12[i] = NULL;
	for(i=0;i<InterpNodesNumAxis_2D_Z;i++){
		binghamDataPointer->dB2Bingham_Z_x12[i] = (double*) malloc(InterpNodesNumAxis_2D_Z * sizeof(double));
		if(binghamDataPointer->dB2Bingham_Z_x12[i]==NULL){
			binghamDataFree(binghamDataPointer);
			printf("fail to initiate bingham data space: menmory not enough!\n");
			exit(1);
		}
	}
	
	for(i=0; i<InterpNodesNumAxis_2D_Z; i++){
		fread(binghamDataPointer->dB2Bingham_Z_x12[i], sizeof(double), InterpNodesNumAxis_2D_Z, fp);
	}


	binghamDataPointer->binghamVal_Q_x14 = (double **) malloc(InterpNodesNumAxis_2D_Q * sizeof(double * ));
	if(binghamDataPointer->binghamVal_Q_x14==NULL){
		binghamDataFree(binghamDataPointer);
		printf("fail to initiate bingham data space: menmory not enough!\n");
		exit(1);
	}
	for(i=0;i<InterpNodesNumAxis_2D_Q;i++) binghamDataPointer->binghamVal_Q_x14[i] = NULL;
	for(i=0;i<InterpNodesNumAxis_2D_Q;i++){
		binghamDataPointer->binghamVal_Q_x14[i] = (double*) malloc(InterpNodesNumAxis_2D_Q * sizeof(double));
		if(binghamDataPointer->binghamVal_Q_x14[i]==NULL){
			binghamDataFree(binghamDataPointer);
			printf("fail to initiate bingham data space: menmory not enough!\n");
			exit(1);
		}
	}
	
	for(i=0; i<InterpNodesNumAxis_2D_Q; i++){
		fread(binghamDataPointer->binghamVal_Q_x14[i], sizeof(double), InterpNodesNumAxis_2D_Q, fp);
	}

	binghamDataPointer->dB1Bingham_Q_x14 = (double **) malloc(InterpNodesNumAxis_2D_Q * sizeof(double * ));
	if(binghamDataPointer->dB1Bingham_Q_x14==NULL){
		binghamDataFree(binghamDataPointer);
		printf("fail to initiate bingham data space: menmory not enough!\n");
		exit(1);
	}
	for(i=0;i<InterpNodesNumAxis_2D_Q;i++) binghamDataPointer->dB1Bingham_Q_x14[i] = NULL;
	for(i=0;i<InterpNodesNumAxis_2D_Q;i++){
		binghamDataPointer->dB1Bingham_Q_x14[i] = (double*) malloc(InterpNodesNumAxis_2D_Q * sizeof(double));
		if(binghamDataPointer->dB1Bingham_Q_x14[i]==NULL){
			binghamDataFree(binghamDataPointer);
			printf("fail to initiate bingham data space: menmory not enough!\n");
			exit(1);
		}
	}
	
	for(i=0; i<InterpNodesNumAxis_2D_Q; i++){
		fread(binghamDataPointer->dB1Bingham_Q_x14[i], sizeof(double), InterpNodesNumAxis_2D_Q, fp);
	}

	binghamDataPointer->dB2Bingham_Q_x14 = (double **) malloc(InterpNodesNumAxis_2D_Q * sizeof(double * ));
	if(binghamDataPointer->dB2Bingham_Q_x14==NULL){
		binghamDataFree(binghamDataPointer);
		printf("fail to initiate bingham data space: menmory not enough!\n");
		exit(1);
	}
	for(i=0;i<InterpNodesNumAxis_2D_Q;i++) binghamDataPointer->dB2Bingham_Q_x14[i] = NULL;
	for(i=0;i<InterpNodesNumAxis_2D_Q;i++){
		binghamDataPointer->dB2Bingham_Q_x14[i] = (double*) malloc(InterpNodesNumAxis_2D_Q * sizeof(double));
		if(binghamDataPointer->dB2Bingham_Q_x14[i]==NULL){
			binghamDataFree(binghamDataPointer);
			printf("fail to initiate bingham data space: menmory not enough!\n");
			exit(1);
		}
	}
	
	for(i=0; i<InterpNodesNumAxis_2D_Q; i++){
		fread(binghamDataPointer->dB2Bingham_Q_x14[i], sizeof(double), InterpNodesNumAxis_2D_Q, fp);
	}


	binghamDataPointer->binghamVal_Q_x12x22 = (double **) malloc(InterpNodesNumAxis_2D_Q * sizeof(double * ));
	if(binghamDataPointer->binghamVal_Q_x12x22==NULL){
		binghamDataFree(binghamDataPointer);
		printf("fail to initiate bingham data space: menmory not enough!\n");
		exit(1);
	}
	for(i=0;i<InterpNodesNumAxis_2D_Q;i++) binghamDataPointer->binghamVal_Q_x12x22[i] = NULL;
	for(i=0;i<InterpNodesNumAxis_2D_Q;i++){
		binghamDataPointer->binghamVal_Q_x12x22[i] = (double*) malloc(InterpNodesNumAxis_2D_Q * sizeof(double));
		if(binghamDataPointer->binghamVal_Q_x12x22[i]==NULL){
			binghamDataFree(binghamDataPointer);
			printf("fail to initiate bingham data space: menmory not enough!\n");
			exit(1);
		}
	}
	
	for(i=0; i<InterpNodesNumAxis_2D_Q; i++){
		fread(binghamDataPointer->binghamVal_Q_x12x22[i], sizeof(double), InterpNodesNumAxis_2D_Q, fp);
	}

	binghamDataPointer->dB1Bingham_Q_x12x22 = (double **) malloc(InterpNodesNumAxis_2D_Q * sizeof(double * ));
	if(binghamDataPointer->dB1Bingham_Q_x12x22==NULL){
		binghamDataFree(binghamDataPointer);
		printf("fail to initiate bingham data space: menmory not enough!\n");
		exit(1);
	}
	for(i=0;i<InterpNodesNumAxis_2D_Q;i++) binghamDataPointer->dB1Bingham_Q_x12x22[i] = NULL;
	for(i=0;i<InterpNodesNumAxis_2D_Q;i++){
		binghamDataPointer->dB1Bingham_Q_x12x22[i] = (double*) malloc(InterpNodesNumAxis_2D_Q * sizeof(double));
		if(binghamDataPointer->dB1Bingham_Q_x12x22[i]==NULL){
			binghamDataFree(binghamDataPointer);
			printf("fail to initiate bingham data space: menmory not enough!\n");
			exit(1);
		}
	}
	
	for(i=0; i<InterpNodesNumAxis_2D_Q; i++){
		fread(binghamDataPointer->dB1Bingham_Q_x12x22[i], sizeof(double), InterpNodesNumAxis_2D_Q, fp);
	}

	binghamDataPointer->dB2Bingham_Q_x12x22 = (double **) malloc(InterpNodesNumAxis_2D_Q * sizeof(double * ));
	if(binghamDataPointer->dB2Bingham_Q_x12x22==NULL){
		binghamDataFree(binghamDataPointer);
		printf("fail to initiate bingham data space: menmory not enough!\n");
		exit(1);
	}
	for(i=0;i<InterpNodesNumAxis_2D_Q;i++) binghamDataPointer->dB2Bingham_Q_x12x22[i] = NULL;
	for(i=0;i<InterpNodesNumAxis_2D_Q;i++){
		binghamDataPointer->dB2Bingham_Q_x12x22[i] = (double*) malloc(InterpNodesNumAxis_2D_Q * sizeof(double));
		if(binghamDataPointer->dB2Bingham_Q_x12x22[i]==NULL){
			binghamDataFree(binghamDataPointer);
			printf("fail to initiate bingham data space: menmory not enough!\n");
			exit(1);
		}
	}
	
	for(i=0; i<InterpNodesNumAxis_2D_Q; i++){
		fread(binghamDataPointer->dB2Bingham_Q_x12x22[i], sizeof(double), InterpNodesNumAxis_2D_Q, fp);
	}


	binghamDataPointer->finite1D_x0 = (double **) malloc(7 * sizeof(double * ));
	if(binghamDataPointer->finite1D_x0==NULL){
		binghamDataFree(binghamDataPointer);
		printf("fail to initiate bingham data space: menmory not enough!\n");
		exit(1);
	}
	for(i=0;i<7;i++) binghamDataPointer->finite1D_x0[i] = NULL;
	for(i=0;i<7;i++){
		binghamDataPointer->finite1D_x0[i] = (double*) malloc(InterpNodesNumAxis_1D * sizeof(double));
		if(binghamDataPointer->finite1D_x0[i]==NULL){
			binghamDataFree(binghamDataPointer);
			printf("fail to initiate bingham data space: menmory not enough!\n");
			exit(1);
		}
	}
	
	for(i=0; i<7; i++){
		fread(binghamDataPointer->finite1D_x0[i], sizeof(double), InterpNodesNumAxis_1D, fp);
	}


	binghamDataPointer->finite1D_x2 = (double **) malloc(7 * sizeof(double * ));
	if(binghamDataPointer->finite1D_x2==NULL){
		binghamDataFree(binghamDataPointer);
		printf("fail to initiate bingham data space: menmory not enough!\n");
		exit(1);
	}
	for(i=0;i<7;i++) binghamDataPointer->finite1D_x2[i] = NULL;
	for(i=0;i<7;i++){
		binghamDataPointer->finite1D_x2[i] = (double*) malloc(InterpNodesNumAxis_1D * sizeof(double));
		if(binghamDataPointer->finite1D_x2[i]==NULL){
			binghamDataFree(binghamDataPointer);
			printf("fail to initiate bingham data space: menmory not enough!\n");
			exit(1);
		}
	}
	
	for(i=0; i<7; i++){
		fread(binghamDataPointer->finite1D_x2[i], sizeof(double), InterpNodesNumAxis_1D, fp);
	}


	binghamDataPointer->finite1D_x4 = (double **) malloc(7 * sizeof(double * ));
	if(binghamDataPointer->finite1D_x4==NULL){
		binghamDataFree(binghamDataPointer);
		printf("fail to initiate bingham data space: menmory not enough!\n");
		exit(1);
	}
	for(i=0;i<7;i++) binghamDataPointer->finite1D_x4[i] = NULL;
	for(i=0;i<7;i++){
		binghamDataPointer->finite1D_x4[i] = (double*) malloc(InterpNodesNumAxis_1D * sizeof(double));
		if(binghamDataPointer->finite1D_x4[i]==NULL){
			binghamDataFree(binghamDataPointer);
			printf("fail to initiate bingham data space: menmory not enough!\n");
			exit(1);
		}
	}
	
	for(i=0; i<7; i++){
		fread(binghamDataPointer->finite1D_x4[i], sizeof(double), InterpNodesNumAxis_1D, fp);
	}


	fclose(fp);

	return binghamDataPointer;
}


/* function for calculating the moments of Bingham distribution */
double* bingham(double* b, binghamDataSpace* binghamDataPointer, double* Q){
	/* result variables */
	double Z, Z_x12, Z_x22, Q_x12, Q_x22, Q_x14, Q_x24, Q_x12x22;
	/* a matrix for sorting the input array b and store its original order */
	double sortMat[2][3]={{b[0], b[1], b[2]}, {1, 2, 3}};
	/* a array stored the original order of b */
	double BIndex[3];

	/* sort the input array b in a increasing order and record its original order in BIndex */
	sort(sortMat, 2);
	BIndex[0]=sortMat[1][0]; BIndex[1]=sortMat[1][1]; BIndex[2]=sortMat[1][2];

	double b1=sortMat[0][0];
	double b2=sortMat[0][1];
	double b3=sortMat[0][2];

	/* translate the array to make b3=0 */
	double scale=exp(b3);
	b1=b1-b3;
	b2=b2-b3;
	b3=0;

	/* indexes of grid which the target point located in, in calculating the 2D interpolation in section 2.3 in the paper.*/
	int indexB1, indexB2;
	/* Some denotes used in calculation to avoid repeated calculations */
	double b1Cor, coCor, x1, x2;
	if(b1>Divide){
	}
	else if(b2>Divide){
		b1Cor=2*sqrt(-PI/b1);
	}
	else{
		coCor=-2*PI/sqrt(b1*b2);
		x1=1/b1, x2=1/b2;
	}

	/* calculating Z */
	/* calculating in case 2.3 in the paper. */
	if(b1>Divide){
		indexB1=-(int) (b1/AreaSpacing_Z); indexB2=-(int) (b2/AreaSpacing_Z);
		Z=hermiteInterpolate(AreaSpacing_Z, AreaSpacingHalf_Z, AreaSpacingSquare_Z, -indexB1*AreaSpacing_Z-AreaSpacing_Z, -indexB1*AreaSpacing_Z, -indexB2*AreaSpacing_Z-AreaSpacing_Z, -indexB2*AreaSpacing_Z, binghamDataPointer->binghamVal_Z[indexB2+1][indexB1+1], binghamDataPointer->binghamVal_Z[indexB2][indexB1+1], binghamDataPointer->binghamVal_Z[indexB2+1][indexB1], binghamDataPointer->binghamVal_Z[indexB2][indexB1], binghamDataPointer->dB1Bingham_Z[indexB2+1][indexB1+1], binghamDataPointer->dB1Bingham_Z[indexB2][indexB1+1], binghamDataPointer->dB1Bingham_Z[indexB2+1][indexB1], binghamDataPointer->dB1Bingham_Z[indexB2][indexB1], binghamDataPointer->dB2Bingham_Z[indexB2+1][indexB1+1], binghamDataPointer->dB2Bingham_Z[indexB2][indexB1+1], binghamDataPointer->dB2Bingham_Z[indexB2+1][indexB1], binghamDataPointer->dB2Bingham_Z[indexB2][indexB1], b1, b2);
    }
	/* calculating in case 2.2 in the paper. */
	else if (b2>Divide){
		Z=b1Cor*(calcDerivative(b2, 0, binghamDataPointer->finite1D_x0) - calcDerivative(b2, 2, binghamDataPointer->finite1D_x0)/4/b1 + calcDerivative(b2, 4, binghamDataPointer->finite1D_x0)/32/b1/b1 - calcDerivative(b2, 6, binghamDataPointer->finite1D_x0)/384/b1/b1/b1 + calcDerivative(b2, 8, binghamDataPointer->finite1D_x0)/6144/b1/b1/b1/b1 - calcDerivative(b2, 10, binghamDataPointer->finite1D_x0)/122880/b1/b1/b1/b1/b1 + calcDerivative(b2, 12, binghamDataPointer->finite1D_x0)/2949120/b1/b1/b1/b1/b1/b1);
	}
	/* calculating in case 2.1 in the paper. */
	else {
		Z=coCor*(x1*(x1*(x1*(x1*(7.26745605469*x1 + 4.03747558594*x2 - 1.79443359375) + x2*(3.46069335938*x2 - 1.025390625) + 0.5859375) + x2*(x2*(3.46069335938*x2 - 0.9228515625) + 0.3515625) - 0.28125) + x2*(x2*(x2*(4.03747558594*x2 - 1.025390625) + 0.3515625) - 0.1875) + 0.25) + x2*(x2*(x2*(x2*(7.26745605469*x2 - 1.79443359375) + 0.5859375) - 0.28125) + 0.25) - 1.0);
	}

	/* calculating Z(x1^2) and Z(x1^2) */
	/* calculating in case 2.3 in the paper. */
	if(b1>Divide){
		indexB1=-(int) (b1/AreaSpacing_Z); indexB2=-(int) (b2/AreaSpacing_Z);
		Z_x12=hermiteInterpolate(AreaSpacing_Z, AreaSpacingHalf_Z, AreaSpacingSquare_Z, -indexB1*AreaSpacing_Z-AreaSpacing_Z, -indexB1*AreaSpacing_Z, -indexB2*AreaSpacing_Z-AreaSpacing_Z, -indexB2*AreaSpacing_Z, binghamDataPointer->binghamVal_Z_x12[indexB2+1][indexB1+1], binghamDataPointer->binghamVal_Z_x12[indexB2][indexB1+1], binghamDataPointer->binghamVal_Z_x12[indexB2+1][indexB1], binghamDataPointer->binghamVal_Z_x12[indexB2][indexB1], binghamDataPointer->dB1Bingham_Z_x12[indexB2+1][indexB1+1], binghamDataPointer->dB1Bingham_Z_x12[indexB2][indexB1+1], binghamDataPointer->dB1Bingham_Z_x12[indexB2+1][indexB1], binghamDataPointer->dB1Bingham_Z_x12[indexB2][indexB1], binghamDataPointer->dB2Bingham_Z_x12[indexB2+1][indexB1+1], binghamDataPointer->dB2Bingham_Z_x12[indexB2][indexB1+1], binghamDataPointer->dB2Bingham_Z_x12[indexB2+1][indexB1], binghamDataPointer->dB2Bingham_Z_x12[indexB2][indexB1], b1, b2);
		Z_x22=hermiteInterpolate(AreaSpacing_Z, AreaSpacingHalf_Z, AreaSpacingSquare_Z, -indexB2*AreaSpacing_Z-AreaSpacing_Z, -indexB2*AreaSpacing_Z, -indexB1*AreaSpacing_Z-AreaSpacing_Z, -indexB1*AreaSpacing_Z, binghamDataPointer->binghamVal_Z_x12[indexB1+1][indexB2+1], binghamDataPointer->binghamVal_Z_x12[indexB1][indexB2+1], binghamDataPointer->binghamVal_Z_x12[indexB1+1][indexB2], binghamDataPointer->binghamVal_Z_x12[indexB1][indexB2], binghamDataPointer->dB1Bingham_Z_x12[indexB1+1][indexB2+1], binghamDataPointer->dB1Bingham_Z_x12[indexB1][indexB2+1], binghamDataPointer->dB1Bingham_Z_x12[indexB1+1][indexB2], binghamDataPointer->dB1Bingham_Z_x12[indexB1][indexB2], binghamDataPointer->dB2Bingham_Z_x12[indexB1+1][indexB2+1], binghamDataPointer->dB2Bingham_Z_x12[indexB1][indexB2+1], binghamDataPointer->dB2Bingham_Z_x12[indexB1+1][indexB2], binghamDataPointer->dB2Bingham_Z_x12[indexB1][indexB2], b2, b1);
    }
	/* calculating in case 2.2 in the paper. */
	else if (b2>Divide){
		Z_x12=b1Cor*(-calcDerivative(b2, 0, binghamDataPointer->finite1D_x0)/2/b1 + 0.375*calcDerivative(b2, 2, binghamDataPointer->finite1D_x0)/b1/b1 - 0.078125*calcDerivative(b2, 4, binghamDataPointer->finite1D_x0)/b1/b1/b1 + 0.00911458*calcDerivative(b2, 6, binghamDataPointer->finite1D_x0)/b1/b1/b1/b1 - 0.000732422*calcDerivative(b2, 8, binghamDataPointer->finite1D_x0)/b1/b1/b1/b1/b1 + 0.0000447591*calcDerivative(b2, 10, binghamDataPointer->finite1D_x0)/b1/b1/b1/b1/b1/b1);
		Z_x22=b1Cor*(calcDerivative(b2, 0, binghamDataPointer->finite1D_x2) - calcDerivative(b2, 2, binghamDataPointer->finite1D_x2)/4/b1 + calcDerivative(b2, 4, binghamDataPointer->finite1D_x2)/32/b1/b1 - calcDerivative(b2, 6, binghamDataPointer->finite1D_x2)/384/b1/b1/b1 + calcDerivative(b2, 8, binghamDataPointer->finite1D_x2)/6144/b1/b1/b1/b1 - calcDerivative(b2, 10, binghamDataPointer->finite1D_x2)/122880/b1/b1/b1/b1/b1 + calcDerivative(b2, 12, binghamDataPointer->finite1D_x2)/2949120/b1/b1/b1/b1/b1/b1);
	}
	/* calculating in case 2.1 in the paper. */
	else {
		Z_x12=coCor*(-x1*(x1*(x1*(x1*(x1*(39.9710083008*x1 + 18.1686401367*x2 - 8.07495117188) + x2*(12.1124267578*x2 - 3.5888671875) + 2.05078125) + x2*(x2*(8.65173339844*x2 - 2.30712890625) + 0.87890625) - 0.703125) + x2*(x2*(x2*(6.05621337891*x2 - 1.5380859375) + 0.52734375) - 0.28125) + 0.375) + x2*(x2*(x2*(x2*(3.63372802734*x2 - 0.897216796875) + 0.29296875) - 0.140625) + 0.125) - 0.5));
		Z_x22=coCor*(-x2*(x2*(x2*(x2*(x2*(39.9710083008*x2 + 18.1686401367*x1 - 8.07495117188) + x1*(12.1124267578*x1 - 3.5888671875) + 2.05078125) + x1*(x1*(8.65173339844*x1 - 2.30712890625) + 0.87890625) - 0.703125) + x1*(x1*(x1*(6.05621337891*x1 - 1.5380859375) + 0.52734375) - 0.28125) + 0.375) + x1*(x1*(x1*(x1*(3.63372802734*x1 - 0.897216796875) + 0.29296875) - 0.140625) + 0.125) - 0.5));
	}

	/* calculating <x1^2> and <x2^2> */
	Q_x12=Z_x12/Z;
	Q_x22=Z_x22/Z;

	/* calculating <x1^4> and <x2^4> */
	/* calculating in case 2.3 in the paper. */
	if(b1>Divide){
		indexB1=-(int) (b1/AreaSpacing_Q); indexB2=-(int) (b2/AreaSpacing_Q);
		Q_x14=hermiteInterpolate(AreaSpacing_Q, AreaSpacingHalf_Q, AreaSpacingSquare_Q, -indexB1*AreaSpacing_Q-AreaSpacing_Q, -indexB1*AreaSpacing_Q, -indexB2*AreaSpacing_Q-AreaSpacing_Q, -indexB2*AreaSpacing_Q, binghamDataPointer->binghamVal_Q_x14[indexB2+1][indexB1+1], binghamDataPointer->binghamVal_Q_x14[indexB2][indexB1+1], binghamDataPointer->binghamVal_Q_x14[indexB2+1][indexB1], binghamDataPointer->binghamVal_Q_x14[indexB2][indexB1], binghamDataPointer->dB1Bingham_Q_x14[indexB2+1][indexB1+1], binghamDataPointer->dB1Bingham_Q_x14[indexB2][indexB1+1], binghamDataPointer->dB1Bingham_Q_x14[indexB2+1][indexB1], binghamDataPointer->dB1Bingham_Q_x14[indexB2][indexB1], binghamDataPointer->dB2Bingham_Q_x14[indexB2+1][indexB1+1], binghamDataPointer->dB2Bingham_Q_x14[indexB2][indexB1+1], binghamDataPointer->dB2Bingham_Q_x14[indexB2+1][indexB1], binghamDataPointer->dB2Bingham_Q_x14[indexB2][indexB1], b1, b2);
		Q_x24=hermiteInterpolate(AreaSpacing_Q, AreaSpacingHalf_Q, AreaSpacingSquare_Q, -indexB2*AreaSpacing_Q-AreaSpacing_Q, -indexB2*AreaSpacing_Q, -indexB1*AreaSpacing_Q-AreaSpacing_Q, -indexB1*AreaSpacing_Q, binghamDataPointer->binghamVal_Q_x14[indexB1+1][indexB2+1], binghamDataPointer->binghamVal_Q_x14[indexB1][indexB2+1], binghamDataPointer->binghamVal_Q_x14[indexB1+1][indexB2], binghamDataPointer->binghamVal_Q_x14[indexB1][indexB2], binghamDataPointer->dB1Bingham_Q_x14[indexB1+1][indexB2+1], binghamDataPointer->dB1Bingham_Q_x14[indexB1][indexB2+1], binghamDataPointer->dB1Bingham_Q_x14[indexB1+1][indexB2], binghamDataPointer->dB1Bingham_Q_x14[indexB1][indexB2], binghamDataPointer->dB2Bingham_Q_x14[indexB1+1][indexB2+1], binghamDataPointer->dB2Bingham_Q_x14[indexB1][indexB2+1], binghamDataPointer->dB2Bingham_Q_x14[indexB1+1][indexB2], binghamDataPointer->dB2Bingham_Q_x14[indexB1][indexB2], b2, b1);
    }
	/* calculating in case 2.2 in the paper. */
	else if (b2>Divide){
		Q_x14=b1Cor*(0.75*calcDerivative(b2, 0, binghamDataPointer->finite1D_x0)/b1/b1 - 0.9375*calcDerivative(b2, 2, binghamDataPointer->finite1D_x0)/b1/b1/b1 + 0.273438*calcDerivative(b2, 4, binghamDataPointer->finite1D_x0)/b1/b1/b1/b1 - 0.0410156*calcDerivative(b2, 6, binghamDataPointer->finite1D_x0)/b1/b1/b1/b1/b1 + 0.00402832*calcDerivative(b2, 8, binghamDataPointer->finite1D_x0)/b1/b1/b1/b1/b1/b1)/Z;
		Q_x24=b1Cor*(calcDerivative(b2, 0, binghamDataPointer->finite1D_x4) - calcDerivative(b2, 2, binghamDataPointer->finite1D_x4)/4/b1 + calcDerivative(b2, 4, binghamDataPointer->finite1D_x4)/32/b1/b1 - calcDerivative(b2, 6, binghamDataPointer->finite1D_x4)/384/b1/b1/b1 + calcDerivative(b2, 8, binghamDataPointer->finite1D_x4)/6144/b1/b1/b1/b1 - calcDerivative(b2, 10, binghamDataPointer->finite1D_x4)/122880/b1/b1/b1/b1/b1 + calcDerivative(b2, 12, binghamDataPointer->finite1D_x4)/2949120/b1/b1/b1/b1/b1/b1)/Z;
	}
	/* calculating in case 2.1 in the paper. */
	else {
		Q_x14=coCor*(-x1*x1*(x1*(x1*(x1*(44.412231445312*x1 + 16.14990234375*x2 - 9.228515625) + x2*(8.074951171875*x2 - 3.076171875) + 2.4609375) + x2*(x2*(3.84521484375*x2 - 1.318359375) + 0.703125) - 0.9375) + x2*(x2*(x2*(1.3458251953125*x2 - 0.439453125) + 0.2109375) - 0.1875) + 0.75))/Z;
		Q_x24=coCor*(-x2*x2*(x2*(x2*(x2*(16.14990234375*x1 + 44.412231445312*x2 - 9.228515625) + x1*(8.074951171875*x1 - 3.076171875) + 2.4609375) + x1*(x1*(3.84521484375*x1 - 1.318359375) + 0.703125) - 0.9375) + x1*(x1*(x1*(1.3458251953125*x1 - 0.439453125) + 0.2109375) - 0.1875) + 0.75))/Z;
	}

	/* calculating <x1^2*x2^2> */
	/* calculating in case 2.3 in the paper. */
	if(b1>Divide){
		indexB1=-(int) (b1/AreaSpacing_Q); indexB2=-(int) (b2/AreaSpacing_Q);
		Q_x12x22=hermiteInterpolate(AreaSpacing_Q, AreaSpacingHalf_Q, AreaSpacingSquare_Q, -indexB1*AreaSpacing_Q-AreaSpacing_Q, -indexB1*AreaSpacing_Q, -indexB2*AreaSpacing_Q-AreaSpacing_Q, -indexB2*AreaSpacing_Q, binghamDataPointer->binghamVal_Q_x12x22[indexB2+1][indexB1+1], binghamDataPointer->binghamVal_Q_x12x22[indexB2][indexB1+1], binghamDataPointer->binghamVal_Q_x12x22[indexB2+1][indexB1], binghamDataPointer->binghamVal_Q_x12x22[indexB2][indexB1], binghamDataPointer->dB1Bingham_Q_x12x22[indexB2+1][indexB1+1], binghamDataPointer->dB1Bingham_Q_x12x22[indexB2][indexB1+1], binghamDataPointer->dB1Bingham_Q_x12x22[indexB2+1][indexB1], binghamDataPointer->dB1Bingham_Q_x12x22[indexB2][indexB1], binghamDataPointer->dB2Bingham_Q_x12x22[indexB2+1][indexB1+1], binghamDataPointer->dB2Bingham_Q_x12x22[indexB2][indexB1+1], binghamDataPointer->dB2Bingham_Q_x12x22[indexB2+1][indexB1], binghamDataPointer->dB2Bingham_Q_x12x22[indexB2][indexB1], b1, b2);
    }
	/* calculating in case 2.2 in the paper. */
	else if (b2>Divide){
		Q_x12x22=b1Cor*(-calcDerivative(b2, 0, binghamDataPointer->finite1D_x2)/2/b1 + 0.375*calcDerivative(b2, 2, binghamDataPointer->finite1D_x2)/b1/b1 - 0.078125*calcDerivative(b2, 4, binghamDataPointer->finite1D_x2)/b1/b1/b1 + 0.00911458*calcDerivative(b2, 6, binghamDataPointer->finite1D_x2)/b1/b1/b1/b1 - 0.000732422*calcDerivative(b2, 8, binghamDataPointer->finite1D_x2)/b1/b1/b1/b1/b1 + 0.0000447591*calcDerivative(b2, 10, binghamDataPointer->finite1D_x2)/b1/b1/b1/b1/b1/b1)/Z;
	}
	/* calculating in case 2.1 in the paper. */
	else {
		Q_x12x22=coCor*(-x1*(x1*(x1*(x2*(x2*(5.767822265625*x2 - 1.318359375) + 0.3515625) + x1*(4.0374755859375*x1*x2 + x2*(5.38330078125*x2 - 1.025390625))) + x2*(x2*(x2*(5.38330078125*x2 - 1.318359375) + 0.421875) - 0.1875)) + x2*(x2*(x2*(x2*(4.0374755859375*x2 - 1.025390625) + 0.3515625) - 0.1875) + 0.25)))/Z;
	}


	/* reorder the array Q to make it correspond to the original order of b, the final order of Q is {Z, Z(x1^2), Z(x2^2), <x1^2>, <x2^2>, <x1^2*x2^2>, <x1^4>, <x2^4>} */
	if(BIndex[0]==1 && BIndex[1]==2){
		Q[5]=Q_x12x22;
	}
	else if(BIndex[0]==2 && BIndex[1]==1){
		Q[5]=Q_x12x22;
	}
	else if(BIndex[0]==1 && BIndex[1]==3){
		Q[5]=- Q_x14 - Q_x12x22 + Q_x12;
	}
	else if(BIndex[0]==3 && BIndex[1]==1){
		Q[5]=- Q_x12x22 - Q_x24 + Q_x22;
	}
	else if(BIndex[0]==2 && BIndex[1]==3){
		Q[5]=- Q_x14 - Q_x12x22 + Q_x12;
	}
	else if(BIndex[0]==3 && BIndex[1]==2){
		Q[5]=- Q_x12x22 - Q_x24 + Q_x22;
	}

	double resortMat[4][3]={
		{BIndex[0], BIndex[1], BIndex[2]},
		{Z_x12, Z_x22, Z - Z_x12 - Z_x22},
		{Q_x12, Q_x22, 1 - Q_x12 - Q_x22},
		{Q_x14, Q_x24, Q_x14 + 2*Q_x12x22 - 2*Q_x12 + Q_x24 - 2*Q_x22 + 1},
	};
	sort(resortMat, 4);

	Q[0]=scale*Z;
	Q[1]=scale*resortMat[1][0]; Q[2]=scale*resortMat[1][1];
	Q[3]=resortMat[2][0]; Q[4]=resortMat[2][1];
	Q[6]=resortMat[3][0]; Q[7]=resortMat[3][1];

    return Q;

}
