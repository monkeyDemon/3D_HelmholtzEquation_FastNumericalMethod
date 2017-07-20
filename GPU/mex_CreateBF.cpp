#include "mex.h"
#include "AddVectors.h"

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
	if(nrhs!=2)
		mexErrMsgTxt("Invaid number of input arguments");
	if(nlhs!=1)
		mexErrMsgTxt("Invaid number of output arguments");
	if(!mxIsSingle(prhs[0])&&!mxIsSingle(prhs[1]))
		mexErrMsgTxt("input vector data type must be single");

	int numRowsA=(int)mxGetM(prhs[0]);
	int numColsA=(int)mxGetN(prhs[0]);
	int numRowsB=(int)mxGetM(prhs[1]);
	int numColsB=(int)mxGetN(prhs[1]);

	if(numRowsA!=numRowsB || numColsA!=numColsB)
		mexErrMsgTxt("Invalid size. The sizes of two vectors must be same");

	int minSize=(numRowsA<numColsA)?numRowsA:numColsA;
	int maxSize=(numRowsA>numColsA)?numRowsA:numColsA;

	if(minSize!=1)
		mexErrMsgTxt("Invalid size. The vector must be one dimentional");

	float *A=(float*)mxGetData(prhs[0]);
	float *B=(float*)mxGetData(prhs[1]);

	plhs[0]=mxCreateNumericMatrix(numRowsA,numColsB,mxSINGLE_CLASS,mxREAL);
	float *C=(float*)mxGetData(plhs[0]);

	addVectors(A,B,C,maxSize);
}
