#include "mex.h"
#include "prepare.h"

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
	if(nrhs!=1)
		mexErrMsgTxt("Invaid number of input arguments");
	if(nlhs!=1)
		mexErrMsgTxt("Invaid number of output arguments");

	float *inputs=(float*)mxGetData(prhs[0]);
	float size=inputs[0];
	float h=inputs[1];
	float K0=inputs[2];

	//mexPrintf("%f\n",size);
	//mexPrintf("%f\n",h);
	//mexPrintf("%f\n",K0);

	plhs[0]=mxCreateNumericMatrix(size*size*size,1,mxSINGLE_CLASS,mxREAL);
	float *BU=(float*)mxGetData(plhs[0]);

	computeBU(BU,(int)size,h,K0);
}
