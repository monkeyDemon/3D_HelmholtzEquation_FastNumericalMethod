#include<stdio.h>
#include<stdlib.h>
#include "prepare.h"

using namespace std;

int main()
{
	int size=127;
	float h=0.0078125;
	float K0=0;

	float *BU;
	BU= (float*) malloc(size*size*size*sizeof(float));
	
	computeBU(BU,size,h,K0);	
}
