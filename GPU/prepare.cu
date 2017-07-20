#include "prepare.h"
#include "stdio.h"
#include "math.h"
#define PI 3.14159265358979323846

//void computeBU_Z_array(float *BU,int i,int j,int M,int N,int K,float h,float K0);

__host__ __device__ float compute_BoundaryCondition(int M,float h,int i,int j,int k)
{
	float Xstart=0;
	float Ystart=0;
	float Zstart=0;
	
	float x=Xstart+h*i;
	float y=Ystart+h*j;
	float z=Zstart+h*k;
       
	float value=0;	
	
	if(i==0)
	{
	        value=sin(PI*y)*sin(PI*z);
		//printf("x=0,y=%f,z=%f,value:%f\n",y,z,value);
	}
	else if(i==M+1)
	{
	        value=2*sin(PI*y)*sin(PI*z);
		//printf("x=1,y=%f,z=%f,value:%f\n",y,z,value);
	}
	else
	        value=0;
	
	return value;
}


__global__ void computeBU_kernel(float* BU, int size, float h, float K0)
{
	int i=blockIdx.x;
	int j=threadIdx.x;
	
	//printf("%d  %d\n",i,j);

	if(i>=size||j>=size)
		return;
	
	i++;
	j++;	
	
	//if(i==1)
	//	printf("i=%d\n",i);
	//if(i==128)
	//	printf("i=%d\n",i);
	//if(i==129)
	//	printf("!!!i=%d\n",i);

	//computeBU_Z_array(BU,i,j,size,size,size,h,K0);

	int M=size;
	int N=size;
	int K=size;

	int Idx=(i-1)*N*K+ (j-1)*K -1;

	float para1=(1+h*h*K0*K0/12)/(h*h);
	float para2=1/(6*h*h);

	float temp=0;	
	float temp_gama=0;
	
	if( i==1 && j==1  )
	{
		for(int k=1;k<=K;k++)
		{
			if(k==1)
			{
				temp=compute_BoundaryCondition(M,h,0,j,k)+compute_BoundaryCondition(M,h,i,0,k)+compute_BoundaryCondition(M,h,i,j,0);
				BU[Idx+k]=BU[Idx+k]+temp*para1;
				temp=compute_BoundaryCondition(M,h,0,0,k)+compute_BoundaryCondition(M,h,0,j,0)-4*compute_BoundaryCondition(M,h,0,j,k)+compute_BoundaryCondition(M,h,0,j,k+1)+compute_BoundaryCondition(M,h,0,j+1,k)+compute_BoundaryCondition(M,h,i,0,0)-4*compute_BoundaryCondition(M,h,i,0,k)+compute_BoundaryCondition(M,h,i,0,k+1)-4*compute_BoundaryCondition(M,h,i,j,0)+compute_BoundaryCondition(M,h,i,j+1,0)+compute_BoundaryCondition(M,h,i+1,0,k)+compute_BoundaryCondition(M,h,i+1,j,0);
				BU[Idx+k]=BU[Idx+k]+temp*para2;
			}
			else if(k==K)
			{
				temp=compute_BoundaryCondition(M,h,0,j,k)+compute_BoundaryCondition(M,h,i,0,k);
				temp_gama=compute_BoundaryCondition(M,h,i,j,K+1);
				BU[Idx+k]=BU[Idx+k]+temp*para1;
				BU[Idx+k]=BU[Idx+k]+temp_gama*para1;
				temp=compute_BoundaryCondition(M,h,0,0,k)+compute_BoundaryCondition(M,h,0,j,k-1)-4*compute_BoundaryCondition(M,h,0,j,k)+compute_BoundaryCondition(M,h,0,j,k+1)+compute_BoundaryCondition(M,h,0,j+1,k)+compute_BoundaryCondition(M,h,i,0,k-1)-4*compute_BoundaryCondition(M,h,i,0,k)+compute_BoundaryCondition(M,h,i,0,k+1)+compute_BoundaryCondition(M,h,i+1,0,k);
				temp_gama=-4*compute_BoundaryCondition(M,h,i,j,k+1)+compute_BoundaryCondition(M,h,i,j+1,k+1)+compute_BoundaryCondition(M,h,i+1,j,k+1);
				BU[Idx+k]=BU[Idx+k]+temp*para2;
				BU[Idx+k]=BU[Idx+k]+temp_gama*para2;
			}
			else
			{
				temp=compute_BoundaryCondition(M,h,0,j,k)+compute_BoundaryCondition(M,h,i,0,k);
				BU[Idx+k]=BU[Idx+k]+temp*para1;
				temp=compute_BoundaryCondition(M,h,0,0,k)+compute_BoundaryCondition(M,h,0,j,k-1)-4*compute_BoundaryCondition(M,h,0,j,k)+compute_BoundaryCondition(M,h,0,j,k+1)+compute_BoundaryCondition(M,h,0,j+1,k)+compute_BoundaryCondition(M,h,i,0,k-1)-4*compute_BoundaryCondition(M,h,i,0,k)+compute_BoundaryCondition(M,h,i,0,k+1)+compute_BoundaryCondition(M,h,i+1,0,k);
				BU[Idx+k]=BU[Idx+k]+temp*para2;
			}
		}
	}
	else if( i==1 && j==N  )
	{
		for(int k=1;k<=K;k++)
		{
			if(k==1)
			{
				temp=compute_BoundaryCondition(M,h,i-1,j,k)+compute_BoundaryCondition(M,h,i,j+1,k)+compute_BoundaryCondition(M,h,i,j,k-1);
				BU[Idx+k]=BU[Idx+k]+temp*para1;
				temp=compute_BoundaryCondition(M,h,i-1,j-1,k)+compute_BoundaryCondition(M,h,i-1,j,k-1)-4*compute_BoundaryCondition(M,h,i-1,j,k)+compute_BoundaryCondition(M,h,i-1,j,k+1)+compute_BoundaryCondition(M,h,i-1,j+1,k)+compute_BoundaryCondition(M,h,i,j-1,k-1)-4*compute_BoundaryCondition(M,h,i,j,k-1)+compute_BoundaryCondition(M,h,i,j+1,k-1)-4*compute_BoundaryCondition(M,h,i,j+1,k)+compute_BoundaryCondition(M,h,i,j+1,k+1)+compute_BoundaryCondition(M,h,i+1,j,k-1)+compute_BoundaryCondition(M,h,i+1,j+1,k);
				BU[Idx+k]=BU[Idx+k]+temp*para2;
			}
			else if(k==K)
			{
				temp=compute_BoundaryCondition(M,h,i-1,j,k)+compute_BoundaryCondition(M,h,i,j+1,k);
				temp_gama=compute_BoundaryCondition(M,h,i,j,k+1);
				BU[Idx+k]=BU[Idx+k]+temp*para1;
				BU[Idx+k]=BU[Idx+k]+temp_gama*para1;
				temp=compute_BoundaryCondition(M,h,i-1,j-1,k)+compute_BoundaryCondition(M,h,i-1,j,k-1)-4*compute_BoundaryCondition(M,h,i-1,j,k)+compute_BoundaryCondition(M,h,i-1,j,k+1)+compute_BoundaryCondition(M,h,i-1,j+1,k)+compute_BoundaryCondition(M,h,i,j+1,k-1)-4*compute_BoundaryCondition(M,h,i,j+1,k)+compute_BoundaryCondition(M,h,i,j+1,k+1)+compute_BoundaryCondition(M,h,i+1,j+1,k);
				temp_gama=compute_BoundaryCondition(M,h,i,j-1,k+1)-4*compute_BoundaryCondition(M,h,i,j,k+1)+compute_BoundaryCondition(M,h,i+1,j,k+1);
				BU[Idx+k]=BU[Idx+k]+temp*para2;
				BU[Idx+k]=BU[Idx+k]+temp_gama*para2;
			}
			else
			{
				temp=compute_BoundaryCondition(M,h,i-1,j,k)+compute_BoundaryCondition(M,h,i,j+1,k);
				BU[Idx+k]=BU[Idx+k]+temp*para1;
				temp=compute_BoundaryCondition(M,h,i-1,j-1,k)+compute_BoundaryCondition(M,h,i-1,j,k-1)-4*compute_BoundaryCondition(M,h,i-1,j,k)+compute_BoundaryCondition(M,h,i-1,j,k+1)+compute_BoundaryCondition(M,h,i-1,j+1,k)+compute_BoundaryCondition(M,h,i,j+1,k-1)-4*compute_BoundaryCondition(M,h,i,j+1,k)+compute_BoundaryCondition(M,h,i,j+1,k+1)+compute_BoundaryCondition(M,h,i+1,j+1,k);
				BU[Idx+k]=BU[Idx+k]+temp*para2;
			}
		}
	}
	else if( i==1 && j>1 && j<N  )
	{
		for(int k=1;k<=K;k++)
		{
			if(k==1)
			{
				temp=compute_BoundaryCondition(M,h,i-1,j,k)+compute_BoundaryCondition(M,h,i,j,k-1);
				BU[Idx+k]=BU[Idx+k]+temp*para1;
				temp=compute_BoundaryCondition(M,h,i-1,j-1,k)+compute_BoundaryCondition(M,h,i-1,j,k-1)-4*compute_BoundaryCondition(M,h,i-1,j,k)+compute_BoundaryCondition(M,h,i-1,j,k+1)+compute_BoundaryCondition(M,h,i-1,j+1,k)+compute_BoundaryCondition(M,h,i,j-1,k-1)-4*compute_BoundaryCondition(M,h,i,j,k-1)+compute_BoundaryCondition(M,h,i,j+1,k-1)+compute_BoundaryCondition(M,h,i+1,j,k-1);
				BU[Idx+k]=BU[Idx+k]+temp*para2;
			}
			else if(k==K)
			{
				temp=compute_BoundaryCondition(M,h,i-1,j,k);
				temp_gama=compute_BoundaryCondition(M,h,i,j,k+1);
				BU[Idx+k]=BU[Idx+k]+temp*para1;
				BU[Idx+k]=BU[Idx+k]+temp_gama*para1;
				temp=compute_BoundaryCondition(M,h,i-1,j-1,k)+compute_BoundaryCondition(M,h,i-1,j,k-1)-4*compute_BoundaryCondition(M,h,i-1,j,k)+compute_BoundaryCondition(M,h,i-1,j,k+1)+compute_BoundaryCondition(M,h,i-1,j+1,k);
				temp_gama=compute_BoundaryCondition(M,h,i,j-1,k+1)-4*compute_BoundaryCondition(M,h,i,j,k+1)+compute_BoundaryCondition(M,h,i,j+1,k+1)+compute_BoundaryCondition(M,h,i+1,j,k+1);
				BU[Idx+k]=BU[Idx+k]+temp*para2;
				BU[Idx+k]=BU[Idx+k]+temp_gama*para2;
			}
			else
			{
				temp=compute_BoundaryCondition(M,h,i-1,j,k);
				BU[Idx+k]=BU[Idx+k]+temp*para1;
				temp=compute_BoundaryCondition(M,h,i-1,j-1,k)+compute_BoundaryCondition(M,h,i-1,j,k-1)-4*compute_BoundaryCondition(M,h,i-1,j,k)+compute_BoundaryCondition(M,h,i-1,j,k+1)+compute_BoundaryCondition(M,h,i-1,j+1,k);
				BU[Idx+k]=BU[Idx+k]+temp*para2;
			}
		}
	}
	else if( i==M && j==1 )
	{
		for(int k=1;k<=K;k++)
		{
			if(k==1)
			{
				temp=compute_BoundaryCondition(M,h,i+1,j,k)+compute_BoundaryCondition(M,h,i,j-1,k)+compute_BoundaryCondition(M,h,i,j,k-1);
				BU[Idx+k]=BU[Idx+k]+temp*para1;
				temp=compute_BoundaryCondition(M,h,i-1,j-1,k)+compute_BoundaryCondition(M,h,i-1,j,k-1)+compute_BoundaryCondition(M,h,i,j-1,k-1)-4*compute_BoundaryCondition(M,h,i,j-1,k)+compute_BoundaryCondition(M,h,i,j-1,k+1)-4*compute_BoundaryCondition(M,h,i,j,k-1)+compute_BoundaryCondition(M,h,i,j+1,k-1)+compute_BoundaryCondition(M,h,i+1,j-1,k)+compute_BoundaryCondition(M,h,i+1,j,k-1)-4*compute_BoundaryCondition(M,h,i+1,j,k)+compute_BoundaryCondition(M,h,i+1,j,k+1)+compute_BoundaryCondition(M,h,i+1,j+1,k);
				BU[Idx+k]=BU[Idx+k]+temp*para2;
			}
			else if(k==K)
			{
				temp=compute_BoundaryCondition(M,h,i+1,j,k)+compute_BoundaryCondition(M,h,i,j-1,k);
				temp_gama=compute_BoundaryCondition(M,h,i,j,k+1);
				BU[Idx+k]=BU[Idx+k]+temp*para1;
				BU[Idx+k]=BU[Idx+k]+temp_gama*para1;
				temp=compute_BoundaryCondition(M,h,i-1,j-1,k)+compute_BoundaryCondition(M,h,i,j-1,k-1)-4*compute_BoundaryCondition(M,h,i,j-1,k)+compute_BoundaryCondition(M,h,i,j-1,k+1)+compute_BoundaryCondition(M,h,i+1,j-1,k)+compute_BoundaryCondition(M,h,i+1,j,k-1)-4*compute_BoundaryCondition(M,h,i+1,j,k)+compute_BoundaryCondition(M,h,i+1,j,k+1)+compute_BoundaryCondition(M,h,i+1,j+1,k);
				temp_gama=compute_BoundaryCondition(M,h,i-1,j,k+1)-4*compute_BoundaryCondition(M,h,i,j,k+1)+compute_BoundaryCondition(M,h,i,j+1,k+1);
				BU[Idx+k]=BU[Idx+k]+temp*para2;
				BU[Idx+k]=BU[Idx+k]+temp_gama*para2;
			}
			else
			{
				temp=compute_BoundaryCondition(M,h,i+1,j,k)+compute_BoundaryCondition(M,h,i,j-1,k);
				BU[Idx+k]=BU[Idx+k]+temp*para1;
				temp=compute_BoundaryCondition(M,h,i-1,j-1,k)+compute_BoundaryCondition(M,h,i,j-1,k-1)-4*compute_BoundaryCondition(M,h,i,j-1,k)+compute_BoundaryCondition(M,h,i,j-1,k+1)+compute_BoundaryCondition(M,h,i+1,j-1,k)+compute_BoundaryCondition(M,h,i+1,j,k-1)-4*compute_BoundaryCondition(M,h,i+1,j,k)+compute_BoundaryCondition(M,h,i+1,j,k+1)+compute_BoundaryCondition(M,h,i+1,j+1,k);
				BU[Idx+k]=BU[Idx+k]+temp*para2;
			}
		}
	}
	else if( i==M && j==N )
	{
		for(int k=1;k<=K;k++)
		{
			if(k==1)
			{
				temp=compute_BoundaryCondition(M,h,i+1,j,k)+compute_BoundaryCondition(M,h,i,j+1,k)+compute_BoundaryCondition(M,h,i,j,k-1);
				BU[Idx+k]=BU[Idx+k]+temp*para1;
				temp=compute_BoundaryCondition(M,h,i-1,j,k-1)+compute_BoundaryCondition(M,h,i-1,j+1,k)+compute_BoundaryCondition(M,h,i,j-1,k-1)-4*compute_BoundaryCondition(M,h,i,j,k-1)+compute_BoundaryCondition(M,h,i,j+1,k-1)-4*compute_BoundaryCondition(M,h,i,j+1,k)+compute_BoundaryCondition(M,h,i,j+1,k+1)+compute_BoundaryCondition(M,h,i+1,j-1,k)+compute_BoundaryCondition(M,h,i+1,j,k-1)-4*compute_BoundaryCondition(M,h,i+1,j,k)+compute_BoundaryCondition(M,h,i+1,j,k+1)+compute_BoundaryCondition(M,h,i+1,j+1,k);
				BU[Idx+k]=BU[Idx+k]+temp*para2;
			}
			else if(k==K)
			{
				temp=compute_BoundaryCondition(M,h,i+1,j,k)+compute_BoundaryCondition(M,h,i,j+1,k);
				temp_gama=compute_BoundaryCondition(M,h,i,j,k+1);
				BU[Idx+k]=BU[Idx+k]+temp*para1;
				BU[Idx+k]=BU[Idx+k]+temp_gama*para1;
				temp=compute_BoundaryCondition(M,h,i-1,j+1,k)+compute_BoundaryCondition(M,h,i,j+1,k-1)-4*compute_BoundaryCondition(M,h,i,j+1,k)+compute_BoundaryCondition(M,h,i,j+1,k+1)+compute_BoundaryCondition(M,h,i+1,j-1,k)+compute_BoundaryCondition(M,h,i+1,j,k-1)-4*compute_BoundaryCondition(M,h,i+1,j,k)+compute_BoundaryCondition(M,h,i+1,j,k+1)+compute_BoundaryCondition(M,h,i+1,j+1,k);
				temp_gama=compute_BoundaryCondition(M,h,i-1,j,k+1)+compute_BoundaryCondition(M,h,i,j-1,k+1)-4*compute_BoundaryCondition(M,h,i,j,k+1);
				BU[Idx+k]=BU[Idx+k]+temp*para2;
				BU[Idx+k]=BU[Idx+k]+temp_gama*para2;
			}
			else
			{
				temp=compute_BoundaryCondition(M,h,i+1,j,k)+compute_BoundaryCondition(M,h,i,j+1,k);
				BU[Idx+k]=BU[Idx+k]+temp*para1;
				temp=compute_BoundaryCondition(M,h,i-1,j+1,k)+compute_BoundaryCondition(M,h,i,j+1,k-1)-4*compute_BoundaryCondition(M,h,i,j+1,k)+compute_BoundaryCondition(M,h,i,j+1,k+1)+compute_BoundaryCondition(M,h,i+1,j-1,k)+compute_BoundaryCondition(M,h,i+1,j,k-1)-4*compute_BoundaryCondition(M,h,i+1,j,k)+compute_BoundaryCondition(M,h,i+1,j,k+1)+compute_BoundaryCondition(M,h,i+1,j+1,k);
				BU[Idx+k]=BU[Idx+k]+temp*para2;
			}
		}
	}
	else if( i==M && j>1 && j<N  )
	{
		for(int k=1;k<=K;k++)
		{
			if(k==1)
			{
				temp=compute_BoundaryCondition(M,h,i+1,j,k)+compute_BoundaryCondition(M,h,i,j,k-1);
				BU[Idx+k]=BU[Idx+k]+temp*para1;
				temp=compute_BoundaryCondition(M,h,i-1,j,k-1)+compute_BoundaryCondition(M,h,i,j-1,k-1)-4*compute_BoundaryCondition(M,h,i,j,k-1)+compute_BoundaryCondition(M,h,i,j+1,k-1)+compute_BoundaryCondition(M,h,i+1,j-1,k)+compute_BoundaryCondition(M,h,i+1,j,k-1)-4*compute_BoundaryCondition(M,h,i+1,j,k)+compute_BoundaryCondition(M,h,i+1,j,k+1)+compute_BoundaryCondition(M,h,i+1,j+1,k);
				BU[Idx+k]=BU[Idx+k]+temp*para2;
			}
			else if(k==K)
			{
				temp=compute_BoundaryCondition(M,h,i+1,j,k);
				temp_gama=compute_BoundaryCondition(M,h,i,j,k+1);
				BU[Idx+k]=BU[Idx+k]+temp*para1;
				BU[Idx+k]=BU[Idx+k]+temp_gama*para1;
				temp=compute_BoundaryCondition(M,h,i+1,j-1,k)+compute_BoundaryCondition(M,h,i+1,j,k-1)-4*compute_BoundaryCondition(M,h,i+1,j,k)+compute_BoundaryCondition(M,h,i+1,j,k+1)+compute_BoundaryCondition(M,h,i+1,j+1,k);
				temp_gama=compute_BoundaryCondition(M,h,i-1,j,k+1)+compute_BoundaryCondition(M,h,i,j-1,k+1)-4*compute_BoundaryCondition(M,h,i,j,k+1)+compute_BoundaryCondition(M,h,i,j+1,k+1);
				BU[Idx+k]=BU[Idx+k]+temp*para2;
				BU[Idx+k]=BU[Idx+k]+temp_gama*para2;
			}
			else
			{
				temp=compute_BoundaryCondition(M,h,i+1,j,k);
				BU[Idx+k]=BU[Idx+k]+temp*para1;
				temp=compute_BoundaryCondition(M,h,i+1,j-1,k)+compute_BoundaryCondition(M,h,i+1,j,k-1)-4*compute_BoundaryCondition(M,h,i+1,j,k)+compute_BoundaryCondition(M,h,i+1,j,k+1)+compute_BoundaryCondition(M,h,i+1,j+1,k);
				BU[Idx+k]=BU[Idx+k]+temp*para2;
			}
		}
	}
	else if( i>1 && i<M && j==1  )
	{
		for(int k=1;k<=K;k++)
		{	
			if(k==1)
			{
				temp=compute_BoundaryCondition(M,h,i,j-1,k)+compute_BoundaryCondition(M,h,i,j,k-1);
				BU[Idx+k]=BU[Idx+k]+temp*para1;
				temp=compute_BoundaryCondition(M,h,i-1,j-1,k)+compute_BoundaryCondition(M,h,i-1,j,k-1)+compute_BoundaryCondition(M,h,i,j-1,k-1)-4*compute_BoundaryCondition(M,h,i,j-1,k)+compute_BoundaryCondition(M,h,i,j-1,k+1)-4*compute_BoundaryCondition(M,h,i,j,k-1)+compute_BoundaryCondition(M,h,i,j+1,k-1)+compute_BoundaryCondition(M,h,i+1,j-1,k)+compute_BoundaryCondition(M,h,i+1,j,k-1);
				BU[Idx+k]=BU[Idx+k]+temp*para2;
			}
			else if(k==K)
			{
				temp=compute_BoundaryCondition(M,h,i,j-1,k);
				temp_gama=compute_BoundaryCondition(M,h,i,j,k+1);
				BU[Idx+k]=BU[Idx+k]+temp*para1;
				BU[Idx+k]=BU[Idx+k]+temp_gama*para1;
				temp=compute_BoundaryCondition(M,h,i-1,j-1,k)+compute_BoundaryCondition(M,h,i,j-1,k-1)-4*compute_BoundaryCondition(M,h,i,j-1,k)+compute_BoundaryCondition(M,h,i,j-1,k+1)+compute_BoundaryCondition(M,h,i+1,j-1,k);
				temp_gama=compute_BoundaryCondition(M,h,i-1,j,k+1)-4*compute_BoundaryCondition(M,h,i,j,k+1)+compute_BoundaryCondition(M,h,i,j+1,k+1)+compute_BoundaryCondition(M,h,i+1,j,k+1);
				BU[Idx+k]=BU[Idx+k]+temp*para2;
				BU[Idx+k]=BU[Idx+k]+temp_gama*para2;
			}
			else
			{
				temp=compute_BoundaryCondition(M,h,i,j-1,k);
				BU[Idx+k]=BU[Idx+k]+temp*para1;
				temp=compute_BoundaryCondition(M,h,i-1,j-1,k)+compute_BoundaryCondition(M,h,i,j-1,k-1)-4*compute_BoundaryCondition(M,h,i,j-1,k)+compute_BoundaryCondition(M,h,i,j-1,k+1)+compute_BoundaryCondition(M,h,i+1,j-1,k);
				BU[Idx+k]=BU[Idx+k]+temp*para2;
			}
		}
	}
	else if( i>1 && i<M && j==N   )
	{
		for(int k=1;k<=K;k++)
		{
			if(k==1)
			{
				temp=compute_BoundaryCondition(M,h,i,j+1,k)+compute_BoundaryCondition(M,h,i,j,k-1);
				BU[Idx+k]=BU[Idx+k]+temp*para1;
				temp=compute_BoundaryCondition(M,h,i-1,j,k-1)+compute_BoundaryCondition(M,h,i-1,j+1,k)+compute_BoundaryCondition(M,h,i,j-1,k-1)-4*compute_BoundaryCondition(M,h,i,j,k-1)+compute_BoundaryCondition(M,h,i,j+1,k-1)-4*compute_BoundaryCondition(M,h,i,j+1,k)+compute_BoundaryCondition(M,h,i,j+1,k+1)+compute_BoundaryCondition(M,h,i+1,j,k-1)+compute_BoundaryCondition(M,h,i+1,j+1,k);
				BU[Idx+k]=BU[Idx+k]+temp*para2;
			}
			else if(k==K)
			{
				temp=compute_BoundaryCondition(M,h,i,j+1,k);
				temp_gama=compute_BoundaryCondition(M,h,i,j,k+1);
				BU[Idx+k]=BU[Idx+k]+temp*para1;
				BU[Idx+k]=BU[Idx+k]+temp_gama*para1;
				temp=compute_BoundaryCondition(M,h,i-1,j+1,k)+compute_BoundaryCondition(M,h,i,j+1,k-1)-4*compute_BoundaryCondition(M,h,i,j+1,k)+compute_BoundaryCondition(M,h,i,j+1,k+1)+compute_BoundaryCondition(M,h,i+1,j+1,k);
				temp_gama=compute_BoundaryCondition(M,h,i-1,j,k+1)+compute_BoundaryCondition(M,h,i,j-1,k+1)-4*compute_BoundaryCondition(M,h,i,j,k+1)+compute_BoundaryCondition(M,h,i+1,j,k+1);
				BU[Idx+k]=BU[Idx+k]+temp*para2;
				BU[Idx+k]=BU[Idx+k]+temp_gama*para2;
			}
			else
			{
				temp=compute_BoundaryCondition(M,h,i,j+1,k);
				BU[Idx+k]=BU[Idx+k]+temp*para1;
				temp=compute_BoundaryCondition(M,h,i-1,j+1,k)+compute_BoundaryCondition(M,h,i,j+1,k-1)-4*compute_BoundaryCondition(M,h,i,j+1,k)+compute_BoundaryCondition(M,h,i,j+1,k+1)+compute_BoundaryCondition(M,h,i+1,j+1,k);
				BU[Idx+k]=BU[Idx+k]+temp*para2;
			}
		}
	}
	else if( i>1 && i<M && j>1 && j<N )
	{
		for(int k=1;k<=K;k++)
		{
			if(k==1)
			{
				temp=compute_BoundaryCondition(M,h,i,j,k-1);
				BU[Idx+k]=BU[Idx+k]+temp*para1;
				temp=compute_BoundaryCondition(M,h,i-1,j,k-1)+compute_BoundaryCondition(M,h,i,j-1,k-1)-4*compute_BoundaryCondition(M,h,i,j,k-1)+compute_BoundaryCondition(M,h,i,j+1,k-1)+compute_BoundaryCondition(M,h,i+1,j,k-1);
				BU[Idx+k]=BU[Idx+k]+temp*para2;
			}
			else if(k==K)
			{
				temp_gama=compute_BoundaryCondition(M,h,i,j,k+1);
				BU[Idx+k]=BU[Idx+k]+temp_gama*para1;
				temp_gama=compute_BoundaryCondition(M,h,i-1,j,k+1)+compute_BoundaryCondition(M,h,i,j-1,k+1)-4*compute_BoundaryCondition(M,h,i,j,k+1)+compute_BoundaryCondition(M,h,i,j+1,k+1)+compute_BoundaryCondition(M,h,i+1,j,k+1);
				BU[Idx+k]=BU[Idx+k]+temp_gama*para2;
			}
			else
			{
				// do nothing
			}
		}
	}
}


void computeBU(float* BU, int size, float h, float K0)
{
	float *dev_BU=0;

	cudaMalloc(&dev_BU,sizeof(float)*size*size*size);

	cudaMemset(dev_BU,0,sizeof(float)*size*size*size);

	computeBU_kernel<<<size,size>>>(dev_BU,size,h,K0);

	cudaMemcpy(BU,dev_BU,sizeof(float)*size*size*size,cudaMemcpyDeviceToHost);

	cudaFree(dev_BU);
}


//__host__ __device__ void computeBU_Z_array(float *BU,int i,int j,int M,int N,int K,float h,float K0)
//{
//	int Idx=(i-1)*N*K+ (j-1)*K -1;
//
//	float para1=(1+h*h*K0*K0/12)/(h*h);
//	float para2=1/(6*h*h);
//
//	float temp=0;	
//	float temp_gama=0;
//	
//	if( i==1 && j==1  )
//	{
//		for(int k=1;k<=K;k++)
//		{
//			if(k==1)
//			{
//				temp=compute_BoundaryCondition(M,h,0,j,k)+compute_BoundaryCondition(M,h,i,0,k)+compute_BoundaryCondition(M,h,i,j,0);
//				BU[Idx+k]=BU[Idx+k]+temp*para1;
//				temp=compute_BoundaryCondition(M,h,0,0,k)+compute_BoundaryCondition(M,h,0,j,0)-4*compute_BoundaryCondition(M,h,0,j,k)+compute_BoundaryCondition(M,h,0,j,k+1)+compute_BoundaryCondition(M,h,0,j+1,k)+compute_BoundaryCondition(M,h,i,0,0)-4*compute_BoundaryCondition(M,h,i,0,k)+compute_BoundaryCondition(M,h,i,0,k+1)-4*compute_BoundaryCondition(M,h,i,j,0)+compute_BoundaryCondition(M,h,i,j+1,0)+compute_BoundaryCondition(M,h,i+1,0,k)+compute_BoundaryCondition(M,h,i+1,j,0);
//				BU[Idx+k]=BU[Idx+k]+temp*para2;
//			}
//			else if(k==K)
//			{
//				temp=compute_BoundaryCondition(M,h,0,j,k)+compute_BoundaryCondition(M,h,i,0,k);
//				temp_gama=compute_BoundaryCondition(M,h,i,j,K+1);
//				BU[Idx+k]=BU[Idx+k]+temp*para1;
//				BU[Idx+k]=BU[Idx+k]+temp_gama*para1;
//				temp=compute_BoundaryCondition(M,h,0,0,k)+compute_BoundaryCondition(M,h,0,j,k-1)-4*compute_BoundaryCondition(M,h,0,j,k)+compute_BoundaryCondition(M,h,0,j,k+1)+compute_BoundaryCondition(M,h,0,j+1,k)+compute_BoundaryCondition(M,h,i,0,k-1)-4*compute_BoundaryCondition(M,h,i,0,k)+compute_BoundaryCondition(M,h,i,0,k+1)+compute_BoundaryCondition(M,h,i+1,0,k);
//				temp_gama=-4*compute_BoundaryCondition(M,h,i,j,k+1)+compute_BoundaryCondition(M,h,i,j+1,k+1)+compute_BoundaryCondition(M,h,i+1,j,k+1);
//				BU[Idx+k]=BU[Idx+k]+temp*para2;
//				BU[Idx+k]=BU[Idx+k]+temp_gama*para2;
//			}
//			else
//			{
//				temp=compute_BoundaryCondition(M,h,0,j,k)+compute_BoundaryCondition(M,h,i,0,k);
//				BU[Idx+k]=BU[Idx+k]+temp*para1;
//				temp=compute_BoundaryCondition(M,h,0,0,k)+compute_BoundaryCondition(M,h,0,j,k-1)-4*compute_BoundaryCondition(M,h,0,j,k)+compute_BoundaryCondition(M,h,0,j,k+1)+compute_BoundaryCondition(M,h,0,j+1,k)+compute_BoundaryCondition(M,h,i,0,k-1)-4*compute_BoundaryCondition(M,h,i,0,k)+compute_BoundaryCondition(M,h,i,0,k+1)+compute_BoundaryCondition(M,h,i+1,0,k);
//				BU[Idx+k]=BU[Idx+k]+temp*para2;
//			}
//		}
//	}
//	else if( i==1 && j==N  )
//	{
//		for(int k=1;k<=K;k++)
//		{
//			if(k==1)
//			{
//				temp=compute_BoundaryCondition(M,h,i-1,j,k)+compute_BoundaryCondition(M,h,i,j+1,k)+compute_BoundaryCondition(M,h,i,j,k-1);
//				BU[Idx+k]=BU[Idx+k]+temp*para1;
//				temp=compute_BoundaryCondition(M,h,i-1,j-1,k)+compute_BoundaryCondition(M,h,i-1,j,k-1)-4*compute_BoundaryCondition(M,h,i-1,j,k)+compute_BoundaryCondition(M,h,i-1,j,k+1)+compute_BoundaryCondition(M,h,i-1,j+1,k)+compute_BoundaryCondition(M,h,i,j-1,k-1)-4*compute_BoundaryCondition(M,h,i,j,k-1)+compute_BoundaryCondition(M,h,i,j+1,k-1)-4*compute_BoundaryCondition(M,h,i,j+1,k)+compute_BoundaryCondition(M,h,i,j+1,k+1)+compute_BoundaryCondition(M,h,i+1,j,k-1)+compute_BoundaryCondition(M,h,i+1,j+1,k);
//				BU[Idx+k]=BU[Idx+k]+temp*para2;
//			}
//			else if(k==K)
//			{
//				temp=compute_BoundaryCondition(M,h,i-1,j,k)+compute_BoundaryCondition(M,h,i,j+1,k);
//				temp_gama=compute_BoundaryCondition(M,h,i,j,k+1);
//				BU[Idx+k]=BU[Idx+k]+temp*para1;
//				BU[Idx+k]=BU[Idx+k]+temp_gama*para1;
//				temp=compute_BoundaryCondition(M,h,i-1,j-1,k)+compute_BoundaryCondition(M,h,i-1,j,k-1)-4*compute_BoundaryCondition(M,h,i-1,j,k)+compute_BoundaryCondition(M,h,i-1,j,k+1)+compute_BoundaryCondition(M,h,i-1,j+1,k)+compute_BoundaryCondition(M,h,i,j+1,k-1)-4*compute_BoundaryCondition(M,h,i,j+1,k)+compute_BoundaryCondition(M,h,i,j+1,k+1)+compute_BoundaryCondition(M,h,i+1,j+1,k);
//				temp_gama=compute_BoundaryCondition(M,h,i,j-1,k+1)-4*compute_BoundaryCondition(M,h,i,j,k+1)+compute_BoundaryCondition(M,h,i+1,j,k+1);
//				BU[Idx+k]=BU[Idx+k]+temp*para2;
//				BU[Idx+k]=BU[Idx+k]+temp_gama*para2;
//			}
//			else
//			{
//				temp=compute_BoundaryCondition(M,h,i-1,j,k)+compute_BoundaryCondition(M,h,i,j+1,k);
//				BU[Idx+k]=BU[Idx+k]+temp*para1;
//				temp=compute_BoundaryCondition(M,h,i-1,j-1,k)+compute_BoundaryCondition(M,h,i-1,j,k-1)-4*compute_BoundaryCondition(M,h,i-1,j,k)+compute_BoundaryCondition(M,h,i-1,j,k+1)+compute_BoundaryCondition(M,h,i-1,j+1,k)+compute_BoundaryCondition(M,h,i,j+1,k-1)-4*compute_BoundaryCondition(M,h,i,j+1,k)+compute_BoundaryCondition(M,h,i,j+1,k+1)+compute_BoundaryCondition(M,h,i+1,j+1,k);
//				BU[Idx+k]=BU[Idx+k]+temp*para2;
//			}
//		}
//	}
//	else if( i==1 && j>1 && j<N  )
//	{
//		for(int k=1;k<=K;k++)
//		{
//			if(k==1)
//			{
//				temp=compute_BoundaryCondition(M,h,i-1,j,k)+compute_BoundaryCondition(M,h,i,j,k-1);
//				BU[Idx+k]=BU[Idx+k]+temp*para1;
//				temp=compute_BoundaryCondition(M,h,i-1,j-1,k)+compute_BoundaryCondition(M,h,i-1,j,k-1)-4*compute_BoundaryCondition(M,h,i-1,j,k)+compute_BoundaryCondition(M,h,i-1,j,k+1)+compute_BoundaryCondition(M,h,i-1,j+1,k)+compute_BoundaryCondition(M,h,i,j-1,k-1)-4*compute_BoundaryCondition(M,h,i,j,k-1)+compute_BoundaryCondition(M,h,i,j+1,k-1)+compute_BoundaryCondition(M,h,i+1,j,k-1);
//				BU[Idx+k]=BU[Idx+k]+temp*para2;
//			}
//			else if(k==K)
//			{
//				temp=compute_BoundaryCondition(M,h,i-1,j,k);
//				temp_gama=compute_BoundaryCondition(M,h,i,j,k+1);
//				BU[Idx+k]=BU[Idx+k]+temp*para1;
//				BU[Idx+k]=BU[Idx+k]+temp_gama*para1;
//				temp=compute_BoundaryCondition(M,h,i-1,j-1,k)+compute_BoundaryCondition(M,h,i-1,j,k-1)-4*compute_BoundaryCondition(M,h,i-1,j,k)+compute_BoundaryCondition(M,h,i-1,j,k+1)+compute_BoundaryCondition(M,h,i-1,j+1,k);
//				temp_gama=compute_BoundaryCondition(M,h,i,j-1,k+1)-4*compute_BoundaryCondition(M,h,i,j,k+1)+compute_BoundaryCondition(M,h,i,j+1,k+1)+compute_BoundaryCondition(M,h,i+1,j,k+1);
//				BU[Idx+k]=BU[Idx+k]+temp*para2;
//				BU[Idx+k]=BU[Idx+k]+temp_gama*para2;
//			}
//			else
//			{
//				temp=compute_BoundaryCondition(M,h,i-1,j,k);
//				BU[Idx+k]=BU[Idx+k]+temp*para1;
//				temp=compute_BoundaryCondition(M,h,i-1,j-1,k)+compute_BoundaryCondition(M,h,i-1,j,k-1)-4*compute_BoundaryCondition(M,h,i-1,j,k)+compute_BoundaryCondition(M,h,i-1,j,k+1)+compute_BoundaryCondition(M,h,i-1,j+1,k);
//				BU[Idx+k]=BU[Idx+k]+temp*para2;
//			}
//		}
//	}
//	else if( i==M && j==1 )
//	{
//		for(int k=1;k<=K;k++)
//		{
//			if(k==1)
//			{
//				temp=compute_BoundaryCondition(M,h,i+1,j,k)+compute_BoundaryCondition(M,h,i,j-1,k)+compute_BoundaryCondition(M,h,i,j,k-1);
//				BU[Idx+k]=BU[Idx+k]+temp*para1;
//				temp=compute_BoundaryCondition(M,h,i-1,j-1,k)+compute_BoundaryCondition(M,h,i-1,j,k-1)+compute_BoundaryCondition(M,h,i,j-1,k-1)-4*compute_BoundaryCondition(M,h,i,j-1,k)+compute_BoundaryCondition(M,h,i,j-1,k+1)-4*compute_BoundaryCondition(M,h,i,j,k-1)+compute_BoundaryCondition(M,h,i,j+1,k-1)+compute_BoundaryCondition(M,h,i+1,j-1,k)+compute_BoundaryCondition(M,h,i+1,j,k-1)-4*compute_BoundaryCondition(M,h,i+1,j,k)+compute_BoundaryCondition(M,h,i+1,j,k+1)+compute_BoundaryCondition(M,h,i+1,j+1,k);
//				BU[Idx+k]=BU[Idx+k]+temp*para2;
//			}
//			else if(k==K)
//			{
//				temp=compute_BoundaryCondition(M,h,i+1,j,k)+compute_BoundaryCondition(M,h,i,j-1,k);
//				temp_gama=compute_BoundaryCondition(M,h,i,j,k+1);
//				BU[Idx+k]=BU[Idx+k]+temp*para1;
//				BU[Idx+k]=BU[Idx+k]+temp_gama*para1;
//				temp=compute_BoundaryCondition(M,h,i-1,j-1,k)+compute_BoundaryCondition(M,h,i,j-1,k-1)-4*compute_BoundaryCondition(M,h,i,j-1,k)+compute_BoundaryCondition(M,h,i,j-1,k+1)+compute_BoundaryCondition(M,h,i+1,j-1,k)+compute_BoundaryCondition(M,h,i+1,j,k-1)-4*compute_BoundaryCondition(M,h,i+1,j,k)+compute_BoundaryCondition(M,h,i+1,j,k+1)+compute_BoundaryCondition(M,h,i+1,j+1,k);
//				temp_gama=compute_BoundaryCondition(M,h,i-1,j,k+1)-4*compute_BoundaryCondition(M,h,i,j,k+1)+compute_BoundaryCondition(M,h,i,j+1,k+1);
//				BU[Idx+k]=BU[Idx+k]+temp*para2;
//				BU[Idx+k]=BU[Idx+k]+temp_gama*para2;
//			}
//			else
//			{
//				temp=compute_BoundaryCondition(M,h,i+1,j,k)+compute_BoundaryCondition(M,h,i,j-1,k);
//				BU[Idx+k]=BU[Idx+k]+temp*para1;
//				temp=compute_BoundaryCondition(M,h,i-1,j-1,k)+compute_BoundaryCondition(M,h,i,j-1,k-1)-4*compute_BoundaryCondition(M,h,i,j-1,k)+compute_BoundaryCondition(M,h,i,j-1,k+1)+compute_BoundaryCondition(M,h,i+1,j-1,k)+compute_BoundaryCondition(M,h,i+1,j,k-1)-4*compute_BoundaryCondition(M,h,i+1,j,k)+compute_BoundaryCondition(M,h,i+1,j,k+1)+compute_BoundaryCondition(M,h,i+1,j+1,k);
//				BU[Idx+k]=BU[Idx+k]+temp*para2;
//			}
//		}
//	}
//	else if( i==M && j==N )
//	{
//		for(int k=1;k<=K;k++)
//		{
//			if(k==1)
//			{
//				temp=compute_BoundaryCondition(M,h,i+1,j,k)+compute_BoundaryCondition(M,h,i,j+1,k)+compute_BoundaryCondition(M,h,i,j,k-1);
//				BU[Idx+k]=BU[Idx+k]+temp*para1;
//				temp=compute_BoundaryCondition(M,h,i-1,j,k-1)+compute_BoundaryCondition(M,h,i-1,j+1,k)+compute_BoundaryCondition(M,h,i,j-1,k-1)-4*compute_BoundaryCondition(M,h,i,j,k-1)+compute_BoundaryCondition(M,h,i,j+1,k-1)-4*compute_BoundaryCondition(M,h,i,j+1,k)+compute_BoundaryCondition(M,h,i,j+1,k+1)+compute_BoundaryCondition(M,h,i+1,j-1,k)+compute_BoundaryCondition(M,h,i+1,j,k-1)-4*compute_BoundaryCondition(M,h,i+1,j,k)+compute_BoundaryCondition(M,h,i+1,j,k+1)+compute_BoundaryCondition(M,h,i+1,j+1,k);
//				BU[Idx+k]=BU[Idx+k]+temp*para2;
//			}
//			else if(k==K)
//			{
//				temp=compute_BoundaryCondition(M,h,i+1,j,k)+compute_BoundaryCondition(M,h,i,j+1,k);
//				temp_gama=compute_BoundaryCondition(M,h,i,j,k+1);
//				BU[Idx+k]=BU[Idx+k]+temp*para1;
//				BU[Idx+k]=BU[Idx+k]+temp_gama*para1;
//				temp=compute_BoundaryCondition(M,h,i-1,j+1,k)+compute_BoundaryCondition(M,h,i,j+1,k-1)-4*compute_BoundaryCondition(M,h,i,j+1,k)+compute_BoundaryCondition(M,h,i,j+1,k+1)+compute_BoundaryCondition(M,h,i+1,j-1,k)+compute_BoundaryCondition(M,h,i+1,j,k-1)-4*compute_BoundaryCondition(M,h,i+1,j,k)+compute_BoundaryCondition(M,h,i+1,j,k+1)+compute_BoundaryCondition(M,h,i+1,j+1,k);
//				temp_gama=compute_BoundaryCondition(M,h,i-1,j,k+1)+compute_BoundaryCondition(M,h,i,j-1,k+1)-4*compute_BoundaryCondition(M,h,i,j,k+1);
//				BU[Idx+k]=BU[Idx+k]+temp*para2;
//				BU[Idx+k]=BU[Idx+k]+temp_gama*para2;
//			}
//			else
//			{
//				temp=compute_BoundaryCondition(M,h,i+1,j,k)+compute_BoundaryCondition(M,h,i,j+1,k);
//				BU[Idx+k]=BU[Idx+k]+temp*para1;
//				temp=compute_BoundaryCondition(M,h,i-1,j+1,k)+compute_BoundaryCondition(M,h,i,j+1,k-1)-4*compute_BoundaryCondition(M,h,i,j+1,k)+compute_BoundaryCondition(M,h,i,j+1,k+1)+compute_BoundaryCondition(M,h,i+1,j-1,k)+compute_BoundaryCondition(M,h,i+1,j,k-1)-4*compute_BoundaryCondition(M,h,i+1,j,k)+compute_BoundaryCondition(M,h,i+1,j,k+1)+compute_BoundaryCondition(M,h,i+1,j+1,k);
//				BU[Idx+k]=BU[Idx+k]+temp*para2;
//			}
//		}
//	}
//	else if( i==M && j>1 && j<N  )
//	{
//		for(int k=1;k<=K;k++)
//		{
//			if(k==1)
//			{
//				temp=compute_BoundaryCondition(M,h,i+1,j,k)+compute_BoundaryCondition(M,h,i,j,k-1);
//				BU[Idx+k]=BU[Idx+k]+temp*para1;
//				temp=compute_BoundaryCondition(M,h,i-1,j,k-1)+compute_BoundaryCondition(M,h,i,j-1,k-1)-4*compute_BoundaryCondition(M,h,i,j,k-1)+compute_BoundaryCondition(M,h,i,j+1,k-1)+compute_BoundaryCondition(M,h,i+1,j-1,k)+compute_BoundaryCondition(M,h,i+1,j,k-1)-4*compute_BoundaryCondition(M,h,i+1,j,k)+compute_BoundaryCondition(M,h,i+1,j,k+1)+compute_BoundaryCondition(M,h,i+1,j+1,k);
//				BU[Idx+k]=BU[Idx+k]+temp*para2;
//			}
//			else if(k==K)
//			{
//				temp=compute_BoundaryCondition(M,h,i+1,j,k);
//				temp_gama=compute_BoundaryCondition(M,h,i,j,k+1);
//				BU[Idx+k]=BU[Idx+k]+temp*para1;
//				BU[Idx+k]=BU[Idx+k]+temp_gama*para1;
//				temp=compute_BoundaryCondition(M,h,i+1,j-1,k)+compute_BoundaryCondition(M,h,i+1,j,k-1)-4*compute_BoundaryCondition(M,h,i+1,j,k)+compute_BoundaryCondition(M,h,i+1,j,k+1)+compute_BoundaryCondition(M,h,i+1,j+1,k);
//				temp_gama=compute_BoundaryCondition(M,h,i-1,j,k+1)+compute_BoundaryCondition(M,h,i,j-1,k+1)-4*compute_BoundaryCondition(M,h,i,j,k+1)+compute_BoundaryCondition(M,h,i,j+1,k+1);
//				BU[Idx+k]=BU[Idx+k]+temp*para2;
//				BU[Idx+k]=BU[Idx+k]+temp_gama*para2;
//			}
//			else
//			{
//				temp=compute_BoundaryCondition(M,h,i+1,j,k);
//				BU[Idx+k]=BU[Idx+k]+temp*para1;
//				temp=compute_BoundaryCondition(M,h,i+1,j-1,k)+compute_BoundaryCondition(M,h,i+1,j,k-1)-4*compute_BoundaryCondition(M,h,i+1,j,k)+compute_BoundaryCondition(M,h,i+1,j,k+1)+compute_BoundaryCondition(M,h,i+1,j+1,k);
//				BU[Idx+k]=BU[Idx+k]+temp*para2;
//			}
//		}
//	}
//	else if( i>1 && i<M && j==1  )
//	{
//		for(int k=1;k<=K;k++)
//		{	
//			if(k==1)
//			{
//				temp=compute_BoundaryCondition(M,h,i,j-1,k)+compute_BoundaryCondition(M,h,i,j,k-1);
//				BU[Idx+k]=BU[Idx+k]+temp*para1;
//				temp=compute_BoundaryCondition(M,h,i-1,j-1,k)+compute_BoundaryCondition(M,h,i-1,j,k-1)+compute_BoundaryCondition(M,h,i,j-1,k-1)-4*compute_BoundaryCondition(M,h,i,j-1,k)+compute_BoundaryCondition(M,h,i,j-1,k+1)-4*compute_BoundaryCondition(M,h,i,j,k-1)+compute_BoundaryCondition(M,h,i,j+1,k-1)+compute_BoundaryCondition(M,h,i+1,j-1,k)+compute_BoundaryCondition(M,h,i+1,j,k-1);
//				BU[Idx+k]=BU[Idx+k]+temp*para2;
//			}
//			else if(k==K)
//			{
//				temp=compute_BoundaryCondition(M,h,i,j-1,k);
//				temp_gama=compute_BoundaryCondition(M,h,i,j,k+1);
//				BU[Idx+k]=BU[Idx+k]+temp*para1;
//				BU[Idx+k]=BU[Idx+k]+temp_gama*para1;
//				temp=compute_BoundaryCondition(M,h,i-1,j-1,k)+compute_BoundaryCondition(M,h,i,j-1,k-1)-4*compute_BoundaryCondition(M,h,i,j-1,k)+compute_BoundaryCondition(M,h,i,j-1,k+1)+compute_BoundaryCondition(M,h,i+1,j-1,k);
//				temp_gama=compute_BoundaryCondition(M,h,i-1,j,k+1)-4*compute_BoundaryCondition(M,h,i,j,k+1)+compute_BoundaryCondition(M,h,i,j+1,k+1)+compute_BoundaryCondition(M,h,i+1,j,k+1);
//				BU[Idx+k]=BU[Idx+k]+temp*para2;
//				BU[Idx+k]=BU[Idx+k]+temp_gama*para2;
//			}
//			else
//			{
//				temp=compute_BoundaryCondition(M,h,i,j-1,k);
//				BU[Idx+k]=BU[Idx+k]+temp*para1;
//				temp=compute_BoundaryCondition(M,h,i-1,j-1,k)+compute_BoundaryCondition(M,h,i,j-1,k-1)-4*compute_BoundaryCondition(M,h,i,j-1,k)+compute_BoundaryCondition(M,h,i,j-1,k+1)+compute_BoundaryCondition(M,h,i+1,j-1,k);
//				BU[Idx+k]=BU[Idx+k]+temp*para2;
//			}
//		}
//	}
//	else if( i>1 && i<M && j==N   )
//	{
//		for(int k=1;k<=K;k++)
//		{
//			if(k==1)
//			{
//				temp=compute_BoundaryCondition(M,h,i,j+1,k)+compute_BoundaryCondition(M,h,i,j,k-1);
//				BU[Idx+k]=BU[Idx+k]+temp*para1;
//				temp=compute_BoundaryCondition(M,h,i-1,j,k-1)+compute_BoundaryCondition(M,h,i-1,j+1,k)+compute_BoundaryCondition(M,h,i,j-1,k-1)-4*compute_BoundaryCondition(M,h,i,j,k-1)+compute_BoundaryCondition(M,h,i,j+1,k-1)-4*compute_BoundaryCondition(M,h,i,j+1,k)+compute_BoundaryCondition(M,h,i,j+1,k+1)+compute_BoundaryCondition(M,h,i+1,j,k-1)+compute_BoundaryCondition(M,h,i+1,j+1,k);
//				BU[Idx+k]=BU[Idx+k]+temp*para2;
//			}
//			else if(k==K)
//			{
//				temp=compute_BoundaryCondition(M,h,i,j+1,k);
//				temp_gama=compute_BoundaryCondition(M,h,i,j,k+1);
//				BU[Idx+k]=BU[Idx+k]+temp*para1;
//				BU[Idx+k]=BU[Idx+k]+temp_gama*para1;
//				temp=compute_BoundaryCondition(M,h,i-1,j+1,k)+compute_BoundaryCondition(M,h,i,j+1,k-1)-4*compute_BoundaryCondition(M,h,i,j+1,k)+compute_BoundaryCondition(M,h,i,j+1,k+1)+compute_BoundaryCondition(M,h,i+1,j+1,k);
//				temp_gama=compute_BoundaryCondition(M,h,i-1,j,k+1)+compute_BoundaryCondition(M,h,i,j-1,k+1)-4*compute_BoundaryCondition(M,h,i,j,k+1)+compute_BoundaryCondition(M,h,i+1,j,k+1);
//				BU[Idx+k]=BU[Idx+k]+temp*para2;
//				BU[Idx+k]=BU[Idx+k]+temp_gama*para2;
//			}
//			else
//			{
//				temp=compute_BoundaryCondition(M,h,i,j+1,k);
//				BU[Idx+k]=BU[Idx+k]+temp*para1;
//				temp=compute_BoundaryCondition(M,h,i-1,j+1,k)+compute_BoundaryCondition(M,h,i,j+1,k-1)-4*compute_BoundaryCondition(M,h,i,j+1,k)+compute_BoundaryCondition(M,h,i,j+1,k+1)+compute_BoundaryCondition(M,h,i+1,j+1,k);
//				BU[Idx+k]=BU[Idx+k]+temp*para2;
//			}
//		}
//	}
//	else if( i>1 && i<M && j>1 && j<N )
//	{
//		for(int k=1;k<=K;k++)
//		{
//			if(k==1)
//			{
//				temp=compute_BoundaryCondition(M,h,i,j,k-1);
//				BU[Idx+k]=BU[Idx+k]+temp*para1;
//				temp=compute_BoundaryCondition(M,h,i-1,j,k-1)+compute_BoundaryCondition(M,h,i,j-1,k-1)-4*compute_BoundaryCondition(M,h,i,j,k-1)+compute_BoundaryCondition(M,h,i,j+1,k-1)+compute_BoundaryCondition(M,h,i+1,j,k-1);
//				BU[Idx+k]=BU[Idx+k]+temp*para2;
//			}
//			else if(k==K)
//			{
//				temp_gama=compute_BoundaryCondition(M,h,i,j,k+1);
//				BU[Idx+k]=BU[Idx+k]+temp_gama*para1;
//				temp_gama=compute_BoundaryCondition(M,h,i-1,j,k+1)+compute_BoundaryCondition(M,h,i,j-1,k+1)-4*compute_BoundaryCondition(M,h,i,j,k+1)+compute_BoundaryCondition(M,h,i,j+1,k+1)+compute_BoundaryCondition(M,h,i+1,j,k+1);
//				BU[Idx+k]=BU[Idx+k]+temp_gama*para2;
//			}
//			else
//			{
//				// do nothing
//			}
//		}
//	}
//}
