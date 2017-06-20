function [ value ] = compute_BoundaryCondition( i,j,k  )

M=15;
N=15;
K=15;

Xstart=0;
Xend=1;
Ystart=0;
Yend=1;
Zstart=0;
Zend=1;

h=(Xend-Xstart)/(1+M);

x=Xstart+h*i;
y=Ystart+h*j;
z=Zstart+h*k;


if(i==0)
	value=sin(pi*y)*sin(pi*z);
elseif(i==M+1)
	value=2*sin(pi*y)*sin(pi*z);
else
	value=0;
end


