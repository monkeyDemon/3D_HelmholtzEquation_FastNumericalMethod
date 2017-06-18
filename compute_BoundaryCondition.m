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
z=Zstart+h*j;


if( i>=1 && i<=M && j>=1 && j<=N && k>=1 && k<=K)
	a=10086
end



if(i==0)
	value=sin(pi*y)*sin(pi*z);
elseif(i==k)
	value=2*sin(pi*y)*sin(pi*z);
else
	value=0;
end


