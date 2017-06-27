function [ value ] = compute_BoundaryCondition( h,i,j,k  )
% need to modify when equation changed


M=63;

Xstart=0;
Ystart=0;
Zstart=0;

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


