function [ SFlist ] = compute_SourceFunction(x,z,para)
%计算源函数 F x,:,z

% 源函数f(i,j,k)= (-2pi^2)sin(pix)sin(piy)sin(kz)
N = para(1);
h = para(2);
Ystart = para(5);
Yend = para(6);
k = para(9);

vector_y = [Ystart+h : h : Yend-h]';

SFlist = (-2*pi*pi)*sin(pi*x)*sin(k*z)*sin(pi*vector_y);


end

