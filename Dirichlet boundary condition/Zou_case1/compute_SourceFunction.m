function [ SFlist ] = compute_SourceFunction(x,z,para)
%¼ÆËãÔ´º¯Êý F x,:,z

N = para(1);
h = para(2);
Ystart = para(5);
Yend = para(6);
K0 = para(9);

vector_y = [Ystart+h : h : Yend-h]';

SFlist = -1*(3*pi*pi-K0*K0)*sin(pi*x)*sin(pi*z)*sin(pi*vector_y);


end

