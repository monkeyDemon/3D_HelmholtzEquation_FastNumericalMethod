function [ SFlist ] = compute_SourceFunction(x,z,para)
%¼ÆËãÔ´º¯Êı F x,:,z

N = para(1);
h = para(2);
Ystart = para(5);
Yend = para(6);

vector_y = [Ystart+h : h : Yend-h]';

SFlist =zeros(N,1);
%SFlist = sin(vector_y);


end

