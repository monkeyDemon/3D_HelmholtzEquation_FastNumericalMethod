function [ u ] = compute_realU( x,y,z,para )
%计算U的真实值
% need to modify when equation changed

k = para(9);

u = sin(pi*x)*sin(pi*y)*sin(k*z);


end

