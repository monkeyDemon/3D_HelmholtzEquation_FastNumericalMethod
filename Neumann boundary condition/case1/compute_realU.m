function [ u ] = compute_realU( x,y,z,para )
%����U����ʵֵ
% need to modify when equation changed

k = para(9);

u = sin(pi*x)*sin(pi*y)*sin(k*z);


end

