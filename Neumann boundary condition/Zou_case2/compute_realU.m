function [ u ] = compute_realU( x,y,z )
%计算U的真实值
% need to modify when equation changed
u = (sin(pi*x)*sin(pi*y)/sinh(sqrt(2)*pi))*(2*sinh(sqrt(2)*pi*z)+sinh(sqrt(2)*pi*(1-z)));


end

