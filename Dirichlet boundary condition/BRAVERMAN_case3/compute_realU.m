function [ u ] = compute_realU( x,y,z,K0 )
%����U����ʵֵ
% need to modify when equation changed
u = cos(2*x)*cos(sqrt(1+K0*K0)*y)*sinh(sqrt(5)*z)/sinh(sqrt(5)*pi);


end

