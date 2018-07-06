function [ g ] = compute_NeumannCondition( para )
% 计算对应（:,:,L+1）处的纽曼边界条件g

% 本算例的精确解为：
% u(x,y,z) = sin(pix)sin(piy)[2sinh(sqrt(2)*piz)+sinh(sqrt(2)*pi(1-z))]/sinh(sqrt(2)pi)
% 因此，纽曼边界条件g = 偏u/偏z 为：
% u(x,y,z) = sin(pix)sin(piy)[2sqrt(2)pi*cosh(sqrt(2)*piz) -sqrt(2)picosh(sqrt(2)*pi(1-z))]/sinh(sqrt(2)pi)

M = para(1);
h = para(2);
Xstart = para(3);
Xend = para(4);
Ystart = para(5);
Yend = para(6);

vector_x = [Xstart+h : h : Xend-h]';
vector_y = [Ystart+h : h : Yend-h]';

g = zeros(M*M,1);
Z_term = 2*sqrt(2)*pi*cosh(sqrt(2)*pi*1) -sqrt(2)*pi*cosh(sqrt(2)*pi*0);
const_term = Z_term / sinh(sqrt(2)*pi);
YZ_term = sin(pi*vector_y)*const_term;
for i = 1:M
    g((i-1)*M+1:i*M) = sin(pi*vector_x(i))*YZ_term;
end

end

