function [ g ] = compute_NeumannCondition( para )
% 计算对应（:,:,L+1）处的纽曼边界条件g

% 本算例的精确解为：
% u(x,y,z) = sin(pix)sin(piy)sin(kz)
% 因此，纽曼边界条件g = 偏u/偏z(z=1) 为：
% g(x,y,z) = ksin(pix)sin(piy)cos(k)

M = para(1);
h = para(2);
Xstart = para(3);
Xend = para(4);
Ystart = para(5);
Yend = para(6);
k = para(9);

vector_x = [Xstart+h : h : Xend-h]';
vector_y = [Ystart+h : h : Yend-h]';

g = zeros(M*M,1);
YZ_term = k*sin(pi*vector_y)*cos(k);
for i = 1:M
    g((i-1)*M+1:i*M) = sin(pi*vector_x(i))*YZ_term;
end

end

