function [ fz ] = compute_SourceFunction_derivativeZ( para )
% 计算源函数沿z方向偏导数fz
% 对应于(:,:,L+1)这些网格点

% 本算例的源函数为
% f(i,j,k)= (-2pi^2)sin(pix)sin(piy)sin(kz)

M = para(1);
h = para(2);
Xstart = para(3);
Xend = para(4);
Ystart = para(5);
Yend = para(6);
k = para(9);

vector_x = [Xstart+h : h : Xend-h]';
vector_y = [Ystart+h : h : Yend-h]';

fz = zeros(M*M,1);
YZ_term = k*(-2*pi*pi)*sin(pi*vector_y)*cos(k);
for i = 1:M
    fz((i-1)*M+1:i*M) = sin(pi*vector_x(i))*YZ_term;
end

end

