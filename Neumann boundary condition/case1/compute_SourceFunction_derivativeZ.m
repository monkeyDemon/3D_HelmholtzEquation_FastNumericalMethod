function [ fz ] = compute_SourceFunction_derivativeZ( para )
% ����Դ������z����ƫ����fz
% ��Ӧ��(:,:,L+1)��Щ�����

% ��������Դ����Ϊ
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

