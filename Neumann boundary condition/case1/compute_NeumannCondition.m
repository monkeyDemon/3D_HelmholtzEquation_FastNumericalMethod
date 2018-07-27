function [ g ] = compute_NeumannCondition( para )
% �����Ӧ��:,:,L+1������Ŧ���߽�����g

% �������ľ�ȷ��Ϊ��
% u(x,y,z) = sin(pix)sin(piy)sin(kz)
% ��ˣ�Ŧ���߽�����g = ƫu/ƫz(z=1) Ϊ��
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

