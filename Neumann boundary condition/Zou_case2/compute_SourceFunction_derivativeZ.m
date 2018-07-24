function [ fz ] = compute_SourceFunction_derivativeZ( para )
% 计算源函数沿z方向偏导数fz
% 对应于(:,:,L+1)这些网格点

% 本算例的源函数为0
% 因此无需计算fz

M = para(1);
% h = para(2);
% Xstart = para(3);
% Xend = para(4);
% Ystart = para(5);
% Yend = para(6);

fz = zeros(M*M,1);

end

