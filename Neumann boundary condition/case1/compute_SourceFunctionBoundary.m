function [ value,state ] = compute_SourceFunctionBoundary( logo, para )
% 计算源函数在边界处的取值
% need to modify when equation changed

% SourceFunction = (3*pi*pi+K0)*sin(pi*x)*sin(pi*z)*sin(pi*vector_y);
% 故边界项全为0

% M = para(1);
% h = para(2);
% Xstart = para(3);
% Xend = para(4);
% Ystart = para(5);
% Yend = para(6);
% Zstart = para(7);
% Zend = para(8);
% K0 = para(9);
% 
% vector_x = [Xstart+h : h : Xend-h]';
% vector_y = [Ystart+h : h : Yend-h]';
% vector_z = [Zstart+h : h : Zend-h]';


value = 0;
if(strcmp(logo,'BF_top'))
    state = 0;
elseif(strcmp(logo,'BF_bottom'))
    state = 0;
elseif(strcmp(logo,'BF_left'))
    state = 0;
elseif(strcmp(logo,'BF_right'))
    state = 0;
elseif(strcmp(logo,'BF_front'))
    state = 0;
elseif(strcmp(logo,'BF_back'))
    state = 0;
else
   throw('error!'); 
end

