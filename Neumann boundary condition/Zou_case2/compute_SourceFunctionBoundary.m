function [ value,state ] = compute_SourceFunctionBoundary( logo, para )
% 计算源函数在边界处的取值
% need to modify when equation changed


M = para(1);
h = para(2);
Xstart = para(3);
Xend = para(4);
Ystart = para(5);
Yend = para(6);
Zstart = para(7);
Zend = para(8);


vector_x = [Xstart+h : h : Xend-h]';
vector_y = [Ystart+h : h : Yend-h]';
vector_z = [Zstart+h : h : Zend-h]';


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
elseif(strcmp(logo,'BF_top_K+2'))
%     value = zeros(M*M,1);
%     tempY=sin(pi*vector_y);
%     for xindex=1:M
%        value((xindex-1)*M+1:xindex*M) = sin(pi*vector_x(xindex))*tempY; 
%     end
%     state = 1;
    state = 0;
elseif(strcmp(logo,'BF_top_K'))
%     value = zeros(M*M,1);
%     tempY=sin(pi*vector_y);
%     for xindex=1:M
%        value((xindex-1)*M+1:xindex*M) = sin(pi*vector_x(xindex))*tempY; 
%     end
%     state = 1;
    state = 0;
else
   throw('error!'); 
end

