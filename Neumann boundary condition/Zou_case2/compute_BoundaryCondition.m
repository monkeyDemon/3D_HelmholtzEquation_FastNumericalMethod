function [ value,state ] = compute_BoundaryCondition( logo, para )
% 计算dirichelet边界条件
% need to modify when equation changed

% global M;
% global h;
% global Xstart Xend;
% global Ystart Yend;
% global Zstart Zend;

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
if(strcmp(logo,'BP_bottom'))
    value = zeros(M*M,1);
    tempY=sin(pi*vector_y);
    for xindex=1:M
       value((xindex-1)*M+1:xindex*M) = sin(pi*vector_x(xindex))*tempY; 
    end
    state = 1;
elseif(strcmp(logo,'BP_left'))
    state = 0;
elseif(strcmp(logo,'BP_right'))
    state = 0;
elseif(strcmp(logo,'BP_front'))
    state = 0;
elseif(strcmp(logo,'BP_back'))
    state = 0;
elseif(strcmp(logo,'BE_t1'))
    state = 0;
elseif(strcmp(logo,'BE_t2'))
    state = 0;
elseif(strcmp(logo,'BE_t3'))
    state = 0;
elseif(strcmp(logo,'BE_t4'))
    state = 0;
elseif(strcmp(logo,'BE_b1'))
    state = 0;
elseif(strcmp(logo,'BE_b2'))
    state = 0;
elseif(strcmp(logo,'BE_b3'))
    state = 0;
elseif(strcmp(logo,'BE_b4'))
    state = 0;
elseif(strcmp(logo,'BE_l1'))
    state = 0;
elseif(strcmp(logo,'BE_l2'))
    state = 0;
elseif(strcmp(logo,'BE_l3'))
    state = 0;
elseif(strcmp(logo,'BE_l4'))
    state = 0;
else
    state = 0;
end