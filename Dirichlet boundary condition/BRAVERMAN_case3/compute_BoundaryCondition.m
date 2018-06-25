function [ value,state ] = compute_BoundaryCondition( logo, para )
% 计算dirichelet边界条件
% need to modify when equation changed

M = para(1);
h = para(2);
Xstart = para(3);
Xend = para(4);
Ystart = para(5);
Yend = para(6);
Zstart = para(7);
Zend = para(8);
K0 = para(9);


vector_x = [Xstart+h : h : Xend-h]';
vector_y = [Ystart+h : h : Yend-h]';
vector_z = [Zstart+h : h : Zend-h]';


value = 0;
if(strcmp(logo,'BP_top'))
    value = zeros(M*M,1);
    const = sinh(sqrt(5))/sinh(sqrt(5)*pi);
    tempY=cos(sqrt(K0*K0+1)*vector_y);
    for xindex=1:M
       value((xindex-1)*M+1:xindex*M) = cos(2*vector_x(xindex))*tempY*const; 
    end
    state = 1;
elseif(strcmp(logo,'BP_bottom'))
    state = 0;
elseif(strcmp(logo,'BP_left'))
    value = zeros(M*M,1);
    const = 1.0/sinh(sqrt(5)*pi);
    tempZ= sinh(sqrt(5)*vector_z);
    for yindex=1:M
       value((yindex-1)*M+1:yindex*M) = cos(sqrt(K0*K0+1)*vector_y(yindex)) * tempZ * const; 
    end
    state = 1;
elseif(strcmp(logo,'BP_right'))
    value = zeros(M*M,1);
    const = cos(2)/sinh(sqrt(5)*pi);
    tempZ= sinh(sqrt(5)*vector_z);
    for yindex=1:M
       value((yindex-1)*M+1:yindex*M) = cos(sqrt(K0*K0+1)*vector_y(yindex)) * tempZ * const; 
    end
    state = 1;
elseif(strcmp(logo,'BP_front'))
    value = zeros(M*M,1);
    const = 1.0/sinh(sqrt(5)*pi);
    tempZ= sinh(sqrt(5)*vector_z);
    for xindex=1:M
       value((xindex-1)*M+1:xindex*M) = cos(2*vector_x(xindex)) * tempZ * const; 
    end
    state = 1;
elseif(strcmp(logo,'BP_back'))
    value = zeros(M*M,1);
    const = cos(sqrt(K0*K0+1)) / sinh(sqrt(5)*pi);
    tempZ= sinh(sqrt(5)*vector_z);
    for xindex=1:M
       value((xindex-1)*M+1:xindex*M) = cos(2*vector_x(xindex)) * tempZ * const; 
    end
    state = 1;
elseif(strcmp(logo,'BE_t1'))
    tempY = cos(sqrt(K0*K0+1)*0);
    tempZ = sinh(sqrt(5)*1) / sinh(sqrt(5)*pi);
    value = cos(2*vector_x)*tempY*tempZ;
    state = 1;
elseif(strcmp(logo,'BE_t2'))
    tempX = cos(2*1);
    tempZ = sinh(sqrt(5)*1) / sinh(sqrt(5)*pi);
    value = cos(sqrt(K0*K0+1)*vector_y)*tempX*tempZ;
    state = 1;
elseif(strcmp(logo,'BE_t3'))
    tempY = cos(sqrt(K0*K0+1)*1);
    tempZ = sinh(sqrt(5)*1) / sinh(sqrt(5)*pi);
    value = cos(2*vector_x)*tempY*tempZ;
    state = 1;
elseif(strcmp(logo,'BE_t4'))
    tempX = cos(2*0);
    tempZ = sinh(sqrt(5)*1) / sinh(sqrt(5)*pi);
    value = cos(sqrt(K0*K0+1)*vector_y)*tempX*tempZ;
    state = 1;
elseif(strcmp(logo,'BE_b1'))
    state = 0;
elseif(strcmp(logo,'BE_b2'))
    state = 0;
elseif(strcmp(logo,'BE_b3'))
    state = 0;
elseif(strcmp(logo,'BE_b4'))
    state = 0;
elseif(strcmp(logo,'BE_l1'))
    tempX = cos(2*0);
    tempY = cos(sqrt(K0*K0+1)*0);
    value = tempX * tempY * sinh(sqrt(5)*vector_z) / sinh(sqrt(5)*pi);
    state = 1;
elseif(strcmp(logo,'BE_l2'))
    tempX = cos(2*1);
    tempY = cos(sqrt(K0*K0+1)*0);
    value = tempX * tempY * sinh(sqrt(5)*vector_z) / sinh(sqrt(5)*pi);
    state = 1;
elseif(strcmp(logo,'BE_l3'))
    tempX = cos(2*1);
    tempY = cos(sqrt(K0*K0+1)*1);
    value = tempX * tempY * sinh(sqrt(5)*vector_z) / sinh(sqrt(5)*pi);
    state = 1;
elseif(strcmp(logo,'BE_l4'))
    tempX = cos(2*0);
    tempY = cos(sqrt(K0*K0+1)*1);
    value = tempX * tempY * sinh(sqrt(5)*vector_z) / sinh(sqrt(5)*pi);
    state = 1;
else
    state = 0;
end