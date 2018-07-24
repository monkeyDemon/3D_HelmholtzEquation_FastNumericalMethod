function [  ] = plot3Dsolution( para,U)
% ����3D�����ڵ���ֵ��

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


%% ��ʼ��ͼ
figure
% ѡ��ɫ��
colormap(jet(128));
% colormap(hsv);  
% colormap(hsv); 
% colormap(jet); 
% colormap(hot); 
% colormap(cool); 
% colormap(spring); 
% colormap(summer); 
% colormap(autumn); 
% colormap(winter); 
halfSize = floor((M+1)/2);


%% �����������
X = zeros(M,M);
for i =1:M
   X(:,i) = vector_x(i); 
end
Y = zeros(M,M);
for i =1:M
   Y(i,:) = vector_y(i); 
end
Z = ones(M,M)*(Zstart+h);
R = zeros(M,M);
for i =1:M
    for j=1:M
        R(i,j) = U((i-1)*M*M+(j-1)*M+1);
    end
end
surf(X,Y,Z,R);
shading interp;  
colorbar;
hold on 

%% �����������
X = ones(M,halfSize)*(Xstart+h);
Y = zeros(M,halfSize);
for i =1:halfSize
   Y(:,i) = vector_y; 
end
Z = zeros(M,halfSize);
for i =1:halfSize
   Z(:,i) =  vector_z(i);
end
R = zeros(M,halfSize);
for j =1:M
    for k=1:halfSize
        R(j,k) = U((j-1)*M+k);
    end
end
surf(X,Y,Z,R);
shading interp;  
hold on 

X = ones(M-halfSize+1,M-halfSize+1)*(Xstart+h);
Y = zeros(M-halfSize+1,M-halfSize+1);
for i =1:M-halfSize+1
   Y(:,i) = vector_y(halfSize:M); 
end
Z = zeros(M-halfSize+1,M-halfSize+1);
for i =1:M-halfSize+1
   Z(:,i) =  vector_z(halfSize+i-1);
end
R = zeros(M-halfSize+1,M-halfSize+1);
for j =1:M-halfSize+1
    for k=1:M-halfSize+1
        R(j,k) = U((halfSize+j-2)*M+halfSize+k-1);
    end
end
surf(X,Y,Z,R);
shading interp;  
hold on 


%% �����������
X = ones(M,halfSize)*(Xend-h);
Y = zeros(M,halfSize);
for i =1:halfSize
   Y(:,i) = vector_y; 
end
Z = zeros(M,halfSize);
for i =1:halfSize
   Z(:,i) =  vector_z(i);
end
R = zeros(M,halfSize);
for j =1:M
    for k=1:halfSize
        R(j,k) = U((M-1)*M*M+(j-1)*M+k);
    end
end
surf(X,Y,Z,R);
shading interp;  
hold on 

X = ones(M-halfSize+1,M-halfSize+1)*(Xend-h);
Y = zeros(M-halfSize+1,M-halfSize+1);
for i =1:M-halfSize+1
   Y(:,i) = vector_y(halfSize:M); 
end
Z = zeros(M-halfSize+1,M-halfSize+1);
for i =1:M-halfSize+1
   Z(:,i) =  vector_z(halfSize+i-1);
end
R = zeros(M-halfSize+1,M-halfSize+1);
for j =1:M-halfSize+1
    for k=1:M-halfSize+1
        R(j,k) = U((M-1)*M*M+(halfSize+j-2)*M+halfSize+k-1);
    end
end
surf(X,Y,Z,R);
shading interp;  
hold on 


%% ����y=0����
X = zeros(halfSize,M);
for i =1:M
   X(:,i) = vector_x(i); 
end
Y = ones(halfSize,M)*(Ystart+h);
Z = zeros(halfSize,M);
for i=1:M
   Z(:,i) = vector_z(1:halfSize); 
end
R = zeros(halfSize,M);
for i =1:M
    for k=1:halfSize
        R(k,i) = U((i-1)*M*M+k);
    end
end
surf(X,Y,Z,R);
shading interp;  
hold on 

%% ����y=0.5�İ�ƽ��
X = zeros(M,M);
for i =1:M
   X(:,i) = vector_x(i); 
end
Y = ones(M,M)*(Ystart+halfSize*h);
Z = zeros(M,M);
for i=1:M
   Z(:,i) = vector_z; 
end
R = zeros(M,M);
for i =1:M
    for k=1:M
        R(k,i) = U((i-1)*M*M+(halfSize-1)*M+k);
    end
end
surf(X,Y,Z,R);
shading interp;  
hold on 

%% ����z=0.5�İ�ƽ��
X = zeros(halfSize,M);
for i =1:M
   X(:,i) = vector_x(i); 
end
Y = zeros(halfSize,M);
for i=1:halfSize
   Y(i,:) = vector_y(i); 
end
Z = ones(halfSize,M)*(Zstart+halfSize*h);
R = zeros(halfSize,M);
for i =1:M
    for j=1:halfSize
        R(j,i) = U((i-1)*M*M+(j-1)*M+halfSize);
    end
end
surf(X,Y,Z,R);
shading interp;  
hold on 


%% ����z=1�İ�ƽ��
X = zeros(halfSize,M);
for i =1:M
   X(:,i) = vector_x(i); 
end
Y = zeros(halfSize,M);
for i=1:halfSize
   Y(i,:) = vector_y(halfSize+i-1); 
end
Z = ones(halfSize,M)*(Zend-h);
R = zeros(halfSize,M);
for i =1:M
    for j=1:halfSize
        R(j,i) = U((i-1)*M*M+(halfSize+j-2)*M+M);
    end
end
surf(X,Y,Z,R);
shading interp;  
hold off


end

