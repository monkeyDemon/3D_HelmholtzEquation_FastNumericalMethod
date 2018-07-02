% author: AnSheng
% this code can sovle the 3 dimentional Helmholtz function with dirichlet boundary condition
% use a 4 order fast numerical method

% this demo is the save memory version
% only solve one plane in three-dimensional space
% the spatial complexity is O(M^2)

% if you want to run our code on your numerical examples
% these functions need to be modified:
% threeD4Order_dirichlet_finegrid.m 
% compute_BoundaryCondition.m            the function to compute boundary value of u
% compute_SourceFunctionBoundary.m       the function to compute boundary value of source function f 
% compute_SourceFunction.m               the function to compute f(accurately saying is calculate f(x,:,z))
% compute_realU.m                        the function to compute real solution of u
% and you should change these variables according to the real situation:
% M
% Xstart,Ystart,Zstart 
% Xend,  Yend,  Zend
% hasSourceFunction
% extractZ


% -----------------------------------------------------------------------------------------
% here the save memory version code is below:
% -----------------------------------------------------------------------------------------
clc,clear
close all
format long;
warning off; %忽略解方程时的精度警告



%-----------------need to change on different numerical examples-------------------
K0=0;

M=15;

Xstart=0;
Xend=1;
Ystart=0;
Yend=1;
Zstart=0;
Zend=1;

hasSourceFunction = 0;  % 标识是否存在源函数，即f是否为0
extractZ=(M+1)/2;% set which layer we extract ,then we can compare with the true value of u

%-----------------need to change on different numerical examples-------------------


% in order to simplified the derivation process, in our code we also have hx=hy=hz
hx=(Xend-Xstart)/(M+1);
hy=(Yend-Ystart)/(M+1);
hz=(Zend-Zstart)/(M+1);
h=hx; 


    
tic

Lambda=2*(M+1)*sin((1:M)*pi/(2*(M+1)))/(Xend-Xstart);
Lambda=-Lambda.*Lambda;

Miu=2*(M+1)*sin((1:M)*pi/(2*(M+1)))/(Yend-Ystart);
Miu=-Miu.*Miu;

P=2*(M+1)*sin((1:M)*pi/(2*(M+1)))/(Zend-Zstart);
P=-P.*P;

% compute FFT matrix
SM=FstMatrix(M);

I = [1 : M, 2 : M, 1 : M - 1]' ;
J = [1 : M, 1 : M - 1, 2 : M]' ;
S = [-2 * ones(M, 1) ; ones(2 * M - 2, 1)] ;
AK = sparse(I, J, S, M, M) / h / h ;

aN1 = sparse (1, 1, 1, M, 1)/(h*h);
aN2 = sparse (M, 1, 1, M, 1)/(h*h);

IK = sparse(1 : M, 1 : M, 1);





%--------------------------start the timer-----------------------------------

% 利用dirichlet边界条件计算出各种边界面/边处的值
para = [M, h, Xstart, Xend, Ystart, Yend, Zstart, Zend, K0];
[U_top,U_top_state] = compute_BoundaryCondition( 'BP_top',para );             % compute U::K+1
[U_bottom,U_bottom_state] = compute_BoundaryCondition( 'BP_bottom',para );    % compute U::0
[U_left,U_left_state] = compute_BoundaryCondition( 'BP_left',para );          % compute U0::
[U_right,U_right_state] = compute_BoundaryCondition( 'BP_right',para );       % compute UM+1::
[U_front,U_front_state] = compute_BoundaryCondition( 'BP_front',para );       % compute U:0:
[U_back,U_back_state] = compute_BoundaryCondition( 'BP_back',para );          % compute U:N+1:

[U_t1,U_t1_state] = compute_BoundaryCondition( 'BE_t1',para );                % compute U : 0 K+1
[U_t2,U_t2_state] = compute_BoundaryCondition( 'BE_t2',para );                % compute U M+1 : K+1
[U_t3,U_t3_state] = compute_BoundaryCondition( 'BE_t3',para );                % compute U : N+1 K+1
[U_t4,U_t4_state] = compute_BoundaryCondition( 'BE_t4',para );                % compute U 0 : K+1
[U_b1,U_b1_state] = compute_BoundaryCondition( 'BE_b1',para );                % compute U : 0 0
[U_b2,U_b2_state] = compute_BoundaryCondition( 'BE_b2',para );                % compute U M+1 : 0
[U_b3,U_b3_state] = compute_BoundaryCondition( 'BE_b3',para );                % compute U : N+1 0
[U_b4,U_b4_state] = compute_BoundaryCondition( 'BE_b4',para );                % compute 0 : 0
[U_l1,U_l1_state] = compute_BoundaryCondition( 'BE_l1',para );                % compute 0 0 :
[U_l2,U_l2_state] = compute_BoundaryCondition( 'BE_l2',para );                % compute M+1 0 :
[U_l3,U_l3_state] = compute_BoundaryCondition( 'BE_l3',para );                % compute M+1 N+1 :
[U_l4,U_l4_state] = compute_BoundaryCondition( 'BE_l4',para );                % compute 0 N+1 :

[F_top,F_top_state] = compute_SourceFunctionBoundary( 'BF_top',para );             % compute F::K+1
[F_bottom,F_bottom_state] = compute_SourceFunctionBoundary( 'BF_bottom',para );    % compute F::0
[F_left,F_left_state] = compute_SourceFunctionBoundary( 'BF_left',para );          % compute F0::
[F_right,F_right_state] = compute_SourceFunctionBoundary( 'BF_right',para );       % compute FM+1::
[F_front,F_front_state] = compute_SourceFunctionBoundary( 'BF_front',para );       % compute F:0:
[F_back,F_back_state] = compute_SourceFunctionBoundary( 'BF_back',para );          % compute F:N+1:



% U_top_bar = （SM×SN）U::K+1 
U_top_bar = zeros(M*M,1);
if(U_top_state==1) %另一种写法，更快
    for j = 1:M
        SNU = SM * U_top((j-1)*M+1:(j-1)*M+M);
        for i = 1:M
            U_top_bar((i-1)*M+1:(i-1)*M+M) = U_top_bar((i-1)*M+1:(i-1)*M+M) + SM(i,j) * SNU;
        end
    end
end
% if(U_top_state==1)
%     for i=1:M
%         temp=zeros(N,1);
%         for j=1:M
%             temp= temp + SN*U_top((j-1)*M+1:(j-1)*M+N)*SM(i,j);   
%         end
%         U_top_bar((i-1)*M+1:(i-1)*M+N) = temp;
%     end
% end
U_bottom_bar = zeros(M*M,1);
if(U_bottom_state==1) %另一种写法，更快
    for j = 1:M
        SNU = SM * U_bottom((j-1)*M+1:(j-1)*M+M);
        for i = 1:M
            U_bottom_bar((i-1)*M+1:(i-1)*M+M) = U_bottom_bar((i-1)*M+1:(i-1)*M+M) + SM(i,j) * SNU;
        end
    end
end
% if(U_bottom_state==1)
%     for i=1:M
%         temp=zeros(N,1);
%         for j=1:M
%         temp= temp + SN*U_bottom((j-1)*M+1:(j-1)*M+N)*SM(i,j);   
%         end
%         U_bottom_bar((i-1)*M+1:(i-1)*M+N) = temp;
%     end
% end


% U_left_bar_AK = (aM1×SN×AK) * U0::
% U_left_bar_IK = (aM1×SN×IK) * U0::
U_left_bar_AK = zeros(M*M,1);
U_left_bar_IK = zeros(M*M,1);
if(U_left_state==1) %另一种写法，更快
    for j =1:M
        AKU = AK * U_left((j-1)*M+1:(j-1)*M+M);
        IKU = U_left((j-1)*M+1:(j-1)*M+M);
        for i = 1:M
            U_left_bar_AK((i-1)*M+1:(i-1)*M+M) = U_left_bar_AK((i-1)*M+1:(i-1)*M+M) + AKU * SM(i,j);
            U_left_bar_IK((i-1)*M+1:(i-1)*M+M) = U_left_bar_IK((i-1)*M+1:(i-1)*M+M) + IKU * SM(i,j);
        end
    end
end
% if(U_left_state==1)
%     for i=1:N
%         temp_AK=zeros(K,1);
%         temp_IK=zeros(K,1);
%         for j=1:N
%             temp_AK= temp_AK + AK*U_left((j-1)*N+1:(j-1)*N+K)*SN(i,j);   
%             temp_IK= temp_IK + IK*U_left((j-1)*N+1:(j-1)*N+K)*SN(i,j);   
%         end
%         U_left_bar_AK((i-1)*N+1:(i-1)*N+K) = temp_AK;
%         U_left_bar_IK((i-1)*N+1:(i-1)*N+K) = temp_IK;
%     end
% end
U_right_bar_AK = zeros(M*M,1);
U_right_bar_IK = zeros(M*M,1);
if(U_right_state==1) %另一种写法，更快
    for j =1:M
        AKU = AK * U_right((j-1)*M+1:(j-1)*M+M);
        IKU = U_right((j-1)*M+1:(j-1)*M+M);
        for i = 1:M
            U_right_bar_AK((i-1)*M+1:(i-1)*M+M) = U_right_bar_AK((i-1)*M+1:(i-1)*M+M) + AKU * SM(i,j);
            U_right_bar_IK((i-1)*M+1:(i-1)*M+M) = U_right_bar_IK((i-1)*M+1:(i-1)*M+M) + IKU * SM(i,j);
        end
    end
end
% if(U_right_state==1)
%     for i=1:N
%         temp_AK=zeros(K,1);
%         temp_IK=zeros(K,1);
%         for j=1:N
%             temp_AK= temp_AK + AK*U_right((j-1)*N+1:(j-1)*N+K)*SN(i,j);   
%             temp_IK= temp_IK + IK*U_right((j-1)*N+1:(j-1)*N+K)*SN(i,j);   
%         end
%         U_right_bar_AK((i-1)*N+1:(i-1)*N+K) = temp_AK;
%         U_right_bar_IK((i-1)*N+1:(i-1)*N+K) = temp_IK;
%     end
% end


% U_front_bar = （SM×aN1×IK）U:0:
U_front_bar = zeros(M*M,1);
if(U_front_state==1)
    for j=1:M
       temp = U_front(j:M:(M-1)*M+j);
       SMU = SM * temp;
       for i=1:M
           U_front_bar((i-1)*M+j) = SMU(i);
       end   
    end
end

U_back_bar = zeros(M*M,1);
if(U_back_state==1)
    for j=1:M
       temp = U_back(j:M:(M-1)*M+j);
       SMU = SM * temp;
       for i=1:M
           U_back_bar((i-1)*M+j) = SMU(i);
       end   
    end
end


U_t1_bar = zeros(M,1); % U_t1_bar = (SM×aN1)U: 0 K+1
if(U_t1_state==1)
    U_t1_bar = SM * U_t1;
end
U_t2_bar = zeros(M,1); % U_t2_bar = (aM2×SN)U M+1 : K+1
if(U_t2_state==1)
    U_t2_bar = SM * U_t2;
end
U_t3_bar = zeros(M,1); % U_t3_bar = (SM×aN2)U: N+1 K+1
if(U_t3_state==1)
    U_t3_bar = SM * U_t3;
end
U_t4_bar = zeros(M,1); % U_t4_bar = (aM1×SN)U 0 : K+1
if(U_t4_state==1)
    U_t4_bar = SM * U_t4;
end
U_b1_bar = zeros(M,1); % U_b1_bar = (SM×aN1)U: 0 0
if(U_b1_state==1)
    U_b1_bar = SM * U_b1;
end
U_b2_bar = zeros(M,1); % U_b2_bar = (aM2×SN)U M+1 : 0
if(U_b2_state==1)
    U_b2_bar = SM * U_b2;
end
U_b3_bar = zeros(M,1); % U_b3_bar = (SM×aN2)U: N+1 0
if(U_b3_state==1)
    U_b3_bar = SM * U_b3;
end
U_b4_bar = zeros(M,1); % U_b4_bar = (aM1×SN)U 0 : 0
if(U_b4_state==1)
    U_b4_bar = SM * U_b4;
end
U_l1_bar = zeros(M,1); % U_l1_bar = (aM1×aN1×IK)U 0 0 : 
if(U_l1_state==1)
    U_l1_bar = U_l1;
end
U_l2_bar = zeros(M,1); % U_l2_bar = (aM2×aN1×IK)U M+1 0 : 
if(U_l2_state==1)
    U_l2_bar = U_l2;
end
U_l3_bar = zeros(M,1); % U_l3_bar = (aM2×aN2×IK)U M+1 N+1 : 
if(U_l3_state==1)
    U_l3_bar = U_l3;
end
U_l4_bar = zeros(M,1); % U_l4_bar = (aM1×aN2×IK)U 0 N+1 : 
if(U_l4_state==1)
    U_l4_bar = U_l4;
end



% F_top_bar = （SM×SN）F::K+1 
F_top_bar = zeros(M*M,1);
if(F_top_state==1) %另一种写法，更快
    for j = 1:M
        SNU = SM * F_top((j-1)*M+1:(j-1)*M+M);
        for i = 1:M
            F_top_bar((i-1)*M+1:(i-1)*M+M) = F_top_bar((i-1)*M+1:(i-1)*M+M) + SM(i,j) * SNU;
        end
    end
end
F_bottom_bar = zeros(M*M,1);
if(F_bottom_state==1) %另一种写法，更快
    for j = 1:M
        SNU = SM * F_bottom((j-1)*M+1:(j-1)*M+M);
        for i = 1:M
            F_bottom_bar((i-1)*M+1:(i-1)*M+M) = F_bottom_bar((i-1)*M+1:(i-1)*M+M) + SM(i,j) * SNU;
        end
    end
end
% F_left_bar = (aM1×SN×IK) * F0::
F_left_bar = zeros(M*M,1);
if(F_left_state==1) %另一种写法，更快
    for j =1:M
        IKU = F_left((j-1)*M+1:(j-1)*M+M);
        for i = 1:M
            F_left_bar((i-1)*M+1:(i-1)*M+M) = F_left_bar((i-1)*M+1:(i-1)*M+M) + IKU * SM(i,j);
        end
    end
end
F_right_bar = zeros(M*M,1);
if(F_right_state==1) %另一种写法，更快
    for j =1:M
        IKU = F_right((j-1)*M+1:(j-1)*M+M);
        for i = 1:M
            F_right_bar((i-1)*M+1:(i-1)*M+M) = F_right_bar((i-1)*M+1:(i-1)*M+M) + IKU * SM(i,j);
        end
    end
end
% F_front_bar = （SM×aN1×IK）F:0:
F_front_bar = zeros(M*M,1);
if(F_front_state==1)
    for j=1:M
       temp = F_front(j:M:(M-1)*M+j);
       SMU = SM * temp;
       for i=1:M
           F_front_bar((i-1)*M+j) = SMU(i);
       end   
    end
end
F_back_bar = zeros(M*M,1);
if(F_back_state==1)
    for j=1:M
       temp = F_back(j:M:(M-1)*M+j);
       SMU = SM * temp;
       for i=1:M
           F_back_bar((i-1)*M+j) = SMU(i);
       end   
    end
end



p1 = 1+h*h*K0*K0/12;
p2 = h*h/6;
p4 = K0*K0;

% storage cavity interface
U_gamma_bar = zeros(M*M,1);
for i=1:M
	
	for j=1:M
         % calculation the left term of the equation------------------------------------------------------------------
        %% compute Hij
         temp = (Lambda(i)+Miu(j))*IK + AK;
         Hij = p1*temp + p2*(Lambda(i)*Miu(j)*IK + (Miu(j)+Lambda(i))*AK) + p4*IK;
        
        %% compute Fij_ba & Qij (coefficient matrix of F_bar)----------------------------------
         Fij_ba = zeros(M,1);
         if hasSourceFunction==1
             Qij = temp * h*h/12 + IK;
             for kk = 1:M
                 tempkk = 0;
                 tempSNj = SM(j,:);
                 for ii =1:M
                     x=Xstart + h*ii;
                     z=Zstart + h*kk;
                     tempkk = tempkk + SM(i,ii) * tempSNj * compute_SourceFunction(x,z,para);
                 end
                 Fij_ba(kk) = tempkk;
             end
             Fij_ba = Qij * Fij_ba;
         end
         
        %% compute BUij_bar---------------------------------------------------
         BU_bar = zeros(M,1);
        
         if(U_top_state==1)
            BP_top_bar = zeros(M,1);
            BP_top_bar(M) = ((Lambda(i)+Miu(j))*p2 + p1)*U_top_bar((i-1)*M+j);
            BU_bar = BU_bar + BP_top_bar;
         end
         if(U_bottom_state==1)
            BP_bottom_bar = zeros(M,1);
            BP_bottom_bar(1) = ((Lambda(i)+Miu(j))*p2 + p1)*U_bottom_bar((i-1)*M+j);
            BU_bar = BU_bar + BP_bottom_bar;
         end
         if(U_left_state==1)
            BP_left_bar = SM(i,1) * p2 * U_left_bar_AK((j-1)*M+1:(j*M));
            BP_left_bar = BP_left_bar + SM(i,1) * (Miu(j)*p2 + p1) * U_left_bar_IK((j-1)*M+1:(j*M));
            BU_bar = BU_bar + BP_left_bar;
         end
         if(U_right_state==1)
            BP_right_bar = SM(i,M)* p2 * U_right_bar_AK((j-1)*M+1:(j*M));
            BP_right_bar = BP_right_bar + SM(i,M)*(Miu(j)*p2 + p1)*U_right_bar_IK((j-1)*M+1:(j*M));
            BU_bar = BU_bar + BP_right_bar;
         end
         if(U_front_state==1)
            BP_front_bar = (p2*Lambda(i) + p1) * SM(j,1) * U_front_bar((i-1)*M+1:(i*M));
            BP_front_bar = BP_front_bar + p2 * SM(j,1) * AK * U_front_bar((i-1)*M+1:(i*M));
            BU_bar = BU_bar + BP_front_bar;
         end
         if(U_back_state==1)
            BP_back_bar = (p2*Lambda(i) + p1) * SM(j,M) * U_back_bar((i-1)*M+1:(i*M));
            BP_back_bar = BP_back_bar + p2 * SM(j,M) * AK * U_back_bar((i-1)*M+1:(i*M));
            BU_bar = BU_bar + BP_back_bar;
         end
        
        
         if(U_t1_state==1)
             BE_t1_bar = zeros(M,1);
             BE_t1_bar(M) = p2 * SM(j,1) * U_t1_bar(i);
             BU_bar = BU_bar + BE_t1_bar;
         end
         if(U_t2_state==1)
             BE_t2_bar = zeros(M,1);
             BE_t2_bar(M) = p2 * SM(i,M) * U_t2_bar(j);
             BU_bar = BU_bar + BE_t2_bar;
         end
         if(U_t3_state==1)
             BE_t3_bar = zeros(M,1);
             BE_t3_bar(M) = p2 * SM(j,M) * U_t3_bar(i);
             BU_bar = BU_bar + BE_t3_bar;
         end
         if(U_t4_state==1)
             BE_t4_bar = zeros(M,1);
             BE_t4_bar(M) = p2 * SM(i,1) * U_t4_bar(j);
             BU_bar = BU_bar + BE_t4_bar;
         end
         if(U_b1_state==1)
             BE_b1_bar = zeros(M,1);
             BE_b1_bar(1) = p2 * SM(j,1) * U_b1_bar(i);
             BU_bar = BU_bar + BE_b1_bar;
         end
         if(U_b2_state==1)
             BE_b2_bar = zeros(M,1);
             BE_b2_bar(1) = p2 * SM(i,M) * U_b2_bar(j);
             BU_bar = BU_bar + BE_b2_bar;
         end
         if(U_b3_state==1)
             BE_b3_bar = zeros(M,1);
             BE_b3_bar(1) = p2 * SM(j,M) * U_b3_bar(i);
             BU_bar = BU_bar + BE_b3_bar;
         end
         if(U_b4_state==1)
             BE_b4_bar = zeros(M,1);
             BE_b4_bar(1) = p2 * SM(i,1) * U_b4_bar(j);
             BU_bar = BU_bar + BE_b4_bar;
         end
         if(U_l1_state==1)
             BE_l1_bar = p2 * SM(i,1) * SM(j,1) * U_l1_bar;
             BU_bar = BU_bar + BE_l1_bar;
         end
         if(U_l2_state==1)
             BE_l2_bar = p2 * SM(i,M) * SM(j,1) * U_l2_bar;
             BU_bar = BU_bar + BE_l2_bar;
         end
         if(U_l3_state==1)
             BE_l3_bar = p2 * SM(i,M) * SM(j,M) * U_l3_bar;
             BU_bar = BU_bar + BE_l3_bar;
         end
         if(U_l4_state==1)
             BE_l4_bar = p2 * SM(i,1) * SM(j,M) * U_l4_bar;
             BU_bar = BU_bar + BE_l4_bar;
         end
        
        %% compute BFij_bar---------------------------------------------------
         BF_bar = zeros(M,1);
         
         if(F_top_state==1)
            BF_top_bar = zeros(M,1);
            BF_top_bar(M) = p1 * F_top_bar((i-1)*M+j);
            BF_bar = BF_bar + BF_top_bar;
         end
         if(F_bottom_state==1)
            BF_bottom_bar = zeros(M,1);
            BF_bottom_bar(1) = p1 * F_bottom_bar((i-1)*M+j);
            BF_bar = BF_bar + BF_bottom_bar;
         end
         if(F_left_state==1)
            BF_left_bar = SM(i,1) * p1 * F_left_bar((j-1)*M+1:(j*M));
            BF_bar = BF_bar + BF_left_bar;
         end
         if(F_right_state==1)
            BF_right_bar = SM(i,M) * p1 * F_right_bar((j-1)*M+1:(j*M));
            BF_bar = BF_bar + BF_right_bar;
         end
        if(F_front_state==1)
            BF_front_bar = p1 * SM(j,1) * F_front_bar((i-1)*M+1:(i*M));
            BF_bar = BF_bar + BF_front_bar;
         end
         if(F_back_state==1)
            BF_back_bar = p1 * SM(j,M) * F_back_bar((i-1)*M+1:(i*M));
            BF_bar = BF_bar + BF_back_bar;
         end
        
		%% compute F^
		 right_term =Fij_ba +BF_bar -BU_bar/(h*h);

		 % solve the linear equation set : get the value of Uij:_ba
		 % compute Uij:_ba
         % BarV storage the value of Uij:_ba
		 [BarV,flag,relres,iter]=bicgstab(Hij,right_term,1e-14,1800);
		
         U_gamma_bar((i-1)*M+j) = BarV(extractZ);
	end
end

%--------------------------end the timer-----------------------------------


% U = SMN * U_bar
U = zeros(M*M,1);
for i=1:M,
    for j=1:M
        temp =0;
        for k=1:M
            U((i-1)*M+j)= U((i-1)*M+j)+(SM(i,k)*SM(j,:)) * U_gamma_bar((k-1)*M+1:k*M);
        end
    end
    
end
toc


% get and draw the cavity interface
Cavity_interface=zeros(M,M);
for i=1:M
	for j=1:M
		Cavity_interface(j,i)=U((i-1)*M+j);
	end
end
plot2Dsolution( para, extractZ, Cavity_interface);


% draw the real value of U
real_U_gama=zeros(M,M);
for i=1:M
	for j=1:M
		x=h*i;
		y=h*j;
		z=h*extractZ;
		real_U_gama(j,i)=compute_realU( x,y,z );	
	end
end
plot2Dsolution( para, extractZ, real_U_gama);


% compute error
ErrorM_abs=abs(Cavity_interface-real_U_gama);
tempE=sum(sum(ErrorM_abs.*ErrorM_abs));
tempE=tempE*(Xend-Xstart)*(Yend-Ystart)/(M*M);
e2=sqrt(tempE)

% % compute U_ba error
% ErrorM=real_U_ba-U_ba;
% ErrorM(1:10)
% tempE=sum(abs(ErrorM));
% tempE=tempE*(Xend-Xstart)*(Yend-Ystart)*(Zend-Zstart)/(M*N*K);
% e2=sqrt(tempE)
