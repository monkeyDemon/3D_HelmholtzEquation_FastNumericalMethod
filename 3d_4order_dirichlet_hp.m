% author: AnSheng
% this code can sovle the 3 dimentional Helmholtz fuction with dirichlet boundary condition
% use a 4 order fast numerical method

% this demo has two versions: hp and sm
% hp means the code solving with the fast performance as it can ,but it can't deal with the probelm with large number grid division M,N,K
% sm means the code solving as the mode of saving memory,but it will run more slowly 

% if you want to run our code on your numerical examples
% these functions need to be modified:
% 3d_4order_dirichlet_sm.m         the main function of save memory version
% 3d_4order_dirichlet_hp.m         the main function of high performance version
% compute_BoundaryCondition.m      the function to compute boundary value of u
% compute_F_ijk.m                  the function to compute f
% you should change M,N,K and Xstart,Ystart,Zstart according to the real situation

% here the sm version code is below:
% -----------------------------------------------------------------------------------------

k0=pi;
%epr=1+i;
epr=0;
K0=k0^2*epr;

%-----------------need to change on different numerical examples-------------------
M=31;
N=31;
K=31;

Xstart=0;
Xend=1;
Ystart=0;
Yend=1;
Zstart=0;
Zend=1;
%-----------------need to change on different numerical examples-------------------


% in order to simplified the derivation process, in our code we also have hx=hy=hz
hx=(Xend-Xstart)/(M+1);
hy=(Yend-Ystart)/(N+1);
hz=(Zend-Zstart)/(K+1);
h=hx; 

% draw the real value of U
extractZ=K;% set which layer we extract ,then we can compare with the true value of u
real_U_gama=zeros(M,N);
for i=1:M
	for j=1:N
		x=h*i;
		y=h*j;
		z=h*extractZ;
		real_U_gama(i,j)=(sin(pi*y)*sin(pi*z)/sinh(sqrt(2)*pi))*(2*sinh(sqrt(2)*pi*x)+sinh(sqrt(2)*pi*(1-x)));	
	end
end
figure
mesh(Xstart+hx:hx:Xend-hx,Ystart+hy:hy:Yend-hy,real_U_gama);

%--------------------------start the timer-----------------------------------
tic

Lambda=2*(M+1)*sin((1:M)*pi/(2*(M+1)));
Lambda=-Lambda.*Lambda;

Miu=2*(N+1)*sin((1:N)*pi/(2*(N+1)));
Miu=-Miu.*Miu;

P=2*(K+1)*sin((1:K)*pi/(2*(K+1)));
P=-P.*P;

% compute FFT matrix
SM=FstMatrix(M);
SN=FstMatrix(N);
SK=FstMatrix(K);

% compute SM°ËSN°ËSK
SMNK=zeros(M*N*K,M*N*K);
for i=1:M*N
	for j=1:M*N
		ni=mod(i,N);
        if ni==0
            ni=N;
        end
        mi=(i-ni)/N+1;
        nj=mod(j,N);
        if nj==0
            nj=N;
        end
        mj=(j-nj)/N+1;

        SMNK((i-1)*K+1:(i-1)*K+K,(j-1)*K+1:(j-1)*K+K)=SM(mi,mj)*SN(ni,nj)*SK;
    end
end


% storage B_ba
B=zeros(M*N*K,1);
B_ba=zeros(M*N*K,1);
for bi=1:M
	for bj=1:N
		B((bi-1)*N*K+(bj-1)*K+1 : (bi-1)*N*K+(bj-1)*K+K)=computeBU_3d_4o_di(bi,bj,M,N,K,h,K0);
    end
end
B_ba=SMNK*B;
clear B;

% compute F_ij_ba
F_ba=zeros(M*N*K,1);
F=zeros(M*N*K,1);
for bi=1:M
	for bj=1:N
        F((bi-1)*N*K+(bj-1)*K+1 : (bi-1)*N*K+(bj-1)*K+K)=compute_F_col(bi,bj,M,N,K,h);
	end
end
F_ba=SMNK*F;
clear F;

% compute BF_ij_ba
BF=zeros(M*N*K,1);
BF_ba=zeros(M*N*K,1);
for bi=1:M
	for bj=1:N
        BF((bi-1)*N*K+(bj-1)*K+1 : (bi-1)*N*K+(bj-1)*K+K)=computeBF_3d_4o_di(bi,bj,M,N,K,h,K0);
	end
end
BF_ba=SMNK*BF;
clear BF;


% storage cavity interface
U_ba=zeros(M*N*K,1);
U=zeros(M*N*K,1);
for i=1:M
	i
	for j=1:N
		% compute linear equation set :similar as (4-2-9).Zou P32,see my derivation for the detail
		
		% calculation the left term of the equation------------------------------------------------------------------
		% compute Hij
		Hij=sparse(1:K,1:K,(Lambda(i)+Miu(j)+P)*(1+h*h*K0*K0/12)+(h*h/6)*(Lambda(i)*Miu(j)+(Miu(j)+Lambda(i))*P)+K0^2);

		% compute B_ij_ba
		B_ij_ba=B_ba( (i-1)*N*K+(j-1)*K+1 :(i-1)*N*K+(j-1)*K+K );


		% calculation the right term of the equation------------------------------------------------------------------
		% compute F_ij_ba
		F_ij_ba=F_ba( (i-1)*N*K+(j-1)*K+1 :(i-1)*N*K+(j-1)*K+K );


		% compute BF_ij_ba
		BF_ij_ba=BF_ba( (i-1)*N*K+(j-1)*K+1 :(i-1)*N*K+(j-1)*K+K );


		% compute F_ij_all (all the right term of the equation)
		HF_ij=sparse(1:K,1:K,(Lambda(i)+Miu(j)+P)*h*h/12);
		F_ij_all=F_ij_ba+HF_ij*F_ij_ba+BF_ij_ba;
		
		%compute F^
		F_jian=F_ij_all-B_ij_ba;

		% solve the linear equation set : get the value of Uij:_ba
		% compute Uij:_ba
        % BarV storage the value of Uij:_ba
		[BarV,flag,relres,iter]=bicgstab(Hij,F_jian,1e-14,1200);
		
        U_ba( (i-1)*N*K+(j-1)*K+1 :(i-1)*N*K+(j-1)*K+K )=BarV;
        
		j
	end
end

U=SMNK*U_ba;
clear SMNK;
clear U_ba;

Cavity_interface=zeros(M,N);
% get the cavity interface
for i=1:M
	for j=1:N
		Cavity_interface(i,j)=U((i-1)*N*K+(j-1)*K+extractZ);
	end
end


Cavity_interface_mag=abs(Cavity_interface);

% compute error
ErrorM=Cavity_interface_mag-real_U_gama;
tempE=sum(sum(abs(ErrorM)));
tempE=tempE*(Xend-Xstart)*(Yend-Ystart)/(M*N);
e2=sqrt(tempE)

toc
%--------------------------end the timer-----------------------------------


figure
mesh(Xstart+hx:hx:Xend-hx,Ystart+hy:hy:Yend-hy,Cavity_interface_mag);

pause

