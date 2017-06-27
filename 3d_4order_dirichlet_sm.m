% author: AnSheng
% this code is aimed to sovle the 3 dimentional Helmholtz fuction with dirichlet boundary condition
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

% here the hp version code is below:
% -----------------------------------------------------------------------------------------

k0=pi;
%epr=1+i;
epr=0;
K0=k0^2*epr;

%-----------------need to change on different numerical examples-------------------
M=63;
N=63;
K=63;

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
		x=hx*i;
		y=hy*j;
		z=hz*extractZ;
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

% storage cavity interface
Cavity_interface=zeros(M,N);
for i=1:M
	i
	for j=1:N
		% compute linear equation set :similar as (4-2-9).Zou P32,see my derivation for the detail
		
		% calculation the left term of the equation------------------------------------------------------------------
		% compute Hij
		Hij=sparse(1:K,1:K,(Lambda(i)+Miu(j)+P)*(1+h*h*K0*K0/12)+(h*h/6)*(Lambda(i)*Miu(j)+(Miu(j)+Lambda(i))*P)+K0^2);

		% compute B_ij_ba
		B_ij_ba=zeros(K,1);% bound of i,j,: after fft
		B_ij=zeros(K,1);
		for bi=1:M
			for bj=1:N
				B_ij=computeBU_3d_4o_di(bi,bj,M,N,K,h,K0);% bound of bi,bj,: before fft
				B_ij_ba=B_ij_ba+SM(i,bi)*SN(j,bj)*SK*B_ij;	
			end
		end


		% calculation the right term of the equation------------------------------------------------------------------
		% compute F_ij_ba
		F_ij_ba=zeros(K,1);% F's bound of i,j,: after fft
		F_ij=zeros(K,1);
		for bi=1:M
			for bj=1:N
				F_ij=compute_F_col(bi,bj,M,N,K,h);% bound of bi,bj,: before fft
				F_ij_ba=F_ij_ba+SM(i,bi)*SN(j,bj)*SK*F_ij;	
			end
		end

		% compute BF_ij_ba
		BF_ij=zeros(K,1);% bound of bi,bj,: before fft
		BF_ij_ba=zeros(K,1);
		for bi=1:M
			for bj=1:N
				BF_ij=computeBF_3d_4o_di(bi,bj,M,N,K,h,K0);% bound of bi,bj,: before fft
				BF_ij_ba=BF_ij_ba+SM(i,bi)*SN(j,bj)*SK*BF_ij;	
			end
		end

		% compute F_ij_all (all the right term of the equation)
		HF_ij=sparse(1:K,1:K,(Lambda(i)+Miu(j)+P)*h*h/12);
		F_ij_all=F_ij_ba+HF_ij*F_ij_ba+BF_ij_ba;
		
		%compute F^
		F_jian=F_ij_all-B_ij_ba;

		% solve the linear equation set : get the value of Uij:_ba
		% compute Uij:_ba
		[BarV,flag,relres,iter]=bicgstab(Hij,F_jian,1e-14,1200);

		% BarV storage the value of Uij:_ba


		% save in cavity interface
		for bi=1:M
			for bj=1:N
				Cavity_interface(bi,bj)=Cavity_interface(bi,bj)+SM(bi,i)*SN(bj,j)*SK(extractZ,:)*BarV;
			end
		end

		j
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

