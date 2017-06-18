
k0=pi;
%epr=1+i;
epr=0;
K0=k0^2*epr;
M=15;
N=15;
K=15;
h=1/(M+1);

extractZ=5;% set which layer we extract ,then we can compare with the true value of u
real_U_gama=zeros(M,N);
for i=1:M
	for j=1:N
		x=h*i;
		y=h*j;
		z=h*extractZ;
		real_U_gama(i,j)=(sin(pi*y)*sin(pi*z)/sin(sqrt(2)*pi))*(2*sinh(sqrt(2)*pi*x)+sinh(sqrt(2)*pi*(1-x)));	
	end
end
figure
mesh(h:h:1-h,h:h:1-h,real_U_gama);


Lambda=2*(M+1)*sin((1:M)*pi/(2*(M+1)));
Lambda=-Lambda.*Lambda;

Miu=2*(N+1)*sin((1:N)*pi/(2*(N+1)));
Miu=-Miu.*Miu;

P=2*(K+1)*sin((1:K)*pi/(2*(K+1)));
P=-P.*P;
	
% compute P
%Pmatrix=sparse(1:K,1:K,P);

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
				F_ij=compute_F_col(bi,bj,M,N,K);% bound of bi,bj,: before fft
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

%Cavity_interface_mag=abs(Cavity_interface);
Cavity_interface

figure
%mesh(h:h:1-h,h:h:1-h,Cavity_interface_mag);
mesh(h:h:1-h,h:h:1-h,Cavity_interface);

pause

% compute U (Transform U_ba to U use FFT)

