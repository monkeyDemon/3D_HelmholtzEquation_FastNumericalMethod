
k0=pi;
epr=1+i;
K0=k0^2*epr;
M=32;
N=32;
K=32;
h=1/(M+1);

Lambda=2*(M+1)*sin((1:M)*pi/(2*(M+1)));
Lambda=-Lambda.*Lambda;

Miu=2*(N+1)*sin((1:N)*pi/(2*(N+1)));
Miu=-Miu.*Miu;

P=2*(K+1)*sin((1:K)*pi/(2*(K+1)));
P=-P.*P;
	
% compute P
Pmatrix=sparse(1:K,1:K,P);
% storage cavity interface
Cavity_interface=zeros(M,N);
for i=1:M
	for j=1:N
		% compute (4-2-9) JingWen.Zou P32
		
		% compute P (we have done before)

		% compute Hj
		temp1=sparse(1:K,1:K,Lambda(i)+Miu(j))+Pmatrix;
		temp2=sparse(1:K,1:K,(h*h/6)*(1-(K0^2)*(h^2)/30)*Lambda(i)*Miu(j))+(Miu(j)+Lambda(i))*Pmatrix;
		temp3=(h^4)/30*Lambda(i)*Miu(j)*Pmatrix;
		temp4=sparse(1:K,1:K,-(K0^2)*(1+(h^2*K0^2)/12+(h^4*K0^4)/360));
		Hj=temp1+temp2+temp3+temp4;

		% compute F_ba
		F=zeros(K,1);

		% compute B_ba
		B=ones(K,1);
		
		%compute F^
		F_j=F-B;

		% compute Uij:_ba
		[BarV,flag,relres,iter]=bicgstab(Hj,F_j,1e-14,1200);

		% BarV storage the value of Uij:_ba
		% save in cavity interface
		Cavity_interface(i,j)=BarV(K);
		i
	end
end

Cavity_interface_mag=abs(Cavity_interface);
figure
mesh(h:h:1-h,h:h:1-h,Cavity_interface_mag);

pause

% compute U (Transform U_ba to U use FFT)
SM=FstMatrix(M);
SN=FstMatrix(N);
SK=FstMatrix(K);

