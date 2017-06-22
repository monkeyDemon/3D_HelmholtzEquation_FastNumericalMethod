function[ F_ij ]=compute_F_col(i,j,M,N,K,h)

F_ij=zeros(K,1);
for k=1:K
	F_ij(k)=compute_F_ijk(h,i,j,k);
end

