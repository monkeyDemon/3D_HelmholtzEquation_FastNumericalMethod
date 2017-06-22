function [ BF ] = computeBF_3d_4o_di( i,j,M,N,K,h,K0 )
%Calculate the bound of F_ij:

BF=zeros(K,1);

para1=(1+h*h*K0*K0/12)/(h*h);

if( i==1 && j==1  )
	for k=1:K
		if(k==1)
			temp=compute_F_ijk(h,0,j,k)+compute_F_ijk(h,i,0,k)+compute_F_ijk(h,i,j,0);
			BF(k)=BF(k)+temp*para1;
		elseif(k==K)
			temp=compute_F_ijk(h,0,j,k)+compute_F_ijk(h,i,0,k)+compute_F_ijk(h,i,j,K+1);
			BF(k)=BF(k)+temp*para1;
		else
			temp=compute_F_ijk(h,0,j,k)+compute_F_ijk(h,i,0,k);
			BF(k)=BF(k)+temp*para1;
		end
	end
elseif( i==1 && j==N  )
	for k=1:K
		if(k==1)
			temp=compute_F_ijk(h,i-1,j,k)+compute_F_ijk(h,i,j+1,k)+compute_F_ijk(h,i,j,k-1);
			BF(k)=BF(k)+temp*para1;
		elseif(k==K)
			temp=compute_F_ijk(h,i-1,j,k)+compute_F_ijk(h,i,j+1,k)+compute_F_ijk(h,i,j,k+1);
			BF(k)=BF(k)+temp*para1;
		else
			temp=compute_F_ijk(h,i-1,j,k)+compute_F_ijk(h,i,j+1,k);
			BF(k)=BF(k)+temp*para1;
		end
	end
elseif( i==1 && j>1 && j<N  )
	for k=1:K
		if(k==1)
			temp=compute_F_ijk(h,i-1,j,k)+compute_F_ijk(h,i,j,k-1);
			BF(k)=BF(k)+temp*para1;
		elseif(k==K)
			temp=compute_F_ijk(h,i-1,j,k)+compute_F_ijk(h,i,j,k+1);
			BF(k)=BF(k)+temp*para1;
		else
			temp=compute_F_ijk(h,i-1,j,k);
			BF(k)=BF(k)+temp*para1;
		end
	end
elseif( i==M && j==1 )
	for k=1:K
		if(k==1)
			temp=compute_F_ijk(h,i+1,j,k)+compute_F_ijk(h,i,j-1,k)+compute_F_ijk(h,i,j,k-1);
			BF(k)=BF(k)+temp*para1;
		elseif(k==K)
			temp=compute_F_ijk(h,i+1,j,k)+compute_F_ijk(h,i,j-1,k)+compute_F_ijk(h,i,j,k+1);
			BF(k)=BF(k)+temp*para1;
		else
			temp=compute_F_ijk(h,i+1,j,k)+compute_F_ijk(h,i,j-1,k);
			BF(k)=BF(k)+temp*para1;
		end
	end
elseif( i==M && j==N )
	for k=1:K
		if(k==1)
			temp=compute_F_ijk(h,i+1,j,k)+compute_F_ijk(h,i,j+1,k)+compute_F_ijk(h,i,j,k-1);
			BF(k)=BF(k)+temp*para1;
		elseif(k==K)
			temp=compute_F_ijk(h,i+1,j,k)+compute_F_ijk(h,i,j+1,k)+compute_F_ijk(h,i,j,k+1);
			BF(k)=BF(k)+temp*para1;
		else
			temp=compute_F_ijk(h,i+1,j,k)+compute_F_ijk(h,i,j+1,k);
			BF(k)=BF(k)+temp*para1;
		end
	end
elseif( i==M && j>1 && j<N  )
	for k=1:K
		if(k==1)
			temp=compute_F_ijk(h,i+1,j,k)+compute_F_ijk(h,i,j,k-1);
			BF(k)=BF(k)+temp*para1;
		elseif(k==K)
			temp=compute_F_ijk(h,i+1,j,k)+compute_F_ijk(h,i,j,k+1);
			BF(k)=BF(k)+temp*para1;
		else
			temp=compute_F_ijk(h,i+1,j,k);
			BF(k)=BF(k)+temp*para1;
		end
	end
elseif( i>1 && i<M && j==1  )
	for k=1:K
		if(k==1)
			temp=compute_F_ijk(h,i,j-1,k)+compute_F_ijk(h,i,j,k-1);
			BF(k)=BF(k)+temp*para1;
		elseif(k==K)
			temp=compute_F_ijk(h,i,j-1,k)+compute_F_ijk(h,i,j,k+1);
			BF(k)=BF(k)+temp*para1;
		else
			temp=compute_F_ijk(h,i,j-1,k);
			BF(k)=BF(k)+temp*para1;
		end
	end
elseif( i>1 && i<M && j==N   )
	for k=1:K
		if(k==1)
			temp=compute_F_ijk(h,i,j+1,k)+compute_F_ijk(h,i,j,k-1);
			BF(k)=BF(k)+temp*para1;
		elseif(k==K)
			temp=compute_F_ijk(h,i,j+1,k)+compute_F_ijk(h,i,j,k+1);
			BF(k)=BF(k)+temp*para1;
		else
			temp=compute_F_ijk(h,i,j+1,k);
			BF(k)=BF(k)+temp*para1;
		end
	end
elseif( i>1 && i<M && j>1 && j<N )
	for k=1:K
		if(k==1)
			temp=compute_F_ijk(h,i,j,k-1);
			BF(k)=BF(k)+temp*para1;
		elseif(k==K)
			temp=compute_F_ijk(h,i,j,k+1);
			BF(k)=BF(k)+temp*para1;
		else
			% do nothing
		end
	end
end

%
%for k=1:K
%	% find the bound of Z 
%
%	index=(i-1)*M+(j-1)*N+k;
%
%	if k==1
%		% case one: (k-1 = 0)
%		BU_normal(index)=BU_normal(index)+compute_BoundaryCondition(i,j,0)*para1;
%		temp=-4*compute_BoundaryCondition(i,j,0)+compute_BoundaryCondition(i,j-1,0)+compute_BoundaryCondition(i,j+1,0)+compute_BoundaryCondition(i-1,j,0)+compute_BoundaryCondition(i+1,j,0);
%		BU_normal(index)=BU_normal(index)+temp*para2;
%	elseif k==K
%		% case one: (k+1= K+1)
%		BU_gama(index)=BU_gama(index)+compute_BoundaryCondition(i,j,K+1)*para1;
%		temp=-4*compute_BoundaryCondition(i,j,K+1)+compute_BoundaryCondition(i,j-1,K+1)+compute_BoundaryCondition(i,j+1,K+1)+compute_BoundaryCondition(i-1,j,K+1)+compute_BoundaryCondition(i+1,j,K+1);
%		BU_gama(index)=BU_gama(index)+temp*para2;
%	else
%		% others: not bound of Z,but we should check X Y axis
%	end
%end
%
%
%% we have look over all the bound of Z
%% next we need to check X and Y
%if i=1
%	% case one: (i-1 = 0)
%	for k=1:K
%		index=(i-1)*M+(j-1)*N+k;
%		if(k==1)
%			BU_normal(index)=BU_normal(index)+compute_BoundaryCondition(0,j,k)*para1;
%			temp=-4*compute_BoundaryCondition(0,j,k)+compute_BoundaryCondition(0,j-1,k)+compute_BoundaryCondition(0,j+1,k)+compute_BoundaryCondition(0,j,k+1);
%			BU_normal(index)=BU_normal(index)+temp*para2;
%		elseif(k==K)
%			BU_normal(index)=BU_normal(index)+compute_BoundaryCondition(0,j,k)*para1;
%			temp=-4*compute_BoundaryCondition(0,j,k)+compute_BoundaryCondition(0,j-1,k)+compute_BoundaryCondition(0,j+1,k)+compute_BoundaryCondition(0,j,k-1);
%			BU_normal(index)=BU_normal(index)+temp*para2;
%		else
%			BU_normal(index)=BU_normal(index)+compute_BoundaryCondition(0,j,k)*para1;
%			temp=-4*compute_BoundaryCondition(0,j,k)+compute_BoundaryCondition(0,j-1,k)+compute_BoundaryCondition(0,j+1,k)+compute_BoundaryCondition(0,j,k-1)+compute_BoundaryCondition(0,j,k+1);
%			BU_normal(index)=BU_normal(index)+temp*para2;
%		end
%	end
%elseif i=M
%	% case one: (i+1= M+1)
%	for k=1:K
%		index=(i-1)*M+(j-1)*N+k;
%		if(k==1)
%			BU_normal(index)=BU_normal(index)+compute_BoundaryCondition(M+1,j,k)*para1;
%			temp=-4*compute_BoundaryCondition(M+1,j,k)+compute_BoundaryCondition(M+1,j-1,k)+compute_BoundaryCondition(M+1,j+1,k)+compute_BoundaryCondition(M+1,j,k+1);
%			BU_normal(index)=BU_normal(index)+temp*para2;
%		elseif(k==K)
%			BU_normal(index)=BU_normal(index)+compute_BoundaryCondition(M+1,j,k)*para1;
%			temp=-4*compute_BoundaryCondition(M+1,j,k)+compute_BoundaryCondition(M+1,j-1,k)+compute_BoundaryCondition(M+1,j+1,k)+compute_BoundaryCondition(M+1,j,k-1);
%			BU_normal(index)=BU_normal(index)+temp*para2;
%		else
%			BU_normal(index)=BU_normal(index)+compute_BoundaryCondition(M+1,j,k)*para1;
%			temp=-4*compute_BoundaryCondition(M+1,j,k)+compute_BoundaryCondition(M+1,j-1,k)+compute_BoundaryCondition(M+1,j+1,k)+compute_BoundaryCondition(M+1,j,k-1)+compute_BoundaryCondition(M+1,j,k+1);
%			BU_normal(index)=BU_normal(index)+temp*para2;
%		end
%	end
%else
%	% others: not bound of X and Z,but we stil need to check Y axis
%	% remember: i=2:M-1 for this part, we don't need to consider i
%	if j=1
%		% case one: (j-1 = 0)
%		for k=1:K
%			index=(i-1)*M+(j-1)*N+k;
%			if(k==1)
%				BU_normal(index)=BU_normal(index)+compute_BoundaryCondition(i,0,k)*para1;
%				temp=-4*compute_BoundaryCondition(i,0,k)+compute_BoundaryCondition(i,0,k+1)+compute_BoundaryCondition(i-1,0,k)+compute_BoundaryCondition(i+1,0,k);
%				BU_normal(index)=BU_normal(index)+temp*para2;
%			elseif(k==K)
%				BU_normal(index)=BU_normal(index)+compute_BoundaryCondition(i,0,k)*para1;
%				temp=-4*compute_BoundaryCondition(i,0,k)+compute_BoundaryCondition(i,0,k-1)+compute_BoundaryCondition(i-1,0,k)+compute_BoundaryCondition(i+1,0,k);
%				BU_normal(index)=BU_normal(index)+temp*para2;
%			else
%				BU_normal(index)=BU_normal(index)+compute_BoundaryCondition(i,0,k)*para1;
%				temp=-4*compute_BoundaryCondition(i,0,k)+compute_BoundaryCondition(i,0,k-1)+compute_BoundaryCondition(i,0,k+1)+compute_BoundaryCondition(i-1,0,k)+compute_BoundaryCondition(i+1,0,k);
%				BU_normal(index)=BU_normal(index)+temp*para2;
%			end
%		end
%	elseif j=N
%		% case one: (j+1= N+1)
%		for k=1:K
%			index=(i-1)*M+(j-1)*N+k;
%			if(k==1)
%				BU_normal(index)=BU_normal(index)+compute_BoundaryCondition(i,N+1,k)*para1;
%				temp=-4*compute_BoundaryCondition(i,N+1,k)+compute_BoundaryCondition(i,N+1,k+1)+compute_BoundaryCondition(i-1,N+1,k)+compute_BoundaryCondition(i+1,N+1,k);
%				BU_normal(index)=BU_normal(index)+temp*para2;
%			elseif(k==K)
%				BU_normal(index)=BU_normal(index)+compute_BoundaryCondition(i,N+1,k)*para1;
%				temp=-4*compute_BoundaryCondition(i,N+1,k)+compute_BoundaryCondition(i,N+1,k-1)+compute_BoundaryCondition(i-1,N+1,k)+compute_BoundaryCondition(i+1,N+1,k);
%				BU_normal(index)=BU_normal(index)+temp*para2;
%			else
%				BU_normal(index)=BU_normal(index)+compute_BoundaryCondition(i,N+1,k)*para1;
%				temp=-4*compute_BoundaryCondition(i,N+1,k)+compute_BoundaryCondition(i,N+1,k-1)+compute_BoundaryCondition(i,N+1,k+1)+compute_BoundaryCondition(i-1,N+1,k)+compute_BoundaryCondition(i+1,N+1,k);
%				BU_normal(index)=BU_normal(index)+temp*para2;
%			end
%		end
%	else
%		% for thie part: i=2:M-1 j=2:N-1 
%		% will not produce bound ,do nothing	
%	end
%end

end
