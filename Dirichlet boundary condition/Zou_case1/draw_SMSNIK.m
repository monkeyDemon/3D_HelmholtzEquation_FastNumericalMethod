clc,clear
close all

M=3;
N=3;
K=3;

% compute FFT matrix
SM=FstMatrix(M);
SN=FstMatrix(N);
SK=FstMatrix(K);

SMNK=zeros(M*N*K,M*N*K);
IK = eye(K);
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

        SMNK((i-1)*K+1:(i-1)*K+K,(j-1)*K+1:(j-1)*K+K)=SM(mi,mj)*SN(ni,nj)*IK;
    end
end
spy(SMNK);
hold on
xlabel('')
set(gca,'xtick',[])
set(gca,'ytick',[])
title('structure of $ S_M \bigotimes S_N \bigotimes I_K $','interpreter','latex')

% 画横向分割线
Bold = 1;
plot([0,28],[0,0],'k','linewidth',1.5);
plot([0,28],[28,28],'k','linewidth',1.5);
for t = 3.5:3:24.5
    leftx = 0;
    lefty = t;
    rightx = 28;
    righty = t;
    if Bold==3
        plot([leftx,rightx],[lefty,righty],'k','linewidth',1.5);
        Bold = 0;
    else
        plot([leftx,rightx],[lefty,righty],'k','linewidth',0.5);
    end
    Bold = Bold+1;
end

% 画纵向分割线
Bold = 1;
plot([0,0],[0,28],'k','linewidth',1.5);
plot([28,28],[0,28],'k','linewidth',1.5);
for t = 3.5:3:24.5
    leftx = t;
    lefty = 0;
    rightx = t;
    righty = 28;
    if Bold==3
        plot([leftx,rightx],[lefty,righty],'k','linewidth',1.5);
        Bold = 0;
    else
        plot([leftx,rightx],[lefty,righty],'k','linewidth',0.5);
    end
    Bold = Bold+1;
end