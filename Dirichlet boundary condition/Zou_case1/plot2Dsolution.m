function [  ] = plot2Dsolution( para, extractZ, U_2D)
% 绘制3D区域内某个2D平面的数值解

M = para(1);
h = para(2);
Xstart = para(3);
Xend = para(4);
Ystart = para(5);
Yend = para(6);
Zstart = para(7);
Zend = para(8);


figure
[Xmesh,Ymesh] = meshgrid( Xstart+h:h:Xend-h, Ystart+h:h:Yend-h);
Zmesh = ones(M,M)*h*extractZ;
surf(Xmesh,Ymesh,Zmesh,U_2D);
xlabel('x');ylabel('y');zlabel('z');
colorbar;
colormap(jet(128));
shading interp;  
set(gca, 'xlim', [Xstart,Xend]);
set(gca, 'ylim', [Ystart,Yend]);
set(gca, 'zlim', [Zstart,Zend]);


end

