[X,Y] = meshgrid(-0.5:0.2:1,-0.5:0.2:1);
[X2,Y2] = meshgrid(-1:0.1:1.5,-1:0.1:1.5);
%%
P = sqrt(6-X.^2-Y.^2)-0.5;
[U,V,W] = surfnorm(X,Y,P);

%%
figure(9)
clf
set(gcf, 'PaperSize',[10 8])
set(gcf, 'PaperPosition', [0 0 10 8])
quiver3(X,Y,P,U,V,W,'k');
hold on
mesh(X2,Y2,3+exp(-(X2.^2+Y2.^2)))
hold off
axis tight
grid off
axis off

saveas(gcf,'./plots/Surface_P.pdf')
%%
[X,Y] = meshgrid(-0.7:0.2:1,-1:0.2:1);
[X2,Y2] = meshgrid(-1:0.1:1.5,-1:0.1:1.5);
S = 0.01*sin(5*X);

figure(9)
clf
set(gcf, 'PaperSize',[10 8])
set(gcf, 'PaperPosition', [0 0 10 8])
quiver3(X,Y,zeros(size(X)),zeros(size(X)),zeros(size(X)),S,'k');
hold on
mesh(X2,Y2,1+0.05*sin(5*X2))
hold off
axis tight
grid off
axis off

saveas(gcf,'./plots/Surface_S.pdf')