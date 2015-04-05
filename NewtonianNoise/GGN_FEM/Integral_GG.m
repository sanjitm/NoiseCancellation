set(0,'DefaultAxesFontSize',22);
set(0,'DefaultTextFontSize',22);

%%
x = linspace(0,3,200);

dipole = 1/3+cos(2*pi*x)./(2*pi*x).^2-sin(2*pi*x)./(2*pi*x).^3;

s2 = 0.5*(sin(2*pi*x)./(2*pi*x)+2*cos(2*pi*x)./(2*pi*x).^2 ...
    -2*sin(2*pi*x)./(2*pi*x).^3);

figure(14)
set(gcf, 'PaperSize',[8 5])
set(gcf, 'PaperPosition', [0 0 8 5])
plot(x,dipole,'k',x,dipole-s2,'k--',...
    output.dist_GG/1000,output.conv_GG(:,1)/nr(1),'ko','MarkerSize',5,'LineWidth',2)
grid
xlabel('Distance [r/\lambda]')
ylabel('Gravity-Gradient Integral')
legend('Spherical volume','Convergence, infinite volume',...
    'Simulation result','Location','SouthEast')

%saveas(gcf,'./plots/GG_Integrals.pdf')

%%

