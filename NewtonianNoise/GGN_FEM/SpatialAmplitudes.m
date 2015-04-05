set(0,'DefaultAxesFontSize',22);
set(0,'DefaultTextFontSize',22);

k = linspace(0,3,200);
k0=1;
r0=1;

da = -2*pi*i./k.*(2*k0./(k.^2-k0^2)+1./k.*log(abs((k+k0)./(k-k0))))...
    .*(cos(k*r0)./(k*r0)-sin(k*r0)./(k*r0).^2);

figure(1)
set(gcf, 'PaperSize',[8 5])
set(gcf, 'PaperPosition', [0 0 8 5])
plot(k,abs(da),'LineWidth',2)
grid
set(gca,'ylim',[0 350])
xlabel('Wave Number [k/k_0]')
ylabel('Spatial Amplitude Spectrum')

%saveas(gcf,'./plots/SpatialSpectrum.eps')
