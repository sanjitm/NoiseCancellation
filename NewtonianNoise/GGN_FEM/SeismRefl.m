set(0,'DefaultAxesFontSize',20);
set(0,'DefaultTextFontSize',20);

c_P = 1;
c_S = 1732/3000;

p1 = linspace(1e-3,1/c_S,1000);
p2 = linspace(1e-3,1/c_S,1000);

q1_P = sqrt(1/c_P^2-p1.^2);
q1_S = sqrt(1/c_S^2-p1.^2);
q2_P = sqrt(1/c_P^2-p2.^2);
q2_S = sqrt(1/c_S^2-p2.^2);


r_in = ((2*p1.^2-1/c_S^2).^2-4*p1.^2.*q1_P.*q1_S)...
    ./(4*p1.^2.*q1_P.*q1_S+(2*p1.^2-1/c_S^2).^2);

r_cross = -i*(4*p2.*sqrt(q2_P.*q2_S).*(2*p2.^2-1/c_S^2))...
    ./(4*p2.^2.*q2_P.*q2_S+(2*p2.^2-1/c_S^2).^2);

%%
figure(1)
clf
set(gcf, 'PaperSize',[10 8])
set(gcf, 'PaperPosition', [0 0 10 8])
plot(p1,real(r_in),p2,real(r_cross),'--','LineWidth',2)
grid
set(gca,'xlim',[0 1/c_S])
xlabel('Horizontal Slowness [1/c_P]')
ylabel('Re(Reflection Coefficient)')
legend('In mode','Cross mode','Location','SouthWest')
%saveas(gcf,'./plots/Refl_real.pdf')
%saveas(gcf,'c:/MyStuff/publish/paperGGN/Refl_real.eps')

figure(2)
clf
set(gcf, 'PaperSize',[10 8])
set(gcf, 'PaperPosition', [0 0 10 8])
plot(p1,imag(r_in),p2,imag(r_cross),'--','LineWidth',2)
grid
set(gca,'xlim',[0 1/c_S])
xlabel('Horizontal Slowness [1/c_P]')
ylabel('Im(Reflection Coefficient)')
legend('In mode','Cross mode','Location','SouthWest')
%saveas(gcf,'./plots/Refl_imag.pdf')
%saveas(gcf,'c:/MyStuff/publish/paperGGN/Refl_imag.eps')

figure(3)
clf
set(gcf, 'PaperSize',[10 8])
set(gcf, 'PaperPosition', [0 0 10 8])
[AX,H1,H2] = plotyy(p1,abs(r_in),p1,180/pi*angle(r_in));
set(get(AX(1),'Ylabel'),'String','Reflection Coefficient [abs]') 
set(AX(1),'XLim',[0 1/c_S]) 
set(AX(2),'XLim',[0 1/c_S]) 
set(AX(1),'YLim',[0 1]) 
set(AX(1),'YTick',[0 0.2 0.4 0.6 0.8 1]) 
set(AX(2),'YLim',[-180 180]) 
set(AX(2),'YTick',[-180 -90 0 90 180]) 
set(get(AX(2),'Ylabel'),'String','Reflection Coefficient [arg]') 
xlabel('Horizontal Slowness [1/c_P]')
set(H1,'LineStyle','-','LineWidth',2)
set(H2,'LineStyle','-.','LineWidth',2)
grid
%set(gca,'xlim',[0 1/c_S])
title('P-P and SV-SV Reflection')
%saveas(gcf,'./plots/Refl_IM.pdf')
%saveas(gcf,'c:/MyStuff/publish/paperGGN/Refl_IM.eps')

figure(4)
clf
set(gcf, 'PaperSize',[10 8])
set(gcf, 'PaperPosition', [0 0 10 8])
[AX,H1,H2] = plotyy(p2,abs(r_cross),p2,180/pi*angle(r_cross));
set(get(AX(1),'Ylabel'),'String','Reflection Coefficient [abs]') 
set(AX(1),'XLim',[0 1/c_S]) 
set(AX(2),'XLim',[0 1/c_S]) 
set(AX(1),'YLim',[0 sqrt(2)])
set(AX(1),'YTick',[0 0.2 0.4 0.6 0.8 1 1.2 1.4]) 
set(AX(2),'YLim',[-180 180])  
set(AX(2),'YTick',[-180 -90 0 90 180]) 
set(get(AX(2),'Ylabel'),'String','Reflection Coefficient [arg]') 
xlabel('Horizontal Slowness [1/c_P]')
set(H1,'LineStyle','-','LineWidth',2)
set(H2,'LineStyle','-.','LineWidth',2)
grid
%set(gca,'xlim',[0 1/c_S])
title('P-SV Conversion Coefficient')
%saveas(gcf,'./plots/Refl_CM.pdf')
%saveas(gcf,'c:/MyStuff/publish/paperGGN/Refl_CM.eps')