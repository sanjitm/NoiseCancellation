%%
%(surface, 1000m P wavelength, seismic=1)
%GGN from Rayleigh waves: 2.2e-12
%(1250m depth, 1000m P wavelength, seismic=1)
%GGN from Rayleigh waves: 3e-15
%(2500m depth, 1000m P wavelength, seismic=1)
%GGN from Rayleigh waves: 4.6e-16

%GGN from body waves: 1.4e-12
%%
clear all
tic
%--------------------------------------------------------
% seismic waves
%--------------------------------------------------------
model.freq = 3; %frequency of seismic waves
model.sn_iso = 1;%5e-10/model.freq^1.7; %seismic displacement [m/rHz]
model.sn_scatter = 1;%1e-10/model.freq^1.7; 
model.sn_surface = 1;%3e-9/model.freq^1.7;
model.c_P = 3000; %speed of P waves [m/s]
model.Poisson = 0.25; %Poisson ratio usually between 0.2 and 0.3
model.rho0 = 2500; %mean rock density [kg/m^3]
model.cavity = 2; %radius of excision of rock body around test mass [m]
model.n_P = 50; %number of P waves per field
model.n_S = 0;%100; %number of S waves per field
model.n_R = 0;%10; %number of Rayleigh waves per field
model.n_SC = 0; %number of localized sources, i.e. spherical waves
%--------------------------------------------------------
% grid
%--------------------------------------------------------
%points at which to calculate GG
r0 = [0 0 0];
%r0 = [zeros(40,1) zeros(40,1) linspace(1000,-400,40)'];
%model.radius = model.c_P/model.freq*ones(1,3); %radius of FEM grid [m]
model.radius = 2*[1000 1000 1000]; %radius of FEM grid [m]
%2D cut through grid; one non-zero number only; value defines cut location
model.slice_2D = [0 0 0.5]; 
model.coh_level = [0.5];
model.surface = 2*1000;%model.radius(3); %surface level [m]
model.ng = [100 100 round(100*(model.surface+model.radius(1))/(2*model.radius(1)))];
%--------------------------------------------------------
% auxiliary
%--------------------------------------------------------
model.samples = 1; %number of fields calculated to derive standard deviations etc

model = addWaveModes(model);
model.isoField = [];

for k = 1:length(r0(:,1))
    plots.plcnt = 1;
    model.r0 = r0(k,:);
    [output, model, plots] = gravityGradient(model,plots);
    if model.samples > 1
        stddevs(k,:) = model.a_std;
    end
end

%model.a_std/(2*pi*model.freq)^2/2000
disp(['Calculation time of simulation: ' num2str(toc)])
%%
% figure(plots.plcnt)
% plots.plcnt = plots.plcnt+1;
% clf
% set(gcf, 'PaperSize',[10 8])
% set(gcf, 'PaperPosition', [0 0 10 8])
% h=semilogy(-(r0(:,3)-max(r0(:,3))),stddevs,'LineWidth',2);
% set(h,{'LineStyle'},{'-';':';'--'})
% xlabel('Depth')
% ylabel('GGN [(m/s^2)/\surd Hz]')
% grid
% axis tight
% legend('\deltag_x','\deltag_y','\deltag_z') 
% set(gca,'xtick',[0 250 500 750 1000 1250])
% set(gca,'xticklabel',{' ',' ',' '})
% text(-10,7e-10,'0')
% text(450,7e-10,'\lambda_P/2')
% text(980,7e-10,'\lambda_P')
% saveas(gcf,'C:/MyStuff/publish/paperGGN/DepthGG.eps')
% saveas(gcf,'./plots/DepthGG.pdf')
%%
% set(gca,'xtick',[-2000 -1000 0 1000 2000])
% set(gca,'xticklabel',{' ',' ',' '})
% text(-2100,-1.1,'-2\lambda_P')
% text(-1060,-1.1,'-\lambda_P')
% text(-40,-1.1,'0')
% text(940,-1.1,'\lambda_P')
% text(1900,-1.1,'2\lambda_P')
% saveas(gcf,'./plots/Coh_zz.pdf')
% saveas(gcf,'C:/MyStuff/publish/paperGGN/Coh_zz.eps')
