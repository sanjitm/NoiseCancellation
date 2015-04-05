%%
%(surface, 1000m P wavelength, seismic=1)
%GGN from Rayleigh waves: 2.2e-12
%(1250m depth, 1000m P wavelength, seismic=1)
%GGN from Rayleigh waves: 3e-15
%(2500m depth, 1000m P wavelength, seismic=1)
%GGN from Rayleigh waves: 4.6e-16

%GGN from body waves: 1.4e-12
%%
set(0,'DefaultAxesFontSize',22);
set(0,'DefaultTextFontSize',22);

clear all
close all
%%
tic
%--------------------------------------------------------
% seismic waves
%--------------------------------------------------------
model.freq = 10; %frequency of seismic waves
model.sn_iso = 5e-10/model.freq^1.7; %seismic displacement [m/rHz]
model.sn_scatter = 1; %1e-10/model.freq^1.7; 
model.sn_surface = 1; %3e-9/model.freq^1.7;
model.c_P = 3000; %speed of P waves [m/s]
model.Poisson = 0.25; %Poisson ratio usually between 0.2 and 0.3
model.rho0 = 2500; %mean rock density [kg/m^3]
model.cavity = 2; %radius of excision of rock body around test mass [m]
model.localMin = 40; %size of local source (minimum radius) [m]
model.n_P = 5; %number of P waves per field
model.n_S = 10; %number of S waves per field
model.n_R = 0; %10; %number of Rayleigh waves per field
model.n_SC_P = 0;%3087; %number of spherical P waves
model.n_SC_S = 0;%6174; %number of "spherical" S waves
%--------------------------------------------------------
% grid
%--------------------------------------------------------
%points at which to calculate GG
model.r0 = [0 0 0];
model.radius = 5*[100 100 100]; %radius of FEM grid [m]
%2D cut through grid; one non-zero number only; value defines cut location
model.slice_2D = [0 0 0.5]; 
model.coh_level = [0.5];
model.surface = model.radius(3); %surface level [m]
model.ng = 1*[100 100 round(100*(model.surface+model.radius(3))/(2*model.radius(3)))];
%--------------------------------------------------------
% auxiliary
%--------------------------------------------------------
model.samples = 1; %number of fields calculated to derive standard deviations etc

model = addWaveModes(model);
model.displacementField = [];

km = 1;
limit = zeros(km,1);
centSph = zeros(km,3);
for k = 1:km
    plots.plcnt = 1;
    
    [output, model, plots] = gravityGradient(model,plots);
    if model.samples > 1
        stddevs(k,:) = model.a_std;
    end
    
    nr = 8*pi*6.673e-11*model.rho0*output.disp0;
    
    limit(k) = output.conv_GG(end,1)/nr(1);
    if model.n_SC_P+model.n_SC_S>0
        centSph(k,:) = model.loc_SC_P(1,:);
    end
end
%%
figure(1),saveas(gcf,'./plots/Displacement.png')
figure(2),saveas(gcf,'./plots/Density.png')
%%
plots.plcnt = 4;
figure(plots.plcnt)
plots.plcnt = plots.plcnt+1;
clf
set(gcf, 'PaperSize',[10 8])
set(gcf, 'PaperPosition', [0 0 10 8])
n = length(output.conv_GG(:,1));
avc = NaN*ones(n-1,n);
for k1 = 1:(n-1)
    for k2 = 1:(n-k1)
        avc(k1,k2) = abs(mean(output.conv_GG(k1:(k1+k2),1))/output.conv_GG(end,1)-1);
    end
end

contourf(output.dist_GG/1000,output.dist_GG(2:end)/1000,log10(avc),50);
shading flat
xlabel('Start Distance [ \lambda_P]')
ylabel('Averaging Length [ \lambda_P]')
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String','Relative Error [log]')
axis([0 2 0.4 1])
set(gca,'clim',[-5 -1])
saveas(gcf,'C:\MyStuff\publish\paperGGN\Averaging.eps','epsc')
%%
figure(4)
set(gcf, 'PaperSize',[8 5])
set(gcf, 'PaperPosition', [0 0 8 5])
dist = sqrt(sum(centSph.^2,2));
plot(dist/1000,(limit*3-1)./dist*min(dist),'ko','MarkerFaceColor',[0 0 0])
grid
set(gca,'ylim',[-0.2 0.2])
xlabel('Distance [\lambda]')
ylabel('Relative Deviation')
%saveas(gcf,'./plots/SphericalTest.pdf')
%%
av = mean(output.conv_GG,1);
nfft = length(output.conv_GG(:,1));
bins = floor(nfft/2);
dl = 2*model.radius./model.ng;
spatial = fft(output.conv_GG-av(ones(nfft,1),:),[],1)*dl(1);
spatial = abs(spatial(1:bins,:)).^2./(2*model.radius(ones(bins,1),:));

figure(plots.plcnt)
plots.plcnt = plots.plcnt+1;
k_val = 2*pi*linspace(0,1/(2*dl(1)),bins);
h = semilogx(k_val,spatial,'LineWidth',2);
set(h,{'LineStyle'},{'-';':';'--'})
grid
axis tight
legend('\deltag_x','\deltag_y','\deltag_z')
xlabel('Wavenumber [1/m]')
ylabel('Spatial Spectrum')

figure(plots.plcnt)
plots.plcnt = plots.plcnt+1;
lambda = 2*pi./k_val;
h = semilogx(lambda(end:-1:2),spatial(end:-1:2,:),'LineWidth',2);
set(h,{'LineStyle'},{'-';':';'--'})
grid
axis tight
legend('\deltag_x','\deltag_y','\deltag_z','Location','NorthWest')
xlabel('Wavelength [m]')
ylabel('Spatial Spectrum')


mean(output.conv_GG(logical((output.dist_GG>1000).*(output.dist_GG<1500)),:)-100,1)
mean(output.conv_GG(logical((output.dist_GG>1000).*(output.dist_GG<2000)),:)-100,1)

disp(['Calculation time of simulation: ' num2str(toc)])
%%
figure(2)
y_pos = -8;
xlabel('Distance')
set(gca,'xtick',[0 500 1000 1500 2000 2500 3000])
set(gca,'xticklabel',{' ',' ',' ',' ',' '})
text(0,y_pos,'0')
text(970,y_pos,'\lambda_P')
text(1930,y_pos,'2\cdot\lambda_P')
text(2930,y_pos,'3\cdot\lambda_P')
legend('Location','SouthEast')

figure(3)
y_pos = -11;
xlabel('Distance')
set(gca,'ylim',[-10 10])
set(gca,'xtick',[0 500 1000 1500 2000 2500 3000])
set(gca,'xticklabel',{' ',' ',' ',' ',' '})
text(0,y_pos,'0')
text(970,y_pos,'\lambda_P')
text(1930,y_pos,'2\cdot\lambda_P')
text(2930,y_pos,'3\cdot\lambda_P')
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
