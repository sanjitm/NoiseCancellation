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
%close all
%--------------------------------------------------------
% seismic waves
%--------------------------------------------------------
model.freq = 3; %frequency of seismic waves
model.sn_iso = 1; %5e-10/model.freq^1.7; %seismic displacement [m/rHz]
model.sn_scatter = 1;%1e-10/model.freq^1.7; 
model.sn_surface = 1;%3e-9/model.freq^1.7;
model.c_P = 3000; %speed of P waves [m/s]
model.Poisson = 0.25; %Poisson ratio usually between 0.2 and 0.3
model.rho0 = 2500; %mean rock density [kg/m^3]
model.cavity = 10; %radius of excision of rock body around test mass [m]
model.n_P = 50; %number of P waves per field
model.n_S = 100; %number of S waves per field
model.n_R = 0;%10; %number of Rayleigh waves per field
model.n_SC_P = 0; %number of spherical pressure waves
model.n_SC_S = 0; %number of spherical shear waves
%--------------------------------------------------------
% grid
%--------------------------------------------------------
%points at which to calculate GG
model.r0 = [0 0 0];
model.radius = 1500*ones(1,3); %radius of FEM grid [m]
model.ng = 70*ones(1,3);
%2D cut through grid; one non-zero number only; value defines cut location
model.slice_2D = [0 0 0.5]; 
model.coh_level = [];
model.surface = 2000; %model.radius(3); %surface level [m]
%--------------------------------------------------------
% auxiliary
%--------------------------------------------------------
model.samples = 1; %number of fields calculated to derive standard deviations etc
model = addWaveModes(model);

%%
% --------- define coordinate grids and original displacement field
[r, model] = makeGrid(model);
r = reshape(r,prod(model.ng),3);

orig_iso = displaceIsotropically(model,r);
orig_ng = model.ng;
orig_radius = model.radius;

r = reshape(r,orig_ng(1),orig_ng(2),orig_ng(3),3);

model.isoField = orig_iso;
model.localField = [];
model.surfaceField = [];

%%
% ----- STEP 1: fine grid simulation (reality)
plots.plcnt = 1;
[out_1, model, plots] = gravityGradient(model,plots);
%%
% ----- STEP 2: coarse grid simulation (seismometer data)
factor = 5;
model.isoField = reshape(orig_iso,orig_ng(1),orig_ng(2),orig_ng(3),3);
model.isoField = reduceGrid(model.isoField,model,factor,'uniform');
r_ru = reduceGrid(r,model,factor,'uniform');
disp(['Number of seismometers (uniform density): ' num2str(size(r_ru(:,:,:,1)))])
model.ng = size(squeeze(r_ru(:,:,:,1)));
model.isoField = reshape(model.isoField,prod(model.ng),3);
model.radius = max(reshape(r_ru,prod(model.ng),3),[],1);

[out_2, model, plots] = gravityGradient(model,plots);
%%
% ----- STEP 3: interpolated field to see how well one can recover reality
reduced = reshape(model.isoField,model.ng(1),model.ng(2),model.ng(3),3);
first = interp3(squeeze(reduced(:,:,:,1)),2,'cubic');
model.isoField = zeros([size(first) 3]);
model.isoField(:,:,:,1) = first;
model.isoField(:,:,:,2) = interp3(squeeze(reduced(:,:,:,2)),2,'cubic');
model.isoField(:,:,:,3) = interp3(squeeze(reduced(:,:,:,3)),2,'cubic');
clear reduced
model.isoField = reshape(model.isoField,numel(first),3);
model.ng = size(first);
clear first

[out_3, model, plots] = gravityGradient(model,plots);

%%
% ----- STEP 4: irregular seismometer array (quadratic density)
factor = 5;
model.ng = orig_ng;
model.isoField = reshape(orig_iso,orig_ng(1),orig_ng(2),orig_ng(3),3);
optimized = reduceGrid(model.isoField,model,factor,'quadratic');
r_rq = reduceGrid(r,model,factor,'quadratic');
disp(['Number of seismometers (quadratic density): ' num2str(size(r_rq(:,:,:,1)))])
model.ng = size(squeeze(model.isoField(:,:,:,1)));

model.isoField = zeros([orig_ng(:)' 3]);
model.isoField(:,:,:,1) = interp3(squeeze(r_rq(:,:,:,1)),squeeze(r_rq(:,:,:,2)),squeeze(r_rq(:,:,:,3)),...
    squeeze(optimized(:,:,:,1)),...
    squeeze(r(:,:,:,1)),squeeze(r(:,:,:,2)),squeeze(r(:,:,:,3)),'spline');
model.isoField(:,:,:,2) = interp3(squeeze(r_rq(:,:,:,1)),squeeze(r_rq(:,:,:,2)),squeeze(r_rq(:,:,:,3)),...
    squeeze(optimized(:,:,:,2)),...
    squeeze(r(:,:,:,1)),squeeze(r(:,:,:,2)),squeeze(r(:,:,:,3)),'spline');
model.isoField(:,:,:,3) = interp3(squeeze(r_rq(:,:,:,1)),squeeze(r_rq(:,:,:,2)),squeeze(r_rq(:,:,:,3)),...
    squeeze(optimized(:,:,:,3)),...
    squeeze(r(:,:,:,1)),squeeze(r(:,:,:,2)),squeeze(r(:,:,:,3)),'spline');
model.ng = orig_ng;
model.isoField = reshape(model.isoField,prod(model.ng),3);

model.radius = orig_radius;
[out_4, model, plots] = gravityGradient(model,plots);

%%
% ----- STEP 5: irregular seismometer array (cubic density)
factor = 5;
model.ng = orig_ng;
model.isoField = reshape(orig_iso,orig_ng(1),orig_ng(2),orig_ng(3),3);
optimized = reduceGrid(model.isoField,model,factor,'cubic');
r_rc = reduceGrid(r,model,factor,'cubic');
disp(['Number of seismometers (cubic density): ' num2str(size(r_rc(:,:,:,1)))])
model.ng = size(squeeze(model.isoField(:,:,:,1)));

model.isoField = zeros([orig_ng(:)' 3]);
model.isoField(:,:,:,1) = interp3(squeeze(r_rc(:,:,:,1)),squeeze(r_rc(:,:,:,2)),squeeze(r_rc(:,:,:,3)),...
    squeeze(optimized(:,:,:,1)),...
    squeeze(r(:,:,:,1)),squeeze(r(:,:,:,2)),squeeze(r(:,:,:,3)),'spline');
model.isoField(:,:,:,2) = interp3(squeeze(r_rc(:,:,:,1)),squeeze(r_rc(:,:,:,2)),squeeze(r_rc(:,:,:,3)),...
    squeeze(optimized(:,:,:,2)),...
    squeeze(r(:,:,:,1)),squeeze(r(:,:,:,2)),squeeze(r(:,:,:,3)),'spline');
model.isoField(:,:,:,3) = interp3(squeeze(r_rc(:,:,:,1)),squeeze(r_rc(:,:,:,2)),squeeze(r_rc(:,:,:,3)),...
    squeeze(optimized(:,:,:,3)),...
    squeeze(r(:,:,:,1)),squeeze(r(:,:,:,2)),squeeze(r(:,:,:,3)),'spline');
model.ng = orig_ng;
model.isoField = reshape(model.isoField,prod(model.ng),3);

model.radius = orig_radius;
[out_5, model, plots] = gravityGradient(model,plots);
%%
c1 = out_1.conv_GG;
c2 = zeros(size(c1));
c2(:,1) = interp1(out_2.dist_GG,out_2.conv_GG(:,1),out_1.dist_GG,'cubic');
c2(:,2) = interp1(out_2.dist_GG,out_2.conv_GG(:,2),out_1.dist_GG,'cubic');
c2(:,3) = interp1(out_2.dist_GG,out_2.conv_GG(:,3),out_1.dist_GG,'cubic');

ratio = 100*abs((c1-c2)./c1);

figure(plots.plcnt)
plots.plcnt = plots.plcnt+1;
clf
set(gcf, 'PaperSize',[10 8])
set(gcf, 'PaperPosition', [0 0 10 8])
h = plot(out_1.dist_GG(2:end),ratio(2:end,:),'LineWidth',2);
set(h,{'LineStyle'},{'-';':';'--'})
xlabel('Distance [m]')
ylabel('Error, uniform [%]')
legend('\deltag_x','\deltag_y','\deltag_z')
grid
axis tight
set(gca,'ylim',[0 100])

%%
c1 = out_1.conv_GG;
c2 = zeros(size(c1));
c2(:,1) = interp1(out_3.dist_GG,out_3.conv_GG(:,1),out_1.dist_GG,'cubic');
c2(:,2) = interp1(out_3.dist_GG,out_3.conv_GG(:,2),out_1.dist_GG,'cubic');
c2(:,3) = interp1(out_3.dist_GG,out_3.conv_GG(:,3),out_1.dist_GG,'cubic');

ratio = 100*abs((c1-c2)./c1);

figure(plots.plcnt)
plots.plcnt = plots.plcnt+1;
clf
set(gcf, 'PaperSize',[10 8])
set(gcf, 'PaperPosition', [0 0 10 8])
h = plot(out_1.dist_GG(2:end),ratio(2:end,:),'LineWidth',2);
set(h,{'LineStyle'},{'-';':';'--'})
xlabel('Distance [m]')
ylabel('Error, uniform, interp. [%]')
legend('\deltag_x','\deltag_y','\deltag_z')
grid
axis tight
set(gca,'ylim',[0 100])

%%
c1 = out_1.conv_GG;
c2 = zeros(size(c1));
c2(:,1) = interp1(out_4.dist_GG,out_4.conv_GG(:,1),out_1.dist_GG,'cubic');
c2(:,2) = interp1(out_4.dist_GG,out_4.conv_GG(:,2),out_1.dist_GG,'cubic');
c2(:,3) = interp1(out_4.dist_GG,out_4.conv_GG(:,3),out_1.dist_GG,'cubic');

ratio = 100*abs((c1-c2)./c1);

figure(plots.plcnt)
plots.plcnt = plots.plcnt+1;
clf
set(gcf, 'PaperSize',[10 8])
set(gcf, 'PaperPosition', [0 0 10 8])
h = plot(out_1.dist_GG(2:end),ratio(2:end,:),'LineWidth',2);
set(h,{'LineStyle'},{'-';':';'--'})
xlabel('Distance [m]')
ylabel('Error, quadratic, interp. [%]')
legend('\deltag_x','\deltag_y','\deltag_z')
grid
axis tight
set(gca,'ylim',[0 100])

%%
c1 = out_1.conv_GG;
c2 = zeros(size(c1));
c2(:,1) = interp1(out_5.dist_GG,out_5.conv_GG(:,1),out_1.dist_GG,'cubic');
c2(:,2) = interp1(out_5.dist_GG,out_5.conv_GG(:,2),out_1.dist_GG,'cubic');
c2(:,3) = interp1(out_5.dist_GG,out_5.conv_GG(:,3),out_1.dist_GG,'cubic');

ratio = 100*abs((c1-c2)./c1);

figure(plots.plcnt)
plots.plcnt = plots.plcnt+1;
clf
set(gcf, 'PaperSize',[10 8])
set(gcf, 'PaperPosition', [0 0 10 8])
h = plot(out_1.dist_GG(2:end),ratio(2:end,:),'LineWidth',2);
set(h,{'LineStyle'},{'-';':';'--'})
xlabel('Distance [m]')
ylabel('Error, cubic, interp. [%]')
legend('\deltag_x','\deltag_y','\deltag_z')
grid
axis tight
set(gca,'ylim',[0 100])

%%
figure(20)
clf
set(gcf, 'PaperSize',[10 8])
set(gcf, 'PaperPosition', [0 0 10 8])
plot(1:length(r_ru(1,:,1)),r_ru(1,:,1),'v',...
    1:length(r_rq(1,:,1)),r_rq(1,:,1),'+',...
    1:length(r_rc(1,:,1)),r_rc(1,:,1),'o','MarkerSize',10)
grid
axis([0 16 -1550 1550])
xlabel('Seismometer Count (1D)')
ylabel('Distance [m]')
legend('Linear','Squared','Cubed','Location','NorthWest')
%%
figure(16)
xlabel('Distance')
set(gca,'xtick',[0 500 1000 1500])
set(gca,'xticklabel',{' ',' ',' ',' '})
text(-10,-5,'0')
text(450,-5,'\lambda_P/2')
text(980,-5,'\lambda_P')
text(1400,-5,'1.5\cdot\lambda_P')
saveas(gcf,'./plots/ErrorUniform.pdf')
saveas(gcf,'./plots/ErrorUniform.eps')
figure(17)
xlabel('Distance')
set(gca,'xtick',[0 500 1000 1500])
set(gca,'xticklabel',{' ',' ',' ',' '})
text(-10,-5,'0')
text(450,-5,'\lambda_P/2')
text(980,-5,'\lambda_P')
text(1400,-5,'1.5\cdot\lambda_P')
saveas(gcf,'./plots/ErrorUniInterp.pdf')
saveas(gcf,'./plots/ErrorUniInterp.eps')
figure(18)
xlabel('Distance')
set(gca,'xtick',[0 500 1000 1500])
set(gca,'xticklabel',{' ',' ',' ',' '})
text(-10,-5,'0')
text(450,-5,'\lambda_P/2')
text(980,-5,'\lambda_P')
text(1400,-5,'1.5\cdot\lambda_P')
saveas(gcf,'./plots/ErrorQuadInterp.pdf')
saveas(gcf,'./plots/ErrorQuadInterp.eps')
figure(19)
xlabel('Distance')
set(gca,'xtick',[0 500 1000 1500])
set(gca,'xticklabel',{' ',' ',' ',' '})
text(-10,-5,'0')
text(450,-5,'\lambda_P/2')
text(980,-5,'\lambda_P')
text(1400,-5,'1.5\cdot\lambda_P')
saveas(gcf,'./plots/ErrorCubicInterp.pdf')
saveas(gcf,'./plots/ErrorCubicInterp.eps')