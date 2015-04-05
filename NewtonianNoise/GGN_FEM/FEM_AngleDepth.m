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
model.n_P = 0; %number of P waves per field
model.n_S = 0; %number of S waves per field
model.n_R = 0; %number of Rayleigh waves per field
model.n_SC_P = 0; %number of spherical pressure waves
model.n_SC_S = 0; %number of spherical shear waves
%--------------------------------------------------------
% grid
%--------------------------------------------------------
model.radius = 3000*ones(1,3); %radius of FEM grid [m]
model.ng = [101 101 51];
%2D cut through grid; one non-zero number only; value defines cut location
model.slice_2D = [0 0 0.5]; 
model.coh_level = [];
model.surface = 0; %surface level [m]
%--------------------------------------------------------
% auxiliary
%--------------------------------------------------------
model.samples = 1; %number of fields calculated to derive standard deviations etc
model = addWaveModes(model);
nr = 8/3*pi*6.673e-11*model.rho0;

%%
% --------- define coordinate grids and original displacement field
[r, model] = makeGrid(model);
iz = 26:51;
n_depth = length(iz);
r0 = [zeros(n_depth,1) zeros(n_depth,1) squeeze(r(1,1,iz,3))];
r = reshape(r,prod(model.ng),3);

%%
n_angle = 51;
angles = linspace(0,90,n_angle);
err_rel = zeros(n_angle,n_depth,3);
err_abs = zeros(n_angle,n_depth,3);
xi_abs = zeros(n_angle,n_depth,3);
for k = 1:n_depth
    model.r0 = r0(k,:);
    disp(['Depth = ' num2str(model.r0(3)) 'm'])
    for a = 1:n_angle
        model.displacementField = injectPlaneWave(angles(a),'P',model,r);

        plots.plcnt = 1;        
        [output, model, plots] = gravityGradient(model,plots);

        err_rel(a,k,:) = (output.conv_GG(end,:)-nr*output.disp0)./output.conv_GG(end,:);
        err_abs(a,k,:) = (output.conv_GG(end,:)-nr*output.disp0)/nr;
        gg_abs(a,k,:) = output.conv_GG(end,:)/nr;
    end
end

%%
x = r0(:,3)/1000;
y = angles;
[X,Y] = meshgrid(x,y);
dim = 3;

plots.plcnt = 4;
figure(plots.plcnt)
plots.plcnt = plots.plcnt+1;
set(gcf, 'PaperSize',[10 6])
set(gcf, 'PaperPosition', [0 0 10 6])
pcolor(X,Y,log10(abs(squeeze(err_rel(:,:,dim)))))
%contourf(x,y,log10(abs(squeeze(err_rel(:,:,dim)))),20);
shading interp
xlabel('Depth [\lambda]')
ylabel('Inclination angle [deg]')
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String','log10(err_{rel}(g_z))'); 
set(gca,'clim',[-3 1])

figure(plots.plcnt)
plots.plcnt = plots.plcnt+1;
set(gcf, 'PaperSize',[10 6])
set(gcf, 'PaperPosition', [0 0 10 6])
pcolor(X,Y,log10(abs(squeeze(err_abs(:,:,dim)))))
%contourf(x,y,log10(abs(squeeze(err_abs(:,:,dim)))))
shading interp
xlabel('Depth [\lambda]')
ylabel('Inclination angle [deg]')
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String','log10(err_{abs}(g_z)/(8\pi G\rho_0\xi_0/3))');
set(gca,'clim',[-3 0])

figure(plots.plcnt)
plots.plcnt = plots.plcnt+1;
set(gcf, 'PaperSize',[10 6])
set(gcf, 'PaperPosition', [0 0 10 6])
pcolor(X,Y,log10(abs(squeeze(gg_abs(:,:,dim)))))
%contourf(x,y,log10(abs(squeeze(gg_abs(:,:,dim)))),100)
shading interp
xlabel('Depth [\lambda]')
ylabel('Inclination angle [deg]')
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String','log10(g_z/(8\pi G\rho_0\xi_0/3))'); 
set(gca,'clim',[-3 0])

%%
% figure(4),saveas(gcf,'./plots/AngleDepth_Rel_H_P.pdf')
% figure(5),saveas(gcf,'./plots/AngleDepth_Abs_H_P.pdf')
% figure(6),saveas(gcf,'./plots/AngleDepth_g_H_P.pdf')
%%
% figure(4),saveas(gcf,'./plots/AngleDepth_Rel_Z_P.pdf')
% figure(5),saveas(gcf,'./plots/AngleDepth_Abs_Z_P.pdf')
% figure(6),saveas(gcf,'./plots/AngleDepth_g_Z_P.pdf')