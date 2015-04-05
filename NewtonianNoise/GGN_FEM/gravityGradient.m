function [data, model, plots] = gravityGradient(model,plots)

set(0,'DefaultAxesFontSize',24);
set(0,'DefaultTextFontSize',24);

G = 6.673e-11; %Newton's constant

r0 = model.r0;
cavity = model.cavity;
radius = model.radius;

[r, model] = makeGrid(model);

ng = model.ng;

x = squeeze(r(1,:,1,1))';
y = squeeze(r(:,1,1,2));
z = squeeze(r(1,1,:,3));

r = reshape(r,prod(ng),3);

dV = prod(max(r,[],1)-min(r,[],1))/prod(ng);
rm = sqrt(sum((r-r0(ones(prod(ng),1),:)).^2,2));

cohi = floor(model.coh_level*ng(3));

line = zeros(length(cohi),model.samples,ng(1),3);
a_gradient = zeros(model.samples,3);

%scattering centers fixed during simulation
model.loc_SC_P = max(floor(ng(1)*ng(2)*ng(3)*rand(model.n_SC_P,1)),1);
model.loc_SC_S = max(floor(ng(1)*ng(2)*ng(3)*rand(model.n_SC_S,1)),1);

[min_rm, r0i] = min(rm);

for s = 1:model.samples
    
    if isempty(model.displacementField)
        xi = displaceIsotropically(model,r);
    else
        xi = model.displacementField;
    end
    
    a_field = G*dV*model.rho0./rm(:,ones(1,3)).^3.*(...
        -3*(sum((r-r0(ones(prod(ng),1),:)).*xi,2)*ones(1,3))...
        .*(r-r0(ones(prod(ng),1),:))./rm(:,ones(1,3)).^2+xi);
    
    a_field(rm<model.cavity,:) = 0;
    a_gradient(s,:) = sum(a_field,1);

    data.disp0 = (xi(r0i+1,:)+xi(r0i-1,:))/2;
    
    xi = reshape(xi,ng(1),ng(2),ng(3),3);
    
    for k = 1:length(cohi)
        line(k,s,:,:) = squeeze(xi(:,floor(ng(2)/2),cohi(k),:));
    end
end


model.loc_SC_P = r(model.loc_SC_P,:);
model.loc_SC_S = r(model.loc_SC_S,:);

if model.samples > 1
    curves = {};
    for k = 1:length(cohi)  
        points_x = squeeze(line(k,:,:,1));
        point_x = squeeze(line(k,:,floor(ng(1)/2),1))'*ones(1,ng(1));
        points_z = squeeze(line(k,:,:,3));
        point_z = squeeze(line(k,:,floor(ng(1)/2),3))'*ones(1,ng(1));

        stds_x = sqrt(mean(points_x.^2,1));
        std_x = sqrt(mean(point_x.^2,1));

        stds_z = sqrt(mean(points_z.^2,1));
        std_z = sqrt(mean(point_z.^2,1));

        corr_xx = mean(points_x.*point_x,1);
        corr_zz = mean(points_z.*point_z,1);

        model.coh_xx(k,:) = corr_xx./(stds_x.*std_x);
        model.coh_zz(k,:) = corr_zz./(stds_z.*std_z);

        curves = {curves{:}, [num2str(10*round(z(cohi(k))/10)) 'm']};
    end
    model.a_std = sqrt(squeeze(mean(a_gradient.^2,1)));
end

data.conv_GG = zeros(floor(max(ng)/2),3);
data.contr_GG = zeros(floor(max(ng)/2),3);
data.dist_GG = linspace(0,max(radius),floor(max(ng)/2));

for k = 2:floor(ng/2)
    ai = (rm>cavity).*(rm>data.dist_GG(k-1)).*(rm<=data.dist_GG(k));
    data.contr_GG(k,:) = sum(a_field.*ai(:,ones(1,3)),1);
end
data.conv_GG = cumsum(data.contr_GG,1);
data.beyond_GG = cumsum(data.contr_GG(end:-1:1,:),1);
n = length(data.conv_GG(:,1));

xi = reshape(xi,ng(1),ng(2),ng(3),3);
r = reshape(r,ng(1),ng(2),ng(3),3);

n_cuts = length(model.slice_2D(:,1));
n1 = 0;
n2 = 0;
n3 = 0;
for k = 1:n_cuts
    model.slice_2D(k,:);
    dim = find(model.slice_2D(k,:));
    if dim == 1
        n1 = n1+1;
        fac = floor(model.slice_2D(k,dim)*ng(dim));
        ivec1 = 1:3:ng(2); %pick quiver indices
        ivec2 = 1:3:ng(3); %pick quiver indices
        data(k).coor_2D = squeeze(r(:,fac,:,2:3));
        data(k).xi_2D = squeeze(xi(:,fac,:,2:3));
        model.cut_x(n1) = x(fac);
        lab1 = 'North - South [m]';
        lab2 = 'Up - Down [m]';
    elseif dim == 2
        n2 = n2+1;
        fac = floor(model.slice_2D(k,dim)*ng(dim));
        ivec1 = 1:3:ng(1); %pick quiver indices
        ivec2 = 1:3:ng(3); %pick quiver indices
        data(k).coor_2D = squeeze(r(fac,:,:,[1 3]));
        data(k).xi_2D = squeeze(xi(fac,:,:,[1 3]));
        model.cut_y(n2) = y(fac);
        lab1 = 'East - West [m]';
        lab2 = 'Up - Down [m]';
    else
        n3 = n3+1;
        fac = floor(model.slice_2D(k,dim)*ng(dim));
        ivec1 = 1:3:ng(1); %pick quiver indices
        ivec2 = 1:3:ng(2); %pick quiver indices
        data(k).coor_2D = squeeze(r(:,:,fac,1:2));
        data(k).xi_2D = squeeze(xi(:,:,fac,1:2));
        model.cut_z(n3) = z(fac);
        lab1 = 'East - West [m]';
        lab2 = 'North - South [m]';
    end

    figure(plots.plcnt)
    plots.plcnt = plots.plcnt+1;
    set(gcf, 'PaperSize',[10 8])
    set(gcf, 'PaperPosition', [0 0 10 8])
    quiver(data(k).coor_2D(ivec1,ivec2,1),data(k).coor_2D(ivec1,ivec2,2),...
        data(k).xi_2D(ivec1,ivec2,1),data(k).xi_2D(ivec1,ivec2,2),'k')
    axis tight
    set(gca,'DataAspectRatio',[1 1 1])
    xlabel(lab1)
    ylabel(lab2)
    grid
    
    figure(plots.plcnt)
    plots.plcnt = plots.plcnt+1;
    set(gcf, 'PaperSize',[12 8])
    set(gcf, 'PaperPosition', [0 0 10 8])
    drho = -model.rho0*divergence(squeeze(r(:,:,:,1)),squeeze(r(:,:,:,2)),squeeze(r(:,:,:,3)),...
        squeeze(xi(:,:,:,1)),squeeze(xi(:,:,:,2)),squeeze(xi(:,:,:,3)));
    contourf(data(k).coor_2D(:,:,1),data(k).coor_2D(:,:,2),drho(:,:,fac),50)
    shading flat
    xlabel(lab1)
    ylabel(lab2)
    grid
    t = colorbar('peer',gca);
    set(get(t,'ylabel'),'String','Density change [kg/m^3]')
end

if model.samples > 1
    figure(plots.plcnt)
    plots.plcnt = plots.plcnt+1;
    set(gcf, 'PaperSize',[10 8])
    set(gcf, 'PaperPosition', [0 0 10 8])
    plot(linspace(-radius(1),radius(1),ng(1)),model.coh_xx,'LineWidth',2)
    axis tight
    set(gca,'ylim',[-1 1])
    grid
    xlabel('Distance [m]')
    ylabel('Coherence \langle\xi_x,\xi_x\rangle')
    legend(curves)

    figure(plots.plcnt)
    plots.plcnt = plots.plcnt+1;
    set(gcf, 'PaperSize',[10 8])
    set(gcf, 'PaperPosition', [0 0 10 8])
    plot(linspace(-radius(1),radius(1),ng(1)),model.coh_zz,'LineWidth',2)
    axis tight
    set(gca,'ylim',[-1 1])
    grid
    xlabel('Distance [m]')
    ylabel('Coherence \langle\xi_z,\xi_z\rangle')
    legend(curves)
end

figure(plots.plcnt)
plots.plcnt = plots.plcnt+1;
set(gcf, 'PaperSize',[8 6])
set(gcf, 'PaperPosition', [0 0 8 6])
tmp = 100*data.conv_GG./data.conv_GG(n*ones(n,1),:);
h = plot(data.dist_GG(2:end),tmp(2:end,:),'LineWidth',2);
set(h,{'LineStyle'},{'-';':';'--'})
axis tight
xlabel('Distance [m]')
ylabel('Gravity Gradient [%]')
legend('\deltag_x','\deltag_y','\deltag_z')
grid

figure(plots.plcnt)
plots.plcnt = plots.plcnt+1;
set(gcf, 'PaperSize',[8 6])
set(gcf, 'PaperPosition', [0 0 8 6])
tmp = 100*data.beyond_GG(end:-1:1,:)./data.conv_GG(n*ones(n,1),:);
h = plot(data.dist_GG(2:end),tmp(2:end,:),'LineWidth',2);
set(h,{'LineStyle'},{'-';':';'--'})
axis tight
xlabel('Distance [m]')
ylabel('Gravity Gradient [%]')
legend('\deltag_x','\deltag_y','\deltag_z')
grid
