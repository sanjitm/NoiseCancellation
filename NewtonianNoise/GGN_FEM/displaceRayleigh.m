function xi = displaceRayleigh(model,r)

z_sur = model.surface; %z-coordinate of surface
ng = model.ng; %number of grid points
xi = zeros(ng(1)*ng(2)*ng(3),3); %displacement field
xi_sur = zeros(ng(1)*ng(2),3); %surface displacement only required for normalization
n_R = model.n_R; %number of plane Rayleigh waves
if n_R > 0
    l_R = model.lambda_R; %Rayleigh wave length determined previously
    l_S = model.lambda_S; %shear wave length
    l_P = model.lambda_P; %compressional wave length
    
    kh = 2*pi/l_R; %horizontal wave number
    kv_S = 2*pi*sqrt(1/l_R^2-1/l_S^2); %vertical shear-wave number
    kv_P = 2*pi*sqrt(1/l_R^2-1/l_P^2); %vertical compressional-wave number
    
    zeta = sqrt(kv_P/kv_S); %parameter used in Rayleigh-wave field

    phi = 2*pi*rand(n_R,1); %uniform distribution of phi
    p0 = 2*pi*rand(n_R,1); %initial phases

    k_R = [cos(phi) sin(phi) zeros(n_R,1)]; %wave vector
    
    %function argument r comes as 2D vector is converted into 4D:
    % 3D for location + 1D for components of displacement vector
    r = reshape(r,ng(1),ng(2),ng(3),3); 
    
    r_sur = squeeze(r(:,:,1,:)); %coordinates of surface grid
    r_sur(:,:,3) = z_sur;
    
    %now r can be converted back into the more suitable 2D form
    r = reshape(r,ng(1)*ng(2)*ng(3),3);
    r_sur = reshape(r_sur,ng(1)*ng(2),3);
    
    %vector of distances of grid points to the surface
    rm_R = abs(r(:,3)-z_sur(ones(ng(1)*ng(2)*ng(3),1)));

    for k = 1:n_R     
        kk = k_R(ones(ng(1)*ng(2)*ng(3),1)*k,:); %assign wavevector to each grid point

        phase = p0(k)+sum(kh*kk.*r,2); %phase of the wave at each grid point

        %phase of the wave at the surface
        phase_sur = p0(k)+sum(kh*k_R(ones(ng(1)*ng(2),1)*k,:).*r_sur,2);

        %horizontal and vertical displacement at each grid point
        xi_h = (kh*exp(-kv_P*rm_R)-zeta*kv_S*exp(-kv_S*rm_R)).*sin(phase);
        xi_v = -(kv_P*exp(-kv_P*rm_R)-zeta*kh*exp(-kv_S*rm_R)).*cos(phase);
        
        %horizontal and vertical surface displacement at each grid point
        xi_sur_h = (kh-zeta*kv_S)*sin(phase_sur);
        xi_sur_v = -(kv_P-zeta*kh)*cos(phase_sur);
        
        %sum contributions from all n_R Rayleigh waves
        xi = xi+xi_h*k_R(k,:)+xi_v*[0 0 1];
        xi_sur = xi_sur+xi_sur_h*k_R(k,:)+xi_sur_v*[0 0 1];
    end

    std_xi = std(xi_sur,0,1); %calculate normalization factor
    %normalize dispalcement field using model amplitude spectrum
    xi = xi./std_xi(ones(ng(1)*ng(2)*ng(3),1),:)*model.sn_surface;
end

%determine distance of grid points to the test mass at r0
rm = sqrt(sum((r-model.r0(ones(prod(model.ng),1),:)).^2,2)); 
%set displacement to zero at grid points that would lie within the cavity
%around the test mass
xi(rm<model.cavity,:) = 0;