function xi = injectPlaneWave(angle,type,model,r)

ng = model.ng;
gridsize = ng(1)*ng(2)*ng(3);
xi = zeros(gridsize,3);

if strcmp(type,'P')
    lambda_P = model.lambda_P;

    phi_P = 2*pi*rand(1);
    theta_P = angle*pi/180;

    k_P = [sin(theta_P)*sin(phi_P) sin(theta_P)*cos(phi_P) cos(theta_P)];

    p0_P = 2*pi*rand(1); %initial phase
    
    kk = k_P(ones(gridsize,1),:);
    phase = p0_P+2*pi/lambda_P*sum(kk.*r,2);
    xi = kk.*cos(phase(:,ones(1,3)));
end

if strcmp(type,'S')
    lambda_S = model.lambda_S;

    phi_S = 2*pi*rand(1);
    theta_S = angle*pi/180;

    psi_S = 2*pi*rand(1); %polarization angle of shear wave

    e1 = [cos(theta_S)*sin(phi_S) cos(theta_S)*cos(phi_S) -sin(theta_S)];
    e2 = [cos(phi_S) -sin(phi_S) zeros(n_S,1)];

    pol_S = cos(psi_S(:,ones(1,3))).*e1+sin(psi_S(:,ones(1,3))).*e2;

    k_S = [sin(theta_S)*sin(phi_S) sin(theta_S)*cos(phi_S) cos(theta_S)];

    p0_S = 2*pi*rand(1);

    kk = k_S(ones(gridsize,1),:);
    phase = p0_S+2*pi/lambda_S*sum(kk.*r,2);
    xi = pol_S(ones(gridsize,1),:).*cos(phase(:,ones(1,3)));
end

rm = sqrt(sum((r-model.r0(ones(prod(model.ng),1),:)).^2,2)); 
xi(rm<model.cavity,:) = 0;
