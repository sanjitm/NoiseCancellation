function xi = displaceIsotropically(model,r)

%z_sur = model.surface;
ng = model.ng;
gridsize = ng(1)*ng(2)*ng(3);
xi = zeros(gridsize,3);
n_P = model.n_P;
n_S = model.n_S;

if n_P>0
    lambda_P = model.lambda_P;

    phi_P = 2*pi*rand(n_P,1); %uniform distribution of phi
    cost_P = 2*rand(n_P,1)-1; %uniform distribution of cos(theta) with theta in [0,pi]
    %cost_P = rand(n_P,1); %uniform distribution of cos(theta) with theta in [0,pi/2]

    k_P = [sqrt(1-cost_P.^2).*sin(phi_P) sqrt(1-cost_P.^2).*cos(phi_P) cost_P];
    %k_Pr = [sqrt(1-cost_P.^2).*sin(phi_P) sqrt(1-cost_P.^2).*cos(phi_P) -cost_P];

    p0_P = 2*pi*rand(n_P,1); %initial phases
    %p0_Pr = p0_P+2*cost_P*z_sur+pi; %phases of reflected waves
    
    for k = 1:n_P       
        kk = k_P(ones(gridsize,1)*k,:);
        phase = p0_P(k)+2*pi/lambda_P*sum(kk.*r,2);
        xi = xi+model.sn_iso*kk.*cos(phase(:,ones(1,3)));
        
%         %reflected wave P->P 
%         kk = ones(gridsize,1)*k_Pr(k,:);
%         phase = p0_Pr(k)+2*pi/lambda_P*sum(kk.*r,2);
%         xi = xi+kk.*cos(phase(:,ones(1,3)));
    end
end

if n_S>0
    lambda_S = model.lambda_S;

    phi_S = 2*pi*rand(n_S,1); %uniform distribution of phi
    cost_S = 2*rand(n_S,1)-1; %uniform distribution of cos(theta) with theta in [0,pi]
    %cost_S = rand(n_S,1); %uniform distribution of cos(theta) with theta in [0,pi/2]

    psi_S = 2*pi*rand(n_S,1); %polarization angle of shear waves

    e1 = [cost_S.*sin(phi_S) cost_S.*cos(phi_S) -sqrt(1-cost_S.^2)];
    %e1r = [-cost_S.*sin(phi_S) -cost_S.*cos(phi_S) -sqrt(1-cost_S.^2)];
    e2 = [cos(phi_S) -sin(phi_S) zeros(n_S,1)];

    pol_S = cos(psi_S(:,ones(1,3))).*e1+sin(psi_S(:,ones(1,3))).*e2
    %pol_Sr = cos(psi_S(:,ones(1,3))).*e1r+sin(psi_S(:,ones(1,3))).*e2;

    k_S = [sqrt(1-cost_S.^2).*sin(phi_S) sqrt(1-cost_S.^2).*cos(phi_S) cost_S]
    %k_Sr = [sqrt(1-cost_S.^2).*sin(phi_S) sqrt(1-cost_S.^2).*cos(phi_S) -cost_S];

    p0_S = 2*pi*rand(n_S,1);
    %p0_Sr = p0_S+2*cost_S*z_sur;

    for k = 1:n_S
        kk = ones(gridsize,1)*k_S(k,:);
        phase = p0_S(k)+2*pi/lambda_S*sum(kk.*r,2);
        xi = xi+model.sn_iso*pol_S(ones(gridsize,1)*k,:).*cos(phase(:,ones(1,3)));
        
%         %reflected wave SV->SV
%         kk = ones(gridsize,1)*k_Sr(k,:);
%         phase = p0_Sr(k)+2*pi/lambda_S*sum(kk.*r,2);
%         xi = xi+pol_Sr(ones(gridsize,1)*k,:).*cos(phase(:,ones(1,3)));
    end
end

% if n_P+n_S>0
%     std_xi = std(xi,0,1);
%     xi = xi./std_xi(ones(gridsize,1),:)*model.sn_iso;
% end

rm = sqrt(sum((r-model.r0(ones(prod(model.ng),1),:)).^2,2)); 
xi(rm<model.cavity,:) = 0;

% test = reshape(xi,ng(1),ng(2),ng(3),3);
% std(squeeze(std(squeeze(test(:,:,end,:)),0,1)),0,1)
% std(squeeze(std(squeeze(test(:,:,1,:)),0,1)),0,1)