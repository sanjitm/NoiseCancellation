function xi = displaceLocalized(model,r)

ng = model.ng;
gridsize = ng(1)*ng(2)*ng(3);
xi = zeros(gridsize,3);
n_SC_P = model.n_SC_P;
n_SC_S = model.n_SC_S;

if n_SC_P>0
    lambda_P = model.lambda_P;

    loc_SC_P = model.loc_SC_P;

    p0 = 2*pi*rand(n_SC_P,1); %initial phases at scattering center

    for k = 1:n_SC_P
        rm_SC = sqrt(sum((r-r(ones(gridsize,1)*loc_SC_P(k),:)).^2,2));
        k_SC = (r-r(ones(gridsize,1)*loc_SC_P(k),:))./rm_SC(:,ones(1,3));

        phase = p0(k)+2*pi*rm_SC/lambda_P;

        xi0 = cos(phase)./rm_SC;
        xi0 = k_SC.*xi0(:,ones(1,3));
        xi0(isnan(xi0)) = 0;
        sc0i = find(rm_SC<model.localMin);
        xi0(sc0i,:) = 0;

        xi = xi+xi0;
    end
end

if n_SC_S>0
    lambda_S = model.lambda_S;

    loc_SC_S = model.loc_SC_S;

    p0 = 2*pi*rand(n_SC_S,1); %initial phases at scattering center
    psi = 2*pi*rand(n_SC_S,1); %polarization angle of shear waves

    for k = 1:n_SC_S
        rm_SC = sqrt(sum((r-r(ones(gridsize,1)*loc_SC_S(k),:)).^2,2));
        k_SC = (r-r(ones(gridsize,1)*loc_SC_S(k),:))./rm_SC(:,ones(1,3));

        e1(:,1) = -k_SC(:,2);
        e1(:,2) = k_SC(:,1);
        e1(:,3) = zeros(gridsize,1);
        e2 = cross(k_SC,e1,2);
        clear k_SC
        pol_S = cos(psi(k))*e1+sin(psi(k))*e2;
        clear e1 e2
        
        phase = p0(k)+2*pi*rm_SC/lambda_S;

        xi0 = cos(phase)./rm_SC;
        xi0 = pol_S.*xi0(:,ones(1,3));
        xi0(isnan(xi0)) = 0;

        xi = xi+xi0;
    end
end

if (n_SC_P+n_SC_S)>0
    std_xi = std(xi,0,1);
    xi = xi./std_xi(ones(gridsize,1),:)*model.sn_scatter;
end

rm = sqrt(sum((r-model.r0(ones(prod(model.ng),1),:)).^2,2)); 
xi(rm<model.cavity,:) = 0;