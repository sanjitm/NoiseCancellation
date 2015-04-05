function gnd = displace2D(gnd,t,X,Y,waves)

WN = length(waves);

for wi = 1:WN
    f = waves(wi).freq;
    if strcmp(waves(wi).type,'plane')
        for wwi = 1:waves(wi).num
            xi0 = waves(wi).xi0(wwi);
            k = waves(wi).k(:,wwi);
            p0 = waves(wi).p0(wwi);

            gnd = gnd+xi0*cos(2*pi*f*t-k(1)*X-k(2)*Y+p0);
        end
    end
    if strcmp(waves(wi).type,'spherical')
        k0 = waves(wi).k0;
        for wwi = 1:waves(wi).num
            R = sqrt((X-waves(wi).loc(1,wwi)).^2+(Y-waves(wi).loc(2,wwi)).^2);
            xi0 = waves(wi).xi0(wwi);
            p0 = waves(wi).p0(wwi);

            gnd = gnd+xi0*besselj(0,k0*R)*cos(2*pi*f*t+p0);
        end
    end
    if strcmp(waves(wi).type,'wavelet')
        dT = waves(wi).dT;
        for wwi = 1:waves(wi).num
            xi0 = waves(wi).xi0(wwi);
            p0 = waves(wi).p0(wwi);
            k = waves(wi).k(:,wwi);
            X0 = waves(wi).loc0(1,wwi);
            Y0 = waves(wi).loc0(2,wwi);
            
            tau = t-(k(1)*(X-X0)+k(2)*(Y-Y0))/(2*pi*f);

            gnd = gnd+xi0*exp(-tau.^2/(2*dT)^2)...
                .*cos(2*pi*f*tau+p0);
        end
    end
end