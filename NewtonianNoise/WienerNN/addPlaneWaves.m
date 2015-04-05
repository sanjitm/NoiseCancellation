function seism = addPlaneWaves(seism,ff,ss,nn,rr)

i0 = length(seism.waves);

for k = 1:length(ff)
    seism.waves(i0+k).type = 'plane';
    seism.waves(i0+k).freq = ff(k);
    seism.waves(i0+k).speed = ss(k);
    seism.waves(i0+k).num = nn(k);
    seism.waves(i0+k).perc = rr(k);
    
    A = seism.gnd_amp*sqrt(rr(k))/sqrt(nn(k));
       
    %phi = -pi/2+pi*rand(1,nn(k)); %forward propagating
    %phi = pi/2+pi*rand(1,nn(k)); %backward propagating
    phi = 2*pi*rand(1,nn(k)); %isotropic field

    seism.waves(i0+k).xi0 = A*randn(1,nn(k)); %displacement amplitude
    seism.waves(i0+k).p0 = 2*pi*rand(1,nn(k)); %initial phase
    seism.waves(i0+k).k = 2*pi*ff(k)/ss(k)*[cos(phi); sin(phi)]; %wave vector
end
