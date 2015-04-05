function seism = addSphericalWaves(seism,ff,cc,rr,nn,pp)

i0 = length(seism.waves);

for k = 1:length(ff)
    seism.waves(i0+k).type = 'spherical';
    seism.waves(i0+k).freq = ff(k);
    seism.waves(i0+k).speed = cc(k);
    seism.waves(i0+k).num = nn(k);
    seism.waves(i0+k).perc = pp(k);
    seism.waves(i0+k).k0 = 2*pi*ff(k)/cc(k);
    
    A = 10*seism.gnd_amp*sqrt(pp(k))/sqrt(nn(k));

    seism.waves(i0+k).xi0 = A*randn(1,nn(k)); %displacement amplitude
    seism.waves(i0+k).p0 = 2*pi*rand(1,nn(k)); %initial phase
    phi = 2*pi*rand(1,nn(k));
    seism.waves(i0+k).loc = rr(k)*[cos(phi); sin(phi)];
    seism.waves(i0+k).az = phi;
end
