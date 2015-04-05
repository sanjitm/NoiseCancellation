function seism = addWavelets(seism,time,ff,dd,cc,nn,pp)

i0 = length(seism.waves);

for k = 1:length(ff)
    seism.waves(i0+k).type = 'wavelet';
    seism.waves(i0+k).freq = ff(k);
    seism.waves(i0+k).dT = dd(k);
    seism.waves(i0+k).speed = cc(k);
    seism.waves(i0+k).num = nn(k);
    seism.waves(i0+k).perc = pp(k);
    
    %the 10 compensates to get more accurate average amplitudes
    A = 10*seism.gnd_amp*sqrt(pp(k))/sqrt(nn(k));

    %phi = -pi/2+pi*rand(1,nn(k)); %forward propagating
    %phi = pi/2+pi*rand(1,nn(k)); %backward propagating
    phi = 2*pi*rand(1,nn(k)); %isotropic field
    R0 = cc(k)*time.T*rand(1,nn(k)); %initial distance

    seism.waves(i0+k).xi0 = A*randn(1,nn(k)); %displacement amplitude
    seism.waves(i0+k).p0 = 2*pi*rand(1,nn(k)); %initial phase
    seism.waves(i0+k).k = 2*pi*ff(k)/cc(k)*[cos(phi); sin(phi)]; %wave vector
    seism.waves(i0+k).loc0 = -[cos(phi); sin(phi)].*R0(ones(2,1),:);
end