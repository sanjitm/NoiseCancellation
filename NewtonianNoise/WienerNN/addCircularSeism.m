function seism = addCircularSeism(seism,SN)

phi = linspace(0,2*pi,SN+1)';
phi(end) = [];

seism.loc = [seism.loc; seism.L*[cos(phi) sin(phi)]];