function seism = addSpiralSeism(seism,A,SN)

phi = linspace(0,2*pi*A,SN)';
rho = linspace(0,seism.L,SN)';

seism.loc = [cos(phi) sin(phi)].*rho(:,ones(1,2));