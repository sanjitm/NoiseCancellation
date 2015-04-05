function seism = addRandSeism(seism,SN)

L = seism.L;

seism.loc = [seism.loc; -L+2*L*rand(SN,2)];