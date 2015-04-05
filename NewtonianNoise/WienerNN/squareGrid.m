function sgrid = squareGrid2D(sgrid)

L = sgrid.L;
nL = sgrid.nL;

sgrid.dS = (2*L/(nL-1))^2;
sgrid.xx = linspace(-L,L,nL);
sgrid.yy = linspace(-L,L,nL);
[sgrid.X,sgrid.Y] = meshgrid(sgrid.xx,sgrid.yy);
sgrid.R = sqrt(sgrid.X.^2+sgrid.Y.^2+sgrid.hTM^2);
sgrid.cos = sgrid.X./sgrid.R;