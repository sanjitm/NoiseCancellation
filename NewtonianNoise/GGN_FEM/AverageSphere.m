N=100;

x = linspace(0,2,N)';
y = linspace(-1,1,N)';

az = 180*x;
el = 180/pi*acos(y)-90;

[xx, yy] = meshgrid(x,y);
[azz, ell] = meshgrid(az,el);

freq1 = 7;
intensity1 = 0;
for k = 1:100
    angle = 2*pi*rand;
    kvec = freq1*[cos(angle) sin(angle)];
    intensity1 = intensity1 ...
        +amp*cos(2*pi*(kvec(1)*2*xx+kvec(2)*yy+rand));
end

freq2 = 3;
intensity2 = 0;
for k = 1:100
    angle = 2*pi*rand;
    kvec = freq2*[cos(angle) sin(angle)];
    intensity2 = intensity2 ...
        +amp*cos(2*pi*(kvec(1)*xx+kvec(2)*yy+rand));
end

azz = reshape(azz,N*N,1);
ell = reshape(ell,N*N,1);
intensity1 = reshape(intensity1,N*N,1);
intensity2 = reshape(intensity2,N*N,1);

intensity1 = intensity1/max(intensity1);
intensity2 = intensity2/max(intensity2);

figure(1)
set(gcf, 'PaperSize',[10 6])
set(gcf, 'PaperPosition', [0 0 10 6])
subplot(1,2,1)
PlotSphereIntensity(azz,ell,intensity1)
view([-1 1 0])
axis off
subplot(1,2,2)
PlotSphereIntensity(azz,ell,intensity2)
view([-1 1 0])
axis off
%saveas(gcf,'c:/MyStuff/publish/paperGGN/AverageSpheres.eps')