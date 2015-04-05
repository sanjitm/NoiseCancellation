function plotSeismicSpectrum(sim,NN)

global plcnt

addpath ..
[fl, low, fh, high] = NLNM(1); %Peterson noise models

n = length(sim.run);

seismicspectrum = zeros(size(sim.run(1).seismicspectrum));
for k = 1:n
    seismicspectrum = seismicspectrum+sim.run(k).seismicspectrum.^2;
end
seismicspectrum = sqrt(seismicspectrum/n);

T = NN.time.T;
fs = NN.time.fs;
N = fs*T;

ff = linspace(0,fs/2,N/2+1)';
ff(end) = [];

figure(plcnt)
set(gcf, 'PaperSize',[8 6])
set(gcf, 'PaperPosition', [0 0 8 6])
clf
plcnt = plcnt+1;
loglog(ff,seismicspectrum,'LineWidth',2);
hold on
loglog(fl,low,'k-.',fh,high,'k-.','LineWidth',3)
hold off
grid
axis([ff(10) ff(end) 1e-12 1e-7])
xlabel('Frequency [Hz]')
ylabel('Seismic spectrum [m/\surd Hz]')