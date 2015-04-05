function plotSubResiduals(sim,NN)

global plcnt

T = NN.time.T;
fs = NN.time.fs;
N = fs*T;

ff = linspace(0,fs/2,N/2+1)';
ff(end) = [];

AN = length(sim.run(1).results.array);

leg = {'Noise model','Newtonian noise'};
for ai = 1:AN
    leg = [leg{:} strcat(NN.array(ai).name,{', res.'})];
end

figure(plcnt)
set(gcf, 'PaperSize',[8 6])
set(gcf, 'PaperPosition', [0 0 8 6])
clf
plcnt = plcnt+1;
loglog(ff,2*[NN.time.noisemodel sim.NN_ff sim.sub_ff]/4000,'LineWidth',2);
legend(leg,'Location','NorthEast')
grid
axis([NN.time.fc+1 ff(end) 1e-23 1e-20])
xlabel('Frequency [Hz]')
ylabel('Strain sensitivity [1/\surd Hz]')
