clear all
close all
global plcnt

plots.arrayresponse = 0;
plots.spatialspectra = 0;
plots.arraycoherence = 1;
plots.displacementfield = 0;
plots.seismicfieldmovie = 0;
plots.gravityfieldmovie = 0;
plots.timeseries = 0;
plots.noisemodel = 0;

for k = 1:3
    plcnt = 1;
    NN = arrayFilter(plots,k);
    sim.run(k).results = NN.results;
    sim.run(k).NN = NN.time.NN;
    sim.run(k).seismicspectrum = NN.time.seismicspectrum;
end

sim = simProcessing(sim,NN.time);
%%
plotSubFactors(sim,NN)
%%
plotSubResiduals(sim,NN)
%%
plotSeismicSpectrum(sim,NN);
%%
% figure(1)
% saveas(gcf,'./plots/NoiseModel_TM.jpg')
% saveas(gcf,'./plots/NoiseModel_TM.pdf')

% figure(1)
% saveas(gcf,'./plots/Field_SW4f3.jpg')
% saveas(gcf,'./plots/Field_SW4f3.pdf')
%
% figure(2)
% saveas(gcf,'./plots/Data_WA4f2_SW4f1.jpg')
% saveas(gcf,'./plots/Data_WA4f2_SW4f1.pdf')

% figure(7)
% saveas(gcf,'./plots/Map_10Hz_SW4f3.jpg')
% saveas(gcf,'./plots/Map_10Hz_SW4f3.pdf')
% 
% figure(11)
% saveas(gcf,'./plots/Perf_rel_Spiral20_WA4f3.jpg')
% saveas(gcf,'./plots/Perf_rel_Spiral20_WA4f3.pdf')
% 
% figure(13)
% saveas(gcf,'./plots/Perf_abs_Spirals_WA4f3.jpg')
% saveas(gcf,'./plots/Perf_abs_Spirals_WA4f3.pdf')

% figure(5)
% saveas(gcf,'./plots/Residuals_Spirals.jpg')
% saveas(gcf,'./plots/Residuals_Spirals.pdf')
% saveas(gcf,'./plots/Residuals_Spirals.eps')

% figure(5)
% saveas(gcf,'./plots/Residuals_Circles.jpg')
% saveas(gcf,'./plots/Residuals_Circels.pdf')
% saveas(gcf,'./plots/Residuals_Circels.eps')

% figure(6)
% saveas(gcf,'./plots/Perf_abs_WA4f2_SW4f1_acex.jpg')
% saveas(gcf,'./plots/Perf_abs_WA4f2_SW4f1_acex.pdf')
% 
% figure(6)
% saveas(gcf,'./plots/Subtraction_ff.jpg')
% saveas(gcf,'./plots/Subtraction_ff.pdf')