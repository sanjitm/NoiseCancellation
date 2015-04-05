function NN = arrayFilter(plots,cnt)

%2D dynamical simulation of NN Wiener filter

% 1) Simulate seismic noise field
% 2) Simulate test mass displacement from NN and some other random signal
% 3) Calculate Wiener filter
% 4) Test Wiener filter
% 5) Compare with ideal filter 

set(0,'DefaultAxesFontSize',15);
set(0,'DefaultTextFontSize',15);

NN.params.rho0 = 2500; %density of ground
NN.params.G = 6.673e-11; %Newton's constant

%% surface grid
NN.sgrid.L = 50;
NN.sgrid.nL = 201; %number of grid points along each dimension
NN.sgrid.hTM = 1.5;
NN.sgrid = squareGrid(NN.sgrid);
%% seismometers and seismic spectrum
NN.array(1).L = 15;
NN.array(1).name = 'Spiral, N=20, r=15m';
NN.array(1).noise_amp = 2e-9;
NN.array(1).loc = [];
NN.array(1) = addSpiralSeism(NN.array(1),2,20);

NN.array(2).L = 8;
NN.array(2).name = 'Spiral, N=10, r=8m';
NN.array(2).noise_amp = 2e-9;
NN.array(2).loc = [];
NN.array(2) = addSpiralSeism(NN.array(2),2,10);

% Use the result of optimized array calc to make an array.
load ../SensorArrayOptimization/IterationNumber.mat
load(['../SensorArrayOptimization/ArrayOptimum_IterationNum_' num2str(iterationNumber) '.mat']);
NN.array(3).L = 1;  % Value not used, just filler
NN.array(3).name = ['Optimized Array, N=' num2str(FixedParameters(1))];
NN.array(3).noise_amp = 2e-9;  % Value unused at this time
NN.array(3).loc = [u_minimized' .* cos(theta_minimized'), u_minimized' .* sin(theta_minimized')];

% NN.array(1).L = 15;
% NN.array(1).noise_amp = 2e-9;
% NN.array(1).name = 'Circle, N=20, r=15m';
% NN.array(1).loc = [0 0];
% NN.array(1) = addCircularSeism(NN.array(1),19);

% NN.array(2).L = 8;
% NN.array(2).noise_amp = 2e-9;
% NN.array(2).name = 'Circle, N=10, r=8m';
% NN.array(2).loc = [0 0];
% NN.array(2) = addCircularSeism(NN.array(2),9);

% NN.seism.loc = [0 0];
% NN.seism = addRandSeism(NN.seism,9);
%% time series
NN.time.fs = 100;
NN.time.T = 5;
%filter out data below fc
%(only used in simulation, not for the actual filters)
NN.time.fc = 5; 
NN.time = initTimeSeries(NN.time,NN.array,plots);
%% seismic waves
NN.seism.gnd_amp = 1e-7;%5e-9;
NN.seism.waves = [];

%f,c,N,p
% NN.seism = addPlaneWaves(NN.seism,...
%     [10 14 18 22],[200 200 200 200],[3 3 3 3],[0.25 0.25 0.25 0.25]);

%f,c,R,N,p
% NN.seism = addSphericalWaves(NN.seism,...
%     [10 14 18 22],[200 200 200 200],[50 50 50 50],[3 3 3 3],[0.5 0.5 0.5 0.5]); 

%fc,dT,c,N,p
% NN.seism = addWavelets(NN.seism,NN.time,...
%     [10 14 18 22],[0.2 0.2 0.2 0.2],[200 200 200 200],[3 3 3 3],[0.25 0.25 0.25 0.25]); 

NN.seism = addSphericalWaves(NN.seism,...
    [10 14 18 22],[200 200 200 200],[50 50 50 50],[1 1 1 1],[0.4 0.4 0.4 0.4]); 
NN.seism = addPlaneWaves(NN.seism,...
    [10 14 18 22],[600 600 600 600],[1 1 1 1],[0.1 0.1 0.1 0.1]);
NN.seism = addWavelets(NN.seism,NN.time,...
    [10 14 18 22],[0.2 0.2 0.2 0.2],[200 200 200 200],[2 2 2 2],[0.5 0.5 0.5 0.5]); 
%% NN
NN.time = calculateNN(NN.time,NN.params,NN.sgrid,NN.array,NN.seism,plots,cnt); %fig 1
%% coherence
NN.results = calculateCoh(NN.time,NN.sgrid,NN.array,plots); %fig 2
%% Wiener filter
NN.results = calculateWiener(NN.results,NN.time,NN.array,NN.sgrid.hTM,plots); %fig 3
%% spatial spectrum
NN.results = calculateMap(NN.results,NN.time,NN.array,NN.seism.waves,plots,cnt); %fig 4,5
