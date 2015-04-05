clear all
%%
model.freq = 10;
model.c_P = 1000;
model.lambda_P = model.c_P/model.freq;
model.rho0 = 2500;
model.n_SC = 100;
model.n_P = 1;
model.seismicnoise = 1e-9/model.freq^1.7; %slope valid between 1 and 10Hz
model.radius = 2*model.lambda_P*ones(1,3);

%symmetric grid
% model.surface = 2*model.lambda_P;
% model.ng = [100 100 100];

%closer to surface
model.surface = 1;%min(1300,model.radius(1));
model.ng = 10*[10 10 round(10*(model.surface+model.radius(1))/(2*model.radius(1)))];

model.samples = 100; %how many values to calculate to estimate coherence
%%
keep = [];
cavi = linspace(2,10,20);

for  k = 2%cavi(end)
    model.cavity = k
    plots.plcnt = 1;
    [output, plots] = isotropicField(model,plots);
    %[output, plots] = scatteredField(model,plots);
    keep = [keep; output.a_std];
end
%%

set(0,'DefaultAxesFontSize',20);
set(0,'DefaultTextFontSize',20);

figure(5)
set(gcf, 'PaperSize',[10 8])
set(gcf, 'PaperPosition', [0 0 10 8])
plot(cavi/model.lambda_P,keep/keep(1),'LineWidth',2)
grid
axis tight
xlabel('Cavity Size [\lambda]')
ylabel('GG Suppression')

%%
%saveas(gcf,'./plots/Cavity_Iso.pdf','pdf')
saveas(gcf,'./plots/Cavity_Sc.pdf','pdf')