function time = initTimeSeries(time,array,plots)

global plcnt

AN = length(array);

N = time.fs*time.T;
T = time.T;
fs = time.fs;

time.tt = linspace(0,T,N)';

for ai = 1:AN
    aS = array(ai).noise_amp;
    SN = length(array(ai).loc(:,1));
    time.array(ai).noiseS = aS*randn(N,SN); %seismometer noise
end
size(time.array(3).noiseS);
%% simulate something like the low-frequency aLIGO noise
A = detrend(randn(N,1),'linear'); %start with white noise

norm_A = mean(abs(fft(A)/fs/sqrt(T))); %normalize noise power

ff = linspace(0,fs/2,N/2+1)';
ff(end) = [];

%this is the modelled noise spectrum
time.noisemodel = 4000*(5e-23./(ff/10).^14+7e-23./(ff/10).^2+4.5e-24./(ff/40).^(1/50));

[DummyVar,fi] = min(abs(ff-time.fc)); %low-frequency cutoff to get nicer results

%filter white noise by noise model
a = fft(A)/fs/norm_A; %rPSD equivalent amplitude
a = a(1:N/2).*time.noisemodel;
a(1:fi) = 0;
a = [a; 0; conj(flipud(a(2:end)))];

time.noiseTM = ifft(a)*fs; %test mass noise

if plots.noisemodel
    figure(plcnt)
    set(gcf, 'PaperSize',[8 6])
    set(gcf, 'PaperPosition', [0 0 8 6])
    clf
    plcnt = plcnt+1;
    a = fft(time.noiseTM)/fs;
    loglog(ff(fi+1:end),abs(a(fi+1:end/2))/sqrt(T),...
        ff(fi+1:end),abs(time.noisemodel(fi+1:end)),'LineWidth',2)
    grid
    axis tight
    set(gca,'ylim',[1e-21 1e-17])
    xlabel('Frequency [Hz]')
    ylabel('Noise spectrum, TM [m/\surd Hz]')
end