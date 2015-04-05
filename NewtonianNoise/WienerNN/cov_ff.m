function cov_data = cov_ff(data,T,fs,noisemodel)

DN = length(data(1,:));

ff = linspace(0,fs/2,T*fs/2+1);
ff(end) = [];

amps = fft(data)/fs;
amps = amps(1:end/2,:)./noisemodel(:,ones(1,DN));

data = ifft([amps; zeros(1,DN); conj(flipud(amps(2:end,:)))])*fs;

cov_data = cov(data);




