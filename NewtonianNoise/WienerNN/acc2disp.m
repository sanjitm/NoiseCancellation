function data = acc2disp(data,T,fs,fc)

SN = length(data(1,:));

ff = linspace(0,fs/2,T*fs/2+1)';
ff(end) = [];
ff = ff(:,ones(1,SN));

amp = fft(data);
amp = amp(1:(end/2),:);

amp = amp./(2*pi*ff).^2;

amp = [zeros(1,SN); amp(2:end,:); zeros(1,SN); conj(flipud(amp(2:end,:)))];
data = ifft(amp);

[b, a] = butter(5,2*fc/fs,'high');
data = filtfilt(b, a, data);