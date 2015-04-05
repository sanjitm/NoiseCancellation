function S = prefilter(S,T,fs,hTM)

SN = length(S(1,:));
c_ref = 210; %assumed speed of seismic waves

ff = linspace(0,fs/2,T*fs/2+1)';
ff(end) = [];
ff = ff(:,ones(1,SN));

amp = fft(S);
amp = amp(1:(end/2),:);

%take into account acc->disp and exponential suppression
%amp = amp./(2*pi*ff).^2;
amp = amp./(2*pi*ff).^2.*exp(-2*pi*ff*hTM/c_ref);

amp = [zeros(1,SN); amp(2:end,:); zeros(1,SN); conj(flipud(amp(2:end,:)))];
S = ifft(amp);

[b, a] = butter(5,2*5/fs,'high');
for si = 1:SN
    S(:,si) = filtfilt(b, a, S(:,si));
end