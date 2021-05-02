function chirpFFT(psthMatrix,drawStruct)
if size(psthMatrix,1) > 1 
   likelymaxfreqIndex = size(psthMatrix,1);
else
    likelymaxfreqIndex = 1;
end
frequencyStrip = psthMatrix(likelymaxfreqIndex,.35*(10^5):1.35*(10^5));






%FFT STUFF BELOW**%**%**%**%** 
fs=10000; %sampling frequency

L=length(frequencyStrip);
NFFT = 100000;
X = fftshift(fft(frequencyStrip,NFFT));
Pxx=X.*conj(X)/(NFFT*NFFT); %computing power with proper scaling
f = fs*(-NFFT/2:NFFT/2-1)/NFFT; %Frequency Vector

figure(10)
plot(f,abs(X)/(L),'r');
title('Magnitude of FFT');
xlabel('Frequency (Hz)')
ylabel('Magnitude |X(f)|');
xlim([1 drawStruct(size(drawStruct,2)).frequencyMax])

fMag = f';
yMag = (abs(X)/(L))';

% Pxx=X.*conj(X)/(NFFT*NFFT); %computing power with proper scaling
% figure(11)
% plot(f,10*log10(Pxx),'r');
% title('Double Sided - Power Spectral Density');
% xlabel('Frequency (Hz)')
% ylabel('Power Spectral Density- P_{xx} dB/Hz');
% xlim([-100 100])

X = fft(frequencyStrip,NFFT);
X = X(1:NFFT/2+1);%Throw the samples after NFFT/2 for single sided plot
Pxx=X.*conj(X)/(NFFT*NFFT);
f = fs*(0:NFFT/2)/NFFT; %Frequency Vector
figure(12)
plot(f,10*log10(Pxx),'r');
title('Single Sided - Power Spectral Density');
xlabel('Frequency (Hz)')
ylabel('Power Spectral Density- P_{xx} dB/Hz');
xlim([1 drawStruct(size(drawStruct,2)).frequencyMax])

fPower = f';
logY=(10*log10(Pxx))';

save('chirpFFT.mat','fMag','yMag','fPower','logY')
end
