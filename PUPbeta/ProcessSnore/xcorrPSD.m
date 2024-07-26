function [f0, HNR] = xcorrPSD(signal,Fs)

N = length(signal);
window = hamming(N); % hamming(N)
signal = signal(:) .* window(:);
y = fft(signal, N);
absy = abs(y);
frqres = Fs/N;
f20 = (20/frqres)+1;
f500 = 1000/frqres;
r_PSD = xcorr(absy);
r_PSD = r_PSD(ceil(length(r_PSD)/2):end);
[xcorrcoef,idx] = max(r_PSD(f20:f500));
f0 = frqres*(f20-1+idx-1);
HNR = xcorrcoef/(1-xcorrcoef);