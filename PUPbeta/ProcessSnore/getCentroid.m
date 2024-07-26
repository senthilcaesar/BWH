function centroid=getCentroid(buffered,Fs,winsize)

N = floor(winsize/2) + 1;
f = Fs/2*linspace(0,1,N)';

%First apply a hamming window to each signal window
w = hamming(winsize);
Ham = repmat(w,[1,size(buffered,2)]);
buffered = buffered .* Ham;

%Take the FFT of each window
ffts = fft(double(buffered),winsize);
clear buffered;

%Keep only the first N samples (buffered should be a real signal)
ffts = ffts(1:N,:);
mags = abs(ffts);
f=repmat(f,[1 size(mags,2)]);

%Apply the formula to calculate centroid
centroid=sum(f.*mags)./sum(mags);