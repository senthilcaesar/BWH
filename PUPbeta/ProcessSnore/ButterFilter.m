function Y = ButterFilter(X,w1,w2,Fs,n)
%This function bandpass filters signal X in the range of [w1,w2] using Butterworth 
%bandpass filter order 2*n.

Wn = [w1 w2]/(Fs/2);
[b,a] = butter(n,Wn);

Y = filter(b,a,X);
end