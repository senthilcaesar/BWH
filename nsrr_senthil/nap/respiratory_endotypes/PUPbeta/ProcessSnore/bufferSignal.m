function [xin,t] = bufferSignal(x,nx,nwin,noverlap,Fs)
% THIS IS A COPY OF getSTFTColumns --> there was an inconsistency with the
% buffered signal I was making with buffer.m and the buffered signal
% generate by spectrogram. So I decided to copy the spectrogram method.

%getSTFTColumns re-orders input signal into matrix with overlap
%   This function is for internal use only. It may be removed in the future. 
%
%   Copyright 2016 The MathWorks, Inc. 

% Determine the number of columns of the STFT output (i.e., the S output)
ncol = fix((nx-noverlap)/(nwin-noverlap));

coloffsets = (0:(ncol-1))*(nwin-noverlap);
rowindices = (1:nwin)';

% segment x into individual columns with the proper offsets
xin = x(rowindices+coloffsets);

% return time vector whose elements are centered in each segment
t = (coloffsets+(nwin/2)')/Fs;