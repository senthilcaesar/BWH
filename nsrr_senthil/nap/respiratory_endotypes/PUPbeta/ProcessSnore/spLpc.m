% NAME
%   spLpc - The linear predictive coding (one-step finite observation 
%           wiener filter prediction)
% SYNOPSIS
%   [a P e] = spLpc(x, fs, ncoef, show)
% DESCRIPTION
%   Obtain LPC coefficients (AR model)
% INPUTS
%   x        (vector) of size Nx1 which contains signal
%   fs       (scalar) the sampling frequency
%   [ncoef]  (scalar) the number of coefficients. The default uses
%              ncoef = 2 + fs / 1000;
%             as a rule of thumb. 
% OUTPUTS
%   a        (vector) of size ncoefx1 which contains LPC coefficients
%   P        (scalar) variance (power) of the prediction error
%   e        (vector) of size Nx1 which contains residual error signals
% AUTHOR
%   Naotoshi Seo, April 2008
% USES
%   lpc.m (Signal Processing toolbox)
function [ffreq,h,f] = spLpc(x, fs, ncoef, show)
 if ~exist('ncoef', 'var') || isempty(ncoef)
     ncoef = 2 + round(fs / 1000); % rule of thumb for human speech is fs/1000, but I reduced to fs/250
 end
 
 if ~exist('show','var')
     show = 0;
 end
 a = lpc(x, ncoef);
 if nargout > 2
    est_x = filter([0 -a(2:end)],1,x);    % Estimated signal
    e = x - est_x;                        % Residual signal
 end 
 
 % find frequencies by root-solving
 r=roots(a); % find roots of polynomial a
 r=r(imag(r)>0.01); % only look for roots >0Hz up to fs/2
 ffreq=sort(atan2(imag(r),real(r))*fs/(2*pi)); % convert to Hz and sort
 
 % Not sure if this is necessary,  but always want ffreq to be at least 3
 % elements
 if length(ffreq) == 1
     ffreq(2) = nan;
     ffreq(3) = nan;
 elseif length(ffreq) == 2
     ffreq(3) = nan;
 end
 
[h,f]=freqz(1,a,1025,fs);
 
 if show == 1
     figure(84)
     % plot signal
     t=(0:length(x)-1)/fs; % times of sampling instants
     subplot(2,1,1);
     plot(t,x);
     legend('Waveform');
     xlabel('Time (s)');
     ylabel('Amplitude');

     % plot frequency response
     subplot(2,1,2);
     plot(f,20*log10(abs(h)+eps));
     legend('LP Filter');
     xlabel('Frequency (Hz)');
     ylabel('Gain (dB)');
%      
%      % temp
%      dbsspecn = tempPxx + ones(1025,1)*(max(20*log10(abs(h)+eps)) ...
%                    - max(tempPxx));
%      hold on
%      plot(f,dbsspecn);
 end


 
end