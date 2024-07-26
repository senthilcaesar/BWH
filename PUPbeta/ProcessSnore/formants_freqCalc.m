function [formantsfreq] = formants_freqCalc (segments,fs)
F123 = zeros(size(segments,1),3);
for s=1:size(segments,1)
  x1 = segments(s,:)'.*hamming(length(segments(s,:))); %window the breath segment using a Hamming window
  preemph = [1 0.63]; % pre-emphasis  %this is a highpass all-pole (AR(1)) filter 
  x1 = filter(1, preemph, x1);
  A = lpc(x1, 8); %we are obtaining a linear prediction coeffiicients. The model order is 2 times the expected formants plus 2. (2*3+2)
  rts = roots(A);
  rts = rts(imag(rts)>=0); %LPC coefficients are real-valued, the roots occur in complex conjugate pairs.
  angz = atan2(imag(rts),real(rts));
  [frqs,indices] = sort(angz.*(fs/(2*pi))); %angles corresponding to the roots
  bw = -1/2*(fs/(2*pi))*log(abs(rts(indices))); %convert the angular frequencies in rad/samples 
  nn = 1;
  formants = [0 0 0];
  for kk = 1:length(frqs)
    if (frqs(kk) > 90 && bw(kk) <400) % formant criteria
      formants(nn) = frqs(kk);
      nn = nn+1;
    end
  end

  F123(s,1) = formants(1);
  F123(s,2) = formants(2);
  F123(s,3) = formants(3);
end
formantsfreq = F123;