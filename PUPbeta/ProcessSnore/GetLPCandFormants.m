function [formants,h,f] = GetLPCandFormants(x,Fs,ncoef)
x1 = x.*hamming(length(x));
%         preemph = [1 0.63];
%         x1 = filter(1,preemph,x1);
A = lpc(x1,ncoef);
rts = roots(A);
rts = rts(imag(rts)>=0);
angz = atan2(imag(rts),real(rts));
[frqs,indices] = sort(angz.*(Fs/(2*pi)));
bw = -1/2*(Fs/(2*pi))*log(abs(rts(indices)));
nn = 1;
[h,f]=freqz(1,A,1025,Fs);
formants(nn) = nan;
for kk = 1:length(frqs)
    if (frqs(kk) > 90 && bw(kk) <500)
        formants(nn) = frqs(kk);
        nn = nn+1;
    end
end