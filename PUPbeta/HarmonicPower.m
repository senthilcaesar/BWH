function [SNR,diffFromPitch,prop5Hzdiff] = HarmonicPower(signal,Fs,nfft,ShowFigure,f0xcorr)
% f0xcorr = f0_xcorr(sampStart1);
% f0xcorrPSD_ = f0xcorrPSD(sampStart2);
% signal = buffered(:,sampStart1);

% get PSD
N = length(signal);
window = hamming(N); % hamming(N) %%SS says you might also like to do: window = window/rms(window) to remove power attenuation with your windowing
signal = signal(:) .* window(:); 
y = fft(signal, nfft);
absy = abs(y(1:1000/(Fs/N)));

% sample between peaks
pk2pk_i = ceil(f0xcorr/(Fs/N));
%% find peaks NOT USING FOR NOW

% pkwinstart = 1;
% pkwinstop = nan;
% ii = 0;
% maxIdx = [];
% while 1
%     ii = ii + 1; % counter
%     
%     % find peak
%     pkwinstop = (pkwinstart-1) + (pk2pk_i);
%     if length(absy) - pkwinstart < pk2pk_i % end while loop condition
%         break
%     end 
%     
%     
%     [m, idx] = max(absy(pkwinstart:pkwinstop));
%     maxIdx(ii) = idx+pkwinstart-1;
%     
%     % check forward and backwards and make sure not increasing
%     while absy(maxIdx(ii)+1) > absy(maxIdx(ii)) || absy(maxIdx(ii)-1) > absy(maxIdx(ii))
%         if absy(idx+1) > absy(idx) 
%             idx = idx+1;
%             maxIdx(ii) = idx+pkwinstart-1;
%         elseif absy(idx-1) > absy(idx)
%             idx = idx-1;
%             maxIdx(ii) = idx+pkwinstart-1;
%         end
%     end
%     
%     % reset start peak search window
%     pkwinstart = round(maxIdx(ii)+(pk2pk_i/2)-1);
% end
% 
% %% find troughs NOT USING FOR NOW
% trwinstart = maxIdx(1);
% trwinstop = nan;
% ii = 1;
% minIdx = 1;
% while 1
%     ii = ii + 1; % counter
%     
%     % find peak
%     trwinstop = (trwinstart-1) + (pk2pk_i);
%     if length(absy) - trwinstart < pk2pk_i % end while loop condition
%         break
%     end 
%     
%     
%     [m, idx] = min(absy(trwinstart:trwinstop));
%     minIdx(ii) = idx+trwinstart-1;
%     
%     % check forward and backwards and make sure not increasing
%     while absy(minIdx(ii)+1) < absy(minIdx(ii)) || absy(minIdx(ii)-1) < absy(minIdx(ii))
%         if absy(idx+1) > absy(idx) 
%             idx = idx+1;
%             minIdx(ii) = idx+trwinstart-1;
%         elseif absy(idx-1) > absy(idx)
%             idx = idx-1;
%             minIdx(ii) = idx+trwinstart-1;
%         end
%     end
%     
%     % reset start peak search window
%     trwinstart = round(minIdx(ii)+(pk2pk_i/2)-1);
% end

%% Make sure peaks are between troughs
% should probably incorporate a checking system

%% Hilbert envelop
[~,ylower] = envelope(absy,round(pk2pk_i/2),'peak'); % easier

%% plot
if ShowFigure
    f = Fs*(0:N-1)/N;
    f = f(1:1000/(Fs/N));
%     figure(7), plot(f,absy, f(maxIdx), absy(maxIdx),'g-o',...
%        f(minIdx), absy(minIdx),'r-o')
   figure(7), plot(f,absy, f,yupper,'g-',...
       f,ylower,'r')
end

diffFromPitch = mean(abs(f0xcorr - (diff(maxIdx).*5)));\
prop5Hzdiff = sum(abs(f0xcorr - (diff(maxIdx).*5)) > 5)/length(maxIdx);
TotalPower = trapz(f,absy);
Noise = trapz(f,ylower);
SNR = TotalPower/Noise;

