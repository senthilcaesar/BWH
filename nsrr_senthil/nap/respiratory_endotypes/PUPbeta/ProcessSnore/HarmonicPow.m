function SNR = HarmonicPow(signal,Fs,nfft,ShowFigure,f0xcorr)
% f0xcorr = f0_xcorr(sampStart1);
% f0xcorrPSD_ = f0xcorrPSD(sampStart2);
% signal = buffered(:,sampStart1);

% get PSD
N = length(signal);
window = hamming(N); % hamming(N)
signal = signal(:) .* window(:);
y = fft(signal, nfft);

absy = abs(y(1:1000/(Fs/N)));
absy = absy - min(absy);
f = Fs*(0:N-1)/N;
f = f(1:1000/(Fs/N));
    
% sample between peaks
pk2pk_i = ceil(f0xcorr/(Fs/N))+1;
numsamps = 3; % number of samples around harmonic to search for peak = 15Hz
% numsamps = 3 is a large search window, and approaching the max we should use 
% Why? because harmonics can be as low as 30Hz and we  do not want two
% harmonic peaks inside 
% A large window means the peak is likely within the search
% window - in this case we can assume if peak does not demonstrate
% peak-like properties, it is likely halfway up a previos or next peak
% see while statement below for fixing these peaks 
% STOPPED THE ABOVE APPROACH
%% find peak (searches within one sample of harmonic




%% TEMP
totalpeaks = floor(f(end)/f0xcorr);
ii = 0;
maxIdx = [];
minIdx = [];

mini = 1+round(f0xcorr/(Fs/N))-2; % start condition
while (mini+pk2pk_i <= length(absy)) % stop looping when indmin is within 2 breaths from the end
    ii = ii + 1; % counter
    incr = pk2pk_i; % set increment
    
    [maxval,maxi] = max(absy(mini:mini+incr)); % search ahead by incr
    maxi = maxi + mini-1;
    maxIdx(ii) = maxi;
    
    if maxi+incr > length(absy) % ends the loop
        [minval,mini] = min(absy(maxi:end));    
        mini = mini + maxi-1;
        minIdx(ii) = mini;
        break
    end

    [minval,mini] = min(absy(maxi:maxi+incr)); % search ahead by incr    
    mini = mini + maxi-1;
    minIdx(ii) = mini;
end

%% find troughs NOT USING FOR NOW
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

%% Hiblert (lower) envelop
% np = 8;
% [~,ylower] = envelope(absy,np,'peak'); % easier

%% Polynomial (lower) envelope
% minIdx = [1, minIdx];
% p = polyfit(f(minIdx)',absy(minIdx),3);
% ylower = polyval(p,f);

%% plot
if ShowFigure
    figure(8), plot(f,absy, f(maxIdx), absy(maxIdx),'g-o',...
       f(minIdx), absy(minIdx),'r-o')
%     figure(7), plot(f,absy, f(maxIdx), absy(maxIdx),'g-o',...
%         f,ylower)
%    figure(7), plot(f,absy,...
%        f,ylower,'r')
end

%% power in harmonics
harmonicpower = 0;
for hnum = 1:length(maxIdx)
    if maxIdx(hnum) == 200
        harmonicpower = harmonicpower + sum(absy(maxIdx(hnum)-1:maxIdx(hnum)));
    elseif maxIdx(hnum) == 1
        harmonicpower = harmonicpower + sum(absy(maxIdx(hnum):maxIdx(hnum)+1));
    else
        harmonicpower = harmonicpower + sum(absy(maxIdx(hnum)-1:maxIdx(hnum)+1));
    end
end

% diffFromPitch = mean(abs(f0xcorr - (diff(maxIdx).*5)));
% prop5Hzdiff = sum(abs(f0xcorr - (diff(maxIdx).*5)) > 5)/length(maxIdx);
% TotalPower = trapz(f,absy);
TotalPower = sum(absy);
% Noise = trapz(f,ylower);
% Noise = sum(ylower);
% SNR = (TotalPower-Noise)/Noise;
SNR = harmonicpower/TotalPower; % multiply by 10 for scaling

