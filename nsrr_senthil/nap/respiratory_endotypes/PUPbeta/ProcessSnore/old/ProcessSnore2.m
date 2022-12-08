function [SnoreInterpStruct, SnoreStruct] = ProcessSnore2(Snore,Flow,Time,ChannelsList,Channels_Fs,PlotSnore)

disp('Processing Snore Data')
%% Initilize settings
% Fs = 1/Snore.interval;
% SndData = Snore.values;

% Data
Fs = Channels_Fs(strcmp(ChannelsList,'SnoreRaw'));
FsFlow = Channels_Fs(strcmp(ChannelsList,'Flow'));
SndData = Snore;
f = 1025; % highest frequency for spectrogram
StartTime = Time(1);

% Window size
freqRes = 5;
timeRes = 1/freqRes; %Ali says you want this to be 0.1 for snoring
winsize = floor(Fs/freqRes); %window size
percoverlap = 0.75;
noverlap = floor(winsize*percoverlap);  %75% overlap
nfft = max(256,2^nextpow2(winsize));  %the DFT length, which can be longer than the window
ncoefs = 13; %Number of MFCCs 
        
timeResSm = 0.02;
winsizeSm = floor(Fs*timeResSm); % timeRes = 0.025
percoverlapSm = 0.5;
noverlapSm = floor(winsizeSm*percoverlapSm);  % 50% overlap

% winsizeT = 0.5;
% winsizeT50 = 0.1;
% percOverlap = 50;
% FsDownsample=6000;
% FsDS2 = F_samp;
FsBuffer = 1/((1-percoverlap)*winsize/Fs);
FsBufferSm = 1/((1-percoverlapSm)*winsizeSm/Fs);
% FsBuffer50 = 1/(winsizeT50*(percOverlap/100));
SegLen = 10*60; %seconds

WinVarArray = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
ZCrArray = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
timeWin = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
formant1 = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
formant2 = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
formant3 = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
formant1_n = nan(floor(length(SndData)*(FsBuffer/Fs)),1); %testing
formant2_n = nan(floor(length(SndData)*(FsBuffer/Fs)),1); %testing
formant3_n = nan(floor(length(SndData)*(FsBuffer/Fs)),1); %testing
p0to125 = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
p125to250 = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
p250to500 = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
p500to1k = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
p1kto1k5 = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
p1k5to2k = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
p2kto2k5 = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
p2k5to3k = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
ptot = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
pkpower = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
pk99power = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
f0 = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
HNR = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
HNR2 = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
HarmonicPower  = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
centroid = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
timeWin = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
SpectP = nan(floor(length(SndData)*(FsBuffer/Fs)),f);
SpectPabs = nan(floor(length(SndData)*(FsBuffer/Fs)),f);
spectFlux = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
spectRollOff = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
spectSpread = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
spectSkewness = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
spectKurtosis = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
spectFlatness = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
spectSlope = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
spectDecrease = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
spectEntropy = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
spectCrest = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
medianFreq = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
bw = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
flo = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
fhi = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
coeffsAvg = nan(floor(length(SndData)*(FsBuffer/Fs)),ncoefs);
delAvg = nan(floor(length(SndData)*(FsBuffer/Fs)),ncoefs);
ddelAvg = nan(floor(length(SndData)*(FsBuffer/Fs)),ncoefs);
coeffsAvg5k = nan(floor(length(SndData)*(FsBuffer/Fs)),ncoefs); %testing
delAvg5k = nan(floor(length(SndData)*(FsBuffer/Fs)),ncoefs); %testing
ddelAvg5k = nan(floor(length(SndData)*(FsBuffer/Fs)),ncoefs); %testing
lpcenv = nan(floor(length(SndData)*(FsBuffer/Fs)),f);
lpcenv_n = nan(floor(length(SndData)*(FsBuffer/Fs)),f);


% Small window vars
ZCrTransform = nan(floor(length(SndData)*(FsBufferSm/Fs)),1);
timeWinSm = nan(floor(length(SndData)*(FsBufferSm/Fs)),1);
coeffsArray = nan(floor(length(SndData)*(FsBufferSm/Fs)),ncoefs);
delArray = nan(floor(length(SndData)*(FsBufferSm/Fs)),ncoefs);
ddelArray = nan(floor(length(SndData)*(FsBufferSm/Fs)),ncoefs);
coeffsArray5k = nan(floor(length(SndData)*(FsBufferSm/Fs)),ncoefs); %  testing
delArray5k = nan(floor(length(SndData)*(FsBufferSm/Fs)),ncoefs); % testing
ddelArray5k = nan(floor(length(SndData)*(FsBufferSm/Fs)),ncoefs); % testing

%% Process sounds
% isolate 10 min segments
% possibly temporary but might be faster to work in sections
% in which case we loop everything below
for ii = 1:ceil(length(SndData)/(SegLen*Fs))
    %%
    PercAnalysisComplete = 100*ii/ceil(length(SndData)/(SegLen*Fs));
    disp(['Percent complete: ', num2str(round(PercAnalysisComplete)),'%'])
    %%
    start = (ii-1)*SegLen*Fs + 1;
    finish = ii*SegLen*Fs + (percoverlap)*winsize;
    
    if finish > length(SndData) % end condition
        finish = length(SndData);
    end
    
    SndData10Temp = SndData(start:finish);
    t=(start:finish)'/Fs;
    
    startFlow = round((ii-1)*SegLen*FsFlow + 1);
    finishFlow = round(ii*SegLen*FsFlow);
    
    if finishFlow > length(Flow) % end condition
        finishFlow = length(Flow);
    end
    
    Flow10 = Flow(startFlow:finishFlow);
    timeFlow = (startFlow:finishFlow)'/FsFlow;

    %% Filter
    FilterBandwidth = [5 2000];
    SndData10 = SndData10Temp;
%     FilteredSnd = ButterFilter(SndData10,FilterBandwidth(1),FilterBandwidth(2),Fs,5);
    FilteredSnd = SndData10;
    
    %% Downsample TURNED OFF
    if 0
    FilteredSndDS = resample(FilteredSnd,FsDownsample,Fs);
    startDS = (ii-1)*SegLen*FsDownsample + 1; finishDS = ii*SegLen*FsDownsample;
    timeDS(startDS:finishDS) = (startDS:finishDS)'/FsDownsample;
    end

    clear SndData10Temp
    
    %% Buffer data
    % dividing the signal into overlapping segments
%     winsize = round(winsizeT*FsDownsample);
%     winOverlap = floor((percOverlap/100)*winsize);

    %Put signal into windows
    num_windows = floor(length(FilteredSnd)/noverlap) - 1;
    total_length = noverlap * (num_windows + 1);
%     bufferedOld = buffer(FilteredSnd(1:total_length),winsize,noverlap,'nodelay'); % This prevents the buffer function from prepending and appending zeros
    [buffered, t_buffered] = bufferSignal(FilteredSnd,length(FilteredSnd),winsize,noverlap,Fs);
    
    % start condition
    if ii == 1
        startWin = 1;
    else
        startWin = finishWinLastWindow + 1;
    end
    
    finishWin = startWin + size(buffered,2)-1; %floor((ii*SegLen*FsBuffer)-1);
    finishWinLastWindow = finishWin;
    
    if finishWin > length(WinVarArray)
        finishWin = startWin + size(buffered,2)-1;
    end
    
    %%%% Adapted from bufferSignal code
    % Determine the number of columns of the STFT output (i.e., the S output)
    ncol = fix((length(FilteredSnd)-noverlap)/(winsize-noverlap));
    coloffsets = (0:(ncol-1))*(winsize-noverlap);
%     rowindices = (startWin:finishWin)';

    % return time vector whose elements are centered in each segment
    timeWin(startWin:finishWin) = StartTime + (ii-1)*SegLen + (coloffsets+(winsize/2)')/Fs;
    
    %% Buffer data in smaller windows
    % this is for features that need smaller windows (e.g. making a ZCr
    % transformed signal, from which to calculate sample entropy)
       
    % Smaller window buffered signal
    [bufferedSm, ~] = bufferSignal(FilteredSnd,length(FilteredSnd),floor(Fs*timeResSm),floor(Fs*percoverlapSm*timeResSm),Fs);
    
    % start condition
    if ii == 1
        startWinSm = 1;
    else
        startWinSm = finishWinLastWindowSm + 1;
    end
    
    finishWinSm = startWinSm + size(bufferedSm,2)-1; %floor((ii*SegLen*FsBuffer)-1);
    finishWinLastWindowSm = finishWinSm;
    
    if finishWinSm > length(ZCrTransform)
        finishWinSm = startWinSm + size(bufferedSm,2)-1;
    end
    
    %%%% Adapted from bufferSignal code
    % Determine the number of columns of the STFT output (i.e., the S output)
    ncol = fix((length(FilteredSnd)-noverlapSm)/(winsizeSm-noverlapSm));
    coloffsetsSm = (0:(ncol-1))*(winsizeSm-noverlapSm);
%     rowindices = (startWin:finishWin)';

    % return time vector whose elements are centered in each segment
    timeWinSm(startWinSm:finishWinSm) =  StartTime + (ii-1)*SegLen + (coloffsetsSm+(winsizeSm/2)')/Fs;
    
    %% Calculate features
    % Sound amplitude - might be undersampled @ 4Hz at the moment
%     SndIntensity = (FilteredSnd).^2; % converted to intensity
%     SndIntensitydB = 10.*log10(SndIntensity/1e-12); %reference = 1e-12 W/m2
%     SndIntensitydBDS = resample(SndIntensitydB,FsBuffer,Fs); % downsample to F_samp Hz
%     startDS2 = (ii-1)*SegLen*FsDS2 + 1; finishDS2 = ii*SegLen*FsDS2;
    
%     if finishDS2 > length(SndIntensityArray) % end condition
%         finishDS2 = startDS2 + length(SndIntensitydBDS)-1;
%     end
%     timeDS2(startDS2:finishDS2) = (startDS2:finishDS2)'/FsDS2;
    
%     SndIntensityArray(startWin:finishWin) = SndIntensitydB; % FIX THIS LATER
    
    % Windowed Variance, Zero-crossing rate
%     WinVar50 = nanvar(buffered50,1);
    WinVar = nanvar(buffered,1);
    MSs = nanmean(buffered,1);
    ZSs = buffered - MSs;
    ZCr = sum(abs(sign(ZSs(1:end-1,:))-sign(ZSs(2:end,:))))/(2*winsize);
    
%     WinVarArray50(startWin50:finishWin50) = WinVar50';
    WinVarArray(startWin:finishWin) = WinVar';
    ZCrArray(startWin:finishWin) = ZCr';
    
    % Band power (absolute)
    % Peak power is computed below where the PSD is computed, this can be
    % used to normalize
    % tested this and its the exact same answer if I calculate from ps
    % (saved in SpectP so I can do these post process too)
    % e.g Fdiff = diff(F); p0to125 = sum(ps(f0:f125,:),1)*Fdiff(1);
    p0to125(startWin:finishWin) = bandpower(buffered,Fs,[0 125]);
    p125to250(startWin:finishWin) = bandpower(buffered,Fs,[125 250]);
    p250to500(startWin:finishWin) = bandpower(buffered,Fs,[250 500]);
    p500to1k(startWin:finishWin) = bandpower(buffered,Fs,[500 1000]);
    p1kto1k5(startWin:finishWin) = bandpower(buffered,Fs,[1000 1500]);
    p1k5to2k(startWin:finishWin) = bandpower(buffered,Fs,[1500 2000]);
    p2kto2k5(startWin:finishWin) = bandpower(buffered,Fs,[2000 2500]);
    p2k5to3k(startWin:finishWin) = bandpower(buffered,Fs,[2500 3000]);
    ptot(startWin:finishWin) = bandpower(buffered,Fs,[0 3000]);
    
    % Spectral centroid
    centroid(startWin:finishWin) = getCentroid(buffered,Fs,winsize);
    
   %% ZCr tranformed signal and entropy
    % Generate ZCr transformed signal
    MSs2 = nanmean(bufferedSm,1);
    ZSs2 = bufferedSm - MSs2;
    ZCr2 = sum(abs(sign(ZSs2(1:end-1,:))-sign(ZSs2(2:end,:))))/(2*(Fs*0.025));
    ZCrTransform(startWinSm:finishWinSm) = ZCr2;
    
        
    %% Pitch and harmonic ratio
    %%% NEED TO TRY AND ADAPT TO BUFFERED DATA %%%
    % Get pitch from ceptral method
    maxf = 500; % max f to search for harmonics
    minf = 20; % min f to search for harmonics
    ms2=floor(Fs/maxf); % 2ms (500Hz)
    ms20=floor(Fs/minf); % 50ms (20Hz)
    f0_cep = nan(size(buffered,2),1);
    f0_xcorr = nan(size(buffered,2),1);
    f0_xcorrPSD = nan(size(buffered,2),1);
    HNRcorr = nan(size(buffered,2),1);
    HNRcep = nan(size(buffered,2),1);
    HNRcorrPSD = nan(size(buffered,2),1);
    HarmPow = nan(size(buffered,2),1);
    lpcenvtemp = nan(size(buffered,2),f);
    lpcenvtemp_n = nan(size(buffered,2),f);
    
    for winnum = 1:size(buffered,2)
        %% find pitch
        signal = buffered(:,winnum); % sample signal
        
        % cepstral method - tends to fail in finding pitch, which results
        % in weird HNR values
        [c, ~] = spCepstrum(signal, Fs, 'hamming', 0);
        HNRcep(winnum,1) = max(abs(c(ms2:ms20)))/nanmean(abs(c(ms2:ms20)));
        
        [~,idx]=max(abs(c(ms2:ms20)));
        f0_cep(winnum,1) = Fs/(ms2+idx-1);
        
        % Autcorrelation method - this seems to be most effective at
        % finding pitch but fails with HNR calculation in specific
        % circumstances
        % Do autocorrelation
        [r] = spCorrelum(signal, Fs, Fs/20, 0);
        % half is just mirror for real signal
        r = r(floor(length(r)/2):end);
        [Rxx_T0_,idx]=max(r(ms2:ms20)); % Rxx_T0coeff is very similar to below
        
        % find pitch
        f0_xcorr(winnum,1) = Fs/(ms2+idx-1);
        
        % modified ceptstral to get HNR (using pitch found with autocorr)
%         s1 = floor(Fs/f0_xcorr(winnum,1))-5;
%         s2 = floor(Fs/f0_xcorr(winnum,1))+5; 
%         [~,idx]=max(abs(c(s1:s2)));
%         f0_cepMod(winnum,1) = Fs/(s1+idx-1);
%         HNRcepMod(winnum,1) = abs(c(idx))/nanmean(abs(c(ms2:ms20)));
        
        % find harmonic ratio
        % my method
        signal = buffered(:,winnum);
        Rxx_ = xcorr(signal);
        Rxx = Rxx_(floor(length(Rxx_)/2):end); % keep second half
        Rxx_0 = Rxx(1); % power of the entire signal
        [Rxx_T0, ~] = max(Rxx(ms2:ms20));
        HNRcorr(winnum,1) = Rxx_T0/(Rxx_0-Rxx_T0);
               
        % method using autocorrelation on PSD
        [f0_xcorrPSD(winnum,1), HNRcorrPSD(winnum,1)] = xcorrPSD(signal,Fs);
        
        % method computing ration total power to noise floor
        ShowFigure = 0;
        HarmPow(winnum,1) = HarmonicPow(signal,Fs,nfft,ShowFigure,f0_xcorr(winnum));
        
        %% Formants ( move up after pitch detection when testing complete)
        signal = buffered(:, winnum);
        [ffreq,h,f_] = spLpc(signal, Fs, 16, 0);
        formant1(winnum,1) = ffreq(1);
        formant2(winnum,1) = ffreq(2);
        formant3(winnum,1) = ffreq(3);
        lpcenvtemp(winnum,:) = 20*log10(abs(h+eps));
             
        if 0
            figure(24)
            [~,f_1] = min(abs(formant1(winnum,1)-f_));
            [~,f_2] = min(abs(formant2(winnum,1)-f_));
            [~,f_3] = min(abs(formant3(winnum,1)-f_));
            plot(f_,lpcenvtemp(winnum,:),...
                formant1(winnum,1), lpcenvtemp(winnum,f_1),'kx',...
                formant2(winnum,1), lpcenvtemp(winnum,f_2),'ro',...
                formant3(winnum,1), lpcenvtemp(winnum,f_3),'gd');
        end
        
        %%
        % different formant algo
        ncoef = 12;
        formants = nan(1,10);
        [formants_, h,f_n] = GetLPCandFormants(signal,Fs,ncoef);
        formants(1:length(formants_)) = formants_;
        formant1_n(winnum,1) = formants(1);
        formant2_n(winnum,1) = formants(2);
        formant3_n(winnum,1) = formants(3);
        lpcenvtemp_n(winnum,:) = 20*log10(abs(h+eps));

        %plot
        if 0
            figure(25)
            [~,f_1] = min(abs(formants(1)-f_));
            [~,f_2] = min(abs(formants(2)-f_));
            [~,f_3] = min(abs(formants(3)-f_));
            plot(f_,lpcenv,...
                formants(1), lpcenv(f_1),'kx',...
                formants(2), lpcenv(f_2),'ro',...
                formants(3), lpcenv(f_3),'gd',...
                f_,lpcenvtemp(winnum,:),'k');
        end    
    end 

    % Only use pitch data similar between both methods
%     similarIdx = abs(f0_cep - f0_xcorr) < 999; sum(similarIdx); % Value of this parameter not tested
%     f0_xcorr(~similarIdx) = nan;
%     HNRcorr(~similarIdx) = nan; % removes data with low harmonic to noise ratio - not using this right now
%     f0_xcorr(HNRcorr<0.5 | SNR < 2.5) = nan;
    f0(startWin:finishWin) = f0_xcorr;
    HNR(startWin:finishWin) = HNRcorr;
    HNR2(startWin:finishWin) = harmonicRatio(buffered,Fs,'Window',hamming(winsize),'OverlapLength',noverlap);
    HarmPow = smooth(HarmPow,4);
    HarmonicPower(startWin:finishWin) = HarmPow;
    lpcenv(startWin:finishWin,:) = lpcenvtemp;
    lpcenv_n(startWin:finishWin,:) = lpcenvtemp_n;
    
    % Spectrogram
    figure('Visible','off');    
%     spectrogram(FilteredSnd,winsize,noverlap,nfft,Fs,'yaxis');
    [X,F,T,ps] = spectrogram(FilteredSnd,winsize,noverlap,nfft,Fs,'yaxis');
    
    % compute magnitude spectrum - used in some spectral features later
    X = abs(X)*2/nfft;
    X([1 end],:)= X([1 end],:)/sqrt(2); % normalization
    
    % compute PSD in db scale
    P = 10*log10(ps+eps);
    F = F/1000; % convert to kHz
    clim = [min(min(P)) max(max(P))];
    
    % get data to rebuild spectrogram
%     ch = get(gca,'ch');
%     ax = gca;
%     T = get(ch, 'XData') + SegLen*(ii-1);
%     F = get(ch, 'YData');
%     P = get(ch,'CData');
% %     t = linspace(0, length(FilteredSndDS), length(FilteredSndDS));
%     clim = ax.CLim;
    
    SpectP(startWin:finishWin,:) = P'; % store spectral power
    SpectPabs(startWin:finishWin,:) = ps';
    
    % Peak power (to normalize bandpowers)
    pk99power(startWin:finishWin) = prctile(ps,99,1);
    pkpower(startWin:finishWin) = prctile(ps,100,1);
    
    %% MFCC 
    % using smaller window to match best practices - need to do a
    % moving average to match sampling rate of other features
    % note: if maxfreq = 2000, filters are spaced ~47Hz apart
    % if maxfreq = 5000, filters are spaced ~62Hz apart - not clear
    % there should be a preference, went with 2k, can test 5k (which is
    % max for Fs=10k).

    [coeffs,~,~] = melfcc(FilteredSnd, Fs, 'wintime',0.020,'hoptime',0.01,'maxfreq',2000); 
    del = deltas(coeffs);
    % Double deltas are deltas applied twice with a shorter window
    ddel = deltas(deltas(coeffs,5),5);

    coeffsArray(startWinSm:finishWinSm,1:ncoefs) = coeffs';
    delArray(startWinSm:finishWinSm,1:ncoefs) = del';
    ddelArray(startWinSm:finishWinSm,1:ncoefs) = ddel';


    % buffer signal and average to make same resolution as other
    % features
    FsCoeffs = Fs/(0.01/(1/Fs));

    % deal with issues in array size
    [bufferedCoeffs, ~] = bufferSignal(coeffs(1,:),length(coeffs(1,:)),10,5,FsCoeffs); % sample output 
    finishWinBuf = startWin + size(bufferedCoeffs,2)-1;

    % NEED TO FIX THIS
    for coefnum = 1:size(coeffs,1)
        [bufferedCoeffs, ~] = bufferSignal(coeffs(coefnum,:),length(coeffs(coefnum,:)),10,5,FsCoeffs); % still a couple frames long
        coeffsAvg(startWin:finishWinBuf, coefnum) = mean(bufferedCoeffs,1);

        [bufferedDel, ~] = bufferSignal(del(coefnum,:),length(del(coefnum,:)),10,5,FsCoeffs); % still a couple frames long
        delAvg(startWin:finishWinBuf, coefnum) = mean(bufferedDel,1);

        [bufferedDdel, ~] = bufferSignal(ddel(coefnum,:),length(ddel(coefnum,:)),10,5,FsCoeffs); % still a couple frames long
        ddelAvg(startWin:finishWinBuf,coefnum) = mean(bufferedDdel,1);
    end
    %% MFCC 
    % 5k
    [coeffs,~,~] = melfcc(FilteredSnd, Fs, 'wintime',0.020,'hoptime',0.01,'maxfreq',5000); 
    del = deltas(coeffs);
    % Double deltas are deltas applied twice with a shorter window
    ddel = deltas(deltas(coeffs,5),5);

    coeffsArray5k(startWinSm:finishWinSm,1:ncoefs) = coeffs';
    delArray5k(startWinSm:finishWinSm,1:ncoefs) = del';
    ddelArray5k(startWinSm:finishWinSm,1:ncoefs) = ddel';


    % buffer signal and average to make same resolution as other
    % features
    FsCoeffs = Fs/(0.01/(1/Fs));

    % deal with issues in array size
    [bufferedCoeffs, ~] = bufferSignal(coeffs(1,:),length(coeffs(1,:)),10,5,FsCoeffs); % sample output 
    finishWinBuf = startWin + size(bufferedCoeffs,2)-1;

    % NEED TO FIX THIS
    for coefnum = 1:size(coeffs,1)
        [bufferedCoeffs, ~] = bufferSignal(coeffs(coefnum,:),length(coeffs(coefnum,:)),10,5,FsCoeffs); % still a couple frames long
        coeffsAvg5k(startWin:finishWinBuf, coefnum) = mean(bufferedCoeffs,1);

        [bufferedDel, ~] = bufferSignal(del(coefnum,:),length(del(coefnum,:)),10,5,FsCoeffs); % still a couple frames long
        delAvg5k(startWin:finishWinBuf, coefnum) = mean(bufferedDel,1);

        [bufferedDdel, ~] = bufferSignal(ddel(coefnum,:),length(ddel(coefnum,:)),10,5,FsCoeffs); % still a couple frames long
        ddelAvg5k(startWin:finishWinBuf,coefnum) = mean(bufferedDdel,1);
    end

    %% Spectral flux
    % the rate at which spectrum is changing between two consecutive
    % frames. Calculated as the squared difference between normalized
    % spectra of two consecutive frames.
    % algorithm from: https://www.audiocontentanalysis.org/code/audio-features
    % this has a little bug caused by analyzing in segements
    spectFlux(startWin:finishWin) = FeatureSpectralFlux (X, Fs);

    %% Spectral roll off
    % Measures "skewness" of the spectral shape. Concretely, it
    % represenents the frequency at which N% of the power is
    kappa = 0.85;
    spectRollOff(startWin:finishWin) = FeatureSpectralRolloff(X, Fs, kappa);

    %% Spectral spread
    % The standard deviation around the spectral centroid 
    spectSpread(startWin:finishWin) = spectralSpread(ps,F*1000);
    
    %% Spectral skewness 
    % The spectral skewness measures symmetry around the centroid.
    spectSkewness(startWin:finishWin) = spectralSkewness(ps,F*1000);
    
    %% Spectral kurtosis
    % The spectral kurtosis measures the flatness, or non-Gaussianity, of the spectrum around its centroid
    spectKurtosis(startWin:finishWin) = spectralKurtosis(ps,F*1000);
    
    %% Spectral Entropy
    % Because entropy is a measure of disorder, regions of voiced speech have lower entropy compared to regions of unvoiced speech.
    spectEntropy(startWin:finishWin) = spectralEntropy(ps,F*1000);
    
    %% Spectral Flatness
    % Spectral flatness is an indication of the peakiness of the spectrum. 
    % A higher spectral flatness indicates noise, while a lower spectral flatness indicates tonality.
    spectFlatness(startWin:finishWin) = spectralFlatness(ps,F*1000);
    
    %% Spectral Crest
    % Spectral crest is an indication of the peakiness of the spectrum. 
    % A higher spectral crest indicates more tonality, while a lower spectral crest indicates more noise.
    spectCrest(startWin:finishWin) = spectralCrest(ps,F*1000);
    
    %% Spectral Slope
    % Spectral slope measures the amount of decrease of the spectrum:
    % Spectral slope is most pronounced when the energy in the lower formants is much greater than the energy in the higher formants.
    spectSlope(startWin:finishWin) = spectralSlope(ps,F*1000);
    
    %% Spectral Decrease
    % Spectral decrease represents the amount of decrease of the spectrum, 
    % while emphasizing the slopes of the lower frequencies
    spectDecrease(startWin:finishWin) = spectralSlope(ps,F*1000);
    
    %% More spectral features
    % Median frequency
    medianFreq(startWin:finishWin) = medfreq(ps,Fs);
    [bw_,flo_,fhi_,~] = obw(ps,Fs);
    bw(startWin:finishWin) = bw_;
    flo(startWin:finishWin) = flo_;
    fhi(startWin:finishWin) = fhi_;
    %% Additional features to add later
    % Spectral entropy
    % Positive and negative amplitude sum, positive and negative
    % amplitude difference, positive and negative amplitude ratio
    % Power ratio at X Hz (50, 100, 200, 400, 800)

    % Also consider playing around with window size
    % Maybe extract pitch with a small frequency resolution and do the 
    % rest with small time resolution
        %% Plot
    if PlotSnore
        %Plot the spectrogram
        figure(52)
        set(gcf, 'Position', [53 42 1227 605])
        ax1=subplot(7,4,1:8);
        imagesc(T, F, P, clim);
        colormap(flipud(gray))

        set(gca,'YDir','normal', 'ylim',[0 4])
    %     xlabel("Time (s)");
        ylabel("Frequency (kHz)")
        set(gca,'Box','On','FontSize',8,'XTick',[])

        %% overlay processed frequency data   
    %     figure, scatter (f0_cep(similarIdx), f0_xcorr(similarIdx))    
        hold on
%         [f] = PitchSpectralHps (X, Fs);
%         [f] = PitchSpectralAcf (X, Fs);
%         f(HNRcorr < 5) = nan;
%         f0_xcorr_ = f0_xcorr;
%         f0_xcorr_(hr < 0.5) = nan;
        plot(T,f0_xcorr_/1000,'go','MarkerSize',1)
%         plot(T,f/1000,'ro','MarkerSize',1)
%         plot(T,f0_xcorrPSD/1000,'mo','MarkerSize',1)

    %     fmtime = linspace(T(1),T(end),length(formant1(startWin:finishWin)));
    %     plot(fmtime,formant1(startWin:finishWin,1)./1000, 'yo', 'MarkerSize',1)
    %     plot(fmtime,formant2(startWin:finishWin,1)./1000, 'ro', 'MarkerSize',1)
    %     plot(fmtime,formant3(startWin:finishWin,1)./1000, 'mo', 'MarkerSize',1)
        hold off 

        %% plot flow
        ax2=subplot(7,4,9:12);
        plot(timeFlow,Flow10);
        xlabel("Time (s)");
        ylabel("Flow (L/s)")
        set(gca,'Box','On','FontSize',8)
        
        %% plot harmonic ratio
        ax3=subplot(7,4,13:16);
        plot(T',hr,'g'); hold on
%         plot(T',HNRcorr,'r');
%         plot(T',HNRcorrPSD,'m');
%         plot(T',HarmPow,'m');    
        ylim([0 10])
        xlabel("Time (s)");
        ylabel("Harmonic ratio")
        set(gca,'Box','On','FontSize',8)

        linkaxes([ax1,ax2,ax3],'x');

        %% use ginput to mark two spots for closer look
        disp('Script pause. Zoom in to selection area'); pause; % to zoom in
        while 1
            tIn = ginput(2);
            [~, sampStart1] = min(abs(T - tIn(1)));
            [~, sampStart2] = min(abs(T - tIn(2)));

            %% mark selections on spectrogram
            axes(ax1);
            hold on

            ylims = get(gca,'ylim');
            plot([T(sampStart1) T(sampStart1)], [ylims(1) ylims(2)],...
                'Color',[0.75 0 0], 'LineStyle','-')

            ylims = get(gca,'ylim');
            plot([T(sampStart2) T(sampStart2)], [ylims(1) ylims(2)],...
                'Color','m', 'LineStyle','-')   

            hold off

            %% Additional plots
            % PSD plots at select time points (not automated)
            subplot(7,4,17)
            plot(F,P(:,sampStart1),...
                'Color',[0.75 0 0], 'LineStyle','-')
            hold on
            set(gca, 'xlim',[0 1])
            ylims = get(gca,'ylim');
            plot([f0_xcorr(sampStart1)/1000 f0_xcorr(sampStart1)/1000], [ylims(1) ylims(2)],...
                'Color',[0.75 0 0], 'LineStyle','-')
            hold off
            ylabel('Power (dB/Hz)')
            xlabel('Frequency (kHz)')

            subplot(7,4,19)
            plot(F,P(:,sampStart2),...
                'Color','m', 'LineStyle','-')
            hold on
            set(gca, 'xlim',[0 1])
            ylims = get(gca,'ylim');
            plot([f0_xcorr(sampStart2)/1000 f0_xcorr(sampStart2)/1000], [ylims(1) ylims(2)],...
                'Color','m', 'LineStyle','-')
            hold off
            ylabel('Power (dB/Hz)')
            xlabel('Frequency (kHz)')

            % Cepstral / autocorrelation plots
            subplot(7,4,18)
            signal = buffered(:,sampStart1);
%             [c, y] = spCepstrum(signal, Fs, 'hamming', 0);
%             plot(ms2:ms50, c(ms2:ms50),'Color',[0.75 0 0], 'LineStyle','-');
%             ylabel('Power (dB/s)')
%             xlabel('Quefrency (ms)')

            [r] = spCorrelum(signal, Fs, Fs/20, 0);
            r = r(floor(length(r)/2):end);
            plot((1:ms20)./Fs, r(1:ms20),'Color',[0.75 0 0], 'LineStyle','-');
            ylabel('Coeff')
            xlabel('time (ms)')
%             [Rxx_T0_,idx]=max(r(ms2:ms50)); % Rxx_T0coeff is very similar to below

            subplot(7,4,20)
            signal = buffered(:,sampStart2);
%             [c, y] = spCepstrum(signal, Fs, 'hamming', 0);
%             plot(ms2:ms50, c(ms2:ms50),'Color','m', 'LineStyle','-');
%             ylabel('Power (dB/s)')
%             xlabel('Quefrency (ms)')
            [r] = spCorrelum(signal, Fs, Fs/20, 0);
            r = r(floor(length(r)/2):end);
            plot((1:ms20)./Fs, r(1:ms20),'Color','m', 'LineStyle','-');
            ylabel('Coeff')
            xlabel('time (ms)')

            % Formant plot
            signal = buffered(:, sampStart1);
            [~,h,f_] = spLpc(signal, Fs, [], 0);
            subplot(7,4,21)
            plot(f_/1000,20*log10(abs(h)+eps),...
                'Color',[0.75 0 0], 'LineStyle','-')
            set(gca, 'xlim', [0 1])
            ylabel('Gain (dB)')
            xlabel('Frequency (kHz)')

            signal = buffered(:, sampStart2);
            [~,h,f_] = spLpc(signal, Fs, [], 0);
            subplot(7,4,23)
            plot(f_/1000,20*log10(abs(h)+eps),...
                'Color','m', 'LineStyle','-')
            set(gca, 'xlim', [0 1])
            ylabel('Gain (dB)')
            xlabel('Frequency (kHz)')

            % Bandpower plots
            signal = buffered(:, sampStart1);
            pband1 = bandpower(signal,Fs,[50 150]);
            pband2 = bandpower(signal,Fs,[150 300]);
            pband3 = bandpower(signal,Fs,[300 450]);
            pband4 = bandpower(signal,Fs,[450 600]);
            ptottemp = bandpower(signal,Fs,[0 2000]);
            relpband = [pband1/ptottemp;pband2/ptottemp;pband3/ptottemp;pband4/ptottemp];
            subplot(7,4,22)
            bar(1:4,relpband,'r')
            set(gca, 'ylim', [0 1])
            ylabel('% total power')
            xlabel('Band #')

            signal = buffered(:, sampStart2);
            pband1 = bandpower(signal,Fs,[50 150]);
            pband2 = bandpower(signal,Fs,[150 300]);
            pband3 = bandpower(signal,Fs,[300 450]);
            pband4 = bandpower(signal,Fs,[450 600]);
            ptottemp = bandpower(signal,Fs,[0 2000]);
            relpband = [pband1/ptottemp;pband2/ptottemp;pband3/ptottemp;pband4/ptottemp];
            subplot(7,4,24)
            bar(1:4,relpband,'m')
            set(gca, 'ylim', [0 1])
            ylabel('% total power')
            xlabel('Band #')
            
            % MFCC plots
            subplot(7,4,25)
            bar(1:5,coeffsAvg(sampStart1,2:6),'r')
%             set(gca, 'ylim', [0 1])
            ylabel('MFCC')
            xlabel('Coeff #')
            
            subplot(7,4,27)
            bar(1:5,coeffsAvg(sampStart2,2:6),'m')
%             set(gca, 'ylim', [0 1])
            ylabel('MFCC')
            xlabel('Coeff #')
            
            choice = menu('Graph another selection?','Yes','No');
            if choice==2 || choice==0
                break;
            end

            % delete lines
            delete(ax1.Children(1:2))
        end
    end
end
%% Store windowed features into structure
%% Store features in structured array
SnoreStruct = struct();
SnoreStruct.SnoreTime = timeWin;
SnoreStruct.SnoreTimeSm = timeWinSm;
% SnoreStruct.FilteredSnd = FilteredSnd;
% SnoreStruct.timeDS = timeDS;
% SnoreStruct.SndIntensity = SndIntensityArray;
% SnoreStruct.WinVar50 = WinVar50;
SnoreStruct.WinVar = WinVarArray;
SnoreStruct.ZCr = ZCrArray;
SnoreStruct.formant1 = formant1;
SnoreStruct.formant2 = formant2;
SnoreStruct.formant3 = formant3;
SnoreStruct.formant1_n = formant1_n;
SnoreStruct.formant2_n = formant2_n;
SnoreStruct.formant3_n = formant3_n;
SnoreStruct.f0 = f0;
SnoreStruct.p0to125 = p0to125;
SnoreStruct.p125to250 = p125to250;
SnoreStruct.p250to500 = p250to500;
SnoreStruct.p500to1k = p500to1k;
SnoreStruct.p1kto1k5 = p1kto1k5;
SnoreStruct.p1k5to2k = p1k5to2k;
SnoreStruct.p2kto2k5 = p2kto2k5;
SnoreStruct.p2k5to3k = p2k5to3k;
SnoreStruct.ptot = ptot;
SnoreStruct.pkpower = pkpower;
SnoreStruct.pk99power = pk99power;
SnoreStruct.centroid = centroid;
SnoreStruct.HNR = HNR;
SnoreStruct.HNR2 = HNR2;
SnoreStruct.HarmonicPower = HarmonicPower;
SnoreStruct.spectFlux = spectFlux;
SnoreStruct.spectRollOff = spectRollOff;
SnoreStruct.spectSpread = spectSpread;
SnoreStruct.spectSkewness = spectSkewness;
SnoreStruct.spectKurtosis = spectKurtosis;
SnoreStruct.spectFlatness = spectFlatness;
SnoreStruct.spectSlope = spectSlope;
SnoreStruct.spectDecrease = spectDecrease;
SnoreStruct.spectEntropy = spectEntropy;
SnoreStruct.spectCrest = spectCrest;
SnoreStruct.medianFreq = medianFreq;
SnoreStruct.bw = bw;
SnoreStruct.flo = flo;
SnoreStruct.fhi = fhi;
SnoreStruct.ZCrTransform = ZCrTransform;
SnoreStruct.coeffsAvg = coeffsAvg;
SnoreStruct.delAvg = delAvg;
SnoreStruct.ddelAvg = ddelAvg;
SnoreStruct.coeffs = coeffs;
SnoreStruct.del = del;
SnoreStruct.ddel = ddel;

%% Spectrogram
SnoreStruct.Spectrogram.Power = SpectP;
SnoreStruct.Spectrogram.Freq = F;
SnoreStruct.Spectrogram.Time = timeWin;
% SnoreStruct.Spectrogram.ps = ps;
SnoreStruct.Spectrogram.lpcenv = lpcenv;
SnoreStruct.Spectrogram.lpcenv_n = lpcenv_n;


%% Interp windowed features
samplepoints = 1:length(timeWin);
querypoints = linspace(1, length(timeWin), length(Flow));

timeWinInterp = interp1(samplepoints', timeWin, querypoints', 'previous');
WinVarInterp = interp1(samplepoints', WinVarArray, querypoints', 'previous');
ZCrInterp = interp1(samplepoints', ZCrArray, querypoints', 'previous');
formant1Interp = interp1(samplepoints', formant1, querypoints', 'previous');
formant2Interp = interp1(samplepoints', formant2, querypoints', 'previous');
formant3Interp = interp1(samplepoints', formant3, querypoints', 'previous');
formant1nInterp = interp1(samplepoints', formant1_n, querypoints', 'previous');
formant2nInterp = interp1(samplepoints', formant2_n, querypoints', 'previous');
formant3nInterp = interp1(samplepoints', formant3_n, querypoints', 'previous');
p0to125Interp = interp1(samplepoints', p0to125, querypoints', 'previous');
p125to250Interp = interp1(samplepoints', p125to250, querypoints', 'previous');
p250to500Interp = interp1(samplepoints', p250to500, querypoints', 'previous');
p500to1kInterp = interp1(samplepoints', p500to1k, querypoints', 'previous');
p1kto1k5Interp = interp1(samplepoints', p1kto1k5, querypoints', 'previous');
p1k5to2kInterp = interp1(samplepoints', p1k5to2k, querypoints', 'previous');
p2kto2k5Interp = interp1(samplepoints', p2kto2k5, querypoints', 'previous');
p2k5to3kInterp = interp1(samplepoints', p2k5to3k, querypoints', 'previous');
ptotInterp = interp1(samplepoints', ptot, querypoints', 'previous');
pkpowerInterp = interp1(samplepoints', pkpower, querypoints', 'previous');
pk99powerInterp = interp1(samplepoints', pk99power, querypoints', 'previous');
centroidInterp = interp1(samplepoints', centroid, querypoints', 'previous');
f0Interp = interp1(samplepoints', f0, querypoints', 'previous');
HNRInterp = interp1(samplepoints', HNR, querypoints', 'previous');
HNR2Interp = interp1(samplepoints', HNR2, querypoints', 'previous');
HarmonicPowInterp = interp1(samplepoints', HarmonicPower, querypoints', 'previous');
spectFluxInterp = interp1(samplepoints', spectFlux, querypoints', 'previous');
spectRollOffInterp = interp1(samplepoints', spectRollOff, querypoints', 'previous');
spectSkewnessInterp = interp1(samplepoints', spectSkewness, querypoints', 'previous');
spectKurtosisInterp = interp1(samplepoints', spectKurtosis, querypoints', 'previous');
spectFlatnessInterp = interp1(samplepoints', spectFlatness, querypoints', 'previous');
spectSlopeInterp = interp1(samplepoints', spectSlope, querypoints', 'previous');
spectDecreaseInterp = interp1(samplepoints', spectDecrease, querypoints', 'previous');
spectEntropyInterp = interp1(samplepoints', spectEntropy, querypoints', 'previous');
spectCrestInterp = interp1(samplepoints', spectCrest, querypoints', 'previous');
medianFreqInterp = interp1(samplepoints', medianFreq, querypoints', 'previous');
bwInterp = interp1(samplepoints', bw, querypoints', 'previous');
floInterp = interp1(samplepoints', flo, querypoints', 'previous');
fhiInterp = interp1(samplepoints', fhi, querypoints', 'previous');
samplepointsMFCC = 1:length(coeffs);
coeffsInterp = interp1(samplepointsMFCC', coeffs', querypoints', 'previous');
delInterp = interp1(samplepointsMFCC', del', querypoints', 'previous');
ddelInterp = interp1(samplepointsMFCC', ddel', querypoints', 'previous');

%% Store interpolated features in structured array
SnoreInterpStruct = struct();
SnoreInterpStruct.SnoreTime = timeWinInterp;
% SnoreInterpStruct.FilteredSnd = SndData;
% SnoreStruct.timeDS = timeDS;
% SnoreStruct.SndIntensity = SndIntensityArray;
% SnoreStruct.WinVar50 = WinVarInterp50;
SnoreInterpStruct.WinVar = WinVarInterp;
SnoreInterpStruct.ZCr = ZCrInterp;
SnoreInterpStruct.formant1 = formant1Interp;
SnoreInterpStruct.formant2 = formant2Interp;
SnoreInterpStruct.formant3 = formant3Interp;
SnoreInterpStruct.formant1n = formant1nInterp;
SnoreInterpStruct.formant2n = formant2nInterp;
SnoreInterpStruct.formant3n = formant3nInterp;
SnoreInterpStruct.f0 = f0Interp;
SnoreInterpStruct.p0to125 = p0to125Interp;
SnoreInterpStruct.p125to250 = p125to250Interp;
SnoreInterpStruct.p250to500 = p250to500Interp;
SnoreInterpStruct.p500to1k = p500to1kInterp;
SnoreInterpStruct.p1kto1k5 = p1kto1k5Interp;
SnoreInterpStruct.p1k5to2k = p1k5to2kInterp;
SnoreInterpStruct.p2kto2k5 = p2kto2k5Interp;
SnoreInterpStruct.p2k5to3k = p2k5to3kInterp;
SnoreInterpStruct.ptot = ptotInterp;
SnoreInterpStruct.pkpower = pkpowerInterp;
SnoreInterpStruct.pk99power = pk99powerInterp;
SnoreInterpStruct.centroid = centroidInterp;
SnoreInterpStruct.HNR = HNRInterp;
SnoreInterpStruct.HNR2 = HNR2Interp;
SnoreInterpStruct.HarmonicPower = HarmonicPowInterp;
SnoreInterpStruct.spectFlux = spectFluxInterp;
SnoreInterpStruct.spectRollOff = spectRollOffInterp;
SnoreInterpStruct.spectSkewness = spectSkewnessInterp;
SnoreInterpStruct.spectKurtosis = spectKurtosisInterp;
SnoreInterpStruct.spectFlatness = spectFlatnessInterp;
SnoreInterpStruct.spectSlope = spectSlopeInterp;
SnoreInterpStruct.spectDecrease = spectDecreaseInterp;
SnoreInterpStruct.spectEntropy = spectEntropyInterp;
SnoreInterpStruct.spectCrest = spectCrestInterp;
SnoreInterpStruct.medianFreq = medianFreqInterp;
SnoreInterpStruct.bw = bwInterp;
SnoreInterpStruct.flo = floInterp;
SnoreInterpStruct.fhigh = fhiInterp;
SnoreInterpStruct.coeffs = coeffsInterp;
SnoreInterpStruct.del = delInterp;
SnoreInterpStruct.ddel = ddelInterp;

end