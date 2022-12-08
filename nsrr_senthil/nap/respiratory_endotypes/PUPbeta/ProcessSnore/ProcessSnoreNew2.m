function [SnoreStruct] = ProcessSnoreNew2(Snore,Flow,Time,FsFlow,Fs,PlotSnore,windowselect)

disp('Processing Snore Data')
%% Initilize settings
% Fs = 1/Snore.interval;
% SndData = Snore.values;

% Data
% Fs = Channels_Fs(strcmp(ChannelsList,'SnoreRaw'));
% FsFlow = Channels_Fs(strcmp(ChannelsList,'Flow'));
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

SnoreCellArray = cell(2, ceil(length(SndData)/(SegLen*Fs)));
SpectrogramCellArray = cell(1, ceil(length(SndData)/(SegLen*Fs)));
lpcCellArray = cell(2, ceil(length(SndData)/(SegLen*Fs)));

%% Process sounds
% isolate 10 min segments
% possibly temporary but might be faster to work in sections
% in which case we loop everything below

% select windows to analyze
if exist('windowselect','var') && ~isempty(windowselect)
    winnumrange = windowselect;
else
    winnumrange = 1:ceil(length(SndData)/(SegLen*Fs));
end

for ii = winnumrange
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
    FilterBandwidth = [5 2500];
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
%     if ii == 1
%         startWin = 1;
%     else
%         startWin = finishWinLastWindow + 1;
%     end
%     
%     finishWin = startWin + size(buffered,2)-1; %floor((ii*SegLen*FsBuffer)-1);
%     finishWinLastWindow = finishWin;
%     
%     if finishWin > length(WinVarArray)
%         finishWin = startWin + size(buffered,2)-1;
%     end
    
    %%%% Adapted from bufferSignal code
    % Determine the number of columns of the STFT output (i.e., the S output)
    ncol = fix((length(FilteredSnd)-noverlap)/(winsize-noverlap));
    coloffsets = (0:(ncol-1))*(winsize-noverlap);
%     rowindices = (startWin:finishWin)';

    % return time vector whose elements are centered in each segment
    timeWin = (StartTime + (ii-1)*SegLen + (coloffsets+(winsize/2)')/Fs)';

    %% Buffer data in smaller windows
    % this is for features that need smaller windows (e.g. making a ZCr
    % transformed signal, from which to calculate sample entropy)
       
    % Smaller window buffered signal
    [bufferedSm, ~] = bufferSignal(FilteredSnd,length(FilteredSnd),floor(Fs*timeResSm),floor(Fs*percoverlapSm*timeResSm),Fs);
    
    % start condition
%     if ii == 1
%         startWinSm = 1;
%     else
%         startWinSm = finishWinLastWindowSm + 1;
%     end
%     
%     finishWinSm = startWinSm + size(bufferedSm,2)-1; %floor((ii*SegLen*FsBuffer)-1);
%     finishWinLastWindowSm = finishWinSm;
%     
%     if finishWinSm > length(ZCrTransform)
%         finishWinSm = startWinSm + size(bufferedSm,2)-1;
%     end
%     
    %%%% Adapted from bufferSignal code
    % Determine the number of columns of the STFT output (i.e., the S output)
    ncol = fix((length(FilteredSnd)-noverlapSm)/(winsizeSm-noverlapSm));
    coloffsetsSm = (0:(ncol-1))*(winsizeSm-noverlapSm);
%     rowindices = (startWin:finishWin)';

    % return time vector whose elements are centered in each segment
    timeWinSm =  (StartTime + (ii-1)*SegLen + (coloffsetsSm+(winsizeSm/2)')/Fs)';
    
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
    WinVar = nanvar(buffered,1)';
    MSs = nanmean(buffered,1);
    ZSs = buffered - MSs;
    ZCr = (sum(abs(sign(ZSs(1:end-1,:))-sign(ZSs(2:end,:))))/(2*winsize))';
    SnoreDB = 10*log10(WinVar./(0.00002^2));
    
    % Band power (absolute)
    % Peak power is computed below where the PSD is computed, this can be
    % used to normalize
    % tested this and its the exact same answer if I calculate from ps
    % (saved in SpectP so I can do these post process too)
    % e.g Fdiff = diff(F); p0to125 = sum(ps(f0:f125,:),1)*Fdiff(1);
    
    if Fs > 5000
        ptot = bandpower(buffered,Fs,[0 Fs/2])';
        p2kto3k_ = bandpower(buffered,Fs,[2000 3000])';
        p2k5to3k_ = bandpower(buffered,Fs,[2500 3000])';
        p3kto5k_ = bandpower(buffered,Fs,[3000 Fs/2])';
        p3kto5k = log10(p3kto5k_./(0.00002^2))./log10(ptot./(0.00002^2));
    else
        ptot = bandpower(buffered,Fs,[0 2500])';
        p2kto3k_ = nan(size(ptot));
        p2k5to3k_ = nan(size(ptot));
        p3kto5k = nan(size(ptot));
    end
    p2kto3k = log10(p2kto3k_./(0.00002^2))./log10(ptot./(0.00002^2));
    p2k5to3k = log10(p2k5to3k_./(0.00002^2))./log10(ptot./(0.00002^2));
    
    p0to50_ = bandpower(buffered,Fs,[0 50])';
    p0to50 = log10(p0to50_./(0.00002^2))./log10(ptot./(0.00002^2));
    p0to125_ = bandpower(buffered,Fs,[0 125])';
    p0to125 = log10(p0to125_./(0.00002^2))./log10(ptot./(0.00002^2));
    p0to250_ = bandpower(buffered,Fs,[0 250])';
    p0to250 = log10(p0to250_./(0.00002^2))./log10(ptot./(0.00002^2));
    p0to500_ = bandpower(buffered,Fs,[0 500])';
    p0to500 = log10(p0to500_./(0.00002^2))./log10(ptot./(0.00002^2));
    p125to250_ = bandpower(buffered,Fs,[125 250])';
    p125to250 = log10(p125to250_./(0.00002^2))./log10(ptot./(0.00002^2));
    p125to500_ = bandpower(buffered,Fs,[125 500])';
    p125to500 = log10(p125to500_./(0.00002^2))./log10(ptot./(0.00002^2));
    p250to500_ = bandpower(buffered,Fs,[250 500])';
    p250to500 = log10(p250to500_./(0.00002^2))./log10(ptot./(0.00002^2));
    p500to1k_ = bandpower(buffered,Fs,[500 1000])';
    p500to1k = log10(p500to1k_./(0.00002^2))./log10(ptot./(0.00002^2));
    p1kto1k5_ = bandpower(buffered,Fs,[1000 1500])';
    p1kto1k5 = log10(p1kto1k5_./(0.00002^2))./log10(ptot./(0.00002^2));
    p1k5to2k_ = bandpower(buffered,Fs,[1500 2000])';
    p1k5to2k = log10(p1k5to2k_./(0.00002^2))./log10(ptot./(0.00002^2));
    p2kto2k5_ = bandpower(buffered,Fs,[2000 2500])';
    p2kto2k5 = log10(p2kto2k5_./(0.00002^2))./log10(ptot./(0.00002^2));
    
    % compute significant log power ratios
    logratio_125_250_0_125 = 10*log10(p125to250./(0.00002^2))./...
        (10*log10(p0to125./(0.00002^2)));
    logratio_250_500_0_250 = 10*log10(p250to500)./...
        (10*log10((p0to125 + p125to250)./(0.00002^2)));
    logratio_2k_2k5_0_500 = 10*log10(p2kto2k5./0.00002^2)./...
         (10*log10((p0to125 + p125to250 + p250to500)./0.00002^2));
    logratio_2k_2k5_0_250 = 10*log10(p2kto2k5./0.00002^2)./...
        (10*log10((p0to125 + p125to250)./0.00002^2));
    logratio_2k_2k5_250_500 = 10*log10(p2kto2k5./0.00002^2)./...
        (10*log10((p125to250 + p250to500)./0.00002^2));

    %% Spectral centroid
    centroid = getCentroid(buffered,Fs,winsize)';
    
   %% ZCr tranformed signal from which to calculate entropy later
    % Generate ZCr transformed signal
    MSs2 = nanmean(bufferedSm,1);
    ZSs2 = bufferedSm - MSs2;
    ZCr2 = sum(abs(sign(ZSs2(1:end-1,:))-sign(ZSs2(2:end,:))))/(2*(Fs*0.025));
    ZCrTransform = ZCr2';
    
    %% test pwelch to get array size
    pwelchxx = pwelch(buffered(:,1));
    logpwelchxx = 10*log10(pwelchxx);
    L = 15;
    w = bartlett(L);
    pwelchwinsmooth = conv(w,logpwelchxx,'full');  
    psdsmoothtemp = pwelchwinsmooth(L+1:end-L);
    
    clear pwelchxx logpwelchxx L w pwelchwinsmooth
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
    formant1 = nan(size(buffered,2),1);
    formant2 = nan(size(buffered,2),1);
    formant3 = nan(size(buffered,2),1);
    formant1_n = nan(size(buffered,2),1);
    formant2_n = nan(size(buffered,2),1);
    formant3_n = nan(size(buffered,2),1);
    psdsmooth = nan(size(buffered,2),length(psdsmoothtemp));
    
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
        ncoef = 18;
%         ncoef = 2 + round(Fs / 714); % rule of thumb for human speech is fs/1000, but I reduced to fs/250
        [ffreq,h,f_] = spLpc(signal, Fs, ncoef, 0);
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
%         ncoef = 12;
        ncoef = 2 + round(Fs / 1000); % rule of thumb for human speech is fs/1000, but I reduced to fs/250
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
        
        %% Smoothed PSD using pwelch
         pwelchxx = pwelch(signal);
         logpwelchxx = 10*log10(pwelchxx);
         L = 15;
         w = bartlett(L);
         pwelchwinsmooth = conv(w,logpwelchxx,'full');
         psdsmooth(winnum,:) = pwelchwinsmooth(L+1:end-L);
%          figure,plot(pwelchwinsmooth(L+1:end-L)), hold on, plot(10*log10(pwelchxx(:,1)))
        
    end 
%%
    % Only use pitch data similar between both methods
%     similarIdx = abs(f0_cep - f0_xcorr) < 999; sum(similarIdx); % Value of this parameter not tested
%     f0_xcorr(~similarIdx) = nan;
%     HNRcorr(~similarIdx) = nan; % removes data with low harmonic to noise ratio - not using this right now
%     f0_xcorr(HNRcorr<0.5 | SNR < 2.5) = nan;
    f0 = f0_xcorr;
    HNR = HNRcorr;
    HNR2 = harmonicRatio(buffered,Fs,'Window',hamming(winsize),'OverlapLength',noverlap)';
    HarmPow = smooth(HarmPow,4);
    HarmonicPower = HarmPow;
    lpcenv = lpcenvtemp;
    lpcenv_n = lpcenvtemp_n;
    
    % Spectrogram
    figure('Visible','off');    
%     spectrogram(FilteredSnd,winsize,noverlap,nfft,Fs,'yaxis');
    [X,F,T,ps] = spectrogram(FilteredSnd,winsize,noverlap,nfft,Fs,'yaxis');
    T = timeWin; % overwrite T
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
    
    SpectP = P'; % store spectral power
    SpectPabs = ps';
    
    % Peak power (to normalize bandpowers)
    pk99power = prctile(ps,99,1)';
    pkpower = prctile(ps,100,1)';
    
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

    coeffsArray = coeffs';
    delArray = del';
    ddelArray = ddel';


    % buffer signal and average to make same resolution as other
    % features
    FsCoeffs = Fs/(0.01/(1/Fs));

    % deal with issues in array size
    [bufferedCoeffs, ~] = bufferSignal(coeffs(1,:),length(coeffs(1,:)),10,5,FsCoeffs); % sample output 
%     finishWinBuf = startWin + size(bufferedCoeffs,2)-1;

    coeffsAvg = nan(length(timeWin),size(coeffs,1));
    delAvg = nan(length(timeWin),size(coeffs,1));
    ddelAvg = nan(length(timeWin),size(coeffs,1));
    for coefnum = 1:size(coeffs,1)
        [bufferedCoeffs, ~] = bufferSignal(coeffs(coefnum,:),length(coeffs(coefnum,:)),10,5,FsCoeffs); % still a couple frames long
        meanTemp = mean(bufferedCoeffs,1);
        coeffsAvg(:, coefnum) = meanTemp(1:length(coeffsAvg));
        
        [bufferedDel, ~] = bufferSignal(del(coefnum,:),length(del(coefnum,:)),10,5,FsCoeffs); % still a couple frames long
        meanTemp = mean(bufferedDel,1);
        delAvg(:, coefnum) = meanTemp(1:length(coeffsAvg));

        [bufferedDdel, ~] = bufferSignal(ddel(coefnum,:),length(ddel(coefnum,:)),10,5,FsCoeffs); % still a couple frames long
        meanTemp = mean(bufferedDdel,1);
        ddelAvg(:,coefnum) = meanTemp(1:length(coeffsAvg));
    end
    
    if size(coeffsAvg,1) > length(timeWin)
        coeffsAvg = coeffsAvg(1:length(timeWin),:);
        delAvg = coeffsAvg(1:length(timeWin),:);
        ddelAvg = coeffsAvg(1:length(timeWin),:);
    end
    %% MFCC 
    % 5k
    [coeffs,~,~] = melfcc(FilteredSnd, Fs, 'wintime',0.020,'hoptime',0.01,'maxfreq',5000); 
    del = deltas(coeffs);
    % Double deltas are deltas applied twice with a shorter window
    ddel = deltas(deltas(coeffs,5),5);

    coeffsArray5k = coeffs';
    delArray5k = del';
    ddelArray5k = ddel';


    % buffer signal and average to make same resolution as other
    % features
    FsCoeffs = Fs/(0.01/(1/Fs));

    % deal with issues in array size
    [bufferedCoeffs, ~] = bufferSignal(coeffs(1,:),length(coeffs(1,:)),10,5,FsCoeffs); % sample output 
%     finishWinBuf = startWin + size(bufferedCoeffs,2)-1;

    % NEED TO FIX THIS
%     for coefnum = 1:size(coeffs,1)
%         [bufferedCoeffs, ~] = bufferSignal(coeffs(coefnum,:),length(coeffs(coefnum,:)),10,5,FsCoeffs); % still a couple frames long
%         coeffsAvg5k(:, coefnum) = mean(bufferedCoeffs,1);
% 
%         [bufferedDel, ~] = bufferSignal(del(coefnum,:),length(del(coefnum,:)),10,5,FsCoeffs); % still a couple frames long
%         delAvg5k(:, coefnum) = mean(bufferedDel,1);
% 
%         [bufferedDdel, ~] = bufferSignal(ddel(coefnum,:),length(ddel(coefnum,:)),10,5,FsCoeffs); % still a couple frames long
%         ddelAvg5k(:,coefnum) = mean(bufferedDdel,1);
%     end

    %% Spectral flux
    % the rate at which spectrum is changing between two consecutive
    % frames. Calculated as the squared difference between normalized
    % spectra of two consecutive frames.
    % algorithm from: https://www.audiocontentanalysis.org/code/audio-features
    % this has a little bug caused by analyzing in segements
    spectFlux = FeatureSpectralFlux (X, Fs)';

    %% Spectral roll off
    % Measures "skewness" of the spectral shape. Concretely, it
    % represenents the frequency at which N% of the power is
    kappa = 0.85;
    spectRollOff = FeatureSpectralRolloff(X, Fs, kappa)';

    %% Spectral spread
    % The standard deviation around the spectral centroid 
    spectSpread = spectralSpread(ps,F*1000);
    
    %% Spectral skewness 
    % The spectral skewness measures symmetry around the centroid.
    spectSkewness = spectralSkewness(ps,F*1000);
    
    %% Spectral kurtosis
    % The spectral kurtosis measures the flatness, or non-Gaussianity, of the spectrum around its centroid
    spectKurtosis = spectralKurtosis(ps,F*1000);
    
    %% Spectral Entropy
    % Because entropy is a measure of disorder, regions of voiced speech have lower entropy compared to regions of unvoiced speech.
    spectEntropy = spectralEntropy(ps,F*1000);
    
    %% Spectral Flatness
    % Spectral flatness is an indication of the peakiness of the spectrum. 
    % A higher spectral flatness indicates noise, while a lower spectral flatness indicates tonality.
    spectFlatness = spectralFlatness(ps,F*1000);
    
    %% Spectral Crest
    % Spectral crest is an indication of the peakiness of the spectrum. 
    % A higher spectral crest indicates more tonality, while a lower spectral crest indicates more noise.
    spectCrest = spectralCrest(ps,F*1000);
    
    %% Spectral Slope
    % Spectral slope measures the amount of decrease of the spectrum:
    % Spectral slope is most pronounced when the energy in the lower formants is much greater than the energy in the higher formants.
    spectSlope = spectralSlope(ps,F*1000);
    
    %% Spectral Decrease
    % Spectral decrease represents the amount of decrease of the spectrum, 
    % while emphasizing the slopes of the lower frequencies
    spectDecrease = spectralSlope(ps,F*1000);
    
    %% More spectral features
    % Median frequency
    medianFreq = medfreq(ps,Fs)';
    [bw_,flo_,fhi_,~] = obw(ps,Fs);
    bw = bw_';
    flo = flo_';
    fhi = fhi_';
    
    %% Additional features to add later
    % Spectral entropy
    % Positive and negative amplitude sum, positive and negative
    % amplitude difference, positive and negative amplitude ratio
    % Power ratio at X Hz (50, 100, 200, 400, 800)

    % Also consider playing around with window size
    % Maybe extract pitch with a small frequency resolution and do the 
    % rest with small time resolution
    
    %% Put into a table and then into cell array
    SnoreTableWin = table(timeWin,SnoreDB,WinVar,ZCr,formant1,formant2,formant3,...
        formant1_n,formant2_n,formant3_n,f0,p0to125,p0to50,p0to250,p0to500,p125to250,...
        p125to500,p250to500,p500to1k,p1kto1k5,p1k5to2k,p2kto2k5,p2kto3k,p2k5to3k,p3kto5k,...
        logratio_125_250_0_125,logratio_250_500_0_250,logratio_2k_2k5_0_500,...
        logratio_2k_2k5_0_250,logratio_2k_2k5_250_500,...
        ptot,pkpower,pk99power,centroid,HNR,HNR2,HarmonicPower,spectFlux,...
        spectRollOff,spectSpread,spectSkewness,spectKurtosis,spectFlatness,...
        spectSlope,spectDecrease,spectEntropy,spectCrest,medianFreq,...
        bw,flo,fhi,coeffsAvg,delAvg,ddelAvg...
    );

    SnoreTableWinSm = table(timeWinSm,ZCrTransform,coeffsArray,delArray,ddelArray,...
        coeffsArray5k,delArray5k,ddelArray5k);
    
    SnoreCellArray{1,ii} = SnoreTableWin; % store features
    SnoreCellArray{2,ii} = SnoreTableWinSm; % store small window features
    SpectrogramCellArray{1,ii} = SpectP; % store spectrogram
    lpcCellArray{1,ii} = lpcenv; % store lpc
    lpcCellArray{2,ii} = lpcenv_n;
    lpcCellArray{3,ii} = psdsmooth;
    
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
%         plot(T,f0/1000,'go','MarkerSize',1) % UNCOMMENT TO SEE PITCH
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
%         plot(T',hr,'g'); hold on
        plot(T',HNRcorr,'r');
%         plot(T',HNRcorrPSD,'m');
%         plot(T',10*HarmPow,'m');    
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
SnoreStruct.SnoreCellArray = SnoreCellArray;
SnoreStruct.SpectrogramCellArray = SpectrogramCellArray;
SnoreStruct.lpcCellArray = lpcCellArray;
SnoreStruct.Freq = F;
SnoreStruct.lpcFreq = f_;

end