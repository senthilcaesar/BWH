function [SnoreInterpStruct, SnoreStruct] = ProcessSnore2(Snore,Flow,ChannelsList,Channels_Fs,PlotSnore)

%% Initilize settings
% Fs = 1/Snore.interval;
% SndData = Snore.values;

% Data
Fs = Channels_Fs(strcmp(ChannelsList,'Snore'));
FsFlow = Channels_Fs(strcmp(ChannelsList,'Flow'));
SndData = Snore;
f = 1025; % highest frequency for spectrogram

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



% ALL THESE FEATURES TO BE COMPUTED POST PROCESS
WinVarArray = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
ZCrArray = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
timeWin = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
formant1 = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
formant2 = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
formant3 = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
p0to150 = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
p150to300 = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
p300to450 = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
p450to600 = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
p600to750 = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
p750to900 = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
p0to1k = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
p1kto2k = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
ptot = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
f0 = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
HNR = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
HarmonicPower  = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
centroid = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
timeWin = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
SpectP = nan(floor(length(SndData)*(FsBuffer/Fs)),f);
spectFlux = nan(floor(length(SndData)*(FsBuffer/Fs)),f);
spectRollOff = nan(floor(length(SndData)*(FsBuffer/Fs)),f);

% Small window vars
ZCrTransform = nan(floor(length(SndData)*(FsBufferSm/Fs)),1);
timeWinSm = nan(floor(length(SndData)*(FsBufferSm/Fs)),1);
coeffsArray = nan(floor(length(SndData)*(FsBufferSm/Fs)),ncoefs);
delArray = nan(floor(length(SndData)*(FsBufferSm/Fs)),ncoefs);
ddelArray = nan(floor(length(SndData)*(FsBufferSm/Fs)),ncoefs);

%% Process sounds
% isolate 10 min segments
% possibly temporary but might be faster to work in sections
% in which case we loop everything below
for ii = 1:ceil(length(SndData)/(SegLen*Fs))
    start = (ii-1)*SegLen*Fs + 1;
    finish = ii*SegLen*Fs + (percoverlap)*winsize;
    
    if finish > length(SndData) % end condition
        finish = length(SndData);
    end
    
    SndData10Temp = SndData(start:finish);
    t=(start:finish)'/Fs;
    
    startFlow = (ii-1)*SegLen*FsFlow + 1;
    finishFlow = ii*SegLen*FsFlow;
    
    if finishFlow > length(Flow) % end condition
        finishFlow = length(Flow);
    end
    
    Flow10 = Flow(startFlow:finishFlow);
    timeFlow = (startFlow:finishFlow)'/FsFlow;

    %% Filter
    FilterBandwidth = [5 2000];
    SndData10 = SndData10Temp;
    FilteredSnd = ButterFilter(SndData10,FilterBandwidth(1),FilterBandwidth(2),Fs,5);
%     FilteredSnd = SndData10;
    
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
    timeWin(startWin:finishWin) = (ii-1)*SegLen + (coloffsets+(winsize/2)')/Fs;
    
    %% Buffer data in smaller windows
    % this is for features that need smaller windows (e.g. making a ZCr
    % transformed signal, from which to calculate sample entropy)
       
    % Smaller window buffered signal
    [bufferedSm, ~] = bufferSignal(FilteredSnd,length(FilteredSnd),floor(Fs*timeResSm),floor(Fs*percoverlapSm*timeResSm),Fs);
    
    % start condition
    if ii == 1
        startWinSm = 1;
    else
        startWinSm = finishWinLastWindow + 1;
    end
    
    finishWinSm = startWinSm + size(bufferedSm,2)-1; %floor((ii*SegLen*FsBuffer)-1);
    finishWinLastWindow = finishWinSm;
    
    if finishWinSm > length(ZCrTransform)
        finishWinSm = startWinSm + size(bufferedSm,2)-1;
    end
    
    %%%% Adapted from bufferSignal code
    % Determine the number of columns of the STFT output (i.e., the S output)
    ncol = fix((length(FilteredSnd)-noverlapSm)/(winsizeSm-noverlapSm));
    coloffsetsSm = (0:(ncol-1))*(winsizeSm-noverlapSm);
%     rowindices = (startWin:finishWin)';

    % return time vector whose elements are centered in each segment
    timeWinSm(startWinSm:finishWinSm) = (ii-1)*SegLen + (coloffsetsSm+(winsizeSm/2)')/Fs;
    
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
    
    % Band power
    p0to150(startWin:finishWin) = bandpower(buffered,Fs,[0 150]);
    p150to300(startWin:finishWin) = bandpower(buffered,Fs,[150 300]);
    p300to450(startWin:finishWin) = bandpower(buffered,Fs,[300 450]);
    p450to600(startWin:finishWin) = bandpower(buffered,Fs,[450 600]);
    p600to750(startWin:finishWin) = bandpower(buffered,Fs,[600 750]);
    p750to900(startWin:finishWin) = bandpower(buffered,Fs,[750 900]);
    p0to1k(startWin:finishWin) = bandpower(buffered,Fs,[0 1000]);
    p1kto2k(startWin:finishWin) = bandpower(buffered,Fs,[1000 2000]);
    ptot(startWin:finishWin) = bandpower(buffered,Fs,[0 2000]);
    
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
        [ffreq,~,~] = spLpc(signal, Fs, [], 0);
        formant1(winnum,1) = ffreq(1);
        formant2(winnum,1) = ffreq(2);
        formant3(winnum,1) = ffreq(3);
        
    end 
    
    % Only use pitch data similar between both methods
%     similarIdx = abs(f0_cep - f0_xcorr) < 999; sum(similarIdx); % Value of this parameter not tested
%     f0_xcorr(~similarIdx) = nan;
%     HNRcorr(~similarIdx) = nan; % removes data with low harmonic to noise ratio - not using this right now
%     f0_xcorr(HNRcorr<0.5 | SNR < 2.5) = nan;
    f0(startWin:finishWin) = f0_xcorr;
    HNR(startWin:finishWin) = HNRcorr;
    HarmPow = smooth(HarmPow,4);
    HarmonicPower(startWin:finishWin) = HarmPow;
    % Spectrogram
    figure('Visible','off');    
    spectrogram(FilteredSnd,winsize,noverlap,nfft,Fs,'yaxis');

    % get data to rebuild spectrogram
    ch = get(gca,'ch');
    ax = gca;
    T = get(ch, 'XData') + SegLen*(ii-1);
    F = get(ch, 'YData');
    P = get(ch,'CData');
%     t = linspace(0, length(FilteredSndDS), length(FilteredSndDS));
    clim = ax.CLim;
    
    SpectP(startWin:finishWin,:) = P'; % store spectral power
        %% Plot
    if PlotSnore
        %Plot the spectrogram
        figure(52)
        set(gcf, 'Position', [53 42 1227 605])
        ax1=subplot(7,4,1:8);
        imagesc(T, F, P, clim);
        colormap(flipud(gray))

        set(gca,'YDir','normal', 'ylim',[0 1])
    %     xlabel("Time (s)");
        ylabel("Frequency (kHz)")
        set(gca,'Box','On','FontSize',8,'XTick',[])

        %% overlay processed frequency data   
    %     figure, scatter (f0_cep(similarIdx), f0_xcorr(similarIdx))    
        hold on
        plot(T,f0_xcorr/1000,'go','MarkerSize',1)
%         plot(T,f0_cepMod/1000,'ro','MarkerSize',1)
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
        plot(T',HNRcorr,'g'); hold on
%         plot(T',HNRcepMod/10,'r');
%         plot(T',HNRcorrPSD,'m');
        plot(T',HarmPow,'m');    
        ylim([0 10])
        xlabel("Time (s)");
        ylabel("Harmonic ratio")
        set(gca,'Box','On','FontSize',8)

        linkaxes([ax1,ax2,ax3],'x');
        
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
        [bufferedCoeffs, ~] = bufferSignal(coeffs(1,:),length(coeffs(1,:)),10,5,FsCoeffs); % sample output 
        coeffsAvg = nan(size(bufferedCoeffs));
        delAvg = nan(size(bufferedCoeffs));
        ddelAvg = nan(size(bufferedCoeffs));
        
        for coefnum = 1:size(coeffs,1)
            [bufferedCoeffs, ~] = bufferSignal(coeffs(coefnum,:),length(coeffs(coefnum,:)),10,5,FsCoeffs); % still a couple frames long
            coeffsAvg(coefnum,:) = mean(bufferedCoeffs,1);
            
            [bufferedDel, ~] = bufferSignal(del(coefnum,:),length(del(coefnum,:)),10,5,FsCoeffs); % still a couple frames long
            delAvg(coefnum,:) = mean(bufferedDel,1);
            
            [bufferedDdel, ~] = bufferSignal(ddel(coefnum,:),length(ddel(coefnum,:)),10,5,FsCoeffs); % still a couple frames long
            ddelAvg(coefnum,:) = mean(bufferedDdel,1);
        end
        
        
        %% Spectral flux
        % the rate at which spectrum is changing between two consecutive
        % frames. Calculated as the squared difference between normalized
        % spectra of two consecutive frames.
        % algorithm from: https://www.audiocontentanalysis.org/code/audio-features
        spectFlux = FeatureSpectralFlux (P, Fs);
        
        
        %% Spectral roll off
        % Measures "skewness" of the spectral shape. Concretely, it
        % represenents power below 85%
        kappa = 0.85;
        spectRollOff = FeatureSpectralRolloff (P, Fs, kappa);
        
        %% Additional features to add later
        % Spectral entropy
        % Positive and negative amplitude sum, positive and negative
        % amplitude difference, positive and negative amplitude ratio
        % Power ratio at X Hz (50, 100, 200, 400, 800)
        
        % Also consider playing around with window size
        % Maybe extract pitch with a small frequency resolution and do the 
        % rest with small time resolution
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
            [~,h,f] = spLpc(signal, Fs, [], 0);
            subplot(7,4,21)
            plot(f/1000,20*log10(abs(h)+eps),...
                'Color',[0.75 0 0], 'LineStyle','-')
            set(gca, 'xlim', [0 1])
            ylabel('Gain (dB)')
            xlabel('Frequency (kHz)')

            signal = buffered(:, sampStart2);
            [~,h,f] = spLpc(signal, Fs, [], 0);
            subplot(7,4,23)
            plot(f/1000,20*log10(abs(h)+eps),...
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
            bar(1:5,coeffsAvg(2:6,sampStart1),'r')
%             set(gca, 'ylim', [0 1])
            ylabel('MFCC')
            xlabel('Coeff #')
            
            subplot(7,4,27)
            bar(1:5,coeffs(2:6,sampStart2),'m')
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
% SnoreStruct.FilteredSnd = FilteredSnd;
% SnoreStruct.timeDS = timeDS;
% SnoreStruct.SndIntensity = SndIntensityArray;
% SnoreStruct.WinVar50 = WinVar50;
SnoreStruct.WinVar = WinVarArray;
SnoreStruct.ZCr = ZCrArray;
SnoreStruct.formant1 = formant1;
SnoreStruct.formant2 = formant2;
SnoreStruct.formant3 = formant3;
SnoreStruct.f0 = f0;
SnoreStruct.p0to150 = p0to150;
SnoreStruct.p150to300 = p150to300;
SnoreStruct.p300to450 = p300to450;
SnoreStruct.p450to600 = p450to600;
SnoreStruct.p600to750 = p600to750;
SnoreStruct.p750to900 = p750to900;
SnoreStruct.p0to1k = p0to1k;
SnoreStruct.p1kto2k = p1kto2k;
SnoreStruct.ptot = ptot;
SnoreStruct.centroid = centroid;
SnoreStruct.HNR = HNR;
SnoreStruct.HarmonicPower = HarmonicPower;

%% Spectrogram
SnoreStruct.Spectrogram.Power = SpectP;
SnoreStruct.Spectrogram.Freq = F;
SnoreStruct.Spectrogram.Time = timeWin;

%% Interp windowed features
samplepoints = 1:length(timeWin);
querypoints = linspace(1, length(timeWin), length(Flow));

timeWinInterp = interp1(samplepoints', timeWin, querypoints', 'previous');
WinVarInterp = interp1(samplepoints', WinVarArray, querypoints', 'previous');
ZCrInterp = interp1(samplepoints', ZCrArray, querypoints', 'previous');
formant1Interp = interp1(samplepoints', formant1, querypoints', 'previous');
formant2Interp = interp1(samplepoints', formant2, querypoints', 'previous');
formant3Interp = interp1(samplepoints', formant3, querypoints', 'previous');
p0to150Interp = interp1(samplepoints', p0to150, querypoints', 'previous');
p150to300Interp = interp1(samplepoints', p150to300, querypoints', 'previous');
p300to450Interp = interp1(samplepoints', p300to450, querypoints', 'previous');
p450to600Interp = interp1(samplepoints', p450to600, querypoints', 'previous');
p600to750Interp = interp1(samplepoints', p600to750, querypoints', 'previous');
p750to900Interp = interp1(samplepoints', p750to900, querypoints', 'previous');
p0to1kInterp = interp1(samplepoints', p0to1k, querypoints', 'previous');
p1kto2kInterp = interp1(samplepoints', p1kto2k, querypoints', 'previous');
ptotInterp = interp1(samplepoints', ptot, querypoints', 'previous');
centroidInterp = interp1(samplepoints', centroid, querypoints', 'previous');
f0Interp = interp1(samplepoints', f0, querypoints', 'previous');
HNRInterp = interp1(samplepoints', HNR, querypoints', 'previous');
HarmonicPowInterp = interp1(samplepoints', HarmonicPower, querypoints', 'previous');

%% Store interpolated features in structured array
SnoreInterpStruct = struct();
SnoreInterpStruct.SnoreTime = timeWinInterp;
SnoreInterpStruct.FilteredSnd = FilteredSnd;
% SnoreStruct.timeDS = timeDS;
% SnoreStruct.SndIntensity = SndIntensityArray;
% SnoreStruct.WinVar50 = WinVarInterp50;
SnoreInterpStruct.WinVar = WinVarInterp;
SnoreInterpStruct.ZCr = ZCrInterp;
SnoreInterpStruct.formant1 = formant1Interp;
SnoreInterpStruct.formant2 = formant2Interp;
SnoreInterpStruct.formant3 = formant3Interp;
SnoreInterpStruct.f0 = f0Interp;
SnoreInterpStruct.p0to150 = p0to150Interp;
SnoreInterpStruct.p150to300 = p150to300Interp;
SnoreInterpStruct.p300to450 = p300to450Interp;
SnoreInterpStruct.p450to600 = p450to600Interp;
SnoreInterpStruct.p600to750 = p600to750Interp;
SnoreInterpStruct.p750to900 = p750to900Interp;
SnoreInterpStruct.p0to1k = p0to1kInterp;
SnoreInterpStruct.p1kto2k = p1kto2kInterp;
SnoreInterpStruct.ptot = ptotInterp;
SnoreInterpStruct.centroid = centroidInterp;
SnoreInterpStruct.HNR = HNRInterp;
SnoreInterpStruct.HarmonicPower = HarmonicPowInterp;

end