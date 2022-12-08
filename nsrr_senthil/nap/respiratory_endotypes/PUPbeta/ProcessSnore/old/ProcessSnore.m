function SnoreStruct = ProcessSnore(Snore,Flow,ChannelsList,Channels_Fs)

global F_samp

%% Initilize settings
% Fs = 1/Snore.interval;
% SndData = Snore.values;
Fs = Channels_Fs(strcmp(ChannelsList,'Snore'));
FsFlow = Channels_Fs(strcmp(ChannelsList,'Flow'));
SndData = Snore;
winsizeT = 0.5;
winsizeT50 = 0.1;
percOverlap = 50;
FsDownsample=6000;
FsDS2 = F_samp;
FsBuffer = 1/(winsizeT*(percOverlap/100));
FsBuffer50 = 1/(winsizeT50*(percOverlap/100));
SegLen = 10*60; %seconds

SndIntensityArray = nan(floor(length(SndData)*(FsDS2/Fs)),1);
WinVarArray = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
ZCrArray = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
timeWin = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
timeDS2 = nan(floor(length(SndData)*(FsDS2/Fs)),1);
timeDS = nan(size(SndData));
frmnFreq1 = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
frmnFreq2 = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
frmnFreq3 = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
p70to250 = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
p250to500 = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
p500to750 = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
p750to1k = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
p1kto1p5k = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
p1p5kto2k = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
centroid = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
% pitch = nan(floor(length(SndData)*(FsBuffer/Fs)),1);
hr = nan(floor(length(SndData)*(FsBuffer/Fs)),1);

WinVarArray50 = nan(floor(length(SndData)*(FsBuffer50/Fs)),1);
timeWin50 = nan(floor(length(SndData)*(FsBuffer50/Fs)),1);

%% Process sounds
% isolate 10 min segments
% possibly temporary but might be faster to work in sections
% in which case we loop everything below
for ii = 1:ceil(length(SndData)/(SegLen*Fs))
    start = (ii-1)*SegLen*Fs + 1;
    finish = ii*SegLen*Fs;
    
    if finish > length(SndData) % end condition
        finish = length(SndData);
    end
    
    SndData10Temp = SndData(start:finish);
    t=(start:finish)'/Fs;
    
    startFlow = (ii-1)*SegLen*FsFlow + 1;
    finishFlow = ii*SegLen*FsFlow;
    Flow10 = Flow(startFlow:finishFlow);
    timeFlow = (startFlow:finishFlow)'/FsFlow;

    %% Filter and downsample data
    FilterBandwidth = [10 2000];
    SndData10 = SndData10Temp;
    FilteredSnd = ButterFilter(SndData10,FilterBandwidth(1),FilterBandwidth(2),Fs,5);
    FilteredSndDS = resample(FilteredSnd,FsDownsample,Fs);
    startDS = (ii-1)*SegLen*FsDownsample + 1; finishDS = ii*SegLen*FsDownsample;
    timeDS(startDS:finishDS) = (startDS:finishDS)'/FsDownsample;

    clear SndData10Temp
    
    %% Buffer data
    % dividing the signal into overlapping segments
    winsize = round(winsizeT*FsDownsample);
    winOverlap = floor((percOverlap/100)*winsize);

    %Put signal into windows
    num_windows = floor(length(FilteredSndDS)/winOverlap) - 1;
    total_length = winOverlap * (num_windows + 1);
    buffered = buffer(FilteredSndDS(1:total_length),winsize,winOverlap,'nodelay'); % This prevents the buffer function from prepending and appending zeros
        
    startWin = (ii-1)*SegLen*FsBuffer + 1; finishWin = (ii*SegLen*FsBuffer)-1;
    if finishWin > length(WinVarArray)
        finishWin = startWin + size(buffered,2)-1;
    end
    timeWin(startWin:finishWin) = (startWin:finishWin)'/FsBuffer;
    
    %% Buffer at 20Hz
    winsize50 = round(winsizeT50*FsDownsample);
    winOverlap50 = floor((percOverlap/100)*winsize50);

    %Put signal into windows
    num_windows = floor(length(FilteredSndDS)/winOverlap50) - 1;
    total_length = winOverlap50*(num_windows + 1);
    buffered50 = buffer(FilteredSndDS(1:total_length),winsize50,winOverlap50,'nodelay'); % This prevents the buffer function from prepending and appending zeros
        
    startWin50 = (ii-1)*SegLen*FsBuffer50 + 1; finishWin50 = (ii*SegLen*FsBuffer50)-1;
    if finishWin50 > length(WinVarArray50)
        finishWin50 = startWin50 + size(buffered50,2)-1;
    end
    timeWin50(startWin50:finishWin50) = (startWin50:finishWin50)'/FsBuffer50;
    
    %% Calculate features
    % Sound amplitude - might be undersampled @ 4Hz at the moment
    SndIntensity = (FilteredSndDS).^2; % converted to intensity
    SndIntensitydB = 10.*log10(SndIntensity/1e-12); %reference = 1e-12 W/m2
    SndIntensitydBDS = resample(SndIntensitydB,FsDS2,FsDownsample); % downsample to F_samp Hz
    startDS2 = (ii-1)*SegLen*FsDS2 + 1; finishDS2 = ii*SegLen*FsDS2;
    
    if finishDS2 > length(SndIntensityArray) % end condition
        finishDS2 = startDS2 + length(SndIntensitydBDS)-1;
    end
    timeDS2(startDS2:finishDS2) = (startDS2:finishDS2)'/FsDS2;
    
    SndIntensityArray(startDS2:finishDS2) = SndIntensitydBDS;
    
    % Windowed Variance, Zero-crossing rate
    WinVar50 = nanvar(buffered50,1);
    WinVar = nanvar(buffered,1);
    MSs = nanmean(buffered,1);
    ZSs = buffered - MSs;
    ZCr = sum(abs(sign(ZSs(1:end-1,:))-sign(ZSs(2:end,:))))/(2*winsize);
    
    WinVarArray50(startWin50:finishWin50) = WinVar50';
    WinVarArray(startWin:finishWin) = WinVar';
    ZCrArray(startWin:finishWin) = ZCr';
    
    % Formants
    [formantsfreq] = formants_freqCalc (buffered',FsDownsample);
    frmnFreq1(startWin:finishWin) = formantsfreq(:,1)';
    frmnFreq2(startWin:finishWin) = formantsfreq(:,2)';
    frmnFreq3(startWin:finishWin) = formantsfreq(:,3)';
    
    % Band power
    p70to250(startWin:finishWin) = bandpower(buffered,FsDownsample,[70 250]);
    p250to500(startWin:finishWin) = bandpower(buffered,FsDownsample,[250 500]);
    p500to750(startWin:finishWin) = bandpower(buffered,FsDownsample,[500 750]);
    p750to1k(startWin:finishWin) = bandpower(buffered,FsDownsample,[750 1000]);
    p1kto1p5k(startWin:finishWin) = bandpower(buffered,FsDownsample,[1000 1500]);
    p1p5kto2k(startWin:finishWin) = bandpower(buffered,FsDownsample,[1500 2000]);
    
    % Spectral centroid
    centroid(startWin:finishWin) = getCentroid(buffered,FsDownsample,winsize);
    
    % Pitch and harmonic ratio
    [f0,idx] = pitch(FilteredSndDS,FsDownsample);
    
%     hr = harmonicRatio(buffered,FsDownsample, ...
%                    'Window',hann(round(fs.*0.05),'periodic'), ...
%                    'OverlapLength',round(fs.*0.025));
    %%
    % Plot spectrogram and pitch
    freqRes = 5; %this means you can resolve frequencies that are 5 Hz apaprt
    timeRes = 1/freqRes; %Ali says you want this to be 0.1 for snoring
    L = floor(Fs/freqRes); %window size
    noverlap = floor(L*.85);  %85% overlap
    nfft = max(256,2^nextpow2(L));  %the DFT length, which can be longer than the window
    
%     f = 0:1:1000;
    spectrogram(FilteredSnd,L,noverlap,nfft,Fs,'yaxis');
    colormap(flipud(gray))
    
    % get data to rebuild spectrogram
    ch = get(gca,'ch');
    ax = gca;
    T = get(ch, 'XData') + start-1;
    F = get(ch, 'YData');
    P = get(ch,'CData');
%     t = linspace(0, length(FilteredSndDS), length(FilteredSndDS));
    clim = ax.CLim;

    %Plot the spectrogram
    figure('Position', [76 148 1180 470])
    ax1=subplot(3,1,1:2);
    imagesc(T, F, P, clim);

    set(gca,'YDir','normal')
%     xlabel("Time (s)");
    ylabel("Frequency (kHz)")
    set(gca,'Box','On','FontSize',12,'XTick',[])
    
    % overlay processed frequency data
    hold on
    pitchtime = linspace(T(1),T(end),length(f0));
    plot(pitchtime,f0/1000,'k.')
    
    ax2=subplot(3,1,3);
    plot(timeFlow,Flow10);
    xlabel("Time (s)");
    ylabel("Flow (L/s)")
    set(gca,'Box','On','FontSize',12)
    
    linkaxes([ax1,ax2],'x');
end
%% Interp windowed features
samplepoints = 1:length(timeWin);
querypoints = linspace(1, length(timeWin), length(timeDS2));

WinVarInterp = interp1(samplepoints', WinVarArray, querypoints', 'previous');
ZCrInterp = interp1(samplepoints', ZCrArray, querypoints', 'previous');
Frmn1Interp = interp1(samplepoints', frmnFreq1, querypoints', 'previous');
Frmn2Interp = interp1(samplepoints', frmnFreq2, querypoints', 'previous');
Frmn3Interp = interp1(samplepoints', frmnFreq3, querypoints', 'previous');
p70to250Interp = interp1(samplepoints', p70to250, querypoints', 'previous');
p250to500Interp = interp1(samplepoints', p250to500, querypoints', 'previous');
p500to750Interp = interp1(samplepoints', p500to750, querypoints', 'previous');
p750to1kInterp = interp1(samplepoints', p750to1k, querypoints', 'previous');
p1kto1p5kInterp = interp1(samplepoints', p1kto1p5k, querypoints', 'previous');
p1p5kto2kInterp = interp1(samplepoints', p1p5kto2k, querypoints', 'previous');
centroidInterp = interp1(samplepoints', centroid, querypoints', 'previous');

samplepoints50 = 1:length(timeWin50);
querypoints50 = linspace(1, length(timeWin50), length(timeDS2));
WinVarInterp50 = interp1(samplepoints50', WinVarArray50, querypoints50', 'previous');

%% Store features in structured array
SnoreStruct = struct();
SnoreStruct.FilteredSndDS = FilteredSndDS;
SnoreStruct.timeDS = timeDS;
SnoreStruct.SndIntensity = SndIntensityArray;
SnoreStruct.WinVar50 = WinVarInterp50;
SnoreStruct.WinVar = WinVarInterp;
SnoreStruct.ZCr = ZCrInterp;
SnoreStruct.frmnFreq1 = Frmn1Interp;
SnoreStruct.frmnFreq2 = Frmn2Interp;
SnoreStruct.frmnFreq3 = Frmn3Interp;
SnoreStruct.p70to250 = p70to250Interp;
SnoreStruct.p250to500 = p250to500Interp;
SnoreStruct.p500to750 = p500to750Interp;
SnoreStruct.p750to1k = p750to1kInterp;
SnoreStruct.p1kto1p5k = p1kto1p5kInterp;
SnoreStruct.p1p5kto2k = p1p5kto2kInterp;
SnoreStruct.centroid = centroidInterp;
% SnoreStruct.timeDS2 = timeDS2;
% SnoreStruct.timeWin = timeWin;
% SnoreStruct.timeWin = timeWin50;

end