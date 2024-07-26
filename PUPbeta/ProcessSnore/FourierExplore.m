% Exploring the sound data
n=1; m=2; ii=1; SegLen = 10*60; %seconds
%%
for freqRes = 5:5:15
    %%
    % Plot spectrogram and pitch
    % freqRes = 5; %this means you can resolve frequencies that are 5 Hz apaprt
    timeRes = 1/freqRes; %Ali says you want this to be 0.1 for snoring
    L = floor(Fs/freqRes); %window size
    noverlap = floor(L*.85);  %85% overlap
    nfft = max(256,2^nextpow2(L));  %the DFT length, which can be longer than the window

    disp(['window size = ', num2str(L)])
    %%
    % Pitch and harmonic ratio
    [f0,idx] = pitch(FilteredSnd,FsDownsample,...
        'WindowLength',L, ...
        'OverlapLength',noverlap);

    figure('Visible','off');    
    spectrogram(FilteredSnd,L,noverlap,nfft,Fs,'yaxis');

    % get data to rebuild spectrogram
    ch = get(gca,'ch');
    ax = gca;
    T = get(ch, 'XData') + start-1;
    F = get(ch, 'YData');
    P = get(ch,'CData');
    %     t = linspace(0, length(FilteredSndDS), length(FilteredSndDS));
    clim = ax.CLim;

%     %Plot the spectrogram
%     figure('Position', [76 148 1180 470])
%     ax1=subplot(3,1,1:2);
%     imagesc(T, F, P, clim);
%     colormap(flipud(gray))
% 
%     set(gca,'YDir','normal', 'ylim',[0 1])
%     %     xlabel("Time (s)");
%     ylabel("Frequency (kHz)")
%     set(gca,'Box','On','FontSize',12,'XTick',[])

    FsBuffer = 1/(L/Fs*(85/100));
    %Put signal into windows
    num_windows = floor(length(FilteredSnd)/noverlap) - 1;
    total_length = noverlap * (num_windows + 1);
    buffered = buffer(FilteredSnd(1:total_length),winsize,noverlap,'nodelay'); % This prevents the buffer function from prepending and appending zeros

    startWin = (ii-1)*SegLen*FsBuffer + 1; finishWin = (ii*SegLen*FsBuffer)-1;
    if finishWin > length(WinVarArray)
        finishWin = startWin + size(buffered,2)-1;
    end
    timeWin(startWin:finishWin) = (startWin:finishWin)'/FsBuffer;

     %% Window with inharmonic sound
    % plot PSDs for interesting parts of the data
    figure(23)
    subplot(3,2,n);
    tStart = 558.2;
    sampStart = round(tStart*(1/0.03));
    plot(F,P(:,sampStart))
    hold on
    set(gca, 'xlim',[0 1])
    ylims = get(gca,'ylim');
    plot([f0(sampStart)/1000 f0(sampStart)/1000], [ylims(1) ylims(2)])
    hold off
    ylabel('Power (dB/Hz')

    %% Use median of several consecutive windows (i.e. median filtered)
    numcwins = 5;
    sampStart = round(tStart*(1/0.03));

    subplot(3,2,m);
    plot(F,median(P(:,sampStart:sampStart+numcwins),2))
    hold on
    set(gca, 'xlim',[0 1])
    ylims = get(gca,'ylim');
    plot([f0(sampStart)/1000 f0(sampStart)/1000], [ylims(1) ylims(2)])
    hold off
    if freqRes == 15
        xlabel('Frequency (Hz)')
    end
    
    %% Window with harmonic sound
    tStart = 561.9;
    sampStart = round(tStart*(1/0.03));

    figure(24)
    subplot(3,2,n)
    plot(F,P(:,sampStart))
    set(gca, 'xlim',[0 1])
    ylims = get(gca,'ylim');
    hold on
    plot([f0(sampStart)/1000 f0(sampStart)/1000], [ylims(1) ylims(2)])
    ylabel('Power (dB/Hz')
    hold off
    %% Use median of several consecutive windows (i.e. median filtered)
    numcwins = 5;
    sampStart = round(tStart*(1/0.03));

    subplot(3,2,m)
    plot(F,median(P(:,sampStart:sampStart+numcwins),2))
    set(gca, 'xlim',[0 1])
    ylims = get(gca,'ylim');
    hold on
    plot([f0(sampStart)/1000 f0(sampStart)/1000], [ylims(1) ylims(2)])
    hold off
    if freqRes == 15
        xlabel('Frequency (Hz)')
    end
    
    %% count for subplot
    n = n+2; m = m+2;
end