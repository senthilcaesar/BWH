function [StimON, timeBuffered] = FindStimON(CHIN,leftt,rightt,dt,StartTime,stimtype,varargin)

global settings

%% Find HGNS stimulation start and stop from CHIN signal
Fs = 1/(dt);
StimSig = CHIN - mean(CHIN);
FilterBandwidth = [3 64];
FilteredStim = ButterFilter(StimSig,FilterBandwidth(1),FilterBandwidth(2),Fs,5);
Time = linspace(0,round(length(StimSig)/Fs),length(StimSig))';
% figure(22), plot(Time, StimSig, Time, FilteredStim)

%% Window size
freqRes = 20;
timeRes = 1/freqRes; %Ali says you want this to be 0.1 for snoring
winsize = floor(Fs/freqRes); %window size
percoverlap = 0.50;
noverlap = floor(winsize*percoverlap);  %75% overlap

% Buffer data
[buffered, ~] = bufferSignal(FilteredStim,length(FilteredStim),floor(Fs*timeRes),floor(Fs*percoverlap*timeRes),Fs);

% Time vector for buffered data
ncol = fix((length(FilteredStim)-noverlap)/(winsize-noverlap));
coloffsets = (0:(ncol-1))*(winsize-noverlap);
timeBuffered = (StartTime + (coloffsets+(winsize/2)')/Fs)';

% Get RMS
StimRMS = nanvar(buffered,1)';

%%
% Shift signal by offset
intervalBuffer = timeBuffered(2) - timeBuffered(1);
shift = round(noverlap*((1/Fs)/intervalBuffer));
StimRMS(1:shift) = []; nanpad = nan(shift,1);
StimRMS = [StimRMS;nanpad];
% StimRMS_Z = StimRMS - nanmean(StimRMS)./nanstd(StimRMS);

% figure(23), clf, hold on, plot(timeBuffered, StimRMS_Z)

%% Find threshold line
% threshopts = [50:-1:1];
% StimRMS_Z2 = StimRMS_Z;
% for nn = 1:length(threshopts)
%     thresh = prctile(StimRMS_Z,threshopts(nn)); % pick starting point
% %     StimRMS_Z2(StimRMS_Z < thresh) = nan;
%     threshSD(nn) = nanstd(StimRMS_Z2(StimRMS_Z < thresh));
%     figure(23), plot([timeBuffered(1) timeBuffered(end)],[thresh thresh])
% end

if strcmp(stimtype,'on_off')
    %% Use cutoff to find stim on periods (can do this in windowed approach if baseline float around too much
    % get region of the signal that has stim on/off 
    idxwithstimonoff = timeBuffered >= leftt & timeBuffered <= rightt;
    StimRMSsmall = StimRMS(idxwithstimonoff);
    StimRMS2 = StimRMSsmall;
    StimRMS2(StimRMS2 < prctile(StimRMSsmall,80)) = nan;
    timeBufferedSmall = timeBuffered(idxwithstimonoff);
    if rightt - leftt < 7 % condition where there is no good Chin so we have precisely label stim on in notepad
        StimONsmall = false(size(StimRMS2));
        StimONsmall(10:end-9) = true;    
    else % do the more complicated routine which finds stim ON within the given bounds
        % first bit here is to isolate periods of stim as NOT nan
        % will find nans and then do some gap filling
        % then find NOT nans

        % find nan indices
        nanfind = isnan(StimRMS2);

        % find location of consecutive nans
        i1 = find(diff([0 nanfind' 0]) == 1);
        out = i1(find(diff([0 nanfind' 0]) == -1) - i1 >= 1);

        % find length of consecutive nans
        measurements = regionprops(logical(nanfind), 'Area');
        theLengths = [measurements.Area];

        % fill in small gaps of nans 
        smallgaps = find(theLengths < 20);
        for kk = smallgaps
            indices = out(kk):out(kk)+theLengths(kk)-1;
            StimRMS2(indices) = StimRMSsmall(indices);
        end

        % now find consecutive NOT nans (these are periods of stim)
        nanfind = ~isnan(StimRMS2);

        % find location of consecutive nans
        i1 = find(diff([0 nanfind' 0]) == 1);
        out = i1(find(diff([0 nanfind' 0]) == -1) - i1 >= 1);

        % find length of consecutive nans
        measurements = regionprops(logical(nanfind), 'Area');
        theLengths = [measurements.Area];

        % min length = 2s = 2*Fs
        minstimduration = round(3.5/intervalBuffer);
        longstimlengths = theLengths(theLengths >= minstimduration);
        longstimstart = out(theLengths >= minstimduration);
        longstimend = longstimstart + longstimlengths - 1;

        StimONsmall = false(size(StimRMS2));
        for kk = 1:length(longstimstart)
            StimONsmall(longstimstart(kk):longstimend(kk)) = true;
        end

        % Idea: run through again with a step input that looks like what I've
        % generated with StimON

    end
    % Put into the larger StimON array
    StimON = false(length(timeBuffered),1);
    StimON(idxwithstimonoff) = StimONsmall;
    
    %%
    if 1
    % figure(22), clf, plot(Time, StimSig, Time, FilteredStim)
    figure(23), clf, hold on, plot(timeBufferedSmall, StimRMSsmall)
%     figure(24), plot(timeBuffered, StimRMS2)
    figure(23), plot(timeBufferedSmall, StimONsmall.*max(StimRMSsmall))
    end
end