function [StartRowPreced5min,EndRowPreced5min,StartRowControl5min,EndRowControl5min,BrDataTable]...
    = findStableBreathing(Sigs,BrDataTable,Fs)

%%
% Remove duplicates
BrDataTableOriginal = BrDataTable;
BrDataTable(BrDataTable.FDuplicated~=0,:) = [];

events = sum([Sigs.Apneas,Sigs.CApneas,Sigs.Hypopneas],2) > 0;
    
% Find 60 second periods of stable breathing
% i.e. no arousal, no resp events, during sleep
winsizeT = 60;
winsizeT2 = 300;

% Initialize array
RowsWStableBreathing = zeros(size(BrDataTable,1),1);
IdxWStableBreathing = zeros(size(events,1),1);

for rownum = 1:size(BrDataTable,1)
    startRow = rownum;
    endRow = find(BrDataTable.Time_start <= ...
        BrDataTable.Time_start(startRow)+winsizeT, 1, 'last');
    
    % Find start and stop time for timeseries data
    startT = BrDataTable.Time_start(startRow);
    endT = BrDataTable.Time_end(endRow);
    startIdx = startT*Fs;
    endIdx = endT*Fs;
    
    % Check for arousal, events, and sleep within start/end times    
    awakeIdx = Sigs.SleepStage(startIdx:endIdx) >= 4;
    awake = sum(awakeIdx)/(endIdx-startIdx+1) > 0.1;
%     sum(awakeIdx)/(endIdx-startIdx+1) 
    
    if sum(Sigs.Arousal(startIdx:endIdx)) == 0 && ...
            sum(events(startIdx:endIdx)) == 0 && ...
            ~awake
        
        RowsWStableBreathing(rownum:rownum+winsizeT-1) = 1;
        IdxWStableBreathing(startIdx:endIdx) = 1;
        
    end
    
end

% Find 5 mins that precede stable breathing which have at least 3 events
% diffIdxWStable = diff(IdxWStableBreathing);
% StableBrStartIdx = find(diffIdxWStable == 1); % these are one sample behind
% StableBrEndIdx = find(diffIdxWStable == -1); % these are one sample behind
diffRowWStable = diff(RowsWStableBreathing);
StableBrStartRow = find(diffRowWStable == 1); % these are one sample behind
StableBrEndRow = find(diffRowWStable == -1); % these are one sample behind

if length(StableBrStartRow) > length(StableBrEndRow)
    StableBrEndRow(end+1) = size(BrDataTable,1);
end

% Find start/stop of event
EventStart = [diff(BrDataTable.Etype) > 0;0];

% StartPreced5min = nan(size(StableBrStartIdx));
% EndPreced5min = nan(size(StableBrStartIdx));
StartRowPreced5min = nan(size(StableBrStartRow));
EndRowPreced5min = nan(size(StableBrStartRow));
Unstable5min = zeros(length(Sigs.Flow),1);

for ii = 1:length(StableBrStartRow)
    % Find start and end row for preceding 5 mins
    timeDiff = BrDataTable.Time_start(StableBrStartRow(ii)) - winsizeT2;
    [~, startrow] = min(abs(BrDataTable.Time_start - timeDiff));
    endrow = StableBrStartRow(ii)-1;
    
    % check if more than N events
    if sum(EventStart(startrow:endrow)) >= 3
        StartRowPreced5min(ii) = startrow;
        EndRowPreced5min(ii) = endrow;
        
        % convert to indices and plot
        startIdx = BrDataTable.Time_start(startrow)*Fs;
        endIdx = BrDataTable.Time_start(endrow)*Fs;
        Unstable5min(startIdx:endIdx) = 1;  
    end
end

StartRowPreced5min(isnan(StartRowPreced5min)) = [];
EndRowPreced5min(isnan(EndRowPreced5min)) = [];

if 1
    figure(32), plot(Sigs.Time, Sigs.Flow, Sigs.Time, IdxWStableBreathing*-1,...
        Sigs.Time, events, Sigs.Time, Unstable5min)
end

% Generate table with control data (and more for future)
% 5 min window with 50% overlap
EventStartLong = diff(events)>0;
EventsWide = buffer(EventStartLong, 300*Fs, 150*Fs);
NumEvtsPerWin = sum(EventsWide,1);
NumEvtsPerWinShift = circshift(NumEvtsPerWin,-1);

UnstableWide = buffer(Unstable5min,300*Fs,150*Fs);
PercUnstablePerWin = sum(UnstableWide,1)./(300*Fs);

% find windows that DO NOT include unstable preceded by stable and are
% followed by windows with >= 3 events in them

ControlWins = find(NumEvtsPerWin > 3 &...
                   NumEvtsPerWinShift > 3 &...
                   PercUnstablePerWin == 0);

StartRowControl5min = nan(length(ControlWins),1);
EndRowControl5min = nan(length(ControlWins),1);

for jj = 1:length(ControlWins)
    StartT = ControlWins(jj)*150;
    EndT = StartT + 300;
    
    [~,StartRow] = min(abs(BrDataTable.Time_start - StartT));
    [~,EndRow] = min(abs(BrDataTable.Time_start - EndT));
    if BrDataTable.Etype(StartRow) > 0
        StartRow = find(BrDataTable.Etype(1:StartRow) == 0,1,'last') - 10;
    elseif BrDataTable.Etype(StartRow) > 0
        EndRow = EndRow + find(BrDataTable.Etype(EndRow:end) == 0,1,'first') + 10;
    end
    
    StartRowControl5min(jj) = StartRow;
    EndRowControl5min(jj) = EndRow;
end



