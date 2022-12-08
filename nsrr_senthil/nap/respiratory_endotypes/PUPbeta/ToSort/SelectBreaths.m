%% Manual flow shape selection 
% Input subject info here
clear
subID = '1764';
treatment = 1; %Options are 'Baseline' or 'OA'
load('FLBreaths.mat')

%% Load converted file
if treatment
    directory = 'J:\PEOPLE\FACULTY\SANDS\OralApplianceMM2018\Converted\Tx\';
    subjectNames = readtable('ConvertDataSpreadsheet_Tx.xlsx','Sheet', 1,'Range','E3:E35');
    subidx = find(contains(subjectNames.Signals, subID));
    filename = subjectNames.Signals{subidx};
    load([directory, filename(1:end-4),'_XHz.mat'])
    night = 'OA';
else
    directory = 'J:\PEOPLE\FACULTY\SANDS\OralApplianceMM2018\Converted\';
    subjectNames = readtable('ConvertDataSpreadsheet.xlsx','Sheet', 1,'Range','E3:E35');
    subidx = find(contains(subjectNames.Signals, subID));
    filename = subjectNames.Signals{subidx};
    load([directory, filename(1:end-4),'_XHz.mat'])
    night = 'Baseline';
end

% Assign variable names
Time = DataEventHypnog_Mat(:,strcmp(ChannelsList, 'Time'));
EpochsXHz = DataEventHypnog_Mat(:,strcmp(ChannelsList, 'Epochs'));
EventsArXHz = DataEventHypnog_Mat(:,strcmp(ChannelsList, 'EventsAr'));
WakeSleep = DataEventHypnog_Mat(:,strcmp(ChannelsList, 'WakeSleep'));
Pmask = DataEventHypnog_Mat(:,strcmp(ChannelsList, 'Pmask'));
Pos = DataEventHypnog_Mat(:,strcmp(ChannelsList, 'Position'));
Flow = DataEventHypnog_Mat(:,strcmp(ChannelsList, 'Flow'));

%% Find breath shapes - should be able to add and remove breaths
figure(1)
Xminute=10;
Yminute=5;
dsf=5;
X=4;
ax2=[];

if exist('Pepi')
    ax2(4)=subplot(X,1,4); plot(downsample(Time,dsf),downsample(Pepi,dsf));
    X=X+1;
elseif exist('Pes')
    ax2(4)=subplot(X,1,4); plot(downsample(Time,dsf),downsample(Pes,dsf));
    X=X+1;
end

ax2(1)=subplot(X,1,1); plot(downsample(Time,dsf),downsample(EpochsXHz,dsf),'b',downsample(Time,dsf),downsample(EventsArXHz,dsf),'r');
if exist('WakeSleep'), hold('on'); plot(downsample(Time,dsf),downsample(WakeSleep,dsf),'color',[1 0.3 0.3]); hold('off'); end
ax2(2)=subplot(X,1,2); plot(downsample(Time,dsf),downsample(Pmask,dsf));
ax2(3)=subplot(X,1,3); plot(downsample(Time,dsf),downsample(Pos,dsf)); ylim([-1 1]);
ax2(X)=subplot(X,1,X); plot(downsample(Time,dsf),downsample(Flow,dsf)); ylim([-1 1]);

if sum(strcmp(fieldnames(FLBreathsAdd), ['DPW',subID]))
    hold on
    %Plot added breaths
    vals = FLBreathsAdd.(['DPW',subID]).(night);
    for ii = 1:size(vals,1)
        plot(Time(vals(ii,1):vals(ii,2)),...
            Flow(vals(ii,1):vals(ii,2)), 'g')
    end
    
    %Cover-up with removed indices
    vals = FLBreathsRm.(['DPW',subID]).(night);
    for ii = 1:size(vals,1)
        plot(Time(vals(ii,1):vals(ii,2)),...
            Flow(vals(ii,1):vals(ii,2)), 'r')
    end
end

linkaxes(ax2,'x');
for i=1:length(ax2)-1
    set(ax2(i),'Xtick',[],'Xcolor',[1 1 1]);
end
for i=1:length(ax2)
    set(ax2(i),'tickdir','out','box','off');
end

global xvalues yvalues range xvalues_rm yvalues_rm
xvalues=[]; yvalues=[]; xvalues_rm=[]; yvalues_rm=[]; 
range=Xminute*60;
plotwithsliderandselectLR([Time(1) Time(end)]);

%% Adjust timings and store breaths
% convert xvalues to indices
xIdx = nan(size(xvalues));
for ii = 1:size(xvalues,1)
    [~, start] = min(abs(Time - xvalues(ii,1)));
    [~, stop] = min(abs(Time - xvalues(ii,2)));
    xIdx(ii,1:2) = [start stop];
end

xIdx_rm = nan(size(xvalues_rm));
for ii = 1:size(xvalues_rm,1)
    [~, start] = min(abs(Time - xvalues_rm(ii,1)));
    [~, stop] = min(abs(Time - xvalues_rm(ii,2)));
    xIdx_rm(ii,1:2) = [start stop];
end

%add data to structure
if sum(strcmp(fieldnames(FLBreathsAdd), ['DPW',subID]))
   FLBreathsAdd.(['DPW',subID]).(night) = [FLBreathsAdd.(['DPW',subID]).(night);...
                                           xIdx xvalues];
else
    FLBreathsAdd.(['DPW',subID]).(night) = [xIdx xvalues];
end

if sum(strcmp(fieldnames(FLBreathsRm), ['DPW',subID]))
    FLBreathsRm.(['DPW',subID]).(night) = [FLBreathsRm.(['DPW',subID]).(night);...
                                           xIdx_rm xvalues_rm];
else
    FLBreathsRm.(['DPW',subID]).(night) = [xIdx_rm xvalues_rm];
end

%% Identify start and stop of each breath selected and store in Breath table
% load('FLBreaths.mat')
load('BreathTableBig.mat')

BreathDataTableSub = BreathDataTableBig(BreathDataTableBig.Subject == ...
    str2num(subID), :);
BreathDataTableSub.FLselect = zeros(size(BreathDataTableSub,1),1);

% Add breaths
for ii = 1:size(FLBreathsAdd.(['DPW',subID]).(night),1)
    diff1 = BreathDataTableSub.Time_mid - FLBreathsAdd.(['DPW',subID]).(night)(ii,3);
    startRow = sum(diff1 < 0) + 1; %finds the first positive value (when Time_mid first exceeds the start selection)
    FLBreathsAdd.(['DPW',subID]).(night)(ii,5) = startRow;
    
    diff2 = BreathDataTableSub.Time_mid - FLBreathsAdd.(['DPW',subID]).(night)(ii,4);
    stopRow = sum(diff2 < 0); %finds last negative value (when Time_mid closest to, but before stop selection)
    FLBreathsAdd.(['DPW',subID]).(night)(ii,6) = stopRow;
    
    BreathDataTableSub.FLselect(startRow:stopRow) = 1;
end

% Remove breaths
for ii = 1:size(FLBreathsRm.(['DPW',subID]).(night),1)
    diff1 = BreathDataTableSub.Time_mid - FLBreathsRm.(['DPW',subID]).(night)(ii,3);
    startRow = sum(diff1 < 0) + 1; %finds the first positive value (when Time_mid first exceeds the start selection)
    FLBreathsRm.(['DPW',subID]).(night)(ii,5) = startRow;
    
    diff2 = BreathDataTableSub.Time_mid - FLBreathsRm.(['DPW',subID]).(night)(ii,4);
    stopRow = sum(diff2 < 0); %finds last negative value (when Time_mid closest to, but before stop selection)
    FLBreathsRm.(['DPW',subID]).(night)(ii,6) = stopRow;
    
    BreathDataTableSub.FLselect(startRow:stopRow) = 0;
end


BreathDataTableBig.FLselect(BreathDataTableBig.Subject == ...
    str2num(subID), :) = BreathDataTableSub.FLselect;

%% Save FLBreaths and BreathTable
save('FLBreaths.mat', 'FLBreathsAdd', 'FLBreathsRm')
save('BreathTableBig.mat', 'BreathDataTableBig')