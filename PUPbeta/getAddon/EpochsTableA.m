function [Evts]=EpochsTableA(Evts,DataEventHypnog_Mat,ChannelsList,BreathDataTable)
%% Make EpochsT
if ~isfield(Evts,'EpochT')
Epochs = Evts.Hypnogram;
EpochsT = table(Epochs);
EpochsT.EpochStart = Evts.Hypnogram_t;

zerocheck = std(diff(EpochsT.EpochStart))

EpochsT.EpochNum=(1:size(EpochsT,1))';
EpochsT=[EpochsT(:,3) EpochsT(:,1:2)];
Evts.EpochT=EpochsT;
end

%% Make time series of data where breaths were analyzed within BreathDataTable
% Can be improved if we also run analysis with Epochs=8 otherwise we lose lots of windows with epochs = unknown
Time=DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'Time')==1));
dt = Time(2)-Time(1);
SpO2ch=find(strcmp(ChannelsList,'SpO2')==1);
if ~isempty(SpO2ch)
    SpO2=DataEventHypnog_Mat(:,SpO2ch);
else
    SpO2=NaN*Time;
end
SpO2Available = ~isnan(SpO2)*1;

Evts.EpochT.EpochStarti = round((Evts.EpochT.EpochStart-Time(1))/dt + 1);
Evts.EpochT.EpochEndi = round((Evts.EpochT.EpochStart+30-Time(1))/dt);
Evts.EpochT.EpochStarti(Evts.EpochT.EpochStarti<1)=1;
Evts.EpochT.EpochEndi(Evts.EpochT.EpochStarti>length(Time))=length(Time);

Evts.EpochT.SpO2Available = 0*Evts.EpochT.EpochStart;
for i=1:length(Evts.EpochT.EpochStart)
    Evts.EpochT.SpO2Available(i) = mean(SpO2Available(Evts.EpochT.EpochStarti(i):Evts.EpochT.EpochEndi(i)));
end

if exist('BreathDataTable')
    BreathDataAvailable = AnalyzedWindowsXHzFromBDT(BreathDataTable,Time);
    Evts.EpochT.BreathDataAvailable = 0*Evts.EpochT.EpochStart;
    for i=1:length(Evts.EpochT.EpochStart)
        Evts.EpochT.BreathDataAvailable(i) = mean(BreathDataAvailable(Evts.EpochT.EpochStarti(i):Evts.EpochT.EpochEndi(i)));
    end
end

%% Calculate modified TST
TST = sum(Evts.EpochT.Epochs<=3)/2
TSTBreathDataAvailable = sum(Evts.EpochT.Epochs<=3.*Evts.EpochT.BreathDataAvailable)/2 %This is the appropriate denominator for auto AHI
TSTSpO2Available = sum(Evts.EpochT.Epochs<=3.*Evts.EpochT.SpO2Available)/2 %This is the appropriate denominator for auto AHI

% add spo2 analyzed1Hz and repeat-TSTanalyzed2
% check for what happens if spo2 is nan and ar=1 in getAHI
Evts.EpochTInfo.TST=TST;
Evts.EpochTInfo.TSTBreathDataAvailable=TSTBreathDataAvailable;
Evts.EpochTInfo.TSTSpO2Available=TSTSpO2Available;


