function [AHIData,Evts,rowNames] = getAHIEvtSubset(Evts,EpochsXHz,Time,Position,CPAPoff)
global settings
dt = 1/settings.Fs;

%Note "Time" not actually used here currently
%CPAPoff contains ones when there is no CPAP and we want AHI info

%% Setup channel data
% Time = DataEventHypnog_Mat(:,1);
% EpochsChan = find(strcmp(ChannelsList,'Epochs')==1);
% EpochsXHz=DataEventHypnog_Mat(:,EpochsChan);
%Alice: W=11,R=12,N1=13,N2=14,N3=15; Terrill: N3=0,N2=1,N1=2,R=3,W=4

%% check for overlapping events (not updated)
if 0
    List = [6 5 7 15 16];
    if size(DataEventHypnog_Mat,2)<16
        List(5)=[];
    end
    clear overlaptemp;
    for i=1:length(List)
        for j=1:length(List)
            if j>=i
                continue
            end
            overlaptemp(i,j) = dt*sum(DataEventHypnog_Mat(:,ColumnHeads(List(i)))&DataEventHypnog_Mat(:,ColumnHeads(List(j))));
        end
    end
    
    totaloverlapsec = sum(overlaptemp(:));
    if totaloverlapsec>1
        disp(['******* WARNING: respiratory event overlap by ' num2str(totaloverlapsec) ' sec *******']);
    end
end


%% Calculate AHI values

% EventCategoriesLabels= {'Arousal','Obstructive Apnea','Central Apnea','Mixed Apnea','Hypopnea','Unsure'};
% EventCategoriesCodes=[1 2 3 5 4 4];

stageoptions = {'EpochsXHz>=0&EpochsXHz<=3','EpochsXHz>=0&EpochsXHz<3','EpochsXHz==3','EpochsXHz==4','EpochsXHz==2','EpochsXHz==1','EpochsXHz==0'};

%cols: AllSleep,NREM,REM,W,N1,N2,N3
%rows: Duration(min)|ARindex|ApOindex|ApCindex|HypOindex|Mindex|HypCindex|AHI'

issupine = PositionSelector(Position,'Supine');

clear DurationandEventsPerHour_Supine_CPAPoff
rowNames={'Dur','ArI','ApOI','ApCI','HyOI','MAI','HyCI','AHI'}; %j
colNames={'AllSleepSup','NRemSup','RemSup','WSup','N1Sup','N2Sup','N3Sup'}; % k
EventNamesSup={};
for k=1:length(stageoptions)
    criteria = CPAPoff & eval(stageoptions{k}) & issupine==1;
    Duration = sum(criteria)*dt/60;
    EventNamestemp={[colNames{k} rowNames{1}]};
    EventNamesSup=[EventNamesSup;EventNamestemp];
    for j=1:6 %ARindex|ApOindex|ApCindex|HypOindex|Mindex|HypCindex
        EventsPerHour(j) = sum(criteria(Evts.RespT.starttimesi)& Evts.RespT.EventCodes(:,1)==j)/Duration*60;
        EventNamestemp={[colNames{k} rowNames{j+1}]};
        EventNamesSup=[EventNamesSup;EventNamestemp];
    end
    EventsPerHour(7) = sum(EventsPerHour(2:6));
    EventNamestemp={[colNames{k} rowNames{8}]};
    EventNamesSup=[EventNamesSup;EventNamestemp];
    DurationandEventsPerHour_Supine_CPAPoff(:,k) = [Duration,EventsPerHour];
end
EventNamesSup=EventNamesSup';
DurationandEventsPerHour_Supine_CPAPoff_ = reshape(DurationandEventsPerHour_Supine_CPAPoff,1,56);

clear DurationandEventsPerHour_AllPositions_CPAPoff EventNamestemp
rowNames={'Dur','ArI','ApOI','ApCI','HyOI','MAI','HyCI','ahi'}; %j
colNames={'AllSleepAllP','NRemAllP','RemAllP','WAllP','N1AllP','N2AllP','N3AllP'}; % k
EventNamesAllP={};
for k=1:length(stageoptions)
    criteria = CPAPoff&eval(stageoptions{k});
    Duration = sum(criteria)*dt/60;
    EventNamestemp={[colNames{k} rowNames{1}]};
    EventNamesAllP=[EventNamesAllP;EventNamestemp];
    for j=1:6
        EventsPerHour(j) = sum(criteria(Evts.RespT.starttimesi)& Evts.RespT.EventCodes(:,1)==j)/Duration*60;
        EventNamestemp={[colNames{k} rowNames{j+1}]};
        EventNamesAllP=[EventNamesAllP;EventNamestemp];
    end
    EventsPerHour(7) = sum(EventsPerHour(2:6));
    EventNamestemp={[colNames{k} rowNames{8}]};
    EventNamesAllP=[EventNamesAllP;EventNamestemp];
    DurationandEventsPerHour_AllPositions_CPAPoff(:,k) = [Duration,EventsPerHour];
end

AllPositionsAllStatesDurationName={'AllSleepAllPDur_','AllSleepAllPAHI'};
AllPositionsAllStatesDurationAHI = [DurationandEventsPerHour_AllPositions_CPAPoff(1,1) DurationandEventsPerHour_AllPositions_CPAPoff(8,1)];

EventNamesAllP=EventNamesAllP';
AllPositionsAllStatesDurationAHI2 = reshape(DurationandEventsPerHour_AllPositions_CPAPoff,1,56); %%% added 2016-11

clear Duration_Allpositions_CPAPofforon
rowNames={'dur'};
colNames={'AllSleepAllP','NRemAllP','RemAllP','WAllP','N1AllP','N2AllP','N3AllP'}; % k
DurNamesAllP={};
for k=1:length(stageoptions)
    criteria = eval(stageoptions{k});
    Duration_Allpositions_CPAPofforon(k) = sum(criteria)*dt/60;
    EventNamestemp={[colNames{k} rowNames{1}]};
    DurNamesAllP=[DurNamesAllP;EventNamestemp];
end
DurNamesAllP=DurNamesAllP';

%event duration data, medians [all positions, NREM all resp.events is 16/56, i.e. AHIData(144)]
clear AllEventDurations_AllPositions_CPAPoff
rowNames={'Dur','ArDur','ApODur','ApCDur','HyODur','MADur','HyCDur','SumEventDur'}; %j
colNames={'AllSleepAllP','NRemAllP','RemAllP','WAllP','N1AllP','N2AllP','N3AllP'}; % k
EventNamesMedDur={};
for k=1:length(stageoptions)
    criteria = CPAPoff&eval(stageoptions{k});
    for j=1:6
        AllEventDurations_AllPositions_CPAPoff(j) = median(Evts.RespT.EventDuration(criteria(Evts.RespT.starttimesi)& Evts.RespT.EventCodes(:,1)==j));
        EventNamestemp={[colNames{k} rowNames{j+1}]};
        EventNamesMedDur=[EventNamesMedDur;EventNamestemp];
    end
    AllEventDurations_AllPositions_CPAPoff(7) = median(Evts.RespT.EventDuration(criteria(Evts.RespT.starttimesi)& sum((Evts.RespT.EventCodes(:,1)==[2 4 5])')'));
    AllEventDurations_AllPositions_CPAPoff(8) = median(Evts.RespT.EventDuration(criteria(Evts.RespT.starttimesi)& sum((Evts.RespT.EventCodes(:,1)==[2:6])')'));
    EventNamestemp={[colNames{k} 'ObsDur']};
    EventNamesMedDur=[EventNamesMedDur;EventNamestemp];
    EventNamestemp={[colNames{k} 'ApHyDur']};
    EventNamesMedDur=[EventNamesMedDur;EventNamestemp];
    
    AllEventDurations_AllPositions_CPAPoff_(:,k) = [AllEventDurations_AllPositions_CPAPoff'];
end
EventNamesMedDur=EventNamesMedDur';
AllEventDurations_AllPositions_CPAPoff_ = reshape(AllEventDurations_AllPositions_CPAPoff_,1,56);

PercentTSTAllpositionsName={'AllSleepAllPTST','NRemAllPTST','RemAllPTST','WAllPTST','N1AllPTST','N2AllPTST','N3AllPTST'};
PercentTSTAllpositions_CPAPofforon = Duration_Allpositions_CPAPofforon/Duration_Allpositions_CPAPofforon(1)*100;

AHIData = [DurationandEventsPerHour_Supine_CPAPoff_  AllPositionsAllStatesDurationAHI  Duration_Allpositions_CPAPofforon PercentTSTAllpositions_CPAPofforon AllPositionsAllStatesDurationAHI2 AllEventDurations_AllPositions_CPAPoff_];
rowNames=[EventNamesSup AllPositionsAllStatesDurationName DurNamesAllP PercentTSTAllpositionsName EventNamesAllP EventNamesMedDur];

Tahi = array2table(AHIData);
Tahi.Properties.VariableNames = rowNames;

