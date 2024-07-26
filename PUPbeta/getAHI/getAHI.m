function [AHIData,Evts] = getAHI(SigT,CPAPoff,ChannelsList,EventsChannelName)
% Added ChannelsList as input to function---needed for calling in ConvertToMat
% --Please modify in any other functions where getAHI is called --9/29/2020 RMA

global settings 
if ~exist('EventsChannelName')
    EventsChannelName = 'EventsResp';
end
% if isfield(settings,'UseAutoScoredEventsForAHI') && settings.UseAutoScoredEventsForAHI==1
%     EventsChannelName = 'EventsRespAuto';
% end

%CPAPoff contains ones when there is no CPAP and we want AHI info 

%% Setup channel data
Time = SigT.Time;
dt = SigT.Time(2) - SigT.Time(1);

EpochsChan = find(strcmp(ChannelsList,'Epochs')==1);
if ~isempty(EpochsChan)
EpochsXHz=SigT.Epochs;
end

% added code for position codes depending on protocol
PosChan = find(strcmp(ChannelsList,'Position')==1);
if ~isempty(PosChan)
    PositionRaw=SigT.Position;
    Position = PositionTranslator(settings.positioncodes,settings.positioncodesout,PositionRaw);
else  % DLM added this case to handle NaN ColumnHeads(10), 
      % It artificially sets the entire night to supine
    if length(settings.supinepositioncode) > 1
        supinecode = settings.supinepositioncode(2);
    else
        supinecode = settings.supinepositioncode;
    end
    Position=ones(height(SigT),1)*supinecode;
    PositionRaw = Position; %to avoid code to break if no position channel available
end
%hypnogram=15-downsample(num,30); %Alice: W=11,R=12,N1=13,N2=14,N3=15; Terrill: N3=0,N2=1,N1=2,R=3,W=4

if 0 % will break with SigT--change columnheads here
    figure(11); set(gcf,'color',[1 1 1]);
    ax1(1)=subplot(3,1,1); plot(Time,EpochsXHz,Time,SigT(:,ColumnHeads(6))); box('off');
    legend('Hyp','Ar');
    ax1(2)=subplot(3,1,2); plot(Time,SigT(:,ColumnHeads(2))); box('off'); %Flow
    ax1(3)=subplot(3,1,3); plot(Time,SigT(:,ColumnHeads(5)));
    legend('OA','CA','H','MA'); box('off');
    linkaxes(ax1,'x');
end

%% check for overlapping events--will break with SigT--change columnheads here
if 0
List = [6 5 7 15 16];
if size(SigT,2)<16
    List(5)=[];
end

clear overlaptemp;
for i=1:length(List)
    for j=1:length(List)
        if j>=i
            continue
        end
        overlaptemp(i,j) = dt*sum(SigT(:,ColumnHeads(List(i)))&SigT(:,ColumnHeads(List(j))));
    end
end

totaloverlapsec = sum(overlaptemp(:));
if totaloverlapsec>1
    disp(['******* WARNING: respiratory event overlap by ' num2str(totaloverlapsec) ' sec *******']);
end
end

%% Extract events list from matrix
%EventCategories={'Central apnea','Obstructive apnea','Mixed apnea','Hypopnea','Desaturation','µ-arousal'};
%EventCols=[10 11 15 12 13 14];
%ColumnHeads=[1                   3                     2        9           10        11                    12         13       14         4           5       6       7          8     15 ];
%[ColumnHeads(1)=Time ColumnHeads(2)=Thorax 3=Pnasal 4=Hypnogram 5=Central Apneas 6=Obstructive(/Mixed) 7=Hypopnea 8=Desats 9=Arousals 10=RIP_Abdo 11=SpO2 12=C3A2 13=Position ]

%ARindex|ApOindex|ApCindex|HypOindex|Mindex|HypCindex

EventsRespChan = find(strcmp(ChannelsList,EventsChannelName)==1); %default: 'EventsResp'
% SigArray=SigT{:,:}; % table causing prob here
temp = SigT{:,EventsRespChan};
temp(temp>0)=1;
Iresp=1+find(diff(temp)>0);
Iresp2=1+find(diff(temp)<0);
[Iresp,Iresp2] = TidyStartEndEventList(Iresp,Iresp2,length(Time));

alwaysuseoriginalarousalscoring=0;
if alwaysuseoriginalarousalscoring
    EventsArChan = find(strcmp(ChannelsList,'EventsAr')==1);
else
    if isfield(settings,'useWSanalysisToReplaceAr') && settings.useWSanalysisToReplaceAr==1
        EventsArChan = find(strcmp(ChannelsList,'EventsArWS')==1);
    elseif isfield(settings,'useWSanalysisToReplaceAr') && settings.useWSanalysisToReplaceAr==2
        EventsArChan = find(strcmp(ChannelsList,'EventsArWSB')==1);
    else
        EventsArChan = find(strcmp(ChannelsList,'EventsAr')==1);
    end
end

temp = SigT{:,EventsArChan};
temp(temp>0)=1;
Iar=1+find(diff(temp)>0);
Iar2=1+find(diff(temp)<0);
[Iar,Iar2] = TidyStartEndEventList(Iar,Iar2,length(Time));

Evts.starttimesi=[Iresp;Iar];
Evts.starttimes = Time(Evts.starttimesi);
Evts.codes=[SigT{Iresp,EventsRespChan}; SigT{Iar,EventsArChan}];
Evts.duration = [(Iresp2-Iresp)*dt;(Iar2-Iar)*dt];
Evts.state=EpochsXHz(Evts.starttimesi);
Evts.position=Position(Evts.starttimesi);
Evts.positionraw=PositionRaw(Evts.starttimesi); 

%% Calculate AHI values
%Evts.starttimesi=round((Evts.starttimes/dt)+1)

stageoptions = {'EpochsXHz>=0&EpochsXHz<=3','EpochsXHz>=0&EpochsXHz<3','EpochsXHz==3','EpochsXHz==4','EpochsXHz==2','EpochsXHz==1','EpochsXHz==0'};

%cols: AllSleep,NREM,REM,W,N1,N2,N3
%rows: Duration(min)|ARindex|ApOindex|ApCindex|HypOindex|Mindex|HypCindex|AHI'
%issupine = sum((Position==settings.supinepositioncode),2); %2 denotes summing each row
%issupinenrem = eval(stageoptions{2}) & issupine==1;

issupine = PositionSelector(Position,'Supine');

clear DurationandEventsPerHour_Supine_CPAPoff
for k=1:length(stageoptions)
    criteria = CPAPoff & eval(stageoptions{k}) & issupine==1;
    Duration = sum(criteria)*dt/60;
    for j=1:6
        EventsPerHour(j) = sum(criteria(Evts.starttimesi)&Evts.codes(:,1)==j)/Duration*60;
    end
    EventsPerHour(7) = sum(EventsPerHour(2:6));
    DurationandEventsPerHour_Supine_CPAPoff(:,k) = [Duration,EventsPerHour];
end
DurationandEventsPerHour_Supine_CPAPoff_ = reshape(DurationandEventsPerHour_Supine_CPAPoff,1,56);

clear DurationandEventsPerHour_AllPositions_CPAPoff
for k=1:length(stageoptions)
    criteria = CPAPoff&eval(stageoptions{k});
    Duration = sum(criteria)*dt/60;
    for j=1:6
        EventsPerHour(j) = sum(criteria(Evts.starttimesi)&Evts.codes(:,1)==j)/Duration*60;
    end
    EventsPerHour(7) = sum(EventsPerHour(2:6));
    DurationandEventsPerHour_AllPositions_CPAPoff(:,k) = [Duration,EventsPerHour];
end
AllPositionsAllStatesDurationAHI = [DurationandEventsPerHour_AllPositions_CPAPoff(1,1) DurationandEventsPerHour_AllPositions_CPAPoff(8,1)];
AllPositionsAllStatesDurationAHI2 = reshape(DurationandEventsPerHour_AllPositions_CPAPoff,1,56); %%% added 2016-11

clear Duration_Allpositions_CPAPofforon
for k=1:length(stageoptions)
    criteria = eval(stageoptions{k});
    Duration_Allpositions_CPAPofforon(k) = sum(criteria)*dt/60;
end

%event duration data, medians [all positions, NREM all resp.events is 16/56, i.e. AHIData(144)]
clear AllEventDurations_AllPositions_CPAPoff
for k=1:length(stageoptions)
criteria = CPAPoff&eval(stageoptions{k});
for j=1:6
    AllEventDurations_AllPositions_CPAPoff(j) = median(Evts.duration(criteria(Evts.starttimesi)&Evts.codes(:,1)==j));
end
AllEventDurations_AllPositions_CPAPoff(7) = median(Evts.duration(criteria(Evts.starttimesi)&sum((Evts.codes(:,1)==[2 4 5])')'));
AllEventDurations_AllPositions_CPAPoff(8) = median(Evts.duration(criteria(Evts.starttimesi)&sum((Evts.codes(:,1)==[2:6])')'));
AllEventDurations_AllPositions_CPAPoff_(:,k) = [AllEventDurations_AllPositions_CPAPoff'];
end
AllEventDurations_AllPositions_CPAPoff_ = reshape(AllEventDurations_AllPositions_CPAPoff_,1,56);

PercentTSTAllpositions_CPAPofforon = Duration_Allpositions_CPAPofforon/Duration_Allpositions_CPAPofforon(1)*100;  

AHIData = [DurationandEventsPerHour_Supine_CPAPoff_  AllPositionsAllStatesDurationAHI  Duration_Allpositions_CPAPofforon PercentTSTAllpositions_CPAPofforon AllPositionsAllStatesDurationAHI2 AllEventDurations_AllPositions_CPAPoff_];



