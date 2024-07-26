function [AHIData,Evts] = getAHI_(DataEventHypnog_Mat,CPAPoff)
global settings ColumnHeads

%% Setup channel data
Time = DataEventHypnog_Mat(:,1);
dt = 1/settings.Fs;
EpochsXHz=DataEventHypnog_Mat(:,ColumnHeads(4));
% added code for position codes depending on protocol
if ~isnan(ColumnHeads(10))
    PositionRaw=DataEventHypnog_Mat(:,ColumnHeads(10));
    Position = PositionTranslator(settings.positioncodes,settings.positioncodesout,PositionRaw);
    
else  % DLM added this case to handle NaN ColumnHeads(10), 
      % It artificially sets the entire night to supine
    if length(settings.supinepositioncode) > 1
        supinecode = settings.supinepositioncode(2);
    else
        supinecode = settings.supinepositioncode;
    end
    Position=ones(size(DataEventHypnog_Mat,1),1)*supinecode;
end
%hypnogram=15-downsample(num,30); %Alice: W=11,R=12,N1=13,N2=14,N3=15; Terrill: N3=0,N2=1,N1=2,R=3,W=4

if 0
    figure(11); set(gcf,'color',[1 1 1]);
    ax1(1)=subplot(3,1,1); plot(Time,EpochsXHz,Time,DataEventHypnog_Mat(:,ColumnHeads(6))); box('off');
    legend('Hyp','Ar');
    ax1(2)=subplot(3,1,2); plot(Time,DataEventHypnog_Mat(:,ColumnHeads(2))); box('off'); %Flow
    ax1(3)=subplot(3,1,3); plot(Time,DataEventHypnog_Mat(:,ColumnHeads(5)));
    legend('OA','CA','H','MA'); box('off');
    linkaxes(ax1,'x');
end

%% check for overlapping events
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

%% Extract events list from matrix
%EventCategories={'Central apnea','Obstructive apnea','Mixed apnea','Hypopnea','Desaturation','µ-arousal'};
%EventCols=[10 11 15 12 13 14];
%ColumnHeads=[1                   3                     2        9           10        11                    12         13       14         4           5       6       7          8     15 ];
%[ColumnHeads(1)=Time ColumnHeads(2)=Thorax 3=Pnasal 4=Hypnogram 5=Central Apneas 6=Obstructive(/Mixed) 7=Hypopnea 8=Desats 9=Arousals 10=RIP_Abdo 11=SpO2 12=C3A2 13=Position ]

%ARindex|ApOindex|ApCindex|HypOindex|Mindex|HypCindex

Iresp=1+find(diff(DataEventHypnog_Mat(:,ColumnHeads(5)))>0);
Iresp2=1+find(diff(DataEventHypnog_Mat(:,ColumnHeads(5)))<0);
[Iresp,Iresp2] = TidyStartEndEventList(Iresp,Iresp2,length(Time));

Iar=1+find(diff(DataEventHypnog_Mat(:,ColumnHeads(6)))>0);
Iar2=1+find(diff(DataEventHypnog_Mat(:,ColumnHeads(6)))<0);
[Iar,Iar2] = TidyStartEndEventList(Iar,Iar2,length(Time));

Evts.starttimesi=[Iresp;Iar];
Evts.startimes=Time(Evts.starttimesi);
Evts.codes=[DataEventHypnog_Mat(Iresp,ColumnHeads(5)); DataEventHypnog_Mat(Iar,ColumnHeads(6))];
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



