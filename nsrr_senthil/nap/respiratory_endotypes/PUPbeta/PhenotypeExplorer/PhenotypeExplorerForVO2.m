
%% Load options
clear all
%close all
%F:\Work\Projects SS\Overweight and obese nonapneics\Draft and Data\1676N1BariatricTrait
%filenameanddir='J:\PEOPLE\POST DOC\SANDS\Phenotyping and Edi\1722';
%directory = 'E:\Work\Projects SS\Phenotyping and Edi\';
directory = 'J:\PEOPLE\POST DOC\DEMELO\Analysed_Camila\';
cd(directory);
files = {'815SSPhenotypeNava','1722SSPhenotype', '1341Phenotype', '1731Phenotype','1723Phenotype','941SSPhenotype','1657SSPhenotype', '1657TiagabineN1' ,'1334SSPhenotype','1717MuscleStimN2','475MuscleStimN2','1313Phenotype','1729MuscleStimN2','383MetabolismPhenotype','1309MetabolismPhenotype','1469MetabolismPhenotype','1565MetabolismPhenotype', '1263MetabolismPhenotype',  '929SSPhenotype', '1264PhenotypeNAVA', '1343SSPhenotype', '1469Phenotype', '1722TiagabineN1', '1731TiagabineN2', '1738Phenotype', '1742TiagabineN2', '1743Phenotype', '1757Phenotype'}; %'383TiagabineN2'
study=2;
filename = files{study};
%filename = '383TiagabineN2'; %941PhenotypeSS, 815SSPhenotypeNAVA1657/941 %checkedzeroflowactives: 1429, 941, 1264(fixed), 1722(fixed),1657,1343,1723,1469,1429,1309,815,533, 1710,1708
filenameanddir=[directory filename];

%1469 - could not get Ccw
filehandle = matfile(filenameanddir);
w = whos('-file',filenameanddir);

%% Get channels
clear channelnameoptions
if 1 %load these
    channelnameoptions.Evts={'Evts','New_Evts'};
    channelnameoptions.Evts2={'Evts2','New_Evts2'};
    channelnameoptions.Epochs={'Epochs'};
    channelnameoptions.Flow={'Flow','Vflow','VFlow'}; %'Vflow'
    channelnameoptions.Pmask={'Pmask','PMask'};
    channelnameoptions.Position={'Position','Pos','pos','position'};
    channelnameoptions.SaO2={'SaO2','SpO2','Sat','Sao2','Spo2','O2sat','o2sat'};
    channelnameoptions.PCO2={'PCO2','CO2_Ana','CO2_anal','CO2'};
    channelnameoptions.PO2={'PO2','pO2','O2_Ana','O2_anal'};
    %channelnameoptions.RC={'Thorax','RC','Chest','CHEST'};
    %channelnameoptions.ABD={'Abdomen','ABD','Abdom','ABDM','ABDO'};
    channelnameoptions.ThRIP={'ThNoxRIP'};
    channelnameoptions.AbRIP={'AbNoxRIP'};    
    channelnameoptions.kNoxRIP={'kNoxRIP'}; 
    channelnameoptions.kFlow={'kFlow'}; 
    channelnameoptions.alphaFlow={'alphaFlow'}; 
    %channelnameoptions.Pepi={'Pepi','Epi','PEpi','EPI'};
    channelnameoptions.Pes={'Pes'};
    %channelnameoptions.Edi={'Edi'};
    %channelnameoptions.GGpmax={'GGpmax','GGPmax','EMGggpmax'};
    %channelnameoptions.EEG_C3_A2={'EEG_C3_A2','C3_A2'};
    %channelnameoptions.EEG_C4_A1={'EEG_C4_A1','C4_A1'};
    %channelnameoptions.EEG_O2_A1={'EEG_O2_A1','O2_A1'};
    %channelnameoptions.EEG_F3_A2={'EEG_F3_A2','F3_A2'};
end


channelnamestemp=fieldnames(channelnameoptions);

for i=1:length(channelnamestemp)
    temp=eval(['channelnameoptions.' char(channelnamestemp(i))])
    for n=1:length(temp)
        %Does it exist?
        for j=1:length(w)
            foundamatch=strcmp(w(j).name,char(temp(n)));
            if foundamatch
                eval([char(temp(1)) '=filehandle.' char(temp(n))]);
                break
            end
        end
    end
end

%% bugcorrection
try
    alphaFlow.start = ThNoxRIP.start;
    alphaFlow.interval = ThNoxRIP.interval;    
catch me
end


%% Options 2

plotfigs=1;
saveontherun=1;
saveattheend=1;
usePnasalasflow=0;
analyseSaO2=0;
analyseO2=0;

%% Resample flow if at high sampling rate
if Flow.interval<0.008
    dt=0.008;
    %Flow.values=resample(Flow.values,round(dt/Flow.interval),1);
    Flow.values=downsample(Flow.values,round(dt/Flow.interval));
    Flow.interval=0.008;
    Flow.length=length(Flow.values);
end

%% Setup Time array
dt=Flow.interval;
Fs=1/dt;

a=Flow.start;
Time=(a:dt:(a+dt*(Flow.length-1)))';

%% Epochs / Staging channel
EpochsXHz=interp1((Epochs.times+15),double(Epochs.codes(:,1)'),Time,'nearest','extrap');

%% Events info
%Evts list is just called Evts
if exist('New_Evts','var'), Evts=New_Evts; end
if exist('New_Evts2','var'), Evts2=New_Evts2; end

backupEvts = Evts;

%add extra row of zeros with time = inf just in case there is an event without an end in the first of two concatenated sets of events.
Evts.codes(size(Evts.codes,1)+1,:)=[0 0 0 0];
Evts.times(size(Evts.times,1)+1,:)=Inf;
Evts.text=[Evts.text;['                                            ']];

% new code to copy into other versions of code:
while 1
    if Evts.codes(1,1)==0
        Evts.codes(1,:)=[];
        Evts.times(1)=[];
        Evts.length=Evts.length-1;
        Evts.text(1,:)=[];
    else
        break
    end
end


if exist('Evts2','var')
while 1
    if Evts2.codes(1,1)==0
        Evts2.codes(1,:)=[];
        Evts2.times(1)=[];
        Evts2.length=Evts.length-1;
        Evts2.text(1,:)=[];
    else
        break
    end
end
end
if exist('Evts2','var')&&exist('Evts','var')
    %combine events lists
    Evts_backup1=Evts;
    Evts.times=[Evts.times;Evts2.times];
    Evts.codes=[Evts.codes;Evts2.codes];
    Evts.text=[Evts.text;Evts2.text];
    Evts.length=Evts.length+Evts2.length;
end


%delete second row of zeros if there are two concurrent rows of zeros
for i=length(Evts.codes):-1:2
    if Evts.codes(i)==0&&Evts.codes(i-1)==0
        i
        Evts.codes(i,:)=[];
        Evts.times(i)=[];
        Evts.text(i,:)=[];
        Evts.length=Evts.length-1;
        disp('delete_row_of_zeros');
    end
end

%add a row of zeros if there is a missing row of zeros. Assume end of event is the start of the next event (or end of file/Inf). 
for i=length(Evts.codes):-1:1
    if i==length(Evts.codes)
        if Evts.codes(i)>0
            Evts.codes=[Evts.codes;[0 0 0 0]];
            Evts.times=[Evts.times;Inf];
            Evts.text=[Evts.text;['                                            ']];
            Evts.length=Evts.length+1;
            disp('added_row_of_zeros to the end');
        end
    end
    if i>1
        if Evts.codes(i)>0&&Evts.codes(i-1)>0
            i
            Evts.codes=[Evts.codes(1:(i-1),:);[0 0 0 0];Evts.codes(i:end,:)];
            Evts.times=[Evts.times(1:i);Evts.times(i:end)];
            Evts.text=[Evts.text(1:(i-1),:);['                                            '];Evts.text(i:end,:)];
            Evts.length=Evts.length+1;
            disp('added_row_of_zeros');
        end
    end
end

% Make list of event types
clear EventTypeList_text EventTypeList_code tempi
EventTypeList_text=[];
for i=1:2:size(Evts.codes,1)
    codematch=0;
    for j=1:size(EventTypeList_text,1)
        if strcmp(char(EventTypeList_text(j,:)),Evts.text(i,:))
            'found code match'
            codematch=1;
            break
        end
    end
    if codematch==0
        tempi=size(EventTypeList_text,1);
        EventTypeList_code(tempi+1)=Evts.codes(i);
        EventTypeList_text(tempi+1,:)=Evts.text(i,:);
    end
end
EventTypeList_code=EventTypeList_code'
char(EventTypeList_text)

Evts_backup=Evts;

for i=1:size(Evts.codes,1)
    if strncmp(Evts.text(i,:),'AR',2)
        Evts.codes(i,1)=1;
    end
end


%remove zero rows:
Evts.codes(2:2:end,:)=[];

Evts.endtimes=Evts.times(2:2:end);
Evts.starttimes=Evts.times(1:2:end);
Evts.durations=Evts.endtimes-Evts.starttimes;
% Evts.timetoprevarousal=Evts.starttimes-[-Inf;Evts.endtimes(1:end-1)];
% Evts.timetonextarousal=[Evts.starttimes(2:end);Inf]-Evts.starttimes;

% Find sleep state and CPAP level for each event (arousal) %%faster code 2015-04-29
for i=1:length(Evts.starttimes)
    I=round((Evts.starttimes(i)-Time(1))/dt+1);
    Evts.epochs(i,:)=EpochsXHz(I);
    Evts.CPAP(i,:)=Pmask.values(I);
end

% Make arousals events in continuous time
EventsArXHz=zeros(length(Time),1); %Arousals.
for x=1:length(Evts.codes)
    if Evts.codes(x)==1
        lefti=round(Evts.starttimes(x)-a)*Fs+1;
        righti=lefti+round(Evts.durations(x))*Fs;
        if righti==Inf
            righti=length(EventsArXHz);
        end
        %EventsArXHz(Time>=Evts.starttimes(x)&Time<Evts.endtimes(x))=Evts.codes(x); %%slower code because of search
        EventsArXHz(lefti:righti)=Evts.codes(x);
    end
end

% Make respiratory events in continuous time
EventsRespXHz=zeros(length(Time),1); 
for x=1:length(Evts.codes)
    if Evts.codes(x)>1
        lefti=round(Evts.starttimes(x)-a)*Fs+1;
        righti=lefti+round(Evts.durations(x))*Fs;
        if righti==Inf
            righti=length(EventsRespXHz);
        end
        %EventsArXHz(Time>=Evts.starttimes(x)&Time<Evts.endtimes(x))=Evts.codes(x); %%slower code because of search
        EventsRespXHz(lefti:righti)=Evts.codes(x);
    end
end

clear lefti righti

clear EventTypeList_text EventTypeList_code Evts_backup Evts_backup1 Evts2 I backupEvts codematch foundamatch i j n temp tempi x


%% Make all channels the same length as Flow, resample if needed

clear channels_all
for i=4:length(channelnamestemp)
    channels_all{i-3}=eval(['channelnameoptions.' channelnamestemp{i} '{1}']);
end

clear channels
for i=1:length(channels_all)
    if exist(channels_all{i},'var')
        channels{i}=channels_all{i};
    end
end

clear channels_all
%Add code here which resamples if channels are not the same sampling rate
for i=1:length(channels)
    if ~isempty(channels{i})
        dtnew=eval([channels{i} '.interval']);
        if dtnew~=dt
            disp(['resampling: ' channels{i} '.interval = ' num2str(eval([channels{i} '.interval']))]);
            if mod(dt/dtnew,1)~=0 %unfinished coding
                eval([channels{i} '.values = resample(' channels{i} '.values,round(1/dt),round(1/dtnew));' ]);
            else
                %eval([channels{i} '.values = resample(' channels{i} '.values,round(dt/' channels{i} '.interval),1);' ]);
                eval([channels{i} '.values = downsample(' channels{i} '.values,round(dt/dtnew));' ]);
            end
            eval([channels{i} '.interval = dt;']);
        end
    end
end

%Make channels the same length (sometimes these are off by up to 20 samples)
for i=1:length(channels)
    if ~isempty(channels{i})
        while length(eval([channels{i} '.values']))~=length(Time)
            if length(eval([channels{i} '.values']))<length(Time)
                eval([channels{i} '.values(end+1)=' channels{i} '.values(end);']); %add a sample to the end if the channel is too short by 1
            elseif length(eval([channels{i} '.values']))>length(Time)
                eval([channels{i} '.values(end)=[];']); %delete a sample from the end if the channel is too long by 1
            end
        end
    end
end

%Rename channel data: for example, Flow = Flow.values; clear Flow structure variable
for i=1:length(channels)
    if ~isempty(channels{i})
        eval(['temp=' channels{i} '.values;']);
        eval([channels{i} '=temp;']);
        clear temp
    end
end

%channels{length(channels)+1}='Time';
channels{length(channels)+1}='EventsArXHz';
channels{length(channels)+1}='EpochsXHz';

%% Plot preliminary figure to check
figure(1);
set(gcf,'color',[1 1 1]);
nicefig = @() set(gca,'box','off','tickdir','out','xcolor',[1 1 1],'xticklabel',[],'fontsize',8,'fontname','arial narrow');
nicefiglast = @() set(gca,'box','off','tickdir','out','fontsize',8,'fontname','arial narrow');
X=5;
dsf=25;
ax1(1)=subplot(X,1,1); plot(downsample(Time,dsf),downsample(EpochsXHz,dsf),downsample(Time,dsf),downsample(EventsArXHz,dsf),downsample(Time,dsf),downsample(EventsRespXHz,dsf));
ylabel('Hyp/Evnts'); nicefig(); legend('Hyp','Ar','Resp'); legend('boxoff');
ax1(2)=subplot(X,1,2); plot(downsample(Time,dsf),downsample(Flow,dsf)); ylabel('Flow'); nicefig();
ax1(3)=subplot(X,1,3); plot(downsample(Time,dsf),downsample(Pmask,dsf));  ylabel('Pmask'); nicefig();
ax1(4)=subplot(X,1,4); plot(downsample(Time,dsf),downsample(PO2,dsf)); ylabel('PO2'); nicefig();
ax1(6)=subplot(X,1,5); plot(downsample(Time,dsf),downsample(Pes,dsf)); ylabel('Pes'); nicefiglast();

linkaxes(ax1,'x');

%% Position (add import corrections from text here also)

if ~exist('Position')
    Position = 0*Flow;
else
    Position = round(Position*5);
end


%% OSA severity and other PSG information off CPAP
    %filter Pmask
    filter_HFcutoff_butter0 = 1/10; %4Hz=250ms
    filter_order0 = 4;
    [B_butter0,A_butter0] = butter(filter_order0,filter_HFcutoff_butter0/(1/dt/2),'low');
    CPAP = filtfilt(B_butter0,A_butter0,Pmask);
    figure(1);
    ax1(3)=subplot(X,1,3); plot(downsample(Time,dsf),downsample(Pmask,dsf));  ylabel('Pmask'); nicefig();
    hold('on');
    plot(downsample(Time,dsf),downsample(CPAP,dsf),'r');
    
    CPAPrange = [-0.9 0.9];
    CPAPoff = CPAP<CPAPrange(2)&CPAP>CPAPrange(1);
    %plot(downsample(Time,dsf),2*downsample(CPAPoff,dsf),'g');   
    minlength=180; minlengthi = (minlength/dt);
    diffCPAPoff1=diff(CPAPoff);
    CPAPswitch = find(diffCPAPoff1==1|diffCPAPoff1==-1);
    CPAPswitchI = [1;CPAPswitch;length(CPAP)];
    difflengths = diff(CPAPswitchI);
    CPAPoff_ = mod(1+CPAPoff(1)+(1:length(difflengths)),2)'; %
    temp = [difflengths,CPAPoff_];
    Iremove = find(difflengths<minlengthi&CPAPoff_==1);
    for i=1:length(Iremove)
        CPAPoff(CPAPswitchI(Iremove(i)):CPAPswitchI(Iremove(i)+1))=0;
    end
    plot(downsample(Time,dsf),5*downsample(1-CPAPoff,dsf),'g');
    
    Evts.starttimesi=round((Evts.starttimes/dt)+1)
    
    stageoptions = {'EpochsXHz>0&EpochsXHz<=5','EpochsXHz>0&EpochsXHz<5','EpochsXHz==5','EpochsXHz==0','EpochsXHz==1','EpochsXHz==2','EpochsXHz==3'};
    
    %AllSleep,NREM,REM,W,N1,N2,N3
    %Duration|ARindex|ApOindex|ApCindex|HypOindex|Mindex|HypCindex|AHI'
    
    clear DurationandEventsPerHour_Supine_CPAPoff
    for k=1:length(stageoptions)
    criteria = (Position==0|Position==2)&CPAPoff&eval(stageoptions{k});
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
    
    clear Duration_Allpositions_CPAPofforon
    for k=1:length(stageoptions)
    criteria = eval(stageoptions{k});
    Duration_Allpositions_CPAPofforon(k) = sum(criteria)*dt/60;
    end
    
    PercentTSTAllpositions_CPAPofforon = Duration_Allpositions_CPAPofforon/Duration_Allpositions_CPAPofforon(1)*100;  
    
    AllData = [DurationandEventsPerHour_Supine_CPAPoff_  AllPositionsAllStatesDurationAHI  Duration_Allpositions_CPAPofforon PercentTSTAllpositions_CPAPofforon]
    

%% RIP
RIP = exist('ThNoxRIP','var')||exist('VolTh1','var');
if RIP
    if ~exist('VolTh1','var')
VolTh1 = nanmedian(alphaFlow)*(ThNoxRIP-nanmedian(ThNoxRIP));
VolAb1 = nanmedian(alphaFlow)*nanmedian(kFlow)*(AbNoxRIP-nanmedian(AbNoxRIP));
    end
clear AbNoxRIP ThNoxRIP alphaFlow kFlow %don't clear if you want to change these below

% Plot preliminary figure to check
figure(1);
set(gcf,'color',[1 1 1]);
nicefig = @() set(gca,'box','off','tickdir','out','xcolor',[1 1 1],'xticklabel',[],'fontsize',8,'fontname','arial narrow');
nicefiglast = @() set(gca,'box','off','tickdir','out','fontsize',8,'fontname','arial narrow');
X=6;
dsf=50;
ax1(1)=subplot(X,1,1); plot(downsample(Time,dsf),downsample(EpochsXHz,dsf),downsample(Time,dsf),downsample(EventsArXHz,dsf),downsample(Time,dsf),downsample(EventsRespXHz,dsf));
ylabel('Hyp/Evnts'); nicefig();
ax1(2)=subplot(X,1,2); plot(downsample(Time,dsf),downsample(Flow,dsf)); ylabel('Flow'); nicefig();
ax1(3)=subplot(X,1,3); plot(downsample(Time,dsf),downsample(Pmask,dsf));  ylabel('Pmask'); nicefig();
ax1(4)=subplot(X,1,4); plot(downsample(Time,dsf),downsample(VolTh1,dsf)); ylabel('VTh'); nicefig();
ax1(5)=subplot(X,1,5); plot(downsample(Time,dsf),downsample(VolAb1,dsf)); ylabel('VAb'); nicefig();
ax1(6)=subplot(X,1,6); plot(downsample(Time,dsf),downsample(Pes,dsf)); ylabel('Pes'); nicefig();

linkaxes(ax1,'x');
end

%% Filter gently for use in analysis
try %close the preliminary figure to save memory
close (1)
catch me
end

if 1
    filter_HFcutoff_butter0 = 3; %4Hz=250ms
    filter_order0 = 4;
    [B_butter0,A_butter0] = butter(filter_order0,filter_HFcutoff_butter0/(1/dt/2),'low');
    PCO2f = filtfilt(B_butter0,A_butter0,PCO2);
    PO2f = filtfilt(B_butter0,A_butter0,PO2);

    filter_HFcutoff_butter0 = 10;
    filter_order0 = 2;
    [B_butter0,A_butter0] = butter(filter_order0,filter_HFcutoff_butter0/(1/dt/2),'low');
    Flowf = filtfilt(B_butter0,A_butter0,Flow);
    
    filter_HFcutoff_butter0 = 2;
    filter_order0 = 2;
    [B_butter0,A_butter0] = butter(filter_order0,filter_HFcutoff_butter0/(1/dt/2),'low');
    Pesf = filtfilt(B_butter0,A_butter0,Pes);
end

%% Import from Excel
[num,~,~]=xlsread([filename '.xlsx'],1,'B:C');

%% Clear Summary data
delta = 0;
clear VO2summary;

%%  Loop
usefixeddelay=1; %when first running the program set this to zero, then set to 1 and use the average delay for "fixeddelay" ad re-run
    fixeddelay=1.312;
    referenceVO2 = 248.34; %just to plot reference line on graph
    maskdeadspace = 0.00;%maskdeadspace = 0.09;
    pause_to_assess_figures_1=0;
    pause_to_assess_figures_2=0;
    Pesdelay = 0.3;
runall=1;

if runall
    m=1; %start from
    m_max=size(num,1); %%finish
else
    m=31; %<----to run single window change this number and make runall=0; %8, 26
    m_max=m; 
end

%%
while m<=m_max 

L = num(m,1)-delta;
R = num(m,2)+delta;
Li = round((L-Time(1))/dt+1);
Ri = round((R-Time(1))/dt+1);

Flow_ = Flow(Li:Ri);
if exist('Position')
    Position_ = Position(Li:Ri);
else
    Position_ = 0*Flow_;
end
Time_ = Time(Li:Ri);
Pmask_ = Pmask(Li:Ri);
if exist('Pes')
    Pesdelayi = round((Pesdelay/dt));
    Pes_ = Pes(Li-Pesdelayi:Ri-Pesdelayi);
    Pesf_ = Pesf(Li-Pesdelayi:Ri-Pesdelayi);
else
    Pes_ = 0*Flow_;
    Pesf_ = Pes_;
end
EventsRespXHz_ = EventsRespXHz(Li:Ri);
EventsArXHz_ = EventsArXHz(Li:Ri);
EpochsXHz_ = EpochsXHz(Li:Ri);
SaO2_=round(SaO2(Li:Ri));
if exist('VolTh1','var')
    VolRIP_ = VolTh1(Li:Ri)+VolAb1(Li:Ri);
end

    %% Test 0.2 s time constant on flow
    filter_HFcutoff_butter0 = 1/(2*pi*0.1);
    filter_order0 = 1;
    [B_butter0,A_butter0] = butter(filter_order0,filter_HFcutoff_butter0/(1/dt/2),'low');
    Flow_timeconst = filter(B_butter0,A_butter0,Flow_);
    
%% Flow analysis: Calculate end-exp times
useactualzero=0;
if 0
[I,Vdot,VT,Ti,Te,VTi,VTe,FlowC_]=Vflowanalysis2(Flow_,Time_,dt,0,useactualzero,2);
[I3,~,~,~,~,VTi2,VTe2,~]=Vflowanalysis2(Flow_,Time_,dt,0,0,1);
else
[I,Vdot,VT,Ti,Te,VTi,VTe,FlowC_]=Vflowanalysis2(Flow_timeconst,Time_,dt,0,useactualzero,2);
[~,~,~,~,~,VTi2,VTe2,~]=Vflowanalysis2(Flow_timeconst,Time_,dt,0,0,1);  
[I2,Vdotunfiltered,VTunfiltered,Tiunfiltered,Teunfiltered,VTiunfiltered,VTeunfiltered]=Vflowanalysis_knownI(Flow_-mean(Flow_),Time_,dt,0,I);
end
CPAPlevel_est = median([Pmask_(I2.starti);Pmask_(I2.midi)]);
%     if CPAPlevel_est>1
%         useactualzero=0;
%         [I,Vdot,VT,Ti,Te,VTi,VTe]=Vflowanalysis2(Flow_,Time_,dt,0,useactualzero,);
%         CPAPlevel_est = median([Pmask_(I.starti);Pmask_(I.midi)]);
%     end
%% Pes analysis
if exist('Pes')
%find Pcw
voltemp1 = cumsum(Flow_)*dt;
    vollower = interp1(Time_(I2.starti),voltemp1(I2.starti),Time_,'linear');
    vollower(1:I2.starti(1)-1)=voltemp1(I2.starti(1)); %extrapolation using nearest
    vollower(I2.starti(end):end)=voltemp1(I2.starti(end)); %extrapolation using nearest
voltemp2 = voltemp1 - vollower;
Pcw = 5*voltemp2;

% peslower = interp1(Time_(I2.starti),Pesf_(I2.starti),Time_,'linear');
%     peslower(1:I2.starti(1)-1)=Pesf_(I2.starti(1)); %extrapolation using nearest
%     peslower(I2.starti(end):end)=Pesf_(I2.starti(end)); %extrapolation using nearest
%     
%     temp = Pesf_(I2.starti); Pes_baselineest = median(temp);
%     figure(100); plot(Time_,Pesf_,Time_,peslower)
    
    
    %%
    pcrtiles = 0:1:100;
    PesPrctiles = prctile(Pesf_,pcrtiles);
    IX=find(Pesf_(1)<=PesPrctiles,1);
    if 0
        [Vdot_intended,parameters,rsquared]=PesToVflowFit(Flow_,Pesf_,Time_,5);
    else
        parameters = [10 20 0 0 100-pcrtiles(IX)]
    end

[Vdot_intended]=PesToVflowRun(Flow_,Pesf_,Time_,Inf,parameters);
    [Ipes,~,~,Tipes,Tepes,~,~,~]=Vflowanalysis2(Vdot_intended,Time_,dt,1,useactualzero,1);
    figure(42) %Pes breath/baseline detection using low pass filtered trace
    ax42(1)=subplot(2,1,1); plot(Time_,Vdot_intended); hold('on'); 
        plot(Time_(Ipes.starti-1),Vdot_intended(Ipes.starti-1),'r.'); 
        plot(Time_(Ipes.midi),Vdot_intended(Ipes.midi),'k.'); hold('off'); 
    ax42(2)=subplot(2,1,2); plot(Time_,Pesf_); hold('on'); 
        plot(Time_(Ipes.starti-1),Pesf_(Ipes.starti-1),'r.'); 
        plot(Time_(Ipes.midi),Pesf_(Ipes.midi),'k.'); hold('off'); 
    linkaxes(ax42,'x');
    Pes_baselineest = median(Pesf_(Ipes.starti));
    
Pmus = Pes_-Pcw-Pes_baselineest;
%find baseline Pes
%figure(456); 
%plot(Time_,Pes_,Time_,0*Time_+Pes_baselineest);
clear minPesB maxPesB Pesdelta Pmusdelta PTP1B PTP2B Pesdelta2 Pmusdelta2
%Pmusdelta2 uses inspiratory value rather than a constant baseline as a reference point
for i=1:length(I.starti)
    minPesB(i) = min(Pes_(I2.starti(i):I2.endi(i))-Pes_baselineest);
    maxPesB(i) = max(Pes_(I2.starti(i):I2.endi(i))-Pes_baselineest);
    Pesdelta(i) = max(Pes_(I2.starti(i):I2.endi(i)))-min(Pes_(I2.starti(i):I2.endi(i)));
    Pmusdelta(i) = max(Pmus(I2.starti(i):I2.endi(i)))-min(Pmus(I2.starti(i):I2.endi(i)));
    PTP1B(i) = -sum(Pmus(I2.starti(i):I2.midi(i)))*dt;
    temp=Pmus(I2.starti(i):I2.endi(i)); temp(temp>0)=0;
    PTP2B(i) = -sum(temp)*dt;
end

clear minPesB2 maxPesB2 Pesdelta2 Pmusdelta2 PTP1B2 PTP2B2
for i=1:length(Ipes.starti)
    minPesB2(i) = min(Pes_(Ipes.starti(i):Ipes.endi(i))-Pes_baselineest);
    maxPesB2(i) = max(Pes_(Ipes.starti(i):Ipes.endi(i))-Pes_baselineest);
    Pesdelta2(i) = Pes_(Ipes.starti(i))-min(Pes_(Ipes.starti(i):Ipes.endi(i)));
    Pmusdelta2(i) = Pmus(Ipes.starti(i))-min(Pmus(Ipes.starti(i):Ipes.endi(i)));
    temp = Pmus(Ipes.starti(i):Ipes.midi(i)); temp=temp-temp(1);
    PTP1B2(i) = -sum(temp)*dt;
    temp = Pmus(Ipes.starti(i):Ipes.endi(i)); temp=temp-temp(1); temp(temp>0)=0;
    PTP2B2(i) = -sum(temp)*dt;
end

temp = maxPesB>abs(minPesB); %test for artifact: upwards swing from "baseline" is larger than the negative swing from the "baseline"
minPesB(temp)=NaN;
maxPesB(temp)=NaN;
Pmusdelta(temp)=NaN;
Pesdelta(temp)=NaN;
PTP1B(temp)=NaN;
PTP2B(temp)=NaN;

temp2 = maxPesB2>abs(minPesB2); %test for artifact: upwards swing from "baseline" is larger than the downwards swing from the "baseline"
minPesB2(temp2)=NaN;
maxPesB2(temp2)=NaN;
Pesdelta2(temp2)=NaN;
Pmusdelta2(temp2)=NaN;
PTP1B2(temp2)=NaN;
PTP2B2(temp2)=NaN;

Ttot = Ti+Te;
Ttotpes = Tipes+Tepes;

PTP1=60*sum(PTP1B(~isnan(PTP1B)))/sum(Ttot(~isnan(PTP1B)));
PTP2=60*sum(PTP2B(~isnan(PTP2B)))/sum(Ttot(~isnan(PTP2B)));

PTP12=60*sum(PTP1B2(~isnan(PTP1B2)))/sum(Ttotpes(~isnan(PTP1B2)));
PTP22=60*sum(PTP2B2(~isnan(PTP2B2)))/sum(Ttotpes(~isnan(PTP2B2)));

figure(898);
plot(Time_,Pes_-Pes_baselineest,Time_,Pcw,Time_,Pmus);

%Pmus_mean = nanmean(Pmusdelta);
Pmus_mean=sum(Pmusdelta(~isnan(Pmusdelta)).*Ttot(~isnan(Pmusdelta)))/sum(Ttot(~isnan(Pmusdelta)));
Pes_mean=sum(Pesdelta(~isnan(Pesdelta)).*Ttot(~isnan(Pesdelta)))/sum(Ttot(~isnan(Pesdelta)));

Pmus_mean2=sum(Pmusdelta2(~isnan(Pmusdelta2)).*Ttotpes(~isnan(Pmusdelta2)))/sum(Ttotpes(~isnan(Pmusdelta2)));
Pes_mean2=sum(Pesdelta2(~isnan(Pesdelta2)).*Ttotpes(~isnan(Pesdelta2)))/sum(Ttotpes(~isnan(Pesdelta2)));



%% arousal threshold estimation
%first need an arousal signal breath-to-breath
Farousalinspthres = 0.1;
Minlength=3;

%default values
PesAtArmedian = NaN;
PmusAtArmedian = NaN;
PesAtArmedianN = 0;

for k=1:1 %enables break / goto

clear Ar_
for i=1:length(Ipes.starti)
Ar_(i) = mean(EventsArXHz_(Ipes.starti(i):Ipes.midi(i))); %using proportion of inspiration...
end

Ar2_ = Ar_;
Ar2_(Ar_>Farousalinspthres)=1;
Ar2_(Ar_<=Farousalinspthres)=0;

ArOnsetNeg1  = [diff(Ar2_) NaN];
ISleepOn = find(ArOnsetNeg1==-1);
ISleepOff = find(ArOnsetNeg1==1);
if isempty(ISleepOn)||isempty(ISleepOff)
    break
end
    
if ISleepOff(1)<ISleepOn(1)
    ISleepOn = [1 ISleepOn];
end
if ISleepOff(end)<ISleepOn(end)
    ISleepOn(end) = [];
end

LengthSleepBr = ISleepOff-ISleepOn;

ISleepOff(LengthSleepBr<Minlength)=[];

if isempty(ISleepOn)||isempty(ISleepOff)
    break
end

clear PesAtAr PmusAtAr previousorcurrent
for i=1:length(ISleepOff)
    rangetemp = [ISleepOff(i)-1 ISleepOff(i)];
    [PesAtAr(i) previousorcurrent(i)] = max(Pesdelta2(rangetemp));  
    PmusAtAr(i) = max(Pmusdelta2(rangetemp));  
end


figure(56)
ax56(1)=subplot(3,1,1); plot(Time_,EventsArXHz_,'r--'); hold('on');
ax56(2)=subplot(3,1,2); plot(Time_,Pes_,'k'); hold('on');
for i=1:length(PesAtAr)
    itemp = ISleepOff(i)-2+previousorcurrent(i);
    plot([Time_(Ipes.starti(itemp)) Time_(Ipes.endi(itemp))],[Pes_(Ipes.starti(itemp)) Pes_(Ipes.starti(itemp))],'g');
    plot([Time_(Ipes.starti(itemp)) Time_(Ipes.endi(itemp))],-PesAtAr(i)+[Pes_(Ipes.starti(itemp)) Pes_(Ipes.starti(itemp))],'g');
    plot([Time_(Ipes.starti(itemp)) Time_(Ipes.starti(itemp))],[Pes_(Ipes.starti(itemp))-PesAtAr(i) Pes_(Ipes.starti(itemp))],'g');
    plot([Time_(Ipes.endi(itemp)) Time_(Ipes.endi(itemp))],[Pes_(Ipes.starti(itemp))-PesAtAr(i) Pes_(Ipes.starti(itemp))],'g');
end

ax56(3)=subplot(3,1,3); plot(Time_,Flow_,'k'); hold('on');
linkaxes(ax56,'x');


PesAtArmedian = nanmedian(PesAtAr);
PmusAtArmedian = nanmedian(PmusAtAr);
PesAtArmedianN = length(PesAtAr);

end

disp(['arousal threshold: ' num2str(PesAtArmedian)]);

Pmus_mean2_sleep=sum(Pmusdelta2(~isnan(Pmusdelta2)&Ar2_==0).*Ttotpes(~isnan(Pmusdelta2)&Ar2_==0))/sum(Ttotpes(~isnan(Pmusdelta2)&Ar2_==0));
Pes_mean2_sleep=sum(Pesdelta2(~isnan(Pesdelta2)&Ar2_==0).*Ttotpes(~isnan(Pesdelta2)&Ar2_==0))/sum(Ttotpes(~isnan(Pesdelta2)&Ar2_==0));


else %no Pes data available:
    Pmus_mean = NaN;
    PTP1 = NaN;
    PTP2 = NaN;
    PTP12 = NaN;
    PTP22 = NaN;
    Pes_mean= NaN;
    Pmus_mean2= NaN;
    Pes_mean2= NaN;
    Pmus_mean2_sleep= NaN;
    Pes_mean2_sleep= NaN;
    PesAtArmedian= NaN;
    PmusAtArmedian= NaN;
    PesAtArmedianN= NaN;
end
    

%% RIP
if exist('VolTh1','var')
    figure(50)
        plot(Time_,VolRIP_,Time_(I.starti),VolRIP_(I.starti),'r.')

%find vol for each breath:
vollower = interp1(Time_(I.starti),VolRIP_(I.starti),Time_,'linear');
    vollower(1:I.starti(1)-1)=VolRIP_(I.starti(1)); %extrapolation using nearest
    vollower(I.starti(end):end)=VolRIP_(I.starti(end)); %extrapolation using nearest
plot(Time_,VolRIP_,Time_(I.starti),VolRIP_(I.starti),'r.',Time_,vollower)
  
  
VTRIP = VT*0;
for i=1:length(VTRIP)
    range = I.starti(i):I.endi(i);
    VTRIP(i)=max(VolRIP_(range)-vollower(range));
end
VTRIP(VTRIP<0)=0;
VdotRIP = VTRIP./(Tiunfiltered+Teunfiltered);

figure(51)
plot(Vdotunfiltered,VdotRIP,'.',[0 0.5],[0 0.5],'k:'); xlabel('VT'); ylabel('VT RIP');
end
%% find PETO2 and PO2 delays... 1
PO2lag_estimated = 1;
PO2lag_estimatedi = -round(PO2lag_estimated/dt);
Li2 = Li-PO2lag_estimatedi;
Ri2 = Ri-PO2lag_estimatedi;
PO2_=PO2(Li2:Ri2);
PO2f_=PO2f(Li2:Ri2);
PO2_B = PO2f_(I.starti);

%% Find nadir PO2 within each breath (e.g. PETO2)
PO2_Bmin = 0*PO2_B;
PO2_Bmini = 0*PO2_B;
for i=1:length(PO2_B)
    if i<length(PO2_B)
        [PO2_Bmin(i),PO2_Bmini(i)] = min(PO2f_(I.midi(i):I.midi(i+1)));
        PO2_Bmini(i)=PO2_Bmini(i)+I.midi(i)-1;
    else
        [PO2_Bmin(i),PO2_Bmini(i)] = min(PO2f_(I.midi(i):I.endi(i)));
        PO2_Bmini(i)=PO2_Bmini(i)+I.midi(i)-1;
    end
end

figure(2); set(gcf,'color',[1 1 1]);
ax2(1)=subplot(3,1,1); plot(Time_,EpochsXHz_,Time_,EventsArXHz_,Time_,EventsRespXHz_); nicefig(); ylabel('Hyp/Events'); legend('Hyp','Ar','Resp'); legend('boxoff');
ax2(2)=subplot(3,1,2); plot(Time_,Flow_,Time_,FlowC_,Time_,Flow_timeconst,Time_(I.starti),FlowC_(I.starti),'r.',Time_(I.midi),FlowC_(I.midi),'k.'); nicefig(); ylabel('Flow'); legend('original','filtered'); legend('boxoff'); 
ax2(3)=subplot(3,1,3); plot(Time_,PO2_,Time_,PO2f_,Time_(I.starti),PO2f_(I.starti),'r.',Time_(I.midi),PO2f_(I.midi),'k.',Time_(PO2_Bmini),PO2f_(PO2_Bmini),'g.'); nicefiglast(); ylabel('PO2'); legend('original','filtered'); legend('boxoff'); 
 nicefiglast();
linkaxes(ax2,'x');

%% Shift end-tidal to peak value
if 0
    mediandeltai_i=0;
    for i=1:3
        PO2_Bi=I.starti; %location of start of current insp
        
        for j=1:length(PO2_B)
            for k=-round(1/dt)*0.5:round(1/dt)*4 %-round(1/dt)*0.5:round(1/dt)*3.5
                range=(PO2_Bi(j):PO2_Bi(j)+1); %
                if range(1)==1
                    break;
                end
                slope = diff(PO2f_(range)); %mean(gradient(PO2(range)));
                slopeplus1 = diff(PO2f_(range-1));
                slopeless1 = diff(PO2f_(range+1));
                if slope>0
                    direction = -1;
                else
                    direction = 1;
                end
                %if absolute value of gradient is at a local minimum and
                if ((abs(slope)<abs(slopeplus1))&&(abs(slope)<abs(slopeless1)))&&((slope<slopeless1)&&(slope>slopeplus1))
                    break;
                else
                    PO2_Bi(j)=PO2_Bi(j)+direction; %moving towards the peak
                end
            end
        end
        
        PO2_B=PO2f_(PO2_Bi);
        
        if 0
            figure(12); plot(Time_,PO2f_,Time_(I.starti),PO2f_(I.starti),'k.',Time_(PO2_Bi),PO2_B,'r.');
        end
        
        %look at cumulative probability distribution of delays
        if 1
            figure(11); plot(-sum(mediandeltai_i)*dt+PO2lag_estimated-prctile(I.starti-PO2_Bi,1:99)*dt,1:99);
        end
        
        mediandeltai_i(i) = round(median(I.starti-PO2_Bi)); %should be mid?
        %use value at the median shift.
        if 1
            Li2 = Li2-mediandeltai_i(i);
            Ri2 = Ri2-mediandeltai_i(i);
            PO2_=PO2(Li2:Ri2);
            PO2f_=PO2f(Li2:Ri2);
            PO2_B = PO2f_(I.starti);
            PO2_Bi_temp=PO2_Bi;
            PO2_Bi = I.starti; %reset this...
        end
        Vflow_to_PO2_delay_measured = (Li2-Li)*dt;
        
    end %do this twice...
else
    if ~usefixeddelay
    delaytempi = I.endi-PO2_Bmini;
    mediandeltai_i = median(delaytempi);
    else
    mediandeltai_i = -round((fixeddelay-PO2lag_estimated)/dt);    
    end
    Li2 = Li2-mediandeltai_i;
    Ri2 = Ri2-mediandeltai_i;
    PO2_=PO2(Li2:Ri2);
    PO2f_=PO2f(Li2:Ri2);
    PO2_B = PO2f_(I.starti);
%     PO2_Bi_temp=PO2_Bi;
    PO2_Bi = I.starti; %reset this...
    Vflow_to_PO2_delay_measured = (Li2-Li)*dt;
end

%% Find min PO2f
PO2_Bmin = 0*PO2_B;
PO2_Bmini = 0*PO2_B;
PIO2_Bmax = 0*PO2_B;
PIO2_Bmaxi = 0*PO2_B;
for i=1:length(PO2_B)
    if i<length(PO2_B)
        [PO2_Bmin(i),PO2_Bmini(i)] = min(PO2f_(I.midi(i):I.midi(i+1)));
        PO2_Bmini(i)=PO2_Bmini(i)+I.midi(i)-1;
    else
        [PO2_Bmin(i),PO2_Bmini(i)] = min(PO2f_(I.midi(i):I.endi(i)));
        PO2_Bmini(i)=PO2_Bmini(i)+I.midi(i)-1;
    end
    [PIO2_Bmax(i),PIO2_Bmaxi(i)] = max(PO2f_(I.starti(i):I.endi(i)));
    PIO2_Bmaxi(i)=PIO2_Bmaxi(i)+I.starti(i)-1;
end

%%

figure(2); set(gcf,'color',[1 1 1]);
ax2(1)=subplot(3,1,1); plot(Time_,EpochsXHz_,Time_,EventsArXHz_,Time_,EventsRespXHz_); nicefig(); ylabel('Hyp/Events'); legend('Hyp','Ar','Resp'); legend('boxoff');
ax2(2)=subplot(3,1,2); plot(Time_,Flow_,Time_,FlowC_,Time_(I.starti),FlowC_(I.starti),'r.',Time_(I.midi),FlowC_(I.midi),'k.'); nicefig(); ylabel('Flow'); legend('original','filtered'); legend('boxoff'); 
ax2(3)=subplot(3,1,3); plot(Time_,PO2_,Time_,PO2f_,Time_(PO2_Bi),PO2f_(PO2_Bi),'r.',Time_(I.midi),PO2f_(I.midi),'k.',Time_(PO2_Bmini),PO2f_(PO2_Bmini),'g.',Time_(PIO2_Bmaxi),PO2f_(PIO2_Bmaxi),'b.'); nicefiglast(); ylabel('PO2'); legend('original','filtered','end exp (ET)','end insp','min'); legend('boxoff'); nicefiglast();
linkaxes(ax2,'x');

PIO2_B=PO2f_(I.midi); %should be start?


%% Find CO2-to-O2 time lag and select PCO2 range
if 1
    Li3 = Li2-20;
    Ri3 = Ri2-20;
    PCO2f_=PCO2f(Li3:Ri3);
    numLags=round(1/dt*4); %crosscorr looks both directions by numLags
    %skip first and last N elements as might be NaNs:
    Nskip1=500;
    Nskip=500;
    [xcf,lags,bounds] = crosscorr(PO2f_(Nskip1:(length(PO2f_)-Nskip)),PCO2f_(Nskip1:(length(PCO2f_)-Nskip)),numLags);
    figure(16), plot(lags,xcf); ylabel('correlation coefficient PCO2 vs PO2'); xlabel('PCO2 lags PO2 (samples)');
    [tempval,tempi]=min(xcf); %test: tempval is the same as xcf(tempi+numLags);
    tempi = tempi-Nskip1-1;
    %Adjust to match CO2 by tempi
    Li3 = Li3+tempi;
    Ri3 = Ri3+tempi;
    PCO2_=PCO2(Li3:Ri3);
    PCO2f_=PCO2f(Li3:Ri3);    
    PCO2_B = PCO2f_(I.starti);
end

%% Find max PCO2f
PCO2_Bmax = 0*PCO2_B;
PCO2_Bmaxi = 0*PCO2_B;
PICO2_Bmin = 0*PCO2_B;
PICO2_Bmini = 0*PCO2_B;
for i=1:length(PCO2_B)
    if i<length(PCO2_B)
        [PCO2_Bmax(i),PCO2_Bmaxi(i)] = max(PCO2f_(I.midi(i):I.midi(i+1)));
        PCO2_Bmaxi(i)=PCO2_Bmaxi(i)+I.midi(i)-1;
    else
        [PCO2_Bmax(i),PCO2_Bmaxi(i)] = max(PCO2f_(I.midi(i):I.endi(i)));
        PCO2_Bmaxi(i)=PCO2_Bmaxi(i)+I.midi(i)-1;
    end
    [PICO2_Bmin(i),PICO2_Bmini(i)] = min(PCO2f_(I.starti(i):I.endi(i)));
    PICO2_Bmini(i)=PICO2_Bmini(i)+I.starti(i)-1;
end


    
%% Figure
    figure(2); set(gcf,'color',[1 1 1]); X=4;
    ax2(1)=subplot(X,1,1); plot(Time_,EpochsXHz_,Time_,EventsArXHz_,Time_,EventsRespXHz_); nicefig(); ylabel('Hyp/Events'); legend('Hyp','Ar','Resp'); legend('boxoff');
    ax2(2)=subplot(X,1,2); plot(Time_,Flow_,Time_,FlowC_,Time_,Flow_timeconst,Time_(I.starti),FlowC_(I.starti),'r.',Time_(I.midi),FlowC_(I.midi),'k.'); nicefig(); ylabel('Flow'); legend('original','filtered'); legend('boxoff'); 
    ax2(3)=subplot(X,1,3); plot(Time_,PO2_,Time_,PO2f_,Time_(PO2_Bi),PO2f_(PO2_Bi),'r.',Time_(I.midi),PO2f_(I.midi),'k.',Time_(PO2_Bmini),PO2f_(PO2_Bmini),'g.',Time_(PIO2_Bmaxi),PO2f_(PIO2_Bmaxi),'b.'); nicefig; ylabel('PO2'); legend('original','filtered','end exp (ET)','end insp','min'); legend('boxoff'); nicefiglast();
    ax2(4)=subplot(X,1,4); plot(Time_,PCO2_,Time_,PCO2f_,Time_(I.endi),PCO2f_(I.endi),'r.',Time_(I.midi),PCO2f_(I.midi),'k.'); nicefiglast(); ylabel('PCO2'); legend('original','filtered','end exp (ET)','end insp','max'); legend('boxoff'); 
    
linkaxes(ax2,'x');

    
%% Simulate a PO2 waveform

PO2model=0*PO2f_+150;

if 0
PO2model=0*PO2f_+150;
deadspace1 = 0.13;
deadspace2 = 0.13;
PIsignal = 146;
for i=1:length(I.starti)
    VoltempInsp = cumsum(Flow_(I.starti(i):I.midi(i))-1)*dt;
    VoltempExp = -cumsum(Flow_(I.midi(i):I.endi(i)))*dt;
    Itemp=find(VoltempExp>deadspace1,1);
    Itemp2=find(VoltempInsp>deadspace2,1);
    if i>1
        PO2model(I.starti(i):(I.starti(i)+Itemp2-1))=PO2f_(I.endi(i-1));
    end
    PO2model((I.midi(i)+Itemp):(I.endi(i)-1))=PO2f_(I.endi(i));
    PO2model((I.starti(i)+Itemp2):(I.midi(i)+Itemp))=PIsignal;
end

    figure(2)
    plot(Time_,PO2f_,Time_,PO2model,Time_(I.endi),PO2f_(I.endi),'r.',Time_(I.midi),PO2f_(I.midi),'k.'); 

    HF1 = 2; HF2=5;
    PO2modelf=PO2model;
    for i=1:5
    [B_butterX,A_butterX] = butter(1,HF1/(1/dt/2),'low');
    PO2modelf = filter(B_butterX,A_butterX,PO2modelf);
    end
    
    [B_butterY,A_butterY] = butter(1,HF2/(1/dt/2),'low');
    PO2modelf = filtfilt(B_butterY,A_butterY,PO2modelf);
    ax3(2)=subplot(3,1,2); plot(Time_,PO2_,Time_,PO2model,Time_,PO2modelf,Time_(I.endi),PO2_(I.endi),'r.',Time_(I.midi),PO2_(I.midi),'k.'); 
    
    VI_rs=Time_*NaN;
    VI_rs(1:I.starti(1)-1)=Vdot(1);
    for i=I.starti(1):length(Time_);
        VI_rs(i) = Vdot(find(Time_(I.starti)<=Time_(i),1,'last')); 
    end
    
    ax3(3)=subplot(3,1,3); plot(Time_,VI_rs,Time_,Flow_/4); 
    linkaxes(ax3,'x');
    
    VO2modelf = FlowC_.*PO2model;
end
    
    
%%
    for i=1:length(I.starti)
        VoltempInsp = cumsum(FlowC_(I.starti(i):I.midi(i))-1)*dt;
        VoltempExp = -cumsum(FlowC_(I.midi(i):I.endi(i)))*dt;
        Itemp=find(VoltempExp>0.15,1);
        PO2model((I.midi(i)+Itemp):(I.endi(i)-1))=PO2f_(I.endi(i));
        PO2model(I.starti(i):(I.midi(i)+Itemp))=150;
    end


%% Estimate continuous PO2/PCO2 signals (plant gain)
    VTmean = sum(VT.*Ttot)/sum(Ttot);

%remove small breaths
    PCO2_Bmaxi_incl = PCO2_Bmaxi;
    PO2_Bmini_incl = PO2_Bmini;
    criteriatemp = VT<0.100;
    
    PCO2_Bmaxi_incl(criteriatemp) = [];
    PO2_Bmini_incl(criteriatemp) = [];

%remove breaths with unusually high PO2 / low PCO2
    medianPETO2 = median(PO2f_(PO2_Bmini_incl));
    medianPETCO2 = median(PCO2f_(PCO2_Bmaxi_incl));
    deltathres = 15;
    
    criteriatemp = (PO2f_(PO2_Bmini_incl)>medianPETO2+deltathres) & (PCO2f_(PCO2_Bmaxi_incl)<medianPETCO2-deltathres);
    PCO2_Bmaxi_incl(criteriatemp) = [];
    PO2_Bmini_incl(criteriatemp) = [];
    
    
%find continuous PO2/PCO2
    [X2,X2_O2,FVAL,FVAL_O2,meanPCO2,minPCO2,maxPCO2,stdPCO2,meanPO2,minPO2,maxPO2,stdPO2,PCO2est_win,PO2est_win,time_dt_win,figurehandle]...
    =PlantgainO2andCO2(Flow_,PCO2f_,PO2f_,Vdot,VTe,Time_,...
     Ttot,I.starti,Time_(PO2_Bmini_incl),PCO2f_(PCO2_Bmaxi_incl),PO2f_(PO2_Bmini_incl),I.midi);

    

%use values of interpolated PO2/PCO2 for excluded ETO2/ETCO2 data are not available for Xsec
    criteriasummary = VT'<0.100|((PO2f_(PO2_Bmini)>medianPETO2+deltathres) & (PCO2f_(PCO2_Bmaxi)<medianPETCO2-deltathres));

    
    PO2_Bmin_model = PO2_Bmin;
    PCO2_Bmax_model = PCO2_Bmax;
    for i=1:length(criteriasummary)
        if criteriasummary(i)==1
            PCO2_Bmax_model(i)=interp1(time_dt_win,PCO2est_win,Time_(I.endi(i)),'linear'); 
            PO2_Bmin_model(i)=interp1(time_dt_win,PO2est_win,Time_(I.endi(i)),'linear');
        else
            PCO2_Bmax_model(i)=interp1(time_dt_win,PCO2est_win,Time_(PCO2_Bmaxi(i)),'linear'); 
            PO2_Bmin_model(i)=interp1(time_dt_win,PO2est_win,Time_(PO2_Bmini(i)),'linear');
        end
    end
    
    PO2_Bmin_modelcombo = PO2_Bmin;
    PCO2_Bmax_modelcombo = PCO2_Bmax;
    for i=1:length(criteriasummary)
        if (criteriasummary(i)==1||(VT(i)<VTmean&PO2_Bmin_model(i)<PO2_Bmin(i)))&&(PO2_Bmin_model(i)<PO2_Bmin(i))
            PCO2_Bmax_modelcombo(i)=interp1(time_dt_win,PCO2est_win,Time_(I.endi(i)),'linear'); 
            PO2_Bmini(i) = I.endi(i);
            PCO2_Bmaxi(i) = I.endi(i);
            PO2_Bmin_modelcombo(i)=interp1(time_dt_win,PO2est_win,Time_(I.endi(i)),'linear');
            criteriasummary(i)=1;
        else
            criteriasummary(i)=0;
        end           
    end
    

 
figure(345);
set(gcf,'color',[1 1 1]);

    ax345(1)=subplot(3,1,1);plot(Time_,PO2f_,'k',time_dt_win,PO2est_win,'r--',...
        Time_(PO2_Bmini(criteriasummary==0)),PO2_Bmin_model(criteriasummary==0),'r.',...
        Time_(PO2_Bmini(criteriasummary==0)),PO2_Bmin_modelcombo(criteriasummary==0),'k.');
        hold('on')
        plot(Time_(PO2_Bmini(criteriasummary==1)),PO2_Bmin_modelcombo(criteriasummary==1),'ro','markersize',2);
        plot(Time_(PO2_Bmini),PO2_Bmin_modelcombo,'g--');
        plot(Time_(PO2_Bmini),PO2_Bmin,'k:');
        ylim([75 150]);
        hold('off');
        box('off');
 
    %time_dt_win2 = time_dt_win + Time_(PO2_Bmini_incl(1));
    ax345(2)=subplot(3,1,2);plot(Time_,PCO2f_,'k',time_dt_win,PCO2est_win,'r--',...
        Time_(PCO2_Bmaxi(criteriasummary==0)),PCO2_Bmax_model(criteriasummary==0),'r.',...
        Time_(PCO2_Bmaxi(criteriasummary==0)),PCO2_Bmax_modelcombo(criteriasummary==0),'k.')
        hold('on')
        plot(Time_(PCO2_Bmaxi(criteriasummary==1)),PCO2_Bmax_modelcombo(criteriasummary==1),'ro','markersize',2);
        plot(Time_(PCO2_Bmaxi),PCO2_Bmax_modelcombo,'g--'); 
        plot(Time_(PCO2_Bmaxi),PCO2_Bmax,'k:'); 
        ylim([0 50]);
        hold('off');
        box('off');
    ax345(3)=subplot(3,1,3);plot(Time_,Flow_,'k');
    ylim([-2 2]);
    box('off');
    linkaxes(ax345,'x');
    
    tempstr = [];
    for i=1:3-length(num2str(m))
        tempstr = ['0' tempstr];
    end
    saveas(figure(345),[directory 'figures\' filename '_modelPO2_' tempstr num2str(m)]);
    
    if pause_to_assess_figures_1
        pause
    end
    
%% breath-tobreath arousal info for ventilation during sleep

clear Ar_
for i=1:length(I2.starti) %I2 holds the unfiltered flow analysis indices
Ar_(i) = mean(EventsArXHz_(I2.starti(i):I2.midi(i))); %using proportion of inspiration...
end
Ar3_ = Ar_;
Ar3_(Ar_>Farousalinspthres)=1;
Ar3_(Ar_<=Farousalinspthres)=0;

%% Estimate metabolic rate
    PO2_Bmin_shifted = [PO2f_(I.starti(1));PO2_Bmin(1:end-1)];
    PCO2_Bmax_shifted = [PCO2f_(I.starti(1));PCO2_Bmax(1:end-1)];
    AmbTemp = 295;
    AmbHumid = 0.5;
    PvapAmb = (10^(8.07131-1730.63/(233.426+AmbTemp-273.15)));
    
    
    
    %FO2ATPH = FO2dry*(760-PvapAmb*AmbHumid)/760 
    FO2drytoATPHfactor = (760-PvapAmb*AmbHumid)/760;
    FO2drytoBTPSfactor = (760-47*1)/760;
    %Flow_ = Flow_-mean(Flow_(I.starti(i):I.endi(end)))
    
    VO2 = FlowC_.*(PO2f_/713)*1000*60;
    VCO2 = FlowC_.*(PCO2f_/713)*1000*60;
    %VO2 = FlowC_.*(PO2model/713)*1000*60;
    %VO2 = FlowC_.*(PO2_/713)*1000*60;
    %VO2insp = FlowC_.*(150/713)*1000*60;
    VTa = VTeunfiltered-0.15-maskdeadspace;
    VTa(VTa<0)=0.00;
    %VO2mean1 = 273/310*mean(VO2(I.starti(1):I.endi(end)));
    VO2mean1_FlowC_ = mean(FlowC_(I.starti(1):I.endi(end)));
    
    clear VO2mean1B VO2mean1_FlowC_B VO2mean1Bexp VO2inspB VCO2mean1Bexp VCO2inspB
    volumecorrectionB  = VTunfiltered./VT; 
        volumecorrectionB(isnan(volumecorrectionB))=1; 
        volumecorrectionB(volumecorrectionB==Inf)=1;
        
    for i=1:length(I.starti)
        tempflow = FlowC_(I.starti(i):(I.endi(i)-1));%-mean(FlowC_(I.starti(i):I.endi(i)-1)));
        tempPO2 = (PO2f_(I.starti(i):I.endi(i)-1)/713);
        VO2mean1B(i) = 273/310*mean(VO2(I.starti(i):I.endi(i)))*volumecorrectionB(i);   %273/310*mean(tempflow.*tempPO2)*1000*60; %
        VO2mean1Bexp(i) = -mean(VO2(I.midi(i):I.endi(i)))*volumecorrectionB(i);   %273/310*mean(tempflow.*tempPO2)*1000*60; %
        VCO2mean1Bexp(i) = -mean(VCO2(I.midi(i):I.endi(i)))*volumecorrectionB(i);   %273/310*mean(tempflow.*tempPO2)*1000*60; %
        VO2inspB(i) = -mean(FlowC_(I.midi(i):I.endi(i)).*(PIO2_Bmax(i)/713)*1000*60)*volumecorrectionB(i);   %273/310*mean(tempflow.*tempPO2)*1000*60;
        VCO2inspB(i) = -mean(FlowC_(I.midi(i):I.endi(i)).*(PICO2_Bmin(i)/713)*1000*60)*volumecorrectionB(i);   %273/310*mean(tempflow.*tempPO2)*1000*60;
        VO2mean1_FlowC_B(i) = mean(FlowC_(I.starti(i):I.endi(i)-1))*volumecorrectionB(i);
    end
    maskdeadspacevolumebreathedin = maskdeadspace + 0*VT;
    maskdeadspacevolumebreathedin(VTiunfiltered - maskdeadspace<0)=0;
    VO2mean1B(VO2mean1B<0)=0;
    VO2mean1=sum(VO2mean1B.*(Ti+Te))/sum(Ti+Te);
    
    VO2expdata = 273/310*(VO2inspB.*Te./(Ti+Te) - VO2mean1Bexp.*Te./(Ti+Te) - maskdeadspacevolumebreathedin./(Ti+Te).*((PIO2_Bmax'-PO2_Bmin_shifted')/713)*1000*60);
    VCO2expdata = 273/310*(VCO2inspB.*Te./(Ti+Te) - VCO2mean1Bexp.*Te./(Ti+Te) - maskdeadspacevolumebreathedin./(Ti+Te).*((PICO2_Bmin'-PCO2_Bmax_shifted')/713)*1000*60);
    
    VO2expdatamean1=sum(VO2expdata.*(Ti+Te))/sum(Ti+Te);
    VCO2expdatamean1=-sum(VCO2expdata.*(Ti+Te))/sum(Ti+Te);
    
    
    
    clear mixedexpiredB
    for i=1:length(I.starti)
        mixedexpiredB(i)=sum(FlowC_(I.midi(i):I.endi(i)).*PO2f_(I.midi(i):I.endi(i)))/sum(FlowC_(I.midi(i):I.endi(i)));
    end
    mixedexpiredmean = sum(VTeunfiltered.*mixedexpiredB)/sum(VTeunfiltered);
    
    VTaexp = VTeunfiltered-0.15; VTaexp(VTaexp<0)=0;
    VTdexp = 0.15 + 0*VTeunfiltered; VTdexp(VTaexp==0)=VTeunfiltered(VTaexp==0);
    mixedexpiredB3 = (VTdexp.*PIO2_Bmax' + VTaexp.*PO2_Bmin')./(VTdexp+VTaexp);
    mixedexpiredB3((VTdexp+VTaexp)==0)=0;
    mixedexpiredmean3 = sum(VTeunfiltered.*mixedexpiredB3)/sum(VTeunfiltered);
    
    
    %     VTa = VTeunfiltered-0.15-maskdeadspace;
%     VTa(VTa<0)=0.00;
%     VO2_2 = 273/310*1000*(VTa./(Ti+Te).*60.*(150-PO2_Bmin')/760);
    %corrections
    Flow2_=FlowC_;
        Flow2_(Flow2_>0) = 273/AmbTemp*Flow2_(Flow2_>0);
        Flow2_(Flow2_<=0) = 273/310*Flow2_(Flow2_<=0);
    
    PO2f2_=PO2model;
        PO2f2_(Flow2_>0) = FO2drytoATPHfactor*PO2f2_(Flow2_>0);
        PO2f2_(Flow2_<=0) = FO2drytoBTPSfactor*PO2f2_(Flow2_<=0);
    
    VO22 = Flow2_.*(PO2f2_/713)*1000*60;
    VO2mean2 = mean(VO22(I.starti(i):I.endi(end)));
    
    figure(13)
    plot(Time_,VO2,Time_,Time_*0+VO2mean1);
    plot(Time_,VO22,Time_,Time_*0+VO2mean2);
    
    %mean(FlowC_(I.starti(1):I.endi(end)))
    
    clear VO2_B VO2_B2
    for i=1:length(I.starti)
        VO2_B(i) = mean(VO2((I.starti(i)+1):I.endi(i)));
        VO2_B2(i) = mean(VO22((I.starti(i)+1):I.endi(i)));
    end
    VO2_B_mean2 = mean(VO2_B2);
    VO2_B_mean2=sum(VO2_B2.*(Ti+Te))/sum(Ti+Te);
    
    
    VO2_2 = 273/310*1000*(VTa./(Ti+Te).*60.*(150-PO2_Bmin')/760);
    
    Va = VTa./(Tiunfiltered+Teunfiltered);
    meanVa = 60*sum(Va.*(Tiunfiltered+Teunfiltered))/sum(Tiunfiltered+Teunfiltered);
    meanVdot = 60*sum(Vdotunfiltered.*(Tiunfiltered+Teunfiltered))/sum(Tiunfiltered+Teunfiltered);
    
    meanVdot_sleep = 60*sum(Vdotunfiltered(Ar3_==0).*(Tiunfiltered(Ar3_==0)+Teunfiltered(Ar3_==0)))/sum(Tiunfiltered(Ar3_==0)+Teunfiltered(Ar3_==0));
    VO2expdatamean1_sleep=sum(VO2expdata(Ar3_==0).*(Ti(Ar3_==0)+Te(Ar3_==0)))/sum(Ti(Ar3_==0)+Te(Ar3_==0));
    VCO2expdatamean1_sleep=-sum(VCO2expdata(Ar3_==0).*(Ti(Ar3_==0)+Te(Ar3_==0)))/sum(Ti(Ar3_==0)+Te(Ar3_==0));
    
    
    
    if exist('VolTh1','var')
        meanVdotRIP = 60*sum(VdotRIP.*(Tiunfiltered+Teunfiltered))/sum(Tiunfiltered+Teunfiltered);
    else
        meanVdotRIP = NaN;
    end
    
    %repeat but make this time-weighted average
    VO2_2mean=sum(VO2_2.*(Ti+Te))/sum(Ti+Te);
    %VO2_2mean2=mean(VO2_2);
    
    %use actual PIO2 measured (even if not true PIO2)
    VO2_3 = 273/310*1000*(VTa./(Ti+Te).*60.*(PIO2_Bmax-PO2_Bmin)'/760);
    VO2_3mean=sum(VO2_3.*(Ti+Te))/sum(Ti+Te);
    close all
    figure(2); 
        set(gcf,'color',[1 1 1]); 
        X=6;
    ax2(1)=subplot(X,1,1); plot(Time_,EpochsXHz_,Time_,EventsArXHz_,Time_,EventsRespXHz_); nicefig(); ylabel('Hyp/Events'); legend('Hyp','Ar','Resp'); legend('boxoff');
    %ax2(2)=subplot(X,1,2); plot(Time_,Flow_,Time_,FlowC_,Time_(I.starti),FlowC_(I.starti),'r.',Time_(I.midi),FlowC_(I.midi),'k.'); nicefig(); ylabel('Flow'); legend('original','filtered'); legend('boxoff'); 
    ax2(2)=subplot(X,1,2); plot(Time_,Flow_,'k'); nicefig(); ylabel('Flow'); legend('boxoff'); 
    %ax2(3)=subplot(X,1,3); plot(Time_,PO2_,Time_,PO2f_,Time_(PO2_Bi),PO2f_(PO2_Bi),'r.',Time_(I.midi),PO2f_(I.midi),'k.',Time_(PO2_Bmini),PO2f_(PO2_Bmini),'g.',Time_(PIO2_Bmaxi),PO2f_(PIO2_Bmaxi),'b.'); nicefig; ylabel('PO2'); legend('original','filtered','end exp (ET)','end insp','min'); legend('boxoff'); nicefiglast();
    ax2(3)=subplot(X,1,3); plot(Time_,PO2f_,'k',Time_(PO2_Bi),PO2f_(PO2_Bi),'r.',Time_(I.midi),PO2f_(I.midi),'k.',Time_(PO2_Bmini),PO2f_(PO2_Bmini),'g.',Time_(PIO2_Bmaxi),PO2f_(PIO2_Bmaxi),'b.'); nicefig(); ylabel('PO2'); %legend('signal','end exp','min','end insp','max'); legend('boxoff'); 
    ylim([min([79.5 min(PO2f_)-1]) max([151 max(PO2f_)+1])]);
    %ax2(4)=subplot(X,1,4); plot(Time_,PCO2_,Time_,PCO2f_,Time_(I.endi),PCO2f_(I.endi),'r.',Time_(I.midi),PCO2f_(I.midi),'k.'); nicefig(); ylabel('PCO2'); legend('original','filtered','end exp (ET)','end insp','max'); legend('boxoff'); 
    hold('on');
    try
        plot(Time_(PO2_Bmini(criteriasummary==1)),PO2_Bmin_modelcombo(criteriasummary==1),'ro','markersize',2); 
        hold('on');
        plot(time_dt_win,PO2est_win,'r');
        hold('off');
    catch me
        disp(me.message);
    end
    
    ax2(4)=subplot(X,1,4); plot(Time_,round(SaO2_),'r'); nicefig();
    ylim([min([89.5 min(SaO2_)-1]) 100.5]);
    ylabel('SpO2')

    ax2(5)=subplot(X,1,5); plot(Time_,Pes_-Pes_baselineest,'k');  %legend('original','filtered','end exp (ET)','end insp','max'); legend('boxoff'); 
    hold('on');
    try
        for i=1:length(PesAtAr)
        itemp = ISleepOff(i)-2+previousorcurrent(i);
        plot([Time_(Ipes.starti(itemp)) Time_(Ipes.endi(itemp))],[Pesf_(Ipes.starti(itemp)) Pesf_(Ipes.starti(itemp))]-Pes_baselineest,'g');
        plot([Time_(Ipes.starti(itemp)) Time_(Ipes.endi(itemp))],-PesAtAr(i)+[Pesf_(Ipes.starti(itemp)) Pesf_(Ipes.starti(itemp))]-Pes_baselineest,'g');
        plot([Time_(Ipes.starti(itemp)) Time_(Ipes.starti(itemp))],[Pesf_(Ipes.starti(itemp))-PesAtAr(i) Pesf_(Ipes.starti(itemp))]-Pes_baselineest,'g');
        plot([Time_(Ipes.endi(itemp)) Time_(Ipes.endi(itemp))],[Pesf_(Ipes.starti(itemp))-PesAtAr(i) Pesf_(Ipes.starti(itemp))]-Pes_baselineest,'g');
        end
    catch me
    end
    hold('off');
    nicefig(); ylabel('Pes');
    
    ylim([min([-31 min(Pesf_-Pes_baselineest)-1]) max([5 max(Pesf_-Pes_baselineest)+1])]);
    ax2(6)=subplot(X,1,6); 
    %stairs(Time_(I.starti),VO2mean1B,'k'); nicefiglast(); ylabel('VO2'); %legend('original','filtered','end exp (ET)','end insp','max'); legend('boxoff'); 
                           %hold('on'); stairs(Time_(I.starti),VO2mean1+0*VO2mean1B,'k:'); 
                           %stairs(Time_(I.starti),VO2_3,'r'); nicefiglast(); ylabel('VO2'); %legend('original','filtered','end exp (ET)','end insp','max'); legend('boxoff'); 
                           
                           %ylim([0 max(VO2_3)+5]);
    if 1
        stairs(Time_(I.starti),VO2expdata,'b'); nicefiglast(); ylabel('VO2'); %legend('original','filtered','end exp (ET)','end insp','max'); legend('boxoff'); 
        hold('on'); stairs(Time_(I.starti),VO2expdatamean1+0*VO2expdata,'b:'); 
        stairs(Time_(I.starti),referenceVO2+0*VO2_3,'k');            
    end
    if 0
        stairs(Time_(I.starti),VO2mean1B,'b:'); ylabel('VO2'); %legend('original','filtered','end exp (ET)','end insp','max'); legend('boxoff'); 
    end
    ylim([0 max([referenceVO2+20 max(VO2expdata)+1])]);
    
    nicefiglast();
    linkaxes(ax2,'x');
    %set(gcf,'Position',[1100 60 560 900]);
    
    if 1
    for i=1:3%length(I.starti)
        figure(90);        
        subplot(2,1,1); plot((Time_((I.midi(i)):I.endi(i)))-Time_(I.midi(i)),PO2f_((I.midi(i)):I.endi(i))); hold('on');
        subplot(2,1,2); plot(cumsum(-dt*FlowC_((I.midi(i)):I.endi(i))),PO2f_((I.midi(i)):I.endi(i))); hold('on');
        figure(91);        
        subplot(2,1,1); plot((Time_((I.starti(i)):I.midi(i)))-Time_(I.starti(i)),PO2f_((I.starti(i)):I.midi(i))); hold('on');
        subplot(2,1,2); plot(cumsum(dt*FlowC_((I.starti(i)):I.midi(i))),PO2f_((I.starti(i)):I.midi(i))); hold('on');        
    end
    end
    
    %% Summary Info
    try
        II = ~isnan(PO2_Bmin_model);
        PO2_Bmin_model_mean = sum(PO2_Bmin_model(II).*Ttot(II)')/sum(Ttot(II));
        II = ~isnan(PO2_Bmin_modelcombo);
        PO2_Bmin_modelcombo_mean = sum(PO2_Bmin_modelcombo(II).*Ttot(II)')/sum(Ttot(II));
        II = ~isnan(PCO2_Bmax_model);
        PCO2_Bmax_model_mean = sum(PCO2_Bmax_model(II).*Ttot(II)')/sum(Ttot(II));
        II = ~isnan(PCO2_Bmax_modelcombo);
        PCO2_Bmax_modelcombo_mean = sum(PCO2_Bmax_modelcombo(II).*Ttot(II)')/sum(Ttot(II));
    catch me
        disp(me.message);
    end
    
    PO2_Bmin_mean = sum(PO2_Bmin.*Ttot')/sum(Ttot);    
    PCO2_Bmax_mean = sum(PCO2_Bmax.*Ttot')/sum(Ttot);
    
    
    IX=find(diff(EventsRespXHz_)>0);%&(EventsRespXHz_(1:end-1)==0))
    eventslist = EventsRespXHz_(IX+1);
    eventstypes = {'NObsA','NCentralA','NObsH','NMixedA','NCentralH'};
    eventscodes = [2 3 4 5 6];
    for i=1:length(eventscodes)
        eval([eventstypes{i} '=sum(eventslist==eventscodes(i))']);
    end
    Fobstructed1 = sum(EventsRespXHz_==2|EventsRespXHz_==4)/length(EventsRespXHz_);
    Fobstructed2 = sum(EventsRespXHz_==2|EventsRespXHz_==4|EventsRespXHz_==5)/length(EventsRespXHz_);
    Fobstructed = mean([Fobstructed1 Fobstructed2]);
    
    Farousal = sum(EventsArXHz_==1)/length(EventsRespXHz_);
    
    IX=find(diff(EventsArXHz_)>0);%&(EventsRespXHz_(1:end-1)==0))
    arousalslist = EventsArXHz_(IX+1);
    NAr = sum(arousalslist);
    
    temp = VTi2./VTe2-1;
    VTionVTeERR = mean(abs(temp(VTe2>0.05))); %trying to exclude values that given inf when VTe = 0
    
    %SaO2_=round(SaO2(Li:Ri));
    meanSpO2 = mean(SaO2_);
    
    VO2summary(m,:) = [VO2mean1 VO2_3mean Vflow_to_PO2_delay_measured VTionVTeERR ...
        meanSpO2 meanVdot meanVa CPAPlevel_est NObsA NCentralA NObsH NMixedA NCentralH ...
        NAr Fobstructed Farousal PO2_Bmin_mean meanVdotRIP mode(Position_) PCO2_Bmax_mean ...
        VO2expdatamean1 mixedexpiredmean mixedexpiredmean3 VCO2expdatamean1 Pmus_mean2 PTP12 PTP22 ...
        PO2_Bmin_model_mean PO2_Bmin_modelcombo_mean PCO2_Bmax_model_mean PCO2_Bmax_modelcombo_mean ...
        Pes_mean2 PesAtArmedian PmusAtArmedian PesAtArmedianN ...
        meanVdot_sleep VO2expdatamean1_sleep VCO2expdatamean1_sleep Pmus_mean2_sleep Pes_mean2_sleep ... 
        ];


    %% save figure as
    tempstr = [];
    for i=1:3-length(num2str(m))
        tempstr = ['0' tempstr];
    end
    saveas(figure(2),[directory 'figures\' filename '_' tempstr num2str(m)]);

    %% Pause
    if pause_to_assess_figures_2
        pause
    end
    %% save workspace
    %save([directory ' filename '_' tempstr num2str(m)]);
   

    %% Continuous PAO2 trace
    %[IX,VI,VT,Ti,Te]=Vflowanalysis2(Flow,Time,dt,0);
if 0    
    %SaO2_=round(SaO2(Li:Ri));
    
    VI_rs=Time_*NaN;
    VI_rs(1:I.starti(1)-1)=Vdot(1);
    for i=I.starti(1):length(Time_);
        VI_rs(i) = Vdot(find(Time_(I.starti)<=Time_(i),1,'last')); 
    end

[XCR,lags]=crosscorr(VI_rs,SaO2,round(30/dt));
XCR(lags<0)=[];
lags(lags<0)=[];
%figure(101); plot(lags,XCR);
[maxXCR,i]=max(XCR); delay_xcr=lags(i)*dt;
delaymin=round(0.33*delay_xcr);
delaymax=round(0.8*delay_xcr);
if delaymax>=30,
    delaymax=29;
end

%%
figure(1)
ax(1)=subplot(4,1,1); plot(Time_,FlowC_);
ax(2)=subplot(4,1,2); plot(Time_,VI_rs);
ax(3)=subplot(4,1,3); plot(Time_,PO2f_);
iindex = [I.starti;I.endi(end)];
    FVeupnea_est=0;
    Peupnea_est = mean(PO2f_(iindex));
    error1_est = PO2f_(iindex(1))-Peupnea_est;
    k1_est = 50; k2_est = 0.002;
Parameters=[k1_est k2_est error1_est 0.2]; %k tau error1 Fdeadspace Peupnea
lower =[k1_est/10 k2_est/10 -20 0.1];
upper =[k1_est*10 k2_est*10 20 0.4];
ParametersF
%FVeupnea=Parameters(4);%-0.05; Veupnea=mean(VI_rs)*(1-FVeupnea);
%Fdeadspace = Parameters(5)*mean(VI_rs);
%error1=Parameters(3);%0;
%k2=Parameters(2);%40;
%k1=Parameters(1);%15;
iindex = [I.starti;I.endi(end)];
[F,PAO2]=VEtoPAO2f(Parameters,VI_rs,PO2_,dt,iindex,FVeupnea_est,Peupnea_est);

figure(1)
ax(3)=subplot(4,1,3); plot(Time_,PO2f_,Time_,PAO2);

%%
OPTIONS = optimset('Display','off','TolX',1e-1,'TolFun',1e-1,'MaxIter',200,'MaxFunEvals',200,'Algorithm','interior-point');

    clear Parameters1 Parameters2 Fres1 Fres2 delayi
    %FVeupnea=[(FVeupnea_est-0.08):0.01:(FVeupnea_est+0.08)];
    %delay=delaymin:1:delaymax;
    dsf=1;
    j=1;
    for i=1:length(FVeupnea)
        for j=1:length(delay)
        [Parameters1(j,:),Fres1(j),~,~] = fmincon(@(Parameters) VEtoPAO2f(Parameters,downsample(VI_rs,dsf),downsample(PO2_,dsf),dt*dsf,round(downsample(iindex,dsf)),FVeupnea(i)),Parameters,[],[],[],[],lower,upper,[],OPTIONS); %@(Parameters) model_nonlinear_constraints() is now [],[]
        end
        [minFres,j]=min(Fres1);
        Parameters2(i,:)=Parameters1(j,:);
        Fres2(i)=Fres1(j);
    end
    [minFres,i]=min(Fres2);
    figure(99); plot(FVeupnea,Fres2);
    ParametersF=Parameters2(i,:);
    FVeupneaF = FVeupnea(i);
    [F,PAO2]=VEtoPAO2f(ParametersF,VI_rs,PO2_,dt,iindex,FVeupneaF);
    figure(1); ax(3)=subplot(4,1,3); plot(Time_,PO2f_,Time_,PAO2);
    
end
    
    
      m=m+1;  
end
    