global AMasterSpreadsheet settings ChannelsList ChannelsFs
addpath('C:\Users\bwhha\Dropbox (Partners HealthCare)\MEEI DISE\Matlab Scripts')
settings1 = ImportSettingsAnalysis(settings,AMasterSpreadsheet); 
settings=settings1;
[num,patients,~] = xlsread(AMasterSpreadsheet,1,'AD4:AH10003');
NaNlist = (isnan(num)); num(NaNlist) = 0; % set missing or NaN as 0
analyzelist = logical(num(:,2));
settings.invertflowlist = logical(num(:,1));
settings.protocol = patients(:,3);
TotalNumPts=size(patients,1);
settings.patients = patients;

PtRangeTemp = 1:1:TotalNumPts; %normally
PtRange = PtRangeTemp(analyzelist==1);
ptnum = PtRange(1);

%% Load converted data
settings.filename=[char(settings.patients(ptnum,1))]; %seems to be unused. 
% DLM says, settings.filename is used in LGfromFlowBetaPart1 to
% save flowdrive plot, and Part2 to save loop gain plot
directoryn=char(settings.patients(ptnum,2));
MATfilename=[directoryn char(settings.patients(ptnum,1))];
Evts=struct();
temp=[];
temp=load(MATfilename);
DataEventHypnog_Mat=temp.DataEventHypnog_Mat;
Evts=temp.Evts;
%WakeSleepInfo=temp.WakeSleepInfo; % removed after modifying 'Info' in Convert 
ChannelsFs=temp.ChannelsFs;
ChannelsList=temp.ChannelsList;
    
% load snore channels list
SnoreMATfilename = [MATfilename(1:end-4), '_snore'];
temp2=load(SnoreMATfilename);
SnoreChannelsList=temp2.SnoreChannelsList;

%% Sample data
Time = DataEventHypnog_Mat(:,strcmp(ChannelsList,'Time'));
Flow = DataEventHypnog_Mat(:,strcmp(ChannelsList,'Flow'));
Snore = DataEventHypnog_Mat(:,strcmp(ChannelsList,'Snore'));
SnoreDB = DataEventHypnog_Mat(:,strcmp(ChannelsList,'SnoreDB'));
fsFlow = ChannelsFs(strcmp(ChannelsList,'Flow'));
vol=cumsum(Flow(:))*(1/fsFlow);

% need to include a load option
ivmin = [];
ivmax = [];
%% Run GUI
BBStartStopAgain()

