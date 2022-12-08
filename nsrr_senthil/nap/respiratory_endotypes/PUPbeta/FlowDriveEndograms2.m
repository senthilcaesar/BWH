
global SAfontsize settings
subj=49; %59 9
settings.AnalyzedDirectory='E:\Dropbox (Partners HealthCare)\PhenotypeDrive2018\AnalyzedDDOSA\';

% set some defaults in case of encountering an error
Data = NaN; DataCI = NaN; DataN = NaN; AHItotal = NaN; Fstates = NaN;
Veupnea = NaN; medianV = NaN; VEdeciles = NaN; Vdrivedeciles = NaN;
Fsupinetemp = NaN;

AHIdata = NaN;
settings.selectstate=4
 Nbootstrap=50;
 SAfontsize=14
 
%% Hardcoded directory info (and additional settings for Pes)
if exist('settings.directory')
    directory=settings.directory;%'H:\MESA\Analyzed3';
else
    directory=settings.AnalyzedDirectory;%'H:\MESA\Analyzed3';
end

savename =settings.savename;%'MESA2Sept2017';

if ~isfield(settings,'DriveSignal')
    settings.DriveSignal='DeltaEdi';
    disp(['Warning default DriveSignal used: DeltaEdi']);
end
disp(settings.DriveSignal)

scoredarousalsinwake=settings.scoredarousalsinwake;
DriveincmH2OnotEupnea=settings.DriveincmH2OnotEupnea;


%% Default for figures
set(groot,'defaultAxesfontname','arial narrow');
set(groot,'defaultFigureColor',[1 1 1]);
set(groot,'defaultAxesBox','off');
%set(groot,'defaultAxesfontsize',SAfontsize);

%% Load PUPdata

%subj=1;
getCIs=settings.getCIs;
varlist = {'LG1','LGn','delay','VRA','ArThres','Vpassive','Vactive','Vcomp','ArThresActive'};
Nperwindowvars=5;

loadpath=[settings.AnalyzedDirectory settings.savename '_' num2str(subj)];

listOfVariables = who('-file', [loadpath '.mat']);
datatoload = load(loadpath);
if isfield(datatoload,'settings')
    datatoload = rmfield(datatoload,'settings');
end

if sum(strcmp(fieldnames(datatoload),'Evts')) == 0 % 'startsWith' or 'endsWith' or 'contains' are better at string matchy matchy
    if isfield(datatoload,'EvtsData')
        datatoload.Evts = datatoload.EvtsData;
        datatoload = rmfield(datatoload, 'EvtsData');
    elseif isfield(datatoload,'Evtsdata')
        datatoload.Evts = datatoload.Evtsdata;
        datatoload = rmfield(datatoload, 'Evtsdata');
    end
end

%Replaces 'SleepData','LGplusinfo','LG_QualityInfo','DataOut','fitQual',etc with WinT
if ~isfield(datatoload,'WinT')
    WinT = getWinTfromFileW(datatoload);
end

cellfun(@(x,y) assignin('caller',x,y),fieldnames(datatoload),struct2cell(datatoload)) %Shifts all fields from from datatosave to current workspace:

try
    AHIdata=Evts.AHIdata;
end

if iscell(AHIdata)
    AHIdata=AHIdata{1};
end


if ismember('settings', listOfVariables)
    settings2=load(loadpath,'settings');
end

if ismember('settingsAnalyzed', listOfVariables)
    try
        settingstemp=load(loadpath,'settingsAnalyzed'); % DV: THERE IS A BUG HERE - NEED TO RUN ANALYZED BEFORE RUNNING SUMMARY
        settings2.settings = settingstemp.settingsAnalyzed;
    end
end

%%
[BreathDataTable2,~,~,~,~]=GetNonOvlappedVE(BreathDataTable,[]);

BreathDataTable2.Win = floor(BreathDataTable2.UniqueID/1000);

if ~any(BreathDataTable2.Properties.VariableNames=="AR3") && exist('DataOut')==1
    DataOutT = DataOut2T(DataOut);
    BreathDataTable2.AR3 = DataOutT.a;
    BreathDataTable2.Vdr_est = DataOutT.x-1; %1 will be added below
    
end
BreathDataTable2.notAR3 = calculateNotAR(BreathDataTable2.AR3,BreathDataTable2.Win,2);

%Window Edge/Boundary code:
BreathDataTable2.WinBoundary = sum([[diff(BreathDataTable2.Win)~=0;1],[1;diff(BreathDataTable2.Win)~=0]],2);

switch settings.DriveSignal %DeltaPes, DeltaPmus,
    case 'DeltaEdi'
        
        if settings.EdiAutoLinearization==1
            [BreathDataTable2,DriveExponent,Drive] = VdriveFromEdi(BreathDataTable2,'VE',settings.DriveSignal,1);
            BreathDataTable2.DeltaEdi=Drive; %replacement will provide linearized signal.
            if settings.EdimtaCalibration
                BreathDataTable2.DeltaEdi=BreathDataTable2.DeltaEdimtaCal;
            end
        end
    case 'VdriveModel'
        BreathDataTable2.VdriveModel = BreathDataTable2.Vdr_est + 1;
end

%Calculate Drive gradient
BreathDataTable2.DriveSlope=gradient(BreathDataTable2.(settings.DriveSignal));
BreathDataTable2.DriveSlope(BreathDataTable2.WinBoundary==1)=NaN;

%% Add VdriveEdiNorm to BreathDataTable2
%% Comb
if isfield(settings,'comb') && settings.comb(1)
    disp('comb window selection');
    width = settings2.settings.windowlength;
    step=settings2.settings.WindowStep/60;
    
    settings.comb(2);
    winn = (1:height(WinT));
    starttime = (winn - 1).*step;
    endtime=starttime+width;
    temp = [starttime;endtime];
    swap = 1*(settings.comb(2)==2);
    temp2 = mod(temp + settings.comb(3)*swap,2*settings.comb(3))/(2*settings.comb(3));
    criteriacomb = (sum(temp2<=0.5)==2)';
else
    criteriacomb = ones(height(WinT),1);
end

%% Check WinT is complete
if ~any(WinT.Properties.VariableNames=="Nevents")
    disp('Warning WinT is missing essential data, Analysis Failed (suspected issue merging tables with different heights after an error)');
end

%%
if sum(size(AHIdata)==[1,1])==2 && isnan(AHIdata)
    AHIdata=nan(1,184);
end

try
    AHItotal=AHIdata(58);
catch
    AHItotal=NaN;
end

%%


%% Analyse LG, ArThres

clear LG1 LG2 LGn temp temp2 FVAL VRA1 VRA2 ArThres MeanEx

% events, longestwake, position, FremMAX
minNeventsLG = settings.minNeventsLG; %1
maxwakeLG = settings.maxwakeLG; %30