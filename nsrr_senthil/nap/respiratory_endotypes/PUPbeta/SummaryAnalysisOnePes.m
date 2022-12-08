function [Data,DataCI,DataN,varlist,AHItotal,Fstates,...
    Veupnea,medianV,VEdeciles,Vdrivedeciles,AHIdata,Fsupine,x_,VEdecilesUpper,VEdecilesLower,GGdata,GGpdeciles,GGtdeciles,EdecilesMean,TdecilesMean] = SummaryAnalysisOnePes(subj,settings)

global SAfontsize

% set some defaults in case of encountering an error
Data = NaN; DataCI = NaN; DataN = NaN; AHItotal = NaN; Fstates = NaN;
Veupnea = NaN; medianV = NaN; VEdeciles = NaN; Vdrivedeciles = NaN;
Fsupinetemp = NaN;

AHIdata = NaN;

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
set(groot,'defaultAxesfontsize',SAfontsize);

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

% try
%     Fnoise = StoNData.Fnoise;
%     FnoiseT = array2table(Fnoise); FnoiseT.Properties.VariableNames = {'Fnoise1','Fnoise2','Fnoise3'};
%     FnoiseT.Fnoise2Under10percent = FnoiseT.Fnoise2<0.1; % nnz(FnoiseT.Fnoise2Under10percent)
% catch me
%     disp('no StoNData available');
% end

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
    
    % if ~any(BreathDataTable2.Properties.VariableNames=="AR3")
    %     BreathDataTable2.AR3 = BreathDataTable2.ARie;
    %     BreathDataTable2.AR3(BreathDataTable2.Etype>0)=0;
    % end
end
BreathDataTable2.notAR3 = calculateNotAR(BreathDataTable2.AR3,BreathDataTable2.Win,2);
% %
% % BreathDataTable2.notAR = calculateNotAR(BreathDataTable2.AR,BreathDataTable2.Time0,2);
% % BreathDataTable2.notAR3 = calculateNotAR(BreathDataTable2.AR3,BreathDataTable2.Time0,2);

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
state=settings.selectstate; %1=nrem1, 4 nrem, 5 rem, 8 ignore state

%% Analyse LG, ArThres

clear LG1 LG2 LGn temp temp2 FVAL VRA1 VRA2 ArThres MeanEx

% events, longestwake, position, FremMAX
minNeventsLG = settings.minNeventsLG; %1
maxwakeLG = settings.maxwakeLG; %30
%% State, 1=nrem1, 2=nrem2, 3=nrem3, 4 = all nrem, 5 rem, 8 ignore state
%defaults
minFnrem1=-Inf;
minFnrem2=-Inf;
minFnrem3=-Inf;
switch settings.selectstate
    case 1
        maxFREM=0;
        minFrem=-Inf;
        minFnrem1=0.5;%settings.statetolerance;
        AHIstate=AHIdata(88); %not correct
        DurationInState=AHIdata(88-7);
        PlotString = 'N1';
        hypok = [2];
    case 2
        maxFREM=0;
        minFrem=-Inf;
        minFnrem2=0.5; %settings.statetolerance;
        AHIstate=AHIdata(88); %not correct
        DurationInState=AHIdata(88-7);
        PlotString = 'N2';
        hypok = [1];
    case 3
        maxFREM=0;
        minFrem=-Inf;
        minFnrem3=0.5; %settings.statetolerance;
        AHIstate=AHIdata(88); %not correct
        DurationInState=AHIdata(88-7);
        PlotString = 'N3';
        hypok = [0];
    case 4 % at the moment this restricts to any window which does not have any REM.
        maxFREM=0; %1-settings.statetolerance;
        minFrem=-Inf;
        AHIstate=AHIdata(88);
        DurationInState=AHIdata(88-7);
        PlotString = 'NREM';
        hypok = [0 1 2];
    case 5 % at the moment this restricts to a window with >50% of REM, and the rest can be anything.
        maxFREM=Inf;
        minFrem=0.5; %settings.statetolerance;
        AHIstate=AHIdata(96);
        DurationInState=AHIdata(96-7);
        PlotString = 'REM';
        hypok = [3];
    case 8 % at the moment this does not restrict at all, and could include periods which are predominantly wake unless a wake threshold is applied elsewhere.
        maxFREM=Inf;
        minFrem=-Inf;
        AHIstate=AHIdata(80);
        DurationInState=AHIdata(80-7);
        PlotString = 'ALL';
        hypok = [0 1 2 3];
end

%% Find first and last windows containing sleep
TimeOfNightCriteria = ones(height(WinT),1);
if settings.selecttimeofnight
    XTiles=settings.selecttimeofnight_XTiles;
    NthXTile=settings.selecttimeofnight_NthXTile;
    
    FwakePerWindow=WinT.FWake;
    SleepWin1=find(FwakePerWindow<1,1,'first');
    SleepWinN=find(FwakePerWindow<1,1,'last');
    
    rangeWin=round((SleepWinN-SleepWin1+1)/XTiles);
    lowerWinNs=round(SleepWin1+((1:XTiles)-1)*rangeWin);
    upperWinNs=round(lowerWinNs+rangeWin-1);
    
    TimeOfNightCriteria=0*TimeOfNightCriteria;
    TimeOfNightCriteria(lowerWinNs(NthXTile):upperWinNs(NthXTile))=1;
end

%%
usemediannotmeanLG1=1;
% %
% % i=1;
% % clear Tn LG1
% % tempLG1=LGplusinfo(:,7);
% % FisNaN(i)=sum(isnan(tempLG1))/length(tempLG1);
% %
% % temp=LGplusinfo;
% %
% % tempLG1=LGplusinfo(:,7); tempLGn=LGplusinfo(:,5);  tempLG2=LGplusinfo(:,8); tempTn=LGplusinfo(:,6);
% % tempLG0=LGplusinfo(:,2); temptau=LGplusinfo(:,3);  tempdelay=LGplusinfo(:,9);
% % tempVRA1=LGplusinfo(:,10); tempVRA2=LGplusinfo(:,11); tempArThres=LGplusinfo(:,12);
% %
% % temppositiondata=LG_QualityInfo(:,5);
% % minwake=SleepData(:,7);
% % FREM=SleepData(:,6);

% Pos = LG_QualityInfo(:,5);
% if isempty(settings.selectposition)
%     Poscriteria = ones(length(Pos),1);
% else
%     Poscriteria = sum(Pos == settings.selectposition,2)>0;
% end

%% Position

% converting a single cell to table, but if cell array of tables,
% then handle differently, and don't {} convert  (DLM)
if isa(BreathDataTable,'cell') && all(size(BreathDataTable)==([1,1]))
    BreathDataTable=BreathDataTable{1};
end

% old position coding
PosWin = zeros(size(BreathDataTable,2),1);
Fsupinetemp = zeros(size(BreathDataTable,2),1);
% new position coding
try
    positioncodes = PositionCodeLookup(settings.poscodesdatabase,settings.protocol{subj});
catch % RMA added--for running single subject A to Z mode
    positioncodes = PositionCodeLookup(settings.poscodesdatabase,settings.protocol(subj,:));
end

for xx=1:size(BreathDataTable,2)
    if size(BreathDataTable{xx},1)<=1
        PosWin(xx)=NaN;
        Fsupinetemp(xx)=NaN;
        continue
    end
    temp = BreathDataTable{xx}.pos_B;
    temp2 = PositionTranslator(positioncodes,settings.positioncodesout,temp);
    PosWin(xx) = mode(temp2);
    Fsupinetemp(xx) = nanmean(PositionSelector(temp2,'Supine'));
end
Poscriteria = PositionSelector(PosWin,settings.selectposition);
str=[num2str(nnz(Poscriteria)), ' of ', num2str(length(Poscriteria)) ,' windows in selected position'];
disp(str);


%% fractions of each sleep state
%SleepData(winNum+1,:)=[FWAKE FNREM FNREM1 FNREM2 FNREM3 FREM LongestWake];
if 1
    Funknownsleep=round(1-sum([WinT.FWake WinT.FNREM WinT.FREM],2),2);
else
    Funknownsleep=(1-sum([WinT.FWake WinT.FNREM WinT.FREM],2)); %old code, removed 5/7/2021. delete option in future
end
containsunknownsleep=1*(Funknownsleep>0);


%% This criteria defines a minimum proportion of particular non-REM sleep states.
nremXcriteria= WinT.FNREM1>=minFnrem1 &...
    WinT.FNREM2>=minFnrem2 & ...
    WinT.FNREM3>=minFnrem3;


%%

%% criteria for LG and UA phenotype
minNeventsLG = settings.minNeventsLG; %default 1
minNeventsUA = settings.minNeventsUA; %has been default 1, was default 0 for a long time previously
maxwakethresLG = settings.maxwakethresLG; % 30
maxwakethresUA = settings.maxwakethresUA; % 300; more tolerant to more wake for UA measures

try %added SS 11/2/2021
    CPAPoffCriteria = WinT.CPAPoff==1;
catch
    CPAPoffCriteria = zeros(1,height(WinT))+1;
    disp('warning: CPAPoff detection failed, assuming CPAP is off always');
end


criteriaAll = criteriacomb & ...
    Poscriteria==1 & ...
    CPAPoffCriteria==1 & ...
    nremXcriteria==1 & ...
    WinT.FREM<=maxFREM &...
    WinT.FREM>=minFrem & ...
    containsunknownsleep==0 & ...  %recently added contains unknownsleep
    TimeOfNightCriteria==1;
if settings.verbose
    disp(['criteriaAll = ', num2str(nnz(criteriaAll)), ' / ', num2str(length(criteriaAll))]);
end

criteriaUA = criteriaAll & ...
    WinT.Nevents>=minNeventsUA & ... % nnz(isfinite(N_events))
    WinT.LongestWake<=maxwakethresUA;
if settings.verbose
    disp(['criteriaUA = ', num2str(nnz(criteriaUA)), ' / ', num2str(length(criteriaUA))]);
end

% criteria for LG and Arousal Threshold
criteriaLG = criteriaAll & ...
    WinT.Nevents>=minNeventsLG & ...
    WinT.LongestWake<=maxwakethresLG;
if settings.verbose
    disp(['criteriaLG = ', num2str(nnz(criteriaLG)), ' / ', num2str(length(criteriaLG))]);
end

criteriaEither = criteriaUA|criteriaLG; %currently same as criteriaUA i.e. least strict criteria
if settings.verbose
    disp(['criteriaEither = ', num2str(nnz(criteriaEither)), ' / ', num2str(length(criteriaEither))]);
end


%%

% 6 is code for upright, which is a bit weird during sleep, so show a message
PosUpright = sum(PosWin == 6,2)>0; %note: PosWin actually only one column so why sum2?
if sum(PosUpright & criteriaLG)>0; disp(['Caution: ' num2str(sum(PosUpright & criteriaLG)) ' windows of Upright sleep found in this study']); end

NwindowsForLG=sum(criteriaLG);
str=['Using ', num2str(NwindowsForLG), ' of ', num2str(length(criteriaLG)) ,' windows for LG calcs'];
disp(str);

%% Bring in LG/VRA/ArThres window data

clear Tn LG1

%only works for 1 time constant:
alpha=(WinT.LG2./WinT.LG1).^2;
beta=(1-alpha)./(4*alpha-1);
WinT.LG0p17=WinT.LG1.*((1+beta)./(1+beta*(1/6)^2)).^0.5;
WinT.LG0p33=WinT.LG1.*((1+beta)./(1+beta*(1/3)^2)).^0.5;
WinT.LG0p67=WinT.LG1.*((1+beta)./(1+beta*(2/3)^2)).^0.5;
WinT.LG1p33=WinT.LG1.*((1+beta)./(1+beta*(4/3)^2)).^0.5;
WinT.LG1p67=WinT.LG1.*((1+beta)./(1+beta*(5/3)^2)).^0.5;
WinT.LG3=WinT.LG1.*((1+beta)./(1+beta*(3)^2)).^0.5;
WinT.LG4=WinT.LG1.*((1+beta)./(1+beta*(4)^2)).^0.5;

WinT.VRA(WinT.Narousals<1)=NaN;
WinT.VRA2(WinT.Narousals<1)=NaN;

WinT.criteriaLG = criteriaLG;
WinT.Poscriteria = Poscriteria;
WinT.nremXcriteria = nremXcriteria;
WinT.TimeOfNightCriteria = TimeOfNightCriteria;
WinT.criteriacomb = criteriacomb;
WinT.Fsupine = Fsupinetemp;
% try
%     Twin = [Twin FnoiseT];
% catch ErrorConcatenatingWindowTables
%     disp(ErrorConcatenatingWindowTables.message);
% end
%% Calculate summary per criteria
criteriaLGI = find(criteriaLG==1);
[LG1,LGn,delay,VRA,ArThres,LG1N,LGnN,delayN,VRAN,ArThresN,...
    Fstates1,Fsupine1,Fstates2,Fsupine2,...
    LG1values,LGnvalues,delayvalues,VRAvalues,ArThresvalues]=...
    SummaryAnalysisOne_LGmodel(WinT(criteriaLGI,:),varlist);

%%


% nnz(FnoiseT.Fnoise2Under10percent)
% %
% % Fnrem1=SleepData(:,3);
% % Fnrem2=SleepData(:,4);
% % Fnrem3=SleepData(:,5);
% % nremXcriteria=Fnrem1>=minFnrem1&Fnrem2>=minFnrem2&Fnrem3>=minFnrem3;
% %
% % minwake(length(tempLG1)+1:end)=[];
% % FREM(length(tempLG1)+1:end)=[];
% % nremXcriteria(length(tempLG1)+1:end)=[];
% % TimeOfNightCriteria(length(tempLG1)+1:end)=[];


% % %SleepData(winNum+1,:)=[FWAKE FNREM FNREM1 FNREM2 FNREM3 FREM LongestWake];
% %
% % alpha=(tempLG2./tempLG1).^2;
% % beta=(1-alpha)./(4*alpha-1);
% % tempLG3min=tempLG1.*((1+beta)./(1+beta*(1/3)^2)).^0.5; %error, needs fixing...
% % tempLG6min=tempLG1.*((1+beta)./(1+beta*(1/6)^2)).^0.5; %error, needs fixing...
% % tempLG90s=tempLG1.*((1+beta)./(1+beta*(2/3)^2)).^0.5; %error, needs fixing...

%GetRsquared:
%temp3=fitQual;
% %
% % clear OneMinusRsq
% % for ii=1:size(fitQual,2)
% %     temp3pt1=fitQual{ii};
% %     if length(temp3pt1)>1
% %         OneMinusRsq(ii,:)=1-temp3pt1(2); %Element #2 of FitQual
% %     else
% %         OneMinusRsq(ii,:)=NaN;
% %     end
% % end
% % % Rsq=temp3(:,1);
% % FVAL=LGplusinfo(:,13);
% %
% % MeanE=LG_QualityInfo(:,3); %is actual scored events. LG_QualityInfo=[N_arousals_used N_events mean(E) mean(E1) position_mode N_position_changes Percent_position EXITFLAG];
% % SStot=FVAL./OneMinusRsq;
% %
% % N_events=LG_QualityInfo(:,2); tempN_arousals=LG_QualityInfo(:,1);
% %
% % criteria=CPAPData(:,1)==1&Poscriteria==1&nremXcriteria==1&TimeOfNightCriteria==1&minwake<=maxwakeLG&N_events>=minNeventsLG&FREM<=maxFREM&FREM>=minFrem; %(temppositiondata==0|temppositiondata==x(4))&
% %
% % %CPAPData(winNum,:)=[CPAPoff_ CPAPmedian CPAPstd CPAPabs95];
% %
% % % New:
% % tempVRA1(tempN_arousals<1)=NaN;
% % tempVRA2(tempN_arousals<1)=NaN;
% %
% %
% % LG1(i,1)=prctile(tempLG1(~isnan(tempLG1)&criteria),50);
% % LG1values=tempLG1(~isnan(tempLG1)&criteria);
% %
% % LG1_SD(i,1)=std(tempLG1(~isnan(tempLG1)&criteria));
% % LG1_N(i,1)=sum(~isnan(tempLG1)&criteria); % no_of_LGmeasures(i,1)=sum(~isnan(tempLG1)&criteria);
% % LG1_SE(i,1)=LG1_SD(i,1)/LG1_N(i,1)^0.5;
% % LG1_N_nocriteria(i,1)=sum(~isnan(tempLG1)); % no_of_LGmeasures(i,1)=sum(~isnan(tempLG1)&criteria);
% %
% % meanSE=mean(LG1_SE);
% % F_incl=LG1_N./LG1_N_nocriteria;
% % meanFincl=mean(F_incl);
% % minLGN=min(LG1_N);
% %
% % %minmeanmedianNmeasures=[min(no_of_LGmeasures) mean(no_of_LGmeasures) median(no_of_LGmeasures)]
% % % LGn_=tempLGn(criteria);
% % % VRA1_=tempVRA1(criteria);
% % % VRA2_=tempVRA2(criteria);
% % MeanEx(i,1)=prctile(MeanE(~isnan(MeanE)&criteria),50);
% % % LG1(i,1)=prctile(tempLG1(criteria),50)
% % LG2(i,1)=prctile(tempLG2(~isnan(tempLG2)&criteria),50);
% % LGn(i,1)=prctile(tempLGn(~isnan(tempLGn)&criteria),50);
% % LGnvalues=tempLGn(~isnan(tempLGn)&criteria);
% % Tn(i,1)=prctile(tempTn(~isnan(tempTn)&criteria),50);
% % VRA1(i,1)=100*prctile(tempVRA1(~isnan(tempVRA1)&criteria),50);
% % VRA1values=100*tempVRA1(~isnan(tempVRA1)&criteria);
% % VRA2(i,1)=100*prctile(tempVRA2(~isnan(tempVRA1)&criteria),50);
% % ArThres(i,1)=100*prctile(tempArThres(~isnan(tempArThres)&criteria),50);
% % ArThresvalues=100*tempArThres(~isnan(tempArThres)&criteria);
% % % Rsq_median(i,1)=prctile(Rsq(criteria),50)
% % LG0direct(i,1)=prctile(tempLG0(~isnan(tempLG0)&criteria),50);
% % tau(i,1)=prctile(temptau(~isnan(temptau)&criteria),50);
% % delay(i,1)=prctile(tempdelay(~isnan(tempdelay)&criteria),50);
% % delayvalues=tempdelay(~isnan(tempdelay)&criteria);
% % LG3min(i,1)=prctile(tempLG3min(~isnan(tempLG3min)&criteria),50);
% % LG6min(i,1)=prctile(tempLG6min(~isnan(tempLG6min)&criteria),50);
% % LG90s(i,1)=prctile(tempLG90s(~isnan(tempLG90s)&criteria),50);


%limitVRA1=1.59;
%VRA1_boundaryN(1,i,1)=sum(tempVRA1(~isnan(tempVRA1)&criteria)>limitVRA1);

%%
% % Fstates1=[mean(Fnrem1(criteria)) mean(Fnrem2(criteria)) mean(Fnrem3(criteria)) mean(FREM(criteria))]';
% % Fstates2=[mean(Fnrem1(criteria&tempN_arousals>0)) mean(Fnrem2(criteria&tempN_arousals>0)) mean(Fnrem3(criteria&tempN_arousals>0)) mean(FREM(criteria&tempN_arousals>0))]';
% %
% % tempN_arousals=LG_QualityInfo(:,1);
% % Fsupine1 = nanmean(Fsupinetemp);
% % Fsupine2 = nanmean(Fsupinetemp(tempN_arousals>0));

%% Plots LG, ArThres
if settings.plotfigs
    figure(1);
    
    subplot(2,6,3);
    xedges=0:0.1:1.2;
    h=histogram(LG1values,xedges,'normalization','probability','EdgeAlpha',0);
    box('off');
    xlim([min(xedges) max(xedges)]);
    %ylim([0 max([0.2;h.Values(:)])]);
    xlabel('LG1')
    hold('on');
    
    subplot(2,6,4);
    xedges=0:0.1:1.2;
    h=histogram(LGnvalues,xedges,'normalization','probability','EdgeAlpha',0);
    box('off');
    xlim([min(xedges) max(xedges)]);
    %ylim([0 max([0.2;h.Values(:)])]);
    xlabel('LGn');
    hold('on');
    
    subplot(2,6,5);
    xedges=0:2:30;
    h=histogram(delayvalues,xedges,'normalization','probability','EdgeAlpha',0);
    box('off');
    xlim([min(xedges) max(xedges)]);
    %ylim([0 max([0.2;h.Values(:)])]);
    xlabel('delay')
    hold('on');
    
    subplot(2,6,6);
    xedges=0:10:100;
    h=histogram(VRAvalues,xedges,'normalization','probability','EdgeAlpha',0);
    box('off');
    xlim([min(xedges) max(xedges)]);
    %ylim([0 max([0.2;h.Values(:)])]);
    xlabel('VRA')
    hold('on');
    
    plotarthresmodel=0;
    if plotarthresmodel
        subplot(2,6,9);
        xedges=(0.8:0.1:3)*100;
        h=histogram(ArThresvalues,xedges,'normalization','probability','EdgeAlpha',0);
        box('off');
        xlim([min(xedges) max(xedges)]);
        %ylim([0 max([0.2;h.Values(:)])]);
        xlabel('ArThres')
        hold('on');
    end
end

%get DataN
for i=1:Nperwindowvars
    eval([varlist{i} 'N=length(' varlist{i} 'values);']);
end

if getCIs
    Nbootstrap=50;
    %default
    for i=1:length(varlist)
        eval([varlist{i} 'CI=[-Inf;Inf];']);
    end
    for i=1:Nperwindowvars
        try
            eval([varlist{i} 'CI=bootci(Nbootstrap,@median,' varlist{i} 'values);']);
        catch me
        end
    end
    
    if 0&settings.plotfigs
        i=1;
        subplot(2,6,3);
        ylims=get(gca,'YLim');
        hold('on');
        plot(eval(varlist{i})*[1 1],ylims,'k');
        %plot(eval([varlist{i} 'CI(1)'])*[1 1],ylims,'color',[0.5 0.5 0.5]);
        %plot(eval([varlist{i} 'CI(2)'])*[1 1],ylims,'color',[0.5 0.5 0.5]);
        
        i=2;
        subplot(2,6,4);
        ylims=get(gca,'YLim');
        hold('on');
        plot(eval(varlist{i})*[1 1],ylims,'k');
        %plot(eval([varlist{i} 'CI(1)'])*[1 1],ylims,'color',[0.5 0.5 0.5]);
        %plot(eval([varlist{i} 'CI(2)'])*[1 1],ylims,'color',[0.5 0.5 0.5]);
        
        i=3;
        subplot(2,6,5);
        ylims=get(gca,'YLim');
        hold('on');
        plot(eval(varlist{i})*[1 1],ylims,'k');
        %plot(eval([varlist{i} 'CI(1)'])*[1 1],ylims,'color',[0.5 0.5 0.5]);
        %plot(eval([varlist{i} 'CI(2)'])*[1 1],ylims,'color',[0.5 0.5 0.5]);
        
        i=4;
        subplot(2,6,6);
        ylims=get(gca,'YLim');
        hold('on');
        plot(eval(varlist{i})*[1 1],ylims,'k');
        %plot(eval([varlist{i} 'CI(1)'])*[1 1],ylims,'color',[0.5 0.5 0.5]);
        %plot(eval([varlist{i} 'CI(2)'])*[1 1],ylims,'color',[0.5 0.5 0.5]);
        
        if plotarthresmodel
            i=5;
            subplot(2,6,9);
            ylims=get(gca,'YLim');
            hold('on');
            plot(eval(varlist{i})*[1 1],ylims,'k');
            %plot(eval([varlist{i} 'CI(1)'])*[1 1],ylims,'color',[0.5 0.5 0.5]);
            %plot(eval([varlist{i} 'CI(2)'])*[1 1],ylims,'color',[0.5 0.5 0.5]);
        end
    end
else
    Nbootstrap=0;
end
%% Summary
%ArThres(ArThres==0)=NaN;
ArThresPSG1=ArThres;

%%
Fstates3=[mean(WinT.FNREM1(criteriaUA)) mean(WinT.FNREM2(criteriaUA)) mean(WinT.FNREM3(criteriaUA)) mean(WinT.FREM(criteriaUA))]';
Fsupine3 = nanmean(Fsupinetemp(criteriaUA));

NwindowsForUA=sum(criteriaUA);
str=['Using ', num2str(NwindowsForUA), ' of ', num2str(length(criteriaUA)) ,' windows for UA calcs'];
disp(str);

for ii=(Nperwindowvars+1):length(varlist)
    eval([varlist{ii} 'N=NwindowsForUA;']);
end
fontsize_= SAfontsize;

%% UA phenotype using model drive
if 0
    subplot(1,3,1);
    %make subset table
    criteriaUAI = find(criteriaUA==1);
    criteriaRow = sum((BreathDataTable2.Win == criteriaUAI'),2)>0;
    BreathDataTable3 = BreathDataTable2(criteriaRow==1,:);
    
    [Vpassive,Vactive,VEdeciles,Vdrivedeciles,VpassiveCI,VactiveCI,Vcomp,VcompCI]=SummaryAnalysisOne_UAmodel(length(criteriaUAI),BreathDataTable3,ArThres,varlist,settings.plotfigs,settings);
end

%% UA Phenotype, gold standard drive

try
    
    ArThresActive=NaN;
    
    if NwindowsForUA>0
        if settings.plotfigs
            plot(1,1,'marker','none'); hold('on');
        end
        hold('on');
        %%[16:Pes, 17:Pmus, 18:VEpes, 19:Edi, 20:VEedi, 21:GGpeak 24:FlowEdi, 23:FlowPes]
        dividedrivebyttot=0;
        
        %criteria = N_events>=minNevents & minwake<maxwake2 & FREM<=maxFREM ; %& (~isnan(Pos)&supine==1)
        Nwindows=sum(criteriaUA);
        
        %% Wake data, used to normalize drive to eupneic flow (later)
        %calculate "g" = nanmedian(y(a1)./x(a1)); y = flow, x=drive, a1=wake breaths
        BreathDataTable3=BreathDataTable2;
        
        x=BreathDataTable3.(settings.DriveSignal);
        y = BreathDataTable3.VE;
        a = BreathDataTable3.AR3;
        hyp = BreathDataTable3.hypnog_B;
        win = floor(BreathDataTable3.UniqueID/1000);
        
        %need y, x, a1, win
        if settings.scoredarousalsinwake==0
            a=(a==1)|(hyp==4);
        end
        [ad]=howfarawayfromsleep(a,win);
        %find breaths during wakefulness and arousals
        minNwakebreaths = 50;
        hh=hist(ad,[1:11]); hh(end)=NaN; th=find(hh>minNwakebreaths,1,'last');
        threshold = min([4 th]);
        a1 = ad>threshold;
        g = nanmedian(y(a1)./x(a1));
        x_wake=nanmedian(x(a1));
        y_wake=nanmedian(y(a1));
        G_N = length(y(a1)./x(a1));
        
        %% GET BREATH-BY-BREATH SLEEP DATA
        
        
        win1 = floor(BreathDataTable2.UniqueID/1000);
        criteriaRow = sum((win1 == (find(criteriaUA==1))'),2)>0;
        BreathDataTable3 = BreathDataTable2(criteriaRow==1,:);
        t=BreathDataTable3.Time_start;
        x=BreathDataTable3.(settings.DriveSignal);
        xslope = BreathDataTable3.DriveSlope;
        
        y = BreathDataTable3.VE;
        a = BreathDataTable3.AR3;
        veup = BreathDataTable3.Veup;
        hyp = BreathDataTable3.hypnog_B;
        win = floor(BreathDataTable3.UniqueID/1000);
        t2=BreathDataTable3.Time_end - BreathDataTable3.Time_start;
        e=BreathDataTable3.E1;
        Etype=BreathDataTable3.Etype;
        nota = BreathDataTable3.notAR3;
        GGp = BreathDataTable3.GGpeak;
        GGt = BreathDataTable3.GGtonic;
        REMphasic = BreathDataTable3.REMphasic;
        REMtonic = BreathDataTable3.REMtonic;
        WS=BreathDataTable3.WakeSleep
        
        StableBreathing = BreathDataTable3.StableBreathing;
        
       
        if settings.constantEupnea %the single-value for Veupnea will stay constant regardless of selected states
            criteriaEupnea = WinT.LongestWake<=settings.maxwakeUA & WinT.FREM<=0; %use all data except wake and REM; should adapt with window length
            win1 = floor(BreathDataTable2.UniqueID/1000);
            criteriaRow = sum((win1 == (find(criteriaEupnea==1))'),2)>0;
            veup1 = BreathDataTable2{criteriaRow==1,'Veup'};
        end
        
        if dividedrivebyttot
            x=x./t2;
        end
        
        %normalize data
        if ~settings.normalizeusingconstantEupnea %eupnea stays constant (does not adapt to signal size), window by window
            y = y./veup; %%window-based normalization. note Y is unnormalized VE from LGfromflow
            y_ = nanmedian(veup);
        else
            if settings.constantEupnea
                y = y./nanmedian(veup1); %study-based normalization, but removes windows predominately wake/REM. note Y is unnormalized VE from LGfromflow
                y_ = nanmedian(veup1);
            else
                y = y./nanmedian(veup); %never actually used this, study-based normalization of all windows Veupnea
                y_ = nanmedian(veup);
            end
        end
        Veupnea = y_*60; %
       
        
        if 0 %historical way (uses selected criteria data to generate eupneic drive, seems odd)
            win1 = floor(BreathDataTable2.UniqueID/1000);
            criteriaRow = sum((win1 == (find(criteria==1))'),2)>0;
            BreathDataTable4 = BreathDataTable2(criteriaRow==1,:);
            y_n = nanmedian(BreathDataTable4.Veup);
            x_ = y_n/g; %xeup: drive eupneic value, veup divided by wake scale factor
            %this would artificially raise drive in %eupnea if V was lower during selected periods
        else %use the actual eupneic flow we just calculated
            y_n = y_;
            x_ = y_/g; %xeup: drive eupneic value, veup divided by wake scale factor
        end
        
        if ~exist('EupneicDriveIsUnity')
            EupneicDriveIsUnity=0;
        end
        if settings.DriveSignal=="VdriveModel"
            EupneicDriveIsUnity=1;
        end
        if EupneicDriveIsUnity==1
            x_=1;
        end
        
        x_n=x_
        
        %Arousal threshold
        
        %First find arthres for Vactive/Vcomp:
        settings2=[3 0 1 0 1 0 1 0.2];
        %1:ignorefirstXbreathsaftersleeponset [4]
        %2:swapbreathsforhigherdrive [N]
        %3:Nbreathsattributabletoarousal [N]
        %4:increasingdriveonly [N]
        %5:deleteifbeloweupnea [Y]
        %6:setaseupneaifbeloweupnea [N]
        %7: use median drive not ROC curves
        %8: not used, was relevant for ROC.
        
        criteriabreath = sum(hyp==hypok,2)>0; % criteria for sleep stages in breath analysis.
        xforAr = x;
        if 1 %new
            xforAr(criteriabreath==0)=NaN;
        end
        
        if ~isfield(settings,'ArThresOKBelowEupnea')
            ArThresOKBelowEupnea=1; %swapped to default 1, 20210712
        else
            ArThresOKBelowEupnea=settings.ArThresOKBelowEupnea;
        end
        if ArThresOKBelowEupnea
            settings2(5)=0;
        end
        
        [ArThres,arthres_N,~,ArThresvalues] = ArThresNew(a',xforAr',x_*1.05,win,settings2); %x_*1.05 means throw out xforAr rows if values are less than x_*1.05 [if settings2(5)==1]
        if 1
            [~,~,~,ArThresvalues_] = ArThresNew(a',xforAr',x_,win,[3 0 1 0 0 0 1 0.2]);
        end
        
        if length(ArThresvalues)>1 && Nbootstrap>0
            ArThresCI=bootci(Nbootstrap,@median,ArThresvalues)/x_*100;
        else
            ArThresCI=[NaN;NaN];
        end
        ArThresN = length(ArThresvalues);
        
        ArThresActive=ArThres;
        if ArThresActive<(1.05*x_),ArThresActive=1.05*x_; end
        %if ArThres(n)<(1.00*x_),ArThres(n)=1.00*x_; end %irrelevant since min value is 1.05Eupnea
        
        
        %% Remove wakefulness breaths for further analysis
        if scoredarousalsinwake
            Ikeep = (nota==1);
        else
            Ikeep = (nota==1&hyp~=4);
        end
        
        x = x(Ikeep);
        y = y(Ikeep);
        veup = veup(Ikeep);
        GGp = GGp(Ikeep);
        GGt = GGt(Ikeep);
        hyp = hyp(Ikeep);
        try
            xslope = xslope(Ikeep);
            Etype = Etype(Ikeep);
            StableBreathing = StableBreathing(Ikeep);
            WS=WS(Ikeep)
        catch me
        end
        
        try
            REMtonic = REMtonic(Ikeep);
            REMphasic = REMphasic(Ikeep);
        catch me
            disp('no REM tonic phasic data')
        end
        
        criteriabreath = ones(length(x),1)==1; %default
        
        if settings.breathlevelinclusion %excludes individual breaths if not in appropriate state as selected; new default
            %hypok = [0 1 2 3];
            criteriabreath = sum(hyp==hypok,2)>0; % criteria for sleep stages in breath analysis.
        end
        
        %settings.TonicREMonly=1;
        %settings.PhasicREMonly=1;
        disp(['N breaths: ' num2str(sum(criteriabreath))]);
        if state==5
            if isfield(settings,'PhasicREMonly') && settings.PhasicREMonly==1 %Phasic only
                NbreathsPhasic = sum(REMphasic==1)
                NbreathsTotal = sum(criteriabreath==1)
                criteriabreath = criteriabreath & (REMphasic==1);
            end
            
            if isfield(settings,'TonicREMonly') && settings.TonicREMonly==1 %Tonic REM only
                NbreathsTonic = sum(REMtonic==1)
                NbreathsTotal = sum(criteriabreath==1)
                criteriabreath = criteriabreath & (REMtonic==1);
            end
        end
        
        % fun new options :)
        if isfield(settings,'EndogramDriveslopecriteria')
            if settings.EndogramDriveslopecriteria>0
                criteriabreath = criteriabreath & (xslope>0);
            elseif settings.EndogramDriveslopecriteria<0
                criteriabreath = criteriabreath & (xslope<0);
            end
        end
        
        if isfield(settings,'EndogramEventcriteria')
            if settings.EndogramEventcriteria==1 %event only
                criteriabreath = criteriabreath & (Etype>0);
            elseif settings.EndogramEventcriteria==-1 %non event only
                criteriabreath = criteriabreath & (Etype==0);
            elseif settings.EndogramEventcriteria==2 %apnea
                criteriabreath = criteriabreath & (Etype==2);
            elseif settings.EndogramEventcriteria==3 %hypopnea
                criteriabreath = criteriabreath & (Etype==4);
            end
        end
        
        if isfield(settings,'StableBreathingcriteria')
            if settings.StableBreathingcriteria==0
                criteriabreath = criteriabreath & (StableBreathing==settings.StableBreathingcriteria);
            elseif settings.StableBreathingcriteria>0
                criteriabreath = criteriabreath & (StableBreathing>=settings.StableBreathingcriteria);
            end
        end
        
              if isfield(settings,'WScriteria')
            if settings.WScriteria==1
                criteriabreath = criteriabreath & (WS>=0.230);
            elseif settings.WScriteria==2
                criteriabreath = criteriabreath & (WS>=0.0123 & WS<0.230);
            elseif settings.WScriteria==3
                criteriabreath = criteriabreath & (WS<=0.0123 & WS>0);
            end
        end
        
        
        
        
        disp(['N breaths: ' num2str(sum(criteriabreath))]);
        x = x(criteriabreath);
        y = y(criteriabreath);
        GGp = GGp(criteriabreath);
        GGt = GGt(criteriabreath);
        veup = veup(criteriabreath);
        
        
        medianV = 100*median(y);
        
        x=x/x_;
        ArThres=ArThres/x_*100; %output to Table
        ArThresActive=ArThresActive/x_;
        
        if ~settings.normalizeusingconstantEupnea
            veonvdrive=y./x*100;
            veonvwake=y.*veup./y_wake*100;
            driveondrivewake=x.*x_/x_wake*100;
        else
            veonvdrive=y./x*100;
            veonvwake=y.*nanmedian(veup1)./y_wake*100;
            driveondrivewake=x.*x_/x_wake*100;
        end
        
        if length(x)<20
            length(x)
            error('not enough breaths left to analyze')
        end
        %NcileLimits=[5 95];
        %settings.Nciles=131;
        
        %NcileLimits=[2 80];
        %settings.Nciles=119;
        %
        %%   Endogram Analysis and Plot
        
        if ~isfield(settings,'SummaryEndogramFilt121')
            filt121 = 1;  %default on, since 20201027
        else
            filt121 = settings.SummaryEndogramFilt121;
        end
        
        ploton=settings.plotfigs;
        if ploton==1
            figure(1);
            subplot(1,3,1);
            hold on
        end
        
        if ~isfield(settings,'NcileLimits')
            NcileLimits=[];
        else
            NcileLimits=settings.NcileLimits;
        end
        
        if ~DriveincmH2OnotEupnea
            if settings.plotfigs
                plot([100 100],[0 150],'--','color',[0.7 0.7 0.7]);
                hold('on');
                plot([0 3*100],[100 100],'--','color',[0.7 0.7 0.7]);
            end
            [Vpassive,Vactive,VEdeciles,Vdrivedeciles,VpassiveCI,VactiveCI,VEdecilesUpper,VEdecilesLower]=VEVdriveAnalysis(100*x,100*y,100*ArThresActive,ploton,100,settings.Nciles,Nbootstrap,NcileLimits,filt121);
        else
            if settings.plotfigs
                plot([x_ x_],[0 150],'--','color',[0.7 0.7 0.7]);
                hold('on');
                plot([0 3*x_],[100 100],'--','color',[0.7 0.7 0.7]);
            end
            [Vpassive,Vactive,VEdeciles,Vdrivedeciles,VpassiveCI,VactiveCI,VEdecilesUpper,VEdecilesLower]=VEVdriveAnalysis(x*x_,100*y,ArThresActive*x_,ploton,x_,settings.Nciles,Nbootstrap,NcileLimits,filt121);
        end
        VpassiveCI=VpassiveCI(:);
        VactiveCI=VactiveCI(:);
        
        
        if settings.plotfigs
            set(gca,'fontsize',fontsize_)
            ylim([0 max([prctile(100*y,95) 100])]);
            ylabel('Ventilation, %eupnea');
            if ~DriveincmH2OnotEupnea
                xlim([0 max([prctile(100*x,95) 300])]);
                xlabel('Drive, %eupnea');
            else
                xlim([0 max([prctile(x*x_,95) 8*x_])]);
                xlabel('Drive, absolute units');
            end
            %hold('off');
            %try title(['AHI=',num2str(round(AHIstate),'%u')]); end
            
            subplot(2,6,9);
            xedges3=0.7:0.1:3;
            
            dx=0.1;
            h=histogram(100*ArThresvalues/x_,100*xedges3,'normalization','probability','EdgeAlpha',0);
            
            box('off');
            xlim(100*[min(xedges3) max(xedges3)]);
            %ylim([0 max([0.2;h.Values(:)])]);
            xlabel('ArThres');
            hold('on');
            
            
            subplot(2,6,10);
            xedges=0:0.1:1.4;
            xedges2=0:0.1:4.0;
            dx=0.1;
            h=histogram(100*y(x>(1-dx)&x<(1+dx)),100*xedges,'normalization','probability','EdgeAlpha',0);
            box('off');
            xlim(100*[min(xedges) max(xedges)]);
            ylim([0 max([0.2;h.Values(:)])]);
            xlabel('Vpassive');
            hold('on');
            
            
            subplot(2,6,11);
            h=histogram(100*y(x>(ArThresActive-dx)&x<(ArThresActive+dx)),100*xedges,'normalization','probability','EdgeAlpha',0);
            box('off');
            xlim(100*[min(xedges) max(xedges)]);
            ylim([0 max([0.2;h.Values(:)])]);
            xlabel('Vactive');
            hold('on');
            
            
        end
        if 0 & settings.plotfigs
            
            i=5;
            subplot(2,6,9);
            ylims=get(gca,'YLim');
            hold('on');
            plot(eval(varlist{i})*[1 1],ylims,'k');
            plot(eval([varlist{i} 'CI(1)'])*[1 1],ylims,'color',[0.5 0.5 0.5]);
            plot(eval([varlist{i} 'CI(2)'])*[1 1],ylims,'color',[0.5 0.5 0.5]);
            
            i=6;
            subplot(2,6,10);
            ylims=get(gca,'YLim');
            hold('on');
            plot(eval(varlist{i})*[1 1],ylims,'k');
            plot(eval([varlist{i} 'CI(1)'])*[1 1],ylims,'color',[0.5 0.5 0.5]);
            plot(eval([varlist{i} 'CI(2)'])*[1 1],ylims,'color',[0.5 0.5 0.5]);
            
            i=7;
            subplot(2,6,11);
            ylims=get(gca,'YLim');
            hold('on');
            plot(eval(varlist{i})*[1 1],ylims,'k');
            plot(eval([varlist{i} 'CI(1)'])*[1 1],ylims,'color',[0.5 0.5 0.5]);
            plot(eval([varlist{i} 'CI(2)'])*[1 1],ylims,'color',[0.5 0.5 0.5]);
            
            %     subplot(2,6,12);
            %     ylims=get(gca,'YLim');
            %     hold('on');
            %     plot(eval(varlist{i})*[1 1],ylims,'k');
            %     plot(eval([varlist{i} 'CI(1)'])*[1 1],ylims,'color',[0.5 0.5 0.5]);
            %     plot(eval([varlist{i} 'CI(2)'])*[1 1],ylims,'color',[0.5 0.5 0.5]);
        end
        Vcomp = Vactive-Vpassive;
        VcompCI = [VactiveCI(1)-VpassiveCI(2);VactiveCI(2)-VpassiveCI(1)];
        VcompCI = sort(VcompCI);
        
        figure(2); clf(2);
        varlist2 = {'veonvdrive','veonvwake','driveondrivewake'};
        xlimupper = [140 140 300];
        
        for i=1:length(varlist2)
            if 0
                h=histogram(eval(varlist2{i}),0:10:xlimupper(i),'normalization','probability','EdgeAlpha',0);
            else
                Data = eval(varlist2{i});
                dStep=10;
                Centers=5:dStep:(xlimupper(i)-dStep/2);
                Edges=(Centers(1)-dStep/2):dStep:(Centers(end)+dStep/2); Edges(end)=Inf;
                [h1,edges] = histcounts(Data,Edges);
                bar(Centers,h1/sum(h1),'EdgeAlpha',0,'FaceAlpha',0.6,'BarWidth',1);
                set(gca,'xtick',[0:50:300],'tickdir','out');
                
            end
            hold('on');
            eval([varlist2{i} '(isnan(' varlist2{i} '))=[];']);
            if Nbootstrap>0
                eval([varlist2{i} 'CI=bootci(Nbootstrap,@median,' varlist2{i} ');']);
            end
        end
        box('off');
        xlim([0 max(xlimupper)]);
        
        xlabel('Value, %');
        ylabel('Proportion of Sleep');
        
        ylims=get(gca,'YLim');
        i=1;
        plot(eval(['nanmedian(' varlist2{i} ')'])*[1 1],ylims,'k','color',[0 0.45 0.74]);
        i=2;
        plot(eval(['nanmedian(' varlist2{i} ')'])*[1 1],ylims,'k','color',[0.85 0.33 0.1]);
        i=3;
        plot(eval(['nanmedian(' varlist2{i} ')'])*[1 1],ylims,'k','color',[0.93 0.69 0.13]);
        
        h=legend('Actual Airflow / Intended Airflow','Actual Airflow / Wake Baseline','Intended Airflow / Wake Baseline');
        
        %% GG
        try
            
            ploton=settings.plotfigs;
            figure(4);
            
            
            if ~DriveincmH2OnotEupnea
                if settings.plotfigs
                    plot([100 100],[0 prctile(GGp,95)],'--','color',[0.7 0.7 0.7]);
                    hold('on');
                    %plot([0 3*100],[100 100],'--','color',[0.7 0.7 0.7]);
                end
                [GGppassive,GGpactive,GGpdeciles,GGpdeciles_drive,GGppassiveCI,GGpactiveCI,GGpdecilesUpper,GGpdecilesLower,GGpdecilesMean]=VEVdriveAnalysis(100*x,GGp,100*ArThresActive,ploton,100,settings.Nciles,Nbootstrap,settings.NcileLimits,filt121);
            else
                if settings.plotfigs
                    plot([x_ x_],[0 150],'--','color',[0.7 0.7 0.7]);
                    hold('on');
                    plot([0 3*x_],[100 100],'--','color',[0.7 0.7 0.7]);
                end
                [GGppassive,GGpactive,GGpdeciles,GGpdeciles_drive,GGppassiveCI,GGpactiveCI,GGpdecilesUpper,GGpdecilesLower,GGpdecilesMean]=VEVdriveAnalysis(x*x_,GGp,ArThresActive*x_,ploton,x_,settings.Nciles,Nbootstrap,settings.NcileLimits,filt121);
            end
            ylim([0 max([GGpdecilesUpper;5])]);
            set(gca,'fontsize',fontsize_);
            
            %ylim([0 max([prctile(100*y,95) 100])]);
            ylabel('GG peak, %max');
            if ~DriveincmH2OnotEupnea
                xlim([0 max([prctile(100*x,95) 300])]);
                xlabel('Drive, %eupnea');
            else
                xlim([0 max([prctile(x*x_,95) 8*x_])]);
                xlabel('Drive, absolute units');
            end
            
            hold('on');
            if ~DriveincmH2OnotEupnea
                if settings.plotfigs
                    plot([100 100],[0 prctile(GGp,95)],'--','color',[0.7 0.7 0.7]);
                    hold('on');
                    %plot([0 3*100],[100 100],'--','color',[0.7 0.7 0.7]);
                end
                [GGtpassive,GGtactive,GGtdeciles,GGtdeciles_drive,GGtpassiveCI,GGtactiveCI,GGtdecilesUpper,GGtdecilesLower]=VEVdriveAnalysis(100*x,GGt,100*ArThresActive,ploton,100,settings.Nciles,Nbootstrap,settings.NcileLimits,filt121);
            else
                if settings.plotfigs
                    plot([x_ x_],[0 150],'--','color',[0.7 0.7 0.7]);
                    hold('on');
                    plot([0 3*x_],[100 100],'--','color',[0.7 0.7 0.7]);
                end
                [GGtpassive,GGtactive,GGtdeciles,GGtdeciles_drive,GGtpassiveCI,GGtactiveCI,GGtdecilesUpper,GGtdecilesLower]=VEVdriveAnalysis(x*x_,GGt,ArThresActive*x_,ploton,x_,settings.Nciles,Nbootstrap,settings.NcileLimits,filt121);
            end
            
            ylim([0 max([GGpdecilesUpper;5])]);
            
            GGppassiveCI=GGppassiveCI(:);
            GGtactiveCI=GGtactiveCI(:);
            
            set(gca,'fontsize',fontsize_);
            
            %ylim([0 max([prctile(100*y,95) 100])]);
            ylabel('GG, %max');
            if ~DriveincmH2OnotEupnea
                xlim([0 max([prctile(100*x,95) 300])]);
                xlabel('Drive/Effort, %eupnea');
            else
                xlim([0 max([prctile(x*x_,95) 8*x_])]);
                xlabel('Drive/Effort, %absolute');
            end
            GGdata = [GGppassive;GGpactive;GGtpassive;GGtactive];
            
        catch me
            %disp(me.message);
            disp('GG analysis failed. Possibly no GG data available')
            
            GGdata = [NaN;NaN;NaN;NaN];
            GGpdeciles = NaN*ones(settings.Nciles,1);
            GGtdeciles = NaN*ones(settings.Nciles,1);
            close(4);
        end
        
        
        
        %% Event probability
                
        e=e(Ikeep);
        e=e(criteriabreath);
        [Epassive,Eactive,Edeciles,Edeciles_drive,EpassiveCI,EactiveCI,EdecilesUpper,EdecilesLower,EdecilesMean]=VEVdriveAnalysis(x*x_,e,ArThresActive*x_,0,x_,settings.Nciles,Nbootstrap,NcileLimits);
        
        if 1 %add colorbar to individual plot
            figure(1)
            subplot(1,3,1);
            DispEventLikelihoodOnEndogram(Vdrivedeciles,VEdeciles,1-EdecilesMean,0)
        end
        
        decileTime=t(Ikeep);
        decileTime=decileTime(criteriabreath);
        [Tpassive,Tactive,Tdeciles,Tdeciles_drive,TpassiveCI,TactiveCI,TdecilesUpper,TdecilesLower,TdecilesMean]=VEVdriveAnalysis(x*x_,decileTime,ArThresActive*x_,0,x_,settings.Nciles,Nbootstrap,NcileLimits);
        
        
    else
        Vpassive=NaN;
        Vactive=NaN;
        Vcomp=NaN;
        VpassiveCI=[-Inf;Inf];
        VactiveCI=[-Inf;Inf];
        VcompCI=[-Inf;Inf];
        GGdata = [NaN;NaN;NaN;NaN];
        VEdeciles = NaN*ones(settings.Nciles,1);
        Vdrivedeciles = NaN*ones(settings.Nciles,1);
        GGpdeciles = NaN*ones(settings.Nciles,1);
        GGtdeciles = NaN*ones(settings.Nciles,1);
        VEdecilesUpper= NaN*ones(settings.Nciles,1);
        VEdecilesLower= NaN*ones(settings.Nciles,1);
        x_ = NaN;
        y_ = NaN;
        x_wake = NaN;
        y_wake = NaN;
        Veupnea = NaN;
        medianV = NaN;
        EdecilesMean = NaN*ones(settings.Nciles,1);
        TdecilesMean = NaN*ones(settings.Nciles,1);
    end
    
catch me
    disp('failed UA analysis');
    disp(me.message);
    Vpassive=NaN;
    Vactive=NaN;
    Vcomp=NaN;
    VpassiveCI=[-Inf;Inf];
    VactiveCI=[-Inf;Inf];
    VcompCI=[-Inf;Inf];
    GGdata = [NaN;NaN;NaN;NaN];
    VEdeciles = NaN*ones(settings.Nciles,1);
    Vdrivedeciles = NaN*ones(settings.Nciles,1);
    GGpdeciles = NaN*ones(settings.Nciles,1);
    GGtdeciles = NaN*ones(settings.Nciles,1);
    VEdecilesUpper= NaN*ones(settings.Nciles,1);
    VEdecilesLower= NaN*ones(settings.Nciles,1);
    x_ = NaN;
    y_ = NaN;
    x_wake = NaN;
    y_wake = NaN;
    Veupnea = NaN;
    medianV = NaN;
    EdecilesMean = NaN*ones(settings.Nciles,1);
    TdecilesMean = NaN*ones(settings.Nciles,1);
    
end
%%
disp(['Eupneic VE, Vdrive sleep: ' num2str([y_ x_ ]) '; wake: ' num2str([y_wake x_wake]) ' (raw)']);

%% Save
criteriaset=[1 1 1 2 3 3 3 3 2];
clear Data DataCI Fstates Fsupine
for i=1:length(varlist)
    Data(i) = eval(varlist{i});
    if settings.getCIs
        DataCI(:,i) = eval([varlist{i} 'CI']);
    else
        DataCI(:,i) = NaN;
    end
    DataN(i) = eval([varlist{i} 'N']);
    Fstates(:,i) = eval(['Fstates' num2str(criteriaset(i))]);
    Fsupine(:,i) = eval(['Fsupine' num2str(criteriaset(i))]);
end
%pause

