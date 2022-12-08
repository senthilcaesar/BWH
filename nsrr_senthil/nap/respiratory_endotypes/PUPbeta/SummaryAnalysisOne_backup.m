function [Data,DataCI,DataN,varlist,AHItotal,Fstates,...
    Veupnea,medianV,VEdeciles,Vdrivedeciles,AHIdata,Fsupine,Twin,BootT,Evts] = SummaryAnalysisOne_backup(subj,settings)

global SAfontsize

% set some defaults in case of encountering an error
Data = NaN; DataCI = NaN; DataN = NaN; AHItotal = NaN; Fstates = NaN;
Veupnea = NaN; medianV = NaN; VEdeciles = NaN; Vdrivedeciles = NaN;
Fsupinetemp = NaN;

AHIdata = NaN; Twin = NaN;


%% Default for figures
set(groot,'defaultAxesfontname','arial narrow');
set(groot,'defaultFigureColor',[1 1 1]);
set(groot,'defaultAxesBox','off');
set(groot,'defaultAxesfontsize',SAfontsize);

%% Load PUPdata
varlist = {'LG1','LGn','delay','VRA','ArThres','Vpassive','Vactive','Vcomp'};
Nperwindowvars=5;

loadpath=[settings.AnalyzedDirectory settings.savename '_' num2str(subj)];
%load(loadpath,'AHIdata','SleepData','LGplusinfo','LG_QualityInfo','DataOut','fitQual','BreathDataTable'); %EvtsData contains rounded position codes, PositionData contains main positioncode for each window
% datatoload = load(loadpath,'AHIdata','SleepData','LGplusinfo','LG_QualityInfo','DataOut','fitQual','BreathDataTable','Evts','StoNData');
% datatoload = load(loadpath,'SleepData','LGplusinfo','LG_QualityInfo','DataOut','fitQual','BreathDataTable','StoNData');


% check to see if Evts is present--applicable to newly analyzed files where
% Evts contain AHIdata and all other info and EvtsData is not saved..

listOfVariables = who('-file', [loadpath '.mat']);
if ismember('Evts', listOfVariables) % returns true
    try
        datatoload = load(loadpath,'Evts','SleepData','LGplusinfo','LG_QualityInfo','DataOut','fitQual','BreathDataTable','StoNData');
    end
elseif ismember ('EvtsData', listOfVariables)
    try
        datatoload = load(loadpath,'AHIdata','SleepData','LGplusinfo','LG_QualityInfo','DataOut','fitQual','BreathDataTable','EvtsData','StoNData');
    end
else % for cfs--its Evtsdata
    try
        datatoload = load(loadpath,'AHIdata','SleepData','LGplusinfo','LG_QualityInfo','DataOut','fitQual','BreathDataTable','Evtsdata','StoNData');
    end
end

if sum(strcmp(fieldnames(datatoload),'Evts')) == 0
    if isfield(datatoload,'EvtsData')
        datatoload.Evts = datatoload.EvtsData;
        datatoload = rmfield(datatoload, 'EvtsData');
    elseif isfield(datatoload,'Evtsdata')
        datatoload.Evts = datatoload.Evtsdata;
        datatoload = rmfield(datatoload, 'Evtsdata');
    end
end

datatoloadlist = fieldnames(datatoload);
for i=1:length(datatoloadlist)
    if isfield(datatoload,datatoloadlist{i}) && iscell(datatoload.(datatoloadlist{i})) && size(datatoload.(datatoloadlist{i}),1)==1 ...
            && size(datatoload.(datatoloadlist{i}),2)==1  %exist, be a 1x1 cell, contains a cell
        datatoload.(datatoloadlist{i}) = datatoload.(datatoloadlist{i}){1}; 
    end
end      

cellfun(@(x,y) assignin('caller',x,y),fieldnames(datatoload),struct2cell(datatoload)) %Shifts all fields from from datatosave to current workspace:

try
    AHIdata=Evts.AHIdata;
end

if iscell(AHIdata)
    AHIdata=AHIdata{1};
end
                    
try
    Fnoise = StoNData.Fnoise;
    FnoiseT = array2table(Fnoise); FnoiseT.Properties.VariableNames = {'Fnoise1','Fnoise2','Fnoise3'};
    FnoiseT.Fnoise2Under10percent = FnoiseT.Fnoise2<0.1; % nnz(FnoiseT.Fnoise2Under10percent)
catch me
    disp('no StoNData available');
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

%% Comb
if settings.comb(1)
    disp('comb window selection');
    width = settings2.settings.windowlength;
    step=settings2.settings.WindowStep/60;
    
    settings.comb(2);
    winn = (1:size(LG_QualityInfo,1));
    starttime = (winn - 1).*step;
    endtime=starttime+width;
    temp = [starttime;endtime];
    swap = 1*(settings.comb(2)==2);
    temp2 = mod(temp + settings.comb(3)*swap,2*settings.comb(3))/(2*settings.comb(3));
    criteriacomb = (sum(temp2<=0.5)==2)';
else
    criteriacomb = ones(size(LG_QualityInfo,1),1);
end


%%
%AHIdata_=AHIdata;
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

clear LG1 LG2 LGn temp temp2 FVAL VRA VRA2 ArThres MeanEx


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
TimeOfNightCriteria = ones(length(SleepData(:,1)),1);
if settings.selecttimeofnight
    XTiles=settings.selecttimeofnight_XTiles;
    NthXTile=settings.selecttimeofnight_NthXTile;
    
    FwakePerWindow=SleepData(:,1);
    SleepWin1=find(FwakePerWindow<1,1,'first');
    SleepWinN=find(FwakePerWindow<1,1,'last');
    
    rangeWin=round((SleepWinN-SleepWin1+1)/XTiles);
    lowerWinNs=round(SleepWin1+((1:XTiles)-1)*rangeWin);
    upperWinNs=round(lowerWinNs+rangeWin-1);
    
    TimeOfNightCriteria=0*TimeOfNightCriteria;
    TimeOfNightCriteria(lowerWinNs(NthXTile):upperWinNs(NthXTile))=1;
end

%% Options
usemediannotmeanLG1=1;
% i=1;

%% Position
% old position data, from LG_QualityInfo
%temppositiondata=LG_QualityInfo(:,5); % this is the mode Pos for window
%temppositiondataNchanges=LG_QualityInfo(:,6);
%temppositiondataPercentInModePos=LG_QualityInfo(:,7);

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
Fwake =SleepData(:,1);
FN1=SleepData(:,3);
FN2=SleepData(:,4);
FN3=SleepData(:,5);
FREM=SleepData(:,6);
longestwake=SleepData(:,7);
Funknownsleep=round(1-sum(SleepData(:,[1 2 6]),2),2); %contains
    
containsunknownsleep=1*(Funknownsleep>0);

%%
% trim any extras, old bugfix
NN = size(LGplusinfo,1);
    longestwake(NN+1:end)=[];
    FREM(NN+1:end)=[];
    Funknownsleep(NN+1:end)=[];
    containsunknownsleep(NN+1:end)=[];
    FN1(NN+1:end)=[];
    FN2(NN+1:end)=[];
    FN3(NN+1:end)=[];
    TimeOfNightCriteria(NN+1:end)=[];
    Poscriteria(NN+1:end)=[];
    PosWin(NN+1:end)=[];
    % following lines added by DLM (found some studies that had long vars:
    % Fwake, Fsupinetemp, FnoiseT
    Fwake(NN+1:end)=[];
    Fsupinetemp(NN+1:end)=[];
    FnoiseT(NN+1:end,:)=[]; 
    
%% This criteria defines a minimum proportion of particular non-REM sleep states.
nremXcriteria= FN1>=minFnrem1 &...
    FN2>=minFnrem2 & ...
    FN2>=minFnrem2;

%% Other critiera Input, used and unused
N_events = LG_QualityInfo(:,2);

%GetRsquared:
%temp3=fitQual;
clear OneMinusRsq
for ii=1:size(fitQual,2)
    temp3pt1=fitQual{ii};
    if length(temp3pt1)>1
        OneMinusRsq(ii,:)=1-temp3pt1(2); %Element #2 of FitQual
    else
        OneMinusRsq(ii,:)=NaN;
    end
end
% Rsq=temp3(:,1);
FVAL=LGplusinfo(:,13);

MeanE=LG_QualityInfo(:,3); %is actual scored events. LG_QualityInfo=[N_arousals_used N_events mean(E) mean(E1) position_mode N_position_changes Percent_position EXITFLAG];
SStot=FVAL./OneMinusRsq;

tempN_arousals=LG_QualityInfo(:,1);

%% criteria for LG and UA phenotype
% reset the minNevents, special case analysis
if isfield(settings,'ResetMinNevents')==1 && settings.ResetMinNevents==1
    minNeventsLG = 0;
    minNeventsUA = 0;
else % normal
    minNeventsLG = 1; %1
    minNeventsUA = 1; %1
end

maxwakethresLG = 30; % shouldn't this come from the spreadsheet
maxwakethresUA = 300; % more tolerant to more wake for UA measures

criteriaAll = criteriacomb & ...
    Poscriteria==1 & ...
    nremXcriteria==1 & ...
    FREM<=maxFREM &...
    FREM>=minFrem & ...
    containsunknownsleep==0 & ...  %recently added containsunknownsleep
    TimeOfNightCriteria==1;
if settings.verbose
    disp(['criteriaAll = ', num2str(nnz(criteriaAll)), ' / ', num2str(length(criteriaAll))]);
end

criteriaUA = criteriaAll & ...
    N_events>=minNeventsUA & ... % nnz(isfinite(N_events))
    longestwake<=maxwakethresUA;
if settings.verbose
    disp(['criteriaUA = ', num2str(nnz(criteriaUA)), ' / ', num2str(length(criteriaUA))]);
end

% criteria for LG and Arousal Threshold
criteriaLG = criteriaAll & ...
    N_events>=minNeventsLG & ...
    longestwake<=maxwakethresLG;  
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
% tempLG1=LGplusinfo(:,7);
% FisNaN(1)=sum(isnan(tempLG1))/length(tempLG1);

% temp=LGplusinfo;

Twin = array2table(LGplusinfo);
Twin.Properties.VariableNames = {'Time','LG0','tau','empty','LGn','Tn','LG1','LG2','delay','VRA','empty2','ArThres','MSE','Ttotmean','TtotSD','TtotMedian','TtotIQR'};
Twin(:,{'empty','empty2','Ttotmean','TtotSD','TtotMedian','TtotIQR'})=[];

% tempLG1=LGplusinfo(:,7); tempLGn=LGplusinfo(:,5);  tempLG2=LGplusinfo(:,8); tempTn=LGplusinfo(:,6);
% tempLG0=LGplusinfo(:,2); temptau=LGplusinfo(:,3);  tempdelay=LGplusinfo(:,9);
% tempVRA1=LGplusinfo(:,10); tempVRA2=LGplusinfo(:,11); tempArThres=LGplusinfo(:,12);

%only works for 1 time constant:
alpha=(Twin.LG2./Twin.LG1).^2;
beta=(1-alpha)./(4*alpha-1);
Twin.LG0p17=Twin.LG1.*((1+beta)./(1+beta*(1/6)^2)).^0.5;
Twin.LG0p33=Twin.LG1.*((1+beta)./(1+beta*(1/3)^2)).^0.5;
Twin.LG0p67=Twin.LG1.*((1+beta)./(1+beta*(2/3)^2)).^0.5;
Twin.LG1p33=Twin.LG1.*((1+beta)./(1+beta*(4/3)^2)).^0.5;
Twin.LG1p67=Twin.LG1.*((1+beta)./(1+beta*(5/3)^2)).^0.5;

Twin.FVAL=FVAL;
Twin.MeanE=MeanE;
Twin.OneMinusRsq=OneMinusRsq;
Twin.Narousals=tempN_arousals;

Twin.VRA(Twin.Narousals<1)=NaN;

Twin.criteriaLG = criteriaLG;
Twin.N_events = N_events;
Twin.Poscriteria = Poscriteria;
Twin.nremXcriteria = nremXcriteria;
Twin.TimeOfNightCriteria = TimeOfNightCriteria;
Twin.longestwake  =longestwake;
Twin.criteriacomb = criteriacomb;
Twin.Fwake = Fwake;
Twin.FN1 = FN1;
Twin.FN2 = FN2;
Twin.FN3 = FN3;
Twin.FREM = FREM;
Twin.Fsupine = Fsupinetemp;
try
    Twin = [Twin FnoiseT];
catch ErrorConcatenatingWindowTables
    disp(ErrorConcatenatingWindowTables.message);
end
% nnz(FnoiseT.Fnoise2Under10percent)
%% Calculate summary per criteria
criteriaLGI = find(criteriaLG==1);
[LG1,LGn,delay,VRA,ArThres,LG1N,LGnN,delayN,VRAN,ArThresN,...
    Fstates1,Fsupine1,Fstates2,Fsupine2,...
    LG1values,LGnvalues,delayvalues,VRAvalues,ArThresvalues]=...
    SummaryAnalysisOne_LGmodel(Twin(criteriaLGI,:),varlist);

%% DLM: option here to replace the LG1 etc outputs with LG0, Tau and Tn
% This is so that we can calculate LG at any frequency, post hoc.
% originally, this was set up as an option to calc summary per criteria 
% for LG4 only (i.e. for infants), but is now being modifed to anyLG
% (the second half of this option is at bottom, where 'Data' is returned)
if isfield(settings,'coreLGvals') && settings.coreLGvals==1
    try
        if 0 % version 1
            %LG0_recalc = LG1.*(1+(2.*pi.*Tau./60).^2).^0.5; % redo LG0 if req'd
            LG4_perWindow = Twin.LG0./(sqrt(1+(((2.*pi.*(Twin.tau./60).*(60./15))).^2)));
            LG4 = nanmedian(LG4_perWindow(criteriaLGI));
        else  % version 2
            if isempty(criteriaLGI)
                LG0 = NaN; Tau = NaN; Tn = NaN; % no data, set as NaN
            else
                LG0 = nanmedian(Twin.LG0(criteriaLGI));
                Tau = nanmedian(Twin.tau(criteriaLGI));
                Tn = nanmedian(Twin.Tn(criteriaLGI));
            end
        end
    catch FaultCollectingCoreLGvals
        LG0 = NaN; Tau = NaN; Tn = NaN; % no data, set as NaN
    end
end


%% Plots LG, ArThres
%******************************************************************
% Plot histograms of phenotypes derived
% from the original LGfromPSG routine. This same plot is continued
% below with confidence intervals (if option selected) and UA anatomy
% phenotypes
%******************************************************************
if settings.plotfigs
    figure(1);
    %fig = gcf;
    %fig.Position = [ 2000 100 1400 420 ]; %%SS removed, causing trouble
    
    subplot(2,6,3);
    xedges=0:0.1:1.2;
    h=histogram(LG1values,xedges,'normalization','probability','EdgeAlpha',0); %facecolor,[0.00 0.45 0.74] --default blue
    box('off');
    xlim([min(xedges) max(xedges)]);
    ylim([0 max([0.2;h.Values(:)])]);
    xlabel('Loop gain, LG1');
    ylabel('Frequency');
    
    subplot(2,6,4);
    xedges=0:0.1:1.2;
    h=histogram(LGnvalues,xedges,'normalization','probability','EdgeAlpha',0);
    box('off');
    xlim([min(xedges) max(xedges)]);
    ylim([0 max([0.2;h.Values(:)])]);
    xlabel('Loop gain, LGn');
    
    subplot(2,6,5);
    xedges=0:2:30;
    h=histogram(delayvalues,xedges,'normalization','probability','EdgeAlpha',0);
    box('off');
    xlim([min(xedges) max(xedges)]);
    ylim([0 max([0.2;h.Values(:)])]);
    xlabel('Delay (s)')
    
    subplot(2,6,6);
    xedges=0:10:100;
    h=histogram(VRAvalues,xedges,'normalization','probability','EdgeAlpha',0);
    box('off');
    xlim([min(xedges) max(xedges)]);
    ylim([0 max([0.2;h.Values(:)])]);
    xlabel('VRA (%)')
    
    subplot(2,6,9);
    xedges=(0.8:0.1:3)*100;
    h=histogram(ArThresvalues,xedges,'normalization','probability','EdgeAlpha',0);
    box('off');
    xlim([min(xedges) max(xedges)]);
    ylim([0 max([0.2;h.Values(:)])]);
    xlabel('Arousal Thres. (%)');
    ylabel('Frequency');
end

%get DataN
% for i=1:Nperwindowvars
%     eval([varlist{i} 'N=length(' varlist{i} 'values);']);
% end

if settings.getCIs
    settings.Nbootstrap=50;
    %default
    for ii=1:length(varlist)
        eval([varlist{ii} 'CI=[-Inf;Inf];']);
    end
    for ii=1:Nperwindowvars
        try
            eval([varlist{ii} 'CI=bootci(settings.Nbootstrap,@median,' varlist{ii} 'values);']);
        catch me
        end
    end
    
    if settings.plotfigs
        figure(1);
        ii=1;
        subplot(2,6,3);
        ylims=get(gca,'YLim');
        hold('on');
        plot(eval(varlist{ii})*[1 1],ylims,'k');
        plot(eval([varlist{ii} 'CI(1)'])*[1 1],ylims,'color',[0.5 0.5 0.5]);
        plot(eval([varlist{ii} 'CI(2)'])*[1 1],ylims,'color',[0.5 0.5 0.5]);
        
        ii=2;
        subplot(2,6,4);
        ylims=get(gca,'YLim');
        hold('on');
        plot(eval(varlist{ii})*[1 1],ylims,'k');
        plot(eval([varlist{ii} 'CI(1)'])*[1 1],ylims,'color',[0.5 0.5 0.5]);
        plot(eval([varlist{ii} 'CI(2)'])*[1 1],ylims,'color',[0.5 0.5 0.5]);
        
        ii=3;
        subplot(2,6,5);
        ylims=get(gca,'YLim');
        hold('on');
        plot(eval(varlist{ii})*[1 1],ylims,'k');
        plot(eval([varlist{ii} 'CI(1)'])*[1 1],ylims,'color',[0.5 0.5 0.5]);
        plot(eval([varlist{ii} 'CI(2)'])*[1 1],ylims,'color',[0.5 0.5 0.5]);
        
        ii=4;
        subplot(2,6,6);
        ylims=get(gca,'YLim');
        hold('on');
        plot(eval(varlist{ii})*[1 1],ylims,'k');
        plot(eval([varlist{ii} 'CI(1)'])*[1 1],ylims,'color',[0.5 0.5 0.5]);
        plot(eval([varlist{ii} 'CI(2)'])*[1 1],ylims,'color',[0.5 0.5 0.5]);
        
        ii=5;
        subplot(2,6,9);
        ylims=get(gca,'YLim');
        hold('on');
        plot(eval(varlist{ii})*[1 1],ylims,'k');
        plot(eval([varlist{ii} 'CI(1)'])*[1 1],ylims,'color',[0.5 0.5 0.5]);
        plot(eval([varlist{ii} 'CI(2)'])*[1 1],ylims,'color',[0.5 0.5 0.5]);
    end
else
    Nbootstrap=0;
end


%% UA phenotype using model drive

NwindowsForUA=sum(criteriaUA);
str=['Using ', num2str(NwindowsForUA), ' of ', num2str(length(criteriaUA)) ,' windows for UA calcs'];
disp(str);

for ii=(Nperwindowvars+1):length(varlist)
    eval([varlist{ii} 'N=NwindowsForUA;']);
end

%move thiese inside SummaryAnalysisOne_UAmodel:
Fstates3=[mean(FN1(criteriaUA)) mean(FN2(criteriaUA)) mean(FN3(criteriaUA)) mean(FREM(criteriaUA))]';
Fsupine3 = nanmean(Fsupinetemp(criteriaUA));

criteriaUAI = find(criteriaUA==1);

[Vpassive,Vactive,VEdeciles,Vdrivedeciles,VpassiveCI,VactiveCI,Vcomp,VcompCI]=SummaryAnalysisOne_UAmodel(DataOut(criteriaUAI),ArThres,varlist,settings.plotfigs,settings);


%% BootStrapping
if isfield(settings,'Boot') && settings.Boot >0  % settings.Boot has to be the num of iterations required
    settings_ = settings;
    settings_.getCIs=0;
    settings_.Nbootstrap=0;
    
    I = find(criteriaEither==1);
    for j=1:settings.Boot
        I2 = datasample(I,length(I));
        criteriaUAIBoot = I2;
        criteriaUAIBoot(criteriaUA(I2)==0)=[];
        criteriaLGIBoot = I2;
        criteriaLGIBoot(criteriaLG(I2)==0)=[];
        
        [LG1Boot(j),LGnBoot(j),delayBoot(j),VRABoot(j),ArThresBoot(j),LG1Boot_N(j),LGnBoot_N(j),delayBoot_N(j),VRABoot_N(j),ArThresBoot_N(j),Fstates1Boot(j,:),Fsupine1Boot(j),Fstates2Boot(j,:),Fsupine2Boot(j)]=SummaryAnalysisOne_LGmodel(Twin(criteriaLGIBoot,:),varlist); %,Fstates1Boot(j),Fsupine1Boot(j),Fstates2Boot(j),Fsupine2Boot(j)
        
        [VpassiveBoot(j),VactiveBoot(j),~,~,~,~,VcompBoot(j)]=SummaryAnalysisOne_UAmodel(DataOut(criteriaUAIBoot),ArThresBoot(j),varlist,0,settings_); %also output Fstates3 Fsupine3
        
        Fstates3Boot(j,:)=[mean(FN1(criteriaUAIBoot)) mean(FN2(criteriaUAIBoot)) mean(FN3(criteriaUAIBoot)) mean(FREM(criteriaUAIBoot))]';
        Fsupine3Boot(j) = nanmean(Fsupinetemp(criteriaUAIBoot));
        
        VpassiveBoot_N(j) = length(criteriaUAIBoot);
        VactiveBoot_N(j) = length(criteriaUAIBoot);
        VcompBoot_N(j) = length(criteriaUAIBoot);
    end
    
    BootT=table(subj*ones(length(LG1Boot(:)),1),LG1Boot(:),LGnBoot(:),delayBoot(:),VRABoot(:),ArThresBoot(:),VpassiveBoot(:),VactiveBoot(:),VcompBoot(:),...
        LG1Boot_N(:),LGnBoot_N(:),delayBoot_N(:),VRABoot_N(:),VpassiveBoot_N(:),VactiveBoot_N(:),VcompBoot_N(:),...
        Fsupine1Boot(:),Fsupine2Boot(:),Fsupine3Boot(:),...
        Fstates1Boot(:,1),Fstates1Boot(:,2),Fstates1Boot(:,3),Fstates1Boot(:,4),...
        Fstates2Boot(:,1),Fstates2Boot(:,2),Fstates2Boot(:,3),Fstates2Boot(:,4),...
        Fstates3Boot(:,1),Fstates3Boot(:,2),Fstates3Boot(:,3),Fstates3Boot(:,4));
    
    BootT.Properties.VariableNames={'Subject','LG1Boot','LGnBoot','delayBoot','VRABoot','ArThresBoot','VpassiveBoot',...
        'VactiveBoot','VcompBoot','LG1Boot_N','LGnBoot_N','delayBoot_N','VRABoot_N','VpassiveBoot_N','VactiveBoot_N','VcompBoot_N'...
        'FsupLoopGain','FsupArTh','FsupUA',...
        'FN1LoopGain','FN2LoopGain','FN3LoopGain','FREMLoopGain',...
        'FN1ArTh','FN2ArTh','FN3ArTh','FREMArTh',...
        'FN1UA','FN2UA','FN3UA','FREMUA'};
    % fitglme(T,'LG1 ~ FN1 + FN3 + FREM + Flateral + (1 | Subj)')
    % fitglme(T,'LG1 ~ FN1*Sex + FN3*Sex + FREM*Sex + Flateral*Sex + (1 | Subj)')
    
else
    BootT=[];
end

%% Save
criteriaset=[1 1 1 2 2 3 3 3];
% this criteriaset variable says:
%   items 1-3 on varlist (i.e. LG1, LGn and delay) use criteria 1,
%   items 4-5 use criteria 2,
%   and items 6-8 use criteria 3.
clear Data DataCI Fstates Fsupine
for ii=1:length(varlist)
    Data(ii) = eval(varlist{ii});
    if settings.getCIs
        DataCI(:,ii) = eval([varlist{ii} 'CI']); %why are these Inf? %% fixed
    else
        DataCI(:,ii) = NaN;
    end
    DataN(ii) = eval([varlist{ii} 'N']);
    Fstates(:,ii) = eval(['Fstates' num2str(criteriaset(ii))]);
    Fsupine(:,ii) = eval(['Fsupine' num2str(criteriaset(ii))]);
end

% returning the following:
% Data,DataCI,DataN,varlist,AHItotal,Fstates,Veupnea,medianV,VEdeciles,Vdrivedeciles,AHIdata,Fsupine

%% and, ONLY if doing it, overwrite LG1, LGn and delay with the core LG vals
% this is fully dodgy, but this hack is really not intended for mainstream
if isfield(settings,'coreLGvals') && settings.coreLGvals==1
    try
        if 0 % version 1
            Data(1) = LG4; % hijack the position for LG1, and replace with LG4
        else % version 2
            % overwrite LG1 and  with LG0 and Tau.
            Data(1) = LG0;  % overwrites LG1
            Data(2) = Tau;  % overwrites LGn
            Data(3) = Tn;   % overwrites delay
            varlist = {'LG0','Tau','Tn','VRA','ArThres','Vpassive','Vactive','Vcomp'};
        end
    catch FaultAddingCoreLGvals
        % write out NaN values
        Data(1) = NaN; Data(2) = NaN; Data(3) = NaN;
    end
end