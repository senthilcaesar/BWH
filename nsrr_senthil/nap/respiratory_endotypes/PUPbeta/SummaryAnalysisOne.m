function [Data,DataCI,DataN,varlist,AHItotal,Fstates,...
    Veupnea,medianV,VEdeciles,Vdrivedeciles,AHIdata,Fsupine,WinT,BootT,Evts,SummaryLargeT] = SummaryAnalysisOne(subj,settings)

global SAfontsize

% set some defaults in case of encountering an error
Data = NaN; DataCI = NaN; DataN = NaN; AHItotal = NaN; Fstates = NaN;
Veupnea = NaN; medianV = NaN; VEdeciles = NaN; Vdrivedeciles = NaN;
Fsupinetemp = NaN;

AHIdata = NaN;

SummaryLargeT = table();
SummaryLargeT.FNoisyLG = NaN;
SummaryLargeT.FnEventsLG = NaN; 
SummaryLargeT.FFailOtherLG = NaN; 
SummaryLargeT.FAnalyzableLG = NaN;
SummaryLargeT.FNoisyUA = NaN;
SummaryLargeT.FnEventsUA = NaN;
SummaryLargeT.FFailOtherUA = NaN; 
SummaryLargeT.FAnalyzableUA = NaN;

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
% if iscell(BreathDataTable) && size(BreathDataTable,1)==1 && size(BreathDataTable,2)==1 && ~isa(BreathDataTable{1},'table') %can't find other code that already does this, harmonize later.
%     BreathDataTable=BreathDataTable{1};
% end

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

if isfield(settings,'PrintSummaryAnalysisSetup') && settings.PrintSummaryAnalysisSetup==1 % debugging the changing selectstate variable
    try disp(['selectstate: ', num2str(settings.selectstate)]); catch; end
    try disp(['selectstateEventAnalysis: ', num2str(settings.selectstateEventAnalysis)]); catch; end
    try disp(['selectposition: ', num2str(settings.selectposition)]); catch; end
    try disp(['analysisRange: ', num2str(settings.analysisRange)]); catch; end
    % try disp(['Mrange: ', num2str(settings.Mrange)]); catch; end
end


%% Analyse LG, ArThres
clear LG1 LG2 LGn temp temp2 FVAL VRA VRA2 ArThres MeanEx


%% State, 1=nrem1, 2=nrem2, 3=nrem3, 4 = all nrem, 5 rem, 8 ignore state
%defaults
minFnrem1=-Inf;
minFnrem2=-Inf;
minFnrem3=-Inf;
minFnrem23=-Inf;
switch settings.selectstate
    case 1
        maxFREM=0;
        minFrem=-Inf;
        minFnrem1=0.5;%settings.statetolerance;
        AHIstate=AHIdata(88); %not correct
        DurationInState=AHIdata(88-7); %not correct
        PlotString = 'N1';
        hypok = [2];
    case 2
        maxFREM=0;
        minFrem=-Inf;
        minFnrem2=0.5; %settings.statetolerance;
        AHIstate=AHIdata(88); %not correct
        DurationInState=AHIdata(88-7); %not correct
        PlotString = 'N2';
        hypok = [1];
    case 3
        maxFREM=0;
        minFrem=-Inf;
        minFnrem3=0.5; %settings.statetolerance;
        AHIstate=AHIdata(88); %not correct
        DurationInState=AHIdata(88-7); %not correct
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
        AHIstate=AHIdata(96); %not correct
        DurationInState=AHIdata(96-7); %not correct
        PlotString = 'REM';
        hypok = [3];
    case -1 % no nREM1 %new option 20211130
        maxFREM=0;
        minFrem=-Inf;
        minFnrem23=0.5; %settings.statetolerance;
        AHIstate=AHIdata(88); %not correct
        DurationInState=AHIdata(88-7); %not correct
        PlotString = 'N2andN3';
        hypok = [0 1];
    case 8 % at the moment this does not restrict at all, and could include periods which are predominantly wake unless a wake threshold is applied elsewhere.
        maxFREM=Inf;
        minFrem=-Inf;
        AHIstate=AHIdata(80); %not correct
        DurationInState=AHIdata(80-7); %not correct
        PlotString = 'ALL';
        hypok = [0 1 2 3];
    case 9 % at the moment this does not restrict at all, and could include periods which are predominantly wake unless a wake threshold is applied elsewhere.
        maxFREM=Inf;
        minFrem=-Inf;
        AHIstate=AHIdata(80); %not correct
        DurationInState=AHIdata(80-7); %not correct
        PlotString = 'ALL';
        hypok = [0 1 2 3 4];
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

%% Options
usemediannotmeanLG1=1;
% i=1;

%% Position
% old position data, from LG_QualityInfo
%temppositiondata=LG_QualityInfo(:,5); % this is the mode Pos for window
%temppositiondataNchanges=LG_QualityInfo(:,6);
%temppositiondataPercentInModePos=LG_QualityInfo(:,7);

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
    WinT.FNREM3>=minFnrem3 & ...
    (WinT.FNREM2 + WinT.FNREM3)>=minFnrem23; %last line is for new option to exclude stage 1 data

%% criteria for LG and UA phenotype
minNeventsLG = settings.minNeventsLG; %default 1
minNeventsUA = settings.minNeventsUA; %has been default 1, was default 0 for a long time previously
maxwakethresLG = settings.maxwakethresLG; % 30
maxwakethresUA = settings.maxwakethresUA; % 300; more tolerant to more wake for UA measures

try %added SS 11/2/2021
    CPAPoffCriteria = WinT.CPAPoff==1;
catch
    CPAPoffCriteria = zeros(1,height(WinT));
    disp('warning: CPAPoff detection failed, assuming CPAP is off always');
end

criteriaAll = criteriacomb & ...
    CPAPoffCriteria==1 & ...        %added SS 11/2/2021
    Poscriteria==1 & ...
    nremXcriteria==1 & ...
    WinT.FREM<=maxFREM & ...
    WinT.FREM>=minFrem & ...
    containsunknownsleep==0 & ...  %recently added contains unknownsleep
    TimeOfNightCriteria==1;
if settings.verbose
    disp(['criteriaAll = ', num2str(nnz(criteriaAll)), ' / ', num2str(length(criteriaAll))]);
end

if isfield (settings,'AutoScoredEventsAllCentral')& settings.AutoScoredEventsAllCentral==1
    WinT.Nevents=ones(height(WinT),1);
    
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

SummaryLargeT.FNoisyLG = sum(criteriaAll & WinT.LongestWake<=maxwakethresLG & WinT.FNoise2>0.1)/sum(criteriaAll & WinT.LongestWake<=maxwakethresLG);
SummaryLargeT.FnEventsLG = sum(criteriaAll & WinT.LongestWake<=maxwakethresLG & WinT.Nevents<minNeventsLG)/sum(criteriaAll & WinT.LongestWake<=maxwakethresLG);
SummaryLargeT.FFailOtherLG = sum(criteriaAll & WinT.LongestWake<=maxwakethresLG & isnan(WinT.Nevents))/sum(criteriaAll & WinT.LongestWake<=maxwakethresLG);
SummaryLargeT.FAnalyzableLG = sum(criteriaAll & WinT.LongestWake<=maxwakethresLG & WinT.Nevents>=minNeventsLG)/sum(criteriaAll & WinT.LongestWake<=maxwakethresLG);
SummaryLargeT.FNoisyUA = sum(criteriaAll & WinT.LongestWake<=maxwakethresUA & WinT.FNoise2>0.1)/sum(criteriaAll & WinT.LongestWake<=maxwakethresUA);
SummaryLargeT.FnEventsUA = sum(criteriaAll & WinT.LongestWake<=maxwakethresUA & WinT.Nevents<minNeventsUA)/sum(criteriaAll & WinT.LongestWake<=maxwakethresUA);
SummaryLargeT.FFailOtherUA = sum(criteriaAll & WinT.LongestWake<=maxwakethresUA & isnan(WinT.Nevents))/sum(criteriaAll & WinT.LongestWake<=maxwakethresUA);
SummaryLargeT.FAnalyzableUA = sum(criteriaAll & WinT.LongestWake<=maxwakethresUA & WinT.Nevents>=minNeventsUA)/sum(criteriaAll & WinT.LongestWake<=maxwakethresUA)


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
% nnz(FnoiseT.Fnoise2Under10percent)

%% Calculate summary per criteria
criteriaLGI = find(criteriaLG==1);
[LG1,LGn,delay,VRA,ArThres,LG1N,LGnN,delayN,VRAN,ArThresN,...
    Fstates1,Fsupine1,Fstates2,Fsupine2,...
    LG1values,LGnvalues,delayvalues,VRAvalues,ArThresvalues]=...
    SummaryAnalysisOne_LGmodel(WinT(criteriaLGI,:),varlist);

%% DLM: option here to replace the LG1 etc outputs with LG0, Tau and Tn
% This is so that we can calculate LG at any frequency, post hoc.
% originally, this was set up as an option to calc summary per criteria 
% for LG4 only (i.e. for infants), but is now being modifed to anyLG
% (the second half of this option is at bottom, where 'Data' is returned)
if isfield(settings,'coreLGvals') && settings.coreLGvals==1
    try
        if 0 % version 1
            %LG0_recalc = LG1.*(1+(2.*pi.*Tau./60).^2).^0.5; % redo LG0 if req'd
            %LG4_perWindow = WinT.LG0./(sqrt(1+(((2.*pi.*(WinT.tau./60).*(60./15))).^2)));
            LG4 = nanmedian(WinT.LG4(criteriaLGI)); %added above
        else  % version 2
            if isempty(criteriaLGI)
                LG0 = NaN; Tau = NaN; Tn = NaN; % no data, set as NaN
            else
                LG0 = nanmedian(WinT.LG0(criteriaLGI));
                Tau = nanmedian(WinT.tau(criteriaLGI));
                Tn = nanmedian(WinT.Tn(criteriaLGI));
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
    try % machine specific figure placement
        if startsWith(getenv('COMPUTERNAME'),'DESKTOP-U4DUN7N')
            fig = gcf; fig.Units = 'Inches'; fig.Position = [ 6 6 14 4.5 ]; 
        end
        
        if startsWith(getenv('COMPUTERNAME'),'MU00188859')
            fig = gcf; fig.Units = 'Inches'; fig.Position = [ 6 6 14 4.5 ]; 
        end
        
    catch
    end
    %fig = gcf;
    %fig.Position = [ 2000 100 1400 420 ]; %%SS removed, causing trouble
    
    subplot(2,6,3);
    hold on %for compare
    xedges=0:0.1:1.2;
    h=histogram(LG1values,xedges,'normalization','probability','EdgeAlpha',0); %facecolor,[0.00 0.45 0.74] --default blue
    box('off');
    xlim([min(xedges) max(xedges)]);
    ylim([0 max([0.2;h.Values(:)])]);
    xlabel('Loop gain, LG1');
    ylabel('Frequency');
    hold off
    
    subplot(2,6,4);
    hold on %for compare
    xedges=0:0.1:1.2;
    h=histogram(LGnvalues,xedges,'normalization','probability','EdgeAlpha',0);
    box('off');
    xlim([min(xedges) max(xedges)]);
    ylim([0 max([0.2;h.Values(:)])]);
    xlabel('Loop gain, LGn');
    hold off
    
    subplot(2,6,5);
    hold on %for compare
    xedges=0:2:30;
    h=histogram(delayvalues,xedges,'normalization','probability','EdgeAlpha',0);
    box('off');
    xlim([min(xedges) max(xedges)]);
    ylim([0 max([0.2;h.Values(:)])]);
    xlabel('Delay (s)')
    hold off
    
    subplot(2,6,6);
    hold on %for compare
    xedges=0:10:100;
    h=histogram(VRAvalues,xedges,'normalization','probability','EdgeAlpha',0);
    box('off');
    xlim([min(xedges) max(xedges)]);
    ylim([0 max([0.2;h.Values(:)])]);
    xlabel('VRA (%)')
    hold off
    
    subplot(2,6,9);
    hold on %for compare
    xedges=(0.8:0.1:3)*100;
    h=histogram(ArThresvalues,xedges,'normalization','probability','EdgeAlpha',0);
    box('off');
    xlim([min(xedges) max(xedges)]);
    ylim([0 max([0.2;h.Values(:)])]);
    xlabel('Arousal Thres. (%)');
    ylabel('Frequency');
    hold off
end

%get DataN
% for i=1:Nperwindowvars
%     eval([varlist{i} 'N=length(' varlist{i} 'values);']);
% end

if settings.getCIs
    Nbootstrap=50;
    %default
    for ii=1:length(varlist)
        eval([varlist{ii} 'CI=[-Inf;Inf];']);
    end
    for ii=1:Nperwindowvars
        try
            eval([varlist{ii} 'CI=bootci(Nbootstrap,@median,' varlist{ii} 'values);']);
        catch me
        end
    end
else
    Nbootstrap=0;
end
    
    if settings.plotfigs
        figure(1);
        ii=1;
        subplot(2,6,3);
        ylims=get(gca,'YLim');
        hold('on');
        plot(eval(varlist{ii})*[1 1],ylims,'k');
%         try
%         plot(eval([varlist{ii} 'CI(1)'])*[1 1],ylims,'color',[0.5 0.5 0.5]);
%         plot(eval([varlist{ii} 'CI(2)'])*[1 1],ylims,'color',[0.5 0.5 0.5]);
%         end
        
        ii=2;
        subplot(2,6,4);
        ylims=get(gca,'YLim');
        hold('on');
        plot(eval(varlist{ii})*[1 1],ylims,'k');
%         try
%         plot(eval([varlist{ii} 'CI(1)'])*[1 1],ylims,'color',[0.5 0.5 0.5]);
%         plot(eval([varlist{ii} 'CI(2)'])*[1 1],ylims,'color',[0.5 0.5 0.5]);
%         end
        ii=3;
        subplot(2,6,5);
        ylims=get(gca,'YLim');
        hold('on');
        plot(eval(varlist{ii})*[1 1],ylims,'k');
%         try
%         plot(eval([varlist{ii} 'CI(1)'])*[1 1],ylims,'color',[0.5 0.5 0.5]);
%         plot(eval([varlist{ii} 'CI(2)'])*[1 1],ylims,'color',[0.5 0.5 0.5]);
%         end
        ii=4;
        subplot(2,6,6);
        ylims=get(gca,'YLim');
        hold('on');
        plot(eval(varlist{ii})*[1 1],ylims,'k');
%         try
%         plot(eval([varlist{ii} 'CI(1)'])*[1 1],ylims,'color',[0.5 0.5 0.5]);
%         plot(eval([varlist{ii} 'CI(2)'])*[1 1],ylims,'color',[0.5 0.5 0.5]);
%         end
        ii=5;
        subplot(2,6,9);
        ylims=get(gca,'YLim');
        hold('on');
        plot(eval(varlist{ii})*[1 1],ylims,'k');
%         try
%         plot(eval([varlist{ii} 'CI(1)'])*[1 1],ylims,'color',[0.5 0.5 0.5]);
%         plot(eval([varlist{ii} 'CI(2)'])*[1 1],ylims,'color',[0.5 0.5 0.5]);
%         end
    end



%% UA phenotype using model drive

NwindowsForUA=sum(criteriaUA);
str=['Using ', num2str(NwindowsForUA), ' of ', num2str(length(criteriaUA)) ,' windows for UA calcs'];
disp(str);

for ii=(Nperwindowvars+1):length(varlist)
    eval([varlist{ii} 'N=NwindowsForUA;']);
end

%move thiese inside SummaryAnalysisOne_UAmodel:
Fstates3=[mean(WinT.FNREM1(criteriaUA)) mean(WinT.FNREM2(criteriaUA)) mean(WinT.FNREM3(criteriaUA)) mean(WinT.FREM(criteriaUA))]';
Fsupine3 = nanmean(Fsupinetemp(criteriaUA));

criteriaUAI = find(criteriaUA==1);

%make subset table
    criteriaRow = sum((BreathDataTable2.Win == criteriaUAI'),2)>0;
    BreathDataTable3 = BreathDataTable2(criteriaRow==1,:);
    
    %note hypok not used in here:
    try
        [Vpassive,Vactive,VEdeciles,Vdrivedeciles,VpassiveCI,VactiveCI,Vcomp,VcompCI]=SummaryAnalysisOne_UAmodel(length(criteriaUAI),BreathDataTable3,ArThres,varlist,settings.plotfigs,settings);
    catch 
        Vpassive=NaN;
        Vactive=NaN;
        Vcomp=NaN;
        VpassiveCI=NaN;
        VactiveCI=NaN;
        VcompCI=NaN;
        VEdeciles=NaN;
        Vdrivedeciles=NaN;
        disp('failed UA endotype')
    end
%% to do
%% Event probability
if 0                
        e=e(Ikeep);
        e=e(criteriabreath);
        [Epassive,Eactive,Edeciles,Edeciles_drive,EpassiveCI,EactiveCI,EdecilesUpper,EdecilesLower,EdecilesMean]=VEVdriveAnalysis(x*x_,e,ArThresActive*x_,0,x_,settings.Nciles,Nbootstrap,NcileLimits);
        
        if 1 %add colorbar to individual plot
            figure(1);
            subplot(1,3,1);
            DispEventLikelihoodOnEndogram(Vdrivedeciles,VEdeciles,1-EdecilesMean,0)
        end
        
        decileTime=t(Ikeep);
        decileTime=decileTime(criteriabreath);
        [Tpassive,Tactive,Tdeciles,Tdeciles_drive,TpassiveCI,TactiveCI,TdecilesUpper,TdecilesLower,TdecilesMean]=VEVdriveAnalysis(x*x_,decileTime,ArThresActive*x_,0,x_,settings.Nciles,Nbootstrap,NcileLimits);
        
end
%% BootStrapping
if isfield(settings,'Boot') && settings.Boot >0  % settings.Boot has to be the num of iterations required
    settings_ = settings;
    settings_.getCIs=0;
    settings_.Nbootstrap=0;
    
    settings.Boot=50;
    
    I = find(criteriaEither==1);
    for j=1:settings.Boot
        I2 = datasample(I,length(I));
        criteriaUAIBoot = I2;
        criteriaUAIBoot(criteriaUA(I2)==0)=[];
        criteriaLGIBoot = I2;
        criteriaLGIBoot(criteriaLG(I2)==0)=[];
        
    %make subset table
    criteriaRowBoot = sum((BreathDataTable2.Win == criteriaUAIBoot'),2)>0;
    BreathDataTable3Boot = BreathDataTable2(criteriaRowBoot==1,:);
    NwindowsForUABoot=sum(criteriaUAIBoot);
        
        [LG1Boot(j),LGnBoot(j),delayBoot(j),VRABoot(j),ArThresBoot(j),LG1Boot_N(j),LGnBoot_N(j),delayBoot_N(j),VRABoot_N(j),ArThresBoot_N(j),Fstates1Boot(j,:),Fsupine1Boot(j),Fstates2Boot(j,:),Fsupine2Boot(j)]=SummaryAnalysisOne_LGmodel(WinT(criteriaLGIBoot,:),varlist); %,Fstates1Boot(j),Fsupine1Boot(j),Fstates2Boot(j),Fsupine2Boot(j)
        
%         [VpassiveBoot(j),VactiveBoot(j),~,~,~,~,VcompBoot(j)]=SummaryAnalysisOne_UAmodel(DataOut(criteriaUAIBoot),ArThresBoot(j),varlist,0,settings_); %also output Fstates3 Fsupine3
        [VpassiveBoot(j),VactiveBoot(j),~,~,~,~,VcompBoot(j)]=SummaryAnalysisOne_UAmodel(NwindowsForUABoot,BreathDataTable3Boot,ArThresBoot(j),varlist,0,settings_);
        Fstates3Boot(j,:)=[mean(WinT.FNREM1(criteriaUAIBoot)) mean(WinT.FNREM2(criteriaUAIBoot)) mean(WinT.FNREM3(criteriaUAIBoot)) mean(WinT.FREM(criteriaUAIBoot))]';
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
% SS agrees, yep, dodgy as
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

