function SummaryAnalysisTable = SummaryAnalysisN(MrangeOverride)

clc
%% load globals
global SAfontsize AMasterSpreadsheet handletext settings


%% default plotting
SAfontsize=12;
set(groot,'defaultAxesfontname','arial narrow');
set(groot,'defaultAxesfontsize',SAfontsize);
set(groot,'defaultFigureColor',[1 1 1]);
set(groot,'defaultAxesBox','off');

%% start processing
displaytext='Starting Summary Analysis';
disp(displaytext); set(handletext,'String',displaytext); drawnow;

%% Import settings
if ~isfield(settings,'ImportedSettingsComplete')
    settings = ImportSettings(settings,AMasterSpreadsheet);
end
if isfield(settings,'SummaryMrange')
    Mrange = settings.SummaryMrange;
end
if isfield(settings,'Mrange')
    %     if ischar(settings.Mrange)
    %         settings.Mrange = eval([settings.Mrange]);
    %     else
    Mrange = settings.Mrange;
    %     end
end
if exist('MrangeOverride')==1
    Mrange = MrangeOverride;
end
if ~exist('Mrange')==1
    disp('Mrange not set, who do you want to Summarize?')
    return
end

% settings.Mrange is used below, and Mrange is also used below, should this just be one thing?

% dlm found that if settings.csv does not exist, then you can get here with
% only Mrange set, but not settings.Mrange set, and this causes crash
% so I'm just testing for that scenario, then manually setting settings.Mrange to keep things happy below
% not sure if this is going to break things elsewhere
if exist('Mrange') && ~isfield(settings,'Mrange')
    settings.Mrange = Mrange;
end


%% Ludo Overwrites Spreadsheet
if isfield(settings,'compareREMvNREM') && settings.compareREMvNREM==1
    settings.Mrange = [2 10 11 16 19 20 25 26 27 30 34 37 40 43 44 46 48 49 50 51 53 55 58 59 61];
    %    settings.Mrange = [41];
    %    settings.Mrange = [10];
    %35,38 - No Edi
end

%% Detect if Mrange is char, if yes -> make array -- edit Sara
if isa(settings.Mrange,'char')
    settings.Mrange = str2num(settings.Mrange);
end

%% Run
cutoffvalues=[0.7 0.4 12 40 115 85 95 10];
posclass=[1 1 1 1 0 0 0 0];
units={'','','(s)','(%eupnea)','(%eupnea)','(%eupnea)','(%eupnea)','(%eupnea)'};
longname={'Loop gain, sensitivity','Loop gain, instability','Chemoreflex delay','Vresponse-to-arousal','Arousal threshold',...
    'Collapsibility, Vpassive','Collapsibility, Vactive','Compensation'};

clear Data DataCI DataN AHItotal success Fstates Veupnea VEdeciles Vdrivedeciles medianV GGdata GGpdeciles GGtdeciles
clear Data_compare Veupnea_compare VEdeciles_compare Vdrivedeciles_compare medianV_compare GGdata_compare GGpdeciles_compare GGtdeciles_compare
clear Data_compare2 Veupnea_compare2 VEdeciles_compare2 Vdrivedeciles_compare2 medianV_compare2 GGdata_compare2 GGpdeciles_compare2 GGtdeciles_compare2

%initializing
Ntraitsout=8 + (settings.method==2);
Data = nan(max(Mrange),Ntraitsout);
Data_compare = nan(max(Mrange),Ntraitsout);
Data_compare2 = nan(max(Mrange),Ntraitsout);
DataN = nan(max(Mrange),Ntraitsout);
DataN_compare = nan(max(Mrange),Ntraitsout);
DataN_compare2 = nan(max(Mrange),Ntraitsout);

GGdata = nan(max(Mrange),4);
GGdata_compare = nan(max(Mrange),4);
GGdata_compare2 = nan(max(Mrange),4);

VEdeciles = nan(settings.Nciles,max(Mrange));
Vdrivedeciles = nan(settings.Nciles,max(Mrange));
VEdeciles_compare = nan(settings.Nciles,max(Mrange));
Vdrivedeciles_compare = nan(settings.Nciles,max(Mrange));
VEdeciles_compare2 = nan(settings.Nciles,max(Mrange));
Vdrivedeciles_compare2 = nan(settings.Nciles,max(Mrange));
GGpdeciles = nan(settings.Nciles,max(Mrange));
GGpdeciles_compare = nan(settings.Nciles,max(Mrange));
GGpdeciles_compare2 = nan(settings.Nciles,max(Mrange));
GGtdeciles = nan(settings.Nciles,max(Mrange));
GGtdeciles_compare = nan(settings.Nciles,max(Mrange));
GGtdeciles_compare2 = nan(settings.Nciles,max(Mrange));

EdecilesMean=nan(settings.Nciles,max(Mrange));
EdecilesMean_compare=nan(settings.Nciles,max(Mrange));
EdecilesMean_compare2=nan(settings.Nciles,max(Mrange));

TdecilesMean=nan(settings.Nciles,max(Mrange));
TdecilesMean_compare=nan(settings.Nciles,max(Mrange));
TdecilesMean_compare2=nan(settings.Nciles,max(Mrange));

medianV = nan(max(Mrange),1);
medianV_compare = nan(max(Mrange),1);
medianV_compare2 = nan(max(Mrange),1);

Veupnea = nan(max(Mrange),1);
Veupnea_compare = nan(max(Mrange),1);
Veupnea_compare2 = nan(max(Mrange),1);

AHIdata = nan(max(Mrange),184);
%AHIdata_compare = nan(max(Mrange),184);

Fsupine = nan(max(Mrange),Ntraitsout);
Fsupine_compare = nan(max(Mrange),Ntraitsout);
Fsupine_compare2 = nan(max(Mrange),Ntraitsout);

AHItotal = nan(max(Mrange),1);

SummaryLargeT = table();
SummaryLargeT.FNoisyLG = nan(max(Mrange),1);
SummaryLargeT.FnEventsLG = nan(max(Mrange),1); 
SummaryLargeT.FFailOtherLG = nan(max(Mrange),1); 
SummaryLargeT.FAnalyzableLG = nan(max(Mrange),1);
SummaryLargeT.FNoisyUA = nan(max(Mrange),1);
SummaryLargeT.FnEventsUA = nan(max(Mrange),1);
SummaryLargeT.FFailOtherUA = nan(max(Mrange),1); 
SummaryLargeT.FAnalyzableUA = nan(max(Mrange),1);

%Fsupine=NaN;

clear BootT

for i=Mrange%1:M
    disp(['Subj: ' num2str(i)])
    try
        figure(4); clf(4);
        figure(2); clf(2);
        figure(1); clf(1);
        disp(settings.method);
        if settings.method==1 %Flow only, model based phenotyping
            %[Data{i},DataCI{i},DataN{i},varlist,AHItotal(i),Fstates{i},Veupnea(i,1),medianV(i,1),~,~,AHIdata(i,1)] = SummaryAnalysisOne(i,settings);
            if settings.runcompare
                settings2 = settings;
                %update with settings, e.g. predefine settings.settingsCompare.selectstate=4; or vary by some other setting e.g. position or approach.
                for currField = fieldnames(settings.settingsCompare)'
                    settings2 = setfield(settings2,currField{:},getfield(settings.settingsCompare,currField{:}));
                end
                [Data_compare(i,:),DataCI_compare{i},DataN_compare(i,:),varlist,AHItotal(i,1),Fstates_compare{i},...
                    Veupnea_compare(i,1),medianV_compare(i,1),VEdeciles_compare(:,i),Vdrivedeciles_compare(:,i),~,Fsupine_compare(i,:),WinTcell_compare{i},~,~,SummaryLargeT(i,:)] = ...
                    SummaryAnalysisOne(i,settings2);                
            end
            
            if isfield(settings,'PrintSummaryAnalysisSetup') && settings.PrintSummaryAnalysisSetup==1 % DLM: debugging the changing selectstate variable
                try disp(['selectstate: ', num2str(settings.selectstate)]); catch; end
                try disp(['selectstateEventAnalysis: ', num2str(settings.selectstateEventAnalysis)]); catch; end
                try disp(['selectposition: ', num2str(settings.selectposition)]); catch; end
                try disp(['analysisRange: ', num2str(settings.analysisRange)]); catch; end
                %try disp(['Mrange: ', num2str(Mrange)]); catch; end
            end
            
            [Data(i,:),DataCI{i},DataN(i,:),varlist,AHItotal(i,1),Fstates{i},...
                Veupnea(i,1),medianV(i,1),VEdeciles(:,i),Vdrivedeciles(:,i),AHIdata(i,:),Fsupine(i,:),WinTcell{i},BootTemp,EvtsTemp,SummaryLargeT(i,:)] = ...
                SummaryAnalysisOne(i,settings);
            
            if isfield(settings, 'Boot') && settings.Boot>0
                if ~isempty(BootTemp)
                    if ~exist('BootT')
                        BootT = BootTemp;
                    else
                        BootT = [BootT;BootTemp];
                    end
                end
            end
            
            
        elseif settings.method==2 %Invasive drive/effort-based phenotyping
            %settings.scoredarousalsinwake=1
            if settings.runcompare
                settings2 = settings;
                %update with settings, e.g. predefine settings.settingsCompare.selectstate=4; or vary by some other setting e.g. position or approach.
                for currField = fieldnames(settings.settingsCompare)'
                    settings2 = setfield(settings2,currField{:},getfield(settings.settingsCompare,currField{:}));
                end
                [Data_compare(i,:),DataCI_compare{i},DataN_compare(i,:),varlist,AHItotal(i,1),Fstates_compare{i},...
                    Veupnea_compare(i,1),medianV_compare(i,1),VEdeciles_compare(:,i),Vdrivedeciles_compare(:,i),~,...
                    Fsupine_compare(i,:),DriveEupnea_compare(i,1),~,~,GGdata_compare(i,:),GGpdeciles_compare(:,i),...
                    GGtdeciles_compare(:,i),EdecilesMean_compare(:,i),TdecilesMean_compare(:,i)] = SummaryAnalysisOnePes(i,settings2);
                
                if settings.runcompare2
                    settings3=settings;
                    for currField = fieldnames(settings.settingsCompare2)'
                        settings3 = setfield(settings2,currField{:},getfield(settings.settingsCompare2,currField{:}));
                    end
                    [Data_compare2(i,:),DataCI_compare2{i},DataN_compare2(i,:),varlist,AHItotal(i,1),Fstates_compare2{i},...
                        Veupnea_compare2(i,1),medianV_compare2(i,1),VEdeciles_compare2(:,i),Vdrivedeciles_compare2(:,i),~,...
                        Fsupine_compare2(i,:),DriveEupnea_compare2(i,1),~,~,GGdata_compare2(i,:),GGpdeciles_compare2(:,i),...
                        GGtdeciles_compare2(:,i),EdecilesMean_compare2(:,i)] = SummaryAnalysisOnePes(i,settings3);
                end
                
            end
            % settings.PhasicREMonly=0;
            % settings.TonicREMonly=1;
            [Data(i,:),DataCI{i},DataN(i,:),varlist,AHItotal(i,1),Fstates{i},...
                Veupnea(i,1),medianV(i,1),VEdeciles(:,i),Vdrivedeciles(:,i),AHIdata(i,:),Fsupine(i,:),DriveEupnea(i,1),~,~,...
                GGdata(i,:),GGpdeciles(:,i),GGtdeciles(:,i),EdecilesMean(:,i),TdecilesMean(:,i)] = ...
                SummaryAnalysisOnePes(i,settings);
            
            if settings.runcompare
                TempDisplay = [Data_compare(i,:);Data(i,:)];
                if isnan(Data_compare(i,1)) | isnan(Data(i,1))
                    compareSuccess=0
                else
                    compareSuccess=1
                end
            end
            try
                run TidyOverlaidSummaryPhenoPlots %should work ok without runcompare
            catch
                'failed TidyOverlaidSummaryPhenoPlots';
            end
        end
        [AllPositionsEventsTable{i},AllPositionsEventsTableN{i},...
            SupineEventsTable{i},EvtDurTable{i}]=...
            AHIdata2Tbls(AHIdata(i,:));
        success(i)=1;
    catch me
        disp(me.message); % getReport(me)
        success(i)=0;
    end
    
    % settings.saveplots % these options are from analysis part of spreadhseet
    % settings.plotfigure
    
    try
        run TidyOverlaidSummaryPhenoPlots %working for SummaryAnalysisOne (OnePes is untested)
    catch
    end
    
    if settings.plotfigs   % plotfigs is an option in summary part of spreadsheet
        %try saveas(1,[settings.SummaryDirectory 'Fig1_' num2str(i) '_s' num2str(settings.selectstate) '_p' num2str(settings.selecttimeofnight_NthXTile) '.png']); end
        try saveas(1,[settings.SummaryDirectory 'Fig1_' num2str(i) '_s' num2str(settings.selectstate) '_p' num2str(settings.selecttimeofnight_NthXTile) '_pos' num2str(settings.selectposition) '.png']); catch end
        
        %try saveas(2,[settings.SummaryDirectory 'Fig2_' num2str(i) '_s' num2str(settings.selectstate) '_p' num2str(settings.selecttimeofnight_NthXTile) '.png']); end
        %try saveas(4,[settings.SummaryDirectory 'Fig4_' num2str(i) '_s' num2str(settings.selectstate) '_p' num2str(settings.selecttimeofnight_NthXTile) '.png']); end
    end
    
    if settings.runcompare
        %     if compareSuccess
        try
            saveas(1,[settings.SummaryDirectory 'Fig1_' settings.Comparename '_' num2str(i) '_s' num2str(settings.selectstate) '_p' num2str(settings.selecttimeofnight_NthXTile) '.png'])
            %    saveas(2,[settings.SummaryDirectory 'Fig2_' settings.Comparename '_' num2str(i) '_s' num2str(settings.selectstate) '_p' num2str(settings.selecttimeofnight_NthXTile) '.png'])
            %    saveas(4,[settings.SummaryDirectory 'Fig4_' settings.Comparename '_' num2str(i) '_s' num2str(settings.selectstate) '_p' num2str(settings.selecttimeofnight_NthXTile) '.png'])
        catch
            disp('failed saving figures')
        end
        %  end
    end
end
I=find(success==1);
Ifail=find(success==0);

%%
if ~exist('varlist')
    disp('Warning: Summary Analysis Failed');
    return
end

%% Find missing arousal scoring
% OLD code
% ArousalDataMissing=zeros(max(Mrange),1);
% for i=Mrange
%     if success(i)==1
%         if DataN{i}(1,4)==0&DataN{i}(1,5)==0&(DataN{i}(1,1)>3)
%             ArousalDataMissing(i)=1;
%         end
%     end
% end

ArousalDataMissing = success(:)==1 & DataN(:,4)==0 & DataN(:,5)==0 & DataN(:,1)>3;
IArousalDataMissing=find(ArousalDataMissing==1);

%% Overwrite, n.b. changes index
if 0
    Data=Data(I);
    DataCI=DataCI(I);
    AHItotal=AHItotal(I);
    DB=DB(I,:);
end

%%
clear DataArray
rowofNaNs = NaN*ones(1,length(varlist));

% for i=1:size(Data,2)
%     DataArray(i,:) = rowofNaNs; %default
%     try
%         DataArray(i,:) = Data{i};
%     catch me
%         DataArray(i,:) = NaN*DataArray(1,:);
%     end
%     DataArray(i,DataN{i}<settings.minwindows)=NaN;
% end
DataArray=Data;
DataArray(DataN<settings.minwindows)=NaN;

if settings.nanifnoarousals
    DataArray(IArousalDataMissing,:)=NaN;
end

if settings.runcompare
    DataArray_compare=Data_compare;
    DataArray_compare(DataN_compare<settings.minwindows)=NaN;
    
    if settings.runcompare2
        DataArray_compare2=Data_compare2;
        DataArray_compare2(DataN_compare2<settings.minwindows)=NaN;
    end
end


%%
% clear DataArrayN
% for i=1:size(Data,2)
%     DataNArray(i,:) = rowofNaNs; %default
%     try
%         DataNArray(i,:) = DataN{i};
%     catch me
%         DataNArray(i,:) = zeros(1,length(DataNArray(1,:)));
%     end
%     %DataNArray(i,DataN{i}<3)=NaN;
% end
DataNArray = DataN;
Temp = zeros(size(DataNArray,1),1); %assume failed analysis had N=0 windows, but analysis not chosen has N=NaN;
Temp(Mrange)=1;
DataNArray(isnan(DataNArray) & Temp==1)=0;
if settings.nanifnoarousals
    DataNArray(IArousalDataMissing,:)=NaN;
end

if settings.runcompare
    DataNArray_compare = DataN_compare;
    Temp = zeros(size(DataNArray_compare,1),1); %assume failed analysis had N=0 windows, but analysis not chosen has N=NaN;
    Temp(Mrange)=1;
    DataNArray_compare(isnan(DataNArray_compare) & Temp==1)=0;
    if settings.nanifnoarousals
        DataNArray_compare(IArousalDataMissing,:)=NaN;
    end
    if settings.runcompare2
        DataNArray_compare2 = DataN_compare2;
        Temp = zeros(size(DataNArray_compare2,1),1); %assume failed analysis had N=0 windows, but analysis not chosen has N=NaN;
        Temp(Mrange)=1;
        DataNArray_compare2(isnan(DataNArray_compare2) & Temp==1)=0;
        if settings.nanifnoarousals
            DataNArray_compare2(IArousalDataMissing,:)=NaN;
        end
    end
end


%%
clear Upper Lower
for i=1:size(Data,1)
    try
        Upper(i,:) = DataCI{i}(2,:);
    catch me
        Upper(i,:) = Inf*DataArray(1,:);
    end
end
for i=1:size(Data,1)
    try
        Lower(i,:) = DataCI{i}(1,:);
    catch me
        Lower(i,:) = -Inf*DataArray(1,:);
    end
end
CIwidth=(Upper-Lower)/2;
%%
clear FstatesArray
for i=1:length(Data)
    for j=1:length(varlist)
        try
            FstatesArray{j}(i,:) = Fstates{i}(:,j)';
        catch me
            FstatesArray{j}(i,:) = [NaN NaN NaN NaN];
        end
    end
end

%%

% clear FsupineArray
% for i=1:length(Data)
%     for j=1:length(varlist)
%         try
%             FsupineArray{j}(i) = Fsupine(i,j);
%         catch me
%             FsupineArray{j}(i) = [NaN NaN NaN NaN];
%         end
%     end
% end
FsupineArray = Fsupine;

%%
LargerData = [DataArray medianV Veupnea];

%% Make Table
SummaryAnalysisTable = array2table([DataArray]);
SummaryAnalysisTable.Properties.VariableNames = varlist;

%check which way up VEdeciles
%vpassive1d
try
    SummaryAnalysisTable.Vmin = VEdeciles(1,1:size(DataArray,1))';
catch
end

% add regular transformation
try
    SummaryAnalysisTable.VpassiveT = fVpassiveT(SummaryAnalysisTable.Vpassive);
    SummaryAnalysisTable.ArThresT = fArThresT(SummaryAnalysisTable.ArThres);
catch
end

% add "A" Transformation
try
    SummaryAnalysisTable.VpassiveTA= fVpassiveTmesa(SummaryAnalysisTable.Vpassive);
    SummaryAnalysisTable.VactiveTA = fVpassiveTmesa(SummaryAnalysisTable.Vactive);
catch
end

SummaryAnalysisTable
SummaryAnalysisTableN = array2table([DataNArray]);
SummaryAnalysisTableN.Properties.VariableNames = varlist;

%% Make Compare Table
if exist('DataArray_compare') && settings.runcompare==1
    SummaryAnalysisTable_compare = array2table([DataArray_compare]);
    SummaryAnalysisTable_compare.Properties.VariableNames = varlist;
    
    %check which way up VEdeciles
    %vpassive1d
    try
        SummaryAnalysisTable_compare.Vmin = VEdeciles_compare(1,1:size(DataArray_compare,1))';
    catch
    end
    
    % add regular transformation
    try
        SummaryAnalysisTable_compare.VpassiveT = fVpassiveT(SummaryAnalysisTable_compare.Vpassive);
        SummaryAnalysisTable_compare.ArThresT = fArThresT(SummaryAnalysisTable_compare.ArThres);
    catch
    end
    
    % add "A" Transformation
    try
        SummaryAnalysisTable_compare.VpassiveTA = fVpassiveTmesa(SummaryAnalysisTable_compare.Vpassive);
        SummaryAnalysisTable_compare.VactiveTA = fVpassiveTmesa(SummaryAnalysisTable_compare.Vactive);
    catch
    end
    
    SummaryAnalysisTable_compare
    SummaryAnalysisTableN_compare = array2table([DataNArray_compare]);
    SummaryAnalysisTableN_compare.Properties.VariableNames = varlist;
end

%% This option overwrites LG1, LG2 and delay with LG0, Tau and Tn
% this is so we can calculate LG at any frequency
% this is all done in SummaryAnalysisOne, but process is interrupted here
% so that we can do post-processing at this point.
if isfield(settings,'coreLGvals') && settings.coreLGvals==1
    keyboard
    SummaryAnalysisTable
    
end

%%
if isfield(settings,'compareREMvNREM') && settings.compareREMvNREM
    save REMvNREMSummaryWorkspace
    %     try
    %     run statsREMvNREM
    %     catch me
    %     end
end

%% Assign in base
if settings.comb(1)
    assignin('base',strcat('SummaryAnalysisTable','_',num2str(settings.comb(2))),SummaryAnalysisTable);
    assignin('base',strcat('SummaryAnalysisTableN','_',num2str(settings.comb(2))),SummaryAnalysisTableN);
    assignin('base',strcat('FstatesArray','_',num2str(settings.comb(2))),FstatesArray);
    assignin('base',strcat('Fsupine','_',num2str(settings.comb(2))),Fsupine);
    assignin('base','AllPositionsEventsTable',AllPositionsEventsTable);
    assignin('base','AllPositionsEventsTableN',AllPositionsEventsTableN);
    assignin('base','SupineEventsTable',SupineEventsTable);
    assignin('base','EvtDurTable',EvtDurTable);
    assignin('base','WinTcell',WinTcell);
    if isfield(settings, 'Boot') && settings.Boot>0
        if exist('BootT')
            assignin('base','BootT',BootT);
        end
    end
else
    assignin('base','SummaryAnalysisTable',SummaryAnalysisTable);
    assignin('base','SummaryAnalysisTableN',SummaryAnalysisTableN);
    assignin('base','VEdeciles',VEdeciles);
    assignin('base','Vdrivedeciles',Vdrivedeciles); %1st data point is 0-5th ctiles; 6th is 0-10th ctiles; 16th is 0-20th ctiles
    assignin('base','FstatesArray',FstatesArray);
    assignin('base','Fsupine',Fsupine);
    assignin('base','AllPositionsEventsTable',AllPositionsEventsTable);
    assignin('base','AllPositionsEventsTableN',AllPositionsEventsTableN);
    assignin('base','SupineEventsTable',SupineEventsTable);
    assignin('base','EvtDurTable',EvtDurTable);
    
    try assignin('base','WinTcell',WindowTable); catch; end
    if isfield(settings, 'Boot') && settings.Boot>0
        if exist('BootT')
            assignin('base','BootT',BootT);
        end
    end
    
end

%% Tidy
clear i I IArousalDataMissing Ifail j posclass success

%% Saving
displaytext='Saving Summary Analysis';
disp(displaytext); set(handletext,'String',displaytext); drawnow;
% default saving
if ~(isfield(settings,'NSRRSave') &&  settings.NSRRSave==1)    
    %savefilename = [settings.SummaryDirectory 'SummaryAnalysis_' '_s' num2str(settings.selectstate) '_p' num2str(settings.selecttimeofnight_NthXTile) '_' datestr(now,'yyyymmdd HHMMSS')];
    savefilename = [settings.SummaryDirectory 'SummaryAnalysis_' '_s' num2str(settings.selectstate) '_p' num2str(settings.selecttimeofnight_NthXTile) '_pos' num2str(settings.selectposition) '_' datestr(now,'yyyymmdd HHMMSS')];
    savefilename(end-6)='T';
    save(savefilename)
    % DLM commented out the next two lines, seems redundant
    %savefilename = [settings.SummaryDirectory 'SummaryAnalysis'];
    %save(savefilename)
    
    % NSRR/A to Z mode saving
elseif isfield(settings,'NSRRSave') &&  settings.NSRRSave==1
    
    Outfoldername=settings.napfolder;
    if ~(exist(Outfoldername, 'dir') == 7)
        mkdir(Outfoldername);
    end
    
    
    ID = settings.nsrrid;
    
    if settings.selectstate==4
        TxtFname=[Outfoldername 'resp_endotype_summary_SS_NREM.txt']
        MatFname = [settings.SummaryDirectory ID '_resp_endotype_summary_SS_NREM_' datestr(now,'yyyymmdd HHMMSS')]
        
    elseif settings.selectstate==5
        TxtFname=[Outfoldername 'resp_endotype_summary_SS_REM.txt']
        MatFname = [settings.SummaryDirectory ID '_resp_endotype_summary_SS_REM_' datestr(now,'yyyymmdd HHMMSS')]
        
        
    elseif settings.selectstate==8
        TxtFname=[Outfoldername 'resp_endotype_summary_SS_ALL.txt']
        MatFname = [settings.SummaryDirectory ID '_resp_endotype_summary_SS_ALL_' datestr(now,'yyyymmdd HHMMSS')]
    end
    
    MatFname(end-6)='T';
    save(MatFname); % saving .mat file
    SummaryAnalysisTableTemp =  round(SummaryAnalysisTable.Variables*100)/100;
    SummaryAnalysisTableTemp=array2table(SummaryAnalysisTableTemp);
    SummaryAnalysisTableTemp.Properties.VariableNames=upper(SummaryAnalysisTable.Properties.VariableNames);
    SummaryAnalysisTableTemp.ID=ID;
    SummaryAnalysisTableTemp=movevars(SummaryAnalysisTableTemp,'ID','Before','LG1');
    writetable(SummaryAnalysisTableTemp, TxtFname,'WriteVariableNames',true,'Delimiter','\t');
end

%%
%% Complete
displaytext='Completed Summary Analysis';
disp(displaytext); set(handletext,'String',displaytext); drawnow;
