clear
%% Get folder info
fpath = 'J:\PEOPLE\FACULTY\SANDS\OralApplianceMM2018\Analyzed';
listing = dir([fpath,'\*_localsignals*']);
numsubs = size(listing,1);
listingCell = struct2cell(listing);
names = listingCell(1,:)';

clear listingCell

%load subject info
subjectInfo = readtable('AnalyzeDataSpreadsheet.xlsx','Sheet', 1,'Range','B2:CT34');

%load finale model table
load('FinalModel_Table.mat')

for subnum = 1:numsubs %excluding 31 bc strange error
    %% Load tables
    if ~exist([fpath, '\MMOAT_', num2str(subnum),'.mat'], 'file')
        continue
    end
    
    load([fpath, '\MMOAT_', num2str(subnum),'.mat'])
    load([fpath, '\MMOAT_', num2str(subnum),'_localsignals.mat'])
    
    %% Check that all tables included match
    [BreathDataTable, BreathFLDataTable, LocalSignals] = ...
        checkTables(BreathDataTable, BreathFLDataTable, LocalSignals);
    
    %% Find dupilcate breaths
    %Expand Tables
    BreathDataTableE = expandtable(BreathDataTable);
    BreathFLDataTableE = expandtable(BreathFLDataTable);
    LocalSignalsE = expandtable(LocalSignals);

    FDuplicated = FindDuplicated(BreathDataTable);

    %% Remove Duplicates and Tidy Table (make this a function later)
    cut = 0.67;
    % Breath data table
    BreathDataTableE_ = RemoveDuplicatesTidyTable(BreathDataTableE,...
                            BreathDataTableE.Time_start,...
                            BreathDataTableE.Time_end,...
                            FDuplicated, cut);

    % MIFL breath data table
    BreathFLDataTableE_ = RemoveDuplicatesTidyTable(BreathFLDataTableE,...
        BreathDataTableE.Time_start,...
        BreathDataTableE.Time_end,...
        FDuplicated, cut);

    %Inspiratory flow shapes
    LocalSignalsETable = array2table(LocalSignalsE); %convert to table for function
    LocalSignalsETable_ = RemoveDuplicatesTidyTable(LocalSignalsETable,...
                            BreathDataTableE.Time_start,...
                            BreathDataTableE.Time_end,...
                            FDuplicated, cut);
                        
    LocalSignalsE_ = table2array(LocalSignalsETable_); %convert back to array
    BreathArray = LocalSignalsE_;
    InspArray = LocalSignalsE_(:,1:250);
    ExpArray = LocalSignalsE_(:,251:500);
    
    %Add subject label column
    Subject = ones(size(BreathDataTableE_,1),1)*str2num(subjectInfo.MATFilename{subnum,1}(1:end-8));
    BreathDataTableE_.Subject = Subject;
    
    %% Compute flow drive ratio
    FlowDrive = computeFlowDrive(FinalModelTable, BreathFLDataTableE_);
    BreathDataTableE_.FlowDrive = FlowDrive;
    
    % Remove FL breaths that are not within a cluster of flow limited breaths
    % To qualify, breaths must be part of a group of at least 4
    win = 4;
    consec = 4;
    increment = 1;
    FLidx = false(size(BreathDataTableE_,1),1);
    for ii = 1:increment:size(BreathDataTableE_,1)-win-1
        FLidxTemp = BreathDataTableE_.FlowDrive(ii:ii+win-1) <= 0.7;
        if sum(FLidxTemp) >= consec
            FLidx(ii:ii+win-1,1) = BreathDataTableE_.FlowDrive(ii:ii+win-1) <= 0.7;
        end
    end

    % Isolate breaths with flow drive < 0.7
    BreathDataTableFLTmp = BreathDataTableE_(FLidx,...
                           BreathDataTableE_.Properties.VariableNames);
    BreathFLDataTableFLTmp = BreathFLDataTableE_(FLidx,...
                             BreathFLDataTableE_.Properties.VariableNames);
    BreathArrayFLTmp = BreathArray(FLidx,:);
    InspArrayFLTmp = InspArray(FLidx,:);
    ExpArrayFLTmp = ExpArray(FLidx,:);
    
    % Concatenate to construct large table with only flow limited breaths
    if ~exist('BreathDataTableFL')
        BreathDataTableFL = BreathDataTableFLTmp;
        BreathFLDataTableFL = BreathFLDataTableFLTmp;
        BreathArrayFL = BreathArrayFLTmp;
        InspArrayFL = InspArrayFLTmp;
        ExpArrayFL = ExpArrayFLTmp;
    else
        BreathDataTableFL = vertcat(BreathDataTableFL, BreathDataTableFLTmp);
        BreathFLDataTableFL = vertcat(BreathFLDataTableFL, BreathFLDataTableFLTmp);
        BreathArrayFL = vertcat(BreathArrayFL, BreathArrayFLTmp);
        InspArrayFL = vertcat(InspArrayFL, InspArrayFLTmp);
        ExpArrayFL = vertcat(ExpArrayFL, ExpArrayFLTmp);
    end
    
    %% Isolate Hypopnea Breaths
    %Generate large arrays with only hypopnea breath
    %Add new column indicating hypopnea number
    %Can later filter by hypopnea number

    HypopShiftDown = circshift(BreathDataTableE_.Etype,1);
    HypopShiftUp = circshift(BreathDataTableE_.Etype,-1);
    HypopShiftDown(isnan(HypopShiftDown)) = 0;
    HypopShiftUp(isnan(HypopShiftUp)) = 0;
    HypopDiffDown = BreathDataTableE_.Etype - HypopShiftDown;
    HypopDiffUp = BreathDataTableE_.Etype - HypopShiftUp;

    HypopFirstIdx = find(HypopDiffDown == 4);
    HypopLastIdx = find(HypopDiffUp == 4);

    %Generate variable indicating hypopnea number to be stored in table
    HypopNum = zeros(size(BreathDataTableE_,1),1);
    NumHypop = zeros(size(BreathDataTableE_,1),1);
    %Loop through hypopnea indices and fill hypopnea number variable
    %Note: 1 is the terminal hypopnea, 2 is second last, etc
    for ii = 1:length(HypopFirstIdx)
        totalHypop = HypopLastIdx(ii) - HypopFirstIdx(ii) + 1;
        HypopNum(HypopFirstIdx(ii):HypopLastIdx(ii),1) = totalHypop:-1:1;
        NumHypop(HypopFirstIdx(ii):HypopLastIdx(ii),1) = totalHypop;
    end

    % Add to table
    BreathDataTableE_.HypopNum = HypopNum;
    BreathDataTableE_.NumHypop = NumHypop;

    % Isolate rows in relevant tables with hypopneas
    BreathDataTableHypTmp = BreathDataTableE_(BreathDataTableE_.HypopNum > 0,...
                           BreathDataTableE_.Properties.VariableNames);

    BreathFLDataTableHypTmp =BreathFLDataTableE_(BreathDataTableE_.HypopNum > 0,...
                           BreathFLDataTableE_.Properties.VariableNames);
    
    BreathArrayHypTmp = BreathArray(BreathDataTableE_.HypopNum > 0,:);
    InspArrayHypTmp = InspArray(BreathDataTableE_.HypopNum > 0,:);
    ExpArrayHypTmp = ExpArray(BreathDataTableE_.HypopNum > 0,:);

    clear HypopShiftDown HypopShiftUp HypopDiffDown HypopDiffUp HypopFirstIdx...
          HypopLastIdx ii totalHypop HypopNum NumHypop

    % Concatenate to construct one large table with only hypopneas
    if ~exist('BreathDataTableHyp')
        BreathDataTableHyp = BreathDataTableHypTmp;
        BreathFLDataTableHyp = BreathFLDataTableHypTmp;
        BreathArrayHyp = BreathArrayHypTmp;
        InspArrayHyp = InspArrayHypTmp;
        ExpArrayHyp = ExpArrayHypTmp;
    else
        BreathDataTableHyp = vertcat(BreathDataTableHyp, BreathDataTableHypTmp);
        BreathFLDataTableHyp = vertcat(BreathFLDataTableHyp, BreathFLDataTableHypTmp);
        BreathArrayHyp = vertcat(BreathArrayHyp, BreathArrayHypTmp);
        InspArrayHyp = vertcat(InspArrayHyp, InspArrayHypTmp);
        ExpArrayHyp = vertcat(ExpArrayHyp, ExpArrayHypTmp);
    end
    
    %% Make big table with all data
    if ~exist('BreathDataTableBig')
        BreathDataTableBig = BreathDataTableE_;
    else
        BreathDataTableBig = vertcat(BreathDataTableBig, BreathDataTableE_);
    end
end  

%% Remove time aberrant breaths
% Remove breaths before/after 2 and 98th percentile
Time_BB = BreathDataTableBig.Time_end - BreathDataTableBig.Time_start;
Time_Lower = prctile(Time_BB, 2);
Time_Upper = prctile(Time_BB, 98);

%Hypopnea
BreathDataTableHyp.Time_BB = BreathDataTableHyp.Time_end - ...
                             BreathDataTableHyp.Time_start;
                         
rmIdxHyp = BreathDataTableHyp.Time_BB > Time_Upper | ...
           BreathDataTableHyp.Time_BB < Time_Lower;

%Flow limited
BreathDataTableFL.Time_BB = BreathDataTableFL.Time_end - ...
                             BreathDataTableFL.Time_start; 
                         
rmIdxFL = BreathDataTableFL.Time_BB > Time_Upper | ...
           BreathDataTableFL.Time_BB < Time_Lower;

%% Normalize expiration (added on because expiration data is interesting)
% Could include this in the analysis script
ExpArrayHyp = normalizeExp(ExpArrayHyp);
ExpArrayFL = normalizeExp(ExpArrayFL);

%% Save tables
%Big
save('BreathTableBig.mat', 'BreathDataTableBig')

%Hypopnea
BreathDataTableFinal = BreathDataTableHyp(~rmIdxHyp,:);
BreathFLDataTableFinal = BreathFLDataTableHyp(~rmIdxHyp,:);
BreathArrayFinal = BreathArrayHyp(~rmIdxHyp,:);
InspArrayFinal = InspArrayHyp(~rmIdxHyp,:);
ExpArrayFinal = ExpArrayHyp(~rmIdxHyp,:);
save('HypopTables.mat', 'BreathDataTableFinal', 'BreathFLDataTableFinal',...
    'BreathArrayFinal', 'InspArrayFinal', 'ExpArrayFinal')

%Flow limited
BreathDataTableFinal = BreathDataTableFL(~rmIdxFL,:);
BreathFLDataTableFinal = BreathFLDataTableFL(~rmIdxFL,:);
BreathArrayFinal = BreathArrayFL(~rmIdxFL,:);
InspArrayFinal = InspArrayFL(~rmIdxFL,:);
ExpArrayFinal = ExpArrayFL(~rmIdxFL,:);
save('FLTables', 'BreathDataTableFinal', 'BreathFLDataTableFinal',...
    'BreathArrayFinal', 'InspArrayFinal', 'ExpArrayFinal')