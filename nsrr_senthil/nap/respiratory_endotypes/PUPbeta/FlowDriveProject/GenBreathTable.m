function GenBreathTable(options)

% Set options
ProjectFolder = options.ProjectFolder;
fpath = options.fpath;
prefix = options.prefix;
PUPversion = options.PUPversion;
SubID = options.SubID;
SubID2 = options.SubID2;
SubIndex = options.SubIndex;
filepath2 = options.filepath2;
processExclude = options.processExclude;

%load finale model table
load([fpath, PUPversion, 'FinalModel_Table.mat'])

for subnum = 1:length(SubIndex)     
    %% Load tables
    if sum(subnum == processExclude) > 0
        disp(['Skipping subject: ', SubID{subnum}, ', Index: ', num2str(SubIndex(subnum))])
        if exist([fpath,'Analyzed\',prefix,'_',num2str(SubIndex(subnum)),'.mat'], 'file')
            disp(['Warning: Skipped subject ', SubID{subnum},' was analyzed'])
            continue
        end
        continue
    end
    
    disp(['Analyzing subject: ', SubID{subnum}, ', Index: ', num2str(SubIndex(subnum))])
    load([fpath,'Analyzed\',prefix,'_',num2str(SubIndex(subnum)),'.mat'])
    
    %% Check that all tables included match
    [BreathDataTable, BreathFLDataTable] = ...
        checkTables(BreathDataTable, BreathFLDataTable);
    
    %% Find dupilcate breaths
    %Expand Tables
    BreathDataTableE = expandtable(BreathDataTable);
    BreathFLDataTableE = expandtable(BreathFLDataTable);

    FDuplicated = FindDuplicated(BreathDataTable);

    %% Remove Duplicates and Tidy Table (make this a function later)
    cut = 0.67;
    
    BreathDataTableE_ = BreathDataTableE;
    BreathDataTableE_(FDuplicated>cut,:) = [];
    
    BreathFLDataTableE_ = BreathFLDataTableE;
    BreathFLDataTableE_(FDuplicated>cut,:) = [];
    
    %% Find duplicate breaths part 2
    FDuplicated2 = FindDuplicatedPt2(BreathDataTableE_);

    BreathDataTableE_(FDuplicated2==1,:) = []; 
    BreathFLDataTableE_(FDuplicated2==1,:) = [];
        
    %Add subject label column
    Subject = repmat({[SubID2{subnum}]}, size(BreathDataTableE_,1), 1); % Made exception for Sara bc of naming inconsistency
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
    
    % Concatenate to construct large table with only flow limited breaths
    if ~exist('BreathDataTableFL')
        BreathDataTableFL = BreathDataTableFLTmp;
        BreathFLDataTableFL = BreathFLDataTableFLTmp;
    else
        BreathDataTableFL = vertcat(BreathDataTableFL, BreathDataTableFLTmp);
        BreathFLDataTableFL = vertcat(BreathFLDataTableFL, BreathFLDataTableFLTmp);
    end
    
    %% Isolate Hypopnea Breaths
    %Generate large arrays with only hypopnea breath
    %Add new column indicating hypopnea number
    %Can later filter by hypopnea number
    Etype4 = BreathDataTableE_.Etype;
    Etype4(Etype4 ~= 4) = 0;
    
    HypopShiftDown = circshift(Etype4,1);
    HypopShiftUp = circshift(Etype4,-1);
    HypopDiffDown = Etype4 - HypopShiftDown;
    HypopDiffUp = Etype4 - HypopShiftUp;

    HypopFirstIdx = find(HypopDiffDown == 4);
    HypopLastIdx = find(HypopDiffUp == 4);

    %Generate variable indicating hypopnea number to be stored in table
    HypopNum = zeros(size(Etype4,1),1);
    NumHypop = zeros(size(Etype4,1),1);
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

    clear HypopShiftDown HypopShiftUp HypopDiffDown HypopDiffUp HypopFirstIdx...
          HypopLastIdx ii totalHypop HypopNum NumHypop

    % Concatenate to construct one large table with only hypopneas
    if ~exist('BreathDataTableHyp')
        BreathDataTableHyp = BreathDataTableHypTmp;
        BreathFLDataTableHyp = BreathFLDataTableHypTmp;
    else
        BreathDataTableHyp = vertcat(BreathDataTableHyp, BreathDataTableHypTmp);
        BreathFLDataTableHyp = vertcat(BreathFLDataTableHyp, BreathFLDataTableHypTmp);
    end
    
    %% Isolate Apnea breaths
    Etype2 = BreathDataTableE_.Etype;
    Etype2(Etype2 ~= 2) = 0;
    
    ApneaShiftDown = circshift(Etype2,1);
    ApneaShiftUp = circshift(Etype2,-1);
    ApneaDiffDown = Etype2 - ApneaShiftDown;
    ApneaDiffUp = Etype2 - ApneaShiftUp;

    ApneaFirstIdx = find(ApneaDiffDown == 2);
    ApneaLastIdx = find(ApneaDiffUp == 2);

    %Generate variable indicating apnea number to be stored in table
    ApneaNum = zeros(size(Etype2,1),1);
    NumApnea = zeros(size(Etype2,1),1);
    %Loop through apnea indices and fill apnea number variable
    %Note: 1 is the terminal hypopnea, 2 is second last, etc
    for ii = 1:length(ApneaFirstIdx)
        totalApnea = ApneaLastIdx(ii) - ApneaFirstIdx(ii) + 1;
        ApneaNum(ApneaFirstIdx(ii):ApneaLastIdx(ii),1) = totalApnea:-1:1;
        NumApnea(ApneaFirstIdx(ii):ApneaLastIdx(ii),1) = totalApnea;
    end

    % Add to table
    BreathDataTableE_.ApneaNum = ApneaNum;
    BreathDataTableE_.NumApnea = NumApnea;
    
    %% Make big table with all data
    if ~exist('BreathDataTableBig')
        BreathDataTableBig = BreathDataTableE_;
        BreathFLDataTableBig = BreathFLDataTableE_;
    else
        BreathDataTableBig = vertcat(BreathDataTableBig, BreathDataTableE_);
        BreathFLDataTableBig = vertcat(BreathFLDataTableBig, BreathFLDataTableE_);
    end
end  

%% Remove time aberrant breaths
% Remove breaths before/after 2 and 98th percentile
if 0
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
 
       
    
%Big FL
BreathDataTableBig.Time_BB = BreathDataTableBig.Time_end - ...
                             BreathDataTableBig.Time_start; 
                         
rmIdxBig = BreathDataTableBig.Time_BB > Time_Upper | ...
           BreathDataTableBig.Time_BB < Time_Lower;
else
    rmIdxHyp = zeros(size(BreathDataTableHyp,1),1);
    rmIdxFL = zeros(size(BreathDataTableFL,1),1);
    rmIdxBig = zeros(size(BreathDataTableBig,1),1);
end

%% Save tables
%Big
BreathDataTableBig = BreathDataTableBig(~rmIdxBig,:);
BreathFLDataTableBig = BreathFLDataTableBig(~rmIdxBig,:);
save([filepath2, ProjectFolder,'\BreathTableBig.mat'],...
    'BreathDataTableBig', 'BreathFLDataTableBig')

%Hypopnea
BreathDataTableFinal = BreathDataTableHyp(~rmIdxHyp,:);
BreathFLDataTableFinal = BreathFLDataTableHyp(~rmIdxHyp,:);
save([filepath2, ProjectFolder,'\HypopTables.mat'], 'BreathDataTableFinal', 'BreathFLDataTableFinal')

%Flow limited
BreathDataTableFinal = BreathDataTableFL(~rmIdxFL,:);
BreathFLDataTableFinal = BreathFLDataTableFL(~rmIdxFL,:);
save([filepath2, ProjectFolder, '\FLTables.mat'], 'BreathDataTableFinal', 'BreathFLDataTableFinal')