function [VIstruct,FeatureTable,BreathTables,AHItable] = ProcessBreathDataTables()

%% global variables and settings
global AMasterSpreadsheet handletext settings n 

%% start processing
% read spreadsheet (files worksheet)
subjectID = settings.MasterWorksheet(:,1);
patients = settings.patients;
[~,~,subjectID2] = xlsread(AMasterSpreadsheet,1,'T4:T10003');

if sum(strcmp(subjectID, '')) > 0 || isempty(subjectID)
    subjectID2_ = subjectID2(settings.analyzelist);
    subjectID = cellfun(@num2str,subjectID2_,'un',0);
end

if strcmp(settings.savename,'SaraODB_OAT')
    [~,subjectID,~] = xlsread(AMasterSpreadsheet,1,'R4:R10003');
end

% Remove dashes from subject ID
subjectIDorig = subjectID;
subjectID = regexprep(subjectIDorig,'-','');

% load('RICCADSA CPAP table.mat')
% Load sample breath data tables
tempi = find(settings.analyzelist,1,'last');
loadpath=[settings.AnalyzedDirectory, 'EventAnalyzed\' settings.savename '_' num2str(tempi)];
load(loadpath,'BreathFLDataTableLong','EAinfo','Evts','EvtsUse'); % Open to get variable names

% CPAPtable(end,:) = [];
FLendIdx = 123;
Varnames = [BreathFLDataTableLong.Properties.VariableNames(1:FLendIdx),...
    fieldnames(EAinfo.EvtFtrs.EventDepth)',fieldnames(EAinfo.EvtFtrs.DesatSlope)',...
    'PRDeltaFromMean','PRDeltaFromMin','PRDeltaFromMeanAr','PRDeltaFromMinAr',...
    'AvgEvDepth','MinEvDepth','EvArea'];
MeanFeatArray = nan(size(patients,1),length(Varnames));
MedFeatArray = nan(size(patients,1),length(Varnames));
SubjectArray = cell(size(patients,1),1);
TableSize = nan(size(patients,1),1);
VIstruct = struct();
AHIdata = nan(size(patients,1),size(Evts.AHIdata2,2)+size(Evts.SpO2,2)+1);
MeanEvens = nan(length(patients), size(BreathFLDataTableLong,2)); 
MeanOdds = nan(length(patients), size(BreathFLDataTableLong,2)); 

%% patient processing <- from here on per patient
for n=1:size(patients,1)
    disp(' '); % add row space for visual clarity in command window
    if settings.analyzelist(n)==0
        displaytext=['Skipping: n=' num2str(n) ', ' char(patients(n,1))];
        disp(displaytext); set(handletext,'String',displaytext); drawnow;
        continue
    end
    displaytext=['Patient ' num2str(n) ': ' char(patients(n,1))];
    disp(displaytext); set(handletext,'String',displaytext); drawnow;

    %% load file
    SubjectArray{n} = subjectID{n,1};

    % Load breath data tables
    loadpath=[settings.AnalyzedDirectory, 'EventAnalyzed\' settings.savename '_' num2str(n)];
    loadpathC = [settings.ConvertedDirectory,settings.patients{n,1}];
    if ~(exist([loadpath,'.mat'],'file')==2)
        disp(['Skipping: n=' num2str(n), ', ' char(patients(n,1)), '. File does not exist'])
        
        % Try to at least get AHI data. Not best way to do this, but only
        % option right now.
        try
        load([settings.AnalyzedDirectory settings.savename '_' num2str(n)],'Evts'); 
        AHIdata(n,:) = [Evts.AHIdata2{2,:} Evts.SpO2{1,:} Evts.HBtotal];
        catch
        end
        
        % skip the rest and move to next subject
        continue
    end
    load(loadpath,'BreathDataTableLong', 'BreathFLDataTableLong','EAinfo','Evts','EvtsUse'); %EvtsData contains rounded position codes, PositionData contains main positioncode for each window
%     load(loadpathC,'Evts');
    
    % remove duplicates (easier to follow the code this way)
    BreathFLDataTableLong2 = BreathFLDataTableLong;
    BreathFLDataTableLong2(BreathDataTableLong.FDuplicated > 0.67...
        | BreathDataTableLong.FDuplicated>0,:) = [];
    
    BreathDataTableLong2 = BreathDataTableLong;
    BreathDataTableLong2(BreathDataTableLong.FDuplicated > 0.67...
        | BreathDataTableLong.FDuplicated>0,:) = [];
    
    TableSize(n) = size(BreathFLDataTableLong2,2);
    
    HypopIdx = BreathDataTableLong2.Etype == 4; 
%     &...
%         (BreathDataTableLong2.pos_B <0.5);
    
    SubjectArrayLong = repmat(SubjectArray(n),size(BreathDataTableLong2,1),1);
    BreathDataTableLong2.Subject = SubjectArrayLong;
  
    %VIarray
    VIsleep = BreathDataTableLong2.VI(BreathDataTableLong2.hypnog_B < 4);
    VIarray = sort(VIsleep);
    Prop = linspace(0,1,length(VIarray));
    VIstruct.(['a',SubjectArray{n}]) = [Prop', VIarray];
    
    % Mean table
    try
        % Desat slope
%         ResatSlpFtrs = nan; %struct2array(EventSignals.ProfileFtrs.ResatSlope);
%         SatAtBsln = EAinfo.EvtFtrs.DesatSlope.SatAtBsln;
       EventDepthFtrs = struct2array(EAinfo.EvtFtrs.EventDepth);
       DesatSlpFtrs = struct2array(EAinfo.EvtFtrs.DesatSlope);
       
       if sum(strcmp('PRDeltaFromMean',Evts.RespT.Properties.VariableNames))
           PRDeltaFromMean = nanmean(Evts.RespT.PRDeltaFromMean);
           PRDeltaFromMin = nanmean(Evts.RespT.PRDeltaFromMin);
           PRDeltaFromMeanAr = nanmean(Evts.ArT.PRDeltaFromMean);
           PRDeltaFromMinAr = nanmean(Evts.ArT.PRDeltaFromMin);
           PRDeltaMean = [PRDeltaFromMean,PRDeltaFromMin,PRDeltaFromMeanAr,PRDeltaFromMinAr];
           
           PRDeltaFromMean = nanmedian(Evts.RespT.PRDeltaFromMean);
           PRDeltaFromMin = nanmedian(Evts.RespT.PRDeltaFromMin);
           PRDeltaFromMeanAr = nanmedian(Evts.ArT.PRDeltaFromMean);
           PRDeltaFromMinAr = nanmedian(Evts.ArT.PRDeltaFromMin);
           PRDeltaMedian = [PRDeltaFromMean,PRDeltaFromMin,PRDeltaFromMeanAr,PRDeltaFromMinAr];
       else
           PRDeltaMean= [nan,nan,nan,nan];
           PRDeltaMedian= [nan,nan,nan,nan];        
       end
       
       % Get event level event depth averages
       if exist('EvtsUse','var')
           AvgEvDepth = nanmean(EvtsUse.RespT.AvgEvDepth);
           MinEvDepth = nanmean(EvtsUse.RespT.MinEvDepth);
           EvArea = nanmean(EvtsUse.RespT.EvArea);
           EvLevelEvDepthMean = [AvgEvDepth,MinEvDepth,EvArea];

           AvgEvDepth = nanmedian(EvtsUse.RespT.AvgEvDepth);
           MinEvDepth = nanmedian(EvtsUse.RespT.MinEvDepth);
           EvArea = nanmedian(EvtsUse.RespT.EvArea);
           EvLevelEvDepthMedian = [AvgEvDepth,MinEvDepth,EvArea];
       end
       
       % Get features of interest from Evts
       
        if ~EAinfo.NoEvents
            % Initialize some features
            MeanFeatArray(n,:) = [nanmean(BreathFLDataTableLong2{HypopIdx,1:FLendIdx},1),...
                EventDepthFtrs,DesatSlpFtrs,PRDeltaMean,EvLevelEvDepthMean];

            % Median table
            MedFeatArray(n,:) = [nanmedian(BreathFLDataTableLong2{HypopIdx,1:FLendIdx},1),...
                EventDepthFtrs,DesatSlpFtrs,PRDeltaMedian,EvLevelEvDepthMedian];
        else
            MeanFeatArray(n,:) = nan(1,length(Varnames));
            MedFeatArray(n,:) = nan(1,length(Varnames));
        end
    catch
        MeanFeatArray(n,:) = nan(1,size(MeanFeatArray,2));
        MedFeatArray(n,:) = nan(1,size(MeanFeatArray,2));
    end
    
    if isfield(settings, 'SaveBreathDataTable') && settings.SaveBreathDataTable == 1
        if ~exist('BigBreathFLDataTable', 'var')
            BigBreathDataTable = BreathDataTableLong2;
            BigBreathFLDataTable = BreathFLDataTableLong2;
        else
            % remove extra columns (temporary fix until all data has been re-converted/analyzed)
            namediffs1 = setdiff(BigBreathDataTable.Properties.VariableNames, BreathDataTableLong2.Properties.VariableNames);
            namediffs2 = setdiff(BreathDataTableLong2.Properties.VariableNames, BigBreathDataTable.Properties.VariableNames);
            
            if ~isempty('namediffs1')
                for jj = 1:length(namediffs1)
                    varname = namediffs1{jj};
                    BigBreathDataTable.(varname) = [];
                end
            end
            
            if ~isempty('namediffs2')
                for jj = 1:length(namediffs2)
                    varname = namediffs2{jj};
                    BreathDataTableLong2.(varname) = [];
                end
            end
 
            BigBreathDataTable = [BigBreathDataTable;BreathDataTableLong2];
            BigBreathFLDataTable = [BigBreathFLDataTable;BreathFLDataTableLong2];
        end
      
    end
    
    if isfield(settings, 'comparewithpnasal') && settings.comparewithpnasal == 1
        if ~exist('EventDepth_Pnasal','var') && isfield(EAinfo.EvtFtrs,'Pnasal')
            EventDepth_Pnasal = nan(size(patients,1),length(fieldnames(EAinfo.EvtFtrs.Pnasal.EventDepth)));
            EventDepth_Pneumo = nan(size(patients,1),length(fieldnames(EAinfo.EvtFtrs.Pneumo.EventDepth)));
        elseif isfield(EAinfo.EvtFtrs,'Pnasal')
            EventDepth_Pnasal(n,:) = struct2array(EAinfo.EvtFtrs.Pnasal.EventDepth);
            EventDepth_Pneumo(n,:) = struct2array(EAinfo.EvtFtrs.Pneumo.EventDepth);
        end
        
    end
    clear BreathDataTableLong BreathFLDataTableLong
    
    % AHI data
    AHIdata(n,:) = [Evts.AHIdata2{2,:} Evts.SpO2{1,:} Evts.HBtotal]; % row 2 is AHI3pa
    
    if isfield(settings, 'flowshaperepeatability') && settings.flowshaperepeatability == 1
        Etype = BreathDataTableLong2.Etype;
        Etype(Etype ~= 4) = 0;
        IdxEvens = false(size(BreathDataTableLong2,1),1);
        IdxOdds = false(size(BreathDataTableLong2,1),1);
        evtstart = find(diff(Etype) > 0);
        evtend = find(diff(Etype) < 0);
        
        if length(evtstart) > length(evtend)
            nanadd = nan(length(evtstart) - length(evtend),1);
            diffevent = [evtend; nanadd] - evtstart;
            if evtstart(end) > evtend(end)
                evtend(end+1) = length(Etype);
            end
        elseif length(evtstart) < length(evtend)
            nanadd = nan(length(evtend) - length(evtstart),1);
            diffevent = evtend - [evtstart;nanadd];
            if evtstart(1) > evtend(1)
                evtstart = [1;evtstart];
            end
        end
        
        for ee = 1:length(evtstart)
           if rem(ee,2) == 0 %even
               IdxEvens(evtstart(ee):evtend(ee)) = true;
           else %odd
               IdxOdds(evtstart(ee):evtend(ee)) = true;
           end
        end
               
        MeanEvens(n,:) = nanmean(BreathFLDataTableLong2{IdxEvens,:},1);
        MeanOdds(n,:) = nanmean(BreathFLDataTableLong2{IdxOdds,:},1);
    end
end

if 1
    FeatureArray = MeanFeatArray;
else
    FeatureArray = MedFeatArray;
end

FeatureTable = array2table(FeatureArray);
FeatureTable.Properties.VariableNames = Varnames;
FeatureTable.SubjectID = SubjectArray;

if isfield(settings, 'comparewithpnasal') && settings.comparewithpnasal == 1
    EventDepthPnasalTbl = array2table(EventDepth_Pnasal);
    EventDepthPnasalTbl.Properties.VariableNames = fieldnames(EAinfo.EvtFtrs.Pnasal.EventDepth);
    EventDepthPnasalTbl.SubjectID = SubjectArray;
    
    EventDepthPneumoTbl = array2table(EventDepth_Pneumo);
    EventDepthPneumoTbl.Properties.VariableNames = fieldnames(EAinfo.EvtFtrs.Pneumo.EventDepth);
    EventDepthPneumoTbl.SubjectID = SubjectArray;
    
    assignin('caller','EventDepthPnasalTbl',EventDepthPnasalTbl)
    assignin('caller','EventDepthPneumoTbl',EventDepthPneumoTbl)

end

if settings.SaveBreathDataTable == 1
    BreathTables = struct();
    BreathTables.BigBreathDataTable = BigBreathDataTable;
    BreathTables.BigBreathFLDataTable = BigBreathFLDataTable;
else
    BreathTables = nan;
end

AHItable = array2table(AHIdata);
AHItable.Properties.VariableNames = [Evts.AHIdata2(2,:).Properties.VariableNames ...
    Evts.SpO2.Properties.VariableNames 'HBtotal'];

if isfield(settings, 'flowshaperepeatability') && settings.flowshaperepeatability == 1
    MeanEvensTbl = array2table(MeanEvens);
    MeanEvensTbl.Properties.VariableNames = BreathFLDataTableLong2.Properties.VariableNames;
    MeanOddsTbl = array2table(MeanOdds);
    MeanOddsTbl.Properties.VariableNames = BreathFLDataTableLong2.Properties.VariableNames;
    assignin('caller','MeanEvensTbl',MeanEvensTbl)
    assignin('caller','MeanOddsTbl',MeanOddsTbl)
end