function [PredSnoreMLR_All, PredSnoreLR_All, PredSnoreFeat_All, success] = RunPredSiteFromSnore(MrangeOverride)

global settings AMasterSpreadsheet ChannelsList ChannelsFs

TotalNumPts=size(settings.patients,1);
PtRangeTemp = 1:1:TotalNumPts; %normally
Mrange = PtRangeTemp(settings.analyzelist==1);

if exist('MrangeOverride','var')
    Mrange = MrangeOverride;
end

% Initialize tables - a bit of headache, but its the most error proof
PredSnoreMLR_All = [];
PredSnoreLR_All = [];
PredSnoreFeat_All = [];

success = ones(length(Mrange),1);
for ptnum = Mrange
    disp(['Analyzing: ' settings.patients{ptnum}(1:end-4)])
    
    % load analyzed data
    directoryA = settings.AnalyzedDirectory;
    filenameA = [directoryA, settings.savename,'_',num2str(ptnum)];
    if ~(exist([filenameA,'.mat'],'file')==2)
        disp('Analyzed data does not exist')
        success(ptnum) = 0;
        continue
    end
    A = load(filenameA);
    
    if ~isfield(A,'SnoreTables')
        disp('SnoreTables do not exist')
        success(ptnum) = 0;
        continue
    end
    BreathSnoreTable = A.SnoreTables.(['BreathSnoreTable',num2str(settings.DBthres)]);
    BreathDataTable = A.BreathDataTable;
    
    if ~istable(BreathDataTable)
        [BreathDataTable2,~,~,~,~]=GetNonOvlappedVE(BreathDataTable,A.BreathFLDataTable);
        dupidx = (BreathDataTable2.FDuplicated>0.67 | BreathDataTable2.FDuplicated2>0.67); 
        BreathDataTable = BreathDataTable2(~dupidx,:);
    end
    
    
    % Pick only include loud snores
    if settings.allsnores
        %%% All snores routine with DISE data %%%
        LoudSnoreIdx = BreathSnoreTable.SnoreDBMedian_i > settings.DBthres;
        lowVIidx = BreathDataTable.VI < 1;
%         FLidx = BreathDataTable.FlowDrive < 0.6;
        SnoreBrIdx = LoudSnoreIdx & lowVIidx;
    else
        %%% DISE routine using only labeled snores %%%
        SnoreBrIdx = ~cellfun(@isempty,BreathDataTable.SiteOfCollapse);
    end
    
    [PredSnoreTblMLR,PredSnoreTblLR] = PredictSiteFromSnore(BreathSnoreTable(SnoreBrIdx,:));
  
    PredSnoreFeatTbl = BreathSnoreTable(SnoreBrIdx,:);
    
    PredSnoreMLR_All = [PredSnoreMLR_All;PredSnoreTblMLR];
    PredSnoreLR_All = [PredSnoreLR_All;PredSnoreTblLR];
    PredSnoreFeat_All = [PredSnoreFeat_All;PredSnoreFeatTbl];
    
    % Also Save In Table that's same dimensions as BreathDataTable
    PredSnoreLR_ = nan(size(BreathDataTable,1),size(PredSnoreTblLR,2));
    PredSnoreLR_(SnoreBrIdx,:) =  PredSnoreTblLR{:,:};
    PredSnoreLR = array2table(PredSnoreLR_);
    PredSnoreLR.Properties.VariableNames = PredSnoreTblLR.Properties.VariableNames;
    A.PredSnoreLR = PredSnoreLR;
    clear PredSnoreLR_ PredSnoreLR
    
    PredSnoreMLR_ = nan(size(BreathDataTable,1),size(PredSnoreTblMLR,2));
    PredSnoreMLR_(SnoreBrIdx,:) =  PredSnoreTblMLR{:,:};
    PredSnoreMLR = array2table(PredSnoreMLR_);
    PredSnoreMLR.Properties.VariableNames = PredSnoreTblMLR.Properties.VariableNames;
    A.PredSnoreMLR = PredSnoreMLR;
    clear PredSnoreMLR_ PredSnoreMLR
    
    save(filenameA,'-struct','A');

end