% Analyze untable breathings segments

folder = 'E:\Dropbox (Partners HealthCare)\PhenotypeDrive2018\Analyzed\UnstableAnalyzed\';
subname = '1897'; % loop once more subjects are included
load([folder 'EvtSigs_' subname '_XHz.mat'])

DataArray = cell(99,3);

% generate table for analysis
count = 0;
StructNames1 = fieldnames(EvtSigsStruct);
for ii = 1:length(StructNames1)
    StructNames2 = fieldnames(EvtSigsStruct.(StructNames1{ii}));
    
    for jj = 1:length(StructNames2)
        count = count+1;
        datastruct = EvtSigsStruct.(StructNames1{ii}).(StructNames2{jj});
        
        DataArray{count,1} = datastruct.FileName;
        DataArray{count,2} = StructNames1{ii};
        DataArray{count,3} = datastruct.EventDepth;
        DataArray{count,4} = datastruct.EventDepth2Min;
        
    end
end

DataTable = cell2table(DataArray);
DataTable.Properties.VariableNames = {'Subject', 'Case', 'EventDepth',...
    'EventDepth2Min'};