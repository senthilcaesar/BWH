%% get arousal info from all files in 1 folder (converted data)

directory = ['G:\Partners Healthcare Dropbox\SATP Group\SaraODB_HGNS2\ATS abstract 2019\Converted', '\']
myFiles = dir(fullfile(directory,'*.mat'));
ArousalTable = []; 

for k = 1:length(myFiles)
    disp(k)
    load([directory, myFiles(k).name])
    getArT
    ArousalTable = [ArousalTable; [array2table({myFiles(k).name}) struct2table(Evts.Info)]];
end
    
    TableVarNames = struct2table(Evts.Info)
    %TableVarNames.Properties.VariableNames{:}
    
    ArousalTable.Properties.VariableNames = {'ID', TableVarNames.Properties.VariableNames{:}}
    writetable(ArousalTable, 'G:\Partners Healthcare Dropbox\SATP Group\SaraODB_HGNS2\ATS abstract 2019\ArousalTable.xlsx')