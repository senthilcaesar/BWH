%% Check if stable breathing data exists in Analyzed data.
%
% Using the same AmasterSpreadsheet that was used during Analysis, this
% script works through the list of Analyzed files,
% checking if Stable Breathing Data exists

close
clear
clc

[file,path] = uigetfile('*.xlsx', 'Select the AMasterSpreadsheet');

% first, set up the job list
AMasterSpreadsheet = [path,file];

% Settings
[~,~,rawSettings] = xlsread(AMasterSpreadsheet,2,'C20:C75');
outputDir = rawSettings{18};
savename = rawSettings{1};

% Files
[~,rawFiles,~] = xlsread(AMasterSpreadsheet,1,'AD4:AE2100');
% trim and tidy
DirNames = rawFiles(:,2);
FileNames = rawFiles(:,1);

errorlog = [];
progressbar('Patient','Breath');
% then work through the list
for i=1:length(FileNames) 
    progressbar(i/(length(FileNames)+1),[]);
    % load the Time and StableBreathing data from converted file
    disp('.'); % add line space in command window for ease of viewing
    AnalyzedFile = [outputDir, 'EventAnalyzed\', savename, '_', num2str(i) ,'.mat'];
    str = ['File: ', savename, '_', num2str(i), ', loading ...'];  disp(str);
    try
        load(AnalyzedFile, 'BreathDataTableFulls');
    catch loaderror
         errorlog{end+1,1} = [savename, '_', num2str(i)];
         errorlog{end,2} = loaderror.message;
    end
    VarNames = BreathDataTableFulls.Properties.VariableNames;
    if isempty(find(strcmp(VarNames,'StableBreathing')==1))
        %log this        
        errorlog{end+1,1} = [savename, '_', num2str(i)];
        errorlog{end,2} = ['No SB in data'];
    end
end


keyboard
progressbar(1,1);
save([outputDir, 'EventAnalyzed\noSB_errorlog.mat'],'errorlog');
disp('All done.')

