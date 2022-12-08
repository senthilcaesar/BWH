function getFDfromAnalyzedMats()
%% this function outputs average FD values.
% it prompts for a folder of analyzed data,
% then traverses across all .mat files within that folder,
% (each file should be the output from one pt analysis)
% generating average FD values for Wake, AllSleep, NREM, and REM sleep.
% 
% this function calls code from PuPBeta, so that will need to be on path
%
%
% ToDo:
%    add FD value for NREM_supine
%
%


PuPDir = 'M:\Dropbox\PUPbeta_git\PUPbeta'; % <--- change this to suit yourself
addpath(genpath(PuPDir));

Results = table(NaN, NaN, NaN, NaN, 'VariableNames', {'Wake','AllSleep','NREM', 'REM'});
Results(1,:) = [];
Names = table({''},'VariableNames', {'ID'});
Names(1,:) = [];

% get directory listing of .mat files
selpath = uigetdir(pwd, 'Select analyzed folder');
% just grab all the zip files in the uiDir and work through each file
dirx = dir(fullfile([selpath filesep '*.mat'])); %dirx(1:2)=[];
for i=1:length(dirx) % i=34
    disp(i); % add line space in command window for ease of viewing
    try
        % open mat file each, in turn
        T = load([dirx(i).folder, filesep, dirx(i).name]);
        
        % if T.BreathDataTable is an array of cells, need to make big table
        if iscell(T.BreathDataTable)
            [BreathDataTable2,~,~,~,~]=GetNonOvlappedVE(T.BreathDataTable,[]);
        end
        
        % get averagew values for major sleep stages
        AllSleepBB = (BreathDataTable2.hypnog_B <4); % index of sleep breaths
        AllSleepFD = nanmedian(BreathDataTable2.FlowDrive(AllSleepBB));
        REMSleepBB = (BreathDataTable2.hypnog_B == 3); % index of REM breaths
        REMSleepFD = nanmedian(BreathDataTable2.FlowDrive(REMSleepBB));
        NREMSleepBB = (BreathDataTable2.hypnog_B < 3); % index of sleep breaths
        NREMSleepFD = nanmedian(BreathDataTable2.FlowDrive(NREMSleepBB));
        WakeBB = (BreathDataTable2.hypnog_B == 4); % index of wake breaths
        WakeFD = nanmedian(BreathDataTable2.FlowDrive(WakeBB));
        
        % save the [filename and/or spreadsheet index, PtAvgFD]
        Results{end+1,:} = [WakeFD,AllSleepFD,NREMSleepFD,REMSleepFD];
        Names{end+1,:}= {dirx(i).name};

    catch runbug
        keyboard
        runbug.getReport
    end
end

% attach Names to front of Results
DataOut = [Names,Results]

% save the output
savename = [dirx(i).folder, filesep, 'FD_Summary.mat'];
save(savename, 'DataOut');

% that's all for now
end
