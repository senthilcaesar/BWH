%% Add stable breathing data to Analyzed data.
%
% Using the same AmasterSpreadsheet that was used during Analysis, this
% script works through the list of converted_XHz data, and Analyzed files,
% retrospectively adding Stable Breathing Data to the Analyzed files.
% (this should be part of the normal Analysis process in future, so as to
% not need this hack process.)
%
% the only "settings" we need from the initial data processing, (StartHere)
% is the name and location of the AMasterSpreadsheet
% so to make it easier, all we are going to ask is the location of the
% AmasterSpreadsheet for the data you wish to work upon

close
clear
clc

[file,path] = uigetfile('*.xlsx', 'Select the AMasterSpreadsheet');
% cd(path);

%% add stable breathing to Analyzed file(s)
% specifically, BreathDataTable
% Process:
% Open the Analyzed file(s) and check if StableBreathing data exists in
% the BreathDataTable. This variable is an array of cells of length windows
% with each cell storing the BreathData for breaths within that window.
% ToDo...? do we really need it at this point?


%% add stable breathing to EventAnalyzed file(s)
% specifically, BreathDataTableFulls
% Process:
% Open the EventAnalyzed file(s) and check if StableBreathing data exists
% in the BreathDataTableFulls. If not, load the Converted file and get
% Time and StableBreathing (this may have just been added above).
% Then go back to Event Analyzed file, and find matches in breath time
% with Time, and then add the corresponding StableBreathing data.

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
    filedir = [DirNames{i}, FileNames{i}]; 
    matObj = matfile(filedir); varlist = who(matObj);
    % before loading the big parts of the file, test here if SB already exists
    ChList = load(filedir, 'ChannelsList');
    if ~isempty(find(strcmp(ChList.ChannelsList,'StableBreathing')==1))
        str = ['File: ', FileNames{i}, ', loading ...'];  disp(str);
        
        % load the converted file
        try
            DataIn = load(filedir);
            TimeChannel = find(strcmp(DataIn.ChannelsList,'Time')==1);
            StableBBChannel = find(strcmp(DataIn.ChannelsList,'StableBreathing')==1);
            TimeDataFull = DataIn.DataEventHypnog_Mat(:,TimeChannel);
            StableBBFull = DataIn.DataEventHypnog_Mat(:,StableBBChannel);
            clearvars DataIn
        catch FaultLoadingConverted
            % log error
            errorlog{end+1,1} = FileNames{i};
            errorlog{end,2} = [];
            errorlog{end,3} = FaultLoadingConverted.message;
        end
        
        % load the analyzed file
        try
            AnalyzedFile = [outputDir, 'EventAnalyzed\', savename, '_', num2str(i) ,'.mat'];
            str = ['File: ', savename, '_', num2str(i), ', loading ...'];  disp(str);
            matObj2 = matfile(AnalyzedFile); varlist2 = who(matObj2);
            load(AnalyzedFile);
        catch FaultLoadingAnalyzed
            % log error
            errorlog{end+1,1} = [];
            errorlog{end,2} = [savename, '_', num2str(i)];
            errorlog{end,3} = FaultLoadingAnalyzed.message;
        end
        
        % confirm it doen not already have StableBreathing before doing it
        VarNames = BreathDataTableFulls.Properties.VariableNames;
        if isempty(find(strcmp(VarNames,'StableBreathing')==1))
            try
                % get the BB start times
                BBtime = BreathDataTableFulls.Time_mid;
                
                % match up breath times (this problem has been solved before...)
                % find the closest TimeDataFull to each BBtime
                % takes a minute or two to run on large files
                disp('Doing breath time matching on full data ...');
                ind = NaN(length(BBtime),1);
                for bb = 1:length(BBtime)
                    progressbar([],bb/length(BBtime));
                    [ d, ind(bb) ] = min( abs( TimeDataFull - BBtime(bb) ) );
                end
                progressbar([],0);
                
                % assign the stable breathing minutes to the matched times
                StableBB = StableBBFull(ind);
                
                % then add to data table
                BreathDataTableFulls.StableBreathing = StableBB;
            catch FaultAddingSB
                % log error
                errorlog{end+1,1} = [];
                errorlog{end,2} = [savename, '_', num2str(i)];
                errorlog{end,3} = FaultAddingSB.message;
                
            end
            
            if 0
               figure(1); clf(figure(1));
               stairs(BreathDataTableFulls.Time_start, ...
               BreathDataTableFulls.WakeSleep+5); hold on
               stairs(BreathDataTableFulls.Time_start, ...
               BreathDataTableFulls.FlowDrive+3);
               stairs(BreathDataTableFulls.Time_start, ...
               BreathDataTableFulls.StableBreathing);              
            end
            
            % then save this update back into event analyzed file
            try
                disp('Saving');
                if 0 % old version, where data loaded into struct
                    save([settings.patients{Pt,2},settings.patients{Pt,1},'.mat'],'-struct','DataIn');
                else % alternate save method, using varlist 
                    save(AnalyzedFile,varlist2{:},'-v7.3');
                    str = ['File: EventAnalyzed\', savename, '_', num2str(i), ', great success']; disp(str);
                end
            catch FaultSavingData
                % log error
                errorlog{end+1,1} = [];
                errorlog{end,2} = [savename, '_', num2str(i)];
                errorlog{end,3} = FaultSavingData.message;
            end
            
        else
            str = ['Event Analyzed file: ', savename, '_', num2str(i) , ' already contains StableBreathing data']; disp(str);
            % no need to log this, it's not an error
        end
    else
        str = ['Converted file: ', FileNames{i}, ', does not have StableBreathing data']; disp(str);
        % log the error
        errorlog{end+1,1} = FileNames{i};
        errorlog{end,2} = [];
        errorlog{end,3} = 'No stable breathing data in Converted_XHz';
    end
end

keyboard
progressbar(1,1);
save([outputDir, 'EventAnalyzed\addSB_errorlog.mat'],'errorlog');
disp('All done.')

