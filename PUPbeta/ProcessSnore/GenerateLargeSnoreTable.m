function GenerateLargeSnoreTable()
clc
diaryfilename = 'AnalysisLog.txt';
analyzefilename = 'AMasterSpreadsheet.xlsx'; % append _OA for oral appliance dataset;

%% test for specific machine - i.e. the following block only runs on DLM work pc
try
[~, str] = system('hostname'); % alternate windows command is getenv('COMPUTERNAME'), mac code is different again
targetstr = 'GS528-3792';
if strncmp(str, targetstr, 10)
    close all
    clear global AMasterSpreadsheet
    clear
    clc
    addpath('C:\Users\uqdmann\Dropbox\PUPbeta_git\');
    cd('C:\Users\uqdmann\Dropbox\PUPbeta_git\FlowAndEdiPesQAO');
    disp('Detected DLM uni computer - overwriting selected variables');
    HowSoonIsNow = datestr(datetime('now'),'yyyymmdd_HHMM');
    diaryfilename = ['AnalysisLog_DLM_',HowSoonIsNow,'.txt'];
    % I use this to run different analysis on versions of matlab
    if verLessThan('matlab', '9.2')   % 2017a is ver 9.2,
        % i.e. matlab 2016b or anything older
        %analyzefilename = 'AMasterSpreadsheet_Pnasal.xlsx'; 
        analyzefilename = 'AMasterSpreadsheet_MESA_dlm.xlsx'; 
    else
        % i.e. matlab 2017a or anything newer
        analyzefilename = 'AMasterSpreadsheet_pneumotach.xlsx';
    end   
end
catch % nothing to report, just move on
end

global particle
particle='\';
if ismac()
    particle='/';
end

%% Turn diary logging on
diary(diaryfilename);
diary on

%% global variables and settings
global AMasterSpreadsheet handletext ChannelsList settings n 
% ChannelsList
%globals include: savename Pnasaldownisinsp LGfromFlowVersion saveplots sqrt_scaling plotfigure usescoredcentralapneas eventsarebreathsfullywithinmargins havescoredcentralhypops manualscoringtouchups maxdelaybreaths Fs exportresultstoxls
if isempty(AMasterSpreadsheet)
    AMasterSpreadsheet = analyzefilename; 
end

%% start processing
t_start = clock;
displaytext='Starting up analysis'; 
disp(displaytext); set(handletext,'String',displaytext); drawnow;

% read spreadsheet (files worksheet)
[num,patients,~] = xlsread(AMasterSpreadsheet,1,'AD4:AH10003');
analyzelist = logical(num(:,2));

% read spreadsheet (options worksheet)
[~,~,raw] = xlsread(AMasterSpreadsheet,2,'C20:C58');
settings.savename = char(raw{1});
settings.OutputDataDirectory = char(raw{18});
if settings.OutputDataDirectory(end)~=particle
    settings.OutputDataDirectory=[settings.OutputDataDirectory particle];
end
try
    displaytext=['Loading saved data from ' settings.savename];
    load([settings.OutputDataDirectory, settings.savename],'settings','AnalysisIndex','LGplusinfo','EventsInfo','LG_QualityInfo','DataOut','ArousalDat','fitQual','MiscEpochData','SleepData','CPAPData','StoNData','AHIdata','SpO2data','BreathDataTable','BreathFLDataTable','LocalSignals','Evtsdata');
catch me
    displaytext=['No saved data file called ' settings.savename];
end
disp(displaytext); set(handletext,'String',displaytext); drawnow;
settings.OutputDataDirectory = char(raw{18}); %reload/overwrite
settings.Fs=raw{13};
settings.protocol = patients(:,3);

settings.rerunspecificwindows=[]; %leave empty, work is needed so that if this is saved, other saved data are preserved

if ~isempty(settings.rerunspecificwindows)
    disp('Warning: rerun specific windows on');
end

%% patient processing <- from here on per patient
for n=1:size(patients,1)
    disp(' '); % add row space for visual clarity in command window
    if analyzelist(n)==0
        displaytext=['Skipping: n=' num2str(n) ', ' char(patients(n,1))];
        disp(displaytext); set(handletext,'String',displaytext); drawnow;
        continue
    end
    displaytext=['Patient ' num2str(n) ': ' char(patients(n,1))];
    disp(displaytext); set(handletext,'String',displaytext); drawnow;

    %% load file
    % Load DataEventHypnogMat
    directoryn=char(patients(n,2));
    if directoryn(end)~=particle
        directoryn=[directoryn particle];
    end
    MATfilename=[directoryn char(patients(n,1))];
    load(MATfilename,'DataEventHypnog_Mat','ChannelsList');
    
    % Load breath data tables
    loadpath=[settings.OutputDataDirectory, 'EventAnalyzed\' settings.savename '_' num2str(n)];
    load(loadpath,'BreathDataTableLong', 'BreathFLDataTableLong','BreathSnoreLong'); %EvtsData contains rounded position codes, PositionData contains main positioncode for each window
    
    % Remove duplicated
    BreathFLDataTableLong(BreathDataTableLong.FDuplicated>0.67 ...
        | BreathDataTableLong.FDuplicated2==1,:) = [];
    BreathSnoreLong(BreathDataTableLong.FDuplicated>0.67 ...
        | BreathDataTableLong.FDuplicated2==1,:) = [];
    BreathDataTableLong(BreathDataTableLong.FDuplicated>0.67 ...
        | BreathDataTableLong.FDuplicated2==1,:) = [];
    
    % Add subjects labels to Breath table
    Subjects = cellstr(repmat(patients{n,1}(1:end-4),size(BreathSnoreLong,1),1));
    BreathDataTableLong.Subject = Subjects;
    
    % Concatenate tables
    if exist('BreathDataTableAll')
        BreathDataTableAll = vertcat(BreathDataTableAll,BreathDataTableLong);
        BreathFLDataTableAll = vertcat(BreathFLDataTableAll,BreathFLDataTableLong);
        BreathSnoreAll = vertcat(BreathSnoreAll,BreathSnoreLong);
    else
        BreathDataTableAll = BreathDataTableLong;
        BreathFLDataTableAll = BreathFLDataTableLong;
        BreathSnoreAll = BreathSnoreLong;
    end
end

savepath=[settings.OutputDataDirectory, 'EventAnalyzed\' settings.savename '_All.mat'];
save(savepath,'BreathDataTableAll', 'BreathFLDataTableAll','BreathSnoreAll');


