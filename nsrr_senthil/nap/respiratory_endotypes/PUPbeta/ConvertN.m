function ConvertN(MrangeOverride)
clc
%% global variables and settings
global AMasterSpreadsheet handletext settings HarmonizedChannelSpreadsheet

%% Turn diary logging on
% should add something to the command window so that we can track back
% things like the machine and time of particular operations
if isfield(settings,'NoKeepConvertLog') && settings.NoKeepConvertLog==1
    % no log file will be recorded
    disp('Convert log file will not be created');
else
    HowSoonIsNow = datestr(datetime('now'),'yyyymmdd_HHMM');
    if ~(exist([settings.workdir,'Logs'], 'dir') == 7)
        mkdir([settings.workdir,'Logs']);
    end
    diaryfilename = [settings.workdir 'Logs' filesep 'ConvertLog_' HowSoonIsNow '.txt'];
    diary(diaryfilename); diary on
    try % write machine name at top of diary log
        if isunix %linux or mac
            disp(['Running on: ',getenv('HOSTNAME')]);
        else % windows
            disp(['Running on: ',getenv('COMPUTERNAME')]);
        end
    catch
        % we failed to determine the machine string, so write nothing
    end
end

%% master spreadsheet
if ~(isfield(settings,'SpreadsheetBypass') && settings.SpreadsheetBypass==1)
    convertfilename = [settings.workdir, filesep, 'AMasterSpreadsheet.xlsx'];
    if isempty(AMasterSpreadsheet)
        AMasterSpreadsheet = convertfilename;
    end
end

%% Run in StartHere now
if ~isfield(settings,'ImportSettingsComplete')
    settings=ImportSettings(settings,AMasterSpreadsheet);
end

%% Start processing
t_start = clock;
displaytext='Starting up the Convert to .mat process';
disp(displaytext); set(handletext,'String',displaytext); drawnow;


%% Which subjects?

if isfield(settings,'Mrange') %If decided above, use this
    Mrange=settings.Mrange;
else
    Mrange=find(settings.ConvertMatFlag==1);
end

if exist('MrangeOverride') %direct function input overrules all
    Mrange=MrangeOverride;
end

%%
errorlist=[];

if isfield(settings,'parallelConvert') && settings.parallelConvert==1 && length(Mrange)>1
    poolobj = gcp('nocreate');
    if ~isempty(poolobj)
        delete(poolobj) % to shutdown all workers, if running already and start fresh.
    end
    settings1=settings;
    AMasterSpreadsheet1=AMasterSpreadsheet;
    HarmonizedChannelSpreadsheet1=HarmonizedChannelSpreadsheet;
    myCluster = parcluster('local');
    if isfield(settings,'NumWorkers')&& settings.NumWorkers>1 % NumWorkers can be set in StartHere
        NumWorkers=settings.NumWorkers;
    else
        NumWorkers=myCluster.NumWorkers-1;
    end
    
    disp(['Parallel Workers used:' num2str(NumWorkers) '/' num2str(myCluster.NumWorkers)])
        
    parfor (ii=1:numel(Mrange),NumWorkers)
        try
            n1=Mrange(ii);
            
            Convert(n1,settings1,AMasterSpreadsheet1,HarmonizedChannelSpreadsheet1);
            
        catch me
            %me.message
            fname=[settings1.Filenames{n1}(1:end-4) '_XHz' '.mat'];
            displaytext=['Error converting: ' fname];
            disp(displaytext);
            % set(handletext,'String',displaytext); drawnow;
            errorlist = [errorlist n1];
        end
    end
    
else % without parallel processing--default
    
    for ii=1:numel(Mrange)
        try
            n1=Mrange(ii);
            Convert(n1);
            
        catch me
            %me.message
            fname=[settings.Filenames{n1}(1:end-4) '_XHz' '.mat'];
            displaytext=['Error converting: ' fname];
            disp(displaytext); set(handletext,'String',displaytext); drawnow;
            errorlist = [errorlist n1];
        end
    end
end


%%
delta_t = etime(clock, t_start); % delta in seconds
D = duration(0,0,round(delta_t),'Format','hh:mm:ss'); % duration in HH:MM:SS
disp(' '); % add row space for visual clarity in command window
displaytext = ['Convert complete. Total time: ', char(D), ' (hh:mm:ss)'];
disp(displaytext); set(handletext,'String',displaytext); drawnow;

if size(errorlist)>0
    displaytext=['Errors: ' num2str(errorlist)];
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
end

diary off

