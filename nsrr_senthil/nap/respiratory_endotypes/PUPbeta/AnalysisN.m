function AnalysisN(MrangeOverride)
global settings AMasterSpreadsheet

%% start processing
t_start = clock;
displaytext='Starting up analysis';
disp(displaytext);

%% Turn diary logging on
if isfield(settings,'NoKeepAnalysisLog') && settings.NoKeepAnalysisLog==1
    % no log file will be recorded
    disp('Analysis log file will not be created');    
else
    HowSoonIsNow = datestr(datetime('now'),'yyyymmdd_HHMM');
    if ~(exist([settings.workdir,'Logs'], 'dir') == 7)
        mkdir([settings.workdir,'Logs']);
    end
    diaryfilename = [settings.workdir filesep 'Logs' filesep 'AnalysisLog_' HowSoonIsNow '.txt'];
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

%% Import Settings & Default settings 

analyzefilename = [settings.workdir, filesep, 'AMasterSpreadsheet.xlsx'];
if isempty(AMasterSpreadsheet) % if this isn't already set
    AMasterSpreadsheet = analyzefilename;
end

if ~isfield(settings,'ImportedSettingsComplete')
settings = ImportSettings(settings,AMasterSpreadsheet); 
end

analyzelist = settings.analyzelist;
patients=settings.patients;
TotalNumPts=size(patients,1);

PtRange = 1:1:TotalNumPts; %normally
PtRange2 = PtRange(analyzelist==1);
Mrange = PtRange2;
if isfield(settings,'Mrange')
    Mrange=settings.Mrange;
end
if exist('MrangeOverride')
    Mrange=MrangeOverride;
end

%%
if isfield(settings,'parallelAnalysis') && settings.parallelAnalysis==1 && length(Mrange)>1
    poolobj = gcp('nocreate');
    if ~isempty(poolobj)
        delete(poolobj) % to shutdown all workers, if running already and start fresh.
    end
    settings1=settings;
    AMasterSpreadsheet1=AMasterSpreadsheet;
    
    myCluster = parcluster('local');
    if isfield(settings,'NumWorkers')&& settings.NumWorkers>1 % NumWorkers can be set in StartHere
        NumWorkers=settings.NumWorkers;
    else
        NumWorkers=myCluster.NumWorkers-1;
    end
    
    disp(['Parallel Workers used:' num2str(NumWorkers) '/' num2str(myCluster.NumWorkers)])
     
        
    parfor (ii=1:numel(Mrange),NumWorkers)
        n1=Mrange(ii);
        Analysis(n1,settings1,AMasterSpreadsheet1);
    end
else
    for ii=1:numel(Mrange)
        n1=Mrange(ii);
        Analysis(n1);
    end
end

%%
%temp1=temp;
delta_t = etime(clock, t_start); % delta in seconds
D = duration(0,0,round(delta_t),'Format','hh:mm:ss'); % duration in HH:MM:SS
displaytext = ['Analysis complete. Total time: ', char(D), ' (hh:mm:ss)'];
disp(displaytext);

diary off

end