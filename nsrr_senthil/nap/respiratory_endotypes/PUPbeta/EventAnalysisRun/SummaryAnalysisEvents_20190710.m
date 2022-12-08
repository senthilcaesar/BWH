function SummaryAnalysisEvents()


diaryfilename = 'SummaryEventAnalysisLog.txt';
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
global AMasterSpreadsheet handletext settings n 
% ChannelsList
%globals include: savename Pnasaldownisinsp LGfromFlowVersion saveplots sqrt_scaling plotfigure usescoredcentralapneas eventsarebreathsfullywithinmargins havescoredcentralhypops manualscoringtouchups maxdelaybreaths Fs exportresultstoxls
if isempty(AMasterSpreadsheet)
    AMasterSpreadsheet = analyzefilename; 
end

%% add path to necessary inbuilt functions from older matlab versions
mydir  = pwd;
    addpath([mydir particle 'Builtin_MatlabFns_confirmed_necessary' particle]);
    %addpath([mydir(1:find(mydir=='\',1,'last')) 'Builtin_MatlabFns\confirmed_necessary\']);


%% start processing
displaytext='Starting up summary'; 
disp(displaytext); set(handletext,'String',displaytext); drawnow;

% read spreadsheet (files worksheet)
[num,patients,~] = xlsread(AMasterSpreadsheet,1,'AD4:AH10003');
analyzelist = logical(num(:,2));

% Read spreadsheet (position codes worksheet)
[~,~,settings.poscodesdatabase] = xlsread(AMasterSpreadsheet,3,'B2:J7');

% read spreadsheet (options worksheet)
[~,~,raw] = xlsread(AMasterSpreadsheet,2,'C20:C52');
settings.savename = char(raw{1});
settings.OutputDataDirectory = char(raw{18});
if settings.OutputDataDirectory(end)~=particle
    settings.OutputDataDirectory=[settings.OutputDataDirectory particle];
end
settings.OutputDataDirectory = char(raw{18}); %reload/overwrite

settings.protocol = patients(:,3);

if ~isempty(settings.rerunspecificwindows)
    disp('Warning: rerun specific windows on');
end

%% patient processing <- from here on per patient
T2 = table();
for n=1:size(patients,1)
    disp(' '); % add row space for visual clarity in command window
    if analyzelist(n)==0
        displaytext=['Skipping: n=' num2str(n) ', ' char(patients(n,1))];
        disp(displaytext); set(handletext,'String',displaytext); drawnow;
        continue
    end
    displaytext=['Patient ' num2str(n) ': ' char(patients(n,1))];
    disp(displaytext); set(handletext,'String',displaytext); drawnow;

    VarToRemove={'Start','End','SyncTime','VE','Arousal','SleepStage','Position',...
        'Dsat','HR','ArousalIntensity','ArousalScore','EventVE','EnTime'};
    try
        %% load file
        loadpath = [settings.OutputDataDirectory 'EventAnalyzed\' settings.savename '_' num2str(n) '.mat'];
        load(loadpath);
        patientName=char(patients(n,1));
        patientName=patientName(1:end-4);
        
        fields = fieldnames(EventSignals);
        clear idx
        for ii=1:length(VarToRemove)
            idx(ii)=find(strcmp(fields,VarToRemove{ii}));
        end
        EventSignals2 = rmfield(EventSignals,fields(idx));
        
        EventSignals2.FileName=patientName;
                           
        T=NestedStruct2table(EventSignals2);
        
        T2(n,:)=T(1,:);
        clear T       
       
        
    catch me
        displaytext=me.message; disp(me.getReport);
        disp(displaytext); set(handletext,'String',displaytext); drawnow;
    end
end

%% Writing the results in a spreadsheet
% T2 = [T2(:,end) T2(:,1:end-1)];
savename1=[settings.OutputDataDirectory 'EventAnalyzed\EventDataSummary.xlsx'];
writetable(T2,savename1)

end


