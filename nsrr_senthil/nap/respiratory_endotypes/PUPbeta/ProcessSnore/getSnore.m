% Computes breath-level snores and saves on patient-level in 
% BreathSnoreTable in the analyzed file
% Must run convert and analyzed first

function getSnore(MrangeOverride) 

global settings AMasterSpreadsheet ChannelsList ChannelsFs

if isfield(settings,'Mrange') %If decided above, use this
    Mrange=settings.Mrange;
else
    TotalNumPts=size(settings.patients,1);
    PtRangeTemp = 1:1:TotalNumPts; %normally
    Mrange = PtRangeTemp(settings.analyzelist==1);
end
SnoreTables = struct();

if isfield(settings,'SiteOfCollapse') && settings.SiteOfCollapse == 1
    Tag = '';
else
    Tag = 'Snore'; % This will have to change for full night PSG data (can make if statement related to System)
end

if exist('MrangeOverride','var')
    Mrange = MrangeOverride;
end
for ptnum = Mrange
        
    %% Load data         
    disp(['Analyzing patient: ', settings.patients{ptnum,1}])
    
    %%% load Source snore data %%%
    directoryS = settings.SourceDirectory;
    SnoreFilename = [directoryS settings.Filenames{ptnum,1}(1:end-4),Tag];
    
    try
        listOfVariables = who('-file',SnoreFilename);
        if ismember('cSnore', listOfVariables) % this is only the case in DISE data where snore is in the spike file
            temp=load(SnoreFilename);
            SnoreTemp=temp.cSnore;
        else % otherwise snore file too big, so it is stored in separate location
            directoryS = settings.SnoreDirectory;
            SnoreFilename = [directoryS settings.Filenames{ptnum,1}(1:end-4),Tag];
            temp=load(SnoreFilename);
            SnoreTemp=temp.cSnore;
        end
    catch
        disp('No snore data')
        continue
    end

    if isfield(settings,'SiteOfCollapse') && settings.SiteOfCollapse == 1
        Snore = SnoreTemp.values;
        FsSnore = 1/SnoreTemp.interval;
    else % This will have to change for Nox data (can make if statement related to System)
        Snore = SnoreTemp.values;
        FsSnore = 1/SnoreTemp.interval;
    end
    clear SnoreTemp

    %%% load converted data %%%
    Evts=struct();
    temp=[];
    directoryC = settings.ConvertedDirectory;
    filenameC=[directoryC settings.patients{ptnum,1} '.mat'];
    if ~isfile(filenameC)
        disp('Converted file does not exist')
        continue
    end
    temp=load(filenameC);
    SigT=temp.SigT;
    Evts=temp.Evts;
    Flow = SigT.Flow;
    FlowTime = SigT.Time;
    FsFlow = 1/(FlowTime(2) - FlowTime(1));
    SnoreDB = SigT.SnoreDB;
    Time = SigT.Time;

    %%% load analyzed data %%%
    directoryA = settings.AnalyzedDirectory;
    filenameA = [directoryA, settings.savename,'_',num2str(ptnum) '.mat'];
    if ~isfile(filenameA)
        disp('Analyzed file does not exist')
        continue
    end
    A = load(filenameA);

    if isfield(A, 'BreathDataTable2')
        BreathDataTable2 = A.BreathDataTable2;
        BreathFLDataTable2 = A.BreathFLDataTable2;
    else
        BreathDataTable = A.BreathDataTable;
        BreathFLDataTable = A.BreathFLDataTable;
        [BreathDataTable2,~,BreathFLDataTable2,~,~]=GetNonOvlappedVE(BreathDataTable,BreathFLDataTable);
    end

    dupidx = (BreathDataTable2.FDuplicated>0.67 | BreathDataTable2.FDuplicated2>0.67); 
    BreathDataTable = BreathDataTable2(~dupidx,:);
    BreathFLDataTable = BreathFLDataTable2(~dupidx,:);

    %% Snore processing: Generates snore features from raw snore data
    % Start with 0db then loop through additional DB threshold options
    PlotSnore=0;
    Fs = 1/(FlowTime(2) - FlowTime(1));
    windowselect = [];
    [SnoreStruct] = ProcessSnoreNew2(Snore,Flow,FlowTime,Fs,FsSnore,PlotSnore,windowselect);
    Freq = SnoreStruct.Freq;

    % Find shape of the SnoreDB 
    if 0 % this didn't really work out
    BBtime = [BreathDataTable.BB_i_start,BreathDataTable.BB_i_mid,BreathDataTable.BB_i_end];
    Apnea_B = false(size(Time));
    SnoreShapes = ComputeSnoreShape(Time, SnoreDB, BBtime, nan, nan, nan, nan, Apnea_B, [0 1 0]);
    SnoreShapes.Properties.VariableNames = strcat(SnoreShapes.Properties.VariableNames,'_S');
    end

    % Subject array
    SubjectID = repmat({settings.patients{ptnum,1}(1:end-4)},length(BreathDataTable.Time_start),1);
    BreathDataTable.SubjectID = SubjectID;
    
    % if it does not exist then make one
    if ~isfield(settings, 'DBthreshOpts') || isempty(settings.DBthreshOpts)
        settings.DBthresOpts = 0;
    end
    
    for kk = 1:length(settings.DBthreshOpts)
%%    
        disp([num2str(settings.DBthreshOpts(kk)),'dB'])
        
        % Compute features based on a DB threshold
        BBtime = [BreathDataTable.Time_start,BreathDataTable.Time_mid,BreathDataTable.Time_end];
        
        [BreathSnoreTable, PSDbreath, EnvelpBreath, EnvelpBreath_n, PwelchSmooth] = ...
        ComputeSnoreFeaturesNew(SnoreStruct,BBtime,settings.DBthreshOpts(kk));
        
        % get max of breath-level envelope
        BreathSnoreTable.EnvelpMax = max(EnvelpBreath,[],2);
        BreathSnoreTable.EnvelpMax_n = max(EnvelpBreath_n,[],2);
        BreathSnoreTable.PwelchMax = max(PwelchSmooth,[],2);
        
        % Break up the MFCC variables
        BreathSnoreTable_ = BreathSnoreTable; 
        Varnames = BreathSnoreTable.Properties.VariableNames;
        coeffsidx = find(contains(Varnames,'coeffsAvg') | contains(Varnames,'delAvg') | contains(Varnames,'ddelAvg'));

        for cc = coeffsidx
            vartemp = BreathSnoreTable{:,cc};
            for cc2 = 1:size(vartemp,2)
                newvarname = [Varnames{cc}(1:end-2),num2str(cc2),Varnames{cc}(end-1:end)];
                BreathSnoreTable_.(newvarname) = vartemp(:,cc2);
            end
            BreathSnoreTable_ = removevars(BreathSnoreTable_,Varnames{cc});
        end
        BreathSnoreTable = BreathSnoreTable_;
        clear BreathSnoreTable_
              
        %% Compute window level features (can only do this on DISE, otherwise prohibitively large)
        if isfield(settings,'includewindowedleveldata') && settings.includewindowedleveldata
            [BreathSnoreTable, PSDbreath, EnvelpBreath, EnvelpBreath_n, PSDsmooth] = ...
                ComputeWinLevelFeatures(SnoreStruct,BBtime,settings.DBthreshOpts(kk));
            
            % get max of breath-level envelope
            WinSnoreTable.EnvelpMax = max(EnvelpWin,[],2);
            WinSnoreTable.EnvelpMax_n = max(EnvelpWin_n,[],2);

            % Break up the MFCC variables
            WinSnoreTable_ = WinSnoreTable; 
            Varnames = WinSnoreTable.Properties.VariableNames;
            coeffsidx = find(contains(Varnames,'coeffsAvg') | contains(Varnames,'delAvg') | contains(Varnames,'ddelAvg'));

            for cc = coeffsidx
                vartemp = WinSnoreTable{:,cc};
                for cc2 = 1:size(vartemp,2)
                    newvarname = [Varnames{cc}(1:end-2),num2str(cc2),Varnames{cc}(end-1:end)];
                    WinSnoreTable_.(newvarname) = vartemp(:,cc2);
                end
                WinSnoreTable_ = removevars(WinSnoreTable_,Varnames{cc});
            end
            WinSnoreTable = WinSnoreTable_;
            clear WinSnoreTable_
        end
        
        %% Store tables in a larger struct
        SnoreTables.(['BreathSnoreTable',num2str(settings.DBthreshOpts(kk))]) = ...
            BreathSnoreTable; 
%         SnoreTables.(['SnoreShapes',num2str(settings.DBthreshOpts(kk))]) = ...
%             SnoreShapes; 
        SnoreTables.(['PSDbreath',num2str(settings.DBthreshOpts(kk))]) = ...
            PSDbreath;
        SnoreTables.(['EnvelpBreath',num2str(settings.DBthreshOpts(kk))]) = ...
            EnvelpBreath;
        SnoreTables.(['EnvelpBreath_n',num2str(settings.DBthreshOpts(kk))]) = ...
            EnvelpBreath_n;
        SnoreTables.(['PwelchSmooth',num2str(settings.DBthreshOpts(kk))]) = ...
            PwelchSmooth;
        SnoreTables.Freq = Freq;

    end
    clear SnoreStruct BreathSnoreTable PSDbreath EnvelpBreath EnvelpBreath_n
   
    %% Save tables for each patients in Analyzed folder
    A.SnoreTables = SnoreTables;
    save(filenameA,'-struct','A')
    
end

end
