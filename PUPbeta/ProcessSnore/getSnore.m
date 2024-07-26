% Computes breath-level snores and saves on patient-level in 
% BreathSnoreTable in the analyzed file
% Must run convert and analyzed first

function getSnore() 

global settings AMasterSpreadsheet ChannelsList ChannelsFs
% settings1 = ImportSettings(settings,AMasterSpreadsheet); 
% settings=settings1;
% settings.DoExtra = 0;
% [num,patients,~] = xlsread(AMasterSpreadsheet,1,'AD4:AH10003');
% [analyzID,~,~] = xlsread(AMasterSpreadsheet,1,'S4:S10003');
% [~,directoryn,~] = xlsread(AMasterSpreadsheet,1,'A24:A27');
% NaNlist = (isnan(num)); num(NaNlist) = 0; % set missing or NaN as 0
% analyzelist = logical(num(:,2));
% settings.invertflowlist = logical(num(:,1));
% settings.protocol = patients(:,3);
TotalNumPts=size(settings.patients,1);
% settings.patients = patients;
analyzID = find(settings.analyzelist);
PtRangeTemp = 1:1:TotalNumPts; %normally
PtRange = PtRangeTemp(settings.analyzelist==1);
% FsSnore = settings.FsSnore;

SnoreTables = struct();

if isfield(settings,'SiteOfCollapse') && settings.SiteOfCollapse == 1
    Tag = '';
else
    Tag = 'Snore'; % This will have to change for full night PSG data (can make if statement related to System)
end

for ptnum = PtRange
    for kk = 1:length(settings.DBthresh)
        disp(['Analyzing patient: ', settings.patients{ptnum,1}, ', DB threshold: ', num2str(settings.DBthresh(kk)),'dB'])
        
        if kk == 1 % load data on first DB threshold loop
            %% Load data
            settings.filename=[char(settings.patients(ptnum,1))]; %seems to be unused. 
            directoryn=char(settings.patients(ptnum,2));
        
            % load Source snore data
            directoryS = settings.SourceDirectory;
            SnoreFilename = [directoryS settings.Filenames{ptnum,1}(1:end-4),Tag];
            temp=load(SnoreFilename);
            SnoreTemp=temp.cSnore;
            
            if isfield(settings,'SiteOfCollapse') && settings.SiteOfCollapse == 1
                Snore = SnoreTemp.values;
                FsSnore = 1/SnoreTemp.interval;
            else % This will have to change for Nox data (can make if statement related to System)
                Snore = SnoreTemp.values;
                FsSnore = 1/SnoreTemp.interval;
            end
            clear SnoreTemp
            
            % load converted data
            directoryC = settings.ConvertedDirectory;
            MATfilename=[directoryC settings.patients{ptnum,1}];
            Evts=struct();
            temp=[];
            temp=load(MATfilename);
            SigT=temp.SigT;
            Evts=temp.Evts;
%             ChannelsFs=temp.ChannelsFs;
%             ChannelsList=temp.ChannelsList;
            Flow = SigT.Flow;
%             SnoreDB = DataEventHypnog_Mat(:,strcmp(ChannelsList,'SnoreDB'));
            FlowTime = SigT.Time;
            FsFlow = 1/(FlowTime(2) - FlowTime(1));
            
            % load analyzed data
            directoryA = settings.AnalyzedDirectory;
            filenameA = [directoryA, settings.savename,'_',num2str(ptnum)];
            A = load(filenameA);
            
            if isfield(A, 'BreathDataTable2')
                BreathDataTable2 = A.BreathDataTable2;
                BreathFLDataTable2 = A.BreathFLDataTable2;
            else
                BreathDataTable = A.BreathDataTable;
                BreathFLDataTable = A.BreathFLDataTable;
                [BreathDataTable2,~,BreathFLDataTable2,~,~]=GetNonOvlappedVE(BreathDataTable,BreathFLDataTable);
            end
            
%             if isfield(settings, 'SiteOfCollapse') && settings.SiteOfCollapse % NEED TO FIX by incorporating VOTE labels in Analyze
%                 BreathDataTable2 = A.BreathDataTable2;
%                 BreathFLDataTable2 = A.BreathFLDataTable2;
%             else
%                 BreathDataTable = A.BreathDataTable;
%                 BreathFLDataTable = A.BreathFLDataTable;
%                 [BreathDataTable2,~,BreathFLDataTable2,~,~]=GetNonOvlappedVE(BreathDataTable,BreathFLDataTable);
%             end
            
            dupidx = (BreathDataTable2.FDuplicated>0.67 | BreathDataTable2.FDuplicated2>0.67); 
            BreathDataTable = BreathDataTable2(~dupidx,:);
            BreathFLDataTable = BreathFLDataTable2(~dupidx,:);
            
            %% Snore processing: Generates snore features from raw snore data
            PlotSnore=0;
%             ChannelsList = [ChannelsList,'SnoreRaw'];
%             ChannelsFs = [ChannelsFs; FsSnore];
            Fs = 1/(FlowTime(2) - FlowTime(1));
            windowselect = [];
            [SnoreStruct] = ProcessSnoreNew2(Snore,Flow,FlowTime,Fs,FsSnore,PlotSnore,windowselect);
            Freq = SnoreStruct.Freq;

            % Subject array
            SubjectID = repmat({settings.patients{ptnum,1}(1:end-4)},length(BreathDataTable.Time_start),1);
            BreathDataTable.SubjectID = SubjectID;

        end
%%        
        % Compute features based on a DB threshold
        BBtime = [BreathDataTable.Time_start,BreathDataTable.Time_mid,BreathDataTable.Time_end];
        
        [BreathSnoreTable, PSDbreath, EnvelpBreath, EnvelpBreath_n, PwelchSmooth] = ...
        ComputeSnoreFeaturesNew(SnoreStruct,BBtime,settings.DBthresh(kk));
        
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
        
        %% Compute airflow features (e.g. average event depth, minimum event ventilation)
        
        %% Compute window level features (can only do this on DISE, otherwise prohibitively large)
        if isfield(settings,'includewindowedleveldata') && settings.includewindowedleveldata
            [BreathSnoreTable, PSDbreath, EnvelpBreath, EnvelpBreath_n, PSDsmooth] = ...
                ComputeWinLevelFeatures(SnoreStruct,BBtime,settings.DBthresh(kk));
            
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
        SnoreTables.(['BreathSnoreTable',num2str(settings.DBthresh(kk))]) = ...
               BreathSnoreTable; 
        SnoreTables.(['PSDbreath',num2str(settings.DBthresh(kk))]) = ...
            PSDbreath;
        SnoreTables.(['EnvelpBreath',num2str(settings.DBthresh(kk))]) = ...
            EnvelpBreath;
        SnoreTables.(['EnvelpBreath_n',num2str(settings.DBthresh(kk))]) = ...
            EnvelpBreath_n;
        SnoreTables.(['PwelchSmooth',num2str(settings.DBthresh(kk))]) = ...
            PwelchSmooth;
        SnoreTables.Freq = Freq;

    end
    clear SnoreStruct BreathSnoreTable PSDbreath EnvelpBreath EnvelpBreath_n
    
    %% Save tables for each patients in Analyzed folder
    A.SnoreTables = SnoreTables;
    save(filenameA,'-struct','A')
    
end

end
