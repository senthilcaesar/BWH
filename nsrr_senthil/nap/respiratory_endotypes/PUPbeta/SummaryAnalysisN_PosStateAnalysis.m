function SummaryAnalysisN_PosStateAnalysis()
clc

%% load globals
global SAfontsize AMasterSpreadsheet handletext settings

%% default plotting
SAfontsize=12;
set(groot,'defaultAxesfontname','arial narrow');
set(groot,'defaultAxesfontsize',SAfontsize);
set(groot,'defaultFigureColor',[1 1 1]);
set(groot,'defaultAxesBox','off');

%% set up spreadsheet if required
summaryanalysisfilename = 'AMasterSpreadsheet.xlsx';
if isempty(AMasterSpreadsheet)
    AMasterSpreadsheet = summaryanalysisfilename; % appending _OA for oral appliance dataset;
end

%% start processing
displaytext='Starting Summary Analysis';
disp(displaytext); set(handletext,'String',displaytext); drawnow;

% read spreadsheet (Files worksheet)
[~,~,raw] = xlsread(AMasterSpreadsheet,2,'C61:C74'); %settings all
% Read spreadsheet (position codes worksheet)
[~,~,settings.poscodesdatabase] = xlsread(AMasterSpreadsheet,3,'B2:J55');
settings.positioncodesout = PositionCodeLookup(settings.poscodesdatabase,'Output');
[~,~,settings.protocol] = xlsread(AMasterSpreadsheet,1,'AF4:AF10003');

settings.directory=char(raw{1});
settings.savename=char(raw{2});
settings.getCIs=logical(raw{3});
settings.plotfigs=logical(raw{4});
settings.selecttimeofnight=logical(raw{5});
settings.selecttimeofnight_XTiles=double(raw{6});
settings.selecttimeofnight_NthXTile=double(raw{7});
settings.selectstate=double(raw{8}); %1=nrem1, etc; 4=nremALL, 5=REM; nonREM states have zero REM; specific states have>50% state X.
% plan to include a tolerance to state option in future release, i.e.
% set the permitted fraction of window to be other than state selected
settings.statetolerance=0.5; % set this manually for now
%Added code for position dependent analysis (position is captured in EvtsData.PositionVariable)
settings.selectposition = char(raw{9});
disp(['Position selected: ', settings.selectposition])
eval(['settings.Mrange=[' raw{10} '];']);
settings.directoryout=char(raw{11});
settings.nanifnoarousals=double(raw{12});

try
    settings.method=double(raw{13});
    method=settings.method;
catch me
    settings.method = 1;
    disp('defaulting to method =1 (i.e. normal with Flow, not Pes method).');
    disp('To avoid seeing this in future, add method=1,2 to row 73 of MasterSpreadsheet--Settings');
end

try
    settings.comb=eval(raw{14});
catch me
    settings.comb=[0 1 30]; %disabled by default
    % first value is used as on/off switch,
    % second value is used to swap comb polarity,
    %   1 is high first, 2 is low first
    % and third value is comb width
end

%% Run
if isfield(settings,'DoPositionAnalysisInSummary') && settings.DoPositionAnalysisInSummary==1
    %% outer loop to run various configuration sets (i.e. position and state selections)
    settings.plotfigs = 0; % forced off for this style of analysis (too many figures)
    for configset=1:9 %3
        switch configset
            case 1
                settings.selectstate = 4; % NREM
                settings.selectposition = 'Supine';
                colout = 'A';
                titleout = {'NonREM Supine'}
            case 2
                settings.selectstate = 4; % NREM
                settings.selectposition = 'Lateral'; % could set as 'NonSupine' to include prone
                colout = 'N';
                titleout = {'NonREM Lateral'}
            case 3
                settings.selectstate = 4; % NREM
                settings.selectposition = 'All';
                colout = 'AA';
                titleout = {'NonREM AllPos'}
            case 4
                settings.selectstate = 8; % All sleep
                settings.selectposition = 'Supine';
                colout = 'AN';
                titleout = {'AllSleep Supine'}
            case 5
                settings.selectstate = 8; % All sleep
                settings.selectposition = 'Lateral'; % could set as 'NonSupine' to include prone
                colout = 'BA';
                titleout = {'AllSleep Lateral'}
            case 6
                settings.selectstate = 8; % All sleep
                settings.selectposition = 'All';
                colout = 'BN';
                titleout = {'AllSleep AllPos'}
            case 7
                settings.selectstate = 5; % REM
                settings.selectposition = 'Supine';
                colout = 'CA';
                titleout = {'REM Supine'}
            case 8
                settings.selectstate = 5; % REM
                settings.selectposition = 'Lateral'; % could set as 'NonSupine' to include prone
                colout = 'CN';
                titleout = {'REM Lateral'}
            case 9
                settings.selectstate = 5; % REM
                settings.selectposition = 'All';
                colout = 'DA';
                titleout = {'REM AllPos'}
        end
        
        % run analysis for the current configuration of state and position
        clear Data DataCI DataN AHItotal success Fstates Veupnea medianV Fsupine
        for i=settings.Mrange%1:M
            disp(['Subj: ' num2str(i)])
            try
                try
                    close(1);
                end
                if settings.method==1
                    [Data{i},DataCI{i},DataN{i},varlist,AHItotal(i),Fstates{i},...
                        Veupnea(i,1),medianV(i,1),~,~,AHIdata(i,1),Fsupine(i,:)] = SummaryAnalysisOne(i);
                    [AllPositionsEventsTable{i},AllPositionsEventsTableN{i},...
                        SupineEventsTable{i},EvtDurTable{i}]=...
                        AHIdata2Tbls(AHIdata(i,1));
                elseif settings.method==2
                    [Data{i},DataCI{i},DataN{i},varlist,AHItotal(i),Fstates{i},...
                        Veupnea(i,1),medianV(i,1),~,~,AHIdata(i,:)] = SummaryAnalysisOnePes(i,settings);
                end
                success(i)=1;
            catch me
                disp(me.message);
                success(i)=0;
            end
            %pause
        end
        I=find(success==1);
        Ifail=find(success==0);
        
        %% Find missing arousal scoring
        ArousalDataMissing=zeros(max(settings.Mrange),1);
        for i=settings.Mrange
            if success(i)==1
                if DataN{i}(1,4)==0 & ...
                        DataN{i}(1,5)==0 & ...
                        (DataN{i}(1,1)>3)
                    ArousalDataMissing(i)=1;
                end
            end
        end
        IArousalDataMissing=find(ArousalDataMissing==1);
        
        %% Overwrite, n.b. changes index
        if 0
            Data=Data(I);
            DataCI=DataCI(I);
            AHItotal=AHItotal(I);
            DB=DB(I,:);
        end
        
        %%
        clear DataArray
        rowofNaNs = NaN*ones(1,length(varlist));
        for i=1:size(Data,2)
            DataArray(i,:) = rowofNaNs; %default
            try
                DataArray(i,:) = Data{i};
            catch me
                DataArray(i,:) = NaN*DataArray(1,:);
            end
            DataArray(i,DataN{i}<3)=NaN;
        end
        if settings.nanifnoarousals
            DataArray(IArousalDataMissing,:)=NaN;
        end
        %%
        clear DataArrayN
        for i=1:size(Data,2)
            DataNArray(i,:) = rowofNaNs; %default
            try
                DataNArray(i,:) = DataN{i};
            catch me
                DataNArray(i,:) = zeros(1,length(DataNArray(1,:)));
            end
            %DataNArray(i,DataN{i}<3)=NaN;
        end
        if settings.nanifnoarousals
            DataNArray(IArousalDataMissing,:)=NaN;
        end
        
        
        %%
        clear Upper Lower
        for i=1:length(Data)
            try
                Upper(i,:) = DataCI{i}(2,:);
            catch me
                Upper(i,:) = Inf*DataArray(1,:);
            end
        end
        for i=1:length(Data)
            try
                Lower(i,:) = DataCI{i}(1,:);
            catch me
                Lower(i,:) = -Inf*DataArray(1,:);
            end
        end
        CIwidth=(Upper-Lower)/2;
        %%
        clear FstatesArray
        
        for i=1:length(Data)
            for j=1:length(varlist)
                try
                    FstatesArray{j}(i,:) = Fstates{i}(:,j)';
                catch me
                    FstatesArray{j}(i,:) = [NaN NaN NaN NaN];
                end
            end
        end
        
        %% Make Table
        SummaryAnalysisTable = array2table([DataArray]);
        SummaryAnalysisTable.Properties.VariableNames = varlist;
        
        SummaryAnalysisTableN = array2table([DataNArray]);
        SummaryAnalysisTableN.Properties.VariableNames = varlist;
        
        %% write out to excel if doing the position/state analysis
        headings = {'ID', varlist{1:3}, 'LG_winNum', varlist{4:5}, 'Ar_winNum', varlist{6:8}, 'UA_winNum'};
        
        DataArray(Ifail,:)=NaN; DataNArray(Ifail,:)=NaN; % fill failed spaces
        dataout = [[1:size(DataArray,1)]',DataArray(:,1:3),DataNArray(:,1),...
            DataArray(:,4:5),DataNArray(:,5),...
            DataArray(:,6:8),DataNArray(:,6)];
        
        xlswrite(AMasterSpreadsheet,titleout,4,[colout num2str(3)]);
        xlswrite(AMasterSpreadsheet,headings,4,[colout num2str(4)]);
        xlswrite(AMasterSpreadsheet,dataout,4,[colout num2str(5)]);
        
        % 1, ID
        % 2, LG1
        % 3, LGn
        % 4, delay
        % 5, LG_winNum
        % 6, VRA,  DataArray(:,4)
        % 7, ArThres
        % 8, Ar_winNum
        % 9, Vpassive
        % 10,Vactive
        % 11,Vcomp
        % 12,UA_winNum
       
        %% Assign in base
        if settings.comb(1)
            assignin('base',strcat('SummaryAnalysisTable','_',num2str(settings.comb(2))),SummaryAnalysisTable);
            assignin('base',strcat('SummaryAnalysisTableN','_',num2str(settings.comb(2))),SummaryAnalysisTableN);
            assignin('base',strcat('FstatesArray','_',num2str(settings.comb(2))),FstatesArray);
            assignin('base',strcat('Fsupine','_',num2str(settings.comb(2))),Fsupine);
            assignin('base','AllPositionsEventsTable',AllPositionsEventsTable);
            assignin('base','AllPositionsEventsTableN',AllPositionsEventsTableN);
            assignin('base','SupineEventsTable',SupineEventsTable);
            assignin('base','EvtDurTable',EvtDurTable);
        else
            assignin('base','SummaryAnalysisTable',SummaryAnalysisTable);
            assignin('base','SummaryAnalysisTableN',SummaryAnalysisTableN);
            assignin('base','FstatesArray',FstatesArray);
            assignin('base','Fsupine',Fsupine);
            assignin('base','AllPositionsEventsTable',AllPositionsEventsTable);
            assignin('base','AllPositionsEventsTableN',AllPositionsEventsTableN);
            assignin('base','SupineEventsTable',SupineEventsTable);
            assignin('base','EvtDurTable',EvtDurTable);
        end
        
        %% Tidy
        clear i I IArousalDataMissing Ifail j success
        
        %% Saving
        displaytext='Saving Summary Analysis';
        disp(displaytext); set(handletext,'String',displaytext); drawnow;
        if 0 % this was the normal way
            savefilename = [settings.directoryout 'SummaryAnalysis_' datestr(now,'yyyymmdd HHMMSS')];
            savefilename(end-6)='T';
            save(savefilename)
        else % this is for position and state analysis
            savefilename = [settings.directoryout  'SummaryAnalysis_' settings.savename '_configset_', num2str(configset)];
            save(savefilename)
        end
        
    end
end

%% Complete
displaytext='Completed Summary Analysis';
disp(displaytext); set(handletext,'String',displaytext); drawnow;
