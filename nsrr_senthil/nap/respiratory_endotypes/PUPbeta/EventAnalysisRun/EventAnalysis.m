function EventAnalysis(MrangeOverride)
%% Was called ...AA
%% Needs to be able to run without "Analyzed" data, i.e. straight from converted file.

%%
clc

%% Turn diary logging on
if 0
diaryfilename = 'EventAnalysisAA_Log.txt';
diary(diaryfilename);
diary on
end

%% global variables and settings
global AMasterSpreadsheet handletext ChannelsList settings n

% if the AMasterSpreadsheet is emtpy, try pointing it to something
if isempty(AMasterSpreadsheet)
    AMasterSpreadsheet = 'AMasterSpreadsheet.xlsx';
end

%% start processing
t_start = clock;
displaytext='Starting EventAnalysis'; 
disp(displaytext); set(handletext,'String',displaytext); drawnow;

%% Import Analysis settings (from Master worksheet table)
% [num,patients,~] = xlsread(AMasterSpreadsheet,1,'AD4:AH10003');
% analyzelist = logical(num(:,2));
% %invertflowlist = logical(num(:,1));

%% Replaced with ImportSettings
settings = ImportSettings(settings,AMasterSpreadsheet);
analyzelist =settings.analyzelist;
patients=settings.patients;

if exist('MrangeOverride')
    analyzelist=analyzelist*0;
    analyzelist(MrangeOverride)=1;
end

%% Special EventAnalysis Settings -- Defaults
%     settingsDefault.selectstate=8; %9 = all events even starting in wake, 8 = all sleep, 4 = NREM only, 5 = REM, 1=N1, 2=n@, 3=N3; %not yet operational
  
  settings.selectstate = settings.selectstateEventAnalysis; %overwrite   
  
  disp(['Warning: settings.selectstate = ' num2str(settings.selectstate) ]);
  disp(['State Selection for Events is now altered by settings.selectstateEventAnalysis in StartHere or in your settings.csv file']);
    
%% Load Decompostion models
Models = load('workspaceWF2','MdlSin','MdlExp','MdlSqu');

%% patient processing <- from here on per patient

runthese = find(analyzelist'==1);
for n=runthese
    displaytext=['Patient ' num2str(n) ': ' char(patients(n,1))];
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
    clear Evts SigT Boxes Ensembles EAinfo datatosave BreathDataTable2
    try
        %% load files
        directoryn=settings.ConvertedDirectory;
        % load the DataHypnogMat from the converted folder
        disp('Loading the Converted file');
        MATfilename=[directoryn char(patients(n,1))];
        
        temp=load(MATfilename);
        if isfield(temp, 'DataEventHypnog_Mat') % convert DataEventHypnog_Mat to SigT
            SigT=array2table(temp.DataEventHypnog_Mat);
            SigT.Properties.VariableNames = temp.ChannelsList;
        else
            SigT=temp.SigT;
            ChannelsFs=temp.ChannelsFs;
            Info=temp.Info;
        end
        %ChannelsList=SigT.Properties.VariableNames;
        clear temp
        
        % load the output mat from the analyzed folder
        disp('Loading the Analyzed file'); % note BreathTables are loaded up at row 39 ?
        
%         if ~isfield(settings, 'OutputDataDirectory') % Added by DLM [then removed by SS :)] I tried running AnalyzeEvents before anything else, and this field wasn't populated.
%             settings.OutputDataDirectory = settings.AnalyzedDirectory;
%         end
        loadpath=[settings.AnalyzedDirectory, settings.savename '_' num2str(n)];
        datatoload = load(loadpath,'BreathDataTable', 'BreathFLDataTable', 'Evts');
        
        datatoloadlist = fieldnames(datatoload);
        for i=1:length(datatoloadlist)
            if isfield(datatoload,datatoloadlist{i}) && iscell(datatoload.(datatoloadlist{i})) && size(datatoload.(datatoloadlist{i}),1)==1 ...
                    && size(datatoload.(datatoloadlist{i}),2)==1  %exist, be a 1x1 cell, contains a cell
                datatoload.(datatoloadlist{i}) = datatoload.(datatoloadlist{i});
            end
        end                  
        cellfun(@(x,y) assignin('caller',x,y),fieldnames(datatoload),struct2cell(datatoload)) %Shifts all fields from from datatosave to current workspace:
        
        % If EvtsData does not exist - open Evts
        if  ~exist('Evts','var')
            try
                temp = load(loadpath,'EvtsData'); % load it from analyzed, old cell array method
                Evts = temp.EvtsData{1}; % take it out of cell array
            end
            if ~exist('Evts','var') %remake the table from the signal if it is fully missing
                load(MATfilename,'Evts'); % load it from Converted
            end
        end

        if ~exist('Evts','var') || ~isfield(Evts,'RespT') %remake the table from the signal if it is fully missing
            Evts.RespT = EventRespSigtoRespT(SigT.EventsResp,SigT.Time);
        end
        
        %make this into a function for general use across multiple
        %functions:

        EvtsUse = Evts; %default
        if isfield(settings,'UseAutoRespInEventAnalysis') && settings.UseAutoRespInEventAnalysis>0
            switch settings.UseAutoRespInEventAnalysis
                case 1
                    EvtsUse = Evts.EvtsAutoRespOnly;
                case 2
                    EvtsUse = Evts.EvtsAuto; %All automated (?), but using best EEG based on human staging 
                case 3
                    EvtsUse = Evts.EvtsAutoB;
                case 4
                    EvtsUse = Evts.EvtsArAuto;
                case 5
                    EvtsUse = Evts.EvtsArAutoB;
            end
        end

        % % EvtsData contains rounded position codes
        % % PositionData contains main positioncode for each window
        
        if isfield(settings,'AnalyzeSnore')==1 && settings.AnalyzeSnore==1
            load(loadpath,'BreathSnoreTable');
            if sum(size(BreathSnoreTable))==2
                BreathSnoreTable=BreathSnoreTable{1};
            end
        end
  

        %% Load VE data from Analyzed Tables and generate NonOvlapped tables
        
        try
            [BreathDataTable2,~,BreathFLDataTable2,~,BreathDataTableFulls]=GetNonOvlappedVE(BreathDataTable,BreathFLDataTable);
        % test BreathDataTable2 and BreathFLDataTable2, if they aren't tables
        % and if the size is [1 1] and the value isnan, then return
        catch %BreathFLDataTable is full of NaN
            disp('Warning: BreathFLDataTable may be empty or full of NaN');
            [BreathDataTable2,~,~,~,BreathDataTableFulls]=GetNonOvlappedVE(BreathDataTable);
        end
        try
            BreathDataTable2 = join(BreathDataTable2,BreathFLDataTable2,'Keys','UniqueID');
            %BreathDataTable2=[BreathDataTable2 BreathFLDataTable2];
        end
        
        
        if ~isa(BreathDataTable2,'table') && all(size(BreathDataTable2)==[1 1]) && isnan(BreathDataTable2) && ...
               ~isa(BreathFLDataTable2,'table') && all(size(BreathFLDataTable2)==[1 1]) && isnan(BreathFLDataTable2)
            disp('Missing BreathData and BreathFLData');
            return
        end
        
        if isfield(settings,'AnalyzeSnore')==1 && settings.AnalyzeSnore==1
            [~,~,BreathSnoreTable2]=GetNonOvlappedVE(BreathDataTable,BreathSnoreTable);
            clear BreathSnoreTable
        end
        
        %% Adjust event start and stop times
        % Adjust start and stop times of events based on ventilation
%         AdjustEventTiming = 1;
        if isfield(settings,'AdjustEventStartStop') &&  settings.AdjustEventStartStop && exist('BreathDataTable2', 'var')
            EvtsUse = AdjustEventStartStop(EvtsUse,BreathDataTable2,settings.Fs);
        end
        
        %% Run EventsResptTable - important if new event times are made could be redundant otherwise 
        try
            EvtsUse = EventRespTable(SigT,EvtsUse,SigT.Properties.VariableNames);
        catch me
            disp('failed to update Events table with new start event times')
        end

       
        %% Generate Event Data
        dt = SigT.Time(2)-SigT.Time(1);
        % Remove SpO2 Artifacts / filter HR
        try
            SigT.SpO2=SpO2ArtifactReject(SigT.SpO2,dt);
        end

        try
            HRtemp = SigT.HR;
            ww = 6;
            temp = movmedian(HRtemp,round(ww/dt),'omitnan');
            %temp2 = 1*(movmean(1*isnan(HRtemp),round(ww/dt))>0.5);
            %temp(temp2==1)=NaN;
            SigT.HRfilt = temp;
        end
        
        ArSignal='EventsAr';
        if isfield(settings,'useWSanalysisToReplaceAr') && settings.useWSanalysisToReplaceAr == 1 %0=use original scoring, 1=use best EEG, 2=use "predicted best" EEG.
            ArSignal='EventsArWS';
        elseif isfield(settings,'useWSanalysisToReplaceAr') && settings.useWSanalysisToReplaceAr == 2
            ArSignal='EventsArWSB';
        end

        List = {'Time','EventsResp','Epochs',ArSignal,'EventsAr','Flow','Position','SpO2','HR','HRfilt','ArIntensity','ArIntensityOrig','WPr','ArPr','PtcCO2','Pulse'}';
        
        
        SignalsInfo = table(List);
        SignalsInfo.Found=zeros(length(List),1);
        for i=1:length(List)
            try
            ch = find(strcmp(SigT.Properties.VariableNames,List{i}));
            if ~isempty(ch)
%                 temp = DataEventHypnog_Mat(:,ch(1));
%                 eval( [SignalsList{i} '= temp;'] );
                SignalsInfo.Found(i)=1;
            end
            catch me
                
            end
        end
        
        % Regenerate EventsResp with EvtsUse
        EventsResp=0*SigT.Time;
        for m=1:length(EvtsUse.RespT.EventStart) %Resp Events
            lefti=round((EvtsUse.RespT.EventStart(m)-SigT.Time(1))*settings.Fs)+1;
            righti=lefti+round((EvtsUse.RespT.EventDuration(m))*settings.Fs);
            if lefti<1, lefti=1; end
            if righti>length(SigT.Time), righti=length(SigT.Time); end
            if EvtsUse.RespT.EventCodes(m)>1
                EventsResp(lefti:righti)=EvtsUse.RespT.EventCodes(m);
            end
        end
        SigT.EventsResp = 1*(EventsResp>0);
        
        SignalsT = SigT(:,SignalsInfo.List(SignalsInfo.Found==1));
        SignalsInfo
        
        try
            logit = @(p) log(p./(1-p));
            logitinverse = @(x) 1./(1+exp(-x));
            SignalsT.WPrLogit = logit(SignalsT.WPr);
            SignalsT.ArPrLogit = logit(SignalsT.ArPr);
        end
        
%             Time = DataEventHypnog_Mat(:,strcmp(ChannelsList,'Time')==1);
%             EventsResp = DataEventHypnog_Mat(:,strcmp(ChannelsList,'EventsResp')==1);
%             Epochs = DataEventHypnog_Mat(:,strcmp(ChannelsList,'Epochs')==1);
%             EventsAr = DataEventHypnog_Mat(:,strcmp(ChannelsList,'EventsAr')==1);
%             Flow = DataEventHypnog_Mat(:,strcmp(ChannelsList,'Flow')==1);
%             Position = DataEventHypnog_Mat(:,strcmp(ChannelsList,'Position')==1);
%             SpO2 = DataEventHypnog_Mat(:,strcmp(ChannelsList,'SpO2')==1);
%             HR = DataEventHypnog_Mat(:,strcmp(ChannelsList,'HR')==1);
%             ArIntensity = DataEventHypnog_Mat(:,strcmp(ChannelsList,'ArIntensity')==1);
%             ArIntensityOrig = DataEventHypnog_Mat(:,strcmp(ChannelsList,'ArIntensityOrig')==1);
%             try WakeSleep = DataEventHypnog_Mat(:,strcmp(ChannelsList,'WakeSleep')==1); catch; end
            

       
        
        %% Compute additional breath table variables
        if isfield(settings,'AddRIPdata') && settings.AddRIPdata==1
            for j=1:height(BreathDataTable2)
                try
                    idx1=find(SigT.Time>BreathDataTable2.Time_start(j),1);
                    idx2=find(SigT.Time>BreathDataTable2.Time_end(j),1)-1;
                    
                    ThoraxMax = max(SigT.Thorax(idx1:idx2));
                    ThoraxMin = min(SigT.Thorax(idx1:idx2));
                    AbMax = max(SigT.Abdomen(idx1:idx2));
                    AbMin = min(SigT.Abdomen(idx1:idx2));
                    
                    BreathDataTable2.DeltaThorax(j)=ThoraxMax-ThoraxMin;
                    
                    BreathDataTable2.DeltaAb(j)=AbMax-AbMin;
                    
                catch
                    BreathDataTable2.DeltaThorax(j)=nan;
                    BreathDataTable2.DeltaAb(j)=nan;
                end
            end
            
            
        end
        
        %% Prepare Drive data
        try
            if sum(isnan(BreathDataTable2.DeltaEdi)) ~= size(BreathDataTable2,1)
                if isfield(settings,'EdiAutoLinearization') && settings.EdiAutoLinearization==1
                    EdiLinearization=1;
                else
                    EdiLinearization=0;
                end
                [BreathDataTable2,DriveExponent] = VdriveFromEdi(BreathDataTable2,'VE','DeltaEdi',EdiLinearization);
            else
                DriveExponent=NaN;
            end
        catch me
            DriveExponent=NaN;
        end
        
        try    
            if sum(isnan(BreathDataTable2.DeltaPmus)) ~= size(BreathDataTable2,1)
                [BreathDataTable2] = VdriveFromEdi(BreathDataTable2,'VE','DeltaPmus',0);
            end 
        end
        try    
            if sum(isnan(BreathDataTable2.DeltaPes)) ~= size(BreathDataTable2,1)
                [BreathDataTable2] = VdriveFromEdi(BreathDataTable2,'VE','DeltaPes',0);
            end 
        end
        
        try    
            if sum(isnan(BreathDataTable2.DeltaPepi)) ~= size(BreathDataTable2,1)
                [BreathDataTable2] = VdriveFromEdi(BreathDataTable2,'VE','DeltaPepi',0);
            end 
        end
        
        try
            BreathDataTable2.Ti=BreathDataTable2.Time_mid-BreathDataTable2.Time_start;
   BreathDataTable2.Te=BreathDataTable2.Time_end-BreathDataTable2.Time_mid;
   BreathDataTable2.Ttot = BreathDataTable2.Ti + BreathDataTable2.Te;
   BreathDataTable2.VT = BreathDataTable2.VE.*BreathDataTable2.Ttot;
        end
        if isfield(settings,'PdriveFits') && settings.PdriveFits==1
     [BreathDataTable2,alpha,scalealpha,scalealpha30,scalealpha10]= PdriveAnalysis(BreathDataTable2,EvtsUse.RespT,settings)
    
        end 
    
        %% Flow Drive-prep
        try
            BreathDataTable2.FlowDrive(BreathDataTable2.FlowDrive>1)=1;
            BreathDataTable2.FlowDrive(BreathDataTable2.FlowDrive<0.1)=0.1;
            BreathDataTable2.FlowDrive(BreathDataTable2.VI<0.1)=0.1;
            BreathDataTable2.FlowDrive(BreathDataTable2.Etype==2)=0.1; %apneas are v.severe FL
            BreathDataTable2.FlowDrive(BreathDataTable2.Etype==3)=NaN; %central apneas are excluded
        end

        %% Criteria -- What events should be used per state. position, other
        
        if isfield(EvtsUse.RespT,'state')
            EvtsUse.RespT.Epochs=EvtsUse.RespT.state; %legacy, delete in 2022
        end
        
        try 
            %set in starthere, for example: settings.selectstate=8
            % settings.selectstate=0
            criteria = CriteriaSelect(n,settings,EvtsUse.RespT.Epochs,EvtsUse.RespT.PositionRaw,EvtsUse.RespT.EventCodes);
            if isfield(settings,'DISEdata') && settings.DISEdata % DV added to accommodate dise data 7/29/21
                criteria(1:end) = true;
            end
            if isfield(settings,'selectscoringcriteria') %options: InclAHI3a, InclAHI3, InclAHI4
                try
                    criteriascoring = getfield(Evts.RespT,settings.selectscoringcriteria);
                    criteria=criteria&(criteriascoring==1);
                end
            end
    
            Ncriteria = sum(criteria);
            disp(['Selected ' num2str(Ncriteria) ' of ' num2str(length(criteria)) ' events (state/pos selection)'])
        catch me
            disp('Warning: event criteria not created');
        end 
        
        %% AA Event Analysis
        try
            clear EAinfo

            EAinfo.FileName=patients{n,1}(1:end-4);
            EAinfo.NoEvents=false;
            EAinfo.DriveExponent = DriveExponent;
            
            
        if isfield(settings,'PdriveFits') && settings.PdriveFits==1
            EAinfo.PdriveExponent = alpha
            EAinfo.PdriveEdiScale = scalealpha
            EAinfo.PdriveEdiScale30 = scalealpha30
            EAinfo.PdriveEdiScale10 = scalealpha10
        end
            
            Nevents = height(EvtsUse.RespT);
            if Nevents > 2 || (isfield(settings,'DISEdata') && settings.DISEdata) % DV added to accommodate dise data 7/29/21
                
                try % try running Scotty's alternate Event Analysis plot
                    
                    if sum(strcmp(BreathDataTable2.Properties.VariableNames,'Vdr_est')) == 1
                        BreathDataTable2.VdriveEst = BreathDataTable2.Vdr_est+1;
                    else
                        BreathDataTable2.VdriveEst = nan(size(BreathDataTable2,1),1);
                    end
                    
                    try
                        maxYFOT=1; 
                        BreathDataTable2.YFOTmeanI(BreathDataTable2.YFOTmeanI>maxYFOT)=maxYFOT;
                        BreathDataTable2.YFOTmeanIW(BreathDataTable2.YFOTmeanIW>maxYFOT)=maxYFOT;
                        BreathDataTable2.YFOTmeanE(BreathDataTable2.YFOTmeanE>maxYFOT)=maxYFOT;
                        BreathDataTable2.YFOTmeanEW(BreathDataTable2.YFOTmeanEW>maxYFOT)=maxYFOT;
                    end
                    
                    [Boxes,Ensembles]=EventAnalysisRun(SignalsT,BreathDataTable2,EvtsUse.RespT,settings.BreathEnsembleMethod);
                    
                    try
                        [BoxesAuto,EnsemblesAuto]=EventAnalysisRun(SignalsT,BreathDataTable2,Evts.EvtsAuto.RespT,settings.BreathEnsembleMethod);
                        datatosave.BoxesAuto = BoxesAuto;
                        datatosave.EnsemblesAuto = EnsemblesAuto;
                    end
%                    settings.PlotEventData2Option=2;
                    EAinfo = PlotEventData2(Boxes,Ensembles,EAinfo,criteria);
                    %savelist = {'EAinfo';'Boxes';'Ensembles'};
                    datatosave.Boxes = Boxes;
                    datatosave.Ensembles = Ensembles;
                    datatosave.EAinfo = EAinfo;
                    
                    %save([settings.OutputDataDirectory 'EventAnalyzed\' settings.savename '_EvtData_' num2str(n)],savelist{:});
                    if 1% settings.saveplots [saveplots is already used for PSG figures, normally off]
                        saveas(gcf,[settings.AnalyzedDirectory 'EventAnalyzed\SS_' settings.savename '_EvtData_' num2str(n)],'jpeg')
                    
                        if  isfield(settings,'EventAnalysisForArousal') & settings.EventAnalysisForArousal==1
                        saveas(gcf,[settings.AnalyzedDirectory 'EventAnalyzed_Arousal\SS_' settings.savename '_EvtData_' num2str(n)],'jpeg')
                        end
                        
                    end
                catch me
                end
                
                %BEING PHASED OUT TO IMPROVE FLEXIBILITY:
%                 try
%                     % Ali's original Event Analysis plot
%                     [EAinfo2]=EventAnalysisRun(EAinfo,Sigs,BreathDataTable2,Models, Fs);
%                     PlotEventData(EAinfo2);
%                 catch me
%                 end
                
                if 1 %[saveplots is already used for PSG figures, normally off]
                    saveas(gcf,[settings.AnalyzedDirectory 'EventAnalyzed\AA_' settings.savename '_EvtData_' num2str(n)],'jpeg')
                end
            else
                EAinfo.NoEvents=true;
            end
            EAinfo.Success=1;
        catch meEA
            disp(meEA.message);  % short version
            %disp(meEA.getReport); % long version
            %EAinfo=NaN;
            EAinfo.Success=0;
        end      
        
        %%
        if EAinfo.NoEvents==1
            disp('No Events: Break')
            continue
        end
        
        %% Additional Analyses
        AvgEvDepth = nan(size(EvtsUse.RespT,1),1);
        MinEvDepth = nan(size(EvtsUse.RespT,1),1);
        EvArea = nan(size(EvtsUse.RespT,1),1);
        Idx = nan(size(EvtsUse.RespT,1),1);
        eventrespend = 101;
        ventilation=Boxes.VIStairs*100;
        RVentLim = 80; % ventilation within which we calc event depth;
        LVentLim = 80; % ventilation within which we calc event depth;
        
        try
        for kk = 1:size(EvtsUse.RespT,1)
            %Find Baseline from eupneic flow in 10 sec before scored event start 
            idxtemp = Boxes.EventsResp(kk,:);
            eventrespstart = eventrespend-find(idxtemp(eventrespend:-1:1)==0,1,'first')+1;
            evtidx = false(1,length(idxtemp));
            evtidx(eventrespstart:eventrespend) = true;
            
            % Average event depth per event SIMPLE
            AvgEvDepth(kk) = 100-nanmean(ventilation(kk,evtidx));
            
            % Average event depth per event COMPLEX: isolates for flow <
            % threshold
            EupBreathTemp = find(ventilation(kk,eventrespend:-1:1) <= RVentLim, 1, 'first');
            EupBreathR_mean = eventrespend - EupBreathTemp+1;
            
            EupBreathTemp = find(ventilation(kk,eventrespstart:1:EupBreathR_mean) < LVentLim, 1, 'first');
            EupBreathL_mean = EupBreathTemp+eventrespstart-1;

            if ~isempty(EupBreathR_mean) && ~isempty(EupBreathL_mean) 
                AvgEvDepth(kk) = 100-nanmean(ventilation(kk,EupBreathL_mean:EupBreathR_mean));

                % check that minimum value is within event indices
%                 [~,minIdx] = min(ventilation(kk,evtidx));
%                 if minIdx < EupBreathR_mean && minIdx > EupBreathL_mean
%                 else
%                     disp('WARNING: Minimum ventilation is outside the bounds of the event in event depth calculation')
%                 end
            else 
                AvgEvDepth(kk) = nan;
            end
            
            % Min event depth per event
            idx =find(Boxes.EventsResp(kk,100:-1:1)<1,1);
            idx(isempty(idx))=100;
            idx(idx<10)=10;
            Idx(kk)=101-idx;

            %Find minimum of event 10 sec
            winlength=9;
            [MinEvDepth(kk) y_idx] = min(mean(buffer(ventilation(kk,Idx(kk)+1:100),winlength+1,winlength,'nodelay')));
%             minpt=Idx(kk)+y_idx;
                        
            % Area under curve of event
            EvArea(kk) = nansum(100-ventilation(kk,evtidx));
            
            % Add event depth features to Evts (using EvtsUse for now but might
            % change later if appropriate) DV 8/2/2021
            EvtsUse.RespT.AvgEvDepth = AvgEvDepth;        
            EvtsUse.RespT.MinEvDepth = MinEvDepth;
            EvtsUse.RespT.EvArea = EvArea;
            
            if 0
                figure(20),clf, plot(Boxes.Time(kk,:), Boxes.VIStairs(kk,:)), hold on
                plot(Boxes.Time(kk,evtidx), Boxes.VIStairs(kk,evtidx))
                plot(Boxes.Time(kk,EupBreathL_mean:EupBreathR_mean), Boxes.VIStairs(kk,EupBreathL_mean:EupBreathR_mean))
            end
            clear idxtemp evtidx
        end
        catch me
            disp('Could not calculate respiratory event properties (e.g event depth, min event VI)')
        end
        %% Add on to get ensemble flowshape event signal; Sara Op de Beeck Project
        try
            [BreathSigT,BreathBox] = EventFlowShape(BreathDataTable2,Evts,SigT,settings);
            datatosave.BreathSigT = BreathSigT;
            datatosave.BreathBox = BreathBox;
            saveas(gcf,[settings.AnalyzedDirectory 'EventAnalyzed\AA2_' settings.savename '_EvtData_' num2str(n)],'jpeg')
        catch
        end
        
        %% Obstructive versus central event features
        try
            I= 1*(EAinfo.EnTime>EAinfo.deltaT2 & EAinfo.EnTime<=0);
            I(I==0)=NaN;
            [EAinfo.OvsC.VINadirV,Imin]=min(Ensembles.VI.*I');
            EAinfo.OvsC.FlowDriveNadirV=Ensembles.FlowDrive(Imin);
            EAinfo.OvsC.FlowDriveBaselineV=interp1(EAinfo.EnTime,Ensembles.FlowDrive,EAinfo.deltaT2);
        catch
            disp(['Failed: OvsC analysis of FlowDriveNadirV']);
        end
      
        %% Event Features by Dan Vena 
        try
            [EAinfo.EvtFtrs]=EventFeaturesRunDV(Boxes,Ensembles,EvtsUse.RespT); % only containts event depth and desat slope - merge with above when above is complete
        catch
            disp(['Failed: EventFeaturesRunDV']);
        end

        %% Waveform Breakdown Analysis
        try
            %just started, work on this.
            WBploton=1;
            % Time=[-100:1:100]'; %should be -100:1:100 but not sure where this is stored, check EAinfo.EnTime, also Ensembles.Time
            % Models = load('workspaceWF2'); % this was done at function start 
            % Models are from importing workspaceWF2.mat
            [EAinfo.EvtFtrs.FsinEst,EAinfo.EvtFtrs.FexpEst,EAinfo.EvtFtrs.FsquEst,EAinfo.EvtFtrs.HarmonicPeaks,EAinfo.EvtFtrs.T1] = ...
                WaveformBreakdown(Ensembles.VI,Ensembles.Time,Models,WBploton);
            
            temp=diff(EvtsUse.RespT.EventEnd);
            temp(temp>1.67*EAinfo.EvtFtrs.T1)=[];
            EAinfo.EvtFtrs.EvtCov=nanstd(temp)./nanmean(temp);
            EAinfo.EvtFtrs.EvtCovN=length(temp);
        catch
            disp('Failed: Waveform Breakdown Analysis');
        end
        
        %% Comparing Pneumotach to Nasal Pressure
        if isfield(settings, 'comparewithpnasal') && settings.comparewithpnasal==1
           pnasalfilelocn = [settings.workdir,'NasalP\Analyzed\EventAnalyzed\',settings.savename '_' num2str(n),'.mat'];
           if isfile(pnasalfilelocn)
               PnasalVars = load(pnasalfilelocn, 'Boxes','Ensembles');
               PnasalNanIdx = isnan(PnasalVars.Boxes.VI); % put nans from pnasal signal into pneumo signal
               PneumoNanIdx = isnan(Boxes.VI); % put nans from pneumo signal into pnasal signal
               PneumoBoxes = Boxes;
               PneumoBoxes.VI(PnasalNanIdx) = nan;
               PnasalVars.Boxes.VI(PneumoNanIdx) = nan;
               try
                   [EAinfo.EvtFtrs.Pnasal]=EventFeaturesRunDV(PnasalVars.Boxes,Ensembles,EvtsUse.RespT);
                   [EAinfo.EvtFtrs.Pneumo]=EventFeaturesRunDV(PneumoBoxes,Ensembles,EvtsUse.RespT);
               catch
               end

           end
        end
        
        %% Space for exploratory detailed arousal threshold analysis.
        if 0
          %%
            settings2=[3 0 1 0 0 0 1 0.2];
            %1:ignorefirstXbreathsaftersleeponset [4]
            %2:swapbreathsforhigherdrive [N]
            %3:Nbreathsattributabletoarousal [N]
            %4:increasingdriveonly [N]
            %5:deleteifbeloweupnea [Y]
            %6:setaseupneaifbeloweupnea [N]
            %7: use median drive not ROC curves
            %8: not used, was relevant for ROC.
            hypok = [0 1 2]; %0=N3,1=N2,2=N1,3=REM,4=W
            criteriabreath = sum(BreathDataTable2.hypnog_B==hypok,2)>0; % criteria for sleep stages in breath analysis.
            
            a = BreathDataTable2.AR3;
            xforAr = BreathDataTable2.DeltaPes;
            win = BreathDataTable2.Time0;
            
            if 1 %new
                xforAr(criteriabreath==0)=NaN;
            end
            
            [ArThres,arthres_N,~,~,IqualifyAT] = ArThresNew(a',xforAr',NaN,win,[3 0 1 0 0 0 1 0.2]);
            
            nanmedian(BreathDataTable2.DeltaPes(IqualifyAT==1));
            sum(IqualifyAT==1);
            Iar = IqualifyAT==1;
            
            Tdata = BreathDataTable2(Iar,{'DeltaPes'});
            
            [Tout,Tset,Xtile]=fNciles(Tdata,Nciles);
            
        end
        
        %% Flow shape analysis
        if ~isempty(BreathFLDataTable2)
        try
            EAinfo.FlowShape = FlowShapeAnalysisRun(BreathDataTable2,BreathFLDataTable2);
        catch me
            disp(me.message);
        end
        end
        
        %% CalculateDesatAndObstructionSeverity
        % added by EAS on 2022-02-18
        if isfield(settings,'CalculateDesatAndObstructionSeverity') && settings.CalculateDesatAndObstructionSeverity == 1
            try
                patientID = num2str(n);
                pupEvents = Evts.Table1;
                pupStartTime = Info.StartTimeInfo.StartTime;
                pupHypnogram = Evts.Hypnogram;
                channelNames = SigT.Properties.VariableNames;
                for i = 1:size(channelNames,2)
                    if (strcmp(channelNames{i}, 'SpO2'))
                        SpO2Index = i;
                        break;
                    end
                end        
                spo2SampleRate = ChannelsFs(SpO2Index);
                spo2 = SigT.SpO2;
                datatosave.DesatAndObstructionSeverityUEF = calcDesatAndObstructionSeverity(patientID, pupEvents, pupStartTime, pupHypnogram, spo2SampleRate, spo2);
            end
        end

        %% CalculateEventAreaDepthDuration
        % added by EAS on 2022-02-18
        if isfield(settings,'CalculateEventAreaDepthDuration') && settings.CalculateEventAreaDepthDuration == 1
            try
                patientNumber = n;
                datatosave.VIEventTables = CalculateEventAreaDepthDuration(BreathDataTableFulls, Evts, patientNumber);
            end
        end
        
        %% set data to save
        datatosave.settings = settings;
        datatosave.EAinfo = EAinfo;
        datatosave.Evts = Evts;
        datatosave.EvtsUse = EvtsUse;
       

        % Save raw breath data tables
        if settings.savelongbreathtables==1 
            datatosave.BreathDataTableLong = BreathDataTable2;
            datatosave.BreathFLDataTableLong = BreathFLDataTable2;
            datatosave.BreathDataTableFulls = BreathDataTableFulls;
            if settings.AnalyzeSnore==1
                datatosave.BreathSnoreLong = BreathSnoreTable2;
            end
            clear BreathDataTable2 BreathFLDataTable2 BreathDataTableFulls...
                BreathDataTable BreathFLDataTable
        end
        
        datatosavenames = fieldnames(datatosave);

        %% save method (append or one file per pt)
        if settings.SavePerSubject==1 % safe, fast, but lots of files
            displaytext=['Save analysis data to ' settings.savename '_' num2str(n) '.mat'];
            disp(displaytext); set(handletext,'String',displaytext); drawnow;
            %cd()
            savename1 = [settings.AnalyzedDirectory 'EventAnalyzed\' settings.savename '_' num2str(n) '.mat'];
            if  isfield(settings,'EventAnalysisForArousal') & settings.EventAnalysisForArousal==1
                     savename1 = [settings.AnalyzedDirectory 'EventAnalyzed_Arousal\' settings.savename '_' num2str(n) '.mat'];
            
            end
            
            
            settings.FileName_In = MATfilename;
            settings.FileName_Out = savename1;
%             for i=2:length(datatosave)
%                 eval([datatosave{i} '=' datatosave{i} ';']);  %'AHIdata=AHIdata(n);'
%             end
%             eval(['save ' savename1 ' ' strjoin(datatosave) ' -v7.3;' ]);
            save(savename1, '-struct', 'datatosave')

            for i=1:length(datatosavenames)
                if ~strcmp(datatosavenames{i},'settings') % save settings but don't clear it
                    eval(['clear ' datatosavenames{i} ';']); %'clear AHIdata;'
                end
            end
        end
        
    catch me
        displaytext=me.message; disp(me.getReport);
        disp(displaytext); set(handletext,'String',displaytext); drawnow;
    end
end


%% wrap it up
delta_t = etime(clock, t_start); % delta in seconds
D = duration(0,0,round(delta_t),'Format','hh:mm:ss'); % duration in HH:MM:SS
displaytext = ['Complete. Total Analysis time: ', char(D), ' (hh:mm:ss)'];
disp(displaytext); set(handletext,'String',displaytext); drawnow;
diary off
end


