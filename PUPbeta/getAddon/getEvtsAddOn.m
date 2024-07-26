function getEvtsAddOn(MrangeOverride)
% RUN StartHere.m first

% this function calculates AHIdata2 (for files converted early which contains only EventsData and
% not Evts or Evts.AHIdata2) and appends it to converted files.

global settings AMasterSpreadsheet

t_start = clock;

%% Load AMasterSpreadsheet

[~,path,~] = xlsread(AMasterSpreadsheet,1,'A25');

[~,~,MasterWorksheet] = xlsread(AMasterSpreadsheet,1,'T4:BF10000');
lastrow = find(1*(~isnan(cell2mat(MasterWorksheet(:,9)))),1,'last');
% lastrow = 62; % data not of the same type error for G:\Partners Healthcare Dropbox\SATP Group\PhenotypeDrive2018\PUPstart

MasterWorksheet(lastrow+1:end,:)=[];

ConvertList = cell2mat(MasterWorksheet(:,9));
patients=(MasterWorksheet(:,11));

NaNlist = (isnan(ConvertList));
ConvertList(NaNlist) = 0; % set missing or NaN as 0



[~,~,raw] = xlsread(AMasterSpreadsheet,2,'C20:C57');
settings.savename = char(raw{1});

if ~isfield(settings,'ImportedSettingsComplete')
    settings = ImportSettings(settings,AMasterSpreadsheet);
end

overwriteflag=1; %redo additional Evts fields even if present, to make sure they are all done


Mrange=1:size(ConvertList,1); %default
if isfield(settings,'Mrange')
    Mrange=settings.Mrange;
end
if exist('MrangeOverride')
    Mrange=MrangeOverride;
end

%%


for n=Mrange
    Flag=0;
    disp(['Processing: ' num2str(n) '/' num2str(Mrange(end))]);
    try
        Convertedfiledir = [path{:} patients{n} '.mat'];
        
        if exist(Convertedfiledir)==2 % if analyzed file exists
            
            % check if EVTS and AHIDATA2 exists in the Convertedfile file
            
            listOfVariables = who('-file',Convertedfiledir);
            if ~ismember('Evts', listOfVariables)   % append  if Evts doesnot exist
                Flag=1;
            else
                EvtsTemp=load(Convertedfiledir,'Evts'); % append  if Evts.AHIdata2 or ArT doesnot exist
                if ~isfield(EvtsTemp.Evts,'AHIdata2') ||  ~isfield(EvtsTemp.Evts,'ArT') || ~isfield(EvtsTemp.Evts,'ArTinfo')
                    Flag=1;
                else
                    disp (['AHIdata2 already present for: ' patients{n} '.mat'])  
                end
                
            end
            
            clear Evts SigT
            if  Flag==1 || overwriteflag==1
                disp (['Calculating AHIdata for : ' patients{n}  '.mat'])
                
                
                % coverted file
                W=load(Convertedfiledir);
                Evts=W.Evts;
                
                
                if isfield(W, 'DataEventHypnog_Mat')
                    SigT=array2table(W.DataEventHypnog_Mat);
                    SigT.Properties.VariableNames = W.ChannelsList;
                    % Field is there.  Remove it.
                    W = rmfield(W, {'DataEventHypnog_Mat','ChannelsList'});
                    
                end
                
                SigT=W.SigT;
                ChannelsList = W.SigT.Properties.VariableNames;
                
                % get CPAP data
                if ~settings.ignoreCPAPdata
                    [CPAPoff,CPAP]=getCPAP(SigT,settings,ChannelsList,0);
                else %a priori knowledge that CPAP is not administered
                    CPAPoff=ones(height(SigT),1);
                    CPAP=zeros(height(SigT),1);
                end
                
                
                
                %% SpO2 Analysis
    
    
    % get Position code
    try
        settings.positioncodes = PositionCodeLookup(settings.poscodesdatabase,settings.protocol{n});
        settings.supinepositioncode = settings.positioncodes(1);
    catch
        % above code failing if settings.protocol is char array.
        settings.positioncodes = PositionCodeLookup(settings.poscodesdatabase,settings.protocol(n,:));
        settings.supinepositioncode = settings.positioncodes(1);
    end
    
    try
        spo2Chan=find(strcmp(ChannelsList,'SpO2')==1);
        if ~isempty(spo2Chan) && sum(isnan(SigT.SpO2)) ~= size(SigT,1)
            [SpO2data,SpO2fixed]=SpO2analysis(SigT,ChannelsList,CPAPoff);
            if 1 %replace with artifact-rejected SpO2
                SigT.SpO2 = SpO2fixed;
                disp('Using fixed Spo2 (artifact-removed SpO2)')
                clear SpO2fixed;
            end
            Evts.SpO2 = SpO2data;
        else
            Evts.SpO2=NaN;
        end
    catch faultSpO2analysis
        disp(faultSpO2analysis.getReport);
        Evts.SpO2=NaN;
    end
    
    
    %% Get Arousal Table
    
    try
        ArChan=find(strcmp(ChannelsList,'EventsAr')==1);
        if ~isempty(ArChan)
            %Regular Scoring
            [Evts.ArT,Evts.ArTinfo]=getArT(SigT,ChannelsList,"EventsAr",Evts.Hypnogram,Evts.Hypnogram_t);
        end
    catch me
        'failed ArT, possibly no scored arousals'
    end
    try
        if (isfield(settings,'WSanalysis') && settings.WSanalysis==1)
            %Autoscoring 1
            [Evts.EvtsArAuto.ArT,Evts.EvtsArAuto.ArTinfo]=getArT(SigT,ChannelsList,"EventsArWS",Evts.Hypnogram,Evts.Hypnogram_t);
            try
                Evts.EvtsArAuto.RespT = Evts.RespT;
                ROImask=CPAPoff; %update
                [Evts.EvtsArAuto]=getAHIAll(SigT,ROImask,Evts.EvtsArAuto,ChannelsList,settings,'EventsResp','EventsArWS');
            end
            
            
            %Autoscoring 2
            [Evts.EvtsArAutoB.ArT,Evts.EvtsArAutoB.ArTinfo]=getArT(SigT,ChannelsList,"EventsArWSB",Evts.Hypnogram,Evts.Hypnogram_t);
            try
                Evts.EvtsArAutoB.RespT = Evts.RespT;
                ROImask=CPAPoff; %update
                [Evts.EvtsArAutoB]=getAHIAll(SigT,ROImask,Evts.EvtsArAutoB,ChannelsList,settings,'EventsResp','EventsArWSB');
            end
        end
        
    catch
        disp('failed creating Autoscored Arousal Table');
    end
    
    % nanmedian(Evts.EvtsArAuto.ArT.WSBalanceMax3(Evts.EvtsArAuto.ArT.AASMarousal==1))
    % nanmedian(Evts.EvtsArAuto.ArT.ArBalanceMax3(Evts.EvtsArAuto.ArT.AASMarousal==1))
    % nanmedian(Evts.EvtsArAuto.ArT.ArInt(Evts.EvtsArAuto.ArT.AASMarousal==1))
    % nanmedian(Evts.EvtsArAuto.ArT.ArIntOr(Evts.EvtsArAuto.ArT.AASMarousal==1))
    
    %% Get AHI, EventsRespT
    
    try
        if  (~isfield(settings, 'SiteOfCollapse') || (isfield(settings, 'SiteOfCollapse') && settings.SiteOfCollapse == 0))
            % skip this if just doing DISE stuff
            % Original Scoring
            disp('Original Scoring:');
            
            %[CPAPoff,CPAP] = getCPAP(SigT,settings,SigT.Properties.VariableNames);
            [Evts]=getAHIAll(SigT,CPAPoff,Evts,SigT.Properties.VariableNames,settings);
            %           Evts.AHIdata=AHIdata;
            %           Evts.AHIdata2=AHIdata2;
            %           Evts.HBtotal=HBtotal;
            %% Get Epoch Table-- epoch level features
            [Evts]=EpochsTable(SigT,Evts,AMasterSpreadsheet,settings,n);
        end
    catch faultAHI
        disp(faultAHI.getReport);
    end
    
    
    %% Ali's Endotypes --delta pulse rate/heart rate for resp and arousal events
    
    try
        Evts=getEvtPRMainFn(Evts,SigT);
        disp('success: Heart rate and/or pulse rate added to Evts.RespT/Evts.ArT');
    catch
        disp('failed adding heart rate and/or pulse rate to Evts.RespT/Evts.ArT');
    end
    
    
    
    %% PSG Report
    
    if isfield(settings,'PsgReport') && settings.PsgReport==1
        settings.onlykeepdataafterCPAPswitchedoff=1;
        fname=extractBefore(char(settings.Filenames(n,1)),'.');
        Evts.PSGreportSummaryT=getPSGreport(SigT,Evts,ChannelsList,settings,fname,CPAPoff);
    end
    %%
                save(Convertedfiledir,'Evts','-append')
                
                disp(['Completed AHIData2/ArT Calculation for study '  patients{n} '.mat']);
            end
        end
    catch
        disp(['Failed Evts Calculation for study '  patients{n} '.mat']);
    end
end

delta_t = etime(clock, t_start); % delta in seconds
D = duration(0,0,round(delta_t),'Format','hh:mm:ss'); % duration in HH:MM:SS
disp(' '); % add row space for visual clarity in command window
displaytext = ['Evts Calculation Complete. Total time: ', char(D), ' (hh:mm:ss)'];
disp(displaytext);
