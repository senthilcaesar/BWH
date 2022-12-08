function AutoScoreRespEvtsAddOn(MrangeOverride)
% RUN StartHere.m first

global settings AMasterSpreadsheet

t_start = clock;

%% Load AMasterSpreadsheet

[~,path,~] = xlsread(AMasterSpreadsheet,1,'A26');

[num,~,MasterWorksheet] = xlsread(AMasterSpreadsheet,1,'AD4:AH10003');
lastrow = find(cell2mat(cellfun(@(x)any(~isnan(x)),MasterWorksheet(:,4),'UniformOutput',false)),1,'last');
% find(1*(~isnan(cell2mat(MasterWorksheet(:,4)))),1,'last');
MasterWorksheet(lastrow+1:end,:)=[];

patients=MasterWorksheet(:,1);
Converteddirlist=MasterWorksheet(:,2);

[~,~,raw] = xlsread(AMasterSpreadsheet,2,'C20:C57');
settings.savename = char(raw{1});

MaxM = size(num,1);

%settings.Mrange=22;

try
    M = max(length(settings.Mrange));
catch
    Mrange=1:MaxM;
    M = max(Mrange);
end

success = zeros(M,1);

if exist('MrangeOverride')
     Mrange = MrangeOverride;
end

%%
for n=Mrange
    clear Evts BreathDataTable temp temp2
    try
        AnalyzedDir = [path{:} settings.savename '_' num2str(n) '.mat'];
        ConvDir=[Converteddirlist{n} patients{n} '.mat'];
        
        if exist(AnalyzedDir)==2 && exist(ConvDir)==2 
            disp(['Processing: ' num2str(n) '/' num2str(M) ': ' settings.savename '_' num2str(n) '.mat']);
            load (AnalyzedDir,'BreathDataTable','BreathFLDataTable');
            load(ConvDir,'DataEventHypnog_Mat','ChannelsList','ChannelsFs','Evts');
            
    
%% AutoScore Events
if 0 %SS turned this off, was causing problems / only analyzing one window for the whole night
    BreathDataTabletemp=BreathDataTable{1,1};
    temp2=BreathDataTabletemp;
    
    for i=1:length(temp2) % empty cells in BreathDataTable causing prob in expandtable.m
        if isempty(temp2{1,i})
            temp2{1,i}=NaN;
        end        
    end
    
    temp{1} = temp2;
    [DataEventHypnog_Mat,ChannelsList,ChannelsFs,Evts] = AutoScoreEvents(temp,DataEventHypnog_Mat,ChannelsList,ChannelsFs,Evts,1);
else
    [DataEventHypnog_Mat,ChannelsList,ChannelsFs,Evts] = AutoScoreEvents(BreathDataTable,DataEventHypnog_Mat,ChannelsList,ChannelsFs,Evts,1);
end

   
    
    plotchanneltext = { ...
        {'Epochs'},{'Flow'},{'EventsResp'},{'EventsRespAuto'},{'VE','VEeupnea'},{'Thorax'},{'Abdomen'},{'SpO2'},{'EventsAr','WakeSleep'}...
        };
%     PlotXHzData(); % causing prob in some files
    
    settings = ImportSettings(settings,AMasterSpreadsheet);

    settings.UseAutoScoredEventsForAHI=1;
    settings.ignoreCPAPdata=1;
     %% Set position codes according to protocol
        settings.positioncodesout = PositionCodeLookup(settings.poscodesdatabase,'Output');
        settings.positioncodes = PositionCodeLookup(settings.poscodesdatabase,settings.protocol(n));
        settings.supinepositioncode = settings.positioncodes(1);
       
    
    try
        timewav=DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'Time')==1));
        settings.Fs=1/(timewav(2)-timewav(1));
        if ~settings.ignoreCPAPdata
            [CPAPoff,CPAP]=getCPAP(DataEventHypnog_Mat);
        else %a priori knowledge that CPAP is not administered
            CPAPoff=ones(size(DataEventHypnog_Mat,1),1);
            CPAP=zeros(size(DataEventHypnog_Mat,1),1);
        end
        
        RespT = EventRespSigtoRespT(DataEventHypnog_Mat(:,find(ChannelsList=="EventsRespAuto")),DataEventHypnog_Mat(:,1));
        BreathDataAvailable = AnalyzedWindowsXHzFromBDT(BreathDataTable,timewav);
        ROImask=(CPAPoff & BreathDataAvailable)*1; % analyzable period per breath analysis and CPAPoff
        
        try
            %running standard version (actual best EEG for arousal analysis based on available and good staging)
            if isfield(Evts,'EvtsAuto') && isfield(Evts.EvtsAuto,'ArT')
%                 
                %Evts.EvtsAuto.ArT;
                %do nothing
            else
                Evts.EvtsAuto.ArT = EventRespSigtoRespT(DataEventHypnog_Mat(:,find(ChannelsList=="EventsArWS")),DataEventHypnog_Mat(:,1));
            end
            Evts.EvtsAuto.RespT = RespT;
            
            [Evts.EvtsAuto]=getAHIAll(DataEventHypnog_Mat,ROImask,Evts.EvtsAuto,ChannelsList,settings,'EventsRespAuto','EventsArWS');
        end
        
        try
            %running B version (estimated best EEG for arousal analysis; for when staging is not available or good)
            if isfield(Evts,'EvtsAutoB') && isfield(Evts.EvtsAutoB,'ArT')
                %Evts.EvtsAuto.ArT;
                %do nothing
            else
                Evts.EvtsAutoB.ArT = EventRespSigtoRespT(DataEventHypnog_Mat(:,find(ChannelsList=="EventsArWSB")),DataEventHypnog_Mat(:,1));
            end
            Evts.EvtsAutoB.RespT = RespT;
            [Evts.EvtsAutoB]=getAHIAll(DataEventHypnog_Mat,ROImask,Evts.EvtsAutoB,ChannelsList,settings,'EventsRespAuto','EventsArWSB');
        end
        
        try
            %running version using manually scored arousals paired with auto scored respiratory events:
            if isfield(Evts,'EvtsAutoRespOnly') && isfield(Evts.EvtsAutoRespOnly,'ArT')
                %Evts.EvtsAuto.ArT;
                %do nothing
            else
                Evts.EvtsAutoRespOnly.ArT = EventRespSigtoRespT(DataEventHypnog_Mat(:,find(ChannelsList=="EventsAr")),DataEventHypnog_Mat(:,1));
            end
            Evts.EvtsAutoRespOnly.RespT = RespT;
            [Evts.EvtsAutoRespOnly]=getAHIAll(DataEventHypnog_Mat,ROImask,Evts.EvtsAutoRespOnly,ChannelsList,settings,'EventsRespAuto','EventsAr');

        end
        
    catch faultAHI
        disp(faultAHI.getReport);
    end
    save(AnalyzedDir,'Evts','-append')
    success(n)=1;
end

% try
%    AHI3paA = Evts.EvtsAutoRespOnly.AHIdata2{1}.AllSleepAllPahi(2)
%    AHI3paB = Evts.EvtsAutoB.AHIdata2{1}.AllSleepAllPahi(2)
% end


    catch
        
        disp(['failed:' num2str(n)]);
        
    end
end

