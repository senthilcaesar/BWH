function getAHI2()
% settings
global settings AMasterSpreadsheet ChannelsList n
settings = ImportSettings(settings, AMasterSpreadsheet);
[num,patients,~] = xlsread(AMasterSpreadsheet,1,'AD4:AH10003');
analyzelist = find(num(:,2))';
settings.protocol = patients(:,3);

% loop through all patients on analyzelist
for n=analyzelist
    
    % load converted
    filedir = [patients{n,2} patients{n,1} '.mat'];
    if exist(filedir)==2
        x = matfile(filedir);
        disp(['Generating AHIdata2 for: ', patients{n,1}])
    else
        disp([filedir ' does not exist'])
        continue
    end
    
    DataEventHypnog_Mat = x.DataEventHypnog_Mat;
    ChannelsList = x.ChannelsList;
    Evts = x.Evts;
    
    % get CPAPoff data for getAHI2 function
    if ~settings.ignoreCPAPdata
%         [CPAPoff,~]=getCPAP(DataEventHypnog_Mat);
        Pmask = DataEventHypnog_Mat(:,strcmp('Pmask',ChannelsList));
        Time = DataEventHypnog_Mat(:,strcmp('Time',ChannelsList)); 
        [~,CPAPoff] = CPAPonoroff(Pmask,Time);
        CPAPoff = TurnCPAPoffForAllDataBeforeCPAPoff2(CPAPoff);
    else %a priori knowledge that CPAP is not administered
        CPAPoff=ones(size(DataEventHypnog_Mat,1),1);
%                 CPAP=zeros(size(DataEventHypnog_Mat,1),1);
    end

    % getAHI2
    % position centric code: generate settings for position - this technically varies per subject
    % so cant do it at the top
    settings.positioncodesout = PositionCodeLookup(settings.poscodesdatabase,'Output');
    settings.positioncodes = PositionCodeLookup(settings.poscodesdatabase,settings.protocol{n});
    settings.supinepositioncode = settings.positioncodes(1);
    [Position,~] = getPos(DataEventHypnog_Mat,ChannelsList,settings);
    
    % run EvtsRespTable if RespT doesnt exist
    if 1 % rerun ALL because random stuff missing%~isfield(Evts,'RespT')
        Evts = EventRespTable(DataEventHypnog_Mat,Evts,ChannelsList);
    end
    
    % the rest
    EpochsXHz=DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'Epochs')==1));
    Time=DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'Time')==1));
    [AHI.Total,Evts,rowNames]=getAHIEvtSubset(Evts,EpochsXHz,Time,Position,CPAPoff);
    AHITable1=array2table(AHI.Total,'VariableNames',rowNames,'RowNames',{'AHITotal'});
    %                 displaytext=['Total AHI: ' num2str(AHI.Total(58),4) ' per hr'];
    %                 disp(displaytext); set(handletext,'String',displaytext); drawnow;
    %                 displaytext=['NREM supine AHI: ' num2str(AHI.Total(16),4) ' per hr'];
    %                 disp(displaytext); set(handletext,'String',displaytext); drawnow;
    %                 
    %Event subset, e.g. 4pc AHI
    desatlist=[3 4]; %note 4percent will not include arousal

    [AHI.ThreePA,Evts,~] = getDesatArSubset(Evts,EpochsXHz,Time,Position,CPAPoff,desatlist(1));
    [AHI.FourP,Evts,rowNames] = getDesatArSubset(Evts,EpochsXHz,Time,Position,CPAPoff,desatlist(2));
    AHITable2=array2table(AHI.ThreePA,'VariableNames',rowNames,'RowNames',{'AHI3PA'});
    AHITable3=array2table(AHI.FourP,'VariableNames',rowNames,'RowNames',{'AHI4P'});
    AHIdata2{1}=[AHITable1;AHITable2;AHITable3];

    % Load analysis file and save all files over again, along with AHdata2 (can't use append
    % because it corrupts the file). Overwripte settings
    loadpath = [settings.OutputDataDirectory settings.savename '_' num2str(n) '.mat'];
    Varlist = who(matfile(loadpath));
    Varlist(strcmp(Varlist,'settings')) = []; % overwrite settings
    Varlist(strcmp(Varlist,'AHIdata2')) = []; % overwrite settings
    load(loadpath, Varlist{:})
    
    % Check for EvtsData
    if 1% rerun all because random stuff missing % ~exist('EvtsData', 'var')
        clear EvtsData
        EvtsData{1} = Evts;
        Varlist = [Varlist;{'EvtsData'}];
    end
    
    VarlistNew = [Varlist;{'AHIdata2'};{'settings'}];
    save(loadpath, VarlistNew{:})
    clear EvtsData
end