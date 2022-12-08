function Data = getData()

%run after SummaryAnalysis
global settings AMasterSpreadsheet ChannelsList

%%
% path = [settings.workdir 'Analyzed'];
[~,path,~] = xlsread(AMasterSpreadsheet,1,'A26');
%range = [1:280];

[num,patients,~] = xlsread(AMasterSpreadsheet,1,'AD4:AH10003');
NaNlist = (isnan(num)); num(NaNlist) = 0; % set missing or NaN as 0
analyzelist = logical(num(:,2));
[~,~,raw] = xlsread(AMasterSpreadsheet,2,'C20:C57');
settings.savename = char(raw{1});

MaxM = size(num,1);

try
    M = max(settings.Mrange);
catch 
    settings.Mrange=1:MaxM;
    M = max(settings.Mrange);
end

%%
clear x T
% T.AHIorig = nan(M,1);
% T.FhypopneasOrig = nan(M,1);
% T.Fhypopneas4 = nan(M,1);
% T.AHI4 = nan(M,1);
% HBtotal = nan(M,1);
success = zeros(M,1);
for i=settings.Mrange
    
    convertedfile = [patients{i,2},patients{i,1},'.mat'];
    
    if exist(convertedfile) == 2
        temp=load(convertedfile);
        DataEventHypnog_Mat=temp.DataEventHypnog_Mat;
%         Evts=temp.Evts;
%         ChannelsFs=temp.ChannelsFs;
        ChannelsList=temp.ChannelsList;
        EpochsXHz = DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'Epochs')==1));
        SaO2XHz = DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'SpO2')==1));
        dt = DataEventHypnog_Mat(2,1)-DataEventHypnog_Mat(1,1);
        % get ODI3/ODI4 in sleep only states
        [ODI3,ODI4]=CalcODI(SaO2XHz,dt,(EpochsXHz ~= 4));
        T.filename(i,1)=patients(i,1);
        T.ODI3(i,1)=ODI3;
        T.ODI4(i,1)=ODI4;        
    end
   
    
    filedir = [path{:} settings.savename '_' num2str(i) '.mat'];
    if exist(filedir)==2
    x = matfile(filedir);

    if sum(strcmp('AHIdata2',fieldnames(x)))==1
        AHIdata2 = x.AHIdata2;
    else
        disp(['AHIdata2 does not exist for patient: ' num2str(i)]);
        continue
    end
    if iscell(AHIdata2)==0 && isnan(AHIdata2)
        disp(['failed: ' num2str(i)]);
        continue
    end
    
    AHIdata2 = AHIdata2{1};
    T.AHI(i,1) = AHIdata2.AllSleepAllPahi(1);
    %T.AHI3pa(i,1) = AHIdata2.AllSleepAllPahi(2);
    T.AHI4(i,1) = AHIdata2.AllSleepAllPahi(3);
    T.Fhypopneas(i,1) = (AHIdata2.AllSleepAllPHyOI(1) + AHIdata2.AllSleepAllPHyCI(1)) / AHIdata2.AllSleepAllPahi(1);
    %T.Fhypopneas3pa(i,1) = (AHIdata2.AllSleepAllPHyOI(1) + AHIdata2.AllSleepAllPHyCI(1)) / AHIdata2.AllSleepAllPahi(1)
    T.Fhypopneas4(i,1) = (AHIdata2.AllSleepAllPHyOI(3) + AHIdata2.AllSleepAllPHyCI(3)) / AHIdata2.AllSleepAllPahi(3);
   
    T.FcentralAHI(i,1) = (AHIdata2.AllSleepAllPApCI(1) + AHIdata2.AllSleepAllPHyCI(1)) / AHIdata2.AllSleepAllPahi(1);
    T.FcentralAHI4(i,1) = (AHIdata2.AllSleepAllPApCI(3) + AHIdata2.AllSleepAllPHyCI(3)) / AHIdata2.AllSleepAllPahi(3);
  
    T.FcentralormixedAHI(i,1) = (AHIdata2.AllSleepAllPApCI(1) + AHIdata2.AllSleepAllPMAI(1) + AHIdata2.AllSleepAllPHyCI(1)) / AHIdata2.AllSleepAllPahi(1);
    T.FcentralormixedAHI4(i,1) = (AHIdata2.AllSleepAllPApCI(3) + AHIdata2.AllSleepAllPMAI(3) + AHIdata2.AllSleepAllPHyCI(3)) / AHIdata2.AllSleepAllPahi(3);
       
    T.FN1(i,1) = AHIdata2.N1AllPDur(1)/AHIdata2.AllSleepAllPDur(1);
%     T.FN2(i,1) = AHIdata2.N2AllPDur(1)/AHIdata2.AllSleepAllPDur(1);
%     T.FN3(i,1) = AHIdata2.N3AllPDur(1)/AHIdata2.AllSleepAllPDur(1);
    T.FN2(i,1) = AHIdata2.N2SupDur(1)/AHIdata2.AllSleepAllPDur(1);
    T.FN3(i,1) = AHIdata2.N3SupDur(1)/AHIdata2.AllSleepAllPDur(1);
    T.FREM(i,1) = AHIdata2.RemAllPDur(1)/AHIdata2.AllSleepAllPDur(1);
    
    if sum(strcmp('EvtsData',fieldnames(x)))==1
        EvtsData = x.EvtsData;
        if iscell(EvtsData)==0 && isnan(EvtsData)
            disp(['failed: ' num2str(i)]);
            continue
        end
    else
        disp('EvtsData does not exist, run getAHI2 to fix this')
%         filedir2 = [patients{i,2} patients{i,1} '.mat'];
%         x = matfile(filedir2);
%         Evts = x.Evts;
%         DataEventHypnog_Mat = x.DataEventHypnog_Mat;
%         if ~isfield('RespT',Evts)
%             Evts = EventRespTable(DataEventHypnog_Mat,Evts,ChannelsList);
%         end
%         EvtsData{1} = Evts;
    end
    
    T.HBtotal(i,1) = nanmean(EvtsData{1}.RespT.HBarea)/60*AHIdata2.AllSleepAllPahi(1);
    T.HBtotal4(i,1) = nanmean(EvtsData{1}.RespT.HBarea(EvtsData{1}.RespT.InclAHI4==1))/60*AHIdata2.AllSleepAllPahi(3);
    
    T.EventDurationMean(i,1) = nanmean(EvtsData{1}.RespT.EventDuration);
    T.EventDurationMean4(i,1) = nanmean(EvtsData{1}.RespT.EventDuration(EvtsData{1}.RespT.InclAHI4==1));
    
    T.Farousal(i,1) = nanmean(EvtsData{1}.RespT.ArE);
    T.Farousal4(i,1) = nanmean(EvtsData{1}.RespT.ArE(EvtsData{1}.RespT.InclAHI4==1));
    
    T.DesatMean(i,1) = nanmean(EvtsData{1}.RespT.SpO2DeltaE);
    T.DesatMean4(i,1) = nanmean(EvtsData{1}.RespT.SpO2DeltaE(EvtsData{1}.RespT.InclAHI4==1));
    
    temp = any(EvtsData{1}.RespT.EventCodes == [4 6:10],2);
    temp2 = isnan(EvtsData{1}.RespT.SpO2DeltaE);
    
    T.FhypopneasUnknownDesat(i,1) = sum(temp & temp2)/sum(temp);
    T.FahiUnknownDesat(i,1) = sum(temp & temp2)/height(EvtsData{1}.RespT);
    
    T.AHI4nrem(i,1) = AHIdata2.NRemAllPahi(3);
    T.AHI4rem(i,1) = AHIdata2.RemAllPahi(3);
        if AHIdata2.RemAllPDur(3)<5
            T.AHI4rem(i,1) = NaN;
        end
    T.AHI4remnrembalance(i,1) = (T.AHI4rem(i,1) - T.AHI4nrem(i,1))./(T.AHI4rem(i,1) + T.AHI4nrem(i,1)); %-1 fully NREM, 0 balanced, +1 = fully REM
    
    T.Flateral(i,1) = 1 - AHIdata2.AllSleepSupDur(3)/AHIdata2.AllSleepAllPDur(3);
        
    success(i)=1;
    disp(['Completed getData for study ' num2str(i)]);
    end
end

%%
T = struct2table(T);
T{success(1:height(T))==0,:}=NaN;

%%
savefilename = [settings.workdir 'Summary\getData_' datestr(now,'yyyymmdd HHMMSS')];
savefilename(end-6)='T';
save(savefilename,'T','-v7.3')
savefilename = [settings.workdir 'Summary\getData'];
save(savefilename,'T','-v7.3')

