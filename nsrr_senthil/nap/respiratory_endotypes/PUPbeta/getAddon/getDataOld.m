function Data = getDataOld()

%run after SummaryAnalysis
global settings AMasterSpreadsheet ChannelsList

t_startGetData = clock;

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
    try
    filedir = [path{:} settings.savename '_' num2str(i) '.mat'];
    if exist(filedir)==2
        
       disp(['Processing: ' num2str(i) '/' num2str(M) ': ' settings.savename '_' num2str(i) '.mat']);
        
        x = matfile(filedir);

        if sum(strcmp('Evts',fieldnames(x)))==1
                Evts = x.Evts;
                AHIdata2 = Evts.AHIdata2;
        else
            if sum(strcmp('AHIdata2',fieldnames(x)))==1
                AHIdata2 = x.AHIdata2;
            else
                disp(['Evts / AHIdata2 does not exist for patient: ' num2str(i)]);
            continue
            end
        end
        if iscell(AHIdata2)==0 && isnan(AHIdata2)
            disp(['failed: ' num2str(i)]);
            continue
        end
        
        AHIdata2 = AHIdata2{1};
        
        
        T.AHI(i,1) = AHIdata2.AllSleepAllPahi(1);
        %T.AHI3pa(i,1) = AHIdata2.AllSleepAllPahi(2);
        T.AHI4(i,1) = AHIdata2.AllSleepAllPahi(3);
        T.AHI3(i,1) = AHIdata2.AllSleepAllPahi(4);
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
        
        
        T.TST(i,1) = AHIdata2.AllSleepAllPDur(1);
        
       
        
        T.AHISupine3pa(i,1) = AHIdata2.AllSleepSupAHI(1);
        T.AHISupine3(i,1) = AHIdata2.AllSleepSupAHI(3);
        T.AHISupine4(i,1) = AHIdata2.AllSleepSupAHI(4);
        if AHIdata2.AllSleepSupDur<5
            T.AHISupine3pa(i,1) = NaN;
            T.AHISupine3(i,1) = NaN;
            T.AHISupine4(i,1) = NaN;            
        end
        
        %T.ArI
        
        if ~exist('Evts')
            if sum(strcmp('EvtsData',fieldnames(x)))==1
                EvtsData = x.EvtsData;
                Evts = EvtsData{1};
                if iscell(EvtsData)==0 && isnan(EvtsData)
                    disp(['failed: ' num2str(i)]);
                    continue
                end
            end
            
        end
        
        T.HBtotal(i,1) = nanmean(Evts.RespT.HBarea)/60*AHIdata2.AllSleepAllPahi(1);
        T.HBtotal4(i,1) = nanmean(Evts.RespT.HBarea(Evts.RespT.InclAHI4==1))/60*AHIdata2.AllSleepAllPahi(3); %for Luigi
        
        T.EventDurationMean(i,1) = nanmean(Evts.RespT.EventDuration);
        T.EventDurationMean4(i,1) = nanmean(Evts.RespT.EventDuration(Evts.RespT.InclAHI4==1));
        
        T.Farousal(i,1) = nanmean(Evts.RespT.ArE);
        T.Farousal4(i,1) = nanmean(Evts.RespT.ArE(Evts.RespT.InclAHI4==1));
        
        T.DesatMean(i,1) = nanmean(Evts.RespT.SpO2DeltaE);
        T.DesatMean4(i,1) = nanmean(Evts.RespT.SpO2DeltaE(Evts.RespT.InclAHI4==1));
        
        temp = any(Evts.RespT.EventCodes == [4 6:10],2);
        temp2 = isnan(Evts.RespT.SpO2DeltaE);
        
        T.FhypopneasUnknownDesat(i,1) = sum(temp & temp2)/sum(temp);
        T.FahiUnknownDesat(i,1) = sum(temp & temp2)/height(Evts.RespT);
        
        T.AHInrem(i,1) = AHIdata2.NRemAllPahi(1);
        T.AHI4nrem(i,1) = AHIdata2.NRemAllPahi(3);
        
        T.AHIrem(i,1) = AHIdata2.RemAllPahi(1);
        T.AHI4rem(i,1) = AHIdata2.RemAllPahi(3);
        
        T.AHIremnrembalance(i,1) = (T.AHIrem(i,1) - T.AHInrem(i,1))./(T.AHIrem(i,1) + T.AHInrem(i,1)); %-1 fully NREM, 0 balanced, +1 = fully REM
        T.AHI4remnrembalance(i,1) = (T.AHI4rem(i,1) - T.AHI4nrem(i,1))./(T.AHI4rem(i,1) + T.AHI4nrem(i,1)); %-1 fully NREM, 0 balanced, +1 = fully REM
        
        
        if AHIdata2.RemAllPDur(1)<5
            T.AHIrem(i,1) = NaN;
            T.AHI4rem(i,1) = NaN;
            T.AHIremnrembalance(i,1) = NaN;
            T.AHI4remnrembalance(i,1) = NaN;
        end
        
        T.Flateral(i,1) = 1 - AHIdata2.AllSleepSupDur(3)/AHIdata2.AllSleepAllPDur(3);
        
        T.ArI(i,1) = NaN;
        T.T90(i,1) = NaN;
        
        
        %default NaN
        T.ODI3(i,1)=NaN;
        T.ODI4(i,1)=NaN;
        T.SpO2mean(i,1)=NaN;
        T.SpO2meanwake(i,1)=NaN;
        T.SpO2nadir(i,1)=NaN;
            
        if ~exist('Evts') || ~isfield(Evts,'SpO2')
            if sum(strcmp('SpO2data',fieldnames(x)))==1
                SpO2data = x.SpO2data;
                if iscell(SpO2data)==0 && isnan(SpO2data)
                    disp(['failed: ' num2str(i)]);
                end
                SpO2data = SpO2data{1};
            else
                disp('SpO2data does not exist');
            end
        else    
            SpO2data = Evts.SpO2;
        end
        
        try
            T.ODI3(i,1)=SpO2data.ODI3;
            T.ODI4(i,1)=SpO2data.ODI4;
            T.SpO2mean(i,1)=SpO2data.SpO2mean;
            T.SpO2meanwake(i,1)=SpO2data.SpO2meanwake;
            T.SpO2nadir(i,1)=SpO2data.SpO2nadir;
            T.T90(i,1) = SpO2data.SpO2below90p_sleep_prct;
        end
        
        success(i)=1;
        disp(['Completed getData for study: ' settings.savename '_' num2str(i) '.mat']);
    end
    catch 
         disp(['Failed getData for study: ' settings.savename '_' num2str(i) '.mat']);
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

delta_tGetData = etime(clock, t_startGetData); % delta in seconds
D = duration(0,0,round(delta_tGetData),'Format','hh:mm:ss'); % duration in HH:MM:SS
disp(' '); % add row space for visual clarity in command window
displaytext = ['GetData Calculation Complete. Total time: ', char(D), ' (hh:mm:ss)'];
disp(displaytext);

