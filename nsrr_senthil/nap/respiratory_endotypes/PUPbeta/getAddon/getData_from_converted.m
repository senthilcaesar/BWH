function Data = getData()

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

%% getconvertedData
[~,~,MasterWorksheet] = xlsread(AMasterSpreadsheet,1,'T4:BF10000');
ConvertMatFlag = cell2mat(MasterWorksheet(:,9));
%last row is based on whether "Convert?" (ConvertMatFlag) has numeric data
lastrow = find(1*(~isnan(ConvertMatFlag)),1,'last');
MasterWorksheet(lastrow+1:end,:)=[];
ConvertMatFlag(lastrow+1:end,:)=[];
Filenames = MasterWorksheet(:,2:8);
%%

MaxM = size(num,1);

try
    M = max(settings.Mrange);
catch
    settings.Mrange=1:MaxM;
    M = max(settings.Mrange);
end

%%
clear x T x2
% T.AHIorig = nan(M,1);
% T.FhypopneasOrig = nan(M,1);
% T.Fhypopneas4 = nan(M,1);
% T.AHI4 = nan(M,1);
% HBtotal = nan(M,1);
success = zeros(M,1);
for i=settings.Mrange
    tempfilename=Filenames{i,1};
    tempfilename(find(tempfilename=='.'):end)=[];
    filedirC = [settings.workdir 'Converted\' tempfilename '_XHz.mat'];
    if exist(filedirC)==2
        try
            clear Evts
            load(filedirC,'Evts');
            if ~exist('Evts')
                disp('no Evts var');
                continue
            end

            if 0 %add this in if we need to pull things from the Analyzed data files, currently we don't; this might add a few seconds
                filedir = [path{:} settings.savename '_' num2str(i) '.mat'];
                if exist(filedir)==2
                    x = matfile(filedir);
                end
            end
            

            AHIdata2 = Evts.AHIdata2{1};
            
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
            
            try %bugfix
                T.FN2(i,1) = AHIdata2.N2AllPDur(1)/AHIdata2.AllSleepAllPDur(1);
                T.FN3(i,1) = AHIdata2.N3AllPDur(1)/AHIdata2.AllSleepAllPDur(1);
            catch
                T.FN2(i,1) = AHIdata2.N2SupAllPDur(1)/AHIdata2.AllSleepAllPDur(1);
                T.FN3(i,1) = AHIdata2.N3SupAllPDur(1)/AHIdata2.AllSleepAllPDur(1);
            end
            
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
            
            T.HBtotal(i,1) = nanmean(Evts.RespT.HBarea)/60*AHIdata2.AllSleepAllPahi(1);
            T.HBtotal4(i,1) = nanmean(Evts.RespT.HBarea(Evts.RespT.InclAHI4==1))/60*AHIdata2.AllSleepAllPahi(3); %for Luigi
            
            try
                I=Evts.RespT.Epochs==3;
            catch
                I=Evts.RespT.state==3; % keeping state here as back-up option.could remove later
            end
            
            T.HBrem(i,1) = nanmean(Evts.RespT.HBarea(I))/60*AHIdata2.RemAllPahi(1);
            T.HBrem4(i,1) = nanmean(Evts.RespT.HBarea(I & Evts.RespT.InclAHI4==1))/60*AHIdata2.RemAllPahi(3);
             try
                I=any(Evts.RespT.Epochs==[0 1 2],2);
            catch
                I=any(Evts.RespT.state==[0 1 2],2);
            end
           
            T.HBnrem(i,1) = nanmean(Evts.RespT.HBarea(I))/60*AHIdata2.NRemAllPahi(1);
            T.HBnrem4(i,1) = nanmean(Evts.RespT.HBarea(I & Evts.RespT.InclAHI4==1))/60*AHIdata2.NRemAllPahi(3);
            
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
            
            
            
            T.Flateral(i,1) = 1 - AHIdata2.AllSleepSupDur(3)/AHIdata2.AllSleepAllPDur(3);
            
            
            T.ArI(i,1) = NaN;
            T.ArInrem(i,1) = NaN;      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            T.ArIrem(i,1) = NaN;      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            T.ArInAASM(i,1) = NaN;
            T.ArInremnAASM(i,1) = NaN;      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            T.ArIremnAASM(i,1) = NaN;      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            try
                I=any(Evts.ArT.Epochs==[0 1 2 3],2) & Evts.ArT.AASMarousal;
                T.ArI(i,1) = sum(I)/(Evts.ArTinfo.TSTmin/60);
                I=any(Evts.ArT.Epochs==[0 1 2],2) & Evts.ArT.AASMarousal;
                T.ArInrem(i,1) = sum(I)/(AHIdata2.NRemAllPDur(1)/60);
                I=any(Evts.ArT.Epochs==[3],2) & Evts.ArT.AASMarousal;
                T.ArIrem(i,1) = sum(I)/(AHIdata2.RemAllPDur(1)/60);
                I=any(Evts.ArT.Epochs==[0 1 2 3],2) & Evts.ArT.nAASMarousal;
                T.ArInAASM(i,1) = sum(I)/(Evts.ArTinfo.TSTmin/60);
                I=any(Evts.ArT.Epochs==[0 1 2],2) & Evts.ArT.nAASMarousal;
                T.ArInremnAASM(i,1) = sum(I)/(AHIdata2.NRemAllPDur(1)/60);
                I=any(Evts.ArT.Epochs==[3],2) & Evts.ArT.nAASMarousal;
                T.ArIremnAASM(i,1) = sum(I)/(AHIdata2.RemAllPDur(1)/60);
            end
            
            if AHIdata2.RemAllPDur(1)<5 %remove all REM metrics if not at least 5 min of REM data available
                T.AHIrem(i,1) = NaN;
                T.AHI4rem(i,1) = NaN;
                T.HBrem(i,1) = NaN;
                T.HBrem4(i,1) = NaN;
                T.AHIremnrembalance(i,1) = NaN;
                T.AHI4remnrembalance(i,1) = NaN;
                T.ArIrem(i,1) = NaN;
                T.ArIremnAASM(i,1) = NaN;
            end
            
            
            %default NaN
            T.ODI3(i,1)=NaN;
            T.ODI4(i,1)=NaN;
            T.SpO2mean(i,1)=NaN;
            T.SpO2meanwake(i,1)=NaN;
            T.SpO2nadir(i,1)=NaN;
            T.T90(i,1) = NaN;
            try
                T.ODI3(i,1)=Evts.SpO2.ODI3;
                T.ODI4(i,1)=Evts.SpO2.ODI4;
                T.SpO2mean(i,1)=Evts.SpO2.SpO2mean;
                T.SpO2meanwake(i,1)=Evts.SpO2.SpO2meanwake;
                T.SpO2nadir(i,1)=Evts.SpO2.SpO2nadir;
                T.T90(i,1) = Evts.SpO2.SpO2below90p_sleep_prct;
            catch
            end
            
            success(i)=1;
            disp(['Completed getData for study: ' settings.savename '_' num2str(i) '.mat']);
        catch
            disp(['Failed getData for study: ' settings.savename '_' num2str(i) '.mat']);
        end
        
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

