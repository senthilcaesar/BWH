function T = getData(MrangeOverride)
% RUN StartHere.m first

% changing to run from converted files-1/15/2021


global settings AMasterSpreadsheet ChannelsList

t_startGetData = clock;

%%

[~,path,~] = xlsread(AMasterSpreadsheet,1,'A26'); % only used incase of analyzed file

[~,~,MasterWorksheet] = xlsread(AMasterSpreadsheet,1,'T4:BF10000');
lastrow = find(1*(~isnan(cell2mat(MasterWorksheet(:,9)))),1,'last');

MasterWorksheet(lastrow+1:end,:)=[];

ConvertList = cell2mat(MasterWorksheet(:,9));

NaNlist = (isnan(ConvertList));
ConvertList(NaNlist) = 0; % set missing or NaN as 0


Filenames = MasterWorksheet(:,2:8);


%%

MaxM = size(ConvertList,1); % right now using all the converted files.

try
    settings.Mrange=1:MaxM;
    M = max(settings.Mrange);
end

%%
if exist('MrangeOverride')
    settings.Mrange = MrangeOverride;
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
    disp(['Processing: ' num2str(i) '/' num2str(settings.Mrange(end))]);

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
            
                        
            loadEvtsFromAnalysis=0;
            if isfield(settings,'useEvtsAutoRespOnlyInGetData') && settings.useEvtsAutoRespOnlyInGetData==1
                loadEvtsFromAnalysis=1;
            end
            if isfield(settings,'UseAutoRespInEventAnalysis') && any(settings.UseAutoRespInEventAnalysis==[1 2 3])
                loadEvtsFromAnalysis=1;
            end
            if loadEvtsFromAnalysis
                filedir = [path{:} settings.savename '_' num2str(i) '.mat'];
                if exist(filedir)==2
                   load(filedir,'Evts');
                end
            end
            
            EvtsTemp = Evts; %Evts is default source of RespT ArT and AHIdata2
            if isfield(settings,'useEvtsAutoRespOnlyInGetData') && settings.useEvtsAutoRespOnlyInGetData==1
                EvtsTemp = Evts.EvtsAutoRespOnly; %must come from analysis not convert
            end
            if isfield(settings,'UseAutoRespInEventAnalysis') && settings.UseAutoRespInEventAnalysis>0 %this setting overrules useEvtsAutoRespOnlyInGetData
                switch settings.UseAutoRespInEventAnalysis
                    case 1
                        EvtsTemp = Evts.EvtsAutoRespOnly; %must come from analysis not convert
                    case 2
                        EvtsTemp = Evts.EvtsAuto; % autoarousal 1 with autoresp. must come from analysis not convert
                    case 3
                        EvtsTemp = Evts.EvtsAutoB; %autoarousal 2 with autoresp. must come from analysis not convert
                    case 4
                        EvtsTemp = Evts.EvtsArAuto; %ok from convert
                    case 5
                        EvtsTemp = Evts.EvtsArAutoB; %ok from convert
                end
            end
            EvtsUse = struct();
            EvtsUse.RespT=EvtsTemp.RespT; EvtsUse.ArT=EvtsTemp.ArT; EvtsUse.AHIdata2=EvtsTemp.AHIdata2;
            clear EvtsTemp;
            
            AHIdata2 = EvtsUse.AHIdata2;
            if ~istable(AHIdata2) %legacy code fix
                AHIdata2 = AHIdata2{1};
            end
            
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
            
            % trying quick fix for prev. created Evts.RespT file which has
            % "state" instead of  "Epochs"---RMA 3/24/21
            if sum(strcmp('state',EvtsUse.RespT.Properties.VariableNames))==1
                EvtsUse.RespT.Epochs=EvtsUse.RespT.state;
            end
            
            
            Iall=EvtsUse.RespT.Epochs<=3;
            if length(Iall)>0 
                T.HBtotal(i,1) = nanmean(EvtsUse.RespT.HBarea(Iall))/60*AHIdata2.AllSleepAllPahi(1);
            else
                T.HBtotal(i,1) = 0;
            end
            Icrit = EvtsUse.RespT.InclAHI4==1 & Iall;
            if length(Icrit)>0 
                T.HBtotal4(i,1) = nanmean(EvtsUse.RespT.HBarea(Icrit))/60*AHIdata2.AllSleepAllPahi(3); %for Luigi
            else
                T.HBtotal4(i,1) = 0;
            end

            Irem=EvtsUse.RespT.Epochs==3;
            if length(Irem)>0 
                T.HBrem(i,1) = nanmean(EvtsUse.RespT.HBarea(Irem))/60*AHIdata2.RemAllPahi(1);
            else 
                T.HBrem(i,1)  = 0; %later overwritten with NaN if insufficient REM duration
            end
            Icrit=EvtsUse.RespT.Epochs==3 & EvtsUse.RespT.InclAHI4==1;
            if length(Icrit)>0 
                T.HBrem4(i,1) = nanmean(EvtsUse.RespT.HBarea(Icrit))/60*AHIdata2.RemAllPahi(3);
            else 
                T.HBrem4(i,1) = 0;
            end
            
            Isup=EvtsUse.RespT.Position==1 & EvtsUse.RespT.Epochs<=3;
            if length(Isup)>0 
                T.HBsup(i,1) = nanmean(EvtsUse.RespT.HBarea(Isup))/60*AHIdata2.AllSleepSupAHI(1);
            else 
                T.HBsup(i,1)  = 0; %later overwritten with NaN if insufficient REM duration
            end
            Icrit=EvtsUse.RespT.Position==1 & EvtsUse.RespT.Epochs<=3 & EvtsUse.RespT.InclAHI4==1;
            if length(Icrit)>0 
                T.HBsup4(i,1) = nanmean(EvtsUse.RespT.HBarea(Icrit))/60*AHIdata2.AllSleepSupAHI(3);
            else 
                T.HBsup4(i,1) = 0;
            end
            
            T.HBnrem(i,1)  = 0;  %default
            T.HBnrem4(i,1) = 0;  %default
            if height(EvtsUse.RespT)>0
                Inrem=any(EvtsUse.RespT.Epochs==[0 1 2],2);
                if length(Inrem)>0 
                    T.HBnrem(i,1) = nanmean(EvtsUse.RespT.HBarea(Inrem))/60*AHIdata2.NRemAllPahi(1);
                end
                Icrit=any(EvtsUse.RespT.Epochs==[0 1 2],2) & EvtsUse.RespT.InclAHI4==1;
                if length(Icrit)>0 
                    T.HBnrem4(i,1) = nanmean(EvtsUse.RespT.HBarea(Icrit))/60*AHIdata2.NRemAllPahi(3);
                end
            else
                Inrem=[];
            end
            
            T.EventDurationMean(i,1) = nanmean(EvtsUse.RespT.EventDuration);
            T.EventDurationMean4(i,1) = nanmean(EvtsUse.RespT.EventDuration(EvtsUse.RespT.InclAHI4==1));
            T.EventDurationMeannrem(i,1) = nanmean(EvtsUse.RespT.EventDuration(Inrem));
            T.EventDurationMean4nrem(i,1) = nanmean(EvtsUse.RespT.EventDuration(Inrem & EvtsUse.RespT.InclAHI4==1));
            T.EventDurationMeanrem(i,1) = nanmean(EvtsUse.RespT.EventDuration(Irem));
            T.EventDurationMean4rem(i,1) = nanmean(EvtsUse.RespT.EventDuration(Irem & EvtsUse.RespT.InclAHI4==1)); 
            
            T.Farousal(i,1) = nanmean(EvtsUse.RespT.ArE);
            T.Farousal4(i,1) = nanmean(EvtsUse.RespT.ArE(EvtsUse.RespT.InclAHI4==1));
            
            T.DesatMean(i,1) = nanmean(EvtsUse.RespT.SpO2DeltaE);
            T.DesatMean4(i,1) = nanmean(EvtsUse.RespT.SpO2DeltaE(EvtsUse.RespT.InclAHI4==1));
            
            T.DesatMeanBaseline(i,1) = nanmean(EvtsUse.RespT.SpO2BaselineE);
            T.DesatMeanBaseline4(i,1) = nanmean(EvtsUse.RespT.SpO2BaselineE(EvtsUse.RespT.InclAHI4==1));
            
            temp = any(EvtsUse.RespT.EventCodes == [4 6:10],2);
            temp2 = isnan(EvtsUse.RespT.SpO2DeltaE);
            
            T.FhypopneasUnknownDesat(i,1) = sum(temp & temp2)/sum(temp);
            T.FahiUnknownDesat(i,1) = sum(temp & temp2)/height(EvtsUse.RespT);
            
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
                I=any(EvtsUse.ArT.Epochs==[0 1 2 3],2) & EvtsUse.ArT.AASMarousal;
                T.ArI(i,1) = sum(I)/(T.TST(i)/60);
                I=any(EvtsUse.ArT.Epochs==[0 1 2],2) & EvtsUse.ArT.AASMarousal;
                T.ArInrem(i,1) = sum(I)/(AHIdata2.NRemAllPDur(1)/60);
                I=any(EvtsUse.ArT.Epochs==[3],2) & EvtsUse.ArT.AASMarousal;
                T.ArIrem(i,1) = sum(I)/(AHIdata2.RemAllPDur(1)/60);
                I=any(EvtsUse.ArT.Epochs==[0 1 2 3],2) & EvtsUse.ArT.nAASMarousal;
                T.ArInAASM(i,1) = sum(I)/(T.TST(i)/60);
                I=any(EvtsUse.ArT.Epochs==[0 1 2],2) & EvtsUse.ArT.nAASMarousal;
                T.ArInremnAASM(i,1) = sum(I)/(AHIdata2.NRemAllPDur(1)/60);
                I=any(EvtsUse.ArT.Epochs==[3],2) & EvtsUse.ArT.nAASMarousal;
                T.ArIremnAASM(i,1) = sum(I)/(AHIdata2.RemAllPDur(1)/60);
            catch
            end
            
            try
                if height(EvtsUse.ArT)>0
                    tempT = EvtsUse.ArT;
                    tempT = tempT(tempT.AASMarousal==1,:);
                                   
                    missingArIntOr=0; missingArInt=0;
                    %avoid breaking analysis if ArIntensity is not run
                    if ~any(tempT.Properties.VariableNames=="ArIntOr")
                        tempT.ArIntOr = nan(height(tempT),1);
                        missingArIntOr=1; 
                    end
                    if ~any(tempT.Properties.VariableNames=="ArInt")
                        tempT.ArInt = nan(height(tempT),1);
                        missingArInt=1;
                    end
                
                    tempT2 = tempT{:,{'WSBalanceMax3','ArBalanceMax3','ArIntOr','ArInt','EventDuration'}};
                    
                    T.ArInt_Med(i,1) = nanmedian(tempT2(:,4));
                    T.ArIntOr_Med(i,1) = nanmedian(tempT2(:,3));
                    T.ArStrength_Med(i,1) = nanmedian(tempT2(:,2));
                    T.WStrength_Med(i,1) = nanmedian(tempT2(:,1));
                    T.ArDuration_Med(i,1) = nanmedian(tempT2(:,5));
                    T.ArInt_Mean(i,1) = nanmean(tempT2(:,4));
                    T.ArIntOr_Mean(i,1) = nanmean(tempT2(:,3));
                    T.ArStrength_Mean(i,1) = nanmean(tempT2(:,2));
                    T.WStrength_Mean(i,1) = nanmean(tempT2(:,1));
                    T.ArDuration_Mean(i,1) = nanmean(tempT2(:,5));
                    
                    T.ArI_All(i,1) = height(tempT)./(AHIdata2.AllSleepAllPDur_(1)/60);
                    T.ArI_ArStrengthOver50p(i,1) = sum(tempT.ArBalanceMax3>0)./(AHIdata2.AllSleepAllPDur_(1)/60); 
                    T.ArI_ArStrengthOver75p(i,1) = sum(tempT.ArBalanceMax3>1.098)./(AHIdata2.AllSleepAllPDur_(1)/60); 
                    T.ArI_ArStrengthOver95p(i,1) = sum(tempT.ArBalanceMax3>2.9444)./(AHIdata2.AllSleepAllPDur_(1)/60);
                    T.ArI_ArStrengthOver99p(i,1) = sum(tempT.ArBalanceMax3>4.5951)./(AHIdata2.AllSleepAllPDur_(1)/60);
                    T.ArI_ArStrengthOver99p9(i,1) = sum(tempT.ArBalanceMax3>6.9068)./(AHIdata2.AllSleepAllPDur_(1)/60);
                    if ~missingArIntOr
                        T.ArI_ArIOrOver1(i,1) = sum(tempT.ArIntOr>1)./(AHIdata2.AllSleepAllPDur_(1)/60);
                        T.ArI_ArIOrOver2(i,1) = sum(tempT.ArIntOr>2)./(AHIdata2.AllSleepAllPDur_(1)/60);
                        T.ArI_ArIOrOver3(i,1) = sum(tempT.ArIntOr>3)./(AHIdata2.AllSleepAllPDur_(1)/60);
                        T.ArI_ArIOrOver4(i,1) = sum(tempT.ArIntOr>4)./(AHIdata2.AllSleepAllPDur_(1)/60);
                        T.ArI_ArIOrOver5(i,1) = sum(tempT.ArIntOr>5)./(AHIdata2.AllSleepAllPDur_(1)/60);
                        T.ArI_ArIOrOver6(i,1) = sum(tempT.ArIntOr>6)./(AHIdata2.AllSleepAllPDur_(1)/60);
                        T.ArI_ArIOrOver7(i,1) = sum(tempT.ArIntOr>7)./(AHIdata2.AllSleepAllPDur_(1)/60);
                    else
                        T.ArI_ArIOrOver1(i,1) = NaN;
                        T.ArI_ArIOrOver2(i,1) = NaN;
                        T.ArI_ArIOrOver3(i,1) = NaN;
                        T.ArI_ArIOrOver4(i,1) = NaN;
                        T.ArI_ArIOrOver5(i,1) = NaN;
                        T.ArI_ArIOrOver6(i,1) = NaN;
                        T.ArI_ArIOrOver7(i,1) = NaN;
                    end
                    
                else
                    error('');
                end
            catch
                T.ArInt_Med(i,1) = NaN;
                T.ArIntOr_Med(i,1) = NaN;
                T.ArStrength_Med(i,1) = NaN;
                T.WStrength_Med(i,1) = NaN;
                T.ArDuration_Med(i,1) = NaN;
                T.ArInt_Mean(i,1) = NaN;
                T.ArIntOr_Mean(i,1) = NaN;
                T.ArStrength_Mean(i,1) = NaN;
                T.WStrength_Mean(i,1) = NaN;
                T.ArDuration_Mean(i,1) = NaN;
                
                T.ArI_All(i,1) = NaN;
                T.ArI_ArStrengthOver50p(i,1) = NaN;
                T.ArI_ArStrengthOver75p(i,1) = NaN;
                T.ArI_ArStrengthOver95p(i,1) = NaN;
                T.ArI_ArStrengthOver99p(i,1) = NaN;
                T.ArI_ArStrengthOver99p9(i,1) = NaN;
                T.ArI_ArIOrOver1(i,1) = NaN;
                T.ArI_ArIOrOver2(i,1) = NaN;
                T.ArI_ArIOrOver3(i,1) = NaN;
                T.ArI_ArIOrOver4(i,1) = NaN;
                T.ArI_ArIOrOver5(i,1) = NaN;
                T.ArI_ArIOrOver6(i,1) = NaN;
                T.ArI_ArIOrOver7(i,1) = NaN;
            end
            
            try
                load(filedirC,'Info');
                T.WSacc(i,1) = max(Info.WSinfo.Acc);
                T.WSaccPred(i,1) = max(Info.WSinfo.AccPred);
            catch
                T.WSacc(i,1) = NaN;
                T.WSaccPred(i,1) = NaN;
            end
            
            try
                T.ReraIndex(i,1) = sum(EvtsUse.RespT.EventCodes==12 & EvtsUse.RespT.Epochs<=3)./(AHIdata2.AllSleepAllPDur_(1)/60); %only useful if RERAs are scored
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
                T.EventDurationMeanrem(i,1) = NaN;
                T.EventDurationMean4rem(i,1) = NaN;
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
            
            try
               if isfield(Evts,'EpochT')
                    temp = Evts.EpochT;
                    Ionset=find(temp.Epochs<4,1,'first'); %points to first sleep
                    Ifirstscoredwake=find(temp.Epochs==4,1,'first'); %points to first wake
                    T.SleepOnsetTime(i,1)=temp.EpochStart(Ionset);
                    T.FirstScoredWake(i,1)=temp.EpochStart(Ifirstscoredwake);
                    Ioffset=find(temp.Epochs<4,1,'last'); %points to last sleep
                    temp2 = temp.Epochs(Ionset:Ioffset);
                    T.WASO(i,1) = sum(temp2==4)/2;
                    temp3 = (temp2==4)*1;
                        temp3(temp2==8) = NaN;
                    T.Nawakenings(i,1) = sum(diff(temp3)==1); %Within sleep; doesn't count final awakening
                    % to get a pseudo sleep onset latency in min use: (T.SleepOnsetTime-T.FirstScoredWake)/60
               else
                   error('')
               end
            catch
                T.SleepOnsetTime(i,1)=NaN;
                T.FirstScoredWake(i,1)=NaN;
                T.WASO(i,1)=NaN;
                T.Nawakenings(i,1)=NaN;
            end
            
            success(i)=1;
            disp(['Completed getData for study: ' tempfilename '_XHz.mat']);
        catch
            disp(['Failed getData for study: ' tempfilename '_XHz.mat']);
        end
        
    end
end

%%
T = struct2table(T);
T{success(1:height(T))==0,:}=NaN;

%% QC edits
minHB=0.01; %used for AHI>0 but AHI<1 and HB=NaN. Repeated for all permutations. 

I = T.AHI==0;
T.HBtotal(I)=0;
I = isnan(T.HBtotal) & T.AHI<1;
T.HBtotal(I)=minHB;

I = T.AHI4==0;
T.HBtotal4(I)=0;
I = isnan(T.HBtotal4) & T.AHI4<1;
T.HBtotal4(I)=minHB;

I = T.AHIrem==0; 
T.HBrem(I)=0;
I = isnan(T.HBrem) & T.AHIrem<1;
T.HBrem(I)=minHB;

I = T.AHI4rem==0; 
T.HBrem4(I)=0;
I = isnan(T.HBrem4) & T.AHI4rem<1;
T.HBrem4(I)=minHB;

I = T.AHISupine3pa==0; 
T.HBsup(I)=0;
I = isnan(T.HBsup) & T.AHISupine3pa<1;
T.HBsup(I)=minHB;

I = T.AHISupine4==0; 
T.HBsup4(I)=0;
I = isnan(T.HBsup4) & T.AHISupine4<1;
T.HBsup4(I)=minHB;

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

