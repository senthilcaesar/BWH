[~,path,~] = xlsread(AMasterSpreadsheet,1,'A26');
global settings
settings.savename = char('WSC');
settings.Mrange=1:1123;

%%
EDcap=90; % was 60
EDNmin=5;

clear AutoBArTinfoT ManualArTinfoT AutoBArTinfo ManualArTinfo AutoArTinfo AutoArTinfoT
RespTAll=[];
AHIT=[];
EvtDurT=[];

for n=1:1123%settings.Mrange
    
    clear AHIManual AHIAutoR AHIAutoRA EvtDurManual EvtDurManualMed EvtDurManualCOV...
        EvtDurManualnREM EvtDurManualnREMCOV EvtDurAutoR EvtDurAutoRMed  EvtDurAutoRCOV...
        EvtDurAutoRnREM EvtDurAutoRnREMCOV EvtDurAutoRA EvtDurAutoRAMed...
        EvtDurAutoRACOV EvtDurAutoRAnREM EvtDurAutoRAnREMCOV tempDur tempAhi);
    
    n
    RsEvT.Sub(n,1)=n;
    
    T.Sub(n,1)=n;
    
    try
        AnalyzedDir = [path{:} settings.savename '_' num2str(n) '.mat']
        clear Evts
        load (AnalyzedDir,'Evts');
        
        
        %RESP
        try
            % Manual Scoring
            clear tempEvts tempRespT tempEDlist
            tempEvts=Evts;
            tempRespT=tempEvts.RespT;
            
            if 1 %keep all events1
                tempRespT.Subj=repmat(n,height(tempRespT),1);
                if n==1
                    RespTAll=tempRespT;
                else
                    RespTAll=[RespTAll;tempRespT];
                end
            end
            
            try
                RsEvT.AHIManual(n,1) = tempEvts.AHIdata2.AllSleepAllPahi(2); % AHI 3pa
            catch
                RsEvT.AHIManual(n,1)= tempEvts.AHIdata2{1}.AllSleepAllPahi(2);
            end
            
            try
                RsEvT.AHIManual3p(n,1) = tempEvts.AHIdata2.AllSleepAllPahi(4); % AHI 3p
            catch
                RsEvT.AHIManual3p(n,1)= tempEvts.AHIdata2{1}.AllSleepAllPahi(4);
            end
            
            
            try
                if sum(ismember(tempRespT.Properties.VariableNames,'state'))
                    tempRespT.Epochs=tempRespT.state;
                end
                criteria=tempRespT.InclAHI3a==1&tempRespT.Epochs<=3; % 3p criteria now
                tempEDlist=tempRespT.EventDuration(criteria);
            end
            tempEDlist(tempEDlist>EDcap)=EDcap; if length(tempEDlist)<EDNmin, tempEDlist=[]; end
            
            RsEvT.EvtDurManual(n,1)=nanmean(tempEDlist);
            RsEvT.EvtDurManualMed(n,1)=nanmedian(tempEDlist);
            RsEvT.EvtDurManualCOV(n,1)=nanstd(tempEDlist)/nanmean(tempEDlist);
            
            % fraction of hypopneas
            RsEvT.NEvtManual(n,1)=size(tempEDlist,1);
            criteria=tempRespT.InclAHI3a==1&tempRespT.Epochs<=3 & tempRespT.EventCodes==4 ; % 3p criteria now
            tempEDlist=tempRespT.EventDuration(criteria);
            RsEvT.NHypManual(n,1)=size(tempEDlist,1);
            criteria=tempRespT.InclAHI3a==1&tempRespT.Epochs<=3 & tempRespT.EventCodes==2 ;
            tempEDlist=tempRespT.EventDuration(criteria);
            RsEvT.NApneaManual(n,1)=size(tempEDlist,1);
            RsEvT.FHypManual(n,1)=RsEvT.NHypManual(n,1)/ RsEvT.NEvtManual(n,1);
            
            try
                criteria=tempRespT.InclAHI3a==1 & tempRespT.Epochs<3;
                tempEDlist=tempRespT.EventDuration(criteria);
            end
            tempEDlist(tempEDlist>EDcap)=EDcap; if length(tempEDlist)<EDNmin, tempEDlist=[]; end
            
            RsEvT.EvtDurManualnREM(n,1)=nanmean(tempEDlist);
            RsEvT.EvtDurManualnREMMed(n,1)=nanmedian(tempEDlist);
            RsEvT.EvtDurManualnREMCOV(n,1)=nanstd(tempEDlist)/nanmean(tempEDlist);
            
            % fraction of hypopneas
            RsEvT.NEvtManualnREM(n,1)=size(tempEDlist,1);
            
            
            try
                criteria=tempRespT.InclAHI3a==1 & tempRespT.Epochs<3 & tempRespT.EventCodes==4 ; % hypopnea
                tempEDlist=tempRespT.EventDuration(criteria);
            end
            tempEDlist(tempEDlist>EDcap)=EDcap; if length(tempEDlist)<EDNmin, tempEDlist=[]; end
            RsEvT.HypDurManualnREM(n,1)=nanmean(tempEDlist);
            RsEvT.HypDurManualnREMMed(n,1)=nanmedian(tempEDlist);
            RsEvT.NHypManualnREM(n,1)=size(tempEDlist,1);
            RsEvT.FHypManualnREM(n,1)=RsEvT.NHypManualnREM(n,1)/ RsEvT.NEvtManualnREM(n,1);
            
            
            try
                criteria=tempRespT.InclAHI3a==1 & tempRespT.Epochs<3 & tempRespT.EventCodes==2 ; % OSA from manual scoring
                tempEDlist=tempRespT.EventDuration(criteria);
            end
            tempEDlist(tempEDlist>EDcap)=EDcap; if length(tempEDlist)<EDNmin, tempEDlist=[]; end
            RsEvT.ApneaDurManualnREM(n,1)=nanmean(tempEDlist);
            RsEvT.ApneaDurManualnREMMed(n,1)=nanmedian(tempEDlist);
            RsEvT.NApneaManualnREM(n,1)=size(tempEDlist,1);
            
            
            %% AutoScoring - Resp Only
            clear tempEvts tempRespT tempEDlist
            tempEvts=Evts.EvtsAutoRespOnly;
            tempRespT=tempEvts.RespT;
            try
                RsEvT.AHIAutoR(n,1) = tempEvts.AHIdata2.AllSleepAllPahi(2);
            catch
                RsEvT.AHIAutoR(n,1)= tempEvts.AHIdata2{1}.AllSleepAllPahi(2);
            end
            
            try
                RsEvT.AHIAutoR3p(n,1) = tempEvts.AHIdata2.AllSleepAllPahi(4);
            catch
                RsEvT.AHIAutoR3p(n,1)= tempEvts.AHIdata2{1}.AllSleepAllPahi(4);
            end
            
            
            try
                if sum(ismember(tempRespT.Properties.VariableNames,'state'))
                    tempRespT.Epochs=tempRespT.state;
                end
                criteria=tempRespT.InclAHI3a==1&tempRespT.Epochs<=3;
                tempEDlist=tempRespT.EventDuration(criteria);
            end
            tempEDlist(tempEDlist>EDcap)=EDcap; if length(tempEDlist)<EDNmin, tempEDlist=[]; end
            
            RsEvT.EvtDurAutoR(n,1)=nanmean(tempEDlist);
            RsEvT.EvtDurAutoRMed(n,1)=nanmedian(tempEDlist);
            RsEvT.EvtDurAutoRCOV(n,1)=nanstd(tempEDlist)/nanmean(tempEDlist);
            % fraction of hypopneas
            RsEvT.NEvtAutoR(n,1)=size(tempEDlist,1);
            criteria=tempRespT.InclAHI3a==1&tempRespT.Epochs<=3 & tempRespT.EventCodes==4 ; % 3p criteria now
            tempEDlist=tempRespT.EventDuration(criteria);
            RsEvT.NHypAutoR(n,1)=size(tempEDlist,1);
            criteria=tempRespT.InclAHI3a==1&tempRespT.Epochs<=3 & tempRespT.EventCodes==2 ;
            tempEDlist=tempRespT.EventDuration(criteria);
            RsEvT.NApneaAutoR(n,1)=size(tempEDlist,1);
            RsEvT.FHypAutoR(n,1)=RsEvT.NHypAutoR(n,1)/ RsEvT.NEvtAutoR(n,1);
            
            try
                criteria=tempRespT.InclAHI3a==1 & tempRespT.Epochs<3;
                tempEDlist=tempRespT.EventDuration(criteria);
            end
            tempEDlist(tempEDlist>EDcap)=EDcap; if length(tempEDlist)<EDNmin, tempEDlist=[]; end
            
            RsEvT.EvtDurAutoRnREM(n,1)=nanmean(tempEDlist);
            RsEvT.EvtDurAutoRnREMCOV(n,1)=nanstd(tempEDlist)/nanmean(tempEDlist);
            RsEvT.EvtDurAutoRnREMMed(n,1)=nanmedian(tempEDlist);
            
            RsEvT.NEvtAutoRnREM(n,1)=size(tempEDlist,1);
            
            try
                criteria=tempRespT.InclAHI3a==1 & tempRespT.Epochs<3 & tempRespT.EventCodes==4 ; % hypopnea
                tempEDlist=tempRespT.EventDuration(criteria);
            end
            tempEDlist(tempEDlist>EDcap)=EDcap; if length(tempEDlist)<EDNmin, tempEDlist=[]; end
            RsEvT.HypDurAutoRnREM(n,1)=nanmean(tempEDlist);
            RsEvT.HypDurAutoRnREMMed(n,1)=nanmedian(tempEDlist);
            RsEvT.NHypAutoRnREM(n,1)=size(tempEDlist,1);
            RsEvT.FHypAutoRnREM(n,1)=RsEvT.NHypAutoRnREM(n,1)/ RsEvT.NEvtAutoRnREM(n,1);
            
            try
                criteria=tempRespT.InclAHI3a==1 & tempRespT.Epochs<3 & tempRespT.EventCodes==2 ; % apnea auto scoring
                tempEDlist=tempRespT.EventDuration(criteria);
            end
            tempEDlist(tempEDlist>EDcap)=EDcap; if length(tempEDlist)<EDNmin, tempEDlist=[]; end
            RsEvT.ApneaDurAutoRnREM(n,1)=nanmean(tempEDlist);
            RsEvT.ApneaDurAutoRnREMMed(n,1)=nanmedian(tempEDlist);
            RsEvT.NApneaAutoRnREM(n,1)=size(tempEDlist,1);
            
            %% AutoScoring - Auto B
            clear tempEvts tempRespT tempEDlist
            tempEvts=Evts.EvtsAutoB;
            tempRespT=tempEvts.RespT;
            try
                RsEvT.AHIAutoRA(n,1) = tempEvts.AHIdata2.AllSleepAllPahi(2);
            catch
                RsEvT.AHIAutoRA(n,1) = tempEvts.AHIdata2{1}.AllSleepAllPahi(2);
            end
            
            try
                RsEvT.AHIAutoRA3p(n,1) = tempEvts.AHIdata2.AllSleepAllPahi(4);
            catch
                RsEvT.AHIAutoRA3p(n,1) = tempEvts.AHIdata2{1}.AllSleepAllPahi(4);
            end
            
            try
                if sum(ismember(tempRespT.Properties.VariableNames,'state'))
                    tempRespT.Epochs=tempRespT.state;
                end
                criteria=tempRespT.InclAHI3a==1&tempRespT.Epochs<=3;
                tempEDlist=tempRespT.EventDuration(criteria);
            end
            tempEDlist(tempEDlist>EDcap)=EDcap; if length(tempEDlist)<EDNmin, tempEDlist=[]; end
            
            RsEvT.EvtDurAutoRA(n,1)=nanmean(tempEDlist);
            RsEvT.EvtDurAutoRAMed(n,1)=nanmedian(tempEDlist);
            RsEvT.EvtDurAutoRACOV(n,1)=nanstd(tempEDlist)/nanmean(tempEDlist);
            % fraction of hypopneas
            RsEvT.NEvtAutoRA(n,1)=size(tempEDlist,1);
            criteria=tempRespT.InclAHI3a==1&tempRespT.Epochs<=3 & tempRespT.EventCodes==4 ; % 3p criteria now
            tempEDlist=tempRespT.EventDuration(criteria);
            RsEvT.NHypAutoRA(n,1)=size(tempEDlist,1);
            criteria=tempRespT.InclAHI3a==1&tempRespT.Epochs<=3 & tempRespT.EventCodes==2 ;
            tempEDlist=tempRespT.EventDuration(criteria);
            RsEvT.NApneaAutoRA(n,1)=size(tempEDlist,1);
            RsEvT.FHypAutoRA(n,1)=RsEvT.NHypAutoRA(n,1)/ RsEvT.NEvtAutoRA(n,1);
            
            try
                criteria=tempRespT.InclAHI3a==1 & tempRespT.Epochs<3;
                tempEDlist=tempRespT.EventDuration(criteria);
            end
            tempEDlist(tempEDlist>EDcap)=EDcap; if length(tempEDlist)<EDNmin, tempEDlist=[]; end
            
            RsEvT.EvtDurAutoRAnREM(n,1)=nanmean(tempEDlist);
            RsEvT.EvtDurAutoRAnREMCOV(n,1)=nanstd(tempEDlist)/nanmean(tempEDlist);
            RsEvT.EvtDurAutoRAnREMMed(n,1)=nanmedian(tempEDlist);
            
            % fraction of hypopneas
            RsEvT.NEvtAutoRAnREM(n,1)=size(tempEDlist,1);
            
            try
                criteria=tempRespT.InclAHI3a==1 & tempRespT.Epochs<3 & tempRespT.EventCodes==4 ; % hypopnea
                tempEDlist=tempRespT.EventDuration(criteria);
            end
            tempEDlist(tempEDlist>EDcap)=EDcap; if length(tempEDlist)<EDNmin, tempEDlist=[]; end
            RsEvT.HypDurAutoRAnREM(n,1)=nanmean(tempEDlist);
            RsEvT.HypDurAutoRAnREMMed(n,1)=nanmedian(tempEDlist);
            RsEvT.NHypAutoRAnREM(n,1)=size(tempEDlist,1);
            RsEvT.FHypAutoRAnREM(n,1)=RsEvT.NHypAutoRAnREM(n,1)/ RsEvT.NEvtAutoRAnREM(n,1);
            
            try
                criteria=tempRespT.InclAHI3a==1 & tempRespT.Epochs<3 & tempRespT.EventCodes==2 ; % apnea auto scoring
                tempEDlist=tempRespT.EventDuration(criteria);
            end
            tempEDlist(tempEDlist>EDcap)=EDcap; if length(tempEDlist)<EDNmin, tempEDlist=[]; end
            RsEvT.ApneaDurAutoRAnREM(n,1)=nanmean(tempEDlist);
            RsEvT.ApneaDurAutoRAnREMMed(n,1)=nanmedian(tempEDlist);
            
            RsEvT.NApneaAutoRAnREM(n,1)=size(tempEDlist,1);
            
        catch
            RsEvT.AHIManual(n,1)=NaN;
            RsEvT.AHIManual3p(n,1)=NaN;
            RsEvT.EvtDurManual(n,1)=NaN;
            RsEvT.EvtDurManualMed(n,1)=NaN;
            RsEvT.EvtDurManualCOV(n,1)=NaN;
            RsEvT.NEvtManual(n,1)=NaN;
            RsEvT.NHypManual(n,1)=NaN;
            RsEvT.NApneaManual(n,1)=NaN;
            RsEvT.FHypManual(n,1)=NaN;
            RsEvT.EvtDurManualnREM(n,1)=NaN;
            RsEvT.EvtDurManualnREMCOV(n,1)=NaN;
            RsEvT.EvtDurManualnREMMed(n,1)=NaN;
            RsEvT.HypDurManualnREM(n,1)=NaN;
            RsEvT.HypDurManualnREMMed(n,1)=NaN;
            RsEvT.ApneaDurManualnREM(n,1)=NaN;
            RsEvT.ApneaDurManualnREMMed(n,1)=NaN;
            RsEvT.NEvtManualnREM(n,1)=NaN;
            RsEvT.NHypManualnREM(n,1)=NaN;
            RsEvT.FHypManualnREM(n,1)=NaN;
            RsEvT.NApneaManualnREM(n,1)=NaN;
            
            RsEvT.AHIAutoR(n,1)=NaN;
            RsEvT.AHIAutoR3p(n,1)=NaN;
            RsEvT.EvtDurAutoR(n,1)=NaN;
            RsEvT.EvtDurAutoRMed(n,1)=NaN;
            RsEvT.EvtDurAutoRCOV(n,1)=NaN;
            RsEvT.NEvtAutoR(n,1)=NaN;
            RsEvT.NHypAutoR(n,1)=NaN;
            RsEvT.NApneaAutoR(n,1)=NaN;
            RsEvT.FHypAutoR(n,1)=NaN;
            RsEvT.EvtDurAutoRnREM(n,1)=NaN;
            RsEvT.EvtDurAutoRnREMCOV(n,1)=NaN;
            RsEvT.EvtDurAutoRnREMMed(n,1)=NaN;
            RsEvT.HypDurAutoRnREM(n,1)=NaN;
            RsEvT.HypDurAutoRnREMMed(n,1)=NaN;
            RsEvT.ApneaDurAutoRnREM(n,1)=NaN;
            RsEvT.ApneaDurAutoRnREMMed(n,1)=NaN;
            RsEvT.NEvtAutoRnREM(n,1)=NaN;
            RsEvT.NHypAutoRnREM(n,1)=NaN;
            RsEvT.FHypAutoRnREM(n,1)=NaN;
            RsEvT.NApneaAutoRnREM(n,1)=NaN;
            
            RsEvT.AHIAutoRA(n,1)=NaN;
            RsEvT.AHIAutoRA3p(n,1)=NaN;
            RsEvT.EvtDurAutoRA(n,1)=NaN;
            RsEvT.EvtDurAutoRAMed(n,1)=NaN;
            RsEvT.EvtDurAutoRACOV(n,1)=NaN;
            RsEvT.NEvtAutoRA(n,1)=NaN;
            RsEvT.NHypAutoRA(n,1)=NaN;
            RsEvT.NApneaAutoRA(n,1)=NaN;
            RsEvT.FHypAutoRA(n,1)=NaN;
            RsEvT.EvtDurAutoRAnREM(n,1)=NaN;
            RsEvT.EvtDurAutoRAnREMCOV(n,1)=NaN;
            RsEvT.EvtDurAutoRAnREMMed(n,1)=NaN;
            RsEvT.HypDurAutoRAnREM(n,1)=NaN;
            RsEvT.HypDurAutoRAnREMMed(n,1)=NaN;
            RsEvT.ApneaDurAutoRAnREM(n,1)=NaN;
            RsEvT.ApneaDurAutoRAnREMMed(n,1)=NaN;
            RsEvT.NEvtAutoRAnREM(n,1)=NaN;
            RsEvT.NHypAutoRAnREM(n,1)=NaN;
            RsEvT.FHypAutoRAnREM(n,1)=NaN;
            RsEvT.NApneaAutoRAnREM(n,1)=NaN;
        end
        
        %% AROUSALS
        % Manual Scoring
        clear tempEvts ManualArTinfo
        if 1
            try
                tempEvts=Evts;
                ManualArTinfo=struct2table(tempEvts.ArTinfo);
                if n==1
                    ManualArTinfoT=ManualArTinfo;
                else
                    ManualArTinfoT(n,:)=ManualArTinfo;
                end
                
                % Auto B
                clear tempEvts AutoBArTinfo
                tempEvts=Evts.EvtsArAutoB;
                AutoBArTinfo=struct2table(tempEvts.ArTinfo);
                if n==1
                    AutoBArTinfoT=AutoBArTinfo;
                else
                    AutoBArTinfoT(n,:)=AutoBArTinfo;
                end
                
                
                % AutoA
                clear tempEvts AutoArTinfo
                tempEvts=Evts.EvtsArAuto;
                AutoArTinfo=struct2table(tempEvts.ArTinfo);
                AutoArTinfo.Sub=n;
                if n==1
                    AutoArTinfoT=AutoArTinfo;
                else
                    AutoArTinfoT(n,:)=AutoArTinfo;
                end
                
            catch
                ManualArTinfoT{n,:}=nan(1,width(ManualArTinfoT));
                AutoBArTinfoT{n,:}=nan(1,width(AutoBArTinfoT));
                AutoArTinfoT{n,:}=nan(1,width(AutoArTinfoT));
            end
            
            % from getData.m
            try
                clear tempT
                tempT = Evts.EvtsArAutoB.ArT{Evts.EvtsArAutoB.ArT.AASMarousal==1,{'WSBalanceBMax3','ArBalanceBMax3','ArIntOr','ArInt','EventDuration'}};
                T.EEGB_Narousals(n,1) = size(tempT,1);
                
                if T.EEGB_Narousals(n,1)>0
                    T.EEGB_ArInt_Med(n,1) = nanmedian(tempT(:,4));
                    T.EEGB_ArIntOr_Med(n,1) = nanmedian(tempT(:,3));
                    T.EEGB_ArStrength_Med(n,1) = nanmedian(tempT(:,2));
                    T.EEGB_WStrength_Med(n,1) = nanmedian(tempT(:,1));
                    T.EEGB_ArDuration_Med(n,1) = nanmedian(tempT(:,5));
                    T.EEGB_ArInt_Mean(n,1) = nanmean(tempT(:,4));
                    T.EEGB_ArIntOr_Mean(n,1) = nanmean(tempT(:,3));
                    T.EEGB_ArStrength_Mean(n,1) = nanmean(tempT(:,2));
                    T.EEGB_WStrength_Mean(n,1) = nanmean(tempT(:,1));
                    T.EEGB_ArDuration_Mean(n,1) = nanmean(tempT(:,5));
                else
                    T.EEGB_ArInt_Med(n,1) = NaN;
                    T.EEGB_ArIntOr_Med(n,1) = NaN;
                    T.EEGB_ArStrength_Med(n,1) = NaN;
                    T.EEGB_WStrength_Med(n,1) = NaN;
                    T.EEGB_ArDuration_Med(n,1) = NaN;
                    T.EEGB_ArInt_Mean(n,1) = NaN;
                    T.EEGB_ArIntOr_Mean(n,1) = NaN;
                    T.EEGB_ArStrength_Mean(n,1) = NaN;
                    T.EEGB_WStrength_Mean(n,1) = NaN;
                    T.EEGB_ArDuration_Mean(n,1) = NaN;
                end
                clear tempT
                tempT = Evts.EvtsArAuto.ArT{Evts.EvtsArAuto.ArT.AASMarousal==1,{'WSBalanceMax3','ArBalanceMax3','ArIntOr','ArInt','EventDuration'}};
                T.EEGA_Narousals(n,1) = size(tempT,1);
                if T.EEGA_Narousals(n,1)>0
                    T.EEGA_ArInt_Med(n,1) = nanmedian(tempT(:,4));
                    T.EEGA_ArIntOr_Med(n,1) = nanmedian(tempT(:,3));
                    T.EEGA_ArStrength_Med(n,1) = nanmedian(tempT(:,2));
                    T.EEGA_WStrength_Med(n,1) = nanmedian(tempT(:,1));
                    T.EEGA_ArDuration_Med(n,1) = nanmedian(tempT(:,5));
                    T.EEGA_ArInt_Mean(n,1) = nanmean(tempT(:,4));
                    T.EEGA_ArIntOr_Mean(n,1) = nanmean(tempT(:,3));
                    T.EEGA_ArStrength_Mean(n,1) = nanmean(tempT(:,2));
                    T.EEGA_WStrength_Mean(n,1) = nanmean(tempT(:,1));
                    T.EEGA_ArDuration_Mean(n,1) = nanmean(tempT(:,5));
                else
                    T.EEGA_ArInt_Med(n,1) = NaN;
                    T.EEGA_ArIntOr_Med(n,1) = NaN;
                    T.EEGA_ArStrength_Med(n,1) = NaN;
                    T.EEGA_WStrength_Med(n,1) = NaN;
                    T.EEGA_ArDuration_Med(n,1) = NaN;
                    T.EEGA_ArInt_Mean(n,1) = NaN;
                    T.EEGA_ArIntOr_Mean(n,1) = NaN;
                    T.EEGA_ArStrength_Mean(n,1) = NaN;
                    T.EEGA_WStrength_Mean(n,1) = NaN;
                    T.EEGA_ArDuration_Mean(n,1) = NaN;
                end
                
                if T.EEGB_Narousals(n,1)>0
                    clear tempT
                    tempT = Evts.EvtsArAuto.ArT(Evts.EvtsArAuto.ArT.AASMarousal==1,:);
                    T.EEGA_ArI_All(n,1) = height(tempT)./(Evts.EvtsArAuto.AHIdata2.AllSleepAllPDur(2)/60); %only useful if RERAs are scored
                    T.EEGA_ArI_ArStrengthOver50p(n,1) = sum(tempT.ArBalanceMax3>0)./(Evts.EvtsArAuto.AHIdata2.AllSleepAllPDur(2)/60); %only useful if RERAs are scored
                    T.EEGA_ArI_ArStrengthOver75p(n,1) = sum(tempT.ArBalanceMax3>1.098)./(Evts.EvtsArAuto.AHIdata2.AllSleepAllPDur(2)/60); %only useful if RERAs are scored
                    T.EEGA_ArI_ArStrengthOver95p(n,1) = sum(tempT.ArBalanceMax3>2.9444)./(Evts.EvtsArAuto.AHIdata2.AllSleepAllPDur(2)/60);
                    T.EEGA_ArI_ArStrengthOver99p(n,1) = sum(tempT.ArBalanceMax3>4.5951)./(Evts.EvtsArAuto.AHIdata2.AllSleepAllPDur(2)/60);
                    T.EEGA_ArI_ArStrengthOver99p9(n,1) = sum(tempT.ArBalanceMax3>6.9068)./(Evts.EvtsArAuto.AHIdata2.AllSleepAllPDur(2)/60);
                    T.EEGA_ArI_WSStrengthOver50p(n,1) = sum(tempT.WSBalanceMax3>0)./(Evts.EvtsArAuto.AHIdata2.AllSleepAllPDur(2)/60); %only useful if RERAs are scored
                    T.EEGA_ArI_WSStrengthOver75p(n,1) = sum(tempT.WSBalanceMax3>1.098)./(Evts.EvtsArAuto.AHIdata2.AllSleepAllPDur(2)/60); %only useful if RERAs are scored
                    T.EEGA_ArI_WSStrengthOver95p(n,1) = sum(tempT.WSBalanceMax3>2.9444)./(Evts.EvtsArAuto.AHIdata2.AllSleepAllPDur(2)/60);
                    T.EEGA_ArI_WSStrengthOver99p(n,1) = sum(tempT.WSBalanceMax3>4.5951)./(Evts.EvtsArAuto.AHIdata2.AllSleepAllPDur(2)/60);
                    T.EEGA_ArI_WSStrengthOver99p9(n,1) = sum(tempT.WSBalanceMax3>6.9068)./(Evts.EvtsArAuto.AHIdata2.AllSleepAllPDur(2)/60);
                    T.EEGA_ArI_ArIOrOver1(n,1) = sum(tempT.ArIntOr>1)./(Evts.EvtsArAuto.AHIdata2.AllSleepAllPDur(2)/60); %only useful if RERAs are scored
                    T.EEGA_ArI_ArIOrOver2(n,1) = sum(tempT.ArIntOr>2)./(Evts.EvtsArAuto.AHIdata2.AllSleepAllPDur(2)/60); %only useful if RERAs are scored
                    T.EEGA_ArI_ArIOrOver3(n,1) = sum(tempT.ArIntOr>3)./(Evts.EvtsArAuto.AHIdata2.AllSleepAllPDur(2)/60); %only useful if RERAs are scored
                    T.EEGA_ArI_ArIOrOver4(n,1) = sum(tempT.ArIntOr>4)./(Evts.EvtsArAuto.AHIdata2.AllSleepAllPDur(2)/60); %only useful if RERAs are scored
                    T.EEGA_ArI_ArIOrOver5(n,1) = sum(tempT.ArIntOr>5)./(Evts.EvtsArAuto.AHIdata2.AllSleepAllPDur(2)/60); %only useful if RERAs are scored
                    T.EEGA_ArI_ArIOrOver6(n,1) = sum(tempT.ArIntOr>6)./(Evts.EvtsArAuto.AHIdata2.AllSleepAllPDur(2)/60); %only useful if RERAs are scored
                    T.EEGA_ArI_ArIOrOver7(n,1) = sum(tempT.ArIntOr>7)./(Evts.EvtsArAuto.AHIdata2.AllSleepAllPDur(2)/60); %only useful if RERAs are scored
                else
                    T.EEGA_ArI_All(n,1) = NaN;
                    T.EEGA_ArI_ArStrengthOver50p(n,1) = NaN;
                    T.EEGA_ArI_ArStrengthOver75p(n,1) = NaN;
                    T.EEGA_ArI_ArStrengthOver95p(n,1) = NaN;
                    T.EEGA_ArI_ArStrengthOver99p(n,1) = NaN;
                    T.EEGA_ArI_ArStrengthOver99p9(n,1) = NaN;
                    T.EEGA_ArI_WSStrengthOver50p(n,1) = NaN;
                    T.EEGA_ArI_WSStrengthOver75p(n,1) = NaN;
                    T.EEGA_ArI_WSStrengthOver95p(n,1) = NaN;
                    T.EEGA_ArI_WSStrengthOver99p(n,1) = NaN;
                    T.EEGA_ArI_WSStrengthOver99p9(n,1)=NaN;
                    T.EEGA_ArI_ArIOrOver1(n,1) = NaN;
                    T.EEGA_ArI_ArIOrOver2(n,1) = NaN;
                    T.EEGA_ArI_ArIOrOver3(n,1) = NaN;
                    T.EEGA_ArI_ArIOrOver4(n,1) = NaN;
                    T.EEGA_ArI_ArIOrOver5(n,1) = NaN;
                    T.EEGA_ArI_ArIOrOver6(n,1) = NaN;
                    T.EEGA_ArI_ArIOrOver7(n,1) = NaN;
                end
            catch
                T.EEGB_Narousals(n,1) = NaN;
                T.EEGB_ArInt_Med(n,1) = NaN;
                T.EEGB_ArIntOr_Med(n,1) = NaN;
                T.EEGB_ArStrength_Med(n,1) = NaN;
                T.EEGB_WStrength_Med(n,1) = NaN;
                T.EEGB_ArDuration_Med(n,1) = NaN;
                T.EEGB_ArInt_Mean(n,1) = NaN;
                T.EEGB_ArIntOr_Mean(n,1) = NaN;
                T.EEGB_ArStrength_Mean(n,1) = NaN;
                T.EEGB_WStrength_Mean(n,1) = NaN;
                T.EEGB_ArDuration_Mean(n,1) = NaN;
                
                T.EEGA_Narousals(n,1) = NaN;
                T.EEGA_ArInt_Med(n,1) = NaN;
                T.EEGA_ArIntOr_Med(n,1) = NaN;
                T.EEGA_ArStrength_Med(n,1) = NaN;
                T.EEGA_WStrength_Med(n,1) = NaN;
                T.EEGA_ArDuration_Med(n,1) = NaN;
                T.EEGA_ArInt_Mean(n,1) = NaN;
                T.EEGA_ArIntOr_Mean(n,1) = NaN;
                T.EEGA_ArStrength_Mean(n,1) = NaN;
                T.EEGA_WStrength_Mean(n,1) = NaN;
                T.EEGA_ArDuration_Mean(n,1) = NaN;
                
                T.EEGA_ArI_All(n,1) = NaN;
                T.EEGA_ArI_ArStrengthOver50p(n,1) = NaN;
                T.EEGA_ArI_ArStrengthOver75p(n,1) = NaN;
                T.EEGA_ArI_ArStrengthOver95p(n,1) = NaN;
                T.EEGA_ArI_ArStrengthOver99p(n,1) = NaN;
                T.EEGA_ArI_ArStrengthOver99p9(n,1) = NaN;
                T.EEGA_ArI_WSStrengthOver50p(n,1) = NaN;
                T.EEGA_ArI_WSStrengthOver75p(n,1) = NaN;
                T.EEGA_ArI_WSStrengthOver95p(n,1) = NaN;
                T.EEGA_ArI_WSStrengthOver99p(n,1) = NaN;
                T.EEGA_ArI_WSStrengthOver99p9(n,1)=NaN;
                T.EEGA_ArI_ArIOrOver1(n,1) = NaN;
                T.EEGA_ArI_ArIOrOver2(n,1) = NaN;
                T.EEGA_ArI_ArIOrOver3(n,1) = NaN;
                T.EEGA_ArI_ArIOrOver4(n,1) = NaN;
                T.EEGA_ArI_ArIOrOver5(n,1) = NaN;
                T.EEGA_ArI_ArIOrOver6(n,1) = NaN;
                T.EEGA_ArI_ArIOrOver7(n,1) = NaN;
                
            end
        end
        
    catch
        disp(['failed ' num2str(n)]);
        
    end
end

RespEventFinalT=struct2table(RsEvT);
ArousalFinalT=struct2table(T);

ManualArTinfoT.Properties.VariableNames=strcat('Manual_',ManualArTinfoT.Properties.VariableNames);
AutoBArTinfoT.Properties.VariableNames=strcat('AutoB_',AutoBArTinfoT.Properties.VariableNames);
AutoArTinfoT.Properties.VariableNames=strcat('Auto',AutoArTinfoT.Properties.VariableNames);
ArTInfoFinalT=[ManualArTinfoT,AutoBArTinfoT,AutoArTinfoT];
save([settings.workdir 'Summary\AutoEvtsTforCompare_RespAr_3p_HypvAp.mat'],'RespEventFinalT','ArousalFinalT',...
    'ArTInfoFinalT','RespTAll');

% save([settings.workdir 'Summary\AutoEvtsTforCompare3pa.mat'],'AHIT','EvtDurT',...
%   'RespTAll');

%% end of analysis.




