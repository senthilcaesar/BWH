function [LG1,LGn,delay,VRA,ArThres,LG1N,LGnN,delayN,VRAN,ArThresN,Fstates1,Fsupine1,Fstates2,Fsupine2,LG1values,LGnvalues,delayvalues,VRAvalues,ArThresvalues]=SummaryAnalysisOne_LGmodel(Twin,varlist)
Nperwindowvars=5;
%Twin1 = Twin(I,:);

LG1(1,1)=prctile(Twin.LG1(~isnan(Twin.LG1)),50); % DLM?: is this not the same as: nanmedian(Twin.LG1)  
LG1values=Twin.LG1(~isnan(Twin.LG1));

% LG1_SD(1,1)=std(tempLG1(~isnan(tempLG1)));
%LG1_N(1,1)=sum(~isnan(Twin.LG1)); % no_of_LGmeasures(1,1)=sum(~isnan(tempLG1)&criteria);
% LG1_SE(1,1)=LG1_SD(1,1)/LG1_N(1,1)^0.5;
% LG1_N_nocriteria(1,1)=sum(~isnan(tempLG1)); % no_of_LGmeasures(1,1)=sum(~isnan(tempLG1)&criteria);

% meanSE=mean(LG1_SE);
% F_incl=LG1_N./LG1_N_nocriteria;
% meanFincl=mean(F_incl);
% minLGN=min(LG1_N);

%minmeanmedianNmeasures=[min(no_of_LGmeasures) mean(no_of_LGmeasures) median(no_of_LGmeasures)]
% LGn_=tempLGn(criteria);
% VRA1_=tempVRA1(criteria);
% VRA2_=tempVRA2(criteria);
% MeanEx(1,1)=prctile(MeanE(~isnan(MeanE)),50);
% LG1(1,1)=prctile(tempLG1(criteriaLG),50)
% LG2(1,1)=prctile(tempLG2(~isnan(tempLG2)),50);
LGn(1,1)=prctile(Twin.LGn(~isnan(Twin.LGn)),50); 
LGnvalues=Twin.LGn(~isnan(Twin.LGn));
% Tn(1,1)=prctile(tempTn(~isnan(tempTn)),50);
VRA(1,1)=100*prctile(Twin.VRA(~isnan(Twin.VRA)),50);
VRAvalues=100*Twin.VRA(~isnan(Twin.VRA));
% VRA2(1,1)=100*prctile(tempVRA2(~isnan(tempVRA1)),50);
ArThres(1,1)=100*prctile(Twin.ArThres(~isnan(Twin.ArThres)),50);
ArThresvalues=100*Twin.ArThres(~isnan(Twin.ArThres));
% Rsq_median(1,1)=prctile(Rsq(criteriaLG),50)
% LG0direct(1,1)=prctile(tempLG0(~isnan(tempLG0)),50);
% tau(1,1)=prctile(temptau(~isnan(temptau)),50);
delay(1,1)=prctile(Twin.delay(~isnan(Twin.delay)),50);
delayvalues=Twin.delay(~isnan(Twin.delay));
% LG3min(1,1)=prctile(tempLG3min(~isnan(tempLG3min)),50);
% LG6min(1,1)=prctile(tempLG6min(~isnan(tempLG6min)),50);
% LG90s(1,1)=prctile(tempLG90s(~isnan(tempLG90s)),50);

%get DataN
for ii=1:Nperwindowvars
    eval([varlist{ii} 'N=length(' varlist{ii} 'values);']);
end

Fstates1=[mean(Twin.FNREM1) mean(Twin.FNREM2) mean(Twin.FNREM3) mean(Twin.FREM)]';
Fstates2=[mean(Twin.FNREM1(Twin.Narousals>0)) mean(Twin.FNREM2(Twin.Narousals>0)) mean(Twin.FNREM3(Twin.Narousals>0)) mean(Twin.FREM(Twin.Narousals>0))]';
Fsupine1 = nanmean(Twin.Fsupine);
Fsupine2 = nanmean(Twin.Fsupine(Twin.Narousals>0));
