clear all; clc; close all;
%% BootTAnalysis
% load BootT table
% remove those with AHI<5
%% Options
OSAAHIthres=5;% 
critwin=3;
addpath(genpath('C:\Users\rma56\Dropbox (Partners HealthCare)\PUPbeta_git\PUPbeta\'));

if 1 % mros visit 1
    % Boot data
    load('C:\Users\rma56\Dropbox (Partners HealthCare)\Reliability Trait Analysis\MrOS\Summary\SummaryAnalysis_V1_AllSleep_AllPos_AllSub_Boot.mat')
    load('C:\Users\rma56\Dropbox (Partners HealthCare)\Reliability Trait Analysis\MrOS\Summary\V1AHIBootTable.mat')
        
    % AHI Data
    Nsrr1=readtable('C:\Users\rma56\Dropbox (Partners HealthCare)\Reliability Trait Analysis\MrOS\NSRR Updated AHI\mros1_nsrr_new_ahis_20200518.csv');
    Ahi3pAr=Nsrr1.POOAHI3;
    

else % mesa
    % Boot data
    load('C:\Users\rma56\Dropbox (Partners HealthCare)\Reliability Trait Analysis\MESA\SummaryAnalysis_AllSleep_AllPos_Boot_.mat')
%     load('C:\Users\rma56\Dropbox (Partners HealthCare)\Reliability Trait Analysis\MESA\SummaryAnalysis_AllSleep_AllPos_Boot.mat')
    load('C:\Users\rma56\Dropbox (Partners HealthCare)\Reliability Trait Analysis\MESA\AHIBootTable.mat')
% AHI data
    MESAdata=load('C:\Users\rma56\Dropbox (Partners HealthCare)\Reliability Trait Analysis\MESA\MESAdata.mat');
    Ahi3pAr =MESAdata.MESAdata.a0h3ai5;
    subid=MESAdata.MESAdata.idno;
    % subid from mesa spreadsheet
    excelSubId=readtable('C:\Users\rma56\Dropbox (Partners HealthCare)\Reliability Trait Analysis\MESA\AMasterSpreadsheet.xlsx','Range','T3:T2063');
    Subiddiff=subid-table2array(excelSubId);
    find(Subiddiff>0) % sub# 40 in excelSubId is repetition
    %code to remove sub# 40??
    I = find(Subiddiff>0)';
    % SubRem = sum(BootT.Subject==I,2)>0;
    % temp = SubRem{:,2:9};
    % temp(SubRem,:)=NaN;
    % BootT{:,2:9} = temp;
    
end

% length(unique(BootT.Subject))

clear betaAll betaStdAll
if 1
    OSA=Ahi3pAr>OSAAHIthres;
    
    %code to limit data based on min number of windows
    DataNtemp=BootT{:,10:16}<critwin; % missing ArTh Num in MrOs bootT table
    DataNtemp = [DataNtemp(:,1:4) DataNtemp(:,4) DataNtemp(:,5:7)];
    temp = BootT{:,2:9};
    temp(DataNtemp)=NaN;
    BootT{:,2:9} = temp;
    
    Irem=find(DataNArray(:,6)<critwin);
    BootTPSGParam(any(Irem==BootTPSGParam.Subject'),:)=[]; %Removing no endotype data per critwin
    BootT(any(Irem==BootT.Subject'),:)=[];
    
    %code to remove those with low AHI here
    I = find(~OSA)';
%     NotOSAline = sum(BootT.Subject==I,2)>0;
%     temp = BootT{:,2:9};
%     temp(NotOSAline,:)=NaN;
%     BootT{:,2:9} = temp;
    
    BootTPSGParam(any(I(:)==BootTPSGParam.Subject'),:)=[]; %Removing nonOSA
    BootT(any(I(:)==BootT.Subject'),:)=[]; %Removing nonOSA
    
%     for ii=1:8
%         SubNtemp(:,ii)=DataNtemp(:,ii);%| NotOSAline;
%         SubN{ii} =unique(BootT.Subject(SubNtemp(:,ii)));
%         SubNN(ii)=length(SubN{ii});
%     end
%     length(unique(BootT.Subject))-max(SubNN)
end

% NperTrait = sum(DataNArray>critwin & OSA);
% NPSGdata= length(unique(BootTPSGParam.Subject))

BootT.FlatLoopGain = 1 - BootT.FsupLoopGain;
BootT.FlatArTh = 1 - BootT.FsupArTh;
BootT.FlatUA = 1 - BootT.FsupUA;


BootT.ArThresBootT = fArThresT(BootT.ArThresBoot);
BootT.VpassiveBootT = fVpassiveT(BootT.VpassiveBoot);

% standardized beta coeffs - 
% Betas are calculated by subtracting the mean from the variable and dividing by its standard deviation
BootTStd=BootT;
for ii=2:36
    BootTStd{:,ii}=(BootTStd{:,ii}-nanmean(BootTStd{:,ii}))./nanstd(BootTStd{:,ii});
end

BootTPSGParam.AHIT = BootTPSGParam.AHI.^0.5;
BootTPSGParam.ArIT = BootTPSGParam.ArI.^0.5;

BootTPSGParamStd=BootTPSGParam;
for ii=2:size(BootTPSGParam,2)
    BootTPSGParamStd{:,ii}=(BootTPSGParam{:,ii}-nanmean(BootTPSGParam{:,ii}))./nanstd(BootTPSGParam{:,ii});
end
% BootTStd.ArThresBootT=(BootTStd.ArThresBootT-nanmean(BootTStd.ArThresBootT))./nanstd(BootTStd.ArThresBootT);
% BootTStd.VpassiveBootT=(BootTStd.VpassiveBootT-nanmean(BootTStd.VpassiveBootT))./nanstd(BootTStd.VpassiveBootT);


j=1
mdl = fitglme(BootT,'LG1Boot ~ FN1LoopGain + FN3LoopGain + FREMLoopGain + FlatLoopGain + (1|Subject)')
betaAll(:,j) = mdl.Coefficients.Estimate(2:end)
mdl2 = fitglme(BootTStd,'LG1Boot ~ FN1LoopGain + FN3LoopGain + FREMLoopGain + FlatLoopGain + (1|Subject)')
betaStdAll(:,j) = mdl.Coefficients.Estimate(2:end)./nanstd(BootT.LG1Boot)

j=2
mdl = fitglme(BootT,'LGnBoot ~ FN1LoopGain + FN3LoopGain + FREMLoopGain + FlatLoopGain + (1|Subject)')
betaAll(:,j) = mdl.Coefficients.Estimate(2:end)
mdl2 = fitglme(BootTStd,'LGnBoot ~ FN1LoopGain + FN3LoopGain + FREMLoopGain + FlatLoopGain + (1|Subject)')
betaStdAll(:,j) = mdl.Coefficients.Estimate(2:end)./nanstd(BootT.LGnBoot)

j=3
mdl = fitglme(BootT,'delayBoot ~ FN1LoopGain + FN3LoopGain + FREMLoopGain + FlatLoopGain + (1|Subject)')
betaAll(:,j) = mdl.Coefficients.Estimate(2:end)
mdl2 = fitglme(BootTStd,'delayBoot ~ FN1LoopGain + FN3LoopGain + FREMLoopGain + FlatLoopGain + (1|Subject)')
betaStdAll(:,j) = mdl.Coefficients.Estimate(2:end)./nanstd(BootT.delayBoot)

j=4
mdl = fitglme(BootT,'VRABoot ~ FN1ArTh + FN3ArTh + FREMArTh + FlatArTh + (1|Subject)')
betaAll(:,j) = mdl.Coefficients.Estimate(2:end)
mdl2 = fitglme(BootTStd,'VRABoot ~ FN1ArTh + FN3ArTh + FREMArTh + FlatArTh + (1|Subject)')
betaStdAll(:,j) = mdl.Coefficients.Estimate(2:end)./nanstd(BootT.VRABoot)


j=5
%mdl1 = fitglme(BootT,'ArThresBoot ~ FN1ArTh + FN3ArTh + FREMArTh + FlatArTh + (1|Subject)')
mdl = fitglme(BootT,'ArThresBootT ~ FN1ArTh + FN3ArTh + FREMArTh + FlatArTh + (1|Subject)')
betaAll(:,j) = mdl.Coefficients.Estimate(2:end)
mdl2 = fitglme(BootTStd,'ArThresBootT ~ FN1ArTh + FN3ArTh + FREMArTh + FlatArTh + (1|Subject)')
betaStdAll(:,j) = mdl2.Coefficients.Estimate(2:end)
betaStdAll(:,j) = mdl.Coefficients.Estimate(2:end)./nanstd(BootT.ArThresBootT)


j=6
% mdl1 = fitglme(BootT,'VpassiveBoot ~ FN1UA + FN3UA + FREMUA + FlatUA + (1|Subject)')
mdl = fitglme(BootT,'VpassiveBootT ~ FN1UA + FN3UA + FREMUA + FlatUA + (1|Subject)')
betaAll(:,j) = mdl.Coefficients.Estimate(2:end)
mdl2 = fitglme(BootTStd,'VpassiveBootT ~ FN1UA + FN3UA + FREMUA + FlatUA + (1|Subject)')
betaStdAll(:,j) = mdl.Coefficients.Estimate(2:end)./nanstd(BootT.VpassiveBootT)

% temp = fVpassiveT(SummaryAnalysisTable.Vpassive);
% temp(SummaryAnalysisTableN.Vpassive<critwin|~OSA)=[];
% 
% temp2 = randomEffects(mdl)
% figure(99);
% plot(temp(),temp2,'.')


j=7
mdl = fitglme(BootT,'VactiveBoot ~ FN1UA + FN3UA + FREMUA + FlatUA + (1|Subject)')
betaAll(:,j) = mdl.Coefficients.Estimate(2:end)
mdl2 = fitglme(BootTStd,'VactiveBoot ~ FN1UA + FN3UA + FREMUA + FlatUA + (1|Subject)')
betaStdAll(:,j) = mdl.Coefficients.Estimate(2:end)./nanstd(BootT.VactiveBoot)


j=8
mdl = fitglme(BootT,'VcompBoot ~ FN1UA + FN3UA + FREMUA + FlatUA + (1|Subject)')
betaAll(:,j) = mdl.Coefficients.Estimate(2:end)
mdl2 = fitglme(BootTStd,'VcompBoot ~ FN1UA + FN3UA + FREMUA + FlatUA + (1|Subject)')
betaStdAll(:,j) = mdl.Coefficients.Estimate(2:end)./nanstd(BootT.VcompBoot)

j=9
mdl = fitglme(BootTPSGParam,'AHIT ~ FN1 + FN3 + FREM + FLateral + (1|Subject)')
mdl.Rsquared.Ordinary
betaAll(:,j) = mdl.Coefficients.Estimate(2:end)
mdl2 = fitglme(BootTPSGParamStd,'AHIT ~ FN1 + FN3 + FREM + FLateral + (1|Subject)')
betaStdAll(:,j) = mdl.Coefficients.Estimate(2:end)./nanstd(BootTPSGParam.AHIT)

j=10
mdl = fitglme(BootTPSGParam,'ArIT ~ FN1 + FN3 + FREM + FLateral + (1|Subject)')
mdl.Rsquared.Ordinary
betaAll(:,j) = mdl.Coefficients.Estimate(2:end)
mdl2 = fitglme(BootTPSGParamStd,'ArIT ~ FN1 + FN3 + FREM + FLateral + (1|Subject)')
betaStdAll(:,j) = mdl2.Coefficients.Estimate(2:end)
betaStdAll(:,j) = mdl.Coefficients.Estimate(2:end)./nanstd(BootTPSGParam.ArIT)




betaNames = {'FN1','FN3','FREM','Flateral'};
varlist={'LG1','LGn','Delay','VRA','ArTh','Vpassive','Vactive','Vcomp','AHI','ArI'};
betaAll = array2table(betaAll);
betaAll.Properties.RowNames = betaNames;
betaAll.Properties.VariableNames = varlist;
betaAll.Properties.VariableNames{5}='ArThresT';
betaAll.Properties.VariableNames{6}='VpassiveT';

betaStdAll = array2table(betaStdAll);
betaStdAll.Properties.VariableNames=betaAll.Properties.VariableNames;
betaStdAll.Properties.RowNames = betaNames;


%% Apply to the data
DataT = array2table(DataArray);
DataT.Properties.VariableNames = varlist;
DataTT = DataT;
DataTT.Properties.VariableNames{5}='ArThresT';
DataTT.Properties.VariableNames{6}='VpassiveT';

%to do make VpassiveT etc.
temp = (DataTT.ArThresT/100)-1;
temp(temp<0)=0;
DataTT.ArThresT=100*((temp.^0.5)+1);


temp = 1-(DataTT.VpassiveT/100);
temp(temp<0)=0;
DataTT.VpassiveT=100*(1-(temp.^0.5));

%make corrected traits
DataTTc = DataTT;

optiontest = 1;
%option 1, apply generic correction to both visit 1 and visit 2.
%Compare. Does correlation improve vs no correction?

%option 2, use Fstates/pos from visit 1 (FstatesArray1,FsupineArray1), and apply correction to visit
%2 data. Apply no correction for visit 1 data.
%Compare. Does correlation improve vs no correction?
%Better than option 1?

%generate generic correction based on this table of expected states/positions
expectedT = repmat([0.07;0.05;0.25;0.5]',size(DataT,1),1);


%for corrected data:
for i=1:size(DataT,1)
    for j=1:length(varlist)
        % i=1 %subject
        % j=1 %LG1
        try
            beta = betaAll{:,j};
            actual = [FstatesArray{j}(i,[1 3 4])';1-FsupineArray(i,1)];
            if optiontest==2
                expected = [FstatesArray1{j}(i,[1 3 4])';1-FsupineArray1(i,1)];
            else
                expected = expectedT(i,:);
            end
            delta = actual(:) - expected(:);
            traiteffect = beta.*delta;
            
            % for viewing:
            %         TraitCorrection = table(beta,actual,expected,delta,traiteffect);
            %         TraitCorrection.Properties.RowNames = betaNames;
            
            DataTTc{i,j} = DataTT{i,j} - sum(traiteffect); %error, need to subtract this.
        catch me
        end
    end
end

%% What are mean states etc from bootstrapped data

BootTmean=array2table(nanmean(BootT{:,:}));
BootTmean.Properties.VariableNames=BootT.Properties.VariableNames;

%% mean states etc from original data
clear FstatesArrayMean
%code to remove those with low AHI and window criteria here
FsupineArray2=FsupineArray;
if 1
I = find(~OSA)';
temp = FsupineArray2;
temp(I,:)=NaN;
FsupineArray2(:,1:8) = temp;

DataNtemp=DataN(:,1:8)<critwin; % missing ArTh Num in MrOs bootT table
temp = FsupineArray2;
temp(DataNtemp)=NaN;
FsupineArray2(:,1:8) = temp;
end
FlatMean=nanmean(1-FsupineArray2);

for jj=1:length(varlist)
    FstatesArraytemp = [FstatesArray{1,jj}];
    FstatesArraytemp(DataNtemp(:,jj),:)=NaN;
    FstatesArraytemp(~OSA)=NaN;
    FstatesArrayMean(jj,:)=nanmean(FstatesArraytemp);
end
FstatesArrayMean=array2table(FstatesArrayMean');
FstatesArrayMean.Properties.VariableNames=varlist;
StageNames = {'FN1','FN2','FN3','FREM'};
FstatesArrayMean.Properties.RowNames=StageNames;
