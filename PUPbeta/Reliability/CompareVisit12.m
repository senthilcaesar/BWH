clear all; close all; clc;

addpath(genpath('C:\Users\rma56\Dropbox (Partners HealthCare)\PUPbeta_git\'));
addpath(genpath('C:\Users\rma56\Dropbox (Partners HealthCare)\Reliability Trait Analysis\MrOS\Summary'));
addpath(genpath('C:\Users\rma56\Dropbox (Partners HealthCare)\Reliability Trait Analysis\MrOS\NSRR Table'));
%% ALL SLEEP TABLE
file1='SummaryAnalysis_V1_AllSleep_AllPos_Boot.mat';
file2='SummaryAnalysis_V2_AllSleep_AllPosUpdated_Boot.mat';
Sum1=load(file1);
Sum2=load(file2);
l1=size(Sum1.SummaryAnalysisTable,1);
l2=size(Sum2.SummaryAnalysisTable,1);
ldel=l2-l1+1; % delete the values corresponding to first visit, n=2096;
SummaryAnalysisTable1_=Sum1.SummaryAnalysisTable;
SummaryAnalysisTable2_=Sum2.SummaryAnalysisTable(ldel:end,:);

SummaryAnalysisTable1=SummaryAnalysisTable1_;
SummaryAnalysisTable2=SummaryAnalysisTable2_;

% transform ArTh
SummaryAnalysisTable1.ArThres=fArThresT(SummaryAnalysisTable1_.ArThres);
SummaryAnalysisTable2.ArThres=fArThresT(SummaryAnalysisTable2_.ArThres);


% transform Vpassive
SummaryAnalysisTable1.Vpassive=fVpassiveT(SummaryAnalysisTable1_.Vpassive);
SummaryAnalysisTable2.Vpassive=fVpassiveT(SummaryAnalysisTable2_.Vpassive);  

% to do: remove pts with less than 3 windows
Nwind1=Sum1.SummaryAnalysisTableN;
Nwind2=Sum2.SummaryAnalysisTableN;


SummaryAnalysisTableN1=Nwind1;
SummaryAnalysisTableN2=Nwind2;

for i=1:8
NwinAvg1(i,1)=nanmean(SummaryAnalysisTableN1.(i));
NwinAvg2(i,1)=nanmean(SummaryAnalysisTableN2.(i));
end

% NREM EXCLUSIVE TABLE
file1nr='SummaryAnalysis_V1_NREM_AllPos.mat';
file2nr='SummaryAnalysis_V2_NREM_AllPosUpdated.mat';
Sum1nr=load(file1nr,'SummaryAnalysisTable');
Sum2nr=load(file2nr,'SummaryAnalysisTable');
l1nr=size(Sum1nr.SummaryAnalysisTable,1);
l2nr=size(Sum2nr.SummaryAnalysisTable,1);
ldelnr=l2nr-l1nr+1; % delete the values corresponding to first visit, n=2096;
SummaryAnalysisTablenr1=Sum1nr.SummaryAnalysisTable;
SummaryAnalysisTablenr2=Sum2nr.SummaryAnalysisTable(ldelnr:end,:);

% transform ArTh
SummaryAnalysisTablenr1.ArThres=fArThresT(SummaryAnalysisTablenr1.ArThres);
SummaryAnalysisTablenr2.ArThres=fArThresT(SummaryAnalysisTablenr2.ArThres);


% transform Vpassive
SummaryAnalysisTablenr1.Vpassive=fVpassiveT(SummaryAnalysisTablenr1.Vpassive);
SummaryAnalysisTablenr2.Vpassive=fVpassiveT(SummaryAnalysisTablenr2.Vpassive);  

% POSITION/STATES
Fsupine_1=Sum1.Fsupine;
Fsupine_2=Sum2.Fsupine;
Fsupine_2=Fsupine_2(ldel:end,:);
Flateral1= 1 - Fsupine_1;
Flateral2=1 - Fsupine_2;

FstatesArray_1=Sum1.FstatesArray;
FstatesArray_2=Sum2.FstatesArray;

for i=1:size(FstatesArray_2,2)
    FstatesArray_2{i}(1:ldel-1,:)=[];
end

for i=1:size(FstatesArray_1,2)
    FstatesArray_1{i}=AddRowsofNaNs(FstatesArray_1{i},size(FstatesArray_1{1,1},1));
    FstatesArray_2{i}=AddRowsofNaNs(FstatesArray_2{i},size(FstatesArray_2{1,1},1));
end

%NUMBER OF WINDOWS IN ALL SLEEP & NREM
SummaryAnalysisTableN_1=Sum1.SummaryAnalysisTableN;
SummaryAnalysisTableN_2=Sum2.SummaryAnalysisTableN;
SummaryAnalysisTableN_2=SummaryAnalysisTableN_2(ldel:end,:);
SummaryAnalysisTable1=AddRowsofNaNs(SummaryAnalysisTable1,l1);
SummaryAnalysisTable2=AddRowsofNaNs(SummaryAnalysisTable2,l1);
SummaryAnalysisTableN_1=AddRowsofZeros(SummaryAnalysisTableN_1,l1);
SummaryAnalysisTableN_2=AddRowsofZeros(SummaryAnalysisTableN_2,l1);


SummaryAnalysisTableNnr_1=load(file1nr,'SummaryAnalysisTableN');
SummaryAnalysisTableNnr_2=load(file2nr,'SummaryAnalysisTableN');
SummaryAnalysisTableNnr_1=SummaryAnalysisTableNnr_1.SummaryAnalysisTableN;
SummaryAnalysisTableNnr_2=SummaryAnalysisTableNnr_2.SummaryAnalysisTableN(ldel:end,:);
SummaryAnalysisTablenr1=AddRowsofNaNs(SummaryAnalysisTablenr1,l1);
SummaryAnalysisTablenr2=AddRowsofNaNs(SummaryAnalysisTablenr2,l1);
SummaryAnalysisTableNnr_1=AddRowsofZeros(SummaryAnalysisTableNnr_1,l1);
SummaryAnalysisTableNnr_2=AddRowsofZeros(SummaryAnalysisTableNnr_2,l1);

AHI1=Sum1.AHItotal;
AHI2=Sum2.AHItotal;
AHI1=AHI1(:);
AHI2=AHI2(ldel:end);
AHI2=AHI2(:);

AHIData1=Sum1.AHIdata;
AHIData2=Sum2.AHIdata;
AHIData2=AHIData2(ldel:end,:);

for jj=1:length(AHIData1)
    try
        TST1(jj,1)=AHIData1{jj}(57);
    catch
        TST1(jj,1)=NaN;
    end
    try
        TST2(jj,1)=AHIData2{jj}(57);
    catch
        TST2(jj,1)=NaN;
    end
end
% TST1 =57
%58-ahi total

%% nsrr table for visit 1 and 2; get AHI from here
load('C:\Users\rma56\Dropbox (Partners HealthCare)\Reliability Trait Analysis\MrOS\NSRR Table\mros_visit1-2.mat');

%AGE
AgeV1=T1_new.vsage1; 
AgeV2=T2_new.vs2age1;

%BMI 
BMIV1=str2double(T1_new.hwbmi);
BMIV2=str2double(T2_new.hwbmi);

% ahi from new updated nsrr tables
Nsrr1=readtable('C:\Users\rma56\Dropbox (Partners HealthCare)\Reliability Trait Analysis\MrOS\NSRR Updated AHI\mros1_nsrr_new_ahis_20200601.csv');
Nsrr2=readtable('C:\Users\rma56\Dropbox (Partners HealthCare)\Reliability Trait Analysis\MrOS\NSRR Updated AHI\mros2_nsrr_new_ahis_20200601.csv');

for ii=1:size(Nsrr1,1)
   indx=find(strcmp(Nsrr2.NSRRID,Nsrr1.NSRRID(ii)));
     if ~isempty(indx)
        Nsrr1.Flag(ii)=1;
    else
        Nsrr1.Flag(ii)=0;
    end
end
% All sleep
Ahi3pArV1=Nsrr1.POOAHI3(Nsrr1.Flag==1); 
Ahi3pArV2=Nsrr2.POOAHI3;

%NREM
% Ahi3pArV1=str2double((Nsrr1.POAHI3N(Nsrr1.Flag==1)));
% Ahi3pArV2=str2double(Nsrr2.POAHI3N);

% %REM
% Ahi3pArV1=str2double(Nsrr1.POAHI3R(Nsrr1.Flag==1));
% Ahi3pArV2=str2double(Nsrr2.POAHI3R);


% Ahi3pArV1=str2double(T1_new.pordi3pa); % overall rdi at 3% desat or arousal
% Ahi0pArV1=str2double(T1_new.pordi0pa); % overall rdi at 0% desat or arousal

% Ahi3pArV2=str2double(T2_new.pordi3pa); 
% Ahi0pArV2=str2double(T2_new.pordi0pa);


corr(BMIV1,BMIV2,'row','Complete')
DeltaBMI =  BMIV2-BMIV1;

%% run SummaryTableV12 for summary



%% ALL SUBJECTS
Subj=(1:1:size(SummaryAnalysisTable1,1))';
Latlist={'FlatLoopGain','FlatArTh','FlatUA'};
FLatV1=table(Flateral1(:,1),Flateral1(:,4),Flateral1(:,6));
FLatV1.Properties.VariableNames = Latlist;
FLatV2=table(Flateral2(:,1),Flateral2(:,4),Flateral2(:,6));
FLatV2.Properties.VariableNames = Latlist;

Stagelist={'FN1LoopGain','FN2LoopGain','FN3LoopGain','FREMLoopGain'};
FStatesLoopGainV1=array2table(FstatesArray_1{1, 1}); 
FStatesLoopGainV1.Properties.VariableNames = Stagelist;
FStatesLoopGainV2=array2table(FstatesArray_2{1, 1});
FStatesLoopGainV2.Properties.VariableNames = Stagelist;

Stagelist={'FN1ArTh','FN2ArTh','FN3ArTh','FREMArTh'};
FStatesArThV1=array2table(FstatesArray_1{1, 4}); 
FStatesArThV1.Properties.VariableNames = Stagelist;
FStatesArThV2=array2table(FstatesArray_2{1, 4});
FStatesArThV2.Properties.VariableNames = Stagelist;

Stagelist={'FN1UA','FN2UA','FN3UA','FREMUA'};
FStatesUAV1=array2table(FstatesArray_1{1, 6}); 
FStatesUAV1.Properties.VariableNames = Stagelist;
FStatesUAV2=array2table(FstatesArray_2{1, 6});
FStatesUAV2.Properties.VariableNames = Stagelist;



%% Data Prep for Regression Models
traitsj=[1 2 3 4 5 6 7 8];
traitslist={'LG1','LGn','Delay','VRA','ArTh','Vpassive','Vactive','Vcomp'}';

thres=3*ones(length(traitsj),1); % minimum 3 windows needed to be considered

% exclude ppl with no trait at all and ahi <5 for table 1 

% from MESA within-night calculations of minimum windows
% thres=[12 24 20 24 16 30 35 37]; old
% thres=[33 29 29 20 20 35 38 38]'; % dummy values for delay,vra & vactive

% decision to be made to use NSRR ahi or our AHI?
% corr(Ahi3pArV1,AHI1,'rows','complete')
% corr(Ahi3pArV2,AHI2,'rows','complete')

OSAAHIthres=5;
OSA = ((Ahi3pArV1>OSAAHIthres)&(Ahi3pArV2>OSAAHIthres));
% OSA = ((AHI1>OSAAHIthres)&(AHI2>OSAAHIthres));

TSTmin=45;
TSTcrit=((TST1>TSTmin)&(TST2>TSTmin));

display(['Num of OSA subj common in V1 & V2: ',num2str(nansum(OSA))])

failedstudy = isnan(AHI1)|isnan(AHI2);
display(['Num of failed study in V1 or V2 based on AHI: ',num2str(nansum(failedstudy))])


temp = isnan(SummaryAnalysisTableN_1{:,:});
temp2 = SummaryAnalysisTableN_1{:,:};
temp2(temp)=0;
temp2(failedstudy,:)=NaN;
SummaryAnalysisTableN_1{:,:} = temp2;

temp = isnan(SummaryAnalysisTableN_2{:,:});
temp2 = SummaryAnalysisTableN_2{:,:};
temp2(temp)=0;
temp2(failedstudy,:)=NaN;
SummaryAnalysisTableN_2{:,:} = temp2;

nanmedian(SummaryAnalysisTableN_1{:,:});
nanmedian(SummaryAnalysisTableN_2{:,:});


%% Get Prediction Intervals & Beta values for Trait2~ Trait1
clear T PredIn N
for i=1:length(traitsj)
    j=traitsj(i);
    N1=SummaryAnalysisTableN_1{:,j}; % visit1
    N2=SummaryAnalysisTableN_2{:,j}; % visit2
    crit = N1>=thres(i) & N2>=thres(i) & OSA; % use only those subjects who has windows greater than threshold
       
    Trait2 = SummaryAnalysisTable2{crit==1,j}; % trait from visit2
    Trait1 = SummaryAnalysisTable1{crit==1,j}; % trait from visit1

    T = table(Trait1, Trait2);
    disp(['trait:', traitslist{j}]);
    mdl1 = fitglm(T,'Trait2 ~ 1')
%     [ypred,yci] = predict(mdl,T);
%      outdata = nanmean(Traitdel) + [0,-1.96, 1.96]*nanstd(Traitdel);
%          pValueAdjFREM(i,1) = mdl.Coefficients.Estimate(1);
     PredIn(i,:)=mdl1.Coefficients.Estimate(1)+[0,-1.96, 1.96]*mdl1.Coefficients.SE(1);
     N(i)=mdl1.NumObservations;
     mdl2=fitglm(T,'Trait2 ~ Trait1')
    BetaTraitsAll(i,1) = mdl2.Coefficients.Estimate(2:end);
    SETraitsAll(i,1)= mdl2.Coefficients.SE(2:end);
    PvalTraitsAll(i,1)=mdl2.Coefficients.pValue(2:end);
%     mdlRsqAll(i,1)= mdl2.Rsquared.Ordinary.^0.5

end

%% Regression Model for NREM exclusive
clear rhoNREMexclusive
traitsj=[6 8 1 2 5];

for i=1:length(traitsj)
    figure(1);
    subplot(2,length(traitsj),i)
    j=traitsj(i);
    N1=SummaryAnalysisTableNnr_1{:,j}; % visit1
    N2=SummaryAnalysisTableNnr_2{:,j}; % visit2
    crit = N1>=thres(j) & N2>=thres(j) & OSA; % use only those subjects who has windows greater than threshold
    
    Nactual(1,i)=nanmedian(N1(crit)); % median # of windows in the actual file
    Nactual(2,i)=nanmedian(N2(crit));
%     
    Trait2 = SummaryAnalysisTablenr2{crit==1,j}; % trait from visit2
    Trait1 = SummaryAnalysisTablenr1{crit==1,j}; % trait from visit1
    
     
    T = table(Trait1,Trait2);
    disp(['trait:', traitslist{j}]);
    mdl = fitglm(T,'Trait2 ~ Trait1')
    
   % pValueAdj(i,1) = mdl.Coefficients.pValue(end);
   % Ypred = predict(mdl,T);
    rhoNREMexclusive(i,1)=corr(T.Trait1,T.Trait2,'Type','Pearson','rows','complete')
   % rhoAdj(i,1)=corr(Ypred,T.Trait2,'Type','Pearson','rows','complete')
   BetaNREMexclusive(i,1) = mdl.Coefficients.Estimate(2:end);
   SENREMexclusive(i,1)= mdl.Coefficients.SE(2:end);
   PvalNREMexclusive(i,1)=mdl.Coefficients.pValue(2:end);
   traitslist1(i,1)=traitslist(j);
end
rhoNREMexclusive1=table(traitslist1,rhoNREMexclusive)

%% Regression Model for NREM vs REM, NREM depths, and Position
traitsj=[6 8 1 2 5];
for i=1:length(traitsj)
    figure(1);
    subplot(2,length(traitsj),i)
    j=traitsj(i);
    N1=SummaryAnalysisTableN_1{:,j}; % visit1
    N2=SummaryAnalysisTableN_2{:,j}; % visit2
    crit = N1>=thres(j) & N2>=thres(j) & OSA; % use only those subjects who has windows greater than threshold
    
    Nactual(1,i)=nanmedian(N1(crit)); % median # of windows in the actual file
    Nactual(2,i)=nanmedian(N2(crit));
    
    
    deltaFlateral = Flateral2(crit,j) - Flateral1(crit,j);
    figure(11); clf(11); histogram(deltaFlateral);
    %     deltaFsup = Fsupine_2(crit,j) - Fsupine_1(crit,j);
    
%     deltaAHI = AHI1(crit) - AHI2(crit);
%     figure(11); clf(11); histogram(deltaAHI);
%     
    deltaFN1 = FstatesArray_2{j}(crit,1) - FstatesArray_1{j}(crit,1); % change in n1
    figure(11); clf(11); histogram(deltaFN1);
    deltaFN3 = FstatesArray_2{j}(crit,3) - FstatesArray_1{j}(crit,3); % change in n3
    figure(11); clf(11); histogram(deltaFN3);
%     deltaFW = (1-sum(FstatesArray_2{j}(crit,:),2)) - (1-sum(FstatesArray_1{j}(crit,:),2)); % change in wake
%     figure(11); clf(11); histogram(deltaFW);
%        
    deltaFREM = FstatesArray_2{j}(crit,4) - FstatesArray_1{j}(crit,4); % change in n3
   % figure(11); clf(11); histogram(deltaFN3);
    
   deltaBMI = DeltaBMI(crit);
   
    Trait2 = SummaryAnalysisTable2{crit==1,j}; % trait from visit2
    Trait1 = SummaryAnalysisTable1{crit==1,j}; % trait from visit1
    
    MeanTrait1(i,1)=nanmean(Trait1);
    SETrait1(i,1)=[nanstd(Trait1)/sqrt(length(Trait1))];

    MeanTrait2(i,1)=nanmean(Trait2);
    SETrait2(i,1)=[nanstd(Trait2)/sqrt(length(Trait2))];

    MeanTraitDel(i,1)=nanmean(Trait2-Trait1);
    SETraitDel(i,1)=[nanstd(Trait2-Trait1)./sqrt(length(Trait2-Trait1))];

       
    T = table(Trait1,Trait2,deltaFlateral,deltaFN1,deltaFN3,deltaFREM,deltaBMI);
    T2 = T; 
    T2.deltaBMI = 0*T2.deltaBMI;
    
    disp(['trait:', traitslist{j}]);
    
    mdl = fitglm(T,'Trait2 ~ Trait1')
    Ypred = predict(mdl,T2);
    rho(i,1)=corr(T.Trait1,T.Trait2,'Type','Pearson','rows','complete')
    Trait1row = find(strcmp(mdl.Coefficients.Properties.RowNames,'Trait1'));
    BetaTrait12(i,1) = mdl.Coefficients.Estimate(Trait1row);
    SETrait12(i,1)= mdl.Coefficients.SE(Trait1row);
    PvalTrait12(i,1)=mdl.Coefficients.pValue(Trait1row);
    
    mdl = fitglm(T,'Trait2 ~ Trait1 + deltaFlateral')
    Posrow = find(strcmp( mdl.Coefficients.Properties.RowNames,'deltaFlateral'));
    pValueAdjPos(i,1) = mdl.Coefficients.pValue(Posrow);
    Ypred = predict(mdl,T);
    rhoAdjPos(i,1)=corr(Ypred,T.Trait2,'Type','Pearson','rows','complete')
    Trait1row = find(strcmp( mdl.Coefficients.Properties.RowNames,'Trait1'));
    BetaAdjPos(i,1) = mdl.Coefficients.Estimate(Trait1row);
    SEAdjPos(i,1)= mdl.Coefficients.SE(Trait1row);
    PvalAdjPos(i,1)=mdl.Coefficients.pValue(Trait1row);
    
    
    mdl = fitglm(T,'Trait2 ~ Trait1 + deltaFlateral+ deltaFREM')
    FREMrow = find(strcmp( mdl.Coefficients.Properties.RowNames,'deltaFREM'));
    pValueAdjFREM(i,1) = mdl.Coefficients.pValue(FREMrow);
    Ypred = predict(mdl,T);
    rhoAdjFREM(i,1)=corr(Ypred,T.Trait2,'Type','Pearson','rows','complete')
    Trait1row = find(strcmp( mdl.Coefficients.Properties.RowNames,'Trait1'))
    BetaAdjFREM(i,1) = mdl.Coefficients.Estimate(Trait1row);
    SEAdjFREM(i,1)= mdl.Coefficients.SE(Trait1row);
    PvalAdjFREM(i,1)=mdl.Coefficients.pValue(Trait1row);
    
    
    mdl = fitglm(T,'Trait2 ~ Trait1 + deltaFlateral+ deltaFREM +deltaFN1 + deltaFN3')
    mdlKeep{i}=mdl;
    N1row = find(strcmp( mdl.Coefficients.Properties.RowNames,'deltaFN1'));
    N3row = find(strcmp( mdl.Coefficients.Properties.RowNames,'deltaFN3'));
    pValueAdjN1(i,1) = mdl.Coefficients.pValue(N1row);
    pValueAdjN3(i,1) = mdl.Coefficients.pValue(N3row);
    Ypred = predict(mdl,T);
    rhoAdjStates(i,1)=corr(Ypred,T.Trait2,'Type','Pearson','rows','complete')
    %mdl.Rsquared.Ordinary.^0.5
    Trait1row = find(strcmp( mdl.Coefficients.Properties.RowNames,'Trait1'))
    BetaAdjN1N3(i,1) = mdl.Coefficients.Estimate(Trait1row);
    SEAdjN1N3(i,1)= mdl.Coefficients.SE(Trait1row);
    PvalAdjN1N3(i,1)=mdl.Coefficients.pValue(Trait1row);
    
    
    BetaAdjAll(:,i) = mdl.Coefficients.Estimate(2:end);
    SEAdjAll(:,i) = mdl.Coefficients.SE(2:end);
   
    %mdl.Rsquared.Ordinary.^0.5
    traitslist1(i,1)=traitslist(j);
end
%%
Trho3 = table(traitslist1,rhoNREMexclusive,rho,rhoAdjPos,pValueAdjPos,rhoAdjFREM,pValueAdjFREM,rhoAdjStates,pValueAdjN1,pValueAdjN3)
BetaTraitsAll=table(traitslist1,BetaTrait12,BetaNREMexclusive,BetaAdjPos,BetaAdjFREM,BetaAdjN1N3);
SETraitsAll=table(traitslist1,SETrait12,SENREMexclusive,SEAdjPos,SEAdjFREM,SEAdjN1N3);
PvalTraitsAll=table(traitslist1,PvalTrait12,PvalNREMexclusive,PvalAdjPos,PvalAdjFREM,PvalAdjN1N3);

BetaAdjAllT=array2table(BetaAdjAll.');
BetaAdjAllT.Properties.VariableNames=mdl.Coefficients.Properties.RowNames(2:end);
BetaAdjAllT.Properties.RowNames=traitslist1;
BetaAdjAllTpart = BetaAdjAllT(:,[2 5 3 4])
%Finding: Minimal effect of adj for position, except Vpassive

SEAdjAllT=array2table(SEAdjAll.');
SEAdjAllT.Properties.VariableNames=mdl.Coefficients.Properties.RowNames(2:end);
SEAdjAllT.Properties.RowNames=traitslist1;
SEAdjAllTpart = SEAdjAllT(:,[2 5 3 4])

MeanTable=table(traitslist1,MeanTrait1,SETrait1,MeanTrait2,SETrait2,MeanTraitDel,SETraitDel)



%% AHI %%
clear rhoNREMexclusiveAHI

Ahi3pArV1=Nsrr1.POOAHI3(Nsrr1.Flag==1);
Ahi3pArV2=Nsrr2.POOAHI3;

traitslist2={'AHI'}

% crit = OSA; % use only those subjects who has windows greater than threshold
crit=OSA & TSTcrit; % including MESA criteria

Trait2 = Ahi3pArV2(crit==1); % trait from visit2
Trait1 = Ahi3pArV1(crit==1); % trait from visit1


T = table(Trait1,Trait2);
disp(['trait:', traitslist{j}]);
mdl = fitglm(T,'Trait2 ~ Trait1')

% pValueAdj(i,1) = mdl.Coefficients.pValue(end);
% Ypred = predict(mdl,T);
rhoNREMexclusive=corr(T.Trait1,T.Trait2,'Type','Pearson','rows','complete')
% rhoAdj(i,1)=corr(Ypred,T.Trait2,'Type','Pearson','rows','complete')
BetaNREMexclusive = mdl.Coefficients.Estimate(2:end);
SENREMexclusive= mdl.Coefficients.SE(2:end);
PvalNREMexclusive=mdl.Coefficients.pValue(2:end);

rhoNREMexclusive1=table(traitslist2,rhoNREMexclusive)

%% AHI %% Regression Model for NREM vs REM, NREM depths, and Position

j=traitsj(3); % using flateral from LG1
crit = OSA ; % use only those subjects who has windows greater than threshold
%  crit=OSA & TSTcrit;

deltaFlateral = Flateral2(crit,j) - Flateral1(crit,j);
     
    deltaFN1 = FstatesArray_2{j}(crit,1) - FstatesArray_1{j}(crit,1); % change in n1
    deltaFN3 = FstatesArray_2{j}(crit,3) - FstatesArray_1{j}(crit,3); % change in n3
    deltaFREM = FstatesArray_2{j}(crit,4) - FstatesArray_1{j}(crit,4); % change in n3
    
   deltaBMI = DeltaBMI(crit);
   
   Trait2 = Ahi3pArV2(crit==1); % trait from visit2
   Trait1 = Ahi3pArV1(crit==1); % trait from visit1
   
   
    MeanTrait1=nanmean(Trait1);
    SETrait1=[nanstd(Trait1)/sqrt(length(Trait1))];

    MeanTrait2=nanmean(Trait2);
    SETrait2=[nanstd(Trait2)/sqrt(length(Trait2))];

    MeanTraitDel=nanmean(Trait2-Trait1);
    SETraitDel=[nanstd(Trait2-Trait1)./sqrt(length(Trait2-Trait1))];

       
    T = table(Trait1,Trait2,deltaFlateral,deltaFN1,deltaFN3,deltaFREM,deltaBMI);
    T2 = T; 
    T2.deltaBMI = 0*T2.deltaBMI;
    
   
    mdl = fitglm(T,'Trait2 ~ Trait1')
    Ypred = predict(mdl,T2);
    rho=corr(T.Trait1,T.Trait2,'Type','Pearson','rows','complete')
    Trait1row = find(strcmp(mdl.Coefficients.Properties.RowNames,'Trait1'));
    BetaTrait12= mdl.Coefficients.Estimate(Trait1row);
    SETrait12= mdl.Coefficients.SE(Trait1row);
    PvalTrait12=mdl.Coefficients.pValue(Trait1row);
    
    mdl = fitglm(T,'Trait2 ~ Trait1 + deltaFlateral')
    Posrow = find(strcmp( mdl.Coefficients.Properties.RowNames,'deltaFlateral'));
    pValueAdjPos = mdl.Coefficients.pValue(Posrow);
    Ypred = predict(mdl,T);
    rhoAdjPos=corr(Ypred,T.Trait2,'Type','Pearson','rows','complete')
    Trait1row = find(strcmp( mdl.Coefficients.Properties.RowNames,'Trait1'));
    BetaAdjPos = mdl.Coefficients.Estimate(Trait1row);
    SEAdjPos= mdl.Coefficients.SE(Trait1row);
    PvalAdjPos=mdl.Coefficients.pValue(Trait1row);
    
    
    mdl = fitglm(T,'Trait2 ~ Trait1 + deltaFlateral+ deltaFREM')
    FREMrow = find(strcmp( mdl.Coefficients.Properties.RowNames,'deltaFREM'));
    pValueAdjFREM= mdl.Coefficients.pValue(FREMrow);
    Ypred = predict(mdl,T);
    rhoAdjFREM=corr(Ypred,T.Trait2,'Type','Pearson','rows','complete')
    Trait1row = find(strcmp( mdl.Coefficients.Properties.RowNames,'Trait1'))
    BetaAdjFREM = mdl.Coefficients.Estimate(Trait1row);
    SEAdjFREM= mdl.Coefficients.SE(Trait1row);
    PvalAdjFREM=mdl.Coefficients.pValue(Trait1row);
    
    
    mdl = fitglm(T,'Trait2 ~ Trait1 + deltaFlateral+ deltaFREM +deltaFN1 + deltaFN3')
    mdlKeep=mdl;
    N1row = find(strcmp( mdl.Coefficients.Properties.RowNames,'deltaFN1'));
    N3row = find(strcmp( mdl.Coefficients.Properties.RowNames,'deltaFN3'));
    pValueAdjN1 = mdl.Coefficients.pValue(N1row);
    pValueAdjN3 = mdl.Coefficients.pValue(N3row);
    Ypred = predict(mdl,T);
    rhoAdjStates=corr(Ypred,T.Trait2,'Type','Pearson','rows','complete')
    %mdl.Rsquared.Ordinary.^0.5
    Trait1row = find(strcmp( mdl.Coefficients.Properties.RowNames,'Trait1'))
    BetaAdjN1N3 = mdl.Coefficients.Estimate(Trait1row);
    SEAdjN1N3= mdl.Coefficients.SE(Trait1row);
    PvalAdjN1N3=mdl.Coefficients.pValue(Trait1row);
    
    SEAdjAll = mdl.Coefficients.SE(2:end);
    BetaAdjAll = mdl.Coefficients.Estimate(2:end);
  
   
  traitslist2={'AHI'}


%%
Trho3 = table(traitslist2,rhoNREMexclusive,rho,rhoAdjPos,pValueAdjPos,rhoAdjFREM,pValueAdjFREM,rhoAdjStates,pValueAdjN1,pValueAdjN3)
BetaTraitsAll=table(traitslist2,BetaTrait12,BetaNREMexclusive,BetaAdjPos,BetaAdjFREM,BetaAdjN1N3);
SETraitsAll=table(traitslist2,SETrait12,SENREMexclusive,SEAdjPos,SEAdjFREM,SEAdjN1N3);
PvalTraitsAll=table(traitslist2,PvalTrait12,PvalNREMexclusive,PvalAdjPos,PvalAdjFREM,PvalAdjN1N3);

BetaAdjAllT=array2table(BetaAdjAll.');
BetaAdjAllT.Properties.VariableNames=mdl.Coefficients.Properties.RowNames(2:end);
BetaAdjAllT.Properties.RowNames=traitslist2;

SEAdjAllT=array2table(SEAdjAll.');
SEAdjAllT.Properties.VariableNames=mdl.Coefficients.Properties.RowNames(2:end);
SEAdjAllT.Properties.RowNames=traitslist2;

%Finding: Minimal effect of adj for position, except Vpassive

MeanTable=table(traitslist2,MeanTrait1,SETrait1,MeanTrait2,SETrait2,MeanTraitDel,SETraitDel)



%% Regression Model for NREM vs REM, NREM depths, and Position - adjusting for BMI

for i=1:length(traitsj)
    figure(1);
    subplot(2,length(traitsj),i)
    j=traitsj(i);
    N1=SummaryAnalysisTableN_1{:,j}; % visit1
    N2=SummaryAnalysisTableN_2{:,j}; % visit2
    crit = N1>=thres(i) & N2>=thres(i) & OSA; % use only those subjects who has windows greater than threshold
    
    Nactual(1,i)=nanmedian(N1(crit)); % median # of windows in the actual file
    Nactual(2,i)=nanmedian(N2(crit));
    
    
    deltaFlateral = Flateral2(crit,j) - Flateral1(crit,j);
    figure(11); clf(11); histogram(deltaFlateral);
    %     deltaFsup = Fsupine_2(crit,j) - Fsupine_1(crit,j);
    
%     deltaAHI = AHI1(crit) - AHI2(crit);
%     figure(11); clf(11); histogram(deltaAHI);
%     
    deltaFN1 = FstatesArray_2{j}(crit,1) - FstatesArray_1{j}(crit,1); % change in n1
    figure(11); clf(11); histogram(deltaFN1);
    deltaFN3 = FstatesArray_2{j}(crit,3) - FstatesArray_1{j}(crit,3); % change in n3
    figure(11); clf(11); histogram(deltaFN3);
%     deltaFW = (1-sum(FstatesArray_2{j}(crit,:),2)) - (1-sum(FstatesArray_1{j}(crit,:),2)); % change in wake
%     figure(11); clf(11); histogram(deltaFW);
%        
    deltaFREM = FstatesArray_2{j}(crit,4) - FstatesArray_1{j}(crit,4); % change in n3
   % figure(11); clf(11); histogram(deltaFN3);
    
   deltaBMI = DeltaBMI(crit);
   
    Trait2 = SummaryAnalysisTable2{crit==1,j}; % trait from visit2
    Trait1 = SummaryAnalysisTable1{crit==1,j}; % trait from visit1
   
    T = table(Trait1,Trait2,deltaFlateral,deltaFN1,deltaFN3,deltaFREM,deltaBMI);
    T2 = T; 
    T2.deltaBMI = 0*T2.deltaBMI;
    
    disp(['trait:', traitslist{j}]);
    
    mdl = fitglm(T,'Trait2 ~ Trait1 + deltaBMI')
    Ypred = predict(mdl,T2);
    BMIrow = find(strcmp( mdl.Coefficients.Properties.RowNames,'deltaBMI'));
    temp = T.Trait2 - mdl.Coefficients.Estimate(BMIrow).*T.deltaBMI;
    %rho(i,1)=corr(Ypred,T.Trait2,'Type','Pearson','rows','complete')
    rho(i,1)=corr(Ypred,temp,'Type','Pearson','rows','complete')
    
    mdl = fitglm(T,'Trait2 ~ Trait1 + deltaFREM + deltaBMI');
    FREMrow = find(strcmp( mdl.Coefficients.Properties.RowNames,'deltaFREM'));
    pValueAdjFREM(i,1) = mdl.Coefficients.pValue(FREMrow);
    Ypred = predict(mdl,T2);
    BMIrow = find(strcmp( mdl.Coefficients.Properties.RowNames,'deltaBMI'));
    temp = T.Trait2 - mdl.Coefficients.Estimate(BMIrow).*T.deltaBMI;
    rhoAdjFREM(i,1)=corr(Ypred,temp,'Type','Pearson','rows','complete');
    
    mdl = fitglm(T,'Trait2 ~ Trait1 + deltaFN1 + deltaFN3 + deltaFREM  + deltaBMI');
    N1row = find(strcmp( mdl.Coefficients.Properties.RowNames,'deltaFN1'));
    N3row = find(strcmp( mdl.Coefficients.Properties.RowNames,'deltaFN3'));

    pValueAdjN1(i,1) = mdl.Coefficients.pValue(N1row);
    pValueAdjN3(i,1) = mdl.Coefficients.pValue(N3row);
    Ypred = predict(mdl,T2);
    BMIrow = find(strcmp( mdl.Coefficients.Properties.RowNames,'deltaBMI'));
    temp = T.Trait2 - mdl.Coefficients.Estimate(BMIrow).*T.deltaBMI;
    rhoAdjStates(i,1)=corr(Ypred,temp,'Type','Pearson','rows','complete');
    %mdl.Rsquared.Ordinary.^0.5
    
    mdl = fitglm(T,'Trait2 ~ Trait1 + deltaFN1 + deltaFN3 + deltaFREM + deltaFlateral  + deltaBMI')
    Posrow = find(strcmp( mdl.Coefficients.Properties.RowNames,'deltaFlateral'));
        
    pValueAdjPos(i,1) = mdl.Coefficients.pValue(Posrow);
    Ypred = predict(mdl,T2);
    BMIrow = find(strcmp( mdl.Coefficients.Properties.RowNames,'deltaBMI'));
    temp = T.Trait2 - mdl.Coefficients.Estimate(BMIrow).*T.deltaBMI;
    rhoAdjPos(i,1)=corr(Ypred,temp,'Type','Pearson','rows','complete');
    %mdl.Rsquared.Ordinary.^0.5
    
end

TrhoBMIadj = table(traitslist,rho,rhoNREMexclusive,rhoAdjFREM,pValueAdjFREM,rhoAdjStates,pValueAdjN1,pValueAdjN3,rhoAdjPos,pValueAdjPos)
Trho3




%% UNUSED CODES  %%%%%%%%%%
%% VISIT 1 as a function of sleep states/position/age/bmi
V1Table=[table(Subj),SummaryAnalysisTable1,table(AgeV1,BMIV1,AHI1,Ahi3pArV1),FStatesLoopGainV1,FStatesArThV1,FStatesUAV1,FLatV1];
traitsj=[1 2 3 4 5 6 7 8];
traitslist={'LG1','LGn','Delay','VRA','ArTh','Vpassive','Vactive','Vcomp'}';

thres=3*ones(length(traitsj)); % minimum 3 windows needed to be considered

OSAAHIthres=5;
OSA =(AHI1>OSAAHIthres);

failedstudy = isnan(AHI1);
temp = isnan(SummaryAnalysisTableN_1{:,:});
temp2 = SummaryAnalysisTableN_1{:,:};
temp2(temp)=0;
temp2(failedstudy,:)=NaN;
SummaryAnalysisTableN_1{:,:} = temp2;

critwin=3;
j=1
traitslist(j)
crit=SummaryAnalysisTableN_1.LG1>critwin & OSA;
crit1(:,j)=crit;
V1Tabletemp=V1Table;

V1Tabletemp{(V1Tabletemp.LG1==0),2}=NaN;
V1Tabletemp{(V1Tabletemp.LGn==0),3}=NaN;
V1Tabletemp{(V1Tabletemp.delay==0),4}=NaN;
V1Tabletemp{(V1Tabletemp.VRA==0),5}=NaN;
V1Tabletemp{(V1Tabletemp.ArThres==0),6}=NaN;
V1Tabletemp{(V1Tabletemp.Vpassive==0),7}=NaN;
V1Tabletemp{(V1Tabletemp.Vactive==0),8}=NaN;
V1Tabletemp{(V1Tabletemp.Vcomp==0),9}=NaN;

V1Tabletemp{~crit,:}=NaN;
mdlAgeLg1V1=fitglme(V1Tabletemp,'LG1 ~ AgeV1+ (1|Subj)');
mdlBMILg1V1=fitglme(V1Tabletemp,'LG1 ~ BMIV1+ (1|Subj)');

mdl = fitglme(V1Tabletemp,'LG1 ~ AgeV1+ BMIV1+ FN1LoopGain + FN3LoopGain + FREMLoopGain + FlatLoopGain + (1|Subj)')
betaAll(:,j) = mdl.Coefficients.Estimate(2:end)
resV1.lg1 = residuals(mdl,'ResidualType','Pearson');
mufitV1.lg1= fitted(mdl);
figure (1); subplot(2,4,j); 
scatter(V1Tabletemp.LG1,resV1.lg1); xlabel('Lg1'); ylabel('residuals');title('LG1');

j=2
traitslist(j)
crit=SummaryAnalysisTableN_1.LGn>critwin & OSA;
crit1(:,j)=crit;
V1Tabletemp=V1Table;
V1Tabletemp{~crit,:}=NaN;
mdlAgeLgn1V1=fitglme(V1Tabletemp,'LGn ~ AgeV1+ (1|Subj)');
mdlBMILgnV1=fitglme(V1Tabletemp,'LGn ~ BMIV1+ (1|Subj)');

mdl = fitglme(V1Tabletemp,'LGn ~ AgeV1+ BMIV1+ FN1LoopGain + FN3LoopGain + FREMLoopGain + FlatLoopGain + (1|Subj)')
betaAll(:,j) = mdl.Coefficients.Estimate(2:end)
resV1.lgn = residuals(mdl,'ResidualType','Pearson');
mufitV1.lgn= fitted(mdl);
subplot(2,4,j);
scatter(V1Tabletemp.LGn,resV1.lgn);xlabel('LGn'); ylabel('residuals');title('LGn');

j=3
traitslist(j)
crit=SummaryAnalysisTableN_1.delay>critwin & OSA;
crit1(:,j)=crit;
V1Tabletemp=V1Table;
V1Tabletemp{~crit,:}=NaN;
mdl = fitglme(V1Tabletemp,'delay ~ AgeV1+ BMIV1+ FN1LoopGain + FN3LoopGain + FREMLoopGain + FlatLoopGain + (1|Subj)')
betaAll(:,j) = mdl.Coefficients.Estimate(2:end)
resV1.delay = residuals(mdl,'ResidualType','Pearson');
mufitV1.delay= fitted(mdl);
subplot(2,4,j);
scatter(V1Tabletemp.delay,resV1.delay);xlabel('Delay'); ylabel('residuals');title('delay');

j=4
traitslist(j)
crit=SummaryAnalysisTableN_1.VRA>critwin & OSA;
crit1(:,j)=crit;
V1Tabletemp=V1Table;
V1Tabletemp{~crit,:}=NaN;
mdl = fitglme(V1Tabletemp,'VRA ~ AgeV1+ BMIV1+ FN1ArTh + FN3ArTh + FREMArTh + FlatArTh + (1|Subj)')
betaAll(:,j) = mdl.Coefficients.Estimate(2:end)
resV1.vra = residuals(mdl,'ResidualType','Pearson');
mufitV1.vra= fitted(mdl);
subplot(2,4,j);
scatter(V1Tabletemp.VRA,resV1.vra);xlabel('VRA'); ylabel('residuals');title('VRA');


j=5
traitslist(j)
crit=SummaryAnalysisTableN_1.ArThres>critwin & OSA;
crit1(:,j)=crit;
V1Tabletemp=V1Table;
V1Tabletemp{~crit,:}=NaN;
mdlAgeArth1V1=fitglme(V1Tabletemp,'ArThres ~ AgeV1+ (1|Subj)');
mdlBMIArthV1=fitglme(V1Tabletemp,'ArThres ~ BMIV1+ (1|Subj)');

mdl = fitglme(V1Tabletemp,'ArThres ~ AgeV1+ BMIV1+ FN1ArTh + FN3ArTh + FREMArTh + FlatArTh + (1|Subj)')
betaAll(:,j) = mdl.Coefficients.Estimate(2:end)
resV1.arth = residuals(mdl,'ResidualType','Pearson');
mufitV1.arth= fitted(mdl);
subplot(2,4,j);
scatter(V1Tabletemp.ArThres,resV1.arth);xlabel('ArThres'); ylabel('residuals');title('ArThres');

j=6
traitslist(j)
crit=SummaryAnalysisTableN_1.Vpassive>critwin & OSA;
crit1(:,j)=crit;
V1Tabletemp=V1Table;
V1Tabletemp{~crit,:}=NaN;
mdl = fitglme(V1Tabletemp,'Vpassive ~ AgeV1+ BMIV1+ FN1UA + FN3UA + FREMUA + FlatUA + (1|Subj)')
betaAll(:,j) = mdl.Coefficients.Estimate(2:end)
resV1.vpassive = residuals(mdl,'ResidualType','Pearson');
mufitV1.vpassive= fitted(mdl);
subplot(2,4,j);
scatter(V1Tabletemp.Vpassive,resV1.vpassive);xlabel('Vpassive'); ylabel('residuals');title('Vpassive');


j=7
traitslist(j)
crit=SummaryAnalysisTableN_1.Vactive>critwin & OSA;
crit1(:,j)=crit;
V1Tabletemp=V1Table;
V1Tabletemp{~crit,:}=NaN;
mdl = fitglme(V1Tabletemp,'Vactive ~ AgeV1+ BMIV1+ FN1UA + FN3UA + FREMUA + FlatUA + (1|Subj)')
betaAll(:,j) = mdl.Coefficients.Estimate(2:end)
resV1.vactive = residuals(mdl,'ResidualType','Pearson');
mufitV1.vactive= fitted(mdl);
subplot(2,4,j);
scatter(V1Tabletemp.Vactive,resV1.vactive);xlabel('Vactive'); ylabel('residuals');title('Vactive');

j=8
traitslist(j)
crit=SummaryAnalysisTableN_1.Vcomp>critwin & OSA;
crit1(:,j)=crit;
V1Tabletemp=V1Table;
V1Tabletemp{~crit,:}=NaN;
mdl = fitglme(V1Tabletemp,'Vcomp ~ AgeV1+ BMIV1+ FN1UA + FN3UA + FREMUA + FlatUA + (1|Subj)')
betaAll(:,j) = mdl.Coefficients.Estimate(2:end)
resV1.vcomp = residuals(mdl,'ResidualType','Pearson');
mufitV1.vcomp= fitted(mdl);
subplot(2,4,j);
scatter(V1Tabletemp.Vcomp,resV1.vcomp);xlabel('Vcomp'); ylabel('residuals');title('Vcomp');



%% VISIT 2 as a function of sleep states/position/age/bmi
clear betaAll;
V2Table=[table(Subj),SummaryAnalysisTable2,table(AgeV2,BMIV2,AHI2,Ahi3pArV2),FStatesLoopGainV2,FStatesArThV2,FStatesUAV2,FLatV2];
traitsj=[1 2 3 4 5 6 7 8];
traitslist={'LG1','LGn','Delay','VRA','ArTh','Vpassive','Vactive','Vcomp'}';

thres=3*ones(length(traitsj)); % minimum 3 windows needed to be considered

OSAAHIthres=5;
OSA =(AHI2>OSAAHIthres);

failedstudy = isnan(AHI2);
clear temp;
temp = isnan(SummaryAnalysisTableN_2{:,:});
temp2 = SummaryAnalysisTableN_2{:,:};
temp2(temp)=0;
temp2(failedstudy,:)=NaN;
SummaryAnalysisTableN_2{:,:} = temp2;

critwin=3;
j=1
traitslist(j)
crit=SummaryAnalysisTableN_2.LG1>critwin & OSA;
crit2(:,j)=crit;
V2Tabletemp=V2Table;

V2Tabletemp{(V2Tabletemp.LG1==0),2}=NaN;
V2Tabletemp{(V2Tabletemp.LGn==0),3}=NaN;
V2Tabletemp{(V2Tabletemp.delay==0),4}=NaN;
V2Tabletemp{(V2Tabletemp.VRA==0),5}=NaN;
V2Tabletemp{(V2Tabletemp.ArThres==0),6}=NaN;
V2Tabletemp{(V2Tabletemp.Vpassive==0),7}=NaN;
V2Tabletemp{(V2Tabletemp.Vactive==0),8}=NaN;
V2Tabletemp{(V2Tabletemp.Vcomp==0),9}=NaN;


V2Tabletemp{~crit,:}=NaN;

mdlAgeLg1V2=fitglme(V2Tabletemp,'LG1 ~ AgeV2+ (1|Subj)');
mdlBMILg1V2=fitglme(V2Tabletemp,'LG1 ~ BMIV2+ (1|Subj)');

mdl = fitglme(V2Tabletemp,'LG1 ~ AgeV2+ BMIV2+ FN1LoopGain + FN3LoopGain + FREMLoopGain + FlatLoopGain + (1|Subj)')
betaAll(:,j) = mdl.Coefficients.Estimate(2:end)
resV2.lg1 = residuals(mdl,'ResidualType','Pearson');
mufitV2.lg1= fitted(mdl);
figure (1); subplot(2,4,j); 
scatter(mufitV2.lg1,resV2.lg1); xlabel('fitted'); ylabel('residuals');title('LG1');

j=2
traitslist(j)
crit=SummaryAnalysisTableN_2.LGn>critwin & OSA;
crit2(:,j)=crit;
V2Tabletemp=V2Table;
V2Tabletemp{~crit,:}=NaN;

mdlAgeLgnV2=fitglme(V2Tabletemp,'LGn ~ AgeV2+ (1|Subj)');
mdlBMILgnV2=fitglme(V2Tabletemp,'LGn ~ BMIV2+ (1|Subj)');


mdl = fitglme(V2Tabletemp,'LGn ~ AgeV2+ BMIV2+ FN1LoopGain + FN3LoopGain + FREMLoopGain + FlatLoopGain + (1|Subj)')
betaAll(:,j) = mdl.Coefficients.Estimate(2:end)
resV2.lgn = residuals(mdl,'ResidualType','Pearson');
mufitV2.lgn= fitted(mdl);
subplot(2,4,j);
scatter(mufitV2.lgn,resV2.lgn);xlabel('fitted'); ylabel('residuals');title('LGn');

j=3
traitslist(j)
crit=SummaryAnalysisTableN_2.delay>critwin & OSA;
crit2(:,j)=crit;
V2Tabletemp=V2Table;
V2Tabletemp{~crit,:}=NaN;
mdl = fitglme(V2Tabletemp,'delay ~ AgeV2+ BMIV2+ FN1LoopGain + FN3LoopGain + FREMLoopGain + FlatLoopGain + (1|Subj)')
betaAll(:,j) = mdl.Coefficients.Estimate(2:end)
resV2.delay = residuals(mdl,'ResidualType','Pearson');
mufitV2.delay= fitted(mdl);
subplot(2,4,j);
scatter(mufitV2.delay,resV2.delay);xlabel('fitted'); ylabel('residuals');title('delay');

j=4
traitslist(j)
crit=SummaryAnalysisTableN_2.VRA>critwin & OSA;
crit2(:,j)=crit;
V2Tabletemp=V2Table;
V2Tabletemp{~crit,:}=NaN;
mdl = fitglme(V2Tabletemp,'VRA ~ AgeV2+ BMIV2+ FN1ArTh + FN3ArTh + FREMArTh + FlatArTh + (1|Subj)')
betaAll(:,j) = mdl.Coefficients.Estimate(2:end)
resV2.vra = residuals(mdl,'ResidualType','Pearson');
mufitV2.vra= fitted(mdl);
subplot(2,4,j);
scatter(mufitV2.vra,resV2.vra);xlabel('fitted'); ylabel('residuals');title('VRA');


j=5
traitslist(j)
crit=SummaryAnalysisTableN_2.ArThres>critwin & OSA;
crit2(:,j)=crit;
V2Tabletemp=V2Table;
V2Tabletemp{~crit,:}=NaN;

mdlAgeArthV2=fitglme(V2Tabletemp,'ArThres ~ AgeV2+ (1|Subj)');
mdlBMIArthV2=fitglme(V2Tabletemp,'ArThres ~ BMIV2+ (1|Subj)');

mdl = fitglme(V2Tabletemp,'ArThres ~ AgeV2+ BMIV2+ FN1ArTh + FN3ArTh + FREMArTh + FlatArTh + (1|Subj)')
betaAll(:,j) = mdl.Coefficients.Estimate(2:end)
resV2.arth = residuals(mdl,'ResidualType','Pearson');
mufitV2.arth= fitted(mdl);
subplot(2,4,j);
scatter(mufitV2.arth,resV2.arth);xlabel('fitted'); ylabel('residuals');title('ArThres');

j=6
traitslist(j)
crit=SummaryAnalysisTableN_2.Vpassive>critwin & OSA;
crit2(:,j)=crit;
V2Tabletemp=V2Table;
V2Tabletemp{~crit,:}=NaN;
mdl = fitglme(V2Tabletemp,'Vpassive ~ AgeV2+ BMIV2+ FN1UA + FN3UA + FREMUA + FlatUA + (1|Subj)')
betaAll(:,j) = mdl.Coefficients.Estimate(2:end)
resV2.vpassive = residuals(mdl,'ResidualType','Pearson');
mufitV2.vpassive= fitted(mdl);
subplot(2,4,j);
scatter(mufitV2.vpassive,resV2.vpassive);xlabel('fitted'); ylabel('residuals');title('Vpassive');


j=7
traitslist(j)
crit=SummaryAnalysisTableN_2.Vactive>critwin & OSA;
crit2(:,j)=crit;
V2Tabletemp=V2Table;
V2Tabletemp{~crit,:}=NaN;
mdl = fitglme(V2Tabletemp,'Vactive ~ AgeV2+ BMIV2+ FN1UA + FN3UA + FREMUA + FlatUA + (1|Subj)')
betaAll(:,j) = mdl.Coefficients.Estimate(2:end)
resV2.vactive = residuals(mdl,'ResidualType','Pearson');
mufitV2.vactive= fitted(mdl);
subplot(2,4,j);
scatter(mufitV2.vactive,resV2.vactive);xlabel('fitted'); ylabel('residuals');title('Vactive');

j=8
traitslist(j)
crit=SummaryAnalysisTableN_2.Vcomp>critwin & OSA;
crit2(:,j)=crit;
V2Tabletemp=V2Table;
V2Tabletemp{~crit,:}=NaN;
mdl = fitglme(V2Tabletemp,'Vcomp ~ AgeV2+ BMIV2+ FN1UA + FN3UA + FREMUA + FlatUA + (1|Subj)')
betaAll(:,j) = mdl.Coefficients.Estimate(2:end)
resV2.vcomp = residuals(mdl,'ResidualType','Pearson');
mufitV2.vcomp= fitted(mdl);
subplot(2,4,j);
scatter(mufitV2.vcomp,resV2.vcomp);xlabel('fitted'); ylabel('residuals');title('Vcomp');

%% residual correlations
traitslist={'lg1','lgn','delay','vra','arth','vpassive','vactive','vcomp'}';
[Corr(1),pval(1)]=corr(resV1.lg1,resV2.lg1,'Rows','pairwise');
[Corr(2),pval(2)]=corr(resV1.lgn,resV2.lgn,'Rows','pairwise');
[Corr(3),pval(3)]=corr(resV1.delay,resV2.delay,'Rows','pairwise');
[Corr(4),pval(4)]=corr(resV1.vra,resV2.vra,'Rows','pairwise');
[Corr(5),pval(5)]=corr(resV1.arth,resV2.arth,'Rows','pairwise');
[Corr(6),pval(6)]=corr(resV1.vpassive,resV2.vpassive,'Rows','pairwise');
[Corr(7),pval(7)]=corr(resV1.vactive,resV2.vactive,'Rows','pairwise');
[Corr(8),pval(8)]=corr(resV1.vcomp,resV2.vcomp,'Rows','pairwise');
Corr
pval

%% percentage change
traitslist={'LG1','LGn','Delay','VRA','ArTh','Vpassive','Vactive','Vcomp'}';
traitsj=[1 2 3 4 5 6 7 8];
clear V1Tabletemp V2Tabletemp
V1Tabletemp=V1Table;
V1Tabletemp{(V1Tabletemp.LG1==0),2}=NaN;
V1Tabletemp{(V1Tabletemp.LGn==0),3}=NaN;
V1Tabletemp{(V1Tabletemp.delay==0),4}=NaN;
V1Tabletemp{(V1Tabletemp.VRA==0),5}=NaN;
V1Tabletemp{(V1Tabletemp.ArThres==0),6}=NaN;
V1Tabletemp{(V1Tabletemp.Vpassive==0),7}=NaN;
V1Tabletemp{(V1Tabletemp.Vactive==0),8}=NaN;
V1Tabletemp{(V1Tabletemp.Vcomp==0),9}=NaN;

V2Tabletemp=V2Table;
V2Tabletemp{(V2Tabletemp.LG1==0),2}=NaN;
V2Tabletemp{(V2Tabletemp.LGn==0),3}=NaN;
V2Tabletemp{(V2Tabletemp.delay==0),4}=NaN;
V2Tabletemp{(V2Tabletemp.VRA==0),5}=NaN;
V2Tabletemp{(V2Tabletemp.ArThres==0),6}=NaN;
V2Tabletemp{(V2Tabletemp.Vpassive==0),7}=NaN;
V2Tabletemp{(V2Tabletemp.Vactive==0),8}=NaN;
V2Tabletemp{(V2Tabletemp.Vcomp==0),9}=NaN;

for ii=1:size(traitsj,2)
    j=traitsj(ii);
    crit1temp=crit1(:,j);
    crit2temp=crit2(:,j);
    critboth=crit1temp & crit2temp;
    V2V1diff{j}=((V2Tabletemp.(j+1)(critboth)-V1Tabletemp.(j+1)(critboth))./V1Tabletemp.(j+1)(critboth))*100;
    MeanV2V1diff(j)=nanmean(V2V1diff{1,j});
    StdV2V1diff(j)=nanstd(V2V1diff{1,j});
    CorrV2V1(j)=corr(V2Tabletemp.(j+1)(critboth),V1Tabletemp.(j+1)(critboth),'rows','pairwise');
   
end
clear critboth crit1temp crit2temp
OSAAHIthres=5;
OSA1 = (Ahi3pArV1>OSAAHIthres)& AHI1>OSAAHIthres;
OSA2=(Ahi3pArV2>OSAAHIthres)& AHI2>OSAAHIthres;


V1Tabletemp{V1Tabletemp.Subj(~OSA1),2:end}=NaN;
V2Tabletemp{V2Tabletemp.Subj(~OSA2),2:end}=NaN;

crit1temp=[logical(ones(size(V1Tabletemp,1),1)),crit1,logical(ones(size(V1Tabletemp,1),4)),...
    repmat(crit1(:,1),1,4),repmat(crit1(:,5),1,4),repmat(crit1(:,6),1,4),...
    crit1(:,1),crit1(:,5),crit1(:,6)];

crit2temp=[logical(ones(size(V2Tabletemp,1),1)),crit2,logical(ones(size(V2Tabletemp,1),4)),...
    repmat(crit2(:,1),1,4),repmat(crit2(:,5),1,4),repmat(crit2(:,6),1,4),...
    crit2(:,1),crit2(:,5),crit2(:,6)];
V2V1Del(:,1)=table2array(V1Tabletemp(:,1));

for jj=2:size(V1Tabletemp,2)
    crit1temp1=crit1temp(:,jj-1);
    crit2temp2=crit2temp(:,jj-1);
    critboth=crit1temp1 & crit2temp2;
     V1temp=V1Tabletemp.(jj);
   V1temp(~critboth)=NaN;
   
    V2temp=V2Tabletemp.(jj);
   V2temp(~critboth)=NaN;
   V2V1Del(:,jj)=V2temp-V1temp;
   
   V2V1Per(:,jj)=V2V1Del(:,jj)./V1temp;
   
end
V2V1PerMean=nanmean(V2V1Per);


V2V1DelTable=array2table(V2V1Del);
V2V1DelTable.Properties.VariableNames=V1Tabletemp.Properties.VariableNames;


%% Delta traits as function of delta state/position etc.
traitsj=[1 2 3 4 5 6 7 8];
traitslist={'LG1','LGn','Delay','VRA','ArTh','Vpassive','Vactive','Vcomp'}';


j=1
traitslist(j)
V2V1DelTabletemp=V2V1DelTable;

mdlAgeLg1del=fitglme(V2V1DelTabletemp,'LG1 ~ AgeV1+ (1|Subj)');
mdlBMILg1del=fitglme(V2V1DelTabletemp,'LG1 ~ BMIV1+ (1|Subj)');

mdl = fitglme(V2V1DelTabletemp,'LG1 ~ AgeV1+ BMIV1+ Ahi3pArV1+ FN1LoopGain + FN3LoopGain + FREMLoopGain + FlatLoopGain + (1|Subj)')
betadel(:,j) = mdl.Coefficients.Estimate(2:end)
resV1.lg1 = residuals(mdl,'ResidualType','Pearson');
mufitV1.lg1= fitted(mdl);
figure (1); subplot(2,4,j); 
scatter(V2V1DelTabletemp.LG1,resV1.lg1); xlabel('Lg1'); ylabel('residuals');title('LG1');

j=2
mdlAgeLgn1V1=fitglme(V2V1DelTabletemp,'LGn ~ AgeV1+ (1|Subj)');
mdlBMILgnV1=fitglme(V2V1DelTabletemp,'LGn ~ BMIV1+ (1|Subj)');

mdl = fitglme(V2V1DelTabletemp,'LGn ~ AgeV1+ BMIV1+ Ahi3pArV1+FN1LoopGain + FN3LoopGain + FREMLoopGain + FlatLoopGain + (1|Subj)')
betadel(:,j) = mdl.Coefficients.Estimate(2:end)
resV1.lgn = residuals(mdl,'ResidualType','Pearson');
mufitV1.lgn= fitted(mdl);
subplot(2,4,j);
scatter(V2V1DelTabletemp.LGn,resV1.lgn);xlabel('LGn'); ylabel('residuals');title('LGn');

j=3
traitslist(j)
mdl = fitglme(V2V1DelTabletemp,'delay ~ AgeV1+ BMIV1+Ahi3pArV1+ FN1LoopGain + FN3LoopGain + FREMLoopGain + FlatLoopGain + (1|Subj)')
betadel(:,j) = mdl.Coefficients.Estimate(2:end)
resV1.delay = residuals(mdl,'ResidualType','Pearson');
mufitV1.delay= fitted(mdl);
subplot(2,4,j);
scatter(V2V1DelTabletemp.delay,resV1.delay);xlabel('Delay'); ylabel('residuals');title('delay');

j=4
traitslist(j)
mdl = fitglme(V2V1DelTabletemp,'VRA ~ AgeV1+ BMIV1+ Ahi3pArV1+FN1ArTh + FN3ArTh + FREMArTh + FlatArTh + (1|Subj)')
betadel(:,j) = mdl.Coefficients.Estimate(2:end)
resV1.vra = residuals(mdl,'ResidualType','Pearson');
mufitV1.vra= fitted(mdl);
subplot(2,4,j);
scatter(V2V1DelTabletemp.VRA,resV1.vra);xlabel('VRA'); ylabel('residuals');title('VRA');


j=5
mdlAgeArth1V1=fitglme(V2V1DelTabletemp,'ArThres ~ AgeV1+ (1|Subj)');
mdlBMIArthV1=fitglme(V2V1DelTabletemp,'ArThres ~ BMIV1+ (1|Subj)');

mdl = fitglme(V2V1DelTabletemp,'ArThres ~ AgeV1+ BMIV1+Ahi3pArV1+FN1ArTh + FN3ArTh + FREMArTh + FlatArTh + (1|Subj)')
betadel(:,j) = mdl.Coefficients.Estimate(2:end)
resV1.arth = residuals(mdl,'ResidualType','Pearson');
mufitV1.arth= fitted(mdl);
subplot(2,4,j);
scatter(V2V1DelTabletemp.ArThres,resV1.arth);xlabel('ArThres'); ylabel('residuals');title('ArThres');

j=6
traitslist(j)
mdl = fitglme(V2V1DelTabletemp,'Vpassive ~ AgeV1+ BMIV1+ Ahi3pArV1+FN1UA + FN3UA + FREMUA + FlatUA + (1|Subj)')
betadel(:,j) = mdl.Coefficients.Estimate(2:end)
resV1.vpassive = residuals(mdl,'ResidualType','Pearson');
mufitV1.vpassive= fitted(mdl);
subplot(2,4,j);
scatter(V2V1DelTabletemp.Vpassive,resV1.vpassive);xlabel('Vpassive'); ylabel('residuals');title('Vpassive');


j=7
traitslist(j)
mdl = fitglme(V2V1DelTabletemp,'Vactive ~ AgeV1+ BMIV1+Ahi3pArV1+FN1UA + FN3UA + FREMUA + FlatUA + (1|Subj)')
betadel(:,j) = mdl.Coefficients.Estimate(2:end)
resV1.vactive = residuals(mdl,'ResidualType','Pearson');
mufitV1.vactive= fitted(mdl);
subplot(2,4,j);
scatter(V2V1DelTabletemp.Vactive,resV1.vactive);xlabel('Vactive'); ylabel('residuals');title('Vactive');

j=8
traitslist(j)
mdl = fitglme(V2V1DelTabletemp,'Vcomp ~ AgeV1+ BMIV1+ Ahi3pArV1+FN1UA + FN3UA + FREMUA + FlatUA + (1|Subj)')
betadel(:,j) = mdl.Coefficients.Estimate(2:end)
resV1.vcomp = residuals(mdl,'ResidualType','Pearson');
mufitV1.vcomp= fitted(mdl);
subplot(2,4,j);
scatter(V2V1DelTabletemp.Vcomp,resV1.vcomp);xlabel('Vcomp'); ylabel('residuals');title('Vcomp');





%% V2 as a function of V1
traitslist={'LG1','LGn','Delay','VRA','ArTh','Vpassive','Vactive','Vcomp'}';
traitsj=[1 2 3 4 5 6 7 8];
for ii=1:size(traitsj,2)
    j=traitsj(ii);
    crit1temp=crit1(:,j);
    crit2temp=crit2(:,j);
    critboth=crit1temp & crit2temp;
    trait1=V1Tabletemp.(j+1)(critboth);
    trait2=V2Tabletemp.(j+1)(critboth);
    V2V1table=table(trait1,trait2);
    mdl=fitglme(V2V1table,'trait2 ~ trait1')
betaV2V1(:,j) = mdl.Coefficients.Estimate(2:end)

end
% mdl1 = fitglm(T,'Trait2 ~ 1')