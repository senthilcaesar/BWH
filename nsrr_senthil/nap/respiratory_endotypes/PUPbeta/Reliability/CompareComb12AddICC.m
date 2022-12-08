clear all; close all;
clc;

where1 = 'C:\Users\rma56\Dropbox (Partners HealthCare)\';
% where1 = 'G:\Partners Healthcare Dropbox\SATP Group\';
% where1 = 'C:\Users\szs88\Partners Healthcare Dropbox\SATP Group\';


addpath(genpath([where1 '\PUPbeta_git\PUPbeta']));
addpath(genpath([where1 'Reliability Trait Analysis\MESA\Reliability']));
addpath(genpath([where1 'MESA Phenotyping Population Factors\DataAnalysis']));

% load('combworkspace30min_2043sub.mat');
% load('SummaryAnalysis_Comb_OddEven.mat');
load('SummaryAnalysis_Comb_OddEven_Sep92021.mat');

load('BigTableSS.mat');
AhiTstWorkspace=load([where1 'Reliability Trait Analysis\MESA\Reliability\AhiTstWorkspace7212020.mat']);
% w = matfile('workspace_NREM_PX_TX.mat');
% w = matfile('SummaryAnalysis_AllSleep_AllPos_Boot.mat');
% w = matfile('SummaryAnalysis_AllSleep_AllPos.mat');
w = matfile('SummaryAnalysis_AllSleep_AllPos_Sep102021.mat');

%%  SET THRESHOLDS, AHI, TST
critwin=3;
traitsj=[1 2 5 6 7 8];
% traitsj=([2 5 6 8]);
upperNtrait=prctile(w.DataNArray,99);
thres=critwin*ones(length(traitsj)); % minimum 3 windows needed to be considered

AHI = BigTableSS.a0h3ai5; % AHI info from NSRR Table
TST = BigTableSS.slpprdp5; % TST from NSRR Table
OSAAHIthres=5;
OSA = AHI>OSAAHIthres;

%%  PSG Features from odd and even windows

PsgOdd=AhiTstWorkspace.CombTPSGParam(strcmp(AhiTstWorkspace.CombTPSGParam.Type,'Odd'),:);
PsgEven=AhiTstWorkspace.CombTPSGParam(strcmp(AhiTstWorkspace.CombTPSGParam.Type,'Even'),:);
PsgTotal=AhiTstWorkspace.CombTPSGParam(strcmp(AhiTstWorkspace.CombTPSGParam.Type,'Total'),:);


%% get MESA data parameters--
load('MESAdata.mat')
MESAdata.RemPer=MESAPSGextra.timeremp; %Calculated - pct time REM
MESAdata.RemPer(~OSA)=NaN;
MESAdata.N1Per=MESAPSGextra.timest1p;  %Calculated - pct time stage 1
MESAdata.N1Per(~OSA)=NaN;
MESAdata.N2Per=MESAPSGextra.timest2p;	%Calculated - pct time stage 2
MESAdata.N2Per(~OSA)=NaN;
MESAdata.N34Per=MESAPSGextra.times34p;	%Calculated - pct time stage 3-4
MESAdata.N34Per(~OSA)=NaN;

%AHIs from NSRR
MESAdata.a0h3aiNREM = MESAdata.oahi3pa_nrem5 + MESAdata.cai0p_nrem5;
MESAdata.a0h3aiREM = MESAdata.oahi3pa_rem5 + MESAdata.cai0p_rem5;


TotalEventN = ...
    MESAPSGextra.OARBP + MESAPSGextra.OAROP + MESAPSGextra.OANBP + MESAPSGextra.OANOP + ...
    MESAPSGextra.CARBP + MESAPSGextra.CAROP + MESAPSGextra.CANBP + MESAPSGextra.CANOP + ...
    MESAPSGextra.Urbpa3 + MESAPSGextra.Uropa3 + MESAPSGextra.Unrbpa3 + MESAPSGextra.unropa3 + ...
    MESAPSGextra.HREMBA3 + MESAPSGextra.HROA3 + MESAPSGextra.HNRBA3 + MESAPSGextra.HNROA3;

MESAdata.Flateral= 1 - (MESAdata.supinep5/100);



%% STARTING ENDOTYPES RELIABILITY

%% Get Fsupine and Number of Windows
if exist('Fsupine_compare')
    Fsupine_1=Fsupine;
    Fsupine_2=Fsupine_compare;
    SummaryAnalysisTable_1=SummaryAnalysisTable;
    SummaryAnalysisTable_2=SummaryAnalysisTable_compare;
    SummaryAnalysisTableN_1=SummaryAnalysisTableN;
    SummaryAnalysisTableN_2=SummaryAnalysisTableN_compare;
   
    clear FstatesArray_1 FstatesArray_2
    for i=1:length(Fstates)
        for j=1:8
            try
                FstatesArray_1{j}(i,:) = Fstates{i}(:,j)';
            catch me
                FstatesArray_1{j}(i,:) = [NaN NaN NaN NaN];
            end
        end
    end
     for i=1:length(Fstates_compare)
        for j=1:8
            try
                FstatesArray_2{j}(i,:) = Fstates_compare{i}(:,j)';
            catch me
                FstatesArray_2{j}(i,:) = [NaN NaN NaN NaN];
            end
        end
    end
end
    
Fsupine_1=AddRowsofNaNs(Fsupine_1,2060);
Fsupine_2=AddRowsofNaNs(Fsupine_2,2060);
for i=1:size(FstatesArray_1,2)
    FstatesArray_1{i}=AddRowsofNaNs(FstatesArray_1{i},2060);
    FstatesArray_2{i}=AddRowsofNaNs(FstatesArray_2{i},2060);
end
SummaryAnalysisTable_1=AddRowsofNaNs(SummaryAnalysisTable_1,2060);
SummaryAnalysisTable_2=AddRowsofNaNs(SummaryAnalysisTable_2,2060);
SummaryAnalysisTableN_1=AddRowsofZeros(SummaryAnalysisTableN_1,2060);
SummaryAnalysisTableN_2=AddRowsofZeros(SummaryAnalysisTableN_2,2060);


%% Set NaNs in N based on AHI

% medianwhole = [53    53    53    40    40    78    78    78] % median # of possible windows based on entire MESA data set
% lower25thwhole = [27   27   27   21   21   44   44   44] % 25th percentile for # of windows


% NaNs in N should be zero

failedstudy = isnan(AHI);
sum(OSA)
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

tempX = w.DataNArray;
temp = isnan(tempX);
tempX(temp)=0;
tempX(failedstudy,:)=NaN;

nanmedian(tempX)
nanmedian(SummaryAnalysisTableN_1{:,:})
nanmedian(SummaryAnalysisTableN_2{:,:})


%repeat for OSA
nanmedian(tempX(OSA,:))
nanmedian(SummaryAnalysisTableN_1{OSA,:})
nanmedian(SummaryAnalysisTableN_2{OSA,:})

%% TABLE 1
if 0
    %thresgood1 = [28    35    34    30    20] %using R=0.7
    
    minwincrit=3; %3; % 3 for normal use and 28 vpassive for r=0.7 use
    
    clear temp temp2
    temp=SummaryAnalysisTableN_1.Vpassive;
    temp(temp<minwincrit)=NaN;
    I1=~isnan(temp);
    temp2=SummaryAnalysisTableN_2.Vpassive;
    temp2(temp2<minwincrit)=NaN;
    I2=~isnan(temp2);
    sum(I1==1 & I2==1)
    I = I1&I2 & OSA;
    sum(I)
    
    if 0
    clear temp temp2 I
    temp = w.DataNArray;
    %     temp2 = w.DataArray;
    temp2=temp(:,6); % take vpassive
    temp2(temp2<minwincrit)=NaN;
    TrueN=~isnan(temp2);
    sum(TrueN==1 & OSA)
    I = TrueN & OSA;
    end
    
%     MESAdata.eventdur=DurEvt_all(:);
%     MESAdata.arousaldur=DurAr_all(:);
%     MESAdata.ODI3=ODI3_all(:);
%     MESAdata.ODI4=ODI4_all(:);
    
    Table1List = {'age5c','gender1','bmi5c',...
        'a0h3ai5','a0h3aiNREM','a0h3aiREM',...
        'ai_all5',...
        'slpprdp5','N1Per','N2Per','N34Per','RemPer',...
        'supinep5'};
    
    Table1Pt = MESAdata(I,Table1List);
    
    Table1Pt.male=Table1Pt.gender1==1;%female=0; male=1
    Table1Pt.female=Table1Pt.gender1==0;
    
    
    % for all vars with 0 and 1
    TableListN={'male','female'};
    
    Table1num=[];
    for ii=1:length(TableListN)
        temp=table2array(Table1Pt(:,TableListN{1,ii}));
        Table1num(1,ii)= sum(temp==1);
        Table1num(2,ii)=size(Table1Pt,1)-Table1num(1,ii);
    end
    OSAseverityN=[];
    temp=Table1Pt.a0h3ai5;
    OSAseverityN(1,:)= sum([temp<5 temp>=5&temp<15 temp>=15&temp<30 temp>=30]); % osa severity
    OSAseverityN(2,:)=size(Table1Pt,1)-OSAseverityN;
    Table1num=array2table(Table1num,'VariableNames',TableListN);
    Table1num=[Table1num array2table(OSAseverityN,'VariableNames',{'normal','mild','moderate','severe'})];
    Table1num.Properties.RowNames={'n','Total#-n'};
    
    
    
    % generate mean and 95CI for tranformed/nontransformed data
    TableListMean={'age5c','bmi5c',...
        'a0h3ai5','a0h3aiNREM','a0h3aiREM','ai_all5',...   
        'slpprdp5','N1Per','N2Per','N34Per','RemPer','supinep5'};
    
   %  'eventdur','arousaldur','ODI3','ODI4',...
    % 1 1 0.25 0.25 ...
    %  0 0 0 0 ...
    %  0 0 0 0  ...
    
    tempfactorx=[1 1 ...
        0.25 0.25 0.75 0.25 ...   
        1 0.5 1 0.5 1 1];
    
    logconv=[0 0 ...
        0 0 0 0 ...
        0 0 0 0 0 0];
    
    nonparametric = [0 0 ...
        0 0 0 0 ...
        0 0 0 0 0 1];
    
    
    Table1mean=[];
    
    if 0 % table 1 based on mean 95% CI
    for ii=1:length(TableListMean)
        tempfactor=tempfactorx(ii);
        logC=logconv(ii);
        if logC==1 && tempfactor==0
            F.Link = @(mu) log10(mu);
            F.Inverse = @(mu) 10.^(mu);
            tempfactor='log'; % for plotting
        else
            F.Link = @(mu) mu.^tempfactor;
            F.Inverse = @(mu) mu.^(1./tempfactor);
        end
        tempdata=table2array(Table1Pt(:,TableListMean{1,ii}));
        tempdatan=sum(~isnan(tempdata));
        
        figure (22); clf(22); subplot(2,1,1); histogram(tempdata);title(TableListMean{1,ii});
        tempdata=F.Link(tempdata);
        tempdata(tempdata==-Inf|tempdata==Inf)=NaN;
        [pValue(ii), ~] = swtest(tempdata)
        subplot(2,1,2); histogram(tempdata);
        title(['factor:',num2str(tempfactor),'; pValue:',num2str(pValue(ii))]);
        %     pause(3);
        if nonparametric(ii)
            outdata = [prctile(tempdata,[50 2.5 97.5])];
        else
            outdata = nanmean(tempdata) + [0,-1.96, 1.96]*nanstd(tempdata);
            outdata = F.Inverse(outdata);
        end
        Table1mean(ii,1:4)=[outdata tempdatan];
    end
    Table1mean=array2table(Table1mean,'VariableNames',{'Mean','LowerCI','UpperCI','N'});
    Table1mean.Properties.RowNames=TableListMean;
    
    else %table 1 based on median IQR
        for ii=1:length(TableListMean)
            
            tempdata=table2array(Table1Pt(:,TableListMean{1,ii}));
            tempdatan=sum(~isnan(tempdata));
            
            figure (22); clf(22); subplot(2,1,1); histogram(tempdata);title(TableListMean{1,ii});
            
            tempdata(tempdata==-Inf|tempdata==Inf)=NaN;
            [pValue(ii), ~] = swtest(tempdata)
            subplot(2,1,2); histogram(tempdata);
            %     pause(3);
            
            outdata = [prctile(tempdata,[50 25 75])];
%             iqr(tempdata)
             Table1mean(ii,1:4)=[outdata tempdatan];
        end
       
   
    Table1mean=array2table(Table1mean,'VariableNames',{'Mean','LowerIQR','UpperIQR','N'});
    Table1mean.Properties.RowNames=TableListMean;
    end
    
end
    




%% Model Endotype Error based on Windows
useconstant=0;
color1 = [0 0 0];
color2 = [0.7 0.4 0.2];
clear MSE x temp thresgood thresgood1 MSEwhole Ngood Ngood2 NgoodP NgoodP2
figure(1);clf (figure(1));
traitsLabelAll={'Loop Gain,LG1','Ventilatory Instability,LGn','Delay','VRA','Arousal Threshold',...
    'Collapsibility','Vactive','Compensation'};

traitsj=[6,8,1,2,5];

% traitsj=[1 2 3 4 5 6 7 8];
traitsLabelj=traitsLabelAll(traitsj);
thres=critwin*ones(1,8);
% allthres=[34 30 30 20 20 28 28 35];
thres=[34 30 30 20 20 28 28 35];

thresgood1 = [28    35    34    30    20] %using R=0.7

for i=1:length(traitsj)
    
    %     subplot(2,length(traitsj),i)
    j=traitsj(i);
    
    N1=SummaryAnalysisTableN_1{:,j}; % comb=1=odd windows
    N2=SummaryAnalysisTableN_2{:,j}; % comb=2=even windows
    
    % use only those subjects who has windows greater than threshold
    crit = N1>=thres(j) & N2>=thres(j)&OSA;
    crit1 = N1>=thres(j) & N2>=thres(j);
       
    
    Nsubj(:,i)=crit;
    Nactual(1,i)=nanmedian(N1(crit)); % median # of odd windows in the actual file
    Nactual(2,i)=nanmedian(N2(crit)); % median # of even windows in the actual file
    
    deltaFN1 = FstatesArray_2{j}(crit,1) - FstatesArray_1{j}(crit,1); % change in n1
    deltaFN3 = FstatesArray_2{j}(crit,3) - FstatesArray_1{j}(crit,3); % change in n3
    deltaFW = (1-sum(FstatesArray_2{j}(crit,:),2)) - (1-sum(FstatesArray_1{j}(crit,:),2)); % change in wake
    deltaFsup = Fsupine_2(crit,j) - Fsupine_1(crit,j); % change in supine position
    Trait2 = SummaryAnalysisTable_2{crit==1,j}; % trait meeting criteria--even windows
    Trait1 = SummaryAnalysisTable_1{crit==1,j}; % trait meeting criteria---odd window
    
    Trait2all = SummaryAnalysisTable_2{:,j}; % all trait irrspective of criteria-- even window
    Trait1all = SummaryAnalysisTable_1{:,j}; %  all trait irrspective of criteria--odd window
    
    
    temp = w.DataNArray;
    TraitTrueN = temp(crit==1,j);
    temp = w.DataArray;
    TraitTrueVal = temp(crit==1,j); %% true value from DataArray meeting criteria
    
    if j==4 % clip VAR
        Trait1(Trait1>90)=90;
        Trait2(Trait2>90)=90;
        Trait1all(Trait1all>90)=90;
        Trait2all(Trait2all>90)=90;
        TraitTrueVal(TraitTrueVal>90)=90;
    end
    
    if j==5 % transform ArTh
        Trait1 = fArThresT(Trait1);
        Trait2 = fArThresT(Trait2);
        Trait1all = fArThresT(Trait1all);
        Trait2all = fArThresT(Trait2all);
        TraitTrueVal = fArThresT(TraitTrueVal);
    end
    
    if j==6 % transform Vpassive
        Trait1 = fVpassiveT(Trait1);
        Trait2 = fVpassiveT(Trait2);
        Trait1all = fVpassiveT(Trait1all);
        Trait2all = fVpassiveT(Trait2all);
        TraitTrueVal = fVpassiveT(TraitTrueVal);
    end
    
    if j==7 % transform Vactive
        Trait1 = fVpassiveT(Trait1);
        Trait2 = fVpassiveT(Trait2);
        Trait1all = fVpassiveT(Trait1all);
        Trait2all = fVpassiveT(Trait2all);
        TraitTrueVal = fVpassiveT(TraitTrueVal);
    end
    
    if j==8 % limit
        VcompLimit=100;
        Trait1all(Trait1all>VcompLimit)=VcompLimit;
        Trait2all(Trait2all>VcompLimit)=VcompLimit;
        Trait2(Trait2>VcompLimit)=VcompLimit;
        Trait1(Trait1>VcompLimit)=VcompLimit;
    end
    
    % MEAN, SEM, DELTA TRAIT & CORRELATION
    TraitMeanOdd(j)=nanmean(Trait1);
    TraitMeanOddSE(j)=[nanstd(Trait1)/sqrt(length(Trait1))];
    TraitMeanOddSD(j)=nanstd(Trait1);
    TraitMeanEven(j)=nanmean(Trait2);
    TraitMeanEvenSE(j)=[nanstd(Trait2)/sqrt(length(Trait2))];
    TraitMeanEvenSD(j)=nanstd(Trait2);
    TraitDel(j)=nanmean(Trait2-Trait1);
    TraitDelSE(j)=[nanstd(Trait2-Trait1)/sqrt(length(Trait2-Trait1))];
    
    [RTrait2(j),pvalueRTrait2(j)] = corr(Trait1,Trait2,'Rows','complete')
    
    % PLOT TRAIT 1 vs TRAIT 2
    
    %     figure(88);
    %     if j==1
    %         subplot(2,3,1)
    %         scatter(Trait1,Trait2,10,'filled','markerfacealpha',0.5);
    %         xlabel('Loop Gain,LG1(odd)'); ylabel('Loop Gain,LG1(even)');
    %     elseif j==2
    %         subplot(2,3,2)
    %         scatter(Trait1,Trait2,10,'filled','markerfacealpha',0.5);
    %         xlabel('Ventilatory Instability(odd)'); ylabel('Ventilatory Instability(even)');
    %     elseif j==5
    %         subplot(2,3,3)
    %         scatter(Trait1,Trait2,10,'filled','markerfacealpha',0.5);
    %         xlabel('Arousal Threshold(odd)'); ylabel('Arousal Threshold(even)');
    %     elseif j==6
    %         subplot(2,3,4)
    %         scatter(Trait1,Trait2,10,'filled','markerfacealpha',0.5);
    %         xlabel('Collapsibility(odd)'); ylabel('Collapsibility(even)');
    %     elseif j==7
    %         subplot(2,3,5)
    %         scatter(Trait1,Trait2,10,'filled','markerfacealpha',0.5);
    %         xlabel('Vactive(odd)'); ylabel('Vactive(even)');
    %     elseif j==8
    %         subplot(2,3,6)
    %         scatter(Trait1,Trait2,10,'filled','markerfacealpha',0.5);
    %         xlabel('Compensation(odd)'); ylabel('Compensation(even)');
    %     end
    
    
    %% uncomment for model % commented out for getting plots and corr coeffs
    if 1
        % Model Trait 2 as a function of trait 1, N1,N3,Wake and supine
        T = table(Trait1,Trait2,deltaFW,deltaFN3,deltaFN1,deltaFsup);
        T.Trait2minusTrait1 = T.Trait2 - T.Trait1;
        
        
%         CenterPoint(i,1)=nanmedian(Trait1);
          CenterPoint(i,1:3)= prctile(Trait1,[33.3333 66.6666 99.9999])
    T.Trait1Class=(1*(T.Trait1<= CenterPoint(i,1))+2*(T.Trait1> CenterPoint(i,1)& T.Trait1<CenterPoint(i,2)))...
       + 3*(T.Trait1>= CenterPoint(i,2)) ; % low =1, med=2, high=3
    
    T.Trait2Class=(1*(T.Trait2<= CenterPoint(i,1))+2*(T.Trait2> CenterPoint(i,1)& T.Trait2<CenterPoint(i,2)))...
       + 3*(T.Trait2>= CenterPoint(i,2)) ; % low =1, med=2, high=3
        
    T.TraitSwitch=T.Trait1Class-T.Trait2Class;
    TraitClassT(1:3,i)=[sum(T.Trait1Class==1);sum( T.Trait1Class==2);sum( T.Trait1Class==3)]; 
    TraitClassT(4:6,i)=[sum(T.Trait2Class==1);sum( T.Trait2Class==2);sum( T.Trait2Class==3)]; 
   
    TraitClassT(7,i)=sum( T.Trait1Class==1 &   T.Trait2Class==1); % v1 low stayed in low
    TraitClassT(8,i)=sum( T.Trait1Class==1 &   T.Trait2Class==2); % v1 low to med
    TraitClassT(9,i)=sum( T.Trait1Class==1 &   T.Trait2Class==3); % v1 low to high
    
    TraitClassT(10,i)=sum( T.Trait1Class==2 &   T.Trait2Class==2); % v1 med stayed in med
    TraitClassT(11,i)=sum( T.Trait1Class==2 &   T.Trait2Class==1); % v1 med to low
    TraitClassT(12,i)=sum( T.Trait1Class==2 &   T.Trait2Class==3); % v1 med to high
    
    TraitClassT(13,i)=sum( T.Trait1Class==3 &   T.Trait2Class==3); % v1 high stayed in high
    TraitClassT(14,i)=sum( T.Trait1Class==3 &   T.Trait2Class==1); % v1 high to low
    TraitClassT(15,i)=sum( T.Trait1Class==3 &   T.Trait2Class==2); % v1 high to med

    TraitClassT(16:18,i)=[ TraitClassT(7,i)/ TraitClassT(1,i);TraitClassT(8,i)/ TraitClassT(1,i);...
        TraitClassT(9,i)/ TraitClassT(1,i);]
    TraitClassT(19:21,i)=[ TraitClassT(10,i)/ TraitClassT(2,i);TraitClassT(11,i)/ TraitClassT(2,i);...
        TraitClassT(12,i)/ TraitClassT(2,i);]
    TraitClassT(22:24,i)=[ TraitClassT(13,i)/ TraitClassT(3,i);TraitClassT(14,i)/ TraitClassT(3,i);...
        TraitClassT(15,i)/ TraitClassT(3,i);]
    
        
        
        
        mdl = fitglm(T,'Trait2 ~ Trait1 + deltaFN3 + deltaFN1 + deltaFW + deltaFsup'); % linear regression model of Trait 2
        Trait2predposstate = predict(mdl,T);
        Trait1adj = (predict(mdl,T) - mdl.Coefficients.Estimate(1))/mdl.Coefficients.Estimate(2); % adjusting for trait 1
        MSE(1,i) = nanmean((Trait1 - Trait2).^2); % mse of trait 1 vs 2
        MSE(2,i) = nanmean((Trait1adj - Trait2).^2); % mse after adjusting for trait 1
               
        % Model Trait 1-2 Error---
        iiii = ~isnan(Trait1) & ~isnan(Trait2); % include only not Nan values
        % find the error between odd and even traits
        MStot = nanmean(((Trait2(iiii) - nanmean(Trait2(iiii))).^2)); %simply variance
        
        if 0 %model abs difference
            Errsq = ((Trait1 - Trait2).^2)/MStot;
        else %model modelerror
            mdlA=fitglm(Trait1,Trait2);
            Trait2pred = predict(mdlA,Trait1);
            Errsq = ((Trait2 - Trait2pred).^2)/MStot;
            Trait2predall = predict(mdlA,Trait1all);
            Errsqall = ((Trait2all - Trait2predall).^2)/MStot;
        end
        Rtemp3(i) = corr(Trait1(iiii),Trait2(iiii),'Rows','complete')
        Rtemp3m(i) = mdlA.Rsquared.Ordinary.^0.5
        Rtemp3mB(i) = (1 - nanmean(Errsq)).^0.5
        
        
        
        %ICCs
        N=sum(iiii==1);
        xbar = mean(mean([Trait1(iiii) Trait2(iiii)])); %1/(2*N).*sum(Trait1(iiii)+Trait2(iiii))
        SSE1 = sum((Trait1(iiii)-xbar).^2);
        SSE2 = sum((Trait2(iiii)-xbar).^2);
        s2 = ((SSE1+SSE2)/2)/N; %meansquarederror
        ErrICC = (Trait1 - xbar).*(Trait2 - xbar);
        ErrICCnorm = ErrICC./s2;
        ICCx(1,i) = mean(ErrICCnorm(iiii))
                
        M = [Trait1(iiii) Trait2(iiii)];
        [n, k] = size(M);

        SStotal = var(M(:)) *(n*k - 1); %seems to be double the SStot, note is ~ var times k (2)
        MSR = var(mean(M, 2)) * k; %variance of mean-rater values, but times 2
        MSW = sum(var(M,[], 2)) / n; %mean variance due to simple measurement differences (half the abs difference, squared)
        MSC = var(mean(M, 1)) * n; %var due to bias e.g. difference in means per rater
        MSE = (SStotal - MSR*(n - 1) - MSC*(k-1))/ ((n - 1) * (k - 1)); %does not increase with bias
        ICCterms = [MSR MSW MSC MSE];
        ICCx(2,i) = (MSR - MSW) / (MSR + (k-1)*MSW); %'1-1'
        ICCx(3,i) = (MSR - MSE) / (MSR + (k-1)*MSE); %'C-1'
        ICCx(4,i) = (MSR - MSE) / (MSR + (k-1)*MSE + k*(MSC-MSE)/n); %'A-1'
    
        %attempt to calculate individual sample normalized error that could
        %be averaged in some way to provide ICC
%         MSWi = var(M,[], 2);
%         ICCtest = (MSR/2 - MSWi/2) / (MSR/2 + MSW/2) %'1-1'
%         ICCtest = (MSR/2) / (MSR/2 + MSW/2) - (MSWi/2) / (MSR/2 + MSW/2)%'1-1'
%         mean(ICCtest)
        
        
        %continue
        [ICCy(1,i), LB(1,i), UB(1,i), F, df1, df2, p(1,i)] = ICC([Trait1(iiii) Trait2(iiii)], '1-1')
        [ICCy(2,i), LB(2,i), UB(2,i), F, df1, df2, p(2,i)] = ICC([Trait1(iiii) Trait2(iiii)], 'C-1')
        [ICCy(3,i), LB(3,i), UB(3,i), F, df1, df2, p(3,i)] = ICC([Trait1(iiii) Trait2(iiii)], 'A-1') %A-1 A-k
%         [r, LB, UB, F, df1, df2, p] = ICC(M, type, alpha, r0)
        %results: ICC provided the same results as the Pearson correlation
        %coefficient, within 0.07% using the consistency method, and within 0.13% using the absolute agreement method.
        
        Rtemp3posstate(i) = corr(Trait1adj(iiii),Trait2(iiii),'Rows','complete')
        Errsqposstate = ((Trait2(iiii) - Trait2predposstate(iiii)).^2)/MStot;
        Rtemp3mBposstate(i) = (1 - nanmean(Errsqposstate)).^0.5

        Errsq((Errsq<0.000001*nanmedian(Errsq)))=0.000001*nanmedian(Errsq);
        
        % calculate mean Trait and number windows
        N_ = (N1(crit) + N2(crit))./2;
        meanTrait = (Trait1 + Trait2)/2;
        N_all = (N1 + N2)./2;
        meanTraitall = (Trait1all + Trait2all)/2;
        
        % for table
        N_AllMean(i)=nanmean(N_);
        N_AllSD(i)=nanstd(N_);
        N_AllSE(i)=[nanstd(N_)/sqrt(length(N_))];
        meanTraitallMean(i)=nanmean(meanTraitall);
        meanTraitallSE(i)=[nanstd(meanTraitall)/sqrt(length(meanTraitall))];
        
        clear F
        F.Link = @(mu) mu.^(-1); F.Derivative= @(mu) (-1)*mu.^((-1)-1); F.Inverse = @(mu) mu.^(1/(-1)); %link
        linkfun=F;
        
        %Nsqrt = N_.^0.5;
        % Predict the error as function of number of windows
        mdlErr = fitglm(N_,Errsq,'Link',linkfun,'Intercept',false);
        %mdlErr = fitglm(Nsqrt,Errsq);
        
        mdlErrKeep2{i}=mdlErr;
        MStot_(i) = MStot;
        SDtrait(i)=MStot^0.5;
        %     mdlErr = fitglm([N_ abs(meanTrait)],Errsq,'Link','log'); % original
        
        
        % Predict the error as function of number of windows and mean
        % trait
        
        Nline = (0:200)';
        %NlineT = 1./(Nline).^1;
        %Npower=1;
        %NT = N_.^Npower;
        if ~useconstant
            %mdlErr = fitglm(Nsqrt,Errsq,'Intercept',false);
            %mdlErr = fitglm(1./NT,Errsq,'Intercept',false);

           % mdlErr = fitglm([1./NT meanTrait-nanmean(meanTrait)],Errsq,'Intercept',false)
            mdlErr = fitglm([N_],Errsq,'Intercept',false,'Link',F)
            %mdlErr = fitglm([N_],Errsq,'Link',linkfun);
            %mdlErr = fitglm([Nsqrt meanTrait],Errsq);
            %             mdlErr = fitglm([N_ meanTrait],Errsq,'Link','log');

           
            mdlErrKeep{i}=mdlErr;
            %ErrsqModel = predict(mdlErr,[NlineT NlineT*0+nanmean(abs(meanTrait))]);
            ErrsqModel = predict(mdlErr,[Nline]);
            %ErrsqModel = predict(mdlErr,[NlineT]);
            %ErrsqModel = predict(mdlErr,[Nline]);
            
            % predict error for number of windows in "Npoints"
            Npoints = [10 60]; % for window length of 10 and 60
            ErrModel(1,i) = ErrsqModel(find(Nline==Npoints(1)))
            ErrModel(2,i) = ErrsqModel(find(Nline==Npoints(2)))
            
            %predict total err SSE (adjusted) for number of windows in whole data set
            ErrsqDataN = predict(mdlErr,[TraitTrueN]);
            % ErrsqDataN = predict(mdlErr,[1./TraitTrueN abs(TraitTrueVal)]);
            %ExclBelow3 = TraitTrueN<3 | isnan(TraitTrueN); % exclude subjects with window # below 3 or NaN
            %ErrsqDataN(ExclBelow3)=NaN;
            TraitTrueNMean(i)=nanmean(TraitTrueN);
            TraitTrueNSD(i)=nanstd(TraitTrueN);
            TraitTrueNSE(i)=[nanstd(TraitTrueN)/sqrt(length(TraitTrueN))];
            TraitTrueValMean(i)=nanmean(TraitTrueVal);
            TraitTrueValSE(i)=[nanstd(TraitTrueVal)/sqrt(length(TraitTrueVal))];
            TraitTrueValSD(i)=nanstd(TraitTrueVal);
           
            
            % MSE of error in the entire cohort
            MSEwhole(1,i) = nanmean(ErrsqDataN);
            
        else %model a flat line just to check R=R
            mdlErrConstant = fitglm([N_*0+1],Errsq,'Link',linkfun); %intercept is allowed in constant model!
            mdlErrKeep{i}=mdlErrConstant;
            ErrsqModel = predict(mdlErrConstant,[Nline*0+1]);
            
            
            % predict error for number of windows in "Npoints"
            Npoints = [10 60]; % for window length of 10 and 60
            ErrModel(1,i) = ErrsqModel(find(Nline==Npoints(1)));
            ErrModel(2,i) = ErrsqModel(find(Nline==Npoints(2)));
            
            
            %predict total err SSE (adjusted) for number of windows in whole data set
            ErrsqDataN = predict(mdlErrConstant,[TraitTrueN*0+1]);
            %ExclBelow3 = w.DataNArray(:,j)<3 | isnan(w.DataNArray(:,j)); % exclude subjects with window # below 3 or NaN
            %ErrsqDataN(ExclBelow3)=NaN;
            
            % MSE of error in the entire cohort
            MSEwhole(1,i) = nanmean(ErrsqDataN);
        end
        
        %repeated here:
        Rtemp3m(i) = mdlA.Rsquared.Ordinary.^0.5
        
        
        % repeat for OSA only; analysis is **unadjusted**
        if exist('temp'); clear temp; end
        if exist('temp2'); clear temp2; end
        temp = w.DataNArray(:,j);
        temp = temp(AHI>OSAAHIthres);
        temp2 = abs(w.DataArray(:,j));
        temp2 = temp2(AHI>OSAAHIthres);
        if ~useconstant
            %ErrsqDataOSAN = predict(mdlErr,[temp temp2]);
            ErrsqDataOSAN = predict(mdlErr,[temp]);
        else
            ErrsqDataOSAN = predict(mdlErrConstant,[temp*0]);
        end
        
        % MSE of error for OSA only
        ExclBelow3 = temp<3 | isnan(temp); % exclude sub with >3 windows
        ErrsqDataOSAN(ExclBelow3)=NaN;
        MSEwholeOSA(1,i) = nanmean(ErrsqDataOSAN);
        
        %finding the number of sub who had enough windows for the analysis
        rsqthresgood = (0.7)^2;
        %try
            thresgood1(i)=Nline(find(ErrsqModel<(1-rsqthresgood),1));
        %catch me
        %    thresgood1(i)=999;
       %     'no data above rsqthres'
        %end
        %         if thresgood1(i)>46
        %             thresgood1(i)=46;
        %         end
        
        %threshold based on percentile of windows in whole cohort.
        % currently not used
        %if 0
        %    thresgood2(i) = prctile(w.DataNArray(:,j),33);
        %    thresgood(i) = thresgood2(i);
        %end
        
        %thresgood(i) = thresgood1(i);
        
        rsqthresgood2 = (0.9)^2;
        %try
            thresgood2(i)=Nline(find(ErrsqModel<(1-rsqthresgood2),1));
        %catch me
        %    thresgood2(i)=999;
        %    'no data above rsqthres2'
        %end
        
        % Fraction of folks with rsquared >0.5 (R>0.7)
        Ngood(i) = sum(TraitTrueN>=thresgood1(i));
        NgoodP(i) = (Ngood(i)/length(TraitTrueN))*100
        
         Ngood2(i) = sum(TraitTrueN>=thresgood2(i));
        NgoodP2(i) = (Ngood2(i)/length(TraitTrueN))*100
        %     NgoodOSA(i) = sum(OSA&~(isnan(AHI)) & TraitTrueN(:,j)>=thresgood(i));
        %     FgoodOSA(i) = NgoodOSA(i)/sum(OSA&~(isnan(AHI)));
        %
        % fraction of controls
        %     NgoodNOSA(i) = sum(~OSA&~(isnan(AHI)) & TraitTrueN(:,j)>=thresgood(i));
        %     FgoodNOSA(i) = NgoodNOSA(i)/sum(~OSA&~(isnan(AHI)));
        
        
        
        % separate plot for error as function of N windows
        %         figure (222);
        %         subplot(3,length(traitsj),i)
        %         tempy=(Errsq).^0.5;
        %         tempy(tempy>2)=2;
        %         plot(N_,tempy,'.');
        %         hold('on');
        %         tempy=(ErrsqModel).^0.5;
        %         tempy(tempy>2)=2;
        %         plot(Nline,tempy,'k','linewidth',2);
        %         ylim([0 2]);
        %         yticks=get(gca,'ytick');
        %         set(gca,'ytick',yticks,'yticklabels',yticks.^2);
        %         xlim([0 upperNtrait(traitsj(i))]);
        %         box('off')
        %         xlabel('Windows'); ylabel('Error');
        %         title(traitsLabelj(i));
        %         set(gca,'FontSize',12)
        
        % plot the OSA vs control for trait 1 vs 2
       
        Ncommon=N1(crit)>thresgood1(i)& N2(crit)>thresgood1(i);
        
        Ncommon2=N1(crit)& N2(crit);
        
        % Bland-Altman
        BAMean=mean([Trait1(Ncommon),Trait2(Ncommon)],2,'omitnan');
        BADiff=Trait2(Ncommon)-Trait1(Ncommon);
        ZeroLine=zeros(length(BAMean),1);
        Meandiff=nanmean(BADiff);
        StdDev1=1.96*nanstd(BADiff);
        
        MeanLine=Meandiff*ones(length(BAMean),1);
        LOA(i,:)=[Meandiff+StdDev1 Meandiff-StdDev1];
        
%         CRcentile3win(i) = prctile(abs(T.Trait2minusTrait1(Ncommon)),95)      
        CRcentileNwin(i) = prctile(abs(T.Trait2minusTrait1(Ncommon)),95)
        % abs guarantees a real result
        bias=0;
        x=T.Trait2minusTrait1(Ncommon);
        CRgaussianNwin(i) = 1.96*(nansum(abs(x - 0).^2) ./ (length(x)-1))^0.5; %0 describes center point for CI %absolute method
        CRgaussianconsistencyNwin(i) = 1.96*nanstd(x);
 
        x=T.Trait2minusTrait1(Ncommon2);
        CRgaussian3win(i) = 1.96*(nansum(abs(x - 0).^2) ./ (length(x)-1))^0.5; %0 describes center point for CI %absolute method
        CRgaussianconsistency3win(i) = 1.96*nanstd(x);

             
        figure(1)
        set(gcf,'Position',get(0,'Screensize'));
        
        % 1st subplot
        % for subjects with >rsq 0.75
        subplot(3,length(traitsj)+1,i);
        scatter(Trait1(~Ncommon),Trait2(~Ncommon),8,color2,'filled','markerfacealpha',0.5);
        hold('on');
        % for subjects with <rsq 0.75
        scatter(Trait1(Ncommon),Trait2(Ncommon),8,color1,'filled','markerfacealpha',0.5);
        
        ylims = get(gca,'ylim');
        rangetemp = max(Trait1) - min(Trait1)
        ylims = [min(Trait1)-0.15*rangetemp max(Trait1)+0.15*rangetemp]
        set(gca,'ylim',ylims);
        
        set(gca,'xlim',ylims);
        yticks = get(gca,'ytick');
        set(gca,'xtick',yticks);
        set(gca,'ytick',yticks);
        set(gca,'FontSize',12);
        hold('on');
                
        xlims = get(gca,'xlim');
        plot(xlims,xlims,'color', [.5 .5 .5])
        if 0 %MrOS method
            plot(xlims,xlims+Meandiff,'b')
            plot(xlims,xlims+Meandiff+StdDev1,'r')
            plot(xlims,xlims+Meandiff-StdDev1,'r')
        else %MESA repeatability
            plot(xlims,xlims,'b')
            plot(xlims,xlims+CRcentileNwin(i),'r')
            plot(xlims,xlims-CRcentileNwin(i),'r')
        end
        hold off
        
%         if j==6
%             set(gca,'yticklabels',fVpassiveBT(yticks),'xticklabels',fVpassiveBT(yticks));
%         end
%         if j==5
%             set(gca,'yticklabels',fArThresT(yticks,1),'xticklabels',fArThresT(yticks,1));
% %             set(gca,'yticklabels',fArThresBT(yticks),'xticklabels',fArThresBT(yticks));
%         end

     
  % bland-altman
     subplot(3,length(traitsj)+1,i+length(traitsj)*2+2);%,i + (length(traitsj)+2)
     scatter(BAMean,BADiff,8,color1,'filled','markerfacealpha',0.5);
     hold on
     
     if 0 %MrOS method
         plot(BAMean,ZeroLine,'color', [.5 .5 .5],'LineWidth',1);
         plot(BAMean,MeanLine,'b','LineWidth',1);
         plot(BAMean,Meandiff+([1 -1]*StdDev1).*(ones(length(BAMean),2)),'r','LineWidth',1);
     else %MESA repeatability
         plot(BAMean,ZeroLine,'color', [.5 .5 .5],'LineWidth',1);
         plot(BAMean,ZeroLine,'b','LineWidth',1);
         plot(BAMean,ZeroLine+([1 -1]*CRcentileNwin(i)).*(ones(length(BAMean),2)),'r','LineWidth',1);
     end 
    
    
     set(gca,'ylim',diff(xlims)/2*[-1 1])
     set(gca,'xlim',xlims)
     set(gca,'xtick',yticks);
     set(gca,'FontSize',12);
     hold off
        
%          if j==6
%             set(gca,'yticklabels',fVpassiveBT(yticks),'xticklabels',fVpassiveBT(yticks));
%         end
%         if j==5
%             set(gca,'yticklabels',fArThresT(yticks,1),'xticklabels',fArThresT(yticks,1));
% %             set(gca,'yticklabels',fArThresBT(yticks),'xticklabels',fArThresBT(yticks));
%         end

        
        % Error vs Windows
        
        % 2nd subplot
        if 0
        subplot(3,length(traitsj)+1,i + 1*(length(traitsj)+1))
        %         tempy=(Errsq(OSA(crit))).^0.5;
        if 1
            tempy=(Errsq).^0.5;
            tempy(tempy>2)=2;
        else
            tempy=log10(Errsq);
            tempy(tempy<-5)=-5;
            tempy(tempy>1.5)=1.5;
        end
        %         plot(N_(OSA(crit)),tempy,'.'); % plotting OSA in blue color
        scatter(N_(~Ncommon),tempy(~Ncommon),8,color2,'filled','markerfacealpha',0.5);
        hold('on');
        scatter(N_(Ncommon),tempy(Ncommon),8,color1,'filled','markerfacealpha',0.5);
        %scatter(N_,tempy,8,color1,'filled','markerfacealpha',0.5);
        %         tempy=(Errsq(~OSA(crit))).^0.5;
        %         tempy(tempy>2)=2;
        %         plot(N_(~OSA(crit)),tempy,'.'); % controls in orange color
        hold('on');
        if 1
            % plot the 1-rsq2 line on OSA vs control plot
            if 1
                tempy=(ErrsqModel).^0.5;
                tempy(tempy>2)=NaN;
            else
                tempy=log10(ErrsqModel);
                tempy(tempy<-5)=NaN;
            end
            Iz = Nline>max(N_);
            plot(Nline(~Iz),tempy(~Iz),'color',[1 0.1 0.1],'linewidth',2);
            plot(Nline(Iz),tempy(Iz),'--','color',[1 0.7 0.7],'linewidth',2);
            %ylim([-3 1.5]);
            
            yticks=get(gca,'ytick');
            set(gca,'ytick',yticks,'yticklabels',yticks.^2);
            xlim([0 150]);
            %xlim([0 upperNtrait(traitsj(i))]);
            box('off')
            set(gca,'FontSize',12);
        end
        end
        
        % Rsqrd for trait 1 vs trait 2
        
        iii = ~isnan(Trait1) & ~isnan(Trait2);
        mdlx = fitglm(Trait1(iii),Trait2(iii))
        RTrait2_(j) = mdlx.Rsquared.Ordinary^0.5;
        
        SSE = nansum(((Trait1 - Trait2).^2)); % sum squared error
        SStot = nansum(((Trait1 - nanmean(Trait1)).^2));
        Rsq1(i) = 1 - SSE/SStot;
        
        
        % RsqEq for ErrModel based on number of windows in "Npoints"
        RsqEq(1,i) = 1 - ErrModel(1,i)/MStot;
        RsqEq(2,i) = 1 - ErrModel(2,i)/MStot;
        
        ErrSDRelative(1,i) = ErrModel(1,i).^0.5/MStot.^0.5;
        ErrSDRelative(2,i) = ErrModel(2,i).^0.5/MStot.^0.5;
        
        %% Histogram
        %
        figure(1);
        % 3rd subplot
        subplot(3,length(traitsj)+1,i + 1*(length(traitsj)+1))
%         subplot(3,length(traitsj)+1,i + (length(traitsj)+1)*2);
        
        
        dStep=5;
        Centers=0:dStep:upperNtrait(traitsj(i));
        Edges=(Centers(1)-dStep/2):dStep:(Centers(end)+dStep/2);
        Edges(end)=Inf;
        temp = w.DataNArray(:,j);
        temp = temp(crit ); %& AHI>15
        
        Ntemp=temp>thresgood1(i);
        [h7,edges] = histcounts(temp,Edges);
        bar(Centers,h7,'facecolor',color2,'EdgeAlpha',0,'BarWidth',1);
        hold on
        [h7,edges] = histcounts(temp(Ntemp),Edges);
        bar(Centers,h7,'facecolor',color1,'EdgeAlpha',0,'BarWidth',1);
        %hold('on');
        %         bar(Centers,h8,'EdgeAlpha',0,'FaceAlpha',0.7,'BarWidth',1);
        box('off')
        set(gca,'FontSize',12);
        
    end
end

traitslist={'Vpassive','Vcomp','LG1','LGn','ArTh'}';
try
    TraitClassTFinal=array2table(TraitClassT,'VariableNames',traitslist');
     TraitClassTFinal.Properties.RowNames={'V1Low','V1Med','V1High','V2Low','V2Med','V2High',...
        'V1Low-Low', 'V1Low-Med','V1Low-High',...
        'V1Med-Med', 'V1Med-Low','V1Med-High',...
        'V1High-High', 'V1High-Low','V1High-Med',...
        'V1Low-LowPer', 'V1Low-MedPer','V1Low-HighPer',...
        'V1Med-MedPer', 'V1Med-LowPer','V1Med-HighPer',...
        'V1High-HighPer', 'V1High-LowPer','V1High-MedPer'}
    
      
end

%
% ErrModel_ = ErrModel.^0.5 % ErrModel is for N=10 and 60 windows
%
% MSEwhole_ = MSEwhole.^0.5
% MSEwholeOSA_ = MSEwholeOSA.^0.5
% temp = w.DataArray;
% OSAstd = nanstd(temp(OSA,traitsj))
%
% temp2 = w.DataNArray;
% temp(temp2<3)=NaN;
% Allstd = nanstd(temp(:,traitsj))
%
% FstdOSA = MSEwholeOSA_ ./ OSAstd
% FvarOSA = MSEwholeOSA_.^2 ./ OSAstd.^2
% (1 - FvarOSA)
% (1 - FvarOSA).^0.5
%
% RTrait2(traitsj)
% Rsq1.^0.5
%
% FvarAll = MSEwhole_.^2 ./ Allstd.^2
% FvarAllx = MSEwhole_.^2 ./ OSAstd.^2
%
% (1 - FvarAll).^0.5
% (1 - FvarAllx).^0.5
%
% for i=1:length(mdlErrKeep)
%     i
%     mdlErrKeep{i}
%     mdlErrKeep2{i};
% end
%
% figure(1);
% subplot(4,4,1);
% axis([0 1 0 1])
% subplot(4,4,8);
% axis([-100 200 -100 200])

%
% 1hr
% MSE =
% 60 min, N=943 (min 30 w)
%     0.0114    0.0051    2.4533  189.6666  106.7375  371.5221
%     0.0113    0.0051    2.4407  170.8697   88.4524  361.8645
% 20 min, N=773 (min 30 w)
%    0.0116    0.0040    1.7693  145.8985   76.5824  339.4778
%     0.0114    0.0039    1.6736  130.2368   69.8442  334.1815
% 20 min, N=969 (min 25 w)
%    0.0120    0.0043    2.8908  184.4956   77.0917  336.6160
%    0.0119    0.0043    2.8116  170.0229   71.1170  330.2710
%
%
%     0.0128    0.0055    3.0505  227.7889  114.3025  359.7883
%     0.0127    0.0054    3.0145  202.9438   93.3785  353.5746


%60
%     0.0065    0.0032    1.6047   83.5396   29.6252   18.4348
%     0.0023    0.0010    0.2634   34.9690   15.7121   18.4348
%     0.0067    0.0032    1.7272   82.5880   33.0852   25.6954
%     0.0023    0.0010    0.2758   32.8536   15.9870   25.6954

%
% ErrModel =
%
%     0.0051    0.0026    1.1947   54.8980   23.7526   17.5982
%     0.0014    0.0006    0.1023   20.0726    8.7060   17.5982
%     0.0054    0.0026    1.2878   74.7426   25.4331   23.0944
%     0.0015    0.0006    0.1036   14.4702    9.3859   23.0944
%

Rcorrected = (1-MSEwhole).^0.5

%% make graphs look pretty
figure(1);
set(gcf,'color',[1 1 1])



%% PSG Parameters Modeling %%
clear temp thresgood thresgood1 

failedstudy = isnan(AHI);
PsgOdd{failedstudy,4:end}=NaN;
PsgEven{failedstudy,4:end}=NaN;
PsgTotal{failedstudy,4:end}=NaN;

useconstant=0;
color1 = [0 0 0];
color2 = [0.7 0.4 0.2];

clear Trait1 Trait2 crit Trait2all Trait1all TraitTrueN TraitTrueVal Trait1adj Trait2pred
figure(2);clf (figure(2));

traitsLabelAll={'AHI'};
traitsj=[1];
traitsLabelj=traitsLabelAll(traitsj);

for i=1:length(traitsj)
    
    %     j=traitsj(i);
   
    
    TST1=PsgOdd.TST;
    TST2=PsgEven.TST;
    
     crit = OSA ;
     crit = OSA & TST1>45 & TST2>45;
     sum(crit)
    
    deltaFN1PSG = PsgEven.FN1(crit) - PsgOdd.FN1(crit); % change in n1
    deltaFN3PSG = PsgEven.FN3(crit) - PsgOdd.FN3(crit); % change in n3
    deltaFsupPSG = (1-PsgEven.FLateral(crit)) -(1- PsgOdd.FLateral(crit)); % change in supine position
    
    Idx=find(strcmp(PsgOdd.Properties.VariableNames,traitsLabelj{i})==1);
    
    Trait2 = PsgEven{crit==1,Idx}; % trait meeting criteria--even windows
    Trait1 = PsgOdd{crit==1,Idx}; % trait meeting criteria---odd window
    
    Trait2all =PsgEven{:,Idx}; % all trait irrspective of criteria-- even window
    Trait1all = PsgOdd{:,Idx};
    
    TraitTrueN = PsgTotal.TST(crit==1);
    TraitTrueVal = table2array(PsgTotal(crit==1,Idx));
    
    TraitMeanOddPSG(i)=nanmean(Trait1);
    TraitMeanOddSEPSG(i)=[nanstd(Trait1)/sqrt(length(Trait1))];
    TraitMeanOddSDPSG(i)=nanstd(Trait1);
    TraitMeanEvenPSG(i)=nanmean(Trait2);
    TraitMeanEvenSEPSG(i)=[nanstd(Trait2)/sqrt(length(Trait2))];
    TraitMeanEvenSDPSG(i)=nanstd(Trait2);
    
    TraitDelPSG(i)=nanmean(Trait2-Trait1);
    TraitDelSEPSG(i)=[nanstd(Trait2-Trait1)/sqrt(length(Trait2-Trait1))];
    
    [RTrait2PSG(i),pvalueRTrait2(i)] = corr(Trait1,Trait2,'Rows','complete')
    
    if 1
        % Model Trait 2 as a function of trait 1, N1,N3 and supine
        TPSG = table(Trait1,Trait2,deltaFN3PSG,deltaFN1PSG,deltaFsupPSG);
        TPSG.Trait2minusTrait1=TPSG.Trait2-TPSG.Trait1;
        
        
        j=6;
       
        %         CenterPoint(j,1)=nanmedian(Trait1);
        CenterPoint(j,1:3)= prctile(Trait1,[33.3333 66.6666 99.9999]);
        TPSG.Trait1Class=(1*(TPSG.Trait1<= CenterPoint(j,1))+2*(TPSG.Trait1> CenterPoint(j,1)& TPSG.Trait1<CenterPoint(j,2)))...
            + 3*(TPSG.Trait1>= CenterPoint(j,2)) ; % low =1, med=2, high=3
        
        TPSG.Trait2Class=(1*(TPSG.Trait2<= CenterPoint(j,1))+2*(TPSG.Trait2> CenterPoint(j,1)& TPSG.Trait2<CenterPoint(j,2)))...
            + 3*(TPSG.Trait2>= CenterPoint(j,2)) ; % low =1, med=2, high=3
        
        TPSG.TraitSwitch=TPSG.Trait1Class-TPSG.Trait2Class;
        TraitClassT(1:3,j)=[sum(TPSG.Trait1Class==1);sum( TPSG.Trait1Class==2);sum( TPSG.Trait1Class==3)];
        TraitClassT(4:6,j)=[sum(TPSG.Trait2Class==1);sum( TPSG.Trait2Class==2);sum( TPSG.Trait2Class==3)];
        
        TraitClassT(7,j)=sum( TPSG.Trait1Class==1 &   TPSG.Trait2Class==1); % v1 low stayed in low
        TraitClassT(8,j)=sum( TPSG.Trait1Class==1 &   TPSG.Trait2Class==2); % v1 low to med
        TraitClassT(9,j)=sum( TPSG.Trait1Class==1 &   TPSG.Trait2Class==3); % v1 low to high
        
        TraitClassT(10,j)=sum( TPSG.Trait1Class==2 &   TPSG.Trait2Class==2); % v1 med stayed in med
        TraitClassT(11,j)=sum( TPSG.Trait1Class==2 &   TPSG.Trait2Class==1); % v1 med to low
        TraitClassT(12,j)=sum( TPSG.Trait1Class==2 &   TPSG.Trait2Class==3); % v1 med to high
        
        TraitClassT(13,j)=sum( TPSG.Trait1Class==3 &   TPSG.Trait2Class==3); % v1 high stayed in high
        TraitClassT(14,j)=sum( TPSG.Trait1Class==3 &   TPSG.Trait2Class==1); % v1 high to low
        TraitClassT(15,j)=sum( TPSG.Trait1Class==3 &   TPSG.Trait2Class==2); % v1 high to med
        
        TraitClassT(16:18,j)=[ TraitClassT(7,j)/ TraitClassT(1,j);TraitClassT(8,j)/ TraitClassT(1,j);...
            TraitClassT(9,j)/ TraitClassT(1,j)];
        TraitClassT(19:21,j)=[ TraitClassT(10,j)/ TraitClassT(2,j);TraitClassT(11,j)/ TraitClassT(2,j);...
            TraitClassT(12,j)/ TraitClassT(2,j)];
        TraitClassT(22:24,j)=[ TraitClassT(13,j)/ TraitClassT(3,j);TraitClassT(14,j)/ TraitClassT(3,j);...
            TraitClassT(15,j)/ TraitClassT(3,j)];
        
        
        
        mdlPSG = fitglm(TPSG,'Trait2 ~ Trait1 + deltaFN3PSG + deltaFN1PSG + deltaFsupPSG'); % linear regression model of Trait 2
        Trait1adj = (predict(mdlPSG,TPSG) - mdlPSG.Coefficients.Estimate(1))/mdlPSG.Coefficients.Estimate(2); % adjusting for trait 1
        MSEPSG(1,i) = nanmean((Trait1 - Trait2).^2); % mse of trait 1 vs 2
        MSEPSG(2,i) = nanmean((Trait1adj - Trait2).^2); % mse after adjusting for trait 1
        
        % Model Trait 1-2 Error---
        iiii = ~isnan(Trait1) & ~isnan(Trait2); % include only not Nan values
        % find the error between odd and even traits
        MStotPSG = nanmean(((Trait2(iiii) - nanmean(Trait2(iiii))).^2)); %simply variance
        
        if 0 %model abs difference
            ErrsqPSG = ((Trait1 - Trait2).^2)/MStotPSG;
            
        else %model modelerror
            mdlAPSG=fitglm(Trait1,Trait2);
            Trait2pred = predict(mdlAPSG,Trait1);
            ErrsqPSG = ((Trait2 - Trait2pred).^2)/MStotPSG;
            Trait2predall = predict(mdlAPSG,Trait1all);
            ErrsqallPSG = ((Trait2all - Trait2predall).^2)/MStotPSG;
        end
        Rtemp3PSG(i) = corr(Trait1(iiii),Trait2(iiii),'Rows','complete')
        Rtemp3mPSG(i) = mdlAPSG.Rsquared.Ordinary.^0.5
        Rtemp3mBPSG(i) = (1 - nanmean(ErrsqPSG)).^0.5
        
        ErrsqPSG((ErrsqPSG<0.000001*nanmedian(ErrsqPSG)))=0.000001*nanmedian(ErrsqPSG);
        
        % calculate mean Trait and TST
        TST_ = (TST1(crit) + TST2(crit))./2;
        meanTrait = (Trait1 + Trait2)/2;
        TST_all = (TST1 + TST2)./2;
        meanTraitall = (Trait1all + Trait2all)/2;
        
        % for table
        TST_AllMean(i)=nanmean(TST_);
        TST_AllSD(i)=nanstd(TST_);
        TST_AllSE(i)=[nanstd(TST_)/sqrt(length(TST_))];
        meanTraitallMeanPSG(i)=nanmean(meanTraitall);
        meanTraitallSDPSG(i)=nanstd(meanTraitall);
        meanTraitallSEPSG(i)=[nanstd(meanTraitall)/sqrt(length(meanTraitall))];
        
        [ICCy(1,6), LB(1,6), UB(1,6), F, df1, df2, p(1,6)] = ICC([Trait1(iiii) Trait2(iiii)], '1-1')
        [ICCy(2,6), LB(2,6), UB(2,6), F, df1, df2, p(2,6)] = ICC([Trait1(iiii) Trait2(iiii)], 'C-1')
        [ICCy(3,6), LB(3,6), UB(3,6), F, df1, df2, p(3,6)] = ICC([Trait1(iiii) Trait2(iiii)], 'A-1') %A-1 A-k
        
        % Predict the error as function of number of TST
        clear F
        F.Link = @(mu) mu.^(-1); F.Derivative= @(mu) (-1)*mu.^((-1)-1); F.Inverse = @(mu) mu.^(1/(-1)); %link
        linkfun=F;
        
        %Nsqrt = N_.^0.5;
        % Predict the error as function of number of windows
        %mdlErr = fitglm(N_,Errsq,'Link',linkfun,'Intercept',false);
        
        mdlErrPSG = fitglm(TST_,ErrsqPSG,'Link',linkfun,'Intercept',false);
        mdlErrKeep2PSG{i}=mdlErrPSG;
        MStot_PSG(i) = MStotPSG;
        SDtraitPSG(i)=MStotPSG^0.5;
        %     mdlErrPSG = fitglm([N_ abs(meanTrait)],ErrsqPSG,'Link','log'); % original
        
        % Predict the error as function of TST and mean
        % trait
        Nline = (0:600)'; % 10 hrs
        if ~useconstant
            %mdlErrPSG = fitglm([TST_ meanTrait],ErrsqPSG,'Link','log');
            mdlErrPSG = fitglm([TST_],ErrsqPSG,'Link',F,'Intercept',false);
            %             mdlErrPSG = fitglm([N_ meanTrait],ErrsqPSG,'Link','log');
            
%             if mdlErrPSG.Coefficients.Estimate(2)>0
%                 mdlErrPSG = fitglm([TST_*0 meanTrait],ErrsqPSG,'Link','log');
%                 %                 mdlErrPSG = fitglm([N_*0 meanTrait],ErrsqPSG,'Link','log');
%             end
            mdlErrKeepPSG{i}=mdlErrPSG;
            %ErrsqModelPSG = predict(mdlErrPSG,[Nline Nline*0+nanmean(abs(meanTrait))]);
            ErrsqModelPSG = predict(mdlErrPSG,[Nline]);
            
            % predict error for number of windows in "Npoints"
            Npoints = [60 360]; % for window length of 1hr and 6 hrs
            ErrModelPSG(1,i) = ErrsqModelPSG(find(Nline==Npoints(1)));
            ErrModelPSG(2,i) = ErrsqModelPSG(find(Nline==Npoints(2)));
            
            %predict total err SSE (adjusted) for number of windows in whole data set
            %ErrsqDataNPSG = predict(mdlErrPSG,[TraitTrueN abs(TraitTrueVal)]);
            ErrsqDataNPSG = predict(mdlErrPSG,[TraitTrueN ]);
            %ExclBelow3 = TraitTrueN<3 | isnan(TraitTrueN); % exclude subjects with window # below 3 or NaN
            %ErrsqDataN(ExclBelow3)=NaN;
            TraitTrueNMeanPSG(i)=nanmean(TraitTrueN);
            TraitTrueNSDPSG(i)=nanstd(TraitTrueN);
            TraitTrueNSEPSG(i)=[nanstd(TraitTrueN)/sqrt(length(TraitTrueN))];
            TraitTrueValMeanPSG(i)=nanmean(TraitTrueVal);
            TraitTrueValSDPSG(i)=nanstd(TraitTrueVal);
            TraitTrueValSEPSG(i)=[nanstd(TraitTrueVal)/sqrt(length(TraitTrueVal))];
            
            % MSE of error in the entire cohort
            MSEwholePSG(1,i) = nanmean(ErrsqDataNPSG);
            
        else %model a flat line just to check R=R [not updated]
            mdlErrConstantPSG = fitglm([TST_*0+1],ErrsqPSG,'Link','log');
            mdlErrKeepPSG{i}=mdlErrConstantPSG;
            ErrsqModelPSG = predict(mdlErrConstantPSG,[Nline*0+1]);
            
            
            % predict error for number of windows in "Npoints"
            Npoints = [60 360]; % for window length of 1hr and 6 hrs
            ErrModelPSG(1,i) = ErrsqModelPSG(find(Nline==Npoints(1)));
            ErrModelPSG(2,i) = ErrsqModelPSG(find(Nline==Npoints(2)));
            
            
            %predict total err SSE (adjusted) for number of windows in whole data set
            ErrsqDataNPSG = predict(mdlErrConstantPSG,[TraitTrueN*0+1]);
            %ExclBelow3 = w.DataNArray(:,j)<3 | isnan(w.DataNArray(:,j)); % exclude subjects with window # below 3 or NaN
            %ErrsqDataNPSG(ExclBelow3)=NaN;
            
            % MSE of error in the entire cohort
            MSEwholePSG(1,i) = nanmean(ErrsqDataNPSG);
        end
        
        Rtemp3mPSG(i) = mdlAPSG.Rsquared.Ordinary.^0.5
        
        
        %finding the number of sub who had enough TST for the analysis
        rsqthresgood = (0.7)^2;
        try
            thresgood1(i)=Nline(find(ErrsqModelPSG<(1-rsqthresgood),1));
        catch me
            thresgood1(i)=999;
            'no data above rsqthres'
        end

        rsqthresgood2 = (0.9)^2;
        try
            thresgood2(i)=Nline(find(ErrsqModelPSG<(1-rsqthresgood2),1));
        catch me
            thresgood2(i)=999;
            'no data above rsqthres'
        end
        
        % Fraction of folks with rsquared >0.5 (R>0.7)
        Ngood(i) = sum(TraitTrueN>=thresgood1(i));
        NgoodP(i) = (Ngood(i)/length(TraitTrueN))*100
        
        Ngood2(i) = sum(TraitTrueN>=thresgood2(i));
        NgoodP2(i) = (Ngood2(i)/length(TraitTrueN))*100
        
        % Rsqrd for trait 1 vs trait 2
        
        iii = ~isnan(Trait1) & ~isnan(Trait2);
        mdlxPSG = fitglm(Trait1(iii),Trait2(iii))
        RTrait2_(j) = mdlxPSG.Rsquared.Ordinary^0.5;
        
        SSEPSG = nansum(((Trait1 - Trait2).^2)); % sum squared error
        SStotPSG = nansum(((Trait1 - nanmean(Trait1)).^2));
        Rsq1PSG(i) = 1 - SSEPSG/SStotPSG;
        
        % RsqEq for ErrModelPSG based on number of windows in "Npoints"
        RsqEqPSG(1,i) = 1 - ErrModelPSG(1,i)/MStotPSG;
        RsqEqPSG(2,i) = 1 - ErrModelPSG(2,i)/MStotPSG;
        
        ErrSDRelativePSG(1,i) = ErrModelPSG(1,i).^0.5/MStotPSG.^0.5;
        ErrSDRelativePSG(2,i) = ErrModelPSG(2,i).^0.5/MStotPSG.^0.5;
        
        
        RcorrectedPSG = (1-MSEwholePSG).^0.5
        
        
%         figure(2); clf (2)
        % 1st subplot
        % for subjects with >rsq 0.75
       
        Ncommon=TST1(crit)>thresgood1(i)& TST2(crit)>thresgood1(i);
        Ncommon2=TST1(crit)& TST2(crit);
        
       
        
        BAMean=mean([Trait1(Ncommon),Trait2(Ncommon)],2,'omitnan');
        BADiff=Trait2(Ncommon)-Trait1(Ncommon);
        ZeroLine=zeros(length(BAMean),1);
        Meandiff=nanmean(BADiff);
        StdDev1=1.96*nanstd(BADiff);
        
        MeanLine=Meandiff*ones(length(BAMean),1);
        
        i=6
        LOA(i,:)=[Meandiff+StdDev1 Meandiff-StdDev1];
                
        CRcentile3win(i) = prctile(abs( TPSG.Trait2minusTrait1(Ncommon2)),95)      
        CRcentileNwin(i) = prctile(abs( TPSG.Trait2minusTrait1(Ncommon)),95)
        
      % abs guarantees a real result
        bias=0;
        x= TPSG.Trait2minusTrait1(Ncommon);
        CRgaussianNwin(i) = 1.96*(nansum(abs(x - 0).^2) ./ (length(x)-1))^0.5; %0 describes center point for CI %absolute method
        CRgaussianconsistencyNwin(i) = 1.96*nanstd(x);
 
        x= TPSG.Trait2minusTrait1(Ncommon2);
        CRgaussian3win(i) = 1.96*(nansum(abs(x - 0).^2) ./ (length(x)-1))^0.5; %0 describes center point for CI %absolute method
        CRgaussianconsistency3win(i) = 1.96*nanstd(x);
 
      
        figure(1);
        subplot(3,6,6);
        scatter(Trait1(~Ncommon),Trait2(~Ncommon),8,color2,'filled','markerfacealpha',0.5);
        hold('on');
        scatter(Trait1(Ncommon),Trait2(Ncommon),8,color1,'filled','markerfacealpha',0.5);
        ylims = get(gca,'ylim');
        
        rangetemp = max(Trait1) - min(Trait1)
        ylims = [min(Trait1)-0.15*rangetemp max(Trait1)+0.15*rangetemp]
        
        set(gca,'ylim',ylims);
        
        set(gca,'xlim',ylims);
        yticks = get(gca,'ytick');
        set(gca,'xtick',yticks);
        set(gca,'ytick',yticks);
        set(gca,'FontSize',12);
        hold('on');
        
        xlims = get(gca,'xlim');
        plot(xlims,xlims,'color', [.5 .5 .5])
        if 0 %MrOS method
            plot(xlims,xlims+Meandiff,'b')
            plot(xlims,xlims+Meandiff+StdDev1,'r')
            plot(xlims,xlims+Meandiff-StdDev1,'r')
        else %MESA repeatability
            plot(xlims,xlims,'b')
            plot(xlims,xlims+CRcentileNwin(i),'r')
            plot(xlims,xlims-CRcentileNwin(i),'r')
        end
        hold off
        
   subplot(3,6,18);%,i + (length(traitsj)+2)
   scatter(BAMean,BADiff,8,color1,'filled','markerfacealpha',0.5);
   hold on
   
   if 0 %MrOS method
       plot(BAMean,ZeroLine,'color', [.5 .5 .5],'LineWidth',1);
       plot(BAMean,MeanLine,'b','LineWidth',1);
       plot(BAMean,Meandiff+([1 -1]*StdDev1).*(ones(length(BAMean),2)),'r','LineWidth',1);
   else %MESA repeatability
       plot(BAMean,ZeroLine,'color', [.5 .5 .5],'LineWidth',1);
       plot(BAMean,ZeroLine,'b','LineWidth',1);
       plot(BAMean,ZeroLine+([1 -1]*CRcentileNwin(i)).*(ones(length(BAMean),2)),'r','LineWidth',1);
   end
   %      xlabel('Mean'); ylabel('Difference');
   %      title(traitslist(traitsj(i)))
   
   
   set(gca,'ylim',diff(xlims)/2*[-1 1])
   set(gca,'xlim',xlims)
   set(gca,'xtick',yticks);
   set(gca,'FontSize',12);
   hold off
   
   figure(1);
   set(gcf,'color',[1 1 1])


        % Error vs Windows
        
        % 2nd subplot
        if 0
        subplot(3,6,12)
        
        if 1
            tempy=(ErrsqPSG).^0.5;
            tempy(tempy>2)=2;
        else
            tempy=log10(ErrsqPSG);
            tempy(tempy<-5)=-5;
            tempy(tempy>1.5)=1.5;
        end
        scatter(TST_(~Ncommon),tempy(~Ncommon),8,color2,'filled','markerfacealpha',0.5);
        hold('on');
        scatter(TST_(Ncommon),tempy(Ncommon),8,color1,'filled','markerfacealpha',0.5);
        hold('on');
        if 1
            % plot the 1-rsq2 line on OSA vs control plot
            if 1
                tempy=(ErrsqModelPSG).^0.5;
                tempy(tempy>2)=NaN;
            else
                tempy=log10(ErrsqModelPSG);
                tempy(tempy<-5)=NaN;
            end
            Iz = Nline>max(TST_);
            plot(Nline(~Iz),tempy(~Iz),'color',[1 0.1 0.1],'linewidth',2);
            plot(Nline(Iz),tempy(Iz),'--','color',[1 0.7 0.7],'linewidth',2);
            %ylim([-3 1.5]);
            
            yticks=get(gca,'ytick');
            set(gca,'ytick',yticks,'yticklabels',yticks.^2);
            xlim([0 650]);
            box('off')
            set(gca,'FontSize',12);
        end
        end
        %% Histogram
        %
%         figure(1);
        % 3rd subplot
        subplot(3,6,12);
        
        
        dStep=15;
        Centers=0:dStep:prctile(PsgTotal{:,11},99);
        Edges=(Centers(1)-dStep/2):dStep:(Centers(end)+dStep/2);
        Edges(end)=Inf;
        temp = table2array(PsgTotal(:,11));
        temp = temp(crit ); %& AHI>5
        
        Ntemp=temp>thresgood1;
        [h7,edges] = histcounts(temp,Edges);
        bar(Centers,h7,'facecolor',color2,'EdgeAlpha',0,'BarWidth',1);
        hold on
        [h7,edges] = histcounts(temp(Ntemp),Edges);
        bar(Centers,h7,'facecolor',color1,'EdgeAlpha',0,'BarWidth',1);
        %hold('on');
        %         bar(Centers,h8,'EdgeAlpha',0,'FaceAlpha',0.7,'BarWidth',1);
        box('off')
        set(gca,'FontSize',12);
        figure(1);
        set(gcf,'color',[1 1 1])
    end
end


try
    TraitClassTFinal=[TraitClassTFinal,table(TraitClassT(:,i),'VariableNames',{'AHI'})];
      
end

%copygraphics(ax,'ContentType','vector')
