clear all;
clc;

addpath(genpath('C:\Users\rma56\Dropbox (Partners HealthCare)\PUPbeta_git\PUPbeta'));
addpath(genpath('C:\Users\rma56\Dropbox (Partners HealthCare)\Reliability Trait Analysis\MESA\Reliability'));
addpath(genpath('C:\Users\rma56\Dropbox (Partners HealthCare)\MESA Phenotyping Population Factors\DataAnalysis'));

load('combworkspace30min_2043sub.mat');
load('BigTableSS.mat');
load('AHI_ODI_workspace');
w = matfile('workspace_NREM_PX_TX.mat');

%% endotypes traits
traitsj=[1 2 3 4 5 6 7 8];
% traitsj=([2 5 6 8]);
upperNtrait=prctile(w.DataNArray,99);
AHI = BigTableSS.a0h3ai5; % AHI info from NSRR Table
TST = BigTableSS.slpprdp5; % TST from NSRR Table
OSAAHIthres=5;
OSA = AHI>OSAAHIthres;

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

%AHIs
MESAdata.a0h3aiNREM = MESAdata.oahi3pa_nrem5 + MESAdata.cai0p_nrem5;
MESAdata.a0h3aiREM = MESAdata.oahi3pa_rem5 + MESAdata.cai0p_rem5;


TotalEventN = ...
MESAPSGextra.OARBP + MESAPSGextra.OAROP + MESAPSGextra.OANBP + MESAPSGextra.OANOP + ...
    MESAPSGextra.CARBP + MESAPSGextra.CAROP + MESAPSGextra.CANBP + MESAPSGextra.CANOP + ...
    MESAPSGextra.Urbpa3 + MESAPSGextra.Uropa3 + MESAPSGextra.Unrbpa3 + MESAPSGextra.unropa3 + ...
    MESAPSGextra.HREMBA3 + MESAPSGextra.HROA3 + MESAPSGextra.HNRBA3 + MESAPSGextra.HNROA3;

MESAdata.Flateral= 1 - (MESAdata.supinep5/100);

%% Cell2mats for AHI, ODI etc.
AHItot_all=cell2mat(AHIall_total)';
AHItot_odd=cell2mat(AHIodd_total)';
AHItot_even=cell2mat(AHIeven_total)';
ODI3_all=cell2mat(ODI3_all);
ODI3_odd=cell2mat(ODI3_odd);
ODI3_even=cell2mat(ODI3_even);
ODI4_all=cell2mat(ODI4_all);
ODI4_odd=cell2mat(ODI4_odd);
ODI4_even=cell2mat(ODI4_even);
MeanOdd=[nanmean(AHItot_odd(OSA)) nanmean(ODI3_odd(OSA)) nanmean(ODI4_odd(OSA))];
MeanOddSE=[nanstd(AHItot_odd(OSA))/sqrt(length((AHItot_odd(OSA)))) nanstd(ODI3_odd(OSA))/sqrt(length((ODI3_odd(OSA))))...
    nanstd(ODI4_odd(OSA))/sqrt(length((ODI4_odd(OSA))))];
MeanEven=[nanmean(AHItot_even(OSA)) nanmean(ODI3_even(OSA)) nanmean(ODI4_even(OSA))];
MeanEvenSE=[nanstd(AHItot_even(OSA))/sqrt(length((AHItot_even(OSA)))) nanstd(ODI3_even(OSA))/sqrt(length((ODI3_even(OSA))))...
    nanstd(ODI4_even(OSA))/sqrt(length((ODI4_even(OSA))))];
%% AHI Event Duration and ODI for odd vs even windows
for i=1:2060
    rowi = 129;
    DurAr_all(i)=AHIData_all{1,i}(rowi);
    DurAr_odd(i)=AHIData_odd{1,i}(rowi);
    DurAr_even(i)=AHIData_even{1,i}(rowi);
    
    rowi = 136;
    DurEvt_all(i)=AHIData_all{1,i}(rowi);
    DurEvt_odd(i)=AHIData_odd{1,i}(rowi);
    DurEvt_even(i)=AHIData_even{1,i}(rowi);
    
    rowi = 73;
    TST_all(i)=AHIData_all{1,i}(rowi);
    TST_odd(i)=AHIData_odd{1,i}(rowi);
    TST_even(i)=AHIData_even{1,i}(rowi);
    
    rowi=74;
    ArI_all(i)=AHIData_all{1,i}(rowi);
    ArI_odd(i)=AHIData_odd{1,i}(rowi);
    ArI_even(i)=AHIData_even{1,i}(rowi);
end

if 0 %sampling rate fudge problem
    AHItot_all = AHItot_all/2;
    AHItot_odd = AHItot_odd/2;
    AHItot_even = AHItot_even/2;
    
    ODI3_odd=ODI3_odd/2;
    ODI3_even=ODI3_even/2;
    ODI3_all = ODI3_all/2;
    ODI4_odd=ODI3_odd/2;
    ODI4_even=ODI4_even/2;
    ODI4_all = ODI4_all/2;
    
    TST_all=TST_all*2;
    TST_odd=TST_odd*2;
    TST_even=TST_even*2;
    
    DurEvt_odd = DurEvt_odd*2;
    DurEvt_even = DurEvt_even*2;
    DurEvt_all = DurEvt_all*2;
    
    DurAr_odd = DurAr_odd*2;
    DurAr_even = DurAr_even*2;
    DurAr_all = DurAr_all*2;
end
MeanOdd=[MeanOdd nanmean(DurAr_odd(OSA)) nanmean(DurEvt_odd(OSA)) nanmean(ArI_odd(OSA))];
MeanOddSE=[MeanOddSE nanstd(DurAr_odd(OSA))/sqrt(length((DurAr_odd(OSA)))) nanstd(DurEvt_odd(OSA))/sqrt(length((DurEvt_odd(OSA))))...
    nanstd(ArI_odd(OSA))/sqrt(length((ArI_odd(OSA))))];

MeanEven=[MeanEven nanmean(DurAr_even(OSA)) nanmean(DurEvt_even(OSA)) nanmean(ArI_even(OSA))];
MeanEvenSE=[MeanEvenSE nanstd(DurAr_even(OSA))/sqrt(length((DurAr_even(OSA)))) nanstd(DurEvt_even(OSA))/sqrt(length((DurEvt_even(OSA))))...
    nanstd(ArI_even(OSA))/sqrt(length((ArI_even(OSA))))];

TST_all=TST_all(OSA);
TST_even=TST_even(OSA);
TST_odd=TST_odd(OSA);
figure(88)
subplot(2,3,1)
scatter(AHItot_odd(OSA),AHItot_even(OSA),10,'filled','markerfacealpha',0.5);
xlabel('AHI_ Odd'); ylabel('AHI_ Even');
subplot(2,3,2)
scatter(DurEvt_odd(OSA), DurEvt_even(OSA),10,'filled','markerfacealpha',0.5);
xlabel('EvtDur_ Odd'); ylabel('EvtDur_ Even');
subplot(2,3,3)
scatter(ArI_odd(OSA), ArI_even(OSA),10,'filled','markerfacealpha',0.5);
xlabel('ArI_ Odd'); ylabel('ArI_ Even');
subplot(2,3,4)
scatter(DurAr_odd(OSA), DurAr_even(OSA),10,'filled','markerfacealpha',0.5);
xlabel('ArDur_ Odd'); ylabel('ArDur_ Even');
subplot(2,3,5)
scatter(ODI3_odd(OSA),ODI3_even(OSA),10,'filled','markerfacealpha',0.5);
xlabel('ODI3_ Odd'); ylabel('ODI3_ Even');
subplot(2,3,6)
scatter(ODI4_odd(OSA),ODI4_even(OSA),10,'filled','markerfacealpha',0.5);
xlabel('ODI4_ Odd'); ylabel('ODI4_ Even');

%% if 0
temp=SummaryAnalysisTable_1.Vpassive;
temp(end+1:2060)=NaN;
temp(temp<3)=NaN;
I1=~isnan(temp);
temp=SummaryAnalysisTable_2.Vpassive;
temp(end+1:2060)=NaN;
temp(temp<3)=NaN;
I2=~isnan(temp);
I = I1&I2 & OSA;

MESAdata.eventdur=DurEvt_all(:);
MESAdata.arousaldur=DurAr_all(:);
MESAdata.ODI3=ODI3_all(:);
MESAdata.ODI4=ODI4_all(:);

Table1List = {'age5c','gender1','bmi5c',...
    'a0h3ai5','a0h3aiNREM','a0h3aiREM',...
    'ai_all5','eventdur','arousaldur','ODI3','ODI4',...
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
    'eventdur','arousaldur','ODI3','ODI4',...
    'slpprdp5','N1Per','N2Per','N34Per','RemPer','supinep5'}; 


tempfactorx=[1 1 ...
            0.25 0.25 0.75 0.25 ... 
            1 1 0.25 0.25 ...
            1 0.5 1 0.5 1 1];
logconv=[0 0 ...
        0 0 0 0 ...
        0 0 0 0 ...
        0 0 0 0 0 0];
       
nonparametric = [0 0 ...
        0 0 0 0 ...
        0 0 0 0  ...
        0 0 0 0 0 1];
      

Table1mean=[];
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

%%
clear MSEsqrt MSEwholeOSA MSEwhole
figure(11); clf(11);
traitname={'AHItot','DurEvt','ArI','DurAr','ODI3','ODI4'};
traitmatrix = [AHItot_all(:) DurEvt_all(:) ArI_all(:) DurAr_all(:) ODI3_all(:) ODI4_all(:)];

J=6;
for j=1:J
    switch j
        case 1
            traitodd = AHItot_odd(OSA);
            traiteven = AHItot_even(OSA);
            traitall = AHItot_all(OSA);
        case 2
            traitodd = DurEvt_odd(:);
            traitodd=traitodd(OSA);
            traiteven = DurEvt_even(:);
            traiteven=traiteven(OSA);
            traitall = DurEvt_all(:);
            traitall=traitall(OSA);
        case 3
            traitodd = ArI_odd(OSA);
            traiteven = ArI_even(OSA);
            traitall = ArI_all(OSA);
        case 4
            traitodd = DurAr_odd(:);
            traitodd=traitodd(OSA);
            traiteven = DurAr_even(:);
             traiteven=traiteven(OSA);
            traitall = DurAr_all(:);
             traitall=traitall(OSA);
        case 5
            traitodd = ODI3_odd(:);
            traitodd=traitodd(OSA);
            traiteven = ODI3_even(:);
             traiteven=traiteven(OSA);
            traitall = ODI3_all(:);
             traitall=traitall(OSA);
        case 6
            traitodd = ODI4_odd(:);
            traitodd=traitodd(OSA);
            traiteven = ODI4_even(:);
              traiteven=traiteven(OSA);
            traitall = ODI4_all(:);
             traitall=traitall(OSA);
    end
    if size(traitodd,1)==1,traitodd=traitodd';end
    if size(traiteven,1)==1,traiteven=traiteven';end
    if size(traitall,1)==1,traitall=traitall';end
    
    TST_avg = ( TST_odd' + TST_even')./2; % avg tst
    deltatrait=traiteven-traitodd; % change in trait between odd and even windows
    DeltatTraitMean(j)=nanmean(deltatrait);
    DeltaTraitMeanSE(j)=[nanstd(deltatrait)/sqrt(length(deltatrait))];
    meanTrait = (traitodd+traiteven)/2; % avg trait value
    
    Errsq_trait = (traitodd -traiteven).^2; % squared error between odd and even windows
    Errsq_trait((Errsq_trait<0.0001*nanmedian(Errsq_trait)))=0.0001*nanmedian(Errsq_trait); % accounting for error due to very small values
    
    mdlErr = fitglm(Errsq_trait,TST_avg,'Link','log') % model the error
    
    mdlErr = fitglm([TST_avg abs(meanTrait)],Errsq_trait,'Link','log')
    Nline = (0:360)'; % for tst
    ErrsqModel = predict(mdlErr,[Nline Nline*0+nanmean(abs(meanTrait))]);
    MStot = nanmean(((traiteven - nanmean(traiteven)).^2)); %simply variance
    
    Rtemp(j) = (1 - nanmean(Errsq_trait)/MStot).^0.5;
    
    Rtemp2(j) = corr(traitodd,traiteven,'Rows','complete')
    
   % temp_as	= nItem*r/(1 + (nItem-1)*r);
    
    subplot(2,J,j)
    scatter(traitodd,traiteven,10,'filled','markerfacealpha',0.5);
    xname=strcat(traitname{j},'_ Odd'); yname=strcat(traitname{j},'_ Even');
    xlabel(xname); ylabel(yname);
    
    subplot(2,J,J+j)
    exponent=0.5;
    plot(TST_avg,(Errsq_trait/MStot).^exponent,'.');
    hold('on');
    plot(Nline,(ErrsqModel/MStot).^exponent,'k','linewidth',2);
    ylim([0 2]);
    yticks=get(gca,'ytick');
    set(gca,'ytick',yticks,'yticklabels',yticks.^(1/exponent));
    box('off');
    xlabel('TST(min)');
    
    ErrsqDataN = predict(mdlErr,[TST_all(:) traitall(:)]); % predict error for the tst in the entire data set
    
    MSEwhole(j) = nanmean(ErrsqDataN);
%     MSEwholeOSA(j) = nanmean(ErrsqDataN(OSA));
    MSEsqrt(j)=MSEwhole(j).^0.5; % error of the trait for data set
    x=[traitodd traiteven];
%     [au(j),as(j)] = CronbachAlpha(x)
    
end

%%
traitsoddeven=[1:6]
MSEwhole_ = MSEwhole.^0.5
%MSEwholeOSA_ = MSEwholeOSA.^0.5
temp = traitmatrix;
OSAstd = nanstd(temp(OSA,traitsoddeven))
Allstd = nanstd(temp(:,traitsoddeven))

% FstdOSA = MSEwholeOSA_ ./ OSAstd
% FvarOSA = MSEwholeOSA_.^2 ./ OSAstd.^2

FvarAll = MSEwhole_.^2 ./ Allstd.^2
FvarAllx = MSEwhole_.^2 ./ OSAstd.^2
traitname


%%
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


%%
%comb filter within-night reliability
figure(2); clf(2);
figure(1); clf(1);
figure(3); clf(3);


% medianwhole = [53    53    53    40    40    78    78    78] % median # of possible windows based on entire MESA data set
% lower25thwhole = [27   27   27   21   21   44   44   44] % 25th percentile for # of windows

thres=3*ones(length(traitsj)); % minimum 3 windows needed to be considered

%% NaNs in N should be zero

failedstudy = isnan(AHI);

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

%%
clear MSE x temp thresgood thresgood1
 clf (figure(88));
 traitsj=[1,2,5,6,8];
for i=1:length(traitsj)
    figure(1);
    subplot(2,length(traitsj),i)
    j=traitsj(i);
    N1=SummaryAnalysisTableN_1{:,j}; % comb=1=odd windows
    N2=SummaryAnalysisTableN_2{:,j}; % comb=2=even windows
    crit = N1>=thres(i) & N2>=thres(i)&OSA; % use only those subjects who has windows greater than threshold
    Nsubj(:,i)=crit;
    Nactual(1,i)=nanmedian(N1(crit)); % median # of odd windows in the actual file
    Nactual(2,i)=nanmedian(N2(crit)); % median # of even windows in the actual file
    
    deltaFN1 = FstatesArray_2{j}(crit,1) - FstatesArray_1{j}(crit,1); % change in n1
    deltaFN3 = FstatesArray_2{j}(crit,3) - FstatesArray_1{j}(crit,3); % change in n3
    deltaFW = (1-sum(FstatesArray_2{j}(crit,:),2)) - (1-sum(FstatesArray_1{j}(crit,:),2)); % change in wake
    deltaFsup = Fsupine_2(crit,j) - Fsupine_1(crit,j); % change in supine position
    Trait2 = SummaryAnalysisTable_2{crit==1,j}; % trait from even window
    Trait1 = SummaryAnalysisTable_1{crit==1,j}; % trait from odd window
        
    if j==5 % transform ArTh
        Trait1(Trait1<100)=100;
        Trait1 = 100*(1+((Trait1/100-1).^0.5));
        Trait2(Trait2<100)=100;
        Trait2 = 100*(1+((Trait2/100-1).^0.5));
    end
    
    if j==6 % transform Vpassive
        Trait1(Trait1>100)=100;
        Trait1 = 100*(1-((1-Trait1/100).^0.5));
        Trait2(Trait2>100)=100;
        Trait2 = 100*(1-((1-Trait2/100).^0.5));
    end
    
%     % CronbachAlpha
%     x=[Trait1 Trait2];
%     [autrait(i),astrait(i)] = CronbachAlpha(x);
    
    TraitMeanOdd(j)=nanmean(Trait1);
    TraitMeanOddSE(j)=[nanstd(Trait1)/sqrt(length(Trait1))];
    TraitMeanEven(j)=nanmean(Trait2);
    TraitMeanEvenSE(j)=[nanstd(Trait2)/sqrt(length(Trait2))];
    TraitDel(j)=nanmean(Trait2-Trait1);
    TraitDelSE(j)=[nanstd(Trait2-Trait1)/sqrt(length(Trait2-Trait1))];
    
    RTrait2(j) = corr(Trait1,Trait2,'Rows','complete')
    
    SSE = nansum(((Trait1 - Trait2).^2)); % sum squared error
    SStot = nansum(((Trait1 - nanmean(Trait1)).^2));
    RTrait2_(j)
    
figure(88);
if j==1
subplot(2,3,1)
scatter(Trait1,Trait2,10,'filled','markerfacealpha',0.5);
xlabel('Loop Gain,LG1(odd)'); ylabel('Loop Gain,LG1(even)');
elseif j==2
subplot(2,3,2)
scatter(Trait1,Trait2,10,'filled','markerfacealpha',0.5);
xlabel('Ventilatory Instability(odd)'); ylabel('Ventilatory Instability(even)');
elseif j==5
subplot(2,3,3)
scatter(Trait1,Trait2,10,'filled','markerfacealpha',0.5);
xlabel('Arousal Threshold(odd)'); ylabel('Arousal Threshold(even)');
elseif j==6
subplot(2,3,4)
scatter(Trait1,Trait2,10,'filled','markerfacealpha',0.5);
xlabel('Collapsibility(odd)'); ylabel('Collapsibility(even)');
elseif j==7
subplot(2,3,5)
scatter(Trait1,Trait2,10,'filled','markerfacealpha',0.5);
xlabel('Vactive(odd)'); ylabel('Vactive(even)');
elseif j==8
subplot(2,3,6)
scatter(Trait1,Trait2,10,'filled','markerfacealpha',0.5);
xlabel('Compensation(odd)'); ylabel('Compensation(even)');
end

%%
    %% uncomment for model % commented out for getting plots and corr coeffs
    if 1
    % Model Trait 2 as a function of trait 1, N1,N3,Wake and supine
    T = table(Trait1,Trait2,deltaFW,deltaFN3,deltaFN1,deltaFsup);
    mdl = fitglm(T,'Trait2 ~ Trait1 + deltaFN3 + deltaFN1 + deltaFW + deltaFsup'); % linear regression model of Trait 2
    Trait1adj = (predict(mdl,T) - mdl.Coefficients.Estimate(1))/mdl.Coefficients.Estimate(2); % adjusting for trait 1
    
    % Calculate MSE of Trait 1 vs 2
    MSE(1,i) = nanmean((Trait1 - Trait2).^2); % mse of trait 1 vs 2
    MSE(2,i) = nanmean((Trait1adj - Trait2).^2); % mse after adjusting for trait 1
    
    % %for plotting Trait 1 vs 2 based on percentiles
    
    %     temp = w.DataNArray(:,j); % # of windows in the MESA-all sub mat file
    %     thres2 = prctile(temp(temp>=3),10) % lowest threshold-10th percentile of all subjects with # windows greater than 3.
    %     thres3 = prctile(temp(temp>=3),33) % 33 percentile
    
    %     if 0
    %         Trait1 = Trait1adj;
    %     end
    % figure(1)
    % sum(crit)
    % %For Trait 1, plot trait 1 and 2 with lowest and highest percentiles. Blue dots
    % %represents the ones between 10th and 33 rd percentiles.
    %
    % %scatter(Trait1,Trait2,10,'filled','markerfacealpha',0.5);
    %
    % scatter(Trait1(N1(crit)<thres2 | N2(crit)<thres2),Trait2(N1(crit)<thres2 | N2(crit)<thres2),10,[0.91 0.41 0.17],'filled','markerfacealpha',0.2);
    % hold('on');
    % %scatter(Trait1(N1(crit)>=thres2 & N2(crit)>=thres2),Trait2(N1(crit)>=thres2 & N2(crit)>=thres2),10,[0.1 0.1 0.8],'filled','markerfacealpha',0.6);
    % scatter(Trait1(N1(crit)>=thres2 & N2(crit)>=thres2 & ~(N1(crit)>=thres3 & N2(crit)>=thres3)),Trait2(N1(crit)>=thres2 & N2(crit)>=thres2 & ~(N1(crit)>=thres3 & N2(crit)>=thres3)),10,[0.1 0.1 0.8],'filled','markerfacealpha',0.5);
    % scatter(Trait1(N1(crit)>=thres3 & N2(crit)>=thres3),Trait2(N1(crit)>=thres3 & N2(crit)>=thres3),10,[0.1 0.8 0.8],'filled','markerfacealpha',0.8);
    %
    % %For Trait 1 adjusted , plot trait 1 and 2 with lowest and highest percentiles. Blue dots
    % %represents the ones between 10th and 33 rd percentiles.
    %
    % ylims = get(gca,'ylim');
    % set(gca,'xlim',ylims);
    % yticks = get(gca,'ytick');
    % set(gca,'xtick',yticks);
    
    % figure(2)
    % subplot(2,length(traitsj),i)
    %
    % scatter(Trait1adj(N1(crit)<thres2 | N2(crit)<thres2),Trait2(N1(crit)<thres2 | N2(crit)<thres2),10,[0.91 0.41 0.17],'filled','markerfacealpha',0.2);
    % hold('on');
    % scatter(Trait1adj(N1(crit)>=thres2 & N2(crit)>=thres2),Trait2(N1(crit)>=thres2 & N2(crit)>=thres2),10,[0.1 0.1 0.8],'filled','markerfacealpha',0.6);
    % scatter(Trait1adj(N1(crit)>=thres3 & N2(crit)>=thres3),Trait2(N1(crit)>=thres3 & N2(crit)>=thres3),10,[0.1 0.8 0.8],'filled','markerfacealpha',0.8);
    
    figure(1);
    subplot(4,length(traitsj),i + 2*length(traitsj))
    
    % find the error between odd and even traits
    MStot = nanmean(((Trait2 - nanmean(Trait2)).^2)); %simply variance
    Errsq = ((Trait1 - Trait2).^2); % squared error between odd and even windows
    Errsq((Errsq<0.0001*nanmedian(Errsq)))=0.0001*nanmedian(Errsq);
    
    N_ = (N1(crit) + N2(crit))./2; % avg # of windows
    meanTrait = (Trait1 + Trait2)/2;
    
    %beta standardized
%     N_=(N_-nanmean(N_))./nanstd(N_);
%     meanTrait=(meanTrait-nanmean(meanTrait))./nanstd(meanTrait);
%     Errsq=(Errsq-nanmean(Errsq))./nanstd(Errsq);
    
    % plotting error for OSA vs Controls
    tempy=(Errsq(OSA(crit))/MStot).^0.5;
    tempy(tempy>2)=2;
    plot(N_(OSA(crit)),tempy,'.'); % plotting OSA in blue color
    hold('on');
    tempy=(Errsq(~OSA(crit))/MStot).^0.5;
    tempy(tempy>2)=2;
    plot(N_(~OSA(crit)),tempy,'.'); % controls in orange color
    hold('on');
    
    % Model the error for odd and even traits
    Nline = (0:200)'; % # windows for prediction;
    mdlErr = fitglm(N_,Errsq,'Link','log');
    mdlErrKeep2{i}=mdlErr;
   
%     mdlErr = fitglm([N_ abs(meanTrait)],Errsq,'Link','log'); % original
    

    % predict error for number of windows in "Nline"
    useconstant=1;
    if ~useconstant
        mdlErr = fitglm([N_ meanTrait],Errsq,'Link','log');
        mdlErrKeep{i}=mdlErr;
        ErrsqModel = predict(mdlErr,[Nline Nline*0+nanmean(abs(meanTrait))]);


        % predict error for number of windows in "Npoints"
        Npoints = [10 60]; % for window length of 10 and 60
        ErrModel(1,i) = ErrsqModel(find(Nline==Npoints(1)));
        ErrModel(2,i) = ErrsqModel(find(Nline==Npoints(2)));


        %predict total err SSE (adjusted) for number of windows in whole data set
        ErrsqDataN = predict(mdlErr,[w.DataNArray(:,j) abs(w.DataArray(:,j))]); % predict error for the total #windows present in the entire data set
        ExclBelow3 = w.DataNArray(:,j)<3 | isnan(w.DataNArray(:,j)); % exclude subjects with window # below 3 or NaN
        ErrsqDataN(ExclBelow3)=NaN;

        % MSE of error in the entire cohort
        MSEwhole(1,i) = nanmean(ErrsqDataN); 
    else %model a flat line just to check R=R
         mdlErr = fitglm([N_*0+1],Errsq,'Link','log');
         mdlErrKeep{i}=mdlErr;
        ErrsqModel = predict(mdlErr,[Nline*0+1]);


        % predict error for number of windows in "Npoints"
        Npoints = [10 60]; % for window length of 10 and 60
        ErrModel(1,i) = ErrsqModel(find(Nline==Npoints(1)));
        ErrModel(2,i) = ErrsqModel(find(Nline==Npoints(2)));


        %predict total err SSE (adjusted) for number of windows in whole data set
        ErrsqDataN = predict(mdlErr,[N_*0+1]); % predict error for the total #windows present in the entire data set
        ExclBelow3 = w.DataNArray(:,j)<3 | isnan(w.DataNArray(:,j)); % exclude subjects with window # below 3 or NaN
        ErrsqDataN(ExclBelow3)=NaN;

        % MSE of error in the entire cohort
        MSEwhole(1,i) = nanmean(ErrsqDataN); 
    end
    
    % repeat for OSA only; analysis is **unadjusted**
    if exist('temp'); clear temp; end
    if exist('temp2'); clear temp2; end
    temp = w.DataNArray(:,j);
    temp = temp(AHI>OSAAHIthres);
    temp2 = abs(w.DataArray(:,j));
    temp2 = temp2(AHI>OSAAHIthres);
    if ~useconstant
        ErrsqDataOSAN = predict(mdlErr,[temp temp2]);
    else
        ErrsqDataOSAN = predict(mdlErr,[temp*0]);
    end
    
    % MSE of error for OSA only
    ExclBelow3 = temp<3 | isnan(temp); % exclude sub with >3 windows
    ErrsqDataOSAN(ExclBelow3)=NaN;
    MSEwholeOSA(1,i) = nanmean(ErrsqDataOSAN);
  
    %finding the number of sub who had enough windows for the analysis
    rsqthresgood = 0.75;
    try
        thresgood1(i)=Nline(find(ErrsqModel/MStot<(1-rsqthresgood),1));
    catch me
        thresgood1(i)=999;
        'no data above rsqthres'
    end
    if thresgood1(i)>46
        thresgood1(i)=46;
    end
    
    %threshold based on percentile of windows in whole cohort.
    % currently not used
    if 0
    thresgood2(i) = prctile(w.DataNArray(:,j),33);
    thresgood(i) = thresgood2(i);
    end
    
    thresgood(i) = thresgood1(i);
    
    % Fraction of OSA with rsquared >0.75
    Ngood(i) = sum(w.DataNArray(:,j)>=thresgood(i));
    NgoodOSA(i) = sum(OSA&~(isnan(AHI)) & w.DataNArray(:,j)>=thresgood(i));
    FgoodOSA(i) = NgoodOSA(i)/sum(OSA&~(isnan(AHI)));
    
    % fraction of controls
    NgoodNOSA(i) = sum(~OSA&~(isnan(AHI)) & w.DataNArray(:,j)>=thresgood(i));
    FgoodNOSA(i) = NgoodNOSA(i)/sum(~OSA&~(isnan(AHI)));
    
    % plot the 1-rsq2 line on OSA vs control plot
    tempy=(ErrsqModel/MStot).^0.5;
    tempy(tempy>2)=2;
    plot(Nline,tempy,'k','linewidth',2);
    ylim([0 2]);
    yticks=get(gca,'ytick');
    set(gca,'ytick',yticks,'yticklabels',yticks.^2);
    xlim([0 upperNtrait(traitsj(i))]);
    box('off')
    
    % plot the OSA vs control for trait 1 vs 2
    figure(1)
    subplot(4,length(traitsj),i);
    hold('off');
    
    % for subjects with >rsq 0.75
    scatter(Trait1(N1(crit)>thresgood(i)&OSA(crit)),Trait2(N1(crit)>thresgood(i)&OSA(crit)),10,[0 0.45 0.74],'filled','markerfacealpha',0.5);
    hold('on');
    scatter(Trait1(N1(crit)>thresgood(i)&~OSA(crit)),Trait2(N1(crit)>thresgood(i)&~OSA(crit)),10,[216 82 24]/256,'filled','markerfacealpha',0.5);
    ylims = get(gca,'ylim');
    set(gca,'xlim',ylims);
    yticks = get(gca,'ytick');
    set(gca,'xtick',yticks);
    
    % for subjects with <rsq 0.75
    subplot(4,length(traitsj),i + length(traitsj)*1);
    hold('off');
    scatter(Trait1(N1(crit)<thresgood(i)&OSA(crit)),Trait2(N1(crit)<thresgood(i)&OSA(crit)),5,[0 0.45 0.74],'filled','markerfacealpha',0.5);
    hold('on');
    scatter(Trait1(N1(crit)<thresgood(i)&~OSA(crit)),Trait2(N1(crit)<thresgood(i)&~OSA(crit)),5,[216 82 24]/256,'filled','markerfacealpha',0.5);
    ylims = get(gca,'ylim');
    set(gca,'xlim',ylims);
    yticks = get(gca,'ytick');
    set(gca,'xtick',yticks);
    
   
    % repeat the above steps for the Trait 1 adjusted...
    if 0
    figure(2);
    subplot(2,length(traitsj),i + length(traitsj))
    
    Errsq = ((Trait1adj - Trait2).^2);
    N_ = (N1(crit) + N2(crit))./2;
    Errsq((Errsq<0.0001*nanmedian(Errsq)))=0.0001*nanmedian(Errsq);
    %plot(N_,log10(Errsq + prctile(Errsq,3)),'.');
    plot(N_,Errsq/MStot,'.');
    hold('on');
        
    mdlErr = fitglm(N_,Errsq,'Link','log')
    meanTrait = (Trait1adj + Trait2)/2;
    mdlErr = fitglm([N_ abs(meanTrait)],Errsq,'Link','log')
    
    ErrsqModel = predict(mdlErr,[Nline Nline*0+nanmean(abs(meanTrait))]); 
    ErrModel(3,i) = ErrsqModel(find(Nline==Npoints(1)));
    ErrModel(4,i) = ErrsqModel(find(Nline==Npoints(2)));
        
    plot(Nline,ErrsqModel/MStot,'r');
    ylim([0 4*max(ErrsqModel)/MStot]);
    ylim([0 1]);
    xlim([0 upperNtrait(traitsj(i))]);
    box('off');
    end 
    
    % Rsqrd for trait 1 vs trait 2
    SSE = nansum(((Trait1 - Trait2).^2)); % sum squared error
    SStot = nansum(((Trait1 - nanmean(Trait1)).^2));
    Rsq1(i) = 1 - SSE/SStot;
    
      
    % RsqEq for ErrModel based on number of windows in "Npoints"
    RsqEq(1,i) = 1 - ErrModel(1,i)/MStot;
    RsqEq(2,i) = 1 - ErrModel(2,i)/MStot;
    
    ErrSDRelative(1,i) = ErrModel(1,i).^0.5/MStot.^0.5;
    ErrSDRelative(2,i) = ErrModel(2,i).^0.5/MStot.^0.5;
    
    %% Histogram for OSA vs Controls
    % OSA blue; controls orange 
    figure(1);
    subplot(4,length(traitsj),i + length(traitsj)*3);
    
    dStep=5;
    Centers=0:dStep:upperNtrait(traitsj(i));
    Edges=(Centers(1)-dStep/2):dStep:(Centers(end)+dStep/2);
    Edges(end)=Inf;
    temp = w.DataNArray(:,j);
    temp = temp(AHI>OSAAHIthres);
    [h7,edges] = histcounts(temp,Edges);
    temp = w.DataNArray(:,j);
    temp = temp(AHI<=OSAAHIthres);
    [h8,edges] = histcounts(temp,Edges);
    bar(Centers,h7,'EdgeAlpha',0,'FaceAlpha',0.3,'BarWidth',1);
    hold('on');
    bar(Centers,h8,'EdgeAlpha',0,'FaceAlpha',0.7,'BarWidth',1);
    box('off')
end
end

ErrModel_ = ErrModel.^0.5 % ErrModel is for N=10 and 60 windows

MSEwhole_ = MSEwhole.^0.5
MSEwholeOSA_ = MSEwholeOSA.^0.5
temp = w.DataArray;
OSAstd = nanstd(temp(OSA,traitsj))

temp2 = w.DataNArray;
temp(temp2<3)=NaN;
Allstd = nanstd(temp(:,traitsj))

FstdOSA = MSEwholeOSA_ ./ OSAstd
FvarOSA = MSEwholeOSA_.^2 ./ OSAstd.^2
(1 - FvarOSA)
(1 - FvarOSA).^0.5

RTrait2(traitsj)
Rsq1.^0.5

FvarAll = MSEwhole_.^2 ./ Allstd.^2
FvarAllx = MSEwhole_.^2 ./ OSAstd.^2

(1 - FvarAll).^0.5
(1 - FvarAllx).^0.5

for i=1:length(mdlErrKeep)
    i
    mdlErrKeep{i}
    mdlErrKeep2{i};
end

figure(1);
subplot(4,4,1);
axis([0 1 0 1])
subplot(4,4,8);
axis([-100 200 -100 200])

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
%% make graphs look pretty
figure(1);
set(gcf,'color',[1 1 1])





