clear all; clc; close all;
%% BootTAnalysis
% load BootT table
% remove those with AHI<5
%% Options
OSAAHIthres=5;% 
critwin=3;
addpath(genpath('C:\Users\rma56\Dropbox (Partners HealthCare)\PUPbeta_git\PUPbeta\'));

%% load mros mesa
% to do : try merging both  mros and mesa

% mros visit 1
    % Boot data
    mrosendotypes=load('C:\Users\rma56\Dropbox (Partners HealthCare)\Reliability Trait Analysis\MrOS\Summary\SummaryAnalysis_V1_AllSleep_AllPos_BootT_Oct092021.mat');
    mrosahi=load('C:\Users\rma56\Dropbox (Partners HealthCare)\Reliability Trait Analysis\MrOS\Summary\AhiTstBootT_Oct2021.mat');
        
    % AHI Data
    mrosNsrr1=readtable('C:\Users\rma56\Dropbox (Partners HealthCare)\Reliability Trait Analysis\MrOS\NSRR Updated AHI\mros1_nsrr_new_ahis_20200601.csv');
    mrosAhi3pAr=mrosNsrr1.POOAHI3;
    

% mesa
    % Boot data
    mesaendotypes=load('C:\Users\rma56\Dropbox (Partners HealthCare)\Reliability Trait Analysis\MESA\Reliability\SummaryAnalysis_AllSleep_AllPos_BootT_Oct082021.mat');
    mesaahi=load('C:\Users\rma56\Dropbox (Partners HealthCare)\Reliability Trait Analysis\MESA\Reliability\AhiTstBootT_Oct2021.mat');
% AHI data
    MESAdata=load('C:\Users\rma56\Dropbox (Partners HealthCare)\Reliability Trait Analysis\MESA\MESAdata.mat');
    mesaAhi3pAr =MESAdata.MESAdata.a0h3ai5;
    subid=MESAdata.MESAdata.idno;
    % subid from mesa spreadsheet
    excelSubId=readtable('C:\Users\rma56\Dropbox (Partners HealthCare)\Reliability Trait Analysis\MESA\AMasterSpreadsheet.xlsx','Range','T3:T2063');
    Subiddiff=subid-table2array(excelSubId);
    find(Subiddiff>0) % sub# 40 in excelSubId is repetition
    %code to remove sub# 40??
    I = find(Subiddiff>0)';
    
% end

% length(unique(BootT.Subject))

%% create combined table
% endotypes
mrosbootT=mrosendotypes.BootT;
mrosbootT.Study=ones(height(mrosbootT),1);
mrosbootT=movevars(mrosbootT,'Study','After','Subject');

mesabootT=mesaendotypes.BootT;
mesabootT.Subject=mesabootT.Subject+2906;
mesabootT.Study=2*ones(height(mesabootT),1);
mesabootT=movevars(mesabootT,'Study','After','Subject');

BootT=[mrosbootT;mesabootT];

clear mrosbootT mesabootT


% psg parameters
mrosbootT=mrosahi.BootTPSGParam;
mrosbootT.Study=ones(height(mrosbootT),1);
mrosbootT=movevars(mrosbootT,'Study','After','Subject');

mesabootT=mesaahi.BootTPSGParam;
mesabootT.Subject=mesabootT.Subject+2906;
mesabootT.Study=2*ones(height(mesabootT),1);
mesabootT=movevars(mesabootT,'Study','After','Subject');

BootTPSGParam=[mrosbootT;mesabootT];

clear mrosbootT mesabootT

mrosDataNArray=[mrosendotypes.DataNArray,ones(length(mrosendotypes.DataNArray),1)];
mesaDataNArray=[mesaendotypes.DataNArray,2*ones(length(mesaendotypes.DataNArray),1)];
DataNArray=[mrosDataNArray;mesaDataNArray];

if 1
    MismatchT=readtable('C:\Users\rma56\Dropbox (Partners HealthCare)\Reliability Trait Analysis\MrOS\nsrr mismatch.xlsx');
    NSRRID=string(MismatchT.NSRRID);
    PUPID=string(upper(extractAfter(MismatchT.pupstart,'visit1-')));
    MismatchT.Compare=zeros(height(MismatchT),1);
    for ii=1:length(NSRRID)
        if sum(strcmp(NSRRID{ii},PUPID))>0
            MismatchT.Compare(ii)=1;
        end
    end
    
    mrosAhi3pAr(MismatchT.Compare==0)=[];
end

Ahi3pAr=[mrosAhi3pAr;mesaAhi3pAr];

%%
clear betaAll betaStdAll
if 1
    OSA=Ahi3pAr>OSAAHIthres;
    OSAMesa=mesaAhi3pAr>OSAAHIthres;
    sum(OSAMesa==1)
    DataNtempMesa=sum(mesaDataNArray(:,1:9)<critwin);
   
    OSAMrOS=mrosAhi3pAr>OSAAHIthres;
    sum(OSAMrOS==1)
    DataNtempMrOS=sum(mrosDataNArray(:,1:9)<critwin);
    
    %code to limit data based on min number of windows
    DataNtemp=BootT{:,10:16}<critwin; % missing ArTh Num in MrOs bootT table
    DataNtemp = [DataNtemp(:,1:4) DataNtemp(:,4) DataNtemp(:,5:7)];
    temp = BootT{:,2:9};
    temp(DataNtemp)=NaN;
    BootT{:,2:9} = temp;
    
    I = find(~OSA)';
    BootTPSGParam(any(I(:)==BootTPSGParam.Subject'),:)=[]; %Removing nonOSA
    BootT(any(I(:)==BootT.Subject'),:)=[]; %Removing nonOSA
    
    Irem=find(DataNArray(:,6)<critwin);
    BootTPSGParam(any(Irem==BootTPSGParam.Subject'),:)=[]; %Removing no endotype data per critwin
    BootT(any(Irem==BootT.Subject'),:)=[];
    
    %code to remove those with low AHI here
    
end

% NperTrait = sum(DataNArray>critwin & OSA);
% NPSGdata= length(unique(BootTPSGParam.Subject))

BootT.FlatLoopGain = 1 - BootT.FsupLoopGain;
BootT.FlatArTh = 1 - BootT.FsupArTh;
BootT.FlatUA = 1 - BootT.FsupUA;

BootT.ArThresBootT = fArThresT(BootT.ArThresBoot);
BootT.VpassiveBootT = fVpassiveT(BootT.VpassiveBoot);

subjs=unique(BootT.Subject);
MeanT=NaN(length(subjs),width(BootT));
for ii=1:length(subjs)
    ii
    MeanT(ii,:)=nanmedian(BootT{BootT.Subject==subjs(ii),:});
end

stds1 = nanstd(MeanT);
stdT=array2table(stds1);
stdT.Properties.VariableNames=BootT.Properties.VariableNames

% standardized beta coeffs - 
% calculated by subtracting the mean from the variable and dividing by its standard deviation
BootTStd=BootT;
for ii=[3:10 36:37]
    BootTStd{:,ii}=(BootTStd{:,ii}-nanmean(BootTStd{:,ii}))./stds1(:,ii);
end

for ii=[18:35]
    BootTStd{:,ii}=BootTStd{:,ii}*100;
end


BootTPSGParam.AHIT = BootTPSGParam.AHI.^0.5;
BootTPSGParam.ArIT = BootTPSGParam.ArI.^0.5;

subjs=unique(BootTPSGParam.Subject);
MeanT2=NaN(length(subjs),width(BootTPSGParam));
for ii=1:length(subjs)
    ii
    MeanT2(ii,:)=nanmedian(BootTPSGParam{BootTPSGParam.Subject==subjs(ii),:});
end
stds2 = nanstd(MeanT2);
std2T=array2table(stds2);
std2T.Properties.VariableNames=BootTPSGParam.Properties.VariableNames


BootTPSGParamStd=BootTPSGParam;
for ii=[8:11]
    BootTPSGParamStd{:,ii}=(BootTPSGParam{:,ii}-nanmean(BootTPSGParam{:,ii}))./stds2(:,ii);
end
for ii=[3:7]
    BootTPSGParamStd{:,ii}=BootTPSGParamStd{:,ii}*100;
end

%%
yvarlist = {'VpassiveBootT','VcompBoot','LG1Boot','LGnBoot','ArThresBootT'}

clear betaStdAll SEStdAll pValueAll
for j=1:length(yvarlist)

mdl2 = fitglme(BootTStd,[yvarlist{j} '~ FN1UA + FN3UA + FREMUA + FlatUA + (1|Subject)']);
betaStdAll(j,:) = mdl2.Coefficients.Estimate(2:end);
SEStdAll(j,:)= mdl2.Coefficients.SE(2:end);
pValueAll(j,:) = mdl2.Coefficients.pValue(2:end);


end

%%
try betaStdAll=betaStdAll{:,:}; end
try SEStdAll=SEStdAll{:,:}; end
try pValueAll=pValueAll{:,:}; end

j=6
mdl2 = fitglme(BootTPSGParamStd,'AHIT ~ FN1 + FN3 + FREM + FLateral + (1|Subject)');
betaStdAll(j,:) = mdl2.Coefficients.Estimate(2:end);
SEStdAll(j,:)= mdl2.Coefficients.SE(2:end);
pValueAll(j,:) = mdl2.Coefficients.pValue(2:end);



betaNames = {'FN1','FN3','FREM','Flateral'};
varlist={'Vpassive','Vcomp','LG1','LGn','ArTh','AHI'};
betaStdAll = array2table(betaStdAll*100);
betaStdAll.Properties.RowNames = varlist;
betaStdAll.Properties.VariableNames = betaNames;
betaStdAll.Properties.RowNames{5}='ArThresT';
betaStdAll.Properties.RowNames{1}='VpassiveT';
betaAllTpart = betaStdAll(:,[4 3 1 2])

SEStdAll = array2table(SEStdAll*100);
SEStdAll.Properties.VariableNames=betaStdAll.Properties.VariableNames;
SEStdAll.Properties.RowNames=betaStdAll.Properties.RowNames;
SEAllTpart = SEStdAll(:,[4 3 1 2])

pValueAll = array2table(pValueAll);
pValueAll.Properties.VariableNames=betaStdAll.Properties.VariableNames;
pValueAll.Properties.RowNames=betaStdAll.Properties.RowNames;
pValueAllTpart = pValueAll(:,[4 3 1 2])


if 0
betaStdAll = array2table(betaStdAll);
betaStdAll.Properties.VariableNames=betaAll.Properties.VariableNames;
betaStdAll.Properties.RowNames = betaNames;
end



%% Apply to the data
cDataT = array2table(DataArray);
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


%% unused
mdl = fitglme(BootT,'delayBoot ~ FN1LoopGain + FN3LoopGain + FREMLoopGain + FlatLoopGain + (1|Subject)')
betaAll(:,j) = mdl.Coefficients.Estimate(2:end)
SEAll(:,j)= mdl.Coefficients.SE(2:end);
mdl2 = fitglme(BootTStd,'delayBoot ~ FN1LoopGain + FN3LoopGain + FREMLoopGain + FlatLoopGain + (1|Subject)')
betaStdAll(:,j) = mdl.Coefficients.Estimate(2:end)./nanstd(BootT.delayBoot)
SEStdAll(:,j)= mdl.Coefficients.SE(2:end)./nanstd(BootT.delayBoot);

mdl = fitglme(BootT,'VRABoot ~ FN1ArTh + FN3ArTh + FREMArTh + FlatArTh + (1|Subject)')
betaAll(:,j) = mdl.Coefficients.Estimate(2:end)
SEAll(:,j)= mdl.Coefficients.SE(2:end);
mdl2 = fitglme(BootTStd,'VRABoot ~ FN1ArTh + FN3ArTh + FREMArTh + FlatArTh + (1|Subject)')
betaStdAll(:,j) = mdl.Coefficients.Estimate(2:end)./nanstd(BootT.VRABoot)
SEStdAll(:,j)= mdl.Coefficients.SE(2:end)./nanstd(BootT.VRABoot);

mdl = fitglme(BootT,'VactiveBoot ~ FN1UA + FN3UA + FREMUA + FlatUA + (1|Subject)')
betaAll(:,j) = mdl.Coefficients.Estimate(2:end)
SEAll(:,j)= mdl.Coefficients.SE(2:end);
mdl2 = fitglme(BootTStd,'VactiveBoot ~ FN1UA + FN3UA + FREMUA + FlatUA + (1|Subject)')
betaStdAll(:,j) = mdl.Coefficients.Estimate(2:end)./nanstd(BootT.VactiveBoot)
SEStdAll(:,j)= mdl.Coefficients.SE(2:end)./nanstd(BootT.VactiveBoot);

mdl = fitglme(BootTPSGParam,'ArIT ~ FN1 + FN3 + FREM + FLateral + (1|Subject)')
mdl.Rsquared.Ordinary
betaAll(:,j) = mdl.Coefficients.Estimate(2:end)
SEAll(:,j)= mdl.Coefficients.SE(2:end);
mdl2 = fitglme(BootTPSGParamStd,'ArIT ~ FN1 + FN3 + FREM + FLateral + (1|Subject)')
betaStdAll(:,j) = mdl2.Coefficients.Estimate(2:end)
betaStdAll(:,j) = mdl.Coefficients.Estimate(2:end)./nanstd(BootTPSGParam.ArIT)
SEStdAll(:,j)= mdl.Coefficients.SE(2:end)./nanstd(BootTPSGParam.ArIT);