clear all;clc;

%% Options
OSAAHIthres=5;% poordi4 
critwin=3;

%% AHI in MrOS

Nsrr1=readtable('C:\Users\rma56\Dropbox (Partners HealthCare)\Reliability Trait Analysis\MrOS\NSRR Updated AHI\mros1_nsrr_new_ahis_20200518.csv');
Nsrr2=readtable('C:\Users\rma56\Dropbox (Partners HealthCare)\Reliability Trait Analysis\MrOS\NSRR Updated AHI\mros2_nsrr_new_ahis_20200518.csv');

for ii=1:size(Nsrr1,1)
   indx=find(strcmp(Nsrr2.NSRRID,Nsrr1.NSRRID(ii)));
     if ~isempty(indx)
        Nsrr1.Flag(ii)=1;
    else
        Nsrr1.Flag(ii)=0;
    end
end
Ahi3pArV1=Nsrr1.POOAHI3(Nsrr1.Flag==1);
Ahi3pArV2=Nsrr2.POOAHI3;

OSAV1=Ahi3pArV1>OSAAHIthres;
OSAV2=Ahi3pArV2>OSAAHIthres;
OSA=Ahi3pArV1>OSAAHIthres & Ahi3pArV2>OSAAHIthres;

%% get the parameters from MESA BootTAnalysis
% using Mesa for getting beta coeffs
TrainBootT=load('C:\Users\rma56\Dropbox (Partners HealthCare)\Reliability Trait Analysis\MESA\SummaryAnalysis_AllSleep_AllPos_Boot.mat');
% TrainBootT=load('C:\Users\rma56\Dropbox (Partners HealthCare)\Reliability Trait Analysis\MrOS\Summary\SummaryAnalysis_V1_AllSleep_AllPos_Boot.mat');

BootT=TrainBootT.BootT;

%code to limit data based on min number of windows
DataNtemp=BootT{:,10:16}<critwin;
DataNtemp = [DataNtemp(:,1:4) DataNtemp(:,4) DataNtemp(:,5:7)];
temp = BootT{:,2:9};
temp(DataNtemp)=NaN;
BootT{:,2:9} = temp;

%code to remove those with low AHI here
I = find(~OSAV1)';
NotOSAline = sum(BootT.Subject==I,2)>0;
temp = BootT{:,2:9};
temp(NotOSAline,:)=NaN;
BootT{:,2:9} = temp;

varlist=TrainBootT.varlist;
BootT.FlatLoopGain = 1 - BootT.FsupLoopGain;
BootT.FlatArTh = 1 - BootT.FsupArTh;
BootT.FlatUA = 1 - BootT.FsupUA;
j=1
TrainMdl1 = fitglme(BootT,'LG1Boot ~ FN1LoopGain + FN3LoopGain + FREMLoopGain + FlatLoopGain + (1|Subject)');
TrainBetaAll(:,j) = TrainMdl1.Coefficients.Estimate(2:end);
j=2
TrainMdl2 = fitglme(BootT,'LGnBoot ~ FN1LoopGain + FN3LoopGain + FREMLoopGain + FlatLoopGain + (1|Subject)');
TrainBetaAll(:,j) = TrainMdl2.Coefficients.Estimate(2:end);
j=3
TrainMdl3 = fitglme(BootT,'delayBoot ~ FN1LoopGain + FN3LoopGain + FREMLoopGain + FlatLoopGain + (1|Subject)');
TrainBetaAll(:,j) = TrainMdl3.Coefficients.Estimate(2:end);
j=4
TrainMdl4 = fitglme(BootT,'VRABoot ~ FN1ArTh + FN3ArTh + FREMArTh + FlatArTh + (1|Subject)');
TrainBetaAll(:,j) = TrainMdl4.Coefficients.Estimate(2:end);

j=5
temp = (BootT.ArThresBoot/100)-1;
temp(temp<0)=0;
BootT.ArThresBootT=100*((temp.^0.5)+1);
%mdl5_ = fitglme(BootT,'ArThresBoot ~ FN1ArTh + FN3ArTh + FREMArTh + FlatArTh + (1|Subject)')
TrainMdl5 = fitglme(BootT,'ArThresBootT ~ FN1ArTh + FN3ArTh + FREMArTh + FlatArTh + (1|Subject)');
TrainBetaAll(:,j) = TrainMdl5.Coefficients.Estimate(2:end);

j=6
temp = 1-(BootT.VpassiveBoot/100);
temp(temp<0)=0;
BootT.VpassiveBootT=100*(1-(temp.^0.5));
% mdl6_ = fitglme(BootT,'VpassiveBoot ~ FN1UA + FN3UA + FREMUA + FlatUA + (1|Subject)')
TrainMdl6 = fitglme(BootT,'VpassiveBootT ~ FN1UA + FN3UA + FREMUA + FlatUA + (1|Subject)');
TrainBetaAll(:,j) = TrainMdl6.Coefficients.Estimate(2:end);

j=7
TrainMdl7 = fitglme(BootT,'VactiveBoot ~ FN1UA + FN3UA + FREMUA + FlatUA + (1|Subject)');
TrainBetaAll(:,j) = TrainMdl7.Coefficients.Estimate(2:end);

j=8
TrainMdl8 = fitglme(BootT,'VcompBoot ~ FN1UA + FN3UA + FREMUA + FlatUA + (1|Subject)');
TrainBetaAll(:,j) = TrainMdl8.Coefficients.Estimate(2:end);

betaNames = {'FN1','FN3','FREM','Flateral'};
TrainBetaAll = array2table(TrainBetaAll);
TrainBetaAll.Properties.RowNames = betaNames;
TrainBetaAll.Properties.VariableNames = varlist;
TrainBetaAll.Properties.VariableNames{5}='ArThresT';
TrainBetaAll.Properties.VariableNames{6}='VpassiveT';
TrainBetaAll


%% Get %state and position from MrOS1 data
MrOs1BootT=load('C:\Users\rma56\Dropbox (Partners HealthCare)\Reliability Trait Analysis\MrOS\Summary\SummaryAnalysis_V1_AllSleep_AllPos_Boot.mat');
DataTV1=array2table(MrOs1BootT.DataArray);

%code to limit data based on min number of windows
DataNtemp=MrOs1BootT.DataNArray<critwin;
temp = DataTV1{:,:};
temp(DataNtemp)=NaN;
DataTV1{:,:} = temp;

%code to remove those with low AHI here
temp = DataTV1{:,:};
temp(~OSAV1,:)=NaN;
DataTV1{:,:} = temp;

DataTV1.Properties.VariableNames = varlist;
clear temp
temp = (DataTV1.ArThres/100)-1;
temp(temp<0)=0;
DataTV1.ArThres=100*((temp.^0.5)+1);
DataTV1.Properties.VariableNames{5}='ArThresT';

temp = 1-(DataTV1.Vpassive/100);
temp(temp<0)=0;
DataTV1.Vpassive=100*(1-(temp.^0.5));
DataTV1.Properties.VariableNames{6}='VpassiveT';

FstatesArray1=MrOs1BootT.FstatesArray;
FsupineArray1=MrOs1BootT.FsupineArray;


%% Apply to the MrOS2 data
MrOs2BootT=load('C:\Users\rma56\Dropbox (Partners HealthCare)\Reliability Trait Analysis\MrOS\Summary\SummaryAnalysis_V2_AllSleep_AllPosUpdated_Boot.mat');
DataArray=MrOs2BootT.DataArray;
DataNtemp=MrOs2BootT.DataNArray<critwin;

if size(DataArray,1)>1026
    DataArray(1:2906,:)=[];
    DataNtemp(1:2906,:)=[];
end
DataT = array2table(DataArray);
%code to limit data based on min number of windows

temp = DataT{:,:};
temp(DataNtemp)=NaN;
DataT{:,:} = temp;

%code to remove those with low AHI here
temp = DataT{:,:};
temp(~OSAV2,:)=NaN;
DataT{:,:} = temp;

varlist=MrOs2BootT.varlist;
FstatesArray=MrOs2BootT.FstatesArray;
if size(FstatesArray{1,1},1)>1026
for i=1:size(FstatesArray,2)
    FstatesArray{i}(1:2906,:)=[];
end
end
FsupineArray=MrOs2BootT.FsupineArray;
if size(FsupineArray,1)>1026
    FsupineArray(1:2906,:)=[];
end

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

optiontest = 2;
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
            beta = TrainBetaAll{:,j};
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


%% correlation between v2 and v1
for ii=1:length(varlist)
    
    %original correlation
    [RHO(ii),PVAL(ii)] = corr(DataTV1.(ii),DataTT.(ii),'Type','Pearson','rows','complete');
    
    [B]=fitglm(DataTV1.(ii),DataTT.(ii));
        B.Rsquared.Ordinary.^0.5;
        
    % correlation after correcting v2
    [RHOcc(ii),PVALcc(ii)] = corr(DataTV1.(ii),DataTTc.(ii),'Type','Pearson','rows','complete');
end
RHO
RHOcc