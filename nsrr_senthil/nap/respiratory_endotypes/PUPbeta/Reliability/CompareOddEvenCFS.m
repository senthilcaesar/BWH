clear all; close all; clc;

addpath(genpath('C:\Users\szs88\Partners Healthcare Dropbox\SATP Group\PUPbeta_git\PUPbeta'));

%% odd and even windows
PsgOdd=load('E:\CFS\Summary\SummaryAnalysis_AllSleepAllPosComb30_Odd.mat');
PsgEven=load('E:\CFS\Summary\SummaryAnalysis_AllSleepAllPosComb30_Even.mat');

%% set thresholds
critwin=3;
traitsj=[1 2 5 6 7 8];

%% Get Number of Windows

SummaryAnalysisTableN_1=PsgOdd.SummaryAnalysisTableN;
temp = isnan(SummaryAnalysisTableN_1{:,:});
temp2 = SummaryAnalysisTableN_1{:,:};
temp2(temp)=0;
SummaryAnalysisTableN_1{:,:} = temp2;

SummaryAnalysisTableN_2=PsgEven.SummaryAnalysisTableN;
temp = isnan(SummaryAnalysisTableN_2{:,:});
temp2 = SummaryAnalysisTableN_2{:,:};
temp2(temp)=0;
SummaryAnalysisTableN_2{:,:} = temp2;

%% traits 
traitsLabelAll={'Loop Gain,LG1','Ventilatory Instability,LGn','Delay','VRA','Arousal Threshold',...
    'Collapsibility','Vactive','Compensation'};

% traitsj=[6,8,1,2,5];

traitsj=[1 2 3 4 5 6 7 8];
traitsLabelj=traitsLabelAll(traitsj);
thres=critwin*ones(1,8);

SummaryAnalysisTable_1=PsgOdd.SummaryAnalysisTable;
SummaryAnalysisTable_2=PsgEven.SummaryAnalysisTable;

for i=1:length(traitsj)
    
    %     subplot(2,length(traitsj),i)
    j=traitsj(i);
    
    N1=SummaryAnalysisTableN_1{:,j}; % comb=1=odd windows
    N2=SummaryAnalysisTableN_2{:,j}; % comb=2=even windows
    
    % use only those subjects who has windows greater than threshold
    crit = N1>=thres(j) & N2>=thres(j);
       
    
    Nsubj(:,i)=crit;
    Nactual(1,i)=nanmedian(N1(crit)); % median # of odd windows in the actual file
    Nactual(2,i)=nanmedian(N2(crit)); % median # of even windows in the actual file
    
    Trait2 = SummaryAnalysisTable_2{crit==1,j}; % trait meeting criteria--even windows
    Trait1 = SummaryAnalysisTable_1{crit==1,j}; % trait meeting criteria---odd window
    
          
    if j==5 % transform ArTh
        Trait1 = fArThresT(Trait1);
        Trait2 = fArThresT(Trait2);
    end
    
    if j==6 % transform Vpassive
        Trait1 = fVpassiveT(Trait1);
        Trait2 = fVpassiveT(Trait2);
    end
    
    if j==7 % transform Vactive
        Trait1 = fVpassiveT(Trait1);
        Trait2 = fVpassiveT(Trait2);
    end
    
    if j==8 % limit
        VcompLimit=100;
        Trait2(Trait2>VcompLimit)=VcompLimit;
        Trait1(Trait1>VcompLimit)=VcompLimit;
    end
    
    % MEAN, SEM, DELTA TRAIT & CORRELATION
    TraitMeanOdd(j)=nanmean(Trait1);
    TraitMeanOddSE(j)=[nanstd(Trait1)/sqrt(length(Trait1))];
    TraitMeanEven(j)=nanmean(Trait2);
    TraitMeanEvenSE(j)=[nanstd(Trait2)/sqrt(length(Trait2))];
    TraitDel(j)=nanmean(Trait2-Trait1);
    TraitDelSE(j)=[nanstd(Trait2-Trait1)/sqrt(length(Trait2-Trait1))];
    
    RTrait2(j) = corr(Trait1,Trait2,'Rows','complete')
end
T=table(traitsLabelAll',TraitMeanOdd(:),TraitMeanOddSE(:),TraitMeanEven(:),...
    TraitMeanEvenSE(:),TraitDel(:),TraitDelSE(:),RTrait2(:));
T.Properties.VariableNames={'Endotype','Odd_Mean','Odd_Se','Even_Mean','Even_Se','Delta_Mean','Delta_Se','Correlation'};

