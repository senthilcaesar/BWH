%% clear workspace 
clear all
close all

%% Options
prefix = '';
savename = 'SummaryAnalysisTable01102018.mat'; %change to OA3_AB from OA3_SS
fnameimport = 'MADData_NoTraits.xlsx';


%% Load matlab data 
load([prefix savename]);

%% load additional ID numbers
[~,~,SubjectID]=xlsread('AMasterSpreadSheet.xlsx',1,'R4:R1003')
SubjectID(size(SummaryAnalysisTable,1)+1:end)=[];
SummaryAnalysisTable.SubjectID = SubjectID;

%% Import clinical variables 
T = readtable(fnameimport);

%%  Merge Tables
T2 = outerjoin(T,SummaryAnalysisTable);

%% Predicting outcomes using multiple traits 

%% Setup Transformed variables
Vpassive = T2.Vpassive;
ArThres = T2.ArThres;
LGn = T2.LGn;
LG1 = T2.LG1;
VRA = T2.VRA;
Vcomp = T2.Vcomp;

%Transform ArThres:
temp=ArThres; temp(temp<100)=100;
ArThres0p5 = 100+(100*((temp-100)/100).^0.5);
%[swtest(ArThres) swtest(ArThres0p5)]

%Transform Vpassive:
temp=Vpassive; temp(Vpassive>100)=100;
Vpassive0p5 = 100-(100*((100-temp)/100).^0.5);
%[swtest(Vpassive) swtest(Vpassive0p5)]

BaselineAHI = T2.BLPSG_AHI;
TreatmentAHI = T2.M3_AHI;
DeltaAHIp = 100*(BaselineAHI-TreatmentAHI)./BaselineAHI;
DeltaAHIf = 100*(BaselineAHI-TreatmentAHI)./(BaselineAHI+TreatmentAHI);
[swtest(DeltaAHIp) swtest(DeltaAHIf)]
if 0
    figure(2)
    subplot(1,2,1); histogram(DeltaAHIp,10);
    subplot(1,2,2); histogram(DeltaAHIf,10);
end


%% Run Atemp.m

DataTable = T2(:,132:139);
load Model
[Pred,Y] = EndotypePredictOAEfficacy(DataTable,Model);

%%
converttop = @(x) 2*x./(x+1);

cutoff=0.4285;
figure(1); clf(1)
subplot(2,1,1)
histogram(DeltaAHIf(Y>cutoff),-100:10:100)
hold('on');
subplot(2,1,2)
histogram(DeltaAHIf(Y<=cutoff),-100:10:100)

converttop([mean(DeltaAHIf(Y>cutoff)) mean(DeltaAHIf(Y<=cutoff))]/100)*100

[~,p]=ttest2(DeltaAHIf(Y>cutoff),DeltaAHIf(Y<=cutoff));









    
    