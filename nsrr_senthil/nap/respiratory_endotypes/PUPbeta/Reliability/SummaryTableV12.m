
OSAAHIthres=5;
OSA = Ahi3pArV1>OSAAHIthres & Ahi3pArV2>OSAAHIthres;
sum(Ahi3pArV1>OSAAHIthres)
sum(Ahi3pArV2>OSAAHIthres)
sum(OSA)

MeanV1=[nanmean(AgeV1(OSA)) nanmean(BMIV1(OSA)) nanmean(Ahi3pArV1(OSA))];
MeanV2=[nanmean(AgeV2(OSA)) nanmean(BMIV2(OSA)) nanmean(Ahi3pArV2(OSA))];

deltaAge=nanmean(AgeV2(OSA)-AgeV1(OSA));
deltaBMI=nanmean(BMIV2(OSA)-BMIV1(OSA));
deltaAHI=nanmean(Ahi3pArV2(OSA)-Ahi3pArV1(OSA));
RTrait2(12)=corr(Ahi3pArV2(OSA),Ahi3pArV1(OSA),'row','Complete');
RTrait2(11)=corr(BMIV2(OSA),BMIV1(OSA),'row','Complete');
RTrait2(10)=corr(AgeV2(OSA),AgeV1(OSA),'row','Complete');

figure(88)
subplot(2,4,7)
scatter(Ahi3pArV1(OSA),Ahi3pArV2(OSA),10,'filled','markerfacealpha',0.5);
xlabel('AHI_ Visit1'); ylabel('AHI_ Visit2');

subplot(2,4,8)
scatter(BMIV1(OSA),BMIV2(OSA),10,'filled','markerfacealpha',0.5);
xlabel('BMI_ Visit1'); ylabel('BMI_ Visit2');


%%
critwin=3;

%% for sas gee analysis
if 0
temp=table2array(SummaryAnalysisTableN_1);
crit1=temp>critwin;
temp=table2array(SummaryAnalysisTableN_2);
crit2=temp>critwin;
CritBoth=crit1&crit2;
CritBothOSA=CritBoth&OSA;
V1Export=table2array(SummaryAnalysisTable1(:,1:8));
V1Export(~CritBothOSA)=NaN;
V2Export=table2array(SummaryAnalysisTable2(:,1:8));
V2Export(~CritBothOSA)=NaN;
for j=1:8
    V1FN1(:,j) = FstatesArray_1{j}(:,1);
    V2FN1(:,j) = FstatesArray_2{j}(:,1);
    V1FN2(:,j) = FstatesArray_1{j}(:,2);
    V2FN2(:,j) = FstatesArray_2{j}(:,2);
    V1FN3(:,j) = FstatesArray_1{j}(:,3);
    V2FN3(:,j) = FstatesArray_2{j}(:,3);
    V1REM(:,j) = FstatesArray_1{j}(:,4);
    V2REM(:,j) = FstatesArray_2{j}(:,4);
    V1Fsup(:,j) = Fsupine_1(:,j); % change in supine position
    V2Fsup(:,j)= Fsupine_2(:,j);
end
V1FN1(~CritBothOSA)=NaN;
V2FN1(~CritBothOSA)=NaN;
V1FN2(~CritBothOSA)=NaN;
V2FN2(~CritBothOSA)=NaN;
V1FN3(~CritBothOSA)=NaN;
V2FN3(~CritBothOSA)=NaN;
V1REM(~CritBothOSA)=NaN;
V2REM(~CritBothOSA)=NaN;
V1Fsup(~CritBothOSA)=NaN;
V2Fsup(~CritBothOSA)=NaN;

BMIV1(~OSA)=NaN;
BMIV2(~OSA)=NaN;
AgeV1(~OSA)=NaN;
AgeV2(~OSA)=NaN;
Ahi3pArV1(~OSA)=NaN;
Ahi3pArV2(~OSA)=NaN;

SubjExp=(1:1026)';
TreatExp1=repmat([1],1026,1);
TreatExp2=repmat([2],1026,1);

VarNames={'Subj','Visit','Trait','FN1','FN2','FN3','REM','Fsup','Age','BMI','AHI'};
outfilen=['C:\Users\rma56\Dropbox (Partners HealthCare)\Reliability Trait Analysis\MrOS\GEE Analysis\TraitsSummary.xlsx'];

LG1Exp=[SubjExp,TreatExp1,V1Export(:,1),V1FN1(:,1),V1FN2(:,1),V1FN3(:,1),V1REM(:,1),V1Fsup(:,1),AgeV1,BMIV1,Ahi3pArV1;...
   SubjExp,TreatExp2,V2Export(:,1),V2FN1(:,1),V2FN2(:,1),V2FN3(:,1),V2REM(:,1),V2Fsup(:,1),AgeV2,BMIV2,Ahi3pArV2];
LG1Exp=sortrows(LG1Exp,1);
LG1Exp(:,isnan(LG1Exp))=[];
LG1Exp=array2table(LG1Exp);
LG1Exp.Properties.VariableNames=VarNames;
writetable(LG1Exp,outfilen,'Sheet','LG1');

LGnExp=[SubjExp,TreatExp1,V1Export(:,2),V1FN1(:,2),V1FN2(:,2),V1FN3(:,2),V1REM(:,2),V1Fsup(:,2),AgeV1,BMIV1,Ahi3pArV1;...
   SubjExp,TreatExp2,V2Export(:,2),V2FN1(:,2),V2FN2(:,2),V2FN3(:,2),V2REM(:,2),V2Fsup(:,2),AgeV2,BMIV2,Ahi3pArV2];
LGnExp=sortrows(LGnExp,1);
LGnExp=array2table(LGnExp);
LGnExp.Properties.VariableNames=VarNames;
writetable(LGnExp,outfilen,'Sheet','LGn');

DelayExp=[SubjExp,TreatExp1,V1Export(:,3),V1FN1(:,3),V1FN2(:,3),V1FN3(:,3),V1REM(:,3),V1Fsup(:,3),AgeV1,BMIV1,Ahi3pArV1;...
    SubjExp,TreatExp2,V2Export(:,3),V2FN1(:,3),V2FN2(:,3),V2FN3(:,3),V2REM(:,3),V2Fsup(:,3),AgeV2,BMIV2,Ahi3pArV2];
DelayExp=sortrows(DelayExp,1);
DelayExp=array2table(DelayExp);
DelayExp.Properties.VariableNames=VarNames;
writetable(DelayExp,outfilen,'Sheet','Delay');

VRAExp=[SubjExp,TreatExp1,V1Export(:,4),V1FN1(:,4),V1FN2(:,4),V1FN3(:,4),V1REM(:,4),V1Fsup(:,4),AgeV1,BMIV1,Ahi3pArV1;...
    SubjExp,TreatExp2,V2Export(:,4),V2FN1(:,4),V2FN2(:,4),V2FN3(:,4),V2REM(:,4),V2Fsup(:,4),AgeV2,BMIV2,Ahi3pArV2];
VRAExp=sortrows(VRAExp,1);
VRAExp=array2table(VRAExp);
VRAExp.Properties.VariableNames=VarNames;
writetable(VRAExp,outfilen,'Sheet','VRA');

ArThExp=[SubjExp,TreatExp1,V1Export(:,5),V1FN1(:,5),V1FN2(:,5),V1FN3(:,5),V1REM(:,5),V1Fsup(:,5),AgeV1,BMIV1,Ahi3pArV1;...
    SubjExp,TreatExp2,V2Export(:,5),V2FN1(:,5),V2FN2(:,5),V2FN3(:,5),V2REM(:,5),V2Fsup(:,5),AgeV2,BMIV2,Ahi3pArV2];
ArThExp=sortrows(ArThExp,1);
ArThExp=array2table(ArThExp);
ArThExp.Properties.VariableNames=VarNames;
writetable(ArThExp,outfilen,'Sheet','ArTh');

VpassiveExp=[SubjExp,TreatExp1,V1Export(:,6),V1FN1(:,6),V1FN2(:,6),V1FN3(:,6),V1REM(:,6),V1Fsup(:,6),AgeV1,BMIV1,Ahi3pArV1;...
    SubjExp,TreatExp2,V2Export(:,6),V2FN1(:,6),V2FN2(:,6),V2FN3(:,6),V2REM(:,6),V2Fsup(:,6),AgeV2,BMIV2,Ahi3pArV2];
VpassiveExp=sortrows(VpassiveExp,1);
VpassiveExp=array2table(VpassiveExp);
VpassiveExp.Properties.VariableNames=VarNames;
writetable(VpassiveExp,outfilen,'Sheet','Vpassive');

VactiveExp=[SubjExp,TreatExp1,V1Export(:,7),V1FN1(:,7),V1FN2(:,7),V1FN3(:,7),V1REM(:,7),V1Fsup(:,7),AgeV1,BMIV1,Ahi3pArV1;...
    SubjExp,TreatExp2,V2Export(:,7),V2FN1(:,7),V2FN2(:,7),V2FN3(:,7),V2REM(:,7),V2Fsup(:,7),AgeV2,BMIV2,Ahi3pArV2];
VactiveExp=sortrows(VactiveExp,1);
VactiveExp=array2table(VactiveExp);
VactiveExp.Properties.VariableNames=VarNames;
writetable(VactiveExp,outfilen,'Sheet','Vactive');

VcompExp=[SubjExp,TreatExp1,V1Export(:,8),V1FN1(:,8),V1FN2(:,8),V1FN3(:,8),V1REM(:,8),V1Fsup(:,8),AgeV1,BMIV1,Ahi3pArV1;...
    SubjExp,TreatExp2,V2Export(:,8),V2FN1(:,8),V2FN2(:,8),V2FN3(:,8),V2REM(:,8),V2Fsup(:,8),AgeV2,BMIV2,Ahi3pArV2];
VcompExp=sortrows(VcompExp,1);
VcompExp=array2table(VcompExp);
VcompExp.Properties.VariableNames=VarNames;
writetable(VcompExp,outfilen,'Sheet','Vcomp');
end


    


%%

traitsj=SummaryAnalysisTable1.Properties.VariableNames;
if  size(traitsj,2)>8
    traitsj(9:end)=[];
end

thres=3*ones(length(traitsj),1);

temp = isnan(SummaryAnalysisTableN_1{:,:});
temp2 = SummaryAnalysisTableN_1{:,:};
temp2(temp)=0;
SummaryAnalysisTableN_1{:,:} = temp2;

temp = isnan(SummaryAnalysisTableN_2{:,:});
temp2 = SummaryAnalysisTableN_2{:,:};
temp2(temp)=0;
SummaryAnalysisTableN_2{:,:} = temp2;

for j=1:length(traitsj)
    N1=SummaryAnalysisTableN_1{:,j}; % comb=1=odd windows
    N2=SummaryAnalysisTableN_2{:,j}; % comb=2=even windows
    crit = N1>=thres(j) & N2>=thres(j)&OSA; % use only those subjects who has windows greater than threshold
    
    Nactual(1,j)=nanmedian(N1(crit)); % median # of odd windows in the actual file
    Nactual(2,j)=nanmedian(N2(crit)); % median # of even windows in the actual file
    
    deltaFN1 = FstatesArray_2{j}(crit,1) - FstatesArray_1{j}(crit,1); % change in n1
    deltaFN3 = FstatesArray_2{j}(crit,3) - FstatesArray_1{j}(crit,3); % change in n3
    deltaFW = (1-sum(FstatesArray_2{j}(crit,:),2)) - (1-sum(FstatesArray_1{j}(crit,:),2)); % change in wake
    deltaFsup = Fsupine_2(crit,j) - Fsupine_1(crit,j); % change in supine position
    Trait2 = SummaryAnalysisTable2{crit==1,j}; % trait from even window
    Trait1 = SummaryAnalysisTable1{crit==1,j}; % trait from odd window
    NSub(j)=size(Trait1,1);
    TraitMeanV1(j)=nanmean(Trait1);
    TraitMeanV2(j)=nanmean(Trait2);
    TraitDel(j)=nanmean(Trait2-Trait1);
    RTrait2(j) = corr(Trait1,Trait2,'Rows','complete')
    
figure(88)
if j==1
subplot(2,4,1)
scatter(Trait1,Trait2,10,'filled','markerfacealpha',0.5);
xlabel('LG1(Visit1)'); ylabel('LG1(Visit2)');
elseif j==2
subplot(2,4,2)
scatter(Trait1,Trait2,10,'filled','markerfacealpha',0.5);
xlabel('LGn(Visit1)'); ylabel('LGn(Visit2)');
elseif j==5
subplot(2,4,3)
scatter(Trait1,Trait2,10,'filled','markerfacealpha',0.5);
xlabel('ArTh(Visit1)'); ylabel('ArTh(Visit2)');
elseif j==6
subplot(2,4,4)
scatter(Trait1,Trait2,10,'filled','markerfacealpha',0.5);
xlabel('Vpassive(Visit1)'); ylabel('Vpassive(Visit2)');
elseif j==7
subplot(2,4,5)
scatter(Trait1,Trait2,10,'filled','markerfacealpha',0.5);
xlabel('Vactive(Visit1)'); ylabel('Vactive(Visit2)');
elseif j==8
subplot(2,4,6)
scatter(Trait1,Trait2,10,'filled','markerfacealpha',0.5);
xlabel('Vcomp(Visit1)'); ylabel('Vcomp(Visit2)');
end
end