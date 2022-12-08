clear all; close all; clc;

whichpc=1;

switch whichpc
    case 1
        userstr = 'C:\Users\rma56\Dropbox (Partners HealthCare)\';
    case 2
        userstr = 'G:\Partners Healthcare Dropbox\SATP Group\';
    case 3
        userstr='C:\Users\szs88\Partners Healthcare Dropbox\SATP Group\';
end

addpath(genpath([userstr 'PUPbeta_git\']));
addpath(genpath([userstr 'Reliability Trait Analysis\MrOS\']));
%addpath(genpath([userstr 'Reliability Trait Analysis\MrOS\Summary']));
%addpath(genpath([userstr 'Reliability Trait Analysis\MrOS\NSRR Table']));

%% find the average years between v1 and v2
Demo1=load([userstr 'Reliability Trait Analysis\MrOS\NSRR Table\mros_visit1-2.mat']);
date1=datetime(Demo1.T1_new.postdydt);
date2=datetime(Demo1.T2_new.postdydt);
diffdate=years(date2-date1);
dtmean=mean(diffdate)
median(diffdate)
rangeyrs=[min(diffdate), max(diffdate),range(diffdate)]
sdyears=std(diffdate)
% cr for subjects lesser and greater than dtmean

%% ALL SLEEP TABLE
file1='SummaryAnalysis_V1_AllSleep_AllPos_Boot.mat'; % boot table has summary as well
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
clear SummaryAnalysisTable1_ SummaryAnalysisTable2_


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

%% tables for demographics etc
[~,~,StudyID] = xlsread('AMasterSpreadsheet_1026visit2.xlsx',1,'T4:T1029');
% from master spreadsheet
nsrrid = {};
for i=1:length(StudyID)
    nsrrid{i,1} = upper(extractAfter(StudyID{i},'visit2-'));
end

DataTableA = table(nsrrid);

%% nsrr table for visit 1 and 2; get AHI from here

% ahi from new updated nsrr tables
Nsrr1=readtable([userstr 'Reliability Trait Analysis\MrOS\NSRR Updated AHI\mros1_nsrr_new_ahis_20200601.csv']);
Nsrr2=readtable([userstr 'Reliability Trait Analysis\MrOS\NSRR Updated AHI\mros2_nsrr_new_ahis_20200601.csv']);
for ii=1:size(Nsrr1,1)
    indx=find(strcmp(Nsrr2.NSRRID,Nsrr1.NSRRID(ii)));
    if ~isempty(indx)
        Nsrr1.Flag(ii)=1;
    else
        Nsrr1.Flag(ii)=0;
    end
end
% All sleep
Ahi3pArV1=Nsrr1.POAHI3A(Nsrr1.Flag==1);
Ahi3pArV2=Nsrr2.POAHI3A;

%NREM
NremAhi3pArV1=str2double((Nsrr1.POAHI3AN(Nsrr1.Flag==1)));
NremAhi3pArV2=str2double(Nsrr2.POAHI3AN);

% REM
RemAhi3pArV1=str2double(Nsrr1.POAHI3AR(Nsrr1.Flag==1));
RemAhi3pArV2=str2double(Nsrr2.POAHI3AR);


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

%% criteria

%%%RUN this to get plots with orange data
thres=3*ones(length(traitsj),1); % minimum 3 windows needed to be considered

%%%RUN this to get actual results for MESA min windows
% from MESA within-night calculations of minimum windows
thres=[34 30 30 20 20 28 28 35]'; % dummy values for delay,vra & vactive

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

%%
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

% %BMI
%     BMIV1=str2double(T1_new.hwbmi);
%     BMIV2=str2double(T2_new.hwbmi);
%     corr(BMIV1,BMIV2,'row','Complete')
%     DeltaBMI =  BMIV2-BMIV1;

%% TABLE 1 (skip)
if 0
    
    % for table1 - use only baseline data all subjects
    AHITable1=table(Nsrr1.NSRRID,Nsrr1.POAHI3A,Nsrr1.POAHI3AN,Nsrr1.POAHI3AR);
    AHITable1.Properties.VariableNames={'nsrrid','Ahi3pAr','NremAhi3pAr','RemAhi3pAr'};
    % age, bmi
    Demo1=readtable('mros-visit1-dataset-0.3.0.csv');
    
    % sleep metrics
    IdLink=readtable('mros_id_linking.csv');
    
    %PSG1 has different id- need to link back to nsrr id using IdLink excel
    %sheet
    PSG1=readtable('mros1psg_20121019.csv');
    I = find(PSG1.Properties.VariableNames=="PPTID")
    PSG1.Properties.VariableNames{I}='id';
    PSG1T=outerjoin(PSG1,IdLink,'MergeKeys', true);
    PSG1T = movevars(PSG1T,'nsrrid','After','id');
    
    nssridT=table(Demo1.nsrrid);
    nssridT.Properties.VariableNames={'nsrrid'};
    
    PSG1FinalT=outerjoin(nssridT,PSG1T);
    
    I = find(PSG1FinalT.Properties.VariableNames=="nsrrid_nssridT")
    PSG1FinalT.Properties.VariableNames{I}='nsrrid';
    PSG1FinalT(strcmp(PSG1FinalT.nsrrid,""),:)=[];
    PSG1FinalT.nsrrid_PSG1T=[];
    clear PSG1 PSG1T IdLink nssridT
    
     
    Visit1Table=outerjoin(AHITable1,Demo1,'MergeKeys', true);
    Visit1Table=outerjoin(Visit1Table,PSG1FinalT,'MergeKeys', true);
    clear AHITable1 Demo1 PSG1FinalT
    
    mrosendotypesV1=load('SummaryAnalysis_V1_AllSleep_AllPos_BootT_Oct092021.mat');
    % get ids from excel sheet
    [~,~,IdV1] = xlsread('AMasterSpreadsheet.xlsx',1,'T4:T2909');
   
    nsrridV1 = {};
    for i=1:length(IdV1)
        nsrridV1{i,1} = upper(extractAfter(IdV1{i},'visit1-'));
    end
    AmasterspreadsheetV1 = table(nsrridV1);
    
    V1endotypeNT=table();
    V1endotypeNT=AmasterspreadsheetV1;
    V1endotypeNT=[AmasterspreadsheetV1,mrosendotypesV1.SummaryAnalysisTableN];
    V1endotypeNT.Properties.VariableNames(1)={'nsrrid'};
    
    % join v1 psg parameters and endotype win numbers
    PSG1FinalTEndo=outerjoin(Visit1Table,V1endotypeNT);
    loc=cellfun('isempty', PSG1FinalTEndo{:,'nsrrid_V1endotypeNT'});
    PSG1FinalTEndo(loc,:)=[];
    I = find(PSG1FinalTEndo.Properties.VariableNames=="nsrrid_Visit1Table")
    PSG1FinalTEndo.Properties.VariableNames{I}='nsrrid';
  
    % calculate patient flow
    
    failedstudy = isnan(PSG1FinalTEndo.Ahi3pAr);
    sum(failedstudy)
    I2= PSG1FinalTEndo.Ahi3pAr<5; % not osa at v1
    sum(I2)
    
    temp=PSG1FinalTEndo.Vpassive;
    temp(temp<3)=NaN;
    temp(failedstudy,:)=NaN;
    I1=~isnan(temp);
    sum(I1)
    I = I1 & ~I2;
    sum(I)
    
    Table1List = {'vsage1','gender','hwbmi',...
        'Ahi3pAr','NremAhi3pAr','RemAhi3pAr',...
        'ai_all',...
        'SLPPRDP','TMSTG1P','TMSTG2P','TMSTG34P','TMREMP',...
        'supinep'};
    %'eventdur','arousaldur','ODI3','ODI4',
    
    TableListMean={'vsage1','hwbmi',...
        'Ahi3pAr','NremAhi3pAr','RemAhi3pAr','ai_all',...
        'SLPPRDP','TMSTG1P','TMSTG2P','TMSTG34P','TMREMP','supinep'};
    %         'eventdur','arousaldur','ODI3','ODI4',...
    
    
    Table1Pt = PSG1FinalTEndo(I,Table1List);
    
   if 1 % find V2 subset patient flow
       AHITable2=table(Nsrr2.NSRRID,Nsrr2.POAHI3A,Nsrr2.POAHI3AN,Nsrr2.POAHI3AR);
       AHITable2.Properties.VariableNames={'nsrrid','Ahi3pAr_2','NremAhi3pAr_2','RemAhi3pAr_2'};
    
       Demo2=readtable('mros-visit2-dataset-0.3.0.csv');
       Demo2T=table(Demo2.nsrrid,Demo2.vs2age1,Demo2.gender,Demo2.hwbmi); % only age, bmi and gender
       Demo2T.Properties.VariableNames={'nsrrid','vs2age1','gender_2','hwbmi_2'};
       Visit2Table=outerjoin(AHITable2,Demo2T,'MergeKeys', true);
      
       
       V2endotypeNT=table();
       AmasterspreadsheetV2=table(nsrrid);
       V2endotypeNT=[AmasterspreadsheetV2,SummaryAnalysisTableN_2];
       V2endotypeNT.Properties.VariableNames={'nsrrid','LG1_2','LGn_2','delay_2','VRA_2','ArThres_2','Vpassive_2','Vactive_2','Vcomp_2'};
         
       PSG2TEndo=outerjoin(Visit2Table,V2endotypeNT,'MergeKeys', true);
       
       % calculate patient flow
       
       failedstudy = isnan(PSG2TEndo.Ahi3pAr_2);
       sum(failedstudy)
       I2=PSG2TEndo.Ahi3pAr_2<5; % not osa at v1
       sum(I2)
       
       temp=PSG2TEndo.Vpassive_2;
       temp(temp<3)=NaN;
       temp(temp<thres(6))=NaN; % mesa criteria
       temp(failedstudy,:)=NaN;
       I1=~isnan(temp);
       sum(I1)
       I = I1 & ~I2;
       sum(I)
       
%        crit = N1>=thres(6) & N2>=thres(6)
       
       clear PSG12Endo I2 I I1
       PSG12Endo=outerjoin(PSG1FinalTEndo,PSG2TEndo); % both 1 and 2
       
       loc=cellfun('isempty', PSG12Endo{:,'nsrrid_PSG2TEndo'});
       PSG12Endo(loc,:)=[];
       
       failedstudy = isnan(PSG12Endo.Ahi3pAr_2)| isnan(PSG12Endo.Ahi3pAr) ;
       sum(failedstudy)
       I2= PSG12Endo.Ahi3pAr<5 | PSG12Endo.Ahi3pAr_2<5; % not osa at v1
       sum(I2)
       
       temp=PSG12Endo.Vpassive;
       temp(temp<thres(6))=NaN;
       temp(failedstudy,:)=NaN;
       I1=~isnan(temp);
       temp=PSG12Endo.Vpassive_2;
       temp(temp<thres(6))=NaN;
       I3=~isnan(temp);
       sum(I1)
       I = I1 & I3& ~I2;
       sum(I)
       
       clear Table1Pt Table1List TableListMean
        
       Table1List = {'vsage1','gender','hwbmi','vs2age1','gender_2','hwbmi_2'...
           'Ahi3pAr','NremAhi3pAr','RemAhi3pAr','Ahi3pAr_2',...
           'ai_all',...
           'SLPPRDP','TMSTG1P','TMSTG2P','TMSTG34P','TMREMP',...
           'supinep'};
       %'eventdur','arousaldur','ODI3','ODI4',
       
    
       TableListMean={'vsage1','hwbmi','vs2age1','hwbmi_2'...
           'Ahi3pAr','NremAhi3pAr','RemAhi3pAr','ai_all','Ahi3pAr_2',...
           'SLPPRDP','TMSTG1P','TMSTG2P','TMSTG34P','TMREMP','supinep'};
       %         'eventdur','arousaldur','ODI3','ODI4',...

       Table1Pt = PSG12Endo(I,Table1List);
   end
    
 
    
   
    clear Visit1Table
    
    Table1Pt.male=Table1Pt.gender==2;%female=0; male=2
    Table1Pt.female=Table1Pt.gender==0;
    
    
    % for all vars with 0 and 1
    TableListN={'male','female'};
    
    Table1num=[];
    for ii=1:length(TableListN)
        temp=table2array(Table1Pt(:,TableListN{1,ii}));
        Table1num(1,ii)= sum(temp==1);
        Table1num(2,ii)=size(Table1Pt,1)-Table1num(1,ii);
    end
    OSAseverityN=[];
    temp=Table1Pt.Ahi3pAr;
    OSAseverityN(1,:)= sum([temp<5 temp>=5&temp<15 temp>=15&temp<30 temp>=30]); % osa severity
    OSAseverityN(2,:)=size(Table1Pt,1)-OSAseverityN;
    Table1num=array2table(Table1num,'VariableNames',TableListN);
    Table1num=[Table1num array2table(OSAseverityN,'VariableNames',{'normal','mild','moderate','severe'})];
    Table1num.Properties.RowNames={'n','Total#-n'};
    
    % generate mean and 95CI for tranformed/nontransformed data
    
    
    tempfactorx=[1 1 ...
        1 1 1 1 ...
        1 1 1 1 1 1];
    %         1 1 0.25 0.25 ...
    
    logconv=[0 0 ...
        0 0 0 0 ...
        0 0 0 0 0 0];
    %         0 0 0 0 ...
    
    
    nonparametric = [1 0 ...
        1 1 1 1 ...
        0 0 0 0 0 1];
    %         0 0 0 0  ...
    
     
    
    
    Table1mean=[];
    if 0
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
        pause
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



%% Get Prediction Intervals & Beta values for Trait2~ Trait1
clear T PredIn N
for i=1:length(traitsj)
    j=traitsj(i);
    N1=SummaryAnalysisTableN_1{:,j}; % visit1
    N2=SummaryAnalysisTableN_2{:,j}; % visit2
    crit = N1>=thres(i) & N2>=thres(i) & OSA; % use only those subjects who has windows greater than threshold
%     crit = N1>=3 & N2>=3 & OSA;
    
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
if 0 % double check num differences
failedstudy = isnan(PSG12Endo.Ahi3pAr_2)| isnan(PSG12Endo.Ahi3pAr) ;
       sum(failedstudy)
       I2= PSG12Endo.Ahi3pAr<5 | PSG12Endo.Ahi3pAr_2<5; % not osa at v1
       sum(~I2)
       
       temp=PSG12Endo.Vpassive;
       temp(temp<thres(6))=NaN;
       temp(failedstudy,:)=NaN;
       I1=~isnan(temp);
       temp=PSG12Endo.Vpassive_2;
       temp(temp<thres(6))=NaN;
       I3=~isnan(temp);
       sum(I1 & I3)
       I = I1 & I3& ~I2;
       sum(I)
       
end   
       
%% MAIN
traitsj=[6 8 1 2 5];
for i=1:length(traitsj)
   
    j=traitsj(i);
    N1=SummaryAnalysisTableN_1{:,j}; % visit1
    N2=SummaryAnalysisTableN_2{:,j}; % visit2
    crit = N1>=thres(j) & N2>=thres(j) & OSA; % use only those subjects who has windows greater than threshold
    
    Nactual(1,i)=nanmedian(N1(crit)); % median # of windows in the actual file
    Nactual(2,i)=nanmedian(N2(crit));
    Ntraits(1,i)=nansum(crit);
    
    deltaFlateral = Flateral2(crit,j) - Flateral1(crit,j);
    deltaTime=diffdate(crit);
    
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
    
    %    deltaBMI = DeltaBMI(crit);
    
    Trait2 = SummaryAnalysisTable2{crit==1,j}; % trait from visit2
    Trait1 = SummaryAnalysisTable1{crit==1,j}; % trait from visit1
    
    
    MeanTrait1(i,1)=nanmean(Trait1);
    SETrait1(i,1)=[nanstd(Trait1)/sqrt(length(Trait1))];
    SDTrait1(i,1)=[nanstd(Trait1)];
    MeanTrait2(i,1)=nanmean(Trait2);
    SETrait2(i,1)=[nanstd(Trait2)/sqrt(length(Trait2))];
    SDTrait2(i,1)=[nanstd(Trait2)];
    Tx = table();
    Tx.Y = [Trait2-Trait1];
    mdlDel=fitglme(Tx,'Y ~ 1');
    
    MeanTraitDel(i,1)=nanmean(Trait2-Trait1);    
    SETraitDel(i,1)=[nanstd(Trait2-Trait1)./sqrt(length(Trait2-Trait1))];
    UpperTraitDel(i,1)=mdlDel.Coefficients.Upper;
    LowerTraitDel(i,1)=mdlDel.Coefficients.Lower;
    pTraitDel(i,1)=mdlDel.Coefficients.pValue;
    
    
    [ICCy(1,i), LB(1,i), UB(1,i), F, df1, df2, p(1,i)] = ICC([Trait1 Trait2], '1-1')
    [ICCy(2,i), LB(2,i), UB(2,i), F, df1, df2, p(2,i)] = ICC([Trait1 Trait2], 'C-1')
    [ICCy(3,i), LB(3,i), UB(3,i), F, df1, df2, p(3,i)] = ICC([Trait1 Trait2], 'A-1') %A-1 A-k
    
    
%     continue
    %     T = table(Trait1,Trait2,deltaFlateral,deltaFN1,deltaFN3,deltaFREM,deltaBMI);
    T = table(Trait1,Trait2,deltaFlateral,deltaFN1,deltaFN3,deltaFREM);
    T2 = T;
    %     T2.deltaBMI = 0*T2.deltaBMI;
    
    disp(['trait:', traitslist{j}]);
    
    T.Trait2minusTrait1 = T.Trait2 - T.Trait1;
    T.deltaTime=deltaTime;
%     CenterPoint(i,1)=nanmedian(Trait1);
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
    
    
    
    mdl = fitglm(T,'Trait2 ~ Trait1')
    Ypred = predict(mdl,T2);
    [rho(i,1),pval(i,1)]=corr(T.Trait1,T.Trait2,'Type','Pearson','rows','complete')
    Trait1row = find(strcmp(mdl.Coefficients.Properties.RowNames,'Trait1'));
    BetaTrait12(i,1) = mdl.Coefficients.Estimate(Trait1row);
    SETrait12(i,1)= mdl.Coefficients.SE(Trait1row);
    PvalTrait12(i,1)=mdl.Coefficients.pValue(Trait1row);
    
%     mdl = fitglm(T,'Trait2 ~ Trait1 + deltaFlateral')
    mdl = fitglm(T,'Trait2minusTrait1 ~ deltaFlateral')
    Posrow = find(strcmp( mdl.Coefficients.Properties.RowNames,'deltaFlateral'));
    pValueAdjPos(i,1) = mdl.Coefficients.pValue(Posrow);
    Ypred = predict(mdl,T);
%     rhoAdjPos(i,1)=corr(Ypred,T.Trait2,'Type','Pearson','rows','complete')
    rhoAdjPos(i,1)=corr(T.Trait1+Ypred,T.Trait2,'Type','Pearson','rows','complete')
%     Trait1row = find(strcmp( mdl.Coefficients.Properties.RowNames,'Trait1'));
%     BetaAdjPos(i,1) = mdl.Coefficients.Estimate(Trait1row);
%     SEAdjPos(i,1)= mdl.Coefficients.SE(Trait1row);
%     PvalAdjPos(i,1)=mdl.Coefficients.pValue(Trait1row);
    
    
%     mdl = fitglm(T,'Trait2 ~ Trait1 + deltaFlateral+ deltaFREM')
    mdl = fitglm(T,'Trait2minusTrait1 ~ deltaFlateral+ deltaFREM')
    FREMrow = find(strcmp( mdl.Coefficients.Properties.RowNames,'deltaFREM'));
    pValueAdjFREM(i,1) = mdl.Coefficients.pValue(FREMrow);
    Ypred = predict(mdl,T);
    corr(T.Trait1,T.Trait2,'Type','Pearson','rows','complete')
    rhoAdjFREM(i,1)=corr(T.Trait1+Ypred,T.Trait2,'Type','Pearson','rows','complete')
%     Trait1row = find(strcmp( mdl.Coefficients.Properties.RowNames,'Trait1'))
%     BetaAdjFREM(i,1) = mdl.Coefficients.Estimate(Trait1row);
%     SEAdjFREM(i,1)= mdl.Coefficients.SE(Trait1row);
%     PvalAdjFREM(i,1)=mdl.Coefficients.pValue(Trait1row);
    rhoAdjFREM_temp=corr(T.Trait1,T.Trait2-Ypred,'Type','Pearson','rows','complete')
    
%     mdl = fitglm(T,'Trait2 ~ Trait1 + deltaFlateral+ deltaFREM +deltaFN1 + deltaFN3')
    mdl = fitglm(T,'Trait2minusTrait1 ~ deltaFlateral+ deltaFREM +deltaFN1 + deltaFN3')
    try mdlKeep{i}=mdl; end
    N1row = find(strcmp( mdl.Coefficients.Properties.RowNames,'deltaFN1'));
    N3row = find(strcmp( mdl.Coefficients.Properties.RowNames,'deltaFN3'));
    pValueAdjN1(i,1) = mdl.Coefficients.pValue(N1row);
    pValueAdjN3(i,1) = mdl.Coefficients.pValue(N3row);
    Ypred = predict(mdl,T); 
%     rhoAdjStates(i,1)=corr(Ypred,T.Trait2,'Type','Pearson','rows','complete')
     rhoAdjStates(i,1)=corr(T.Trait1+Ypred,T.Trait2,'Type','Pearson','rows','complete')
     
     mdl = fitglm(T,'Trait2minusTrait1 ~ deltaFlateral+ deltaFREM +deltaFN1 + deltaFN3 + deltaTime')
     Timerow = find(strcmp( mdl.Coefficients.Properties.RowNames,'deltaTime'));
     pValueAdjTime(i,1) = mdl.Coefficients.pValue(Timerow);
     Ypred = predict(mdl,T);
     rhoAdjTime(i,1)=corr(T.Trait1+Ypred,T.Trait2,'Type','Pearson','rows','complete')
    
    
    
     [ICCyAdj(1,i), LBAdj(1,i), UBAdj(1,i), F, df1, df2, pAdj(1,i)] = ICC([T.Trait1+Ypred,T.Trait2], '1-1');
    [ICCyAdj(2,i), LBAdj(2,i), UBAdj(2,i), F, df1, df2, pAdj(2,i)] = ICC([T.Trait1+Ypred,T.Trait2], 'C-1');    
    [ICCyAdj(3,i), LBAdj(3,i), UBAdj(3,i), F, df1, df2, pAdj(3,i)] = ICC([T.Trait1+Ypred,T.Trait2], 'A-1'); %A-1 A-k
    
    
    biasA = nanmean(T.Trait2minusTrait1)
    stdA = nanstd(T.Trait2minusTrait1)
    LofA_ = 1.96*stdA;
    LofA = biasA + stdA*1.96*[-1 0 1]
    
    
    
    %CR1 = 1.96*rms(T.Trait2minusTrait1)
    %https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0073990
    
    
    
%     T.Error = Ypred - T.Trait1;
%     corr(Ypred,T.Trait2,'Type','Pearson','rows','complete')
%     corr(T.Trait1,T.Trait2,'Type','Pearson','rows','complete')
%     corr(T.Trait1+T.Error,T.Trait2,'Type','Pearson','rows','complete')
%     corr(T.Trait1+T.Error/2,T.Trait2-T.Error/2,'Type','Pearson','rows','complete')
%     corr(T.Trait1,T.Trait2-T.Error,'Type','Pearson','rows','complete')
%     ICC([T.Trait1+T.Error/2,T.Trait2-T.Error/2], 'C-1')
%     ICC([T.Trait1,T.Trait2-T.Error], 'C-1')
%     ICC([T.Trait1+T.Error,T.Trait2], 'C-1')
%     ICC([T.Trait1,T.Trait2-T.Error], 'C-1')
%     ICC([T.Trait1,T.Trait2-T.Error], 'C-1')
%     figure(99)
%     subplot(1,4,1); plot(T.Trait1,T.Trait2,'.')
%     subplot(1,4,2); plot(T.Trait1+T.Error,T.Trait2,'.')
%     subplot(1,4,3); plot(T.Trait1+0.5*T.Error,T.Trait2-0.5*T.Error,'.')
%     subplot(1,4,4); plot(T.Trait1,T.Trait2-T.Error,'.')
%     ICC([T.Trait1,T.Trait2], 'C-1')
%     ICC([T.Trait2,Ypred], 'C-1')
%     SSdiff = sum([T.Trait1-T.Trait2].^2)
%     SSdiff2 = sum([Ypred-T.Trait2].^2)
%     SStot = 0.5*sum([T.Trait1-nanmean(T.Trait1)].^2) + sum([T.Trait2-nanmean(T.Trait2)].^2)
%     SStot2 = 0.5*sum([Ypred-nanmean(Ypred)].^2) + sum([T.Trait2-nanmean(T.Trait2)].^2)
%     
%     
    
    %mdl.Rsquared.Ordinary.^0.5
%     Trait1row = find(strcmp( mdl.Coefficients.Properties.RowNames,'Trait1'))
%     BetaAdjN1N3(i,1) = mdl.Coefficients.Estimate(Trait1row);
%     SEAdjN1N3(i,1)= mdl.Coefficients.SE(Trait1row);
%     PvalAdjN1N3(i,1)=mdl.Coefficients.pValue(Trait1row);
%     
%     
%     BetaAdjAll(:,i) = mdl.Coefficients.Estimate(2:end);
%     SEAdjAll(:,i) = mdl.Coefficients.SE(2:end);
    
    %mdl.Rsquared.Ordinary.^0.5
    traitslist1(i,1)=traitslist(j);
    
    
    %% plots
    thresgood1=[28 35 34 30 20]'; % from mesa criteria
    traitsj=[6 8 1 2 5]; % order in which traits used
    color1 = [0 0 0];
    color2 = [0.7 0.4 0.2];
    Ncommon=N1(crit)>thresgood1(i)& N2(crit)>thresgood1(i);
    
    CRcentile(i) = prctile(abs(T.Trait2minusTrait1(Ncommon)),95)
    
   % T2 = table();
   % T2.Trait = [T.Trait1;T.Trait2];
   % T2.Subj = nominal([[1:length(T.Trait1)]';[1:length(T.Trait1)]']);
   % fitlme(T2,'Trait ~ 1 + (1|Subj) ')
    
        % abs guarantees a real result
        bias=0;
        x=T.Trait2minusTrait1(Ncommon);
        CRgaussian(i) = 1.96*(nansum(abs(x - 0).^2) ./ (length(x)-1))^0.5; %0 describes center point for CI %absolute method
        
        CRgaussianconsistency(i) = 1.96*nanstd(x);
        
        NcommonTime1=Ncommon(deltaTime>median(deltaTime));
        x=T.Trait2minusTrait1(NcommonTime1);
        CRgaussianMedTimeabove(i) = 1.96*(nansum(abs(x - 0).^2) ./ (length(x)-1))^0.5; %0 describes center point for CI %absolute method
        
        NcommonTime2=Ncommon(deltaTime<median(deltaTime));
        x=T.Trait2minusTrait1(NcommonTime2);
        CRgaussianMedTimebelow(i) = 1.96*(nansum(abs(x - 0).^2) ./ (length(x)-1))^0.5; %0 describes center point for CI %absolute method
        
    %T.Trait12mean = (T.Trait1 + T.Trait2)/2;
    %
    %sqrt(2)^2*1.96*nanstd([T.Trait1-T.Trait12mean; T.Trait2-T.Trait12mean])
    
%     ft2 = fittype({'x'});
%     x=Trait1(Ncommon);
%     y=Trait2(Ncommon);
%     idx = ~isnan(x) & ~isnan(y);
%     p2 = fit(x(idx),y(idx),ft2);
%     x_fit = min(x):max(x);
%     y2_fitted = feval(p2,  x_fit);
%     plot(x_fit,y2_fitted,'b--');
    
    
    % xlim ylim for Lgn
    
    figure(1);
    % 2nd subplot
    upperNtrait=prctile(N1,99);
    subplot(3,length(traitsj)+1,i + (length(traitsj)+1));
    
    dStep=5;
    Centers=0:dStep:upperNtrait;
    Edges=(Centers(1)-dStep/2):dStep:(Centers(end)+dStep/2);
    Edges(end)=Inf;
    temp = N1;
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
    hold off
    
     
    % Bland-Altman
     BAMean=mean([Trait1(Ncommon),Trait2(Ncommon)],2,'omitnan');
     BADiff=Trait2(Ncommon)-Trait1(Ncommon);
     ZeroLine=zeros(length(BAMean),1);
     Meandiff=nanmean(BADiff);
     StdDev1=1.96*nanstd(BADiff);
     
     MeanLine=Meandiff*ones(length(BAMean),1);
     LOA(i,:)=[Meandiff+StdDev1 Meandiff-StdDev1];
     
    
     
     
     figure(1)
     set(gcf,'Position',get(0,'Screensize'));
    % 1st subplot
    
    subplot(3,length(traitsj)+1,i);
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
    if 1 %MrOS method
    plot(xlims,xlims+Meandiff,'b')
    plot(xlims,xlims+Meandiff+StdDev1,'r')
    plot(xlims,xlims+Meandiff-StdDev1,'r')
    else %MESA repeatability
    plot(xlims,xlims,'r-')
    plot(xlims,xlims+CRgaussian(i),'r')
    plot(xlims,xlims-CRgaussian(i),'r')        
    end 
    hold off

    
  
     subplot(3,length(traitsj)+1,i+length(traitsj)*2+2);%,i + (length(traitsj)+2)
     scatter(BAMean,BADiff,8,color1,'filled','markerfacealpha',0.5);
     hold on
     
     if 1 %MrOS method
         plot(BAMean,ZeroLine,'color', [.5 .5 .5],'LineWidth',1);
         plot(BAMean,MeanLine,'b','LineWidth',1);
         plot(BAMean,Meandiff+([1 -1]*StdDev1).*(ones(length(BAMean),2)),'r','LineWidth',1);
     else %MESA repeatability
         plot(BAMean,ZeroLine,'k','LineWidth',1);
         plot(BAMean,ZeroLine,'r','LineWidth',1);
         plot(BAMean,ZeroLine+([1 -1]*CRgaussian(i)).*(ones(length(BAMean),2)),'r','LineWidth',1);
     end 
%      xlabel('Mean'); ylabel('Difference');
%      title(traitslist(traitsj(i)))
    
    
     set(gca,'ylim',diff(xlims)/2*[-1 1])
     set(gca,'xlim',xlims)
     set(gca,'xtick',yticks);
     set(gca,'FontSize',12);
     hold off

    
%      BlandAltman(Trait1(Ncommon),Trait2(Ncommon),{traitslist(traitsj(i)),traitslist(traitsj(i))})
     
    
end
%% final tables - endotypes
Trho3 = table(traitslist1,rhoNREMexclusive,rho,rhoAdjPos,pValueAdjPos,rhoAdjFREM,pValueAdjFREM,rhoAdjStates,pValueAdjN1,pValueAdjN3,rhoAdjTime,pValueAdjTime)

try
    CRT=table([CRgaussian;CRgaussianMedTimebelow;CRgaussianMedTimeabove]);
    CRT.Properties.RowNames={'CRgaussian','CRgaussianMedTimebelow','CRgaussianMedTimeabove'};
    CRT
end

try
    TraitClassTFinal=array2table(TraitClassT,'VariableNames',traitslist1');
    TraitClassTFinal.Properties.RowNames={'V1Low','V1Med','V1High','V2Low','V2Med','V2High',...
        'V1Low-Low', 'V1Low-Med','V1Low-High',...
        'V1Med-Med', 'V1Med-Low','V1Med-High',...
        'V1High-High', 'V1High-Low','V1High-Med',...
        'V1Low-LowPer', 'V1Low-MedPer','V1Low-HighPer',...
        'V1Med-MedPer', 'V1Med-LowPer','V1Med-HighPer',...
        'V1High-HighPer', 'V1High-LowPer','V1High-MedPer'}
    
      
end

try
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
end
MeanTable=table(traitslist1,MeanTrait1,SETrait1,SDTrait1,MeanTrait2,SETrait2,SDTrait2,MeanTraitDel,SETraitDel,LowerTraitDel,UpperTraitDel,pTraitDel)

%% AHI %%
clear rhoNREMexclusiveAHI

Ahi3pArV1=Nsrr1.POOAHI3N(Nsrr1.Flag==1); % nrem ahi
Ahi3pArV2=Nsrr2.POOAHI3N;

traitslist2={'AHI'}

crit = OSA; % use only those subjects who has windows greater than threshold
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

clear Ahi3pArV1 Ahi3pArV2

Ahi3pArV1=Nsrr1.POOAHI3(Nsrr1.Flag==1);
Ahi3pArV2=Nsrr2.POOAHI3;


j=traitsj(3); % using flateral from LG1
crit = OSA ; % use only those subjects who has apneas
crit=OSA & TSTcrit;

deltaFlateral = Flateral2(crit,j) - Flateral1(crit,j);

deltaFN1 = FstatesArray_2{j}(crit,1) - FstatesArray_1{j}(crit,1); % change in n1
deltaFN3 = FstatesArray_2{j}(crit,3) - FstatesArray_1{j}(crit,3); % change in n3
deltaFREM = FstatesArray_2{j}(crit,4) - FstatesArray_1{j}(crit,4); % change in n3

%    deltaBMI = DeltaBMI(crit);

Trait2 = Ahi3pArV2(crit==1); % trait from visit2
Trait1 = Ahi3pArV1(crit==1); % trait from visit1


MeanTrait1=nanmean(Trait1);
SETrait1=[nanstd(Trait1)/sqrt(length(Trait1))];
 SDTrait1=[nanstd(Trait1)];
 
MeanTrait2=nanmean(Trait2);
SETrait2=[nanstd(Trait2)/sqrt(length(Trait2))];
 SDTrait2=[nanstd(Trait2)];
 
 Tx = table();
 Tx.Y = [Trait2-Trait1];
 mdlDel=fitglme(Tx,'Y ~ 1');
 
 MeanTraitDel=nanmean(Trait2-Trait1)
 
 SETraitDel=[nanstd(Trait2-Trait1)./sqrt(length(Trait2-Trait1))];
 UpperTraitDel=mdlDel.Coefficients.Upper;
 LowerTraitDel=mdlDel.Coefficients.Lower;
 pTraitDel=mdlDel.Coefficients.pValue;
 
MeanTraitDel=nanmean(Trait2-Trait1);
SETraitDel=[nanstd(Trait2-Trait1)./sqrt(length(Trait2-Trait1))];


T = table(Trait1,Trait2,deltaFlateral,deltaFN1,deltaFN3,deltaFREM);
T2 = T;
T.Trait2minusTrait1 = T.Trait2 - T.Trait1;

i=6;
% CenterPoint(i,1)=nanmedian(Trait1);

CenterPoint(i,1:3)= prctile(Trait1,[33.3333 66.6666 99.9999]);
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



%     T2.deltaBMI = 0*T2.deltaBMI;


mdl = fitglm(T,'Trait2 ~ Trait1')
Ypred = predict(mdl,T2);
rho=corr(T.Trait1,T.Trait2,'Type','Pearson','rows','complete')
Trait1row = find(strcmp(mdl.Coefficients.Properties.RowNames,'Trait1'));
BetaTrait12= mdl.Coefficients.Estimate(Trait1row);
SETrait12= mdl.Coefficients.SE(Trait1row);
PvalTrait12=mdl.Coefficients.pValue(Trait1row);

[ICCy(1,6), LB(1,6), UB(1,6), F, df1, df2, p(1,6)] = ICC([T.Trait1 T.Trait2], '1-1');
[ICCy(2,6), LB(2,6), UB(2,6), F, df1, df2, p(2,6)] = ICC([T.Trait1 T.Trait2], 'C-1');
[ICCy(3,6), LB(3,6), UB(3,6), F, df1, df2, p(3,6)] = ICC([T.Trait1 T.Trait2], 'A-1'); %A-1 A-k

 
%     mdl = fitglm(T,'Trait2 ~ Trait1 + deltaFlateral')
    mdl = fitglm(T,'Trait2minusTrait1 ~ deltaFlateral')
    Posrow = find(strcmp( mdl.Coefficients.Properties.RowNames,'deltaFlateral'));
    pValueAdjPos(6,1) = mdl.Coefficients.pValue(Posrow);
    Ypred = predict(mdl,T);
%     rhoAdjPos(i,1)=corr(Ypred,T.Trait2,'Type','Pearson','rows','complete')
    rhoAdjPos(6,1)=corr(T.Trait1+Ypred,T.Trait2,'Type','Pearson','rows','complete')


%     mdl = fitglm(T,'Trait2 ~ Trait1 + deltaFlateral+ deltaFREM')
    mdl = fitglm(T,'Trait2minusTrait1 ~ deltaFlateral+ deltaFREM')
    FREMrow = find(strcmp( mdl.Coefficients.Properties.RowNames,'deltaFREM'));
    pValueAdjFREM(6,1) = mdl.Coefficients.pValue(FREMrow);
    Ypred = predict(mdl,T);
    rhoAdjFREM(6,1)=corr(T.Trait1+Ypred,T.Trait2,'Type','Pearson','rows','complete')



    %     mdl = fitglm(T,'Trait2 ~ Trait1 + deltaFlateral+ deltaFREM +deltaFN1 + deltaFN3')
    mdl = fitglm(T,'Trait2minusTrait1 ~ deltaFlateral+ deltaFREM +deltaFN1 + deltaFN3')
    try mdlKeep{6}=mdl; end
    N1row = find(strcmp( mdl.Coefficients.Properties.RowNames,'deltaFN1'));
    N3row = find(strcmp( mdl.Coefficients.Properties.RowNames,'deltaFN3'));
    pValueAdjN1(6,1) = mdl.Coefficients.pValue(N1row);
    pValueAdjN3(6,1) = mdl.Coefficients.pValue(N3row);
    Ypred = predict(mdl,T); 
%     rhoAdjStates(i,1)=corr(Ypred,T.Trait2,'Type','Pearson','rows','complete')
     rhoAdjStates(6,1)=corr(T.Trait1+Ypred,T.Trait2,'Type','Pearson','rows','complete')
    
     [ICCyAdj(1,6), LBAdj(1,6), UBAdj(1,6), F, df1, df2, pAdj(1,6)] = ICC([T.Trait1+Ypred,T.Trait2], '1-1');
    [ICCyAdj(2,6), LBAdj(2,6), UBAdj(2,6), F, df1, df2, pAdj(2,6)] = ICC([T.Trait1+Ypred,T.Trait2], 'C-1');    
    [ICCyAdj(3,6), LBAdj(3,6), UBAdj(3,6), F, df1, df2, pAdj(3,6)] = ICC([T.Trait1+Ypred,T.Trait2], 'A-1'); %A-1 A-k
    

    

%mdl.Rsquared.Ordinary.^0.5
% Trait1row = find(strcmp( mdl.Coefficients.Properties.RowNames,'Trait1'))
% BetaAdjN1N3 = mdl.Coefficients.Estimate(Trait1row);
% SEAdjN1N3= mdl.Coefficients.SE(Trait1row);
% PvalAdjN1N3=mdl.Coefficients.pValue(Trait1row);
% 
% SEAdjAll = mdl.Coefficients.SE(2:end);
% BetaAdjAll = mdl.Coefficients.Estimate(2:end);


traitslist2={'AHI'}

%% plot
thresgood2=45;

Ncommon=TST1(crit)>thresgood2& TST2(crit)>thresgood2;

BAMean=mean([Trait1(Ncommon),Trait2(Ncommon)],2,'omitnan');
BADiff=Trait2(Ncommon)-Trait1(Ncommon);
ZeroLine=zeros(length(BAMean),1);
Meandiff=nanmean(BADiff);
StdDev1=1.96*nanstd(BADiff);

MeanLine=Meandiff*ones(length(BAMean),1);

i=6
LOA(i,:)=[Meandiff+StdDev1 Meandiff-StdDev1];


CRcentile(i) = prctile(abs(T.Trait2minusTrait1(Ncommon)),95)

% abs guarantees a real result
bias=0;
x=T.Trait2minusTrait1(Ncommon);
CRgaussian(i) = 1.96*(nansum(abs(x - 0).^2) ./ (length(x)-1))^0.5; %0 describes center point for CI %absolute method

CRgaussianconsistency(i) = 1.96*nanstd(x);
       

figure(1);

% 1st subplot
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
if 1 %MrOS method
    plot(xlims,xlims+Meandiff,'b')
    plot(xlims,xlims+Meandiff+StdDev1,'r')
    plot(xlims,xlims+Meandiff-StdDev1,'r')
else %MESA repeatability
    plot(xlims,xlims,'r-')
    plot(xlims,xlims+CRgaussian(i),'r')
    plot(xlims,xlims-CRgaussian(i),'r')
end
hold off

if 1
subplot(3,6,12);
dStep=15;
Centers=0:dStep:prctile(TST1,99);
Edges=(Centers(1)-dStep/2):dStep:(Centers(end)+dStep/2);
Edges(end)=Inf;
temp = TST1;
temp = temp(crit ); %& AHI>5

Ntemp=temp>thresgood2;
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




subplot(3,6,18);%,i + (length(traitsj)+2)
scatter(BAMean,BADiff,8,color1,'filled','markerfacealpha',0.5);
hold on

if 1 %MrOS method
    plot(BAMean,ZeroLine,'color', [.5 .5 .5],'LineWidth',1);
    plot(BAMean,MeanLine,'b','LineWidth',1);
    plot(BAMean,Meandiff+([1 -1]*StdDev1).*(ones(length(BAMean),2)),'r','LineWidth',1);
else %MESA repeatability
    plot(BAMean,ZeroLine,'k','LineWidth',1);
    plot(BAMean,ZeroLine,'r','LineWidth',1);
    plot(BAMean,ZeroLine+([1 -1]*CRgaussian(i)).*(ones(length(BAMean),2)),'r','LineWidth',1);
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

%% final tables
try
Trho3 = table(traitslist2,rhoNREMexclusive,rho,rhoAdjPos,pValueAdjPos,rhoAdjFREM,pValueAdjFREM,rhoAdjStates,pValueAdjN1,pValueAdjN3)
BetaTraitsAll=table(traitslist2,BetaTrait12,BetaNREMexclusive,BetaAdjPos,BetaAdjFREM,BetaAdjN1N3);
SETraitsAll=table(traitslist2,SETrait12,SENREMexclusive,SEAdjPos,SEAdjFREM,SEAdjN1N3);
PvalTraitsAll=table(traitslist2,PvalTrait12,PvalNREMexclusive,PvalAdjPos,PvalAdjFREM,PvalAdjN1N3);
end

try
    TraitClassTFinal=[TraitClassTFinal,table(TraitClassT(:,i),'VariableNames',{'AHI'})];
      
end

try
BetaAdjAllT=array2table(BetaAdjAll.');
BetaAdjAllT.Properties.VariableNames=mdl.Coefficients.Properties.RowNames(2:end);
BetaAdjAllT.Properties.RowNames=traitslist2;

SEAdjAllT=array2table(SEAdjAll.');
SEAdjAllT.Properties.VariableNames=mdl.Coefficients.Properties.RowNames(2:end);
SEAdjAllT.Properties.RowNames=traitslist2;

BetaAdjAllTpart = BetaAdjAllT(:,[2 5 3 4])
SEAdjAllTpart = SEAdjAllT(:,[2 5 3 4])
end
%Finding: Minimal effect of adj for position, except Vpassive

MeanTable=table(traitslist2,MeanTrait1,SDTrait1,MeanTrait2,SDTrait2,MeanTraitDel,SETraitDel,LowerTraitDel,UpperTraitDel,pTraitDel)

 



%% Unused
% Regression Model for NREM vs REM, NREM depths, and Position - adjusting for BMI

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