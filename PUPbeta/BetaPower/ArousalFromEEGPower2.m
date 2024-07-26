%% Start
addpath(genpath(pwd));
%%
load('BreathDataFullTableBackup2.mat')

%%
logit = @(p) log(p./(1-p));
logitinverse = @(p) 1./(1+exp(-p));
normalci = @(p,n) 1.96*((p.*(1-p)./n).^0.5);
addpath(genpath(pwd))

%% Overwriting directory, lazy

settings.OutputDataDirectory='G:\Dropbox (Personal)\PhenotypeDrive2018\Analyzed\'
    settings.savename='PhenoDrive2019';
    Npatients=50;
    
settings.OutputDataDirectory='G:\Dropbox (Personal)\RICCADSA 2\TraitsAnalysis\Analyzed\'
    settings.savename='RICCADSA'; %note powers were in mV2 and need multipling by 1000^2, i.e. +6 in log scale
    Npatients=222;
    
settings.OutputDataDirectory='G:\Dropbox (Personal)\SaraODB OAT\Analyzed\'
    settings.savename='SaraODB_OAT';
    Npatients=36;
    
settings.OutputDataDirectory='G:\Dropbox (Partners HealthCare)\MAD-OX\Traits\Analyzed\'
    settings.savename='MADOXBaseline'; 
    Npatients=37;
    

%%
clear BreathDataFullTable
%% Loop to make giant data table all patients all breaths, nonoverlapping

tic
for n=1:Npatients
    disp(['start: ' num2str(n)]);
    loadpath=[settings.OutputDataDirectory, settings.savename '_' num2str(n)];

    try
        load(loadpath,'BreathDataTable'); 
        [~,BreathDataTable2]=GetNonOvlappedVE(BreathDataTable);
    catch me
        disp(me.message);
        'no data, skipped'
        continue
    end
    
    clear temp
    for i=1:size(BreathDataTable2,1)
        temp{i,1} = [settings.savename '_' num2str(n)];
    end
    BreathDataTable2.Subj = temp;
     
    if exist('BreathDataFullTable')
        BreathDataFullTable = outerjoin(BreathDataFullTable,BreathDataTable2,'MergeKeys',true);
        %BreathDataFullTable = [BreathDataFullTable;BreathDataTable2];
    else
        BreathDataFullTable = [BreathDataTable2]; %start new table
    end

    disp(['end: ' num2str(n)]);
    toc
end

%% STEP A

BreathDataFullTableBackup=outerjoin(BreathDataFullTableBackup,BreathDataFullTable,'MergeKeys',true);

%%
% BreathDataFullTable1 = BreathDataFullTable;
% BreathDataFullTable2 = BreathDataFullTable;
% BreathDataFullTable=outerjoin(BreathDataFullTable1,BreathDataFullTable2,'MergeKeys',true);
% BreathDataFullTableBackup = BreathDataFullTable;
% BreathDataFullTableBackupShort = BreathDataFullTable;
% 

I = find(BreathDataFullTableBackup.Pbeta==-Inf);
length(I)
BreathDataFullTableBackup(I,:)=[];
I = find(isnan(BreathDataFullTableBackup.Pbeta));
length(I)
BreathDataFullTableBackup(I,:)=[];

I = find(BreathDataFullTable.Pbeta==-Inf);
length(I)
BreathDataFullTable(I,:)=[];
I = find(isnan(BreathDataFullTable.Pbeta));
length(I)
BreathDataFullTable(I,:)=[];

UniqueSubjList = unique(BreathDataFullTableBackup.Subj);

% %% Stats 1
% warning('off')
%

%%

BreathDataFullTable.hypnog_B(BreathDataFullTable.hypnog_B==-5)=NaN

%% Convert mV to uV in RICCADSA
ColsToAmplify={'Pdelta','Pbeta','Palpha','Ptheta','Psigma'};
StudyToAmplify='RICCADSA';
Headers = BreathDataFullTableBackup.Properties.VariableNames;

I=[];
for i=1:length(UniqueSubjList)
    I(i,1)=contains(string(UniqueSubjList{i}),string(StudyToAmplify));
end
I=find(I==1);


I=[];
for i=1:size(BreathDataFullTableBackup,1)
    I(i,1)=contains(string(BreathDataFullTableBackup.Subj{i}),string(StudyToAmplify));
end
I=find(I==1);


J=[];
for j=1:length(Headers)
    J(j,1)=sum(contains(string(Headers{j}),string(ColsToAmplify)))==1;
end
J=find(J==1);

if 1
    BreathDataFullTableBackup{I,J}=BreathDataFullTableBackup{I,J}+6;
end

%% STEP B
BreathDataFullTable=BreathDataFullTableBackup;


%% Remove subjects with negligible data in table
if 0
BreathDataFullTableBackup(strcmp(BreathDataFullTableBackup.Subj,'RICCADSA_114')==1,:)=[];
BreathDataFullTable(strcmp(BreathDataFullTable.Subj,'RICCADSA_114')==1,:)=[];
end

if 0 %attempt to handle duplicate import error
   temp = find(diff(subjrowrangeN)==0)
   temp = temp+1;
   Error = zeros(size(BreathDataFullTable,1),1);
   for i=1:length(temp)
        Error(subjrowrangeI{temp(i)})=1;
   end
   sum(Error)
   BreathDataFullTable(Error==1,:)=[];
   BreathDataFullTableBackup(Error==1,:)=[];
end

%% START Generate Centiles data (few min)
clear temp
NN=size(BreathDataFullTable,1);
    newcollist = {'Ptheta','Pbeta','Palpha','Pdelta'}
    centilexlist = [5 25 50 75 95];
    %centilexlist = [1:20 25 33 50 67 75 90];
    %centilexlist = [1 2 3 4 5 10 15 20 25 33 50 67 75 90];
for i=1:length(newcollist)
    for j=1:length(centilexlist)

newcol = newcollist{i};
    %newcol = 'Pdelta'
    centilex=centilexlist(j);
    newcol_ = [newcol num2str(centilex) 'p']
for n=1:length(UniqueSubjList)
    disp([newcol_ ' ' num2str(n)]);
    SubjTemp = UniqueSubjList{n};
    
    tempdata = eval(['BreathDataFullTable.' newcol]);
    rowrange = strcmp(BreathDataFullTable.Subj,SubjTemp)==1;
    temp = prctile(tempdata(rowrange,:),centilex);
    if any(strcmp(newcol_,BreathDataFullTable.Properties.VariableNames))==0
        nrow = size(BreathDataFullTable,1);
        eval(['BreathDataFullTable.' newcol_ '=zeros(nrow,1)+NaN;']);
    end
    eval(['BreathDataFullTable.' newcol_ '(rowrange)=temp;']);
end

    end
end

%% STEP 2: Make cols of prctile constants referenced to Pbeta5p
%beta + ref + others, no subj

%     newcollist = {'Pbeta','Palpha','Ptheta','Pdelta','Psigma'};
%     centilexlist = [1 2 5 10 15 20 25 33 50 67 75 90];
    
    newcollist = {'Pbeta','Pdelta','Palpha','Ptheta'}%,'Palpha','Pbeta_Pdelta'};
    %newcollist = {'Pbeta','Palpha','Pdelta'}%,'Palpha'};
    centilexlist = [5 25 50 75 95];
    %centilexlist = [5 50 95];
    
    count = 1
    clear constantslist
    constantsliststring=[];
    ref = {'Pbeta5p'}%,'Palpha1p','Palpha5p','Palpha25p','Palpha75p','Palpha10p','Palpha75p','Palpha90p','Ptheta','Pdelta','Psigma'};
        
for i=1:length(newcollist)
    for j=1:length(centilexlist)
        temp=[newcollist{i} num2str(centilexlist(j)) 'p'];
        if strcmp(temp,ref)==1
            continue
        end
        constantslist{count}=temp;
        newcol = constantslist(count);
        count=count+1;
        eval(['BreathDataFullTable.' ref{1} '_' newcol{1} '= BreathDataFullTable.' ref{1} '- BreathDataFullTable.' newcol{1} ';']);
        constantsliststring = [constantsliststring '+ ' ref{1} '_' newcol{1} ' '];
    end
end

clear constantslist2
for i=1:length(constantslist)
    constantslist2{1,i} = [ref{1} '_' constantslist{i}];
end

%overwrite to remove if useful:
%constantslist2 = [];%constantslist2([0]);

%% STEP 3: Reference all powers to Pbeta5p
newcol = {'Pbeta','Palpha','Pdelta','Ptheta','Psigma'};
refs = {'Pbeta5p'};
for j=1:length(refs)
for i=1:length(newcol)
    eval(['BreathDataFullTable.' newcol{i} '_' refs{j} '= BreathDataFullTable.' newcol{i} '- BreathDataFullTable.' refs{j} ';']);
end
end

%% STEP 4: run model with constants
%mdl3b = fitglm(BreathDataFullTable,['ARieF ~ Pbeta_Pbeta5p + Pdelta_Pbeta5p + Palpha_Pbeta5p + Ptheta_Pbeta5p + Psigma_Pbeta5p' constantsliststring],'Distribution','binomial','Link','logit','weights',weights)
%mdl3b = fitglm(BreathDataFullTable,['ARieF ~ Pbeta_Pbeta5p + Pdelta_Pbeta5p + Palpha_Pbeta5p + Psigma_Pbeta5p' constantsliststring],'Distribution','binomial','Link','logit','weights',weights)
%mdl3b = fitglm(BreathDataFullTable,['ARieF ~ Pbeta_Pbeta5p + Pdelta_Pbeta5p + Palpha_Pbeta5p' constantsliststring],'Distribution','binomial','Link','logit','weights',weights)

% Excl arousals
BreathDataFullTable.WakeNoAR = (BreathDataFullTable.hypnog_B==4)*1;
    BreathDataFullTable.WakeNoAR(isnan(BreathDataFullTable.WakeNoAR))=NaN;
    
advancedSWSbalance = 0;    
    if advancedSWSbalance
        temp = logitinverse(StageMeans2)
        temp2 = temp/(max(temp)-min(temp));
        temp2 = temp2-min(temp2) 
        %temp2/temp2(2)*0.67
        BreathDataFullTable.WakeNoAR(BreathDataFullTable.hypnog_B==4) = 1;
        BreathDataFullTable.WakeNoAR(BreathDataFullTable.hypnog_B==2) = 0.10;
        BreathDataFullTable.WakeNoAR(BreathDataFullTable.hypnog_B==3) = 0.07;
        BreathDataFullTable.WakeNoAR(BreathDataFullTable.hypnog_B==1) = 0.02;
        BreathDataFullTable.WakeNoAR(BreathDataFullTable.hypnog_B==0) = 0;
    end
    
    BreathDataFullTable.ExcludeAR = zeros(length(BreathDataFullTable.WakeNoAR),1);
    BreathDataFullTable.ExcludeAR(BreathDataFullTable.hypnog_B==4&BreathDataFullTable.ARieF<0.5)=1;
    BreathDataFullTable.ExcludeAR(BreathDataFullTable.hypnog_B~=4&BreathDataFullTable.ARieF>0.5)=1;
    BreathDataFullTable.ExcludeAR = logical(BreathDataFullTable.ExcludeAR);
 
    % Weights
    balancefactor=0.86;
balance = nanmean(BreathDataFullTable.WakeNoAR(~BreathDataFullTable.ExcludeAR))
weights = zeros(length(BreathDataFullTable.WakeNoAR),1);
weights(BreathDataFullTable.WakeNoAR>0.5)=1-balance*balancefactor;
weights(BreathDataFullTable.WakeNoAR<=0.5)=balance*balancefactor;
sum(weights(BreathDataFullTable.WakeNoAR>0.5&~BreathDataFullTable.ExcludeAR))
sum(weights(BreathDataFullTable.WakeNoAR<=0.5&~BreathDataFullTable.ExcludeAR))

if 1%advancedSWSbalance
    N_SWS = sum(BreathDataFullTable.hypnog_B==0 & ~BreathDataFullTable.ExcludeAR)
    %N_nonSWS = sum(BreathDataFullTable.hypnog_B>0 & BreathDataFullTable.hypnog_B<4 & ~BreathDataFullTable.ExcludeAR)
    N_W = sum(BreathDataFullTable.hypnog_B==4 & ~BreathDataFullTable.ExcludeAR)
    N_R = sum(BreathDataFullTable.hypnog_B==3 & ~BreathDataFullTable.ExcludeAR)
    N_N1 = sum(BreathDataFullTable.hypnog_B==2 & ~BreathDataFullTable.ExcludeAR)
    N_N2 = sum(BreathDataFullTable.hypnog_B==1 & ~BreathDataFullTable.ExcludeAR)
    
    N_total = N_N1 + N_N2 + N_R + N_SWS + N_W;
    F_SWS_actual = N_SWS / N_total;
    F_N1_actual = N_N1 / N_total;
    F_N2_actual = N_N2 / N_total;
    F_R_actual = N_R / N_total;
    F_W_actual = N_W / N_total;
    
    F_W_needed = 0.5;
    F_N1_needed = 0.125;
    F_N2_needed = 0.125;
    F_R_needed = 0.125;
    F_SWS_needed = 0.125;
    
    weights = zeros(length(BreathDataFullTable.WakeNoAR),1);
    weights(BreathDataFullTable.hypnog_B==0)=F_SWS_needed/F_SWS_actual;
    weights(BreathDataFullTable.hypnog_B==1)=F_N2_needed/F_N2_actual;
    weights(BreathDataFullTable.hypnog_B==2)=F_N1_needed/F_N1_actual;
    weights(BreathDataFullTable.hypnog_B==3)=F_R_needed/F_R_actual;
    weights(BreathDataFullTable.hypnog_B==4)=F_W_needed/F_W_actual;
    
end

meanbalancedY = nansum(weights(~BreathDataFullTable.ExcludeAR).*BreathDataFullTable.WakeNoAR(~BreathDataFullTable.ExcludeAR))/nansum(weights(~BreathDataFullTable.ExcludeAR))

if 0
    N_SWS = sum(BreathDataFullTable.hypnog_B==0 & ~BreathDataFullTable.ExcludeAR)
    N_N1 = sum(BreathDataFullTable.hypnog_B==2 & ~BreathDataFullTable.ExcludeAR)
    N_N2 = sum(BreathDataFullTable.hypnog_B==1 & ~BreathDataFullTable.ExcludeAR)
    N_W = sum(BreathDataFullTable.hypnog_B==4 & ~BreathDataFullTable.ExcludeAR)
    N_total = N_nonSWS + N_SWS + N_W + N_N2;
    F_SWS_actual = N_SWS / N_total;
    F_N1_actual = N_N1 / N_total;
    F_N2_actual = N_N2 / N_total;
    F_W_actual = N_W / N_total;
    
    F_W_needed = 0.25;
    F_SWS_needed = 0.25;
    F_N1_needed = 1 - F_W_needed - F_SWS_needed;
    
    weights = zeros(length(BreathDataFullTable.WakeNoAR),1);
    weights(BreathDataFullTable.hypnog_B==0 & ~BreathDataFullTable.ExcludeAR)=F_SWS_needed/F_SWS_actual;
    %weights(BreathDataFullTable.hypnog_B==1 & ~BreathDataFullTable.ExcludeAR)=F_N2_needed/F_N2_actual;
    weights(BreathDataFullTable.hypnog_B==2 & ~BreathDataFullTable.ExcludeAR)=F_N1_needed/F_N1_actual;
    weights(BreathDataFullTable.hypnog_B==4 & ~BreathDataFullTable.ExcludeAR)=F_W_needed/F_W_actual;
   
offset = thresFind; %logit(thresFind), then thresFind
BreathDataFullTable.WSC = NaN*zeros(length(BreathDataFullTable.ExcludeAR),1);
BreathDataFullTable.WSC(BreathDataFullTable.hypnog_B==4)=3 - offset;
BreathDataFullTable.WSC(BreathDataFullTable.hypnog_B==3)=-1.9 - offset;
BreathDataFullTable.WSC(BreathDataFullTable.hypnog_B==2)=-1.6 - offset;
BreathDataFullTable.WSC(BreathDataFullTable.hypnog_B==1)=-3.6 - offset;
BreathDataFullTable.WSC(BreathDataFullTable.hypnog_B==0)=-5 - offset;
mdlWSC = compact(fitglm(BreathDataFullTable,['WSC ~ Ptheta_Pbeta5p*Pdelta_Pbeta5p + Ptheta_Pbeta5p^2 + Pbeta_Pbeta5p + Palpha_Pbeta5p + Pdelta_Pbeta5p + Pdelta_Pbeta5p^2',constantsliststring],'Exclude',BreathDataFullTable.ExcludeAR,'weights',weights))
mdlWSC.Rsquared.Ordinary
end

%%
mdl3b = compact(fitglm(BreathDataFullTable,['WakeNoAR ~ Ptheta_Pbeta5p*Pdelta_Pbeta5p + Ptheta_Pbeta5p^2 + Pbeta_Pbeta5p + Palpha_Pbeta5p + Pdelta_Pbeta5p + Pdelta_Pbeta5p^2',constantsliststring],'Distribution','binomial','Link','logit','Exclude',BreathDataFullTable.ExcludeAR,'weights',weights))
mdl3b.Rsquared.Ordinary
% 
% %weights WSC (see below)
% mdl3bAlt = compact(fitglm(BreathDataFullTable,['WakeNoAR ~ Ptheta_Pbeta5p*Pdelta_Pbeta5p + Ptheta_Pbeta5p^2 + Pbeta_Pbeta5p + Palpha_Pbeta5p + Pdelta_Pbeta5p + Pdelta_Pbeta5p^2',constantsliststring],'Distribution','binomial','Link','logit','Exclude',BreathDataFullTable.ExcludeAR,'weights',weightsWSC))
% mdl3bAlt.Rsquared.Ordinary

%mdl3b = compact(fitglm(BreathDataFullTable,['ARieF ~ Ptheta_Pbeta5p*Pdelta_Pbeta5p + Ptheta_Pbeta5p^2 + Pbeta_Pbeta5p + Palpha_Pbeta5p + Pdelta_Pbeta5p + Pdelta_Pbeta5p^2',constantsliststring],'Distribution','binomial','Link','logit'))

%% STEP 5: convert coeffs for constants to generate single reference value
%mdl3b = fitglm(BreathDataFullTable,['ARieF ~ Pbeta_Pbeta5p^2 + Pdelta_Pbeta5p^2 + Palpha_Pbeta5p^2 + Ptheta_Pbeta5p^2 + Psigma_Psigma5p^2 + Ptheta_Pbeta5p^2 + Pbeta_Pbeta5p*Pdelta_Pbeta5p' constantsliststring],'Distribution','binomial','Link','logit','weights',weights)
%mdl3b=compact(mdl3b);
%mdlX = mdl3bAlt;

mdlX = mdl3b;
constant = 1*(sum((string(mdlX.Coefficients.Properties.RowNames) == string(constantslist2)),2)==1);
constant(1)=NaN;
I1 = find(constant==1);
constantnames = mdlX.Coefficients.Properties.RowNames(I1);
coeffsconstant = mdlX.Coefficients.Estimate(I1);
I2 = find(constant==0);
coeffs1 = mdlX.Coefficients.Estimate(I2);
terms = 0;
temptemp=1/sum(coeffs1)*coeffsconstant;
for i=1:length(I1)
    terms = terms - 1/sum(coeffs1)*coeffsconstant(i)*eval(['BreathDataFullTable.' constantnames{i}]);
end
BreathDataFullTable.ref = BreathDataFullTable.Pbeta5p + terms;

Coeffs = [1-sum(temptemp);temptemp];

clear Constants
for i=1:length(constantnames)
    Constants{i,1} = constantnames{i}(9:end);
end
Constants = [{'Pbeta5p'};Constants];
RefTable = table(Constants,Coeffs)
%How to apply:
ref_=0;
for i=1:length(RefTable.Constants)
    ref_ = ref_ + RefTable.Coeffs(i)*eval(['BreathDataFullTable.' RefTable.Constants{i}]);
end

listlist=[];
newcol = {'Pbeta','Palpha','Pdelta','Ptheta','Psigma'};
refs = {'ref'};
for j=1:length(refs)
for i=1:length(newcol)
    listlist{length(listlist)+1} = [newcol{i} '_' refs{j}];
    eval(['BreathDataFullTable.' newcol{i} '_' refs{j} '= BreathDataFullTable.' newcol{i} '- BreathDataFullTable.' refs{j} ';']);
end
end

%% STEP 5b: rerun using single ref (for presentation, communication), excl arousals

% mdl3bXrefTop4_full_noweights_nobetasq = compact(fitglm(BreathDataFullTable,['ARieF ~ Ptheta_ref*Pdelta_ref + Ptheta_ref^2 + Pbeta_ref + Palpha_ref + Pdelta_ref + Pdelta_ref^2'],'Distribution','binomial','Link','logit'))
% mdl3bXrefTop4_full_noweights_nobetasq.Rsquared.Ordinary

% BreathDataFullTable.WakeNoAR = (BreathDataFullTable.hypnog_B==4)*1;
% ExcludeAR = 0*BreathDataFullTable.WakeNoAR;
%     ExcludeAR(BreathDataFullTable.hypnog_B==4&BreathDataFullTable.ARieF==0)=1;
%     ExcludeAR(BreathDataFullTable.hypnog_B~=4&BreathDataFullTable.ARieF==1)=1;
%     ExcludeAR = logical(ExcludeAR);

%no weights
% mdl3bXrefTop4_full_noweights_nobetasq_W_noweights = compact(fitglm(BreathDataFullTable,['WakeNoAR ~ Ptheta_ref*Pdelta_ref + Ptheta_ref^2 + Pbeta_ref + Palpha_ref + Pdelta_ref + Pdelta_ref^2'],'Distribution','binomial','Link','logit','Exclude',BreathDataFullTable.ExcludeAR))
% mdl3bXrefTop4_full_noweights_nobetasq_W_noweights.Rsquared.Ordinary

%weights
mdl4Alt = compact(fitglm(BreathDataFullTable,['WakeNoAR ~ Ptheta_ref*Pdelta_ref + Ptheta_ref^2 + Pbeta_ref + Palpha_ref + Pdelta_ref + Pdelta_ref^2'],'Distribution','binomial','Link','logit','Exclude',BreathDataFullTable.ExcludeAR,'weights',weights))
mdl4Alt.Rsquared.Ordinary
% 
% mdl4Alt = compact(fitglm(BreathDataFullTable,['WakeNoAR ~ Ptheta_ref*Pdelta_ref + Ptheta_ref^2 + Pbeta_ref + Palpha_ref + Pdelta_ref + Pdelta_ref^2'],'Distribution','binomial','Link','logit','Exclude',BreathDataFullTable.ExcludeAR,'weights',weights))
% mdl4Alt.Rsquared.Ordinary


%% STEP 6: Performance W only
% Incl AR
% clear temp
% [ARieF_pred,ARieF_predCI] = predict(mdlA,BreathDataFullTable);
% temp = [BreathDataFullTable.ARieF ARieF_pred];
% thres = 0.5;
% performance = PredictiveValue(1*(BreathDataFullTable.ARieF>thres),1*(ARieF_pred>thres),BreathDataFullTable.ARieF)
% BreathDataFullTable.ARieF_predModel = ARieF_pred;

% Excl AR
%mdlA = mdl4Alt;
mdlA = mdlX;
%mdlA = mdlWSC;
mdlA.Rsquared.Ordinary;
clear temp
[ARieF_pred,ARieF_predCI] = predict(mdlA,BreathDataFullTable);

[x,y,t,AUC_WS,~] = perfcurve(1*(BreathDataFullTable.WakeNoAR(~BreathDataFullTable.ExcludeAR)>0.5),1*(ARieF_pred(~BreathDataFullTable.ExcludeAR)),1); %need to find the threshold value that gives the OPTROCPT!
[~,I]=max(y+(1-x));
thresFind=mean(t(I:(I+1)))
AUC_WS

BreathDataFullTable.SWSnotN1W = 1*(BreathDataFullTable.hypnog_B==0);
BreathDataFullTable.SWSnotN1W(isnan(BreathDataFullTable.hypnog_B)) = NaN;
BreathDataFullTable.SWSnotN1W(BreathDataFullTable.ExcludeAR==1) = NaN;
BreathDataFullTable.SWSnotN1W(BreathDataFullTable.hypnog_B==1) = NaN;


[x,y,t,AUC_SWS,~] = perfcurve(1*(BreathDataFullTable.SWSnotN1W(~isnan(BreathDataFullTable.SWSnotN1W))>0.5),ARieF_pred(~isnan(BreathDataFullTable.SWSnotN1W)),0); %need to find the threshold value that gives the OPTROCPT!
[~,I]=max(y+(1-x));
thresFindSWS=mean(t(I:(I+1)))
AUC_SWS

if 0
thres = 0.5;
else
thres=thresFind;
end      
performance = PredictiveValue(1*(BreathDataFullTable.WakeNoAR(~BreathDataFullTable.ExcludeAR)>0.5),1*(ARieF_pred(~BreathDataFullTable.ExcludeAR)>thres),BreathDataFullTable.WakeNoAR(~BreathDataFullTable.ExcludeAR))
sensspecave=0.5*(performance.Spec_sem_chance_p(1)+performance.Sens_sem_chance_p(1))

BreathDataFullTable.ARieF_predModel = ARieF_pred;
BreathDataFullTable.ARieF_pred_logit = logit(BreathDataFullTable.ARieF_predModel);
%BreathDataFullTable.ARieF_pred_logit = (BreathDataFullTable.ARieF_predModel);

%adjust score so threshold is 0.5
if 1
    BreathDataFullTable.ARieF_pred_logit = logit(BreathDataFullTable.ARieF_predModel) - logit(thres);
    BreathDataFullTable.ARieF_predModel = logitinverse(BreathDataFullTable.ARieF_pred_logit);
end

%% SKIP: Younes Lookup Table Method, acc=83.94
if 1
PbetaCut = prctile(BreathDataFullTable.Pbeta(~BreathDataFullTable.ExcludeAR),10:10:90);
PalphaCut = prctile(BreathDataFullTable.Palpha(~BreathDataFullTable.ExcludeAR),10:10:90);
PthetaCut = prctile(BreathDataFullTable.Ptheta(~BreathDataFullTable.ExcludeAR),10:10:90);
PdeltaCut = prctile(BreathDataFullTable.Pdelta(~BreathDataFullTable.ExcludeAR),10:10:90);

IDa = [sum(BreathDataFullTable.Pdelta>PdeltaCut,2), ...
                sum(BreathDataFullTable.Ptheta>PthetaCut,2), ...
                sum(BreathDataFullTable.Palpha>PalphaCut,2), ...
                sum(BreathDataFullTable.Pbeta>PbetaCut,2)];

ID = IDa(:,1)*1000 + IDa(:,2)*100 + IDa(:,3)*10 + IDa(:,4)*1;
IDs=(0:9999)';
num=NaN*ones(10000,1);
den=NaN*ones(10000,1);
for i=1:10000
    I=find(ID==IDs(i)&~BreathDataFullTable.ExcludeAR);
    den(i)=length(I);
    num(i)=sum(BreathDataFullTable.WakeNoAR(I));
    %num2(i)=sum(BreathDataFullTable.ARieF_predModel(I));
end

%raw (unbalanced) probability
PrWID_unbalanced = num./den; %no correction for missing lookups
PrWID_unbalanced(den<3)=NaN;

%balanced probabilities (i.e. modified based on prevalence so that cutoff for wake sleep is at Pr=0.5)
WonS_balance =  1/(1./balance - 1);
SonW = 1./PrWID_unbalanced - 1;
SonW_balanced = SonW*WonS_balance;
PrWID_balanced = 1./(1 + SonW_balanced);

%select balanced for further analysis
PrWID = PrWID_balanced;

Tx = table(IDs,PrWID,den);

%PrWIDmodel = num2./den; %no correction for missing lookups
%PrWIDmodel(den<10)=NaN;

IDsa=[floor(rem(IDs,10000)/1000) floor(rem(IDs,1000)/100) floor(rem(IDs,100)/10) floor(rem(IDs,10)/1)];
IDdelta=IDsa(:,1);
IDtheta=IDsa(:,2);
IDalpha=IDsa(:,3);
IDbeta=IDsa(:,4);

% Performance
BreathDataFullTable.ARieF_predID_MY = PrWID(ID+1);   
ARieF_pred = BreathDataFullTable.ARieF_predID_MY;
clear temp
thres = 0.5;
performanceIDMY = PredictiveValue(1*(BreathDataFullTable.WakeNoAR(~BreathDataFullTable.ExcludeAR)>thres),1*(ARieF_pred(~BreathDataFullTable.ExcludeAR)>thres),BreathDataFullTable.WakeNoAR(~BreathDataFullTable.ExcludeAR))

%Model error
PrWIDmodelSSE = nansum((num2.*(PrWIDmodel-PrWID).^2))./nansum(num2)
PrWIDmodelRsq = 1-nansum((num2.*(PrWIDmodel-PrWID).^2))./nansum((num2.*(PrWID-nanmean(PrWID)).^2))
end

%% STEP 7: RERUN_THETA Lookup Table Method, four, normalized by ref (acc=92.5p)
PdeltaCut = prctile(BreathDataFullTable.Pdelta_ref(~BreathDataFullTable.ExcludeAR),10:10:90);
PthetaCut = prctile(BreathDataFullTable.Ptheta_ref(~BreathDataFullTable.ExcludeAR),10:10:90);
PalphaCut = prctile(BreathDataFullTable.Palpha_ref(~BreathDataFullTable.ExcludeAR),10:10:90);
PbetaCut = prctile(BreathDataFullTable.Pbeta_ref(~BreathDataFullTable.ExcludeAR),10:10:90);

IDa = [         sum(BreathDataFullTable.Pdelta_ref>PdeltaCut,2), ...
                sum(BreathDataFullTable.Ptheta_ref>PthetaCut,2), ...
                sum(BreathDataFullTable.Palpha_ref>PalphaCut,2), ...
                sum(BreathDataFullTable.Pbeta_ref>PbetaCut,2)];

ID = IDa(:,1)*1000 + IDa(:,2)*100 + IDa(:,3)*10 + IDa(:,4)*1;
IDs=(0:9999)';

IDsa=[floor(rem(IDs,10000)/1000) floor(rem(IDs,1000)/100) floor(rem(IDs,100)/10) floor(rem(IDs,10)/1)];

num=NaN*ones(10000,1);
num2=NaN*ones(10000,1);
Fnumden=NaN*ones(10000,1);
den=NaN*ones(10000,1);
betaID = NaN*ones(10000,1);
alphaID = NaN*ones(10000,1);
thetaID = NaN*ones(10000,1);
deltaID = NaN*ones(10000,1);
for i=1:10000
    I=find(ID==IDs(i)&~BreathDataFullTable.ExcludeAR);
    den(i)=length(I);
    if 0
        num(i)=sum(BreathDataFullTable.WakeNoAR(I));
        
    else
        Fnumden(i) = nansum(BreathDataFullTable.WakeNoAR(I).*weights(I)) ./ nansum(weights(I));
        num(i)=den(i).*Fnumden(i);
    end
    num2(i)=sum(BreathDataFullTable.ARieF_predModel(I)); %for plots only, comparison with model
    betaID(i) = nanmedian(BreathDataFullTable.Pbeta_ref(I)); %for plots only
    alphaID(i) = nanmedian(BreathDataFullTable.Palpha_ref(I)); %for plots only
    thetaID(i) = nanmedian(BreathDataFullTable.Ptheta_ref(I)); %for plots only
    deltaID(i) = nanmedian(BreathDataFullTable.Pdelta_ref(I)); %for plots only
end

% Keep
DecilesLookupRef.PdeltaCut = PdeltaCut;
DecilesLookupRef.PthetaCut = PthetaCut;
DecilesLookupRef.PalphaCut = PalphaCut;
DecilesLookupRef.PbetaCut = PbetaCut;
DecilesLookupRef.Pr = num./den;
DecilesLookupRef.N = den;
DecilesLookupRef.ID = IDsa;
DecilesLookupRef.IDdelta=IDsa(:,1);
DecilesLookupRef.IDtheta=IDsa(:,2);
DecilesLookupRef.IDalpha=IDsa(:,3);
DecilesLookupRef.IDbeta=IDsa(:,4);

%raw (unbalanced) probability
PrWID_unbalanced = num./den; %no correction for missing lookups
PrWID_unbalanced(den<2)=NaN;

%balanced probabilities (i.e. modified based on prevalence so that cutoff for wake sleep is at Pr=0.5)
% WonS_balance =  1/(1./balance - 1);
% SonW = 1./PrWID_unbalanced - 1;
% SonW_balanced = SonW*WonS_balance;
% PrWID_balanced = 1./(1 + SonW_balanced);
% PrWID = PrWID_balanced;

%select balanced for further analysis
PrWID = PrWID_unbalanced; %%new

NemptyID = sum(isnan(PrWID))
Tx = table(IDs,PrWID,den);

PrWIDmodel = num2./den; %no correction for missing lookups
%PrWIDmodel(den<10)=NaN;

IDsa=[floor(rem(IDs,10000)/1000) floor(rem(IDs,1000)/100) floor(rem(IDs,100)/10) floor(rem(IDs,10)/1)];
IDdelta=IDsa(:,1);
IDtheta=IDsa(:,2);
IDalpha=IDsa(:,3);
IDbeta=IDsa(:,4);

% Performance
BreathDataFullTable.ARieF_predID = PrWID(ID+1);   
ARieF_pred = BreathDataFullTable.ARieF_predID;
clear temp
thres = 0.5;
performanceID = PredictiveValue(1*(BreathDataFullTable.WakeNoAR(~BreathDataFullTable.ExcludeAR)>thres),1*(ARieF_pred(~BreathDataFullTable.ExcludeAR)>thres),BreathDataFullTable.WakeNoAR(~BreathDataFullTable.ExcludeAR))

%ERRORS
PrWIDmodelSSE = nansum((num2.*(PrWIDmodel-PrWID).^2))./nansum(num2)
PrWIDmodelRsq = 1-nansum((num2.*(PrWIDmodel-PrWID).^2))./nansum((num2.*(PrWID-nanmean(PrWID)).^2))

%% STEP 8: RERUN_THETA Plots
%wake vs beta, 9 levels of delta; two levels of theta; three levels of theta
y=[3 8]; %alpha
x=[5]; %theta
figure(1); clf(1);
pauselength = 0.05;
for i=1:length(y)
    for w=[0:9]
        subplot(1,length(y),i);
        I = find(IDtheta==x&IDdelta==w&IDalpha==y(i));
        %I = find(IDdelta==w&IDtheta==x&IDbeta==z);
        IDs(I);
        PrWIDplot = PrWID(I);
        PrWIDplot2 = PrWIDmodel(I);
        figure(1);
        %color_ = [0.5 0.5 z/10];
        %xvals = 0:9;
        xvals = betaID(I);
        plot(xvals,[PrWIDplot],'b.-'); hold('on');
        plot(xvals,[PrWIDplot2],'r.-'); hold('on');
        pause(pauselength);
    end
    %xlabel('beta rank');
    xlabel('beta power');
    text(xvals(1),0.9,['alpha rank = ' num2str(y(i))]);
end

%wake vs beta, 9 levels of delta; mid range level of alpha; three levels of theta 
y=[5]; %alpha
x=[2 5]; %theta
figure(2); clf(2);
for i=1:length(x)
    for w=[0:9]
        subplot(1,length(x),i);
        I = find(IDtheta==x(i)&IDdelta==w&IDalpha==y);
        %I = find(IDdelta==w&IDtheta==x&IDbeta==z);
        IDs(I);
        PrWIDplot = PrWID(I);
        PrWIDplot2 = PrWIDmodel(I);
        figure(2);
        %xvals = 0:9;
        xvals = betaID(I);
        plot(xvals,[PrWIDplot],'b.-'); hold('on');
        plot(xvals,[PrWIDplot2],'r.-'); hold('on');
        pause(pauselength);
    end
    %xlabel('beta rank');
    xlabel('beta power');
    text(xvals(1),0.9,['theta rank = ' num2str(x(i))]);
end

%wake vs theta, 9 levels of beta; mid alpha, three levels of delta;
y=[6]; %alpha
w=[3 7]; %delta
figure(3); clf(3);
for i=1:length(w)
    for z=[0:9]
        subplot(1,length(w),i);
        I = find(IDbeta==z&IDdelta==w(i)&IDalpha==y);
        %I = find(IDdelta==w&IDtheta==x&IDbeta==z);
        IDs(I);
        PrWIDplot = PrWID(I);
        PrWIDplot2 = PrWIDmodel(I);
        %color_ = [0.5 0.5 z/10];
        %xvals = 0:9;
        xvals = thetaID(I);
        plot(xvals,[PrWIDplot],'b.-'); hold('on');
        plot(xvals,[PrWIDplot2],'r.-'); hold('on');
        pause(pauselength);
        
    end
    xlabel('theta power');
    text(xvals(1),0.9,['delta rank = ' num2str(w(i))]);
end

%wake vs delta, 9 levels of beta; mid alpha, three levels of theta;
y=[5]; %alpha
x=[2 8]; %theta
figure(4); clf(4);
for i=1:length(x)
    for z=[0:9]
        subplot(1,length(w),i);
        I = find(IDbeta==z&IDtheta==x(i)&IDalpha==y);
        %I = find(IDdelta==w&IDtheta==x&IDbeta==z);
        IDs(I);
        PrWIDplot = PrWID(I);
        PrWIDplot2 = PrWIDmodel(I);
        figure(4);
        %color_ = [0.5 0.5 z/10];
        %xvals = 0:9;
        xvals = deltaID(I);
        plot(xvals,[PrWIDplot],'b.-'); hold('on');
        plot(xvals,[PrWIDplot2],'r.-'); hold('on');
        pause(pauselength);
        
    end
    xlabel('delta power');
    text(1,0.9,['theta rank = ' num2str(x(i))]);
end


%% plot3d, individual data against beta, alpha, delta

figure(1); clf(1);
ds=200;
plotsize=2;
I = find(BreathDataFullTable.WakeNoAR==1&~BreathDataFullTable.ExcludeAR);
I=I(1:ds:end);
scatter3(BreathDataFullTable.Pdelta_ref(I),BreathDataFullTable.Palpha_ref(I),BreathDataFullTable.Pbeta_ref(I),plotsize,'filled')
hold('on')
%I = find(BreathDataFullTable.ARieF<0.1,Npoints);
I = find(BreathDataFullTable.WakeNoAR==0&~BreathDataFullTable.ExcludeAR);
I=I(1:ds:end);
scatter3(BreathDataFullTable.Pdelta_ref(I),BreathDataFullTable.Palpha_ref(I),BreathDataFullTable.Pbeta_ref(I),plotsize,'filled')
xlabel('delta');
ylabel('alpha');
zlabel('beta');
hold('off')

%% Surface4

Nlines=10;
figure(9); clf(9);

Pdelta_refX=prctile(BreathDataFullTable.Pdelta_ref,1:(99-1)/9:99);
Pbeta_refX=prctile(BreathDataFullTable.Pbeta_ref,1:(99-1)/9:99);
Palpha_refX=prctile(BreathDataFullTable.Palpha_ref,1:(99-1)/9:99);
Ptheta_refX=prctile(BreathDataFullTable.Ptheta_ref,1:(99-1)/9:99);

[Pdelta_ref,Pbeta_ref] = meshgrid(Pdelta_refX,Pbeta_refX);

Palpha_ref = Pbeta_ref*0 + prctile(BreathDataFullTable.Palpha_ref,5);
Ptheta_ref = Pbeta_ref*0 + prctile(BreathDataFullTable.Ptheta_ref,50);
TblNew = table(Palpha_ref(:),Pdelta_ref(:),Pbeta_ref(:),Ptheta_ref(:));
TblNew.Properties.VariableNames = {'Palpha_ref','Pdelta_ref','Pbeta_ref','Ptheta_ref'};
temp = reshape(predict(mdlA,TblNew),[Nlines Nlines]);
figure(9)
subplot(2,2,1)
surf(Pbeta_ref,Pdelta_ref,temp)
xlabel('beta')
ylabel('delta')

Palpha_ref = Pbeta_ref*0 + prctile(BreathDataFullTable.Palpha_ref,95);
Ptheta_ref = Pbeta_ref*0 + prctile(BreathDataFullTable.Ptheta_ref,50);
TblNew = table(Palpha_ref(:),Pdelta_ref(:),Pbeta_ref(:),Ptheta_ref(:));
TblNew.Properties.VariableNames = {'Palpha_ref','Pdelta_ref','Pbeta_ref','Ptheta_ref'};
temp = reshape(predict(mdlA,TblNew),[Nlines Nlines]);
figure(9)
subplot(2,2,2)
surf(Pbeta_ref,Pdelta_ref,temp)
xlabel('beta')
ylabel('delta')

Palpha_ref = Pbeta_ref*0 + prctile(BreathDataFullTable.Palpha_ref,50);
Ptheta_ref = Pbeta_ref*0 + prctile(BreathDataFullTable.Ptheta_ref,5);
TblNew = table(Palpha_ref(:),Pdelta_ref(:),Pbeta_ref(:),Ptheta_ref(:));
TblNew.Properties.VariableNames = {'Palpha_ref','Pdelta_ref','Pbeta_ref','Ptheta_ref'};
temp = reshape(predict(mdlA,TblNew),[Nlines Nlines]);
figure(9)
subplot(2,2,3)
surf(Pbeta_ref,Pdelta_ref,temp)
xlabel('beta')
ylabel('delta')
%
figure(10); clf(10);
[Ptheta_ref,Pbeta_ref] = meshgrid(Ptheta_refX,Pbeta_refX);
Palpha_ref = Pbeta_ref*0 + prctile(BreathDataFullTable.Palpha_ref,50);
Pdelta_ref = Pbeta_ref*0 + prctile(BreathDataFullTable.Pdelta_ref,5);
TblNew = table(Palpha_ref(:),Pdelta_ref(:),Pbeta_ref(:),Ptheta_ref(:));
TblNew.Properties.VariableNames = {'Palpha_ref','Pdelta_ref','Pbeta_ref','Ptheta_ref'};
temp = reshape(predict(mdlA,TblNew),[Nlines Nlines]);

subplot(1,2,1)
surf(Pbeta_ref,Ptheta_ref,temp)
xlabel('beta')
ylabel('theta')

[Ptheta_ref,Pbeta_ref] = meshgrid(Ptheta_refX,Pbeta_refX);
Palpha_ref = Pbeta_ref*0 + prctile(BreathDataFullTable.Palpha_ref,50);
Pdelta_ref = Pbeta_ref*0 + prctile(BreathDataFullTable.Pdelta_ref,95);
TblNew = table(Palpha_ref(:),Pdelta_ref(:),Pbeta_ref(:),Ptheta_ref(:));
TblNew.Properties.VariableNames = {'Palpha_ref','Pdelta_ref','Pbeta_ref','Ptheta_ref'};
temp = reshape(predict(mdlA,TblNew),[Nlines Nlines]);

subplot(1,2,2)
surf(Pbeta_ref,Ptheta_ref,temp)
xlabel('beta')
ylabel('theta')



%% AR detect, STEP 1. Calculate local average for a window, very slow to calculate


UniqueSubjList = unique(BreathDataFullTable.Subj);


BreathDataFullTable.ARieF_pred_localbaseline_1p=BreathDataFullTable.Time0*NaN;
BreathDataFullTable.ARieF_pred_localbaseline_5p=BreathDataFullTable.Time0*NaN;
BreathDataFullTable.ARieF_pred_localbaseline_10p=BreathDataFullTable.Time0*NaN;
BreathDataFullTable.ARieF_pred_localbaseline_25p=BreathDataFullTable.Time0*NaN;
BreathDataFullTable.ARieF_pred_localbaseline_50p=BreathDataFullTable.Time0*NaN;
BreathDataFullTable.ARieF_pred_localbaseline_75p=BreathDataFullTable.Time0*NaN;
BreathDataFullTable.ARieF_pred_localbaseline_95p=BreathDataFullTable.Time0*NaN;
%BreathDataFullTable.ARieF_pred_local2baseline_1p=BreathDataFullTable.Time0*NaN;
BreathDataFullTable.ARieF_pred_localbaseline_min=BreathDataFullTable.Time0*NaN;

BreathDataFullTable.Pdelta_localbaseline_50p=BreathDataFullTable.Time0*NaN;
BreathDataFullTable.Ptheta_localbaseline_50p=BreathDataFullTable.Time0*NaN;
BreathDataFullTable.Pbeta_localbaseline_50p=BreathDataFullTable.Time0*NaN;
BreathDataFullTable.Palpha_localbaseline_50p=BreathDataFullTable.Time0*NaN;

BreathDataFullTable.Pdelta_localbaseline_5p=BreathDataFullTable.Time0*NaN;
BreathDataFullTable.Ptheta_localbaseline_5p=BreathDataFullTable.Time0*NaN;
BreathDataFullTable.Pbeta_localbaseline_5p=BreathDataFullTable.Time0*NaN;
BreathDataFullTable.Palpha_localbaseline_5p=BreathDataFullTable.Time0*NaN;

BreathDataFullTable.Pdelta_localbaseline_25p=BreathDataFullTable.Time0*NaN;
BreathDataFullTable.Ptheta_localbaseline_25p=BreathDataFullTable.Time0*NaN;
BreathDataFullTable.Pbeta_localbaseline_25p=BreathDataFullTable.Time0*NaN;
BreathDataFullTable.Palpha_localbaseline_25p=BreathDataFullTable.Time0*NaN;

BreathDataFullTable.Pdelta_localbaseline_75p=BreathDataFullTable.Time0*NaN;
BreathDataFullTable.Ptheta_localbaseline_75p=BreathDataFullTable.Time0*NaN;
BreathDataFullTable.Pbeta_localbaseline_75p=BreathDataFullTable.Time0*NaN;
BreathDataFullTable.Palpha_localbaseline_75p=BreathDataFullTable.Time0*NaN;

BreathDataFullTable.Pdelta_localbaseline_95p=BreathDataFullTable.Time0*NaN;
BreathDataFullTable.Ptheta_localbaseline_95p=BreathDataFullTable.Time0*NaN;
BreathDataFullTable.Pbeta_localbaseline_95p=BreathDataFullTable.Time0*NaN;
BreathDataFullTable.Palpha_localbaseline_95p=BreathDataFullTable.Time0*NaN;

%ARieF_pred_logitX
BreathDataFullTable.ARieF_pred_localXbaseline_1p=BreathDataFullTable.Time0*NaN;

for i=1:length(UniqueSubjList)
    I = find(string(BreathDataFullTable.Subj)==string(UniqueSubjList{i}));
    tempunique = unique(BreathDataFullTable.Time0(I));
    i
    for j=1:length(tempunique)
        %I2 = find(BreathDataFullTable.Time0==tempunique(j)&string(BreathDataFullTable.Subj)==string(UniqueSubjList{i}));
        I2 = I(find(BreathDataFullTable.Time0(I)==tempunique(j)));
        %BreathDataFullTable.ARieF_pred_local2baseline_1p(I2)=prctile(BreathDataFullTable.ARieF_pred_logit2(I2),1);
        %BreathDataFullTable.ARieF_pred_localXbaseline_1p(I2)=prctile(BreathDataFullTable.ARieF_pred_logitX(I2),1);
%         %BreathDataFullTable.ARieF_pred_localbaseline_min(I2)=min(BreathDataFullTable.ARieF_pred_logit(I2));
%         %BreathDataFullTable.ARieF_pred_localbaseline_5p(I2)=prctile(BreathDataFullTable.ARieF_pred_logit(I2),5);
%         %BreathDataFullTable.ARieF_pred_localbaseline_10p(I2)=prctile(BreathDataFullTable.ARieF_pred_logit(I2),10);
%         %BreathDataFullTable.ARieF_pred_localbaseline_25p(I2)=prctile(BreathDataFullTable.ARieF_pred_logit(I2),25);
%         BreathDataFullTable.ARieF_pred_localbaseline_50p(I2)=prctile(BreathDataFullTable.ARieF_pred_logit(I2),50);
%         BreathDataFullTable.ARieF_pred_localbaseline_75p(I2)=prctile(BreathDataFullTable.ARieF_pred_logit(I2),75);
%         BreathDataFullTable.ARieF_pred_localbaseline_95p(I2)=prctile(BreathDataFullTable.ARieF_pred_logit(I2),95);
% BreathDataFullTable.Pdelta_localbaseline_50p(I2)=prctile(BreathDataFullTable.Pdelta_ref(I2),50);
% BreathDataFullTable.Pbeta_localbaseline_50p(I2)=prctile(BreathDataFullTable.Pbeta_ref(I2),50);
% BreathDataFullTable.Palpha_localbaseline_50p(I2)=prctile(BreathDataFullTable.Palpha_ref(I2),50);
% BreathDataFullTable.Ptheta_localbaseline_50p(I2)=prctile(BreathDataFullTable.Ptheta_ref(I2),50);

BreathDataFullTable.Pdelta_localbaseline_25p(I2)=prctile(BreathDataFullTable.Pdelta_ref(I2),25);
BreathDataFullTable.Pbeta_localbaseline_25p(I2)=prctile(BreathDataFullTable.Pbeta_ref(I2),25);
BreathDataFullTable.Palpha_localbaseline_25p(I2)=prctile(BreathDataFullTable.Palpha_ref(I2),25);
BreathDataFullTable.Ptheta_localbaseline_25p(I2)=prctile(BreathDataFullTable.Ptheta_ref(I2),25);

BreathDataFullTable.Pdelta_localbaseline_5p(I2)=prctile(BreathDataFullTable.Pdelta_ref(I2),5);
BreathDataFullTable.Pbeta_localbaseline_5p(I2)=prctile(BreathDataFullTable.Pbeta_ref(I2),5);
BreathDataFullTable.Palpha_localbaseline_5p(I2)=prctile(BreathDataFullTable.Palpha_ref(I2),5);
BreathDataFullTable.Ptheta_localbaseline_5p(I2)=prctile(BreathDataFullTable.Ptheta_ref(I2),5);

BreathDataFullTable.Pdelta_localbaseline_75p(I2)=prctile(BreathDataFullTable.Pdelta_ref(I2),75);
BreathDataFullTable.Pbeta_localbaseline_75p(I2)=prctile(BreathDataFullTable.Pbeta_ref(I2),75);
BreathDataFullTable.Palpha_localbaseline_75p(I2)=prctile(BreathDataFullTable.Palpha_ref(I2),75);
BreathDataFullTable.Ptheta_localbaseline_75p(I2)=prctile(BreathDataFullTable.Ptheta_ref(I2),75);

    end
end

%% AR detect, STEP 2. Run model

ExcludeW = BreathDataFullTable.hypnog_B==4;% & string(BreathDataFullTable.Subj)~=string(UniqueSubjList(1:10:100));

balanceAR = nanmean(BreathDataFullTable.ARieF(~ExcludeW))*1.5 %adjusted to balance FP and FN using *1.5
weightsAR = 0*BreathDataFullTable.ARieF;
weightsAR(BreathDataFullTable.ARieF>0.5)=1-balanceAR;
weightsAR(BreathDataFullTable.ARieF<=0.5)=balanceAR;

BreathDataFullTable.ARieF_pred_logit_last1 = [NaN;BreathDataFullTable.ARieF_pred_logit(1:end-1)];
BreathDataFullTable.ARieF_pred_logit_last2 = [NaN;NaN;BreathDataFullTable.ARieF_pred_logit(1:end-2)];
BreathDataFullTable.ARieF_pred_logit_last3 = [NaN;NaN;NaN;BreathDataFullTable.ARieF_pred_logit(1:end-3)];
BreathDataFullTable.ARieF_pred_logit_last4 = [NaN;NaN;NaN;NaN;BreathDataFullTable.ARieF_pred_logit(1:end-4)];
BreathDataFullTable.ARieF_pred_logit_last5 = [NaN;NaN;NaN;NaN;NaN;BreathDataFullTable.ARieF_pred_logit(1:end-5)];
BreathDataFullTable.ARieF_pred_logit_last6 = [NaN;NaN;NaN;NaN;NaN;NaN;BreathDataFullTable.ARieF_pred_logit(1:end-6)];
BreathDataFullTable.ARieF_pred_logit_last7 = [NaN;NaN;NaN;NaN;NaN;NaN;NaN;BreathDataFullTable.ARieF_pred_logit(1:end-7)];
BreathDataFullTable.ARieF_pred_logit_last8 = [NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;BreathDataFullTable.ARieF_pred_logit(1:end-8)];
BreathDataFullTable.ARieF_pred_logit_last9 = [NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;BreathDataFullTable.ARieF_pred_logit(1:end-9)];

BreathDataFullTable.ARieF_pred_logit_next1 = [BreathDataFullTable.ARieF_pred_logit(2:end);NaN];
BreathDataFullTable.ARieF_pred_logit_next2 = [BreathDataFullTable.ARieF_pred_logit(3:end);NaN;NaN];
BreathDataFullTable.ARieF_pred_logit_next3 = [BreathDataFullTable.ARieF_pred_logit(4:end);NaN;NaN;NaN];
BreathDataFullTable.ARieF_pred_logit_next4 = [BreathDataFullTable.ARieF_pred_logit(5:end);NaN;NaN;NaN;NaN];
BreathDataFullTable.ARieF_pred_logit_next5 = [BreathDataFullTable.ARieF_pred_logit(6:end);NaN;NaN;NaN;NaN;NaN];
BreathDataFullTable.ARieF_pred_logit_next6 = [BreathDataFullTable.ARieF_pred_logit(7:end);NaN;NaN;NaN;NaN;NaN;NaN];
BreathDataFullTable.ARieF_pred_logit_next7 = [BreathDataFullTable.ARieF_pred_logit(8:end);NaN;NaN;NaN;NaN;NaN;NaN;NaN];
BreathDataFullTable.ARieF_pred_logit_next8 = [BreathDataFullTable.ARieF_pred_logit(9:end);NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN];
BreathDataFullTable.ARieF_pred_logit_next9 = [BreathDataFullTable.ARieF_pred_logit(10:end);NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN];

% retrain
% Pbeta_localbaseline_75p  + Pbeta_localbaseline_25p + Pbeta_localbaseline_50p + Pbeta_localbaseline_5p + Palpha_localbaseline_5p + Palpha_localbaseline_25p + Palpha_localbaseline_50p + Palpha_localbaseline_75p + Pdelta_localbaseline_5p + Pdelta_localbaseline_25p + Pdelta_localbaseline_50p + 
% mdlAR_rebuild = compact(fitglm(BreathDataFullTable,...
%     ['ARieF2 ~ Palpha_ref*Pbeta_ref + Pdelta_ref*Ptheta_ref + Pdelta_ref^2 + Ptheta_ref^2+ Pbeta_ref_last1 + Pbeta_ref_last3 + Pbeta_ref_last5 + Pbeta_ref_last7 + Pbeta_ref_last9 + Pbeta_ref_next1 + Pbeta_ref_next3 + Pbeta_ref_next5 + Pbeta_ref_next7 + Pbeta_ref_next9'],'Distribution','binomial','Link','logit','Exclude',ExcludeW,'weights',weightsAR2))
% mdlAR_rebuild = compact(fitglm(BreathDataFullTable,...
%     ['ARieF ~ Palpha_ref*Pbeta_ref + Pdelta_ref*Ptheta_ref + Pdelta_ref^2 + Ptheta_ref^2+ Pbeta_ref_last1 + Pbeta_ref_last3 + Pbeta_ref_last5 + Pbeta_ref_last7 + Pbeta_ref_last9 + Pbeta_ref_next1 + Pbeta_ref_next3 + Pbeta_ref_next5 + Pbeta_ref_next7 + Pbeta_ref_next9'],'Distribution','binomial','Link','logit','Exclude',ExcludeW,'weights',weightsAR))
% mdlAR_rebuild = compact(fitglm(BreathDataFullTable,...
%     ['ARieF ~ ARieF_pred_logit + Pbeta_ref_last1 + Pbeta_ref_last3 + Pbeta_ref_last5 + Pbeta_ref_last7 + Pbeta_ref_last9 + Pbeta_ref_next1 + Pbeta_ref_next3 + Pbeta_ref_next5 + Pbeta_ref_next7 + Pbeta_ref_next9'],'Distribution','binomial','Link','logit','Exclude',ExcludeW,'weights',weightsAR))
mdlAR_rebuild = compact(fitglm(BreathDataFullTable,...
    ['ARieF ~ ARieF_pred_logit + ARieF_pred_logit_last1 + ARieF_pred_logit_last3 + ARieF_pred_logit_last5 + ARieF_pred_logit_last7 + ARieF_pred_logit_last9 + ARieF_pred_logit_next1 + ARieF_pred_logit_next3  + ARieF_pred_logit_next5  + ARieF_pred_logit_next7  + ARieF_pred_logit_next9'],'Distribution','binomial','Link','logit','Exclude',ExcludeW,'weights',weightsAR))
% mdlAR_rebuild = compact(fitglm(BreathDataFullTable,...
%     ['ARieF ~ ARieF_pred_logit + SWSpred_logit + ARieF_pred_logit_last1 + ARieF_pred_logit_last3 + ARieF_pred_logit_last5 + ARieF_pred_logit_last7 + ARieF_pred_logit_last9 + ARieF_pred_logit_next1 + ARieF_pred_logit_next3  + ARieF_pred_logit_next5  + ARieF_pred_logit_next7  + ARieF_pred_logit_next9'],'Distribution','binomial','Link','logit','Exclude',ExcludeW,'weights',weightsAR))
mdlAR_rebuild.Rsquared.Ordinary

% mdlAR_rebuild = compact(fitglm(BreathDataFullTable,...
%     ['ARieF2 ~ Palpha_ref*Pbeta_ref + Pdelta_ref*Ptheta_ref + Pdelta_ref^2 + Ptheta_ref^2'],'Distribution','binomial','Link','logit','Exclude',ExcludeW,'weights',weightsAR2)) %
% mdlAR_rebuild.Rsquared.Ordinary

newARpred_ = predict(mdlAR_rebuild,BreathDataFullTable);
BreathDataFullTable.newARpredF = newARpred_;
BreathDataFullTable.newARpred_logit = logit(newARpred_);
BreathDataFullTable.newARpred = 1*(newARpred_>0.5);
    BreathDataFullTable.newARpred(isnan(newARpred_)) = NaN;
thres=0.5;
performanceAR = PredictiveValue(1*(BreathDataFullTable.ARieF(~ExcludeW)>0.5),1*(newARpred_(~ExcludeW)>thres),BreathDataFullTable.ARieF(~ExcludeW))
%performanceX = PredictiveValue(1*(BreathDataFullTable.ARieF(~ExcludeW)>0.5),1*(BreathDataFullTable.ARieF_pred_logit(~ExcludeW)>-1.2954),BreathDataFullTable.ARieF(~ExcludeW))
%performanceAR2 = PredictiveValue(1*(BreathDataFullTable.ARieF2(~ExcludeW)>0.5),1*(newARpred_(~ExcludeW)>thres),BreathDataFullTable.ARieF2(~ExcludeW))

I1 = find(BreathDataFullTable.newARpredF>0.5&BreathDataFullTable.ARieF<0.5&~ExcludeW);
I2 = find(BreathDataFullTable.newARpredF<0.5&BreathDataFullTable.ARieF>0.5&~ExcludeW);
I3 = find(BreathDataFullTable.newARpredF<0.5&BreathDataFullTable.ARieF<0.5&~ExcludeW);
I4 = find(BreathDataFullTable.newARpredF>0.5&BreathDataFullTable.ARieF>0.5&~ExcludeW);
I5 = find(BreathDataFullTable.newARpredF>0.5&~ExcludeW);
I6 = find(BreathDataFullTable.ARieF>0.5&~ExcludeW);
VItest = [nanmean(BreathDataFullTable.VI(I1)) nanmean(BreathDataFullTable.VI(I2)) nanmean(BreathDataFullTable.VI(I3)) nanmean(BreathDataFullTable.VI(I4)) nanmean(BreathDataFullTable.VI(I5)) nanmean(BreathDataFullTable.VI(I6))]

[x,y,t,AUC_AR,~] = perfcurve(1*(BreathDataFullTable.ARieF(~ExcludeW)>0.5),1*(newARpred_(~ExcludeW)),1); %need to find the threshold value that gives the OPTROCPT!
[~,I]=max(y+(1-x));
thresFind=mean(t(I:(I+1)))
AUC_AR

%% combined WakeSleepArousal

BreathDataFullTable.WAS = max([BreathDataFullTable.newARpred_logit BreathDataFullTable.ARieF_pred_logit]')';

BreathDataFullTable.WASPr = logitinverse(BreathDataFullTable.WAS);
thres=0.5;
performanceWAS = PredictiveValue(1*(BreathDataFullTable.ARieF>0.5),1*(BreathDataFullTable.WASPr>thres),BreathDataFullTable.ARieF)
%performanceX = PredictiveValue(1*(BreathDataFullTable.ARieF(~ExcludeW)>0.5),1*(BreathDataFullTable.ARieF_pred_logit(~ExcludeW)>-1.2954),BreathDataFullTable.ARieF(~ExcludeW))
%performanceAR2 = PredictiveValue(1*(BreathDataFullTable.ARieF2(~ExcludeW)>0.5),1*(newARpred_(~ExcludeW)>thres),BreathDataFullTable.ARieF2(~ExcludeW))


%% AR hist
% mdl6=mdlAR_rebuild; 
% %[mdl5x.Coefficients.Estimate(2:end) mdlA.Coefficients.Estimate(2:end)]./sum(abs([mdl5x.Coefficients.Estimate(2:end) mdlA.Coefficients.Estimate(2:end)]))
% 
% 
% Imain = find(string(mdl5a2.Coefficients.Properties.RowNames)=="ARieF_pred_logit");
% coefmain = mdl6.Coefficients.Estimate(Imain);
% 
% 
% newARpred_ = predict(mdl6,BreathDataFullTable);
% %newARpred_ = predict(mdl5a_beta,BreathDataFullTable);
% 
% BreathDataFullTable.newARpredF = newARpred_;
% BreathDataFullTable.newARpred_logit = logit(newARpred_);
% BreathDataFullTable.newARpred = 1*(newARpred_>0.5);
%     BreathDataFullTable.newARpred(isnan(newARpred_)) = NaN;
% 
% thres=0.5;
% performance = PredictiveValue(1*(BreathDataFullTable.ARieF(~ExcludeW)>0.5),1*(newARpred_(~ExcludeW)>thres),BreathDataFullTable.ARieF(~ExcludeW))
% %performanceX = PredictiveValue(1*(BreathDataFullTable.ARieF(~ExcludeW)>0.5),1*(BreathDataFullTable.ARieF_pred_logit(~ExcludeW)>-1.2954),BreathDataFullTable.ARieF(~ExcludeW))

%plot histogram
figure(126); clf(126);
dStep=0.1;
Centers=-10:dStep:10;
Edges=(Centers(1)-dStep/2):dStep:(Centers(end)+dStep/2);
ax12(2)=subplot(2,1,1);
[h7,edges] = histcounts(logit(newARpred_(BreathDataFullTable.ARieF>0.99)),Edges);
[h8,edges] = histcounts(logit(newARpred_(BreathDataFullTable.ARieF<0.01)),Edges);
bar(Centers,h7/sum(h7),'EdgeAlpha',0,'FaceAlpha',0.3,'BarWidth',1);
hold('on');
bar(Centers,h8/sum(h8),'EdgeAlpha',0,'FaceAlpha',0.7,'BarWidth',1);


ax12(1)=subplot(2,1,2);
[h7,edges] = histcounts(logit(newARpred_(~ExcludeW&BreathDataFullTable.ARieF>0.5)),Edges);
[h8,edges] = histcounts(logit(newARpred_(~ExcludeW&BreathDataFullTable.ARieF<0.5)),Edges);
bar(Centers,h7/sum(h7),'EdgeAlpha',0,'FaceAlpha',0.3,'BarWidth',1);
hold('on');
bar(Centers,h8/sum(h8),'EdgeAlpha',0,'FaceAlpha',0.7,'BarWidth',1);

linkaxes(ax12,'x')


%% Relationship between scored arousals and flow during sleep
% 
% figure(99)
% scatter(logit(newARpred),BreathDataFullTable.VI,5,'filled','markerfacealpha',0.1)

BreathDataFullTable.ARieF_pred_logitB = (BreathDataFullTable.ARieF_pred_logit>-1.2954)*1;
VInewthres=1.5;

thres1=0.5;thres2=0.5;thres3=0.67;
performanceVI1_flowthres = PredictiveValue(1*(BreathDataFullTable.VI(~ExcludeW)>VInewthres),1*(BreathDataFullTable.ARieF(~ExcludeW)>thres1),BreathDataFullTable.VI(~ExcludeW))
performanceVI2_flowthres = PredictiveValue(1*(BreathDataFullTable.VI(~ExcludeW)>VInewthres),1*(BreathDataFullTable.ARieF_pred_logitB(~ExcludeW)>thres2),BreathDataFullTable.VI(~ExcludeW))
performanceVI3_flowthres = PredictiveValue(1*(BreathDataFullTable.VI(~ExcludeW)>VInewthres),1*(BreathDataFullTable.newARpredF(~ExcludeW)>thres3),BreathDataFullTable.VI(~ExcludeW))

%
mdlxx = compact(fitglm(BreathDataFullTable,'ARieF ~ VI','Exclude',ExcludeW,'Distribution','binomial','Link','logit','weights',weightsAR))
newARpredVI1 = predict(mdlxx,BreathDataFullTable);
thres=0.5
performanceVI1 = PredictiveValue(1*(BreathDataFullTable.ARieF(~ExcludeW)>0.5),1*(newARpredVI1(~ExcludeW)>thres),BreathDataFullTable.ARieF(~ExcludeW))
VIthres = -mdlxx.Coefficients.Estimate(1)/mdlxx.Coefficients.Estimate(2);
performanceVI1_flowthres = PredictiveValue(1*(BreathDataFullTable.VI(~ExcludeW)>1.5),1*(BreathDataFullTable.ARieF(~ExcludeW)>thres),BreathDataFullTable.VI(~ExcludeW))


mdlxx2 = compact(fitglm(BreathDataFullTable,'newARpred ~ VI','Exclude',ExcludeW,'Distribution','binomial','Link','logit','weights',weightsAR))
newARpredVI2 = predict(mdlxx2,BreathDataFullTable);
thres=0.5
performanceVI2 = PredictiveValue(1*(BreathDataFullTable.newARpred(~ExcludeW)>0.5),1*(newARpredVI2(~ExcludeW)>thres),BreathDataFullTable.newARpred(~ExcludeW))
VIthres2 = -mdlxx2.Coefficients.Estimate(1)/mdlxx2.Coefficients.Estimate(2);

performanceVI2_flowthres = PredictiveValue(1*(BreathDataFullTable.VI(~ExcludeW)>1.5),1*(BreathDataFullTable.newARpred(~ExcludeW)>thres),BreathDataFullTable.VI(~ExcludeW))


BreathDataFullTable.ARieF_pred_logitB = (BreathDataFullTable.ARieF_pred_logit>0.5)*1;
mdlxx3 = compact(fitglm(BreathDataFullTable,'ARieF_pred_logitB ~ VI','Exclude',ExcludeW,'Distribution','binomial','Link','logit','weights',weightsAR))
newARpredVI3 = predict(mdlxx3,BreathDataFullTable);
thres=0.5
performanceVI3 = PredictiveValue(1*(BreathDataFullTable.ARieF_pred_logitB(~ExcludeW)>0.5),1*(newARpredVI3(~ExcludeW)>thres),BreathDataFullTable.ARieF_pred_logitB(~ExcludeW))
VIthres3 = -mdlxx3.Coefficients.Estimate(1)/mdlxx3.Coefficients.Estimate(2);

performanceVI3_flowthres = PredictiveValue(1*(BreathDataFullTable.VI(~ExcludeW)>1.2),1*(BreathDataFullTable.ARieF_pred_logitB(~ExcludeW)>thres),BreathDataFullTable.VI(~ExcludeW))


%% Plot variability in Ref
dStep=0.1;
Centers=-1:dStep:2;
Edges=(Centers(1)-dStep/2):dStep:(Centers(end)+dStep/2);
[hall,edges] = histcounts((BreathDataFullTable.ref),Edges);
figure(20)
bar(Centers,hall,'EdgeAlpha',0,'FaceAlpha',0.3,'BarWidth',1);

%% PLOT TRACES
%figure(1); clf(1);
figure()
set(gcf,'color',[1 1 1]);
jj=179500; %952000 shows concordance with flow, better than ar; 456856 shows clearer diffs vss ORP (noweightsversion); 457110 false arousal
           %517110 arousals weak per method; local change detection may add
           % information better. Same with 641910 661910 
           % 651110,652110 low values in wake detection, ORP similar (depending on technique)
           % 10000 is a problem for MY method; 
           % 221000 works well 709910 722910
           % 9000 nice correlation with flow
           % 34000 highlights advantage of knowing underlying model, not just Pr; also MYID fail
           % concordance with flow 42000
           % 170000 is great %335000
           % 123000 great AR but not great W 1100500 24000 1185000 
           % 1204000 179500 (over)scored arousals/overshoots with minimal EEG changes 
           % 3000000 new (underscored)arousals 3120000, no arousals scored in REM 3140000
           % , 3179000, 3182000, 129000, 151000, 34000
           % difficult section: 1203950
           % 179500
           
plotW=1;   
plotARdelta=1;
plotMY=0;  
RefPlot=mode(BreathDataFullTable.ref(jj:jj+Nsampplot));
Nsampplot = 200;
ax1(1)=subplot(3,1,1);
stairs(BreathDataFullTable.Time_start(jj:jj+Nsampplot),BreathDataFullTable.ARieF(jj:jj+Nsampplot),'k');
hold('on');
stairs(BreathDataFullTable.Time_start(jj:jj+Nsampplot),BreathDataFullTable.hypnog_B(jj:jj+Nsampplot)*0.25+1.25,'k');
box('off');
set(gca,'xticklabels',[],'xcolor',[1 1 1])
if plotMY
stairs(BreathDataFullTable.Time_start(jj:jj+Nsampplot),BreathDataFullTable.ARieF_predID(jj:jj+Nsampplot),'g');
stairs(BreathDataFullTable.Time_start(jj:jj+Nsampplot),BreathDataFullTable.ARieF_predID_MY(jj:jj+Nsampplot),'b');
end
if plotW
stairs(BreathDataFullTable.Time_start(jj:jj+Nsampplot),BreathDataFullTable.ARieF_predModel(jj:jj+Nsampplot),'r');
end
if plotARdelta

stairs(BreathDataFullTable.Time_start(jj:jj+Nsampplot),BreathDataFullTable.newARpredF(jj:jj+Nsampplot),'g');
stairs(BreathDataFullTable.Time_start(jj:jj+Nsampplot),(BreathDataFullTable.WASPr(jj:jj+Nsampplot)),'b');
% stairs(BreathDataFullTable.Time_start(jj:jj+Nsampplot),logitinverse(BreathDataFullTable.SWSpred_logit(jj:jj+Nsampplot)),'r:');
end
% stairs(BreathDataFullTable.Time_start(1:jj),ARieF_predCI(1:jj,1),'color',[0.5 0.2 0.2]);
% stairs(BreathDataFullTable.Time_start(1:jj),ARieF_predCI(1:jj,2),'color',[0.5 0.2 0.2]);
hold('off')
ylim([-0.1 2.3])
logit = @(p) log(p./(1-p));
ax1(2)=subplot(3,1,2);
stairs(BreathDataFullTable.Time_start(jj:jj+Nsampplot),BreathDataFullTable.ARieF(jj:jj+Nsampplot),'k');
hold('on');
box('off')
set(gca,'xticklabels',[],'xcolor',[1 1 1])
if plotMY
stairs(BreathDataFullTable.Time_start(jj:jj+Nsampplot),logit(BreathDataFullTable.ARieF_predID(jj:jj+Nsampplot)),'g');
stairs(BreathDataFullTable.Time_start(jj:jj+Nsampplot),logit(BreathDataFullTable.ARieF_predID_MY(jj:jj+Nsampplot)),'b');
end
if plotW
stairs(BreathDataFullTable.Time_start(jj:jj+Nsampplot),logit(BreathDataFullTable.ARieF_predModel(jj:jj+Nsampplot)),'r');
end
if plotARdelta

stairs(BreathDataFullTable.Time_start(jj:jj+Nsampplot),logit(BreathDataFullTable.newARpredF(jj:jj+Nsampplot)),'g');
stairs(BreathDataFullTable.Time_start(jj:jj+Nsampplot),logit(BreathDataFullTable.WASPr(jj:jj+Nsampplot)),'b');
% stairs(BreathDataFullTable.Time_start(jj:jj+Nsampplot),(BreathDataFullTable.SWSpred_logit(jj:jj+Nsampplot)),'r:');
end

ax1(3)=subplot(3,1,3);
stairs(BreathDataFullTable.Time_start(jj:jj+Nsampplot),BreathDataFullTable.VI(jj:jj+Nsampplot),'g-','linewidth',1.5);
hold('on');
if 0
stairs(BreathDataFullTable.Time_start(jj:jj+Nsampplot),BreathDataFullTable.Pdelta_ref(jj:jj+Nsampplot),'k-');
stairs(BreathDataFullTable.Time_start(jj:jj+Nsampplot),BreathDataFullTable.Ptheta_ref(jj:jj+Nsampplot),'b-');
stairs(BreathDataFullTable.Time_start(jj:jj+Nsampplot),BreathDataFullTable.Palpha_ref(jj:jj+Nsampplot),'r-');
stairs(BreathDataFullTable.Time_start(jj:jj+Nsampplot),BreathDataFullTable.Pbeta_ref(jj:jj+Nsampplot),'k-');
end
linkaxes(ax1,'x');
box('off')
%%
save workspace20190610c -v7.3



%% FrWakeSleep (~2 min)

newcol_ = {'FWakeNoAR','FWake','FARinSleep','FARieF_pred_logit','meantemp_ARieF_predModel','meantemp_ARieF_pred_logit'};
subjrowrangeI=[];
subjrowrangeN=[];
for n=1:length(UniqueSubjList)
    n
    
    
    li = SubjStartRow(n);
    if n<length(UniqueSubjList)
        ri = SubjStartRow(n+1)-1;
    else
        ri = length(UniqueSubjList)
    end
    subjrowrange = li:ri;
    
    temp_WakeNoAR = BreathDataFullTable.WakeNoAR(subjrowrange);
    temp_ExcludeAR = BreathDataFullTable.ExcludeAR(subjrowrange);
    temp_Wake = 1*(BreathDataFullTable.hypnog_B(subjrowrange)==4);
        temp_Wake(isnan(BreathDataFullTable.hypnog_B(subjrowrange)))=NaN;
        temp_Wake(BreathDataFullTable.hypnog_B(subjrowrange)>4)=NaN;
    FWakeNoAR = nansum(temp_WakeNoAR(~temp_ExcludeAR)>0.5)/nansum(~isnan(temp_WakeNoAR(~temp_ExcludeAR)));
    temp_AR = BreathDataFullTable.ARieF(subjrowrange);
    temp_ARinSleep = temp_AR>0.5&temp_Wake==0;
    FARinSleep = nansum(temp_ARinSleep==1)/nansum(~isnan(temp_ARinSleep));
    FWake = nansum(temp_Wake>0.5)/nansum(~isnan(temp_Wake));
    
    temp_WakeModel = 1*(BreathDataFullTable.hypnog_B(subjrowrange)==4);
    temp_ARieF_pred_logit = BreathDataFullTable.ARieF_pred_logit(subjrowrange);
    FARieF_pred_logit = nansum(temp_ARieF_pred_logit>0)/nansum(~isnan(temp_ARieF_pred_logit));
    
    temp_ARieF_predModel = BreathDataFullTable.ARieF_predModel(subjrowrange);
    meantemp_ARieF_predModel = nanmean(temp_ARieF_predModel);
    
    temp_ARieF_pred_logit = BreathDataFullTable.ARieF_pred_logit(subjrowrange);
    meantemp_ARieF_pred_logit = nanmean(temp_ARieF_pred_logit);
    %rowrangesubset = strcmp(BreathDataFullTable.Subj,SubjTemp)==1;
    %temp = prctile(tempdata,centilex);
    for j=1:length(newcol_)
    if any(strcmp(newcol_{j},BreathDataFullTable.Properties.VariableNames))==0
        nrow = size(BreathDataFullTable,1);
        eval(['BreathDataFullTable.' newcol_{j} '=zeros(nrow,1)+NaN;']);
    end
    eval(['BreathDataFullTable.' newcol_{j} '(subjrowrange)=' newcol_{j} ';']);
    end
end

%% Compare to ref

figure(121)
scatter(BreathDataFullTable.FWakeNoAR(SubjStartRow),BreathDataFullTable.Pbeta5p(SubjStartRow),10,'filled','markerfacealpha',0.5)
figure(122)
scatter(BreathDataFullTable.FWakeNoAR(SubjStartRow),BreathDataFullTable.ref(SubjStartRow),10,'filled','markerfacealpha',0.5)

figure(13)
scatter(BreathDataFullTable.Pbeta5p(SubjStartRow),BreathDataFullTable.ref(SubjStartRow),10,'filled','markerfacealpha',0.5)
figure(132)
scatter(BreathDataFullTable.Palpha25p(SubjStartRow),BreathDataFullTable.ref(SubjStartRow),10,'filled','markerfacealpha',0.5)

%%
figure(14)
%{'FWakeNoAR','FWake','FARinSleep','FARieF_pred_logit','meantemp_ARieF_predModel'};
subplot(2,2,1);
scatter(BreathDataFullTable.FWake(SubjStartRow),BreathDataFullTable.FARieF_pred_logit(SubjStartRow),10,'filled','markerfacealpha',0.5)

subplot(2,2,2);
scatter(BreathDataFullTable.FWake(SubjStartRow) + BreathDataFullTable.FARinSleep(SubjStartRow),BreathDataFullTable.FARieF_pred_logit(SubjStartRow),10,'filled','markerfacealpha',0.5)

subplot(2,2,3);
scatter(BreathDataFullTable.FWake(SubjStartRow) + BreathDataFullTable.FARinSleep(SubjStartRow),BreathDataFullTable.meantemp_ARieF_pred_logit(SubjStartRow),10,'filled','markerfacealpha',0.5)

%% Predicting the future? Works within N2
figure(16)


delta_=[14];
NbreathsInFuture_= 2 + 0*delta_;

VarX_ = BreathDataFullTable.ARieF_pred_logit;
% VarX_ = logitinverse(BreathDataFullTable.ARieF_predID);
% VarX_ = BreathDataFullTable.ARieF_predID;
% VarX_ = BreathDataFullTable.ARieF_predModel;
% VarX_ = BreathDataFullTable.hypnog_B;

for j=1:length(NbreathsInFuture_)
    NbreathsInFuture=NbreathsInFuture_(j); %15
    delta = delta_(j);
    %BreathDataFullTable.ARieF
    
    I = find((BreathDataFullTable.hypnog_B==1)&(BreathDataFullTable.ARieF==0)==1);
    
    %remove breaths out of range
    I(find((I + NbreathsInFuture + delta)>size(BreathDataFullTable,1)))=[];
    
    %remove breaths if different subjects
    I((find(strcmp(BreathDataFullTable.Subj(I),BreathDataFullTable.Subj(I+NbreathsInFuture))==0)))=[];
    
    if 1 %also remove breaths where the next 1 is arousal from assessment
        I((find(BreathDataFullTable.ARieF(I+1)>0)))=[];
        I((find(BreathDataFullTable.ARieF(I+2)>0)))=[];
        %I((find(BreathDataFullTable.ARieF(I+3)==0)))=[];
    end
    
    %remove breaths if clearly beyond expected level
    I((find(BreathDataFullTable.Time_start(I+NbreathsInFuture)>BreathDataFullTable.Time_start(I)+300)))=[];
    
    VarX = VarX_(I);
    
    currentAR = BreathDataFullTable.ARieF(I);
    futureAR = BreathDataFullTable.ARieF(I+NbreathsInFuture);
    futureAR2 = 1*(sum(BreathDataFullTable.ARieF(I+[NbreathsInFuture:NbreathsInFuture+delta])>0.5,2)>0);
    
    deltaT1 = nanmedian(BreathDataFullTable.Time_start(I+NbreathsInFuture)-BreathDataFullTable.Time_start(I))
    deltaT2 = nanmedian(BreathDataFullTable.Time_end(I+NbreathsInFuture+delta)-BreathDataFullTable.Time_start(I))
    
    PlogitCut = prctile(VarX,5:5:95); %5:5:95
    PlogitCutScore = 1+sum(VarX>PlogitCut,2);
    
    num=NaN*ones(length(PlogitCut)+1,1);
    xval=NaN*ones(length(PlogitCut)+1,1);
    den=NaN*ones(length(PlogitCut)+1,1);
    for i=1:(length(PlogitCut)+1)
        Ix=find(PlogitCutScore==i);
        den(i)=sum(~isnan(futureAR(Ix)));
        num(i)=sum(futureAR2(Ix)>0.5);
        xval(i)=median(VarX(Ix));
        %num2(i)=sum(BreathDataFullTable.ARieF_predModel(I));
    end
    
    normalci_ = normalci(num./den,den);
    
    
    plot(xval,num./den,'k.-','markersize',12);
    hold('on')
    plot(xval,num./den-normalci_,'-','color',[0.5 0.5 0.5]);
    plot(xval,num./den+normalci_,'-','color',[0.5 0.5 0.5]);
    %plot(logitinverse(xval),num./den,'.-') 
end


%%
mdl1test = compact(fitglm(BreathDataFullTable,['WakeNoAR ~ ARieF_pred_logit'],'Distribution','binomial','Link','logit','Exclude',BreathDataFullTable.ExcludeAR)); %'weights',weights
mdl1 = fitglme(BreathDataFullTable,['WakeNoAR ~ ARieF_pred_logit + (1|Subj)'],'Link','logit','Exclude',BreathDataFullTable.ExcludeAR);
mdl1test = compact(fitglm(BreathDataFullTable,['WakeNoAR ~ Pbeta_ref'],'Distribution','binomial','Link','logit','Exclude',BreathDataFullTable.ExcludeAR)); %'weights',weights
mdl1 = fitglme(BreathDataFullTable,['WakeNoAR ~ Pbeta_ref + (1|Subj)'],'Link','logit','Exclude',BreathDataFullTable.ExcludeAR);

%mdlx = fitglme(BreathDataFullTable,'ARieF ~ Pbeta + (1|Subj)','Link','logit')

BreathDataFullTable.ExcludeAR = ExcludeAR;

    %%
    %tempTable1 = BreathDataFullTable(~BreathDataFullTable.ExcludeAR,:);
    tabletemp = tempTable1(1:10:100000,{'WakeNoAR','ARieF_pred_logit','Subj','Pbeta_ref'});
    tabletemp.Subj = nominal(tabletemp.Subj);
mdl1_ = fitglme(tabletemp,['WakeNoAR ~ Pbeta_ref + (1|Subj)'],'Link','logit')

mdl1 = fitglme(tempTable1,['WakeNoAR ~ ARieF_pred_logit + (1|Subj)'],'Link','logit')

%%

% mdl1 = compact(mdl1);
REs = randomEffects(mdl1);

temp = BreathDataFullTable.Subj;
temp2 = ['NaN';BreathDataFullTable.Subj(1:end-1)];
I=find(strcmp(temp,temp2)==0);

BreathDataFullTable.Subj(I);
BreathDataFullTable_ = BreathDataFullTable(I,:);
BreathDataFullTable_.REs = REs;

mdlthres = fitglm(BreathDataFullTable_,'REs ~ Pbeta10p + Pbeta5p + Palpha25p + Palpha5p')
mdlthres.Rsquared.Ordinary
BreathDataFullTable.PowerRef=predict(mdlthres,BreathDataFullTable);


%% IndividualSubject WS Data plots
clear temp
%newcol_ = 'FWakeNoAR';
%subjrowrangeI=[];
%subjrowrangeN=[];
nrows=4;
ncols=8;
clear cutoff_NW_NS
plothist1=1
for n=1:length(UniqueSubjList)
    n
    SubjTemp = UniqueSubjList{n};
    subjrowrange = strcmp(BreathDataFullTable.Subj,SubjTemp)==1;
    %temp_WakeNoAR = BreathDataFullTable.WakeNoAR(strcmp(BreathDataFullTable.Subj,SubjTemp)==1);
    %temp_ExcludeAR = BreathDataFullTable.ExcludeAR(strcmp(BreathDataFullTable.Subj,SubjTemp)==1);
    balanceInd = nanmean(BreathDataFullTable.WakeNoAR(~BreathDataFullTable.ExcludeAR&subjrowrange==1));
    
    weightsInd = 0*BreathDataFullTable.WakeNoAR(subjrowrange==1);
    IW = BreathDataFullTable.WakeNoAR(subjrowrange==1)>0.5 & BreathDataFullTable.ExcludeAR(subjrowrange==1)==0;
    IS = BreathDataFullTable.WakeNoAR(subjrowrange==1)<=0.5 & BreathDataFullTable.ExcludeAR(subjrowrange==1)==0;
    %check
    
%     IW = logical([1 0 0])';
%     IS = logical([0 1 1])';
    %balanceInd_ = sum(IW)/(sum(IW)+sum(IS));
    
    weightsInd(IW)=1-balanceInd;
    weightsInd(IS)=balanceInd;
    
    %check
    sum(weightsInd(IW));
    sum(weightsInd(IS));
        
    mdltemp = compact(fitglm(BreathDataFullTable(subjrowrange==1,:),['WakeNoAR ~ ARieF_pred_logit'],'Distribution','binomial','Link','logit','Exclude',BreathDataFullTable.ExcludeAR(subjrowrange==1),'weights',weightsInd));
    
    cutoff_NW_NS(n,:) = [-mdltemp.Coefficients.Estimate(1)/mdltemp.Coefficients.Estimate(2) sum(IW) sum(IS) balanceInd];
    if plothist1
    %figure(99 + floor((n-1)/(nrows*ncols)))
    figure(99)
    subplot(nrows,ncols,mod(n-1,nrows*ncols)+1)
    dStep=0.25;
    Centers=-10:dStep:10;
    Edges=(Centers(1)-dStep/2):dStep:(Centers(end)+dStep/2);

    [h1,edges] = histcounts(logit(BreathDataFullTable.ARieF_predModel(subjrowrange==1&BreathDataFullTable.hypnog_B==4&BreathDataFullTable.ARieF>0.95)),Edges); 
    [h2,edges] = histcounts(logit(BreathDataFullTable.ARieF_predModel(subjrowrange==1&BreathDataFullTable.hypnog_B<4&BreathDataFullTable.ARieF<0.05)),Edges); 

%stairs(Centers-dStep/2,h1);
%stairs(Centers-dStep/2,h2);
bar(Centers,h1,'EdgeAlpha',0,'FaceAlpha',0.3,'BarWidth',1);
    hold('on');
bar(Centers,h2,'EdgeAlpha',0,'FaceAlpha',0.3,'BarWidth',1);
plot(cutoff_NW_NS(n,1)*[1 1],get(gca,'ylim'),'k')
hold('off');
    end
end
%Compare to ref
% 
% mdltemp1 = compact(fitglm(BreathDataFullTable,['WakeNoAR ~ ARieF_pred_logit'],'Distribution','binomial','Link','logit','Exclude',BreathDataFullTable.ExcludeAR,'weights',weights))
% cutoffalldata = -mdltemp1.Coefficients.Estimate(1)/mdltemp1.Coefficients.Estimate(2)

%% No bias due to non-inclusion of EEG amplitude; but, was there before? would need to rerun analysis using "absolute power" approach

figure(100)
%cutoff_NW_NS(:,1)

temp = BreathDataFullTable.Subj;
temp2 = ['NaN';BreathDataFullTable.Subj(1:end-1)];
SubjStartRow=find(strcmp(temp,temp2)==0);
RefSubj = BreathDataFullTable.ref(SubjStartRow);

plot(RefSubj-mean(RefSubj),cutoff_NW_NS(:,1),'.')

%% Time course of Pwake after sleep onset

logitinverse = @(p) 1./(1+exp(-p));
normalci = @(p,n) 1.96*((p.*(1-p)./n).^0.5);

VarX_ = BreathDataFullTable.ARieF_pred_logit;
% VarX_ = logitinverse(BreathDataFullTable.ARieF_predID);
% VarX_ = BreathDataFullTable.ARieF_predID;
% VarX_ = BreathDataFullTable.ARieF_predModel;
% VarX_ = BreathDataFullTable.hypnog_B;

Exp = 1;
switch Exp
    case 1
    NbreathsInFuture=360;
    NbreathsInPast=10;
    %BreathDataFullTable.ARieF
    
    I=find([NaN;diff(BreathDataFullTable.ARieF)]<-0.5); %find front edge of arousal (could be in wake or sleep)
    %isolate analysis to those in or near wake at sleep onset
    if 1
        temp = 1*(sum(BreathDataFullTable.hypnog_B(I+[-NbreathsInPast:0])==4,2)==0); %remove if no recent wake
        I(temp==1)=[];
        
        temp = 1*(sum(BreathDataFullTable.ARieF(I+[-2:-1])<0.5,2)>0); %remove if any breaths Xback are in sleep
        I(temp==1)=[];
    end
    
    %BreathDataFullTable.ARieF(I(1:10)+[-3:-1])<0.5,@
    
    %remove breaths out of range
        I(find((I + NbreathsInFuture)>size(BreathDataFullTable,1)))=[];
    
    %remove if any arousal in range?
    if 1
        futureAR2 = 1*(sum(BreathDataFullTable.ARieF(I+[0:NbreathsInFuture])>0.5,2)>(NbreathsInFuture*0.1));
        I(futureAR2==1)=[];
    end
    
    case 2
        NbreathsInFuture=60;
        NbreathsInFutureRule=60;
        NbreathsInPast=30; %plot
        NbreathsInPastRule=30;

            I=find([NaN;diff(BreathDataFullTable.ARieF)]<-0.5); %find front edge of arousal (could be in wake or sleep)
        %isolate analysis to those in or near wake at sleep onset
        
        %remove breaths out of range
        I(find((I + NbreathsInFuture)>size(BreathDataFullTable,1)))=[];
    
        if 1
            temp = 1*(sum(BreathDataFullTable.hypnog_B(I+[-NbreathsInPastRule:0])<4,2)==0);
            I(temp==1)=[];
        end

        if 1
            temp = 1*(sum(BreathDataFullTable.ARieF(I+[-NbreathsInPastRule:-10])>0.5,2)>0);
            I(temp==1)=[];
        end

        %any arousal in range?
        if 1
        futureAR2 = 1*(sum(BreathDataFullTable.ARieF(I+[0:NbreathsInFutureRule])>0.5,2)>0);
        I(futureAR2==1)=[];
        end
end
    %I(find(BreathDataFullTable.hypnog_B(I)<4))=[];
   
%     I=find([NaN;diff(BreathDataFullTable.hypnog_B)]==-1 & BreathDataFullTable.hypnog_B==1);
%     I=find([NaN;diff(BreathDataFullTable.hypnog_B)]==-1);
    %remove data if stage is not equal to ?
    if 0
        temp = 1*(sum(BreathDataFullTable.hypnog_B(I+[0:NbreathsInFuture])~=1,2)>0);
        I(temp==1)=[];
    end
    


    %remove breaths if different subjects
    I((find(strcmp(BreathDataFullTable.Subj(I),BreathDataFullTable.Subj(I+NbreathsInFuture))==0)))=[];
    I((find(strcmp(BreathDataFullTable.Subj(I),BreathDataFullTable.Subj(I-NbreathsInPast))==0)))=[];
    
    %I(find(BreathDataFullTable.hypnog_B(I)==3)==1)=[];
    %I(find(BreathDataFullTable.hypnog_B(I)==4)==1)=[];
    
    
    %remove breaths if clearly beyond expected level
    I((find(BreathDataFullTable.Time_start(I+NbreathsInFuture)>BreathDataFullTable.Time_start(I)+NbreathsInFuture*6)))=[];
    I((find(BreathDataFullTable.Time_start(I-NbreathsInPast)<BreathDataFullTable.Time_start(I)-NbreathsInPast*6)))=[];
    
    
    
    
    medians = median(VarX_(I+[-NbreathsInPast:NbreathsInFuture]));
    mediansAR = mean(BreathDataFullTable.ARieF(I+[-NbreathsInPast:NbreathsInFuture]));
    modehyp = mode(BreathDataFullTable.hypnog_B(I+[-NbreathsInPast:NbreathsInFuture]));
    medianSWSpred = mode(BreathDataFullTable.SWSpred_logit(I+[-NbreathsInPast:NbreathsInFuture]));
    
    beta_t = median(BreathDataFullTable.Pbeta_ref(I+[-NbreathsInPast:NbreathsInFuture]));
    alpha_t = median(BreathDataFullTable.Palpha_ref(I+[-NbreathsInPast:NbreathsInFuture]));
    theta_t = median(BreathDataFullTable.Ptheta_ref(I+[-NbreathsInPast:NbreathsInFuture]));
    delta_t = median(BreathDataFullTable.Pdelta_ref(I+[-NbreathsInPast:NbreathsInFuture]));
    sigma_t = median(BreathDataFullTable.Psigma_ref(I+[-NbreathsInPast:NbreathsInFuture]));
    
    medianPredID = nanmedian(BreathDataFullTable.ARieF_predID(I+[-NbreathsInPast:NbreathsInFuture]));
    medianPredID_MY = nanmedian(BreathDataFullTable.ARieF_predID_MY(I+[-NbreathsInPast:NbreathsInFuture]));
    
    Time_ = mean(BreathDataFullTable.Time_start(I+[-NbreathsInPast:NbreathsInFuture]) - BreathDataFullTable.Time_start(I));
    
    meanVI = nanmean(BreathDataFullTable.VI(I+[-NbreathsInPast:NbreathsInFuture]));
    Nsamples = length(I)
    figure(17)
    plot(Time_,medians,'k-','markersize',12);
     hold('on')
    plot(Time_,mediansAR,'b-','markersize',12);
    %plot(Time_,logit(medianPredID),'r-','markersize',12);
    %plot(Time_,logit(medianPredID_MY),'r--','markersize',12);
    
    %plot(Time_,-medianSWSpred,'g-','markersize',12);
    
    plot(Time_,modehyp,'r-','markersize',12);
    %plot(Time_,meanVI*2-3,'k-','markersize',12);
    
%     plot(Time_,beta_t*3-5,'b-','markersize',12);
%     plot(Time_,alpha_t*3-5,'b-','markersize',12);
%     plot(Time_,theta_t*3-5,'b-','markersize',12);
%     plot(Time_,delta_t*3-5,'b-','markersize',12);
    %plot(Time_,sigma_t*3-5,'b:','markersize',12);
  
    %plot(logitinverse(xval),num./den,'.-') 
    
% hold('off')


%% Est SWS, within sleep
BreathDataFullTable.SWS = (BreathDataFullTable.hypnog_B==0)*1;
    BreathDataFullTable.SWS(isnan(BreathDataFullTable.hypnog_B))=NaN;
    
BreathDataFullTable.ExcludeForSWS = zeros(size(BreathDataFullTable,1),1);

%add in these lines to remove wake, or remove arousals
    %BreathDataFullTable.ExcludeForSWS(BreathDataFullTable.hypnog_B==4)=1;
    BreathDataFullTable.ExcludeForSWS(BreathDataFullTable.hypnog_B==4&BreathDataFullTable.ARieF<0.5)=1;
    BreathDataFullTable.ExcludeForSWS(BreathDataFullTable.hypnog_B~=4&BreathDataFullTable.ARieF>0.5)=1;
    BreathDataFullTable.ExcludeForSWS = logical(BreathDataFullTable.ExcludeForSWS);
    
balance = nanmean(BreathDataFullTable.SWS(~BreathDataFullTable.ExcludeForSWS))
weightsSWS = 0*weights;
weightsSWS(BreathDataFullTable.SWS>0.5)=1-balance;
weightsSWS(BreathDataFullTable.SWS<=0.5)=balance;
sum(weightsSWS(BreathDataFullTable.SWS>0.5&~BreathDataFullTable.ExcludeForSWS))
sum(weightsSWS(BreathDataFullTable.SWS<=0.5&~BreathDataFullTable.ExcludeForSWS))    
    
% mdl9 = compact(fitglm(BreathDataFullTable,'SWS ~ Ptheta_ref*Pdelta_ref + Ptheta_ref^2 + Pbeta_ref*Palpha_ref + Pdelta_ref + Pdelta_ref^2','Distribution','binomial','Link','logit','Exclude',BreathDataFullTable.ExcludeForSWS,'weights',weightsSWS))
% mdl9.Rsquared.Ordinary

mdl9 = compact(fitglm(BreathDataFullTable,'SWS ~ ARieF_pred_logit','Distribution','binomial','Link','logit','Exclude',BreathDataFullTable.ExcludeForSWS,'weights',weightsSWS))
mdl9.Rsquared.Ordinary
cutoffSWS_ = -mdl9.Coefficients.Estimate(1)/mdl9.Coefficients.Estimate(2)
% 
% mdl9 = compact(fitglm(BreathDataFullTable,'SWS ~ ARieF_pred_logit + Pbeta_ref + Pdelta_ref','Distribution','binomial','Link','logit','Exclude',BreathDataFullTable.ExcludeForSWS,'weights',weightsSWS))
% mdl9.Rsquared.Ordinary
    
BreathDataFullTable.SWSpred_logit = logit(predict(mdl9,BreathDataFullTable));

%% WeightsWSC and Linear regression version to develop SWS->W continuum "WSC"

%StageMedians_ = [ 4.3469   -2.6775   -3.1186   -4.2288   -5.7585]    W, N1, REM, N2, N3
    
BreathDataFullTable.WSC = NaN*zeros(size(BreathDataFullTable,1),1);
    BreathDataFullTable.WSC(BreathDataFullTable.hypnog_B==4) = 4.3469;
    BreathDataFullTable.WSC(BreathDataFullTable.hypnog_B==3) = -3.1186;
    BreathDataFullTable.WSC(BreathDataFullTable.hypnog_B==2) = -2.6775;
    BreathDataFullTable.WSC(BreathDataFullTable.hypnog_B==1) = -4.2288;
    BreathDataFullTable.WSC(BreathDataFullTable.hypnog_B==0) = -5.7585;
    BreathDataFullTable.WSC(isnan(BreathDataFullTable.hypnog_B))=NaN;
    

BreathDataFullTable.ExcludeAR = zeros(size(BreathDataFullTable,1),1);
    BreathDataFullTable.ExcludeAR(BreathDataFullTable.hypnog_B==4&BreathDataFullTable.ARieF<0.5)=1;
    BreathDataFullTable.ExcludeAR(BreathDataFullTable.hypnog_B~=4&BreathDataFullTable.ARieF>0.5)=1;
    BreathDataFullTable.ExcludeAR = logical(BreathDataFullTable.ExcludeAR);
    
NperStage(1) = sum(BreathDataFullTable.hypnog_B==4&~BreathDataFullTable.ExcludeAR);    
NperStage(2) = sum(BreathDataFullTable.hypnog_B==3&~BreathDataFullTable.ExcludeAR)    
NperStage(3) = sum(BreathDataFullTable.hypnog_B==2&~BreathDataFullTable.ExcludeAR)    
NperStage(4) = sum(BreathDataFullTable.hypnog_B==1&~BreathDataFullTable.ExcludeAR)    
NperStage(5) = sum(BreathDataFullTable.hypnog_B==0&~BreathDataFullTable.ExcludeAR)    
Ntotal = sum(~isnan(BreathDataFullTable.hypnog_B)&~BreathDataFullTable.ExcludeAR)
sum(NperStage)

NperStageF = NperStage/sum(NperStage)

weightsWSC = zeros(size(BreathDataFullTable,1),1);
    weightsWSC(BreathDataFullTable.hypnog_B==4) = 1/NperStageF(1);
    weightsWSC(BreathDataFullTable.hypnog_B==3) = 1/NperStageF(2);
    weightsWSC(BreathDataFullTable.hypnog_B==2) = 1/NperStageF(3);
    weightsWSC(BreathDataFullTable.hypnog_B==1) = 2/NperStageF(4);
    weightsWSC(BreathDataFullTable.hypnog_B==0) = 4/NperStageF(5);
    
sum(weightsWSC(BreathDataFullTable.hypnog_B>=0&BreathDataFullTable.hypnog_B<4))
sum(weightsWSC(BreathDataFullTable.hypnog_B==4))

weightsWSC(BreathDataFullTable.hypnog_B>=0&BreathDataFullTable.hypnog_B<4) = weightsWSC(BreathDataFullTable.hypnog_B>=0&BreathDataFullTable.hypnog_B<4)/nansum(weightsWSC(BreathDataFullTable.hypnog_B>=0&BreathDataFullTable.hypnog_B<4));
weightsWSC(BreathDataFullTable.hypnog_B==4) = weightsWSC(BreathDataFullTable.hypnog_B==4)/nansum(weightsWSC(BreathDataFullTable.hypnog_B==4));

   
sum(weightsWSC(BreathDataFullTable.hypnog_B>=0&BreathDataFullTable.hypnog_B<4))
sum(weightsWSC(BreathDataFullTable.hypnog_B==4))

mdl10 = compact(fitglm(BreathDataFullTable,'WSC ~ Ptheta_ref*Pdelta_ref + Ptheta_ref^2 + Pbeta_ref + Palpha_ref + Pdelta_ref + Pdelta_ref^2','Exclude',BreathDataFullTable.ExcludeAR,'weights',weightsWSC))
mdl10.Rsquared.Ordinary
BreathDataFullTable.WSCpred = predict(mdl10,BreathDataFullTable);
thres = 0.5;
performanceWSC = PredictiveValue(1*(BreathDataFullTable.WakeNoAR(~BreathDataFullTable.ExcludeAR)>thres),1*(BreathDataFullTable.WSCpred(~BreathDataFullTable.ExcludeAR)>thres),BreathDataFullTable.WakeNoAR(~BreathDataFullTable.ExcludeAR))

%use as AR
BreathDataFullTable.WSCpred_last1 = [NaN;BreathDataFullTable.WSCpred(1:end-1)];
BreathDataFullTable.WSCpred_last2 = [NaN;NaN;BreathDataFullTable.WSCpred(1:end-2)];
BreathDataFullTable.WSCpred_last3 = [NaN;NaN;NaN;BreathDataFullTable.WSCpred(1:end-3)];
BreathDataFullTable.WSCpred_last4 = [NaN;NaN;NaN;NaN;BreathDataFullTable.WSCpred(1:end-4)];
BreathDataFullTable.WSCpred_last5 = [NaN;NaN;NaN;NaN;NaN;BreathDataFullTable.WSCpred(1:end-5)];
BreathDataFullTable.WSCpred_last6 = [NaN;NaN;NaN;NaN;NaN;NaN;BreathDataFullTable.WSCpred(1:end-6)];
BreathDataFullTable.WSCpred_last6_min = min([BreathDataFullTable.WSCpred_last1';BreathDataFullTable.WSCpred_last2';BreathDataFullTable.WSCpred_last3';BreathDataFullTable.WSCpred_last4';BreathDataFullTable.WSCpred_last5';BreathDataFullTable.WSCpred_last6'])';

mdlAR = compact(fitglm(BreathDataFullTable,'ARieF ~ WSCpred + WSCpred_last6_min + WSCpred_last3 + WSCpred_last4','Exclude',ExcludeW,'Distribution','binomial','Link','logit','weights',weightsAR))
mdlAR.Rsquared.Ordinary

newARpred_ = predict(mdlAR,BreathDataFullTable);
BreathDataFullTable.newARpredF = newARpred_;
BreathDataFullTable.newARpred_logit = logit(newARpred_);
BreathDataFullTable.newARpred = 1*(newARpred_>0.5);
    BreathDataFullTable.newARpred(isnan(newARpred_)) = NaN;
thres=0.5;
performanceAR = PredictiveValue(1*(BreathDataFullTable.ARieF(~ExcludeW)>0.5),1*(newARpred_(~ExcludeW)>thres),BreathDataFullTable.ARieF(~ExcludeW))
%performanceX = PredictiveValue(1*(BreathDataFullTable.ARieF(~ExcludeW)>0.5),1*(BreathDataFullTable.ARieF_pred_logit(~ExcludeW)>-1.2954),BreathDataFullTable.ARieF(~ExcludeW))

figure(121); clf(121);
dStep=0.1;
Centers=-10:dStep:10;
Edges=(Centers(1)-dStep/2):dStep:(Centers(end)+dStep/2);
subplot(1,1,1);
[h7,edges] = histcounts(logit(newARpred_(~ExcludeW&BreathDataFullTable.ARieF>0.5)),Edges);
[h8,edges] = histcounts(logit(newARpred_(~ExcludeW&BreathDataFullTable.ARieF<0.5)),Edges);
bar(Centers,h7,'EdgeAlpha',0,'FaceAlpha',0.3,'BarWidth',1);
hold('on');
bar(Centers,h8,'EdgeAlpha',0,'FaceAlpha',0.7,'BarWidth',1);




%% Histograms

figure(12); clf(12);
% XvalH = BreathDataFullTable.WSCpred;
% 
% figure(11); clf(11);

% XvalH = logit(BreathDataFullTable.ARieF_predModel); 
 
 XvalH = BreathDataFullTable.WAS; 
 ArTemp = BreathDataFullTable.ARieF;
 XvalH = BreathDataFullTable.ARieF_pred_logit; 

logit = @(p) log(p./(1-p));

dStep=0.1;
Centers=-10:dStep:10;
Edges=(Centers(1)-dStep/2):dStep:(Centers(end)+dStep/2);

[h1,edges] = histcounts((XvalH(BreathDataFullTable.hypnog_B==4&ArTemp>0.95)),Edges); 
[h2,edges] = histcounts((XvalH(BreathDataFullTable.hypnog_B<4&ArTemp>0.95)),Edges); 
[h3,edges] = histcounts((XvalH(BreathDataFullTable.hypnog_B==0&ArTemp<0.05)),Edges); 
[h4,edges] = histcounts((XvalH(BreathDataFullTable.hypnog_B==1&ArTemp<0.05)),Edges); 
[h5,edges] = histcounts((XvalH(BreathDataFullTable.hypnog_B==2&ArTemp<0.05)),Edges); 
[h6,edges] = histcounts((XvalH(BreathDataFullTable.hypnog_B==3&ArTemp<0.05)),Edges);
%[h7,edges] = histcounts(logit(newARpred(BreathDataFullTable.hypnog_B<4&ArTemp>0.95)),Edges);
[h8,edges] = histcounts((XvalH(BreathDataFullTable.hypnog_B==4&ArTemp<0.05)),Edges); 
    

medianPrW = nanmedian(XvalH(BreathDataFullTable.hypnog_B==4&~BreathDataFullTable.ExcludeAR)); 
    meanPrW = nanmean(XvalH(BreathDataFullTable.hypnog_B==4&~BreathDataFullTable.ExcludeAR)); 
medianPrAR = nanmedian(XvalH(BreathDataFullTable.hypnog_B<4&~BreathDataFullTable.ExcludeAR)); 
    meanPrAR = nanmean(XvalH(BreathDataFullTable.hypnog_B<4&~BreathDataFullTable.ExcludeAR)); 
medianPrSOW = nanmedian(XvalH(BreathDataFullTable.hypnog_B==4&~BreathDataFullTable.ExcludeAR)); 
    meanPrSOW = nanmean(XvalH(BreathDataFullTable.hypnog_B==4&~BreathDataFullTable.ExcludeAR)); 
medianPrN3 = nanmedian(XvalH(BreathDataFullTable.hypnog_B==0&~BreathDataFullTable.ExcludeAR)); 
medianPrN2 = nanmedian(XvalH(BreathDataFullTable.hypnog_B==1&~BreathDataFullTable.ExcludeAR)); 
medianPrN1 = nanmedian(XvalH(BreathDataFullTable.hypnog_B==2&~BreathDataFullTable.ExcludeAR)); 
medianPrREM = nanmedian(XvalH(BreathDataFullTable.hypnog_B==3&~BreathDataFullTable.ExcludeAR));

meanPrW = nanmean(XvalH(BreathDataFullTable.hypnog_B==4&~BreathDataFullTable.ExcludeAR)); 
meanPrR = nanmean(XvalH(BreathDataFullTable.hypnog_B==3&~BreathDataFullTable.ExcludeAR)); 
meanPr1 = nanmean(XvalH(BreathDataFullTable.hypnog_B==2&~BreathDataFullTable.ExcludeAR)); 
meanPr2 = nanmean(XvalH(BreathDataFullTable.hypnog_B==1&~BreathDataFullTable.ExcludeAR)); 
meanPr3 = nanmean(XvalH(BreathDataFullTable.hypnog_B==0&~BreathDataFullTable.ExcludeAR)); 

StageMeans2 = [meanPrW meanPr1 meanPrR meanPr2 meanPr3]
logitinverse(StageMeans2)
StageMedians2 = [medianPrW medianPrN1 medianPrREM medianPrN2 medianPrN3]
logitinverse(StageMedians2)

%bar(Centers,h1+h2,'EdgeAlpha',0,'BarWidth',1);

if 0
bar(Centers,h1/sum(h1),'EdgeAlpha',0,'FaceAlpha',0.3,'BarWidth',1);
hold('on');
bar(Centers,h2/sum(h2),'EdgeAlpha',0,'FaceAlpha',0.7,'BarWidth',1);
bar(Centers,h4/sum(h4),'EdgeAlpha',0,'FaceAlpha',0.3,'BarWidth',1);
bar(Centers,h5/sum(h5),'EdgeAlpha',0,'FaceAlpha',0.3,'BarWidth',1);
bar(Centers,h6/sum(h6),'EdgeAlpha',0,'FaceAlpha',0.3,'BarWidth',1);
bar(Centers,h3/sum(h3),'EdgeAlpha',0,'FaceAlpha',0.3,'BarWidth',1);
else
    bar(Centers,h1,'EdgeAlpha',0,'FaceAlpha',0.3,'BarWidth',1);
hold('on');
%stairs(Centers-dStep/2,h1);
stairs(Centers-dStep/2,h2);
bar(Centers,h4,'EdgeAlpha',0,'FaceAlpha',0.3,'BarWidth',1);
bar(Centers,h5,'EdgeAlpha',0,'FaceAlpha',0.3,'BarWidth',1);
bar(Centers,h6,'EdgeAlpha',0,'FaceAlpha',0.3,'BarWidth',1);
bar(Centers,h3,'EdgeAlpha',0,'FaceAlpha',0.3,'BarWidth',1);
%stairs(Centers,h8*10);
%bar(Centers,hall,'EdgeAlpha',0,'FaceAlpha',0.3,'BarWidth',1);
end
%bar(Centers,h7/sum(h7),'EdgeAlpha',0,'FaceAlpha',0.3,'BarWidth',1);
legend('W','AR','N2','N1','REM','N3')
