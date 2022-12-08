function [WhyT,DataTableT] = getCSAvsOSAall(MrangeOverride,reruns)
% MrangeOverride=[1 3:37 39:72]; -- for Oxy study

global settings AMasterSpreadsheet

settings = ImportSettings(settings,AMasterSpreadsheet);

if ~exist('reruns')
    reruns=0;
end

try
    if reruns
        error();
    end
    load([settings.workdir '\Summary\SummaryAnalysis.mat'],'SummaryAnalysisTable');
    
    %adding handling for whether there is data for all "MrangeOverride" patients
    % elements of MrangeOverride should all be present in Summary analysis
    % table--done using ismember function
    
    %add Nan checking
    if height(SummaryAnalysisTable)<size(MrangeOverride,2) || sum(ismember(MrangeOverride,[1:height(SummaryAnalysisTable)]))<size(MrangeOverride,2)
        SummaryAnalysisTable = SummaryAnalysisN(MrangeOverride);
    end
catch 
    SummaryAnalysisTable = SummaryAnalysisN(MrangeOverride); %make this code have an option to avoid saving data (output table only)
end

try
    if reruns
        error();
    end
    load([settings.workdir '\Summary\getFL.mat'],'FLFinalT');
    
    %adding handling for whether there is data for all "MrangeOverride" patients
    if height(FLFinalT)<size(MrangeOverride,2) || sum(ismember(MrangeOverride,[1:height(FLFinalT)]))<size(MrangeOverride,2)
       [FLFinalT,~] = getFL(MrangeOverride,1);
    end

catch 
    [FLFinalT,~] = getFL(MrangeOverride,1);
end

if 0
    EventAnalysis(MrangeOverride);
end

%% IMPORTING EAinfo

% need to have Events Analysis performed prior to this.
tic
if ~exist('EventInfoT')
    clear EventInfoT EventInfoT1
    for i=MrangeOverride
        try
            clear EAinfo EventInfoT1
            load([settings.workdir '\Analyzed\EventAnalyzed\' settings.savename '_' num2str(i) '.mat'],'EAinfo');
            EventInfoT1 = EAinfo.OvsC;
            temp = EAinfo.EvtFtrs;
            temp = rmfield(temp,{'EventDepth','DesatSlope','HarmonicPeaks'});
            EventInfoT1 = mergestructs(EventInfoT1,temp);
            EventInfoT1 = struct2table(EventInfoT1);
            i
            if i==1
                EventInfoT = EventInfoT1;
            else
                EventInfoT(i,:) = EventInfoT1
            end
        catch
%             if exist('EventInfoT')
                EventInfoT{i,:} = nan(1,width(EventInfoT));
%             end
        end
%         try
%             size(EventInfoT) %feedback to user
%         end
    end
end
toc


%% Combine

%DataTable = [SummaryAnalysisTable(MrangeOverride,:) EventInfoT(MrangeOverride,:) FLFinalT];
tempT = [SummaryAnalysisTable EventInfoT FLFinalT];
DataTableT = tempT;
DataTableT{:,:}=NaN;
DataTableT{MrangeOverride,:}=tempT{MrangeOverride,:};

DataTable = [SummaryAnalysisTable(MrangeOverride,:) EventInfoT(MrangeOverride,:) FLFinalT(MrangeOverride,:)];

try
DataTableT.VRAT = ((DataTableT.VRA/100).^0.5)*100;
catch
DataTableT.VRAT = ((DataTableT.VRA1/100).^0.5)*100;    
end

%%
load CSAvsOSA mdlCSA
load CSAvsOSA DataTable2
%%
logitinverse = @(x) 1./(1+exp(-x));
logit = @(p) log(p./(1-p));
mdl=mdlCSA.m3;
DataTableT.CSAlogodds3 = logit(predict(mdl,DataTableT));
mdl=mdlCSA.m4;
DataTableT.CSAlogodds4 = logit(predict(mdl,DataTableT));
mdl=mdlCSA.m5;
DataTableT.CSAlogodds5 = logit(predict(mdl,DataTableT));
mdl=mdlCSA.m6;
DataTableT.CSAlogodds6 = logit(predict(mdl,DataTableT));
mdl=mdlCSA.mtraits;
DataTableT.CSAlogoddstraits = logit(predict(mdl,DataTableT));

% mdl5 is:
% % 
% mdlX1=fitglm(DataTableT,'P ~ LGn  + delay +FlowDriveMedian +FlowDriveNadirV +FlowDriveUnder50p + EvtCov   + VRAT','distribution','binomial')
% mdlX1.Rsquared.Ordinary.^0.5
%fitglm(DataTable2,'CSA ~ EvtCov','distribution','binomial')
%fitglm(DataTable2,'LGn ~ EvtCov')

%%
DataTable = DataTableT(MrangeOverride(end),:);

%getting details for last subject in the list

mdl=mdlCSA.m5;

Labels = mdl.Coefficients.Properties.RowNames
%Labels = Labels(2:end)
Labels = [Labels ; 'Total'];
clear Value GroupMean GroupSD
for i=1:length(Labels)
    
    if i==1
        Beta(i,:) = mdl.Coefficients.Estimate(i);
        Value(i,:) = 1;
        GroupMean(i,:) = 1;%1-(3.3266)/7.207998653673618;
        GroupSD(i,:) = 0;
    elseif i<length(Labels)
        Beta(i,:) = mdl.Coefficients.Estimate(i);
%         Value(i,:) = DataTable{:,Labels{i}}'; % if using single patient
        Value(i,:) = nanmean(DataTable{:,Labels{i}}'); % if using full dataset
        GroupMean(i,:) = nanmean(DataTable2{:,Labels{i}});
        GroupSD(i,:) = nanstd(DataTable2{:,Labels{i}});
    else
        Beta(i,:) = 1;
        Value(i,:) = logit(predict(mdl,DataTable));
        GroupMean(i,:) = nanmean(logit(predict(mdl,DataTable2)));
        GroupSD(i,:) = nanstd(logit(predict(mdl,DataTable2)));
    end
end
if height(DataTable)==1
    WhyT = table(Labels,Value,GroupMean,GroupSD,Beta);
end

WhyT.Diff = (WhyT.Value - WhyT.GroupMean);
WhyT.PartialLogOdds = (WhyT.Value - WhyT.GroupMean).*WhyT.Beta;
temp = GroupMean.*Beta;
WhyT.PartialLogOdds(1) = sum(temp(1:6));
WhyT.PartialLogOdds(end) = sum(WhyT.PartialLogOdds(1:end-1));

writetable(WhyT,[settings.AMasterdir 'CSAvOSAT.xlsx']);


%%
return

%% for oxy study
DataTable.Ypred = logit(predict(mdlCSA.m5,DataTable));
BaselineT=DataTable(1:35,:); % Oxy: using first 36 baseline studies; discarding 2nd subject due to poor FlowQ-
RespNonResp=[ones(9,1);zeros(26,1)]; % first 9 responders
DataTable.RespNonResp=RespNonResp;
varlist={'Ypred'};
clear splithres AUC
for jj=1:length(varlist)
    h= figure(jj);
    set(gcf,'color',[1 1 1]);
    
    Data = DataTable{:,varlist{jj}};
    crit = DataTable2.RespNonResp==1;
    crit2 = DataTable2.RespNonResp==0;
    PlotThreeGroupIndividualData(h,Data,crit,crit2)
    Y = Data;
    labels = crit;
    I = isnan(Y) | (crit==0 & crit2==0);
    Y(I)=[]; labels(I)=[];
    
    [splithres(jj),AUC(jj),~,~,~,~]=ROCAUCSEM(labels,Y);
    %
    ylabel(varlist{jj})
    set(gca,'xtick',[1 2],'xticklabels',{'Responders','Non-responders'},'fontsize',12);
    
    xtickangle(45);
    
end
AUC(jj)

