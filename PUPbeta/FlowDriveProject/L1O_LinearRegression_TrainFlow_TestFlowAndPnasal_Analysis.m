%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
close; clear; clc
mypath = 'C:\Users\uqdmann\Dropbox\QAO\FlowAndEdiPesQAO_CodeForPublication';
addpath(mypath);
cd(mypath);
datadir = 'C:\PSG_Data\FlowDrive\Analyzed\Publication\'; % data
experimentnumber = '__n000'; %
load([datadir, 'FlowDrive_only25Hz_ReNormalized_FeatureSpace_AutoRef2_Edi_Clean_LinRegWorkspace_TrainFlow_TestFlowAndPnasal', experimentnumber, '.mat']);
% probably should just load variables of interest to save time/space/electrons
% 'PtData','Amatrix2','FeatureNames','etc...'

if 0 % this is the final run
load('C:\PSG_Data\FlowDrive\Analyzed\ExclBadR\FlowDrive_only25Hz_ReNormalized_FeatureSpace_AutoRef2_Edi_Clean_FinalRun2.mat');
end

%% Tidy up
clearvars -except ...
'Amatrix2_flow' ...
'Amatrix2_pnasal' ...
'Amatrix_flow' ...
'Amatrix_pnasal' ...
'Bad_Rsq_vals_i' ...
'beta_array_flow' ...
'ErrTrain_array_flow' ...
'FeatureNames' ...
'GS' ...
'Ind' ...
'Isubj' ...
'Labels' ...
'labels_flow' ...
'labels_Step_Subj_flow' ...
'LabelsOrderedOpt' ...
'maxDataSize' ...
'maxG' ...
'numofftrs' ...
'Original_FeatureNames' ...
'predyL1O_array_flow' ...
'predyL1O_array_pnasal' ...
'PtData_flow' ...
'PtData_pnasal' ...
'RemovedFtrs' ...
'RsqTrain_array_flow' ...
'settings' ...
'weights' 

%% Use All breaths for classification performance plots
SleepOnly = 0;  % 0 for All breaths, 1 for sleep only breaths
IncludePnasal = 1;
[BB, BB_] = getAllorSleepOnlyBB(SleepOnly, IncludePnasal, PtData_flow, PtData_pnasal);
str = ['Flow - Using ',num2str(nnz(BB)),' of ', num2str(nnz(~isnan(PtData_flow.PT)))]; disp(str);
str = ['Pnasal - Using ',num2str(nnz(BB_)),' of ', num2str(nnz(~isnan(PtData_pnasal.PT)))]; disp(str);

% Select data to use
PlotUptoNFtrs = min([50 maxDataSize]); % normally 100

ftrnum=min([25 maxDataSize]); % normally 25
ErrorMaxN = size(predyL1O_array_flow,2); % normally 100

% pnasal
predy_pnasal = predyL1O_array_pnasal(BB_,ftrnum);
Yval_pnasal = GS(BB_);
weights_pnasal = weights(BB_);

% flow
predy_flow = predyL1O_array_flow(BB,ftrnum);
Yval_flow = GS(BB);
weights_flow = weights;

savefigas = ''; %   'saveasPNG'; % ''; % 'saveasTIFF'; % options are saveasPNG, saveasFIG, and saveasTIFF
closefigs = 0;

TickFntSz = 12;
LabelFntSz = 18;
FntSz = 18;
sigdig = '%.2f';

%% Process test data. flow and pnasal, for all pts combined
for i=1:size(predyL1O_array_flow,2)% normally size(predyL1O_array,2) PlotUptoNFtrs 
    [RsqL1O_pnasal(i), RL1O_pnasal(i)] = ...
    UnivariateStats(predyL1O_array_pnasal(BB_,i),GS(BB_), weights(BB_)); % new method for R and Rsq
    ErrL1O_pnasal(i) = nanmean(weights(BB_).*abs(predyL1O_array_pnasal((BB_),i)-GS(BB_)));
    [RsqL1O_flow(i), RL1O_flow(i)] = ...
        UnivariateStats(predyL1O_array_flow(BB,i),GS(BB), weights(BB)); % new method for R and Rsq
    ErrL1O_flow(i) = nanmean(weights(BB).*abs(predyL1O_array_flow((BB),i)-GS(BB)));
end

%% For the average of all pts combined
figure(1); clf(figure(1));fig = gcf;
fig.Color = [1 1 1]; fig.Units = 'inches';
fig.Position = [-19  -0.5   16   5];
subplot(1,3,1);
for pt = 1:54
    if ismember(pt, settings.Flow_list)
        %stairs(RsqTrain_array_flow(pt,:).^0.5, 'g'); % is wrong, only testing        
        stairs(RsqTrain_array_flow(pt,:), 'r'); hold on;
        stairs(ErrTrain_array_flow(pt,:), 'k');
    end
end
xlim([0, PlotUptoNFtrs]); ylim([0,1]);
title('Flow - Training'); legend('R squared', 'MAE'); axis square

subplot(1,3,2);
% stairs(RsqL1O_flow.^0.5,'g'); % is wrong, only testing
stairs(RL1O_flow,'b'); hold on;
stairs(RsqL1O_flow,'r');
stairs(ErrL1O_flow, 'k');
xlim([0, PlotUptoNFtrs]); ylim([0,1]);
title('Flow - Testing'); legend('R','R squared','MAE'); axis square

subplot(1,3,3);
%stairs(RL1O_pnasal.^0.5,'g'); % is wrong, only testing
stairs(RL1O_pnasal,'b'); hold on; 
stairs(RsqL1O_pnasal,'r'); 
stairs(ErrL1O_pnasal, 'k');
xlim([0, PlotUptoNFtrs]); ylim([0,1]);
title('Pnasal - Testing'); legend('R','R squared','MAE'); axis square

%% Table E4 - individual patient performance in classification
% i.e. exactly correct, and correct within one severity class
classcutoffs = [0.9 0.7 0.5 0.3];
Nclasses=length(classcutoffs)+1;
ftrs = 25;
ExactClass_perPT_flow=[];
WithinOneClass_perPT_flow=[];
MAE_perPT=[];
for subj=1:54 %length(PT_list) % subj=44
    if ~ismember(subj, settings.Flow_list); continue; end
    Isubj = PtData_flow.PT==subj & BB;
    Yval = GS(Isubj);
    predy = predyL1O_array_flow(Isubj,ftrs);
  
    g = NaN*Yval;
    ghat = NaN*predy;
    g(Yval>classcutoffs(1))=1;
    ghat(predy>classcutoffs(1))=1;
    for i=2:Nclasses
        g(Yval<classcutoffs(i-1))=i;
        ghat(predy<classcutoffs(i-1))=i;
    end
    [C,order] = confusionmat(g,ghat);
    %sumactual=sum(C')';
    sumactual=nansum(C,2);
    sumestimated=nansum(C);
    %C_Factual=C./sum(C')'
    C_Factual=C./nansum(C,2)*100; %rows are actual, cols are estimated
    C_Festimated=C./nansum(C)*100;
    C_Total = (C./nansum(C(:)))*100;
    
    AccA_C_Factual = C_Factual.*diag(ones(1,length(C)));
    AccB_C_Factual = C_Factual.*(diag(ones(1,length(C)-1),1)+diag(ones(1,length(C)-1),-1));
    AccAactual = mean(sum(AccA_C_Factual,2));
    AccBactual = mean(sum(AccA_C_Factual,2)) + mean(sum(AccB_C_Factual,2));
    AccA_C_Festimated = C_Festimated.*diag(ones(1,length(C)));
    AccB_C_Festimated = C_Festimated.*(diag(ones(1,length(C)-1),1)+diag(ones(1,length(C)-1),-1));
    AccAestimated = nanmean(sum(AccA_C_Festimated));
    AccBestimated = nanmean(sum(AccA_C_Festimated)) + nanmean(sum(AccB_C_Festimated));
    ACCs = [AccAactual AccBactual; AccAestimated AccBestimated];
    %second col gives accuracy if accepting next-category error as correct
    
    ExactClass_perPT_flow(subj) = AccAestimated;
    WithinOneClass_perPT_flow(subj) = AccBestimated;
    
end

% PerPtPerformance_flow = table([1:41]',...
%     round(ExactClass_perPT_flow(settings.Flow_list),1)',...
%     round(WithinOneClass_perPT_flow(settings.Flow_list),1)',... %round(MAE_perPT(PT_list),2)', ...
%     'VariableNames', {'PT','Exact','WithinOne'});
PerPtPerformance_flow = table([1:54]',...
    round(ExactClass_perPT_flow,1)',...
    round(WithinOneClass_perPT_flow,1)',... %round(MAE_perPT(PT_list),2)', ...
    'VariableNames', {'PT','Exact','WithinOne'});

ExactClass_perPT_pnasal=[];
WithinOneClass_perPT_pnasal=[];
for subj=1:54 %length(PT_list) % subj=44
    if ~ismember(subj, settings.Pnasal_list); continue; end
    Isubj = PtData_pnasal.PT==subj & BB;
    Yval = GS(Isubj);
    predy = predyL1O_array_pnasal(Isubj,ftrs);
  
    g = NaN*Yval;
    ghat = NaN*predy;
    g(Yval>classcutoffs(1))=1;
    ghat(predy>classcutoffs(1))=1;
    for i=2:Nclasses
        g(Yval<classcutoffs(i-1))=i;
        ghat(predy<classcutoffs(i-1))=i;
    end
    [C,order] = confusionmat(g,ghat);
    %sumactual=sum(C')';
    sumactual=nansum(C,2);
    sumestimated=nansum(C);
    %C_Factual=C./sum(C')'
    C_Factual=C./nansum(C,2)*100; %rows are actual, cols are estimated
    C_Festimated=C./nansum(C)*100;
    C_Total = (C./nansum(C(:)))*100;
    
    AccA_C_Factual = C_Factual.*diag(ones(1,length(C)));
    AccB_C_Factual = C_Factual.*(diag(ones(1,length(C)-1),1)+diag(ones(1,length(C)-1),-1));
    AccAactual = mean(sum(AccA_C_Factual,2));
    AccBactual = mean(sum(AccA_C_Factual,2)) + mean(sum(AccB_C_Factual,2));
    AccA_C_Festimated = C_Festimated.*diag(ones(1,length(C)));
    AccB_C_Festimated = C_Festimated.*(diag(ones(1,length(C)-1),1)+diag(ones(1,length(C)-1),-1));
    AccAestimated = nanmean(sum(AccA_C_Festimated));
    AccBestimated = nanmean(sum(AccA_C_Festimated)) + nanmean(sum(AccB_C_Festimated));
    ACCs = [AccAactual AccBactual; AccAestimated AccBestimated];
    %second col gives accuracy if accepting next-category error as correct
    
    ExactClass_perPT_pnasal(subj) = AccAestimated;
    WithinOneClass_perPT_pnasal(subj) = AccBestimated;
    
end

% PerPtPerformance_pnasal = table([1:17]',...
%     round(ExactClass_perPT_pnasal(settings.Pnasal_list),1)',...
%     round(WithinOneClass_perPT_pnasal(settings.Pnasal_list),1)',... %round(MAE_perPT(PT_list),2)', ...
%     'VariableNames', {'PT','Exact','WithinOne'})
PerPtPerformance_pnasal = table([1:54]',...
    round(ExactClass_perPT_pnasal,1)',...
    round(WithinOneClass_perPT_pnasal,1)',... %round(MAE_perPT(PT_list),2)', ...
    'VariableNames', {'PT','Exact','WithinOne'});

PerPtPerformance = outerjoin(PerPtPerformance_flow, PerPtPerformance_pnasal, 'key', 'PT', 'mergekeys', 1);

exclude = PerPtPerformance.Exact_PerPtPerformance_flow==0;
PerPtPerformance(exclude,:)=[];
PerPtPerformance.PT = [1:1:41]';
Table_E4 = PerPtPerformance;

%% Figure E3 - flow
figure(83); clf(figure(83)); fig = gcf;
fig.Color = [1 1 1]; fig.Units = 'inches';
fig.Position = [22 1.5   16   5];

% performance vs number of features
subplot(1,3,1);
%stairs([RL1O_flow(1:PlotUptoNFtrs);RsqL1O_flow(1:PlotUptoNFtrs);ErrL1O_flow(1:PlotUptoNFtrs)]'); hold on;
stairs([1-RsqL1O_flow(1:PlotUptoNFtrs);ErrL1O_flow(1:PlotUptoNFtrs)]'); hold on;
%currentYlim=ylim();
%plot([ftrnum, ftrnum], [0, currentYlim(2)],'k-', 'linewidth', 1.5);
plot([ftrnum, ftrnum], [0, 0.5],'k-', 'linewidth', 1.5);
ax = gca; set(ax,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz);
xlabel('Number of Features', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
ylabel('Performance (Testing)', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
ylim([0 0.7]); xlim([0 50]); axis square
if 1
    %text(30, 0.7, {'correlation coefficient'}, 'fontname','arial narrow','FontSize',10, 'FontWeight', 'normal');
    text(30, 0.45, {'1 - coefficient of determination'}, 'fontname','arial narrow','FontSize',10, 'FontWeight', 'normal');
    text(30, 0.23, {'mean absolute error'}, 'fontname','arial narrow','FontSize',10, 'FontWeight', 'normal');
end

% proportion of total error
subplot(1,3,2);
a_out = mean(ErrTrain_array_flow(settings.Flow_list,1:ErrorMaxN),1);
a_diff = diff(a_out);
err100 = a_out(ErrorMaxN);
p_of100 = 100.*(err100 ./ a_out(1:ErrorMaxN));
plot(100-p_of100, 'k-'); hold on;
currentYlim=ylim();
plot([ftrnum, ftrnum], [0, currentYlim(2)/2],'k-', 'linewidth', 1.5);
ax = gca; set(ax,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz);
xlabel('Number of Features', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
ylabel('Proportion of absolute error (%)', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
ylim([0 30]); xlim([0 50]); axis square

% confusion mat
subplot(1,3,3);
customcmap = GetCustomColorMap('gray'); % 'SS'
g = NaN*Yval_flow;
ghat = NaN*predy_flow;
classes = 5;
switch classes
    case 5
        classcutoffs = [0.9 0.7 0.5 0.3];
        labeltext={'Normal','Mild','Moderate','Severe','V.Severe'};
    case 4
        classcutoffs = [0.9 0.5 0.3];
        labeltext={'Normal','Mild','Moderate','Severe'};
    case 3
        classcutoffs = [0.9 0.5];
        labeltext={'Normal','Intermediate','Flow Limited'};
    otherwise
        classcutoffs = [0.7]; %normal mild moderate severe v.severe / normal borderline mild moderate severe
        labeltext={'Normal','Flow Limited'};
end
Nclasses=length(classcutoffs)+1;
g(Yval_flow>classcutoffs(1))=1;
ghat(predy_flow>classcutoffs(1))=1;
for i=2:Nclasses
    g(Yval_flow<classcutoffs(i-1))=i;
    ghat(predy_flow<classcutoffs(i-1))=i;
end
[C,order] = confusionmat(g,ghat);
%sumactual=sum(C')';
sumactual=sum(C,2);
sumestimated=sum(C);
%C_Factual=C./sum(C')'
C_Factual=C./sum(C,2)*100; %rows are actual, cols are estimated
C_Festimated=C./sum(C)*100;
C_Total = (C./sum(C(:)))*100;

if 1
    AccA_C_Factual = C_Factual.*diag(ones(1,length(C)));
    AccB_C_Factual = C_Factual.*(diag(ones(1,length(C)-1),1)+diag(ones(1,length(C)-1),-1));
    AccAactual = mean(sum(AccA_C_Factual,2));
    AccBactual = mean(sum(AccA_C_Factual,2)) + mean(sum(AccB_C_Factual,2));
    AccA_C_Festimated = C_Festimated.*diag(ones(1,length(C)));
    AccB_C_Festimated = C_Festimated.*(diag(ones(1,length(C)-1),1)+diag(ones(1,length(C)-1),-1));
    AccAestimated = mean(sum(AccA_C_Festimated));
    AccBestimated = mean(sum(AccA_C_Festimated)) + mean(sum(AccB_C_Festimated));
    ACCs = [AccAactual AccBactual; AccAestimated AccBestimated]
    % first row is actual, second row is estimate - use second row
    %second col gives accuracy if accepting next-category error as correct
end

x = order'; y = order';
if 1; C1 = C_Festimated; else; C1 = C_Factual; end
%C1 = C_Total;
xSplit = diff(x)/2;                 % Find edge points
ySplit = diff(y)/2;
xEdges = [x(1)-xSplit(1) x(2:end)-xSplit x(end)+xSplit(end)];
yEdges = [y(1)-ySplit(1) y(2:end)-ySplit y(end)+ySplit(end)];
[XGrid, YGrid] = meshgrid(xEdges,yEdges);
YGrid = flipud(YGrid);              % To match expected behavior
XGrid = fliplr(XGrid);
C2 = [[C1 zeros(size(C1,1),1)] ; zeros(1,size(C1,2)+1)];% Last row/col ignored
pcolor(XGrid,YGrid,(1-(1-C2/100).^2)*100)
hold on                             % Plot original data points
[X,Y] = meshgrid(x,y);
set(gcf,'colormap',customcmap);
C1 = flipud(C1);
C1 = fliplr(C1);
C1 = C1';
for i=1:size(C1,1)
    for j=1:size(C1,2)
        if C1(i,j)<((max(max(C1))-min(min(C1)))/2+min(min(C1)))
            textcolor=[1 1 1];
        else
            textcolor=[0 0 0];
        end
        text(x(i),y(j),num2str(round(C1(i,j))),'color',textcolor,'horizontalalignment','center','fontname','arial narrow')
    end
end
ax = gca; set(ax,'fontname','arial narrow','FontSize', TickFntSz);
yticks(y); yticklabels(gca,fliplr(labeltext));
xticks(x); xticklabels(gca,fliplr(labeltext));
xlabel('Flow Shape Classification', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
ylabel('Gold Standard Classification', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
axis square

% % Add labels A B etc to plot space
% subplot(1,2,1); hold on;
% text(-15, 0.99, 'A', 'FontSize', 20, 'FontWeight', 'Bold');
% subplot(1,2,2); hold on;
% text(-15, 0.99, 'B', 'FontSize', 20, 'FontWeight', 'Bold');

str = ['..\Figures\Figure_E3',settings.experimentnumber];

switch savefigas
    case 'saveasTIFF'; print(fig, str, '-dtiff', '-r1000');
    case 'saveasPNG'; saveas(fig, str, 'png');
    case 'saveasFIG'; savefig(str);
end
if closefigs; close(fig); end

SummaryStats_FigE3 = [RL1O_flow(ftrnum), RsqL1O_flow(ftrnum),  AccAestimated, AccBestimated];

%% Figure E4 - pnasal
figure(84); clf(figure(84)); fig = gcf;
fig.Color = [1 1 1]; fig.Units = 'inches';
fig.Position = [22  0.5   12   4.5];

% performance vs number of features
subplot(1,2,1); % pnasal
%stairs([RL1O_pnasal(1:PlotUptoNFtrs);RsqL1O_pnasal(1:PlotUptoNFtrs);ErrL1O_pnasal(1:PlotUptoNFtrs)]'); hold on;
stairs([1-RsqL1O_pnasal(1:PlotUptoNFtrs);ErrL1O_pnasal(1:PlotUptoNFtrs)]'); hold on;
%currentYlim=ylim();
%plot([ftrnum, ftrnum], [0, currentYlim(2)],'k-', 'linewidth', 1.5);
plot([ftrnum, ftrnum], [0, 0.6],'k-', 'linewidth', 1.5);
ax = gca; set(ax,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz);
xlabel('Number of Features', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
ylabel('Performance (Testing)', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow'); 
ylim([0 0.7]); xlim([0 50]); axis square;
%text(30, 0.62, {'correlation coefficient'}, 'fontname','arial narrow','FontSize',10, 'FontWeight', 'normal');
text(30, 0.55, {'1 - coefficient of determination'}, 'fontname','arial narrow','FontSize',10, 'FontWeight', 'normal');
text(30, 0.23, {'mean absolute error'}, 'fontname','arial narrow','FontSize',10, 'FontWeight', 'normal');


% % proportion of total error
% subplot(1,3,2);
% % nothing to show, no training error for pnasal, as no training in pnasal
% ax = gca; set(ax,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz);
% xlabel('Number of Features', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
% ylabel('Proportion of absolute error (%)', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
% ylim([0 30]); xlim([0 50]); axis square
% 
% text(15, 15, {'no training data for pnasal'}, 'fontname','arial narrow','FontSize',10, 'FontWeight', 'normal');

% confusion mat
subplot(1,2,2);
customcmap = GetCustomColorMap('gray'); % 'SS'
g = NaN*Yval_pnasal;
ghat = NaN*predy_pnasal;
classes = 5;
switch classes
    case 5
        classcutoffs = [0.9 0.7 0.5 0.3];
        labeltext={'Normal','Mild','Moderate','Severe','V.Severe'};
    case 4
        classcutoffs = [0.9 0.5 0.3];
        labeltext={'Normal','Mild','Moderate','Severe'};
    case 3
        classcutoffs = [0.9 0.5];
        labeltext={'Normal','Intermediate','Flow Limited'};
    otherwise
        classcutoffs = [0.7]; %normal mild moderate severe v.severe / normal borderline mild moderate severe
        labeltext={'Normal','Flow Limited'};
end
Nclasses=length(classcutoffs)+1;
g(Yval_pnasal>classcutoffs(1))=1;
ghat(predy_pnasal>classcutoffs(1))=1;
for i=2:Nclasses
    g(Yval_pnasal<classcutoffs(i-1))=i;
    ghat(predy_pnasal<classcutoffs(i-1))=i;
end
[C,order] = confusionmat(g,ghat);
%sumactual=sum(C')';
sumactual=sum(C,2);
sumestimated=sum(C);
%C_Factual=C./sum(C')'
C_Factual=C./sum(C,2)*100; %rows are actual, cols are estimated
C_Festimated=C./sum(C)*100;
C_Total = (C./sum(C(:)))*100;
 
if 1
    AccA_C_Factual = C_Factual.*diag(ones(1,length(C)));
    AccB_C_Factual = C_Factual.*(diag(ones(1,length(C)-1),1)+diag(ones(1,length(C)-1),-1));
    AccAactual = mean(sum(AccA_C_Factual,2));
    AccBactual = mean(sum(AccA_C_Factual,2)) + mean(sum(AccB_C_Factual,2));
    AccA_C_Festimated = C_Festimated.*diag(ones(1,length(C)));
    AccB_C_Festimated = C_Festimated.*(diag(ones(1,length(C)-1),1)+diag(ones(1,length(C)-1),-1));
    AccAestimated = mean(sum(AccA_C_Festimated));
    AccBestimated = mean(sum(AccA_C_Festimated)) + mean(sum(AccB_C_Festimated));
    ACCs = [AccAactual AccBactual; AccAestimated AccBestimated]
    % first row is actual, second row is estimate - use second row
    %second col gives accuracy if accepting next-category error as correct
end

x = order'; y = order';
if 1; C1 = C_Festimated; else; C1 = C_Factual; end
%C1 = C_Total;
xSplit = diff(x)/2;                 % Find edge points
ySplit = diff(y)/2;
xEdges = [x(1)-xSplit(1) x(2:end)-xSplit x(end)+xSplit(end)];
yEdges = [y(1)-ySplit(1) y(2:end)-ySplit y(end)+ySplit(end)];
[XGrid, YGrid] = meshgrid(xEdges,yEdges);
YGrid = flipud(YGrid);              % To match expected behavior
XGrid = fliplr(XGrid);
C2 = [[C1 zeros(size(C1,1),1)] ; zeros(1,size(C1,2)+1)];% Last row/col ignored
pcolor(XGrid,YGrid,(1-(1-C2/100).^2)*100)
hold on                             % Plot original data points
[X,Y] = meshgrid(x,y);
set(gcf,'colormap',customcmap);
C1 = flipud(C1);
C1 = fliplr(C1);
C1 = C1';
for i=1:size(C1,1)
    for j=1:size(C1,2)
        if C1(i,j)<((max(max(C1))-min(min(C1)))/2+min(min(C1)))
            textcolor=[1 1 1];
        else
            textcolor=[0 0 0];
        end
        text(x(i),y(j),num2str(round(C1(i,j))),'color',textcolor,'horizontalalignment','center','fontname','arial narrow')
    end
end
ax = gca; set(ax,'fontname','arial narrow','FontSize', TickFntSz);
yticks(y); yticklabels(gca,fliplr(labeltext));
xticks(x); xticklabels(gca,fliplr(labeltext));
xlabel('Flow Shape Classification', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
ylabel('Gold Standard Classification', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
axis square;

str = ['..\Figures\Figure_E4',settings.experimentnumber]; 

switch savefigas     
    case 'saveasTIFF'; print(fig, str, '-dtiff', '-r1000');     
    case 'saveasPNG'; saveas(fig, str, 'png');      
    case 'saveasFIG'; savefig(str);  
end 
if closefigs; close(fig); end

SummaryStats_FigE4 = [RL1O_pnasal(ftrnum), RsqL1O_pnasal(ftrnum),  AccAestimated, AccBestimated];

%% Add apnea and low flow breaths back in for everything hereafter
RemovedBB_pnasal = load([settings.datadir, 'PnasalDrive_only25Hz_ReNormalized_FeatureSpace_AutoRef2_Edi_Clean.mat'],'RemovedBB_Apnoea','RemovedBB_LowFlow');
RemovedBB_flow = load([settings.datadir, 'FlowDrive_only25Hz_ReNormalized_FeatureSpace_AutoRef2_Edi_Clean.mat'],'RemovedBB_Apnoea','RemovedBB_LowFlow');

numBBinTest_flow = NaN(54,6); numBBinTest_pnasal = NaN(54,6);
PredY_pnasal = []; Gtest_pnasal = []; PTarray_pnasal=[];
PredY_flow = []; Gtest_flow = []; PTarray_flow=[];
Gtest_flow_avg = NaN(54,1); PredY_flow_avg = NaN(54,1);
Gtest_pnasal_avg = NaN(54,1); PredY_pnasal_avg = NaN(54,1);
Gtest_flow_avg_NoAp = NaN(54,1); Gtest_pnasal_avg_NoAp = NaN(54,1);
for subj=1:54
    if ismember(subj, settings.Pnasal_list) % add in pnasal RemBB
        Isubj=(PtData_pnasal.PT==subj); % find all the breaths that belong to this pt, this is just for counting
        numBBinTest_pnasal(subj,1) = nnz(Isubj); % the total number of breaths for this pt
        Isubj=(PtData_pnasal.PT==subj)&(PtData_pnasal.Hypnog<4)&(PtData_pnasal.Ar==0); % find the sleep only breaths that belong to this pt
        numBBinTest_pnasal(subj,2) = nnz(Isubj); % the number of sleep breaths for this pt
        Gtest_pt_pnasal = GS(Isubj); Gtest_pnasal_avg_NoAp(subj) = median(Gtest_pt_pnasal);
        PredY_pt_pnasal = predyL1O_array_pnasal(Isubj,ftrnum);
        % add apnea breaths as 0.1 to both Gtest_pt and PredY_pt
        indAP = find(RemovedBB_pnasal.RemovedBB_Apnoea.Pt == subj);
        indLF = find(RemovedBB_pnasal.RemovedBB_LowFlow.Pt == subj);
        numApBB = RemovedBB_pnasal.RemovedBB_Apnoea{indAP,2};
        if isempty(numApBB); numApBB=0; end
        numBBinTest_pnasal(subj,3) = numApBB;
        numLFBB = RemovedBB_pnasal.RemovedBB_LowFlow{indLF,2};
        if isempty(numLFBB); numLFBB=0; end
        numBBinTest_pnasal(subj,4) = numLFBB;
        numRemBB = numApBB + numLFBB;
        PredY_pt_pnasal = [PredY_pt_pnasal;ones(numRemBB,1)*0.09];
        numBBinTest_pnasal(subj,5) = length(PredY_pt_pnasal);
        Gtest_pt_pnasal = [Gtest_pt_pnasal;ones(numRemBB,1)*0.09];
        numBBinTest_pnasal(subj,6) = length(Gtest_pt_pnasal);
        Gtest_pnasal_avg(subj) = median(Gtest_pt_pnasal);
        PredY_pnasal_avg(subj) = median(PredY_pt_pnasal);
        PTarray_pnasal = [PTarray_pnasal; ones(length(PredY_pt_pnasal),1)*subj]; % same as PtData.PT but includes apnea breaths
        PredY_pnasal = [PredY_pnasal; PredY_pt_pnasal];
        Gtest_pnasal = [Gtest_pnasal; Gtest_pt_pnasal];
    end
    
    if ismember(subj, settings.Flow_list) % add in flow RemBB
        Isubj=(PtData_flow.PT==subj); % find all the breaths that belong to this pt, this is just for counting
        numBBinTest_flow(subj,1) = nnz(Isubj); % the total number of breaths for this pt
        Isubj=(PtData_flow.PT==subj)&(PtData_flow.Hypnog<4)&(PtData_flow.Ar==0); % find the sleep only breaths that belong to this pt
        numBBinTest_flow(subj,2) = nnz(Isubj); % the number of sleep breaths for this pt
        Gtest_pt_flow = GS(Isubj); Gtest_flow_avg_NoAp(subj) = median(Gtest_pt_flow);
        PredY_pt_flow = predyL1O_array_flow(Isubj,ftrnum);
        % add apneoa breaths as 0.1 to both Gtest_pt and PredY_pt
        indAP = find(RemovedBB_flow.RemovedBB_Apnoea.Pt == subj);
        indLF = find(RemovedBB_flow.RemovedBB_LowFlow.Pt == subj);
        numApBB = RemovedBB_flow.RemovedBB_Apnoea{indAP,2};
        if isempty(numApBB); numApBB=0; end
        numBBinTest_flow(subj,3) = numApBB;
        numLFBB = RemovedBB_flow.RemovedBB_LowFlow{indLF,2};
        if isempty(numLFBB); numLFBB=0; end
        numBBinTest_flow(subj,4) = numLFBB;
        numRemBB = numApBB + numLFBB;
        PredY_pt_flow = [PredY_pt_flow;ones(numRemBB,1)*0.09];
        numBBinTest_flow(subj,5) = length(PredY_pt_flow);
        Gtest_pt_flow = [Gtest_pt_flow;ones(numRemBB,1)*0.09];
        numBBinTest_flow(subj,6) = length(Gtest_pt_flow);        
        Gtest_flow_avg(subj) = median(Gtest_pt_flow);
        PredY_flow_avg(subj) = median(PredY_pt_flow);
        PTarray_flow = [PTarray_flow; ones(length(PredY_pt_flow),1)*subj]; % same as PtData.PT but includes apnea breaths
        PredY_flow = [PredY_flow; PredY_pt_flow];
        Gtest_flow = [Gtest_flow; Gtest_pt_flow];
    end
end

% tidy up table, remove nan rows
% rows in numBBinTest are each patient
% cols in numBBinTest are:
% Total BB, Sleep BB, Ap BB, Low flow BB, Total BB to use (with AP)
numBBinTest_flow = numBBinTest_flow(settings.Flow_list,:);
numBBinTest_pnasal = numBBinTest_pnasal(settings.Pnasal_list,:);
Gtest_flow_avg_NoAp=Gtest_flow_avg_NoAp(settings.Flow_list);
Gtest_flow_avg=Gtest_flow_avg(settings.Flow_list);
PredY_flow_avg=PredY_flow_avg(settings.Flow_list);
Gtest_pnasal_avg=Gtest_pnasal_avg(settings.Pnasal_list);
PredY_pnasal_avg=PredY_pnasal_avg(settings.Pnasal_list);
Gtest_pnasal_avg_NoAp=Gtest_pnasal_avg_NoAp(settings.Pnasal_list);

%% get hypnogram data for breaths
% for subj=1:54
%      if ismember(subj, settings.Flow_list) % add in flow RemBB
%         Isubj=(PtData_flow.PT==subj);
%      end
%      
%      
% end
% Hypnog_flow = []; 
% Hypnog_pnasal = [];


%% checking and cleaning
if 0
    if length(PredY_pnasal) ~= length(Gtest_pnasal); keyboard; end              % check same length
    if nnz(isnan(PredY_pnasal))>0 || nnz(isnan(GS))>0; keyboard; end            % check for any NaN's
    if nnz(PredY_pnasal<0)>0; Predy(PredY_pnasal<0)=0; end                      % force lower limit
    if nnz(PredY_pnasal>maxG)>0; Predy(PredY_pnasal>maxG)=1.5; end              % force upper limit
end

%% Figure 3 - flow
figure(3); clf(figure(3)); fig = gcf;
fig.Color = [1 1 1]; fig.Units = 'inches';
fig.Position = [22   5  12   4.5];

% scatter with box overlay
subplot(1,2,1);
scatter(100*predy_flow,100*Yval_flow,2,'filled','markerfacealpha',0.06, 'markerfacecolor', [0 0 0]); hold on;
dx=0.2; xbins=[0 0.3:dx:0.9 1.5];
clear medianX medianY upperIQRY lowerIQRY upperIQRY2 lowerIQRY2
for i=1:length(xbins)-1
    Ix=predy_flow>xbins(i)&predy_flow<xbins(i+1);
    medianX(i)=prctile(predy_flow(Ix),50);
    medianY(i)=prctile(Yval_flow(Ix),50);
    upperIQRY(i)=prctile(Yval_flow(Ix),75);
    lowerIQRY(i)=prctile(Yval_flow(Ix),25);
    upperIQRY2(i)=prctile(Yval_flow(Ix),90);
    lowerIQRY2(i)=prctile(Yval_flow(Ix),10);
end
plotbarwitherrorsmedian(100*medianY,100*medianX,100*upperIQRY2,100*lowerIQRY2,zeros(length(medianY),3),NaN,100*0.005,100*0.01);
plotbarwitherrorsmedian(100*medianY,100*medianX,100*upperIQRY,100*lowerIQRY,zeros(length(medianY),3),NaN,100*dx/2,100*0.01);
xlim([0 150]); 
ax = gca;set(ax,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz);
ax.XTick=[0:25:150]; xticklabels(ax, {'0', '25', '50', '75', '100', '125'});
ax.YTick=[0:25:150]; yticklabels(ax, {'0', '25', '50', '75', '100', '125'});
xlabel('Flow Shape Predicted {\itflow:drive} (%)', 'FontSize', FntSz, 'FontName', 'Arial Narrow');
ylabel('Gold Standard {\itflow:drive} (%)', 'FontSize', FntSz, 'FontName', 'Arial Narrow');
%[F6B_r_1,~]=corr(predy_flow,Yval_flow);  % don't use corr if we have weights
%[F6B_Rsq_2,~,~,~,~]=glmfitFast(predy_flow,Yval_flow,weights_flow,1);
[F4A_Rsq, F4A_r] = UnivariateStats(predy_flow,Yval_flow,weights_flow);
[F4A_Rsq_uw, F4A_r_uw] = UnivariateStats(predy_flow,Yval_flow);
text(110, 15, ['R^2 = ', num2str(F4A_Rsq,sigdig)]); 
%text(110, 5, ['r = ', num2str(round(F4A_r,2))]); 
text(-45, 148, 'A', 'FontSize', 20, 'FontWeight', 'Bold'); % Add label to plot space
axis square

subplot(1,2,2);
[F4B_r, ~] = plotregressionwithSEM(PredY_flow_avg.*100, Gtest_flow_avg.*100);
%[F6C_r, ~] =corr(PredY_flow_avg,Gtest_flow_avg); % don't need corr if have plotreg
w = ones(length(Gtest_flow_avg),1);
[F4B_Rsq,~,~,~,~]=glmfitFast(PredY_flow_avg, Gtest_flow_avg,w,1);
%[F6C_Rsq_2, F6C_r_2] = UnivariateStats(PredY_flow_avg,Gtest_flow_avg); % don't use without weights
xlim([-5 110]); ylim([-5 110]);
ax = gca; set(ax,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz);
ax.XTick=[0:25:125];
ax.YTick=[0:25:125];
xlabel({'Flow Shape Predicted {\itflow:drive}','median (%)'}, 'FontSize', FntSz, 'FontName', 'Arial Narrow');
ylabel({'Gold Standard {\itflow:drive}','median (%)'}, 'FontSize', FntSz, 'FontName', 'Arial Narrow'); 
axis square
text(80, 7, ['R^2 = ', num2str(F4B_Rsq,sigdig)]); 
%text(80, 0, ['r = ', num2str(round(F4B_r,2))]); 
text(-45, 108, 'B', 'FontSize', 20, 'FontWeight', 'Bold'); % Add label to plot space

str = ['..\Figures\Figure_3',settings.experimentnumber];

switch savefigas
    case 'saveasTIFF'; print(fig, str, '-dtiff', '-r1000');
    case 'saveasPNG'; saveas(fig, str, 'png');
    case 'saveasFIG'; savefig(str);
end
if closefigs; close(fig); end

SummaryStats_Fig3 = [NaN, NaN, F4A_r, F4A_Rsq, F4B_r, F4B_Rsq];

% per breath
%100*mean(ErrTrain_array_flow(settings.Flow_list,ftrnum)) % training mean absolute error
[F4A_MAE] = CalcMAE(predy_flow,Yval_flow); % all breaths
% 100*mean(abs(PredY_flow-Gtest_flow)) % sleep only breaths, including Ap.

% per patient
% 100*mean(abs(PredY_flow_avg-Gtest_flow_avg)) 
[F4B_MAE] = CalcMAE(PredY_flow_avg, Gtest_flow_avg);


%% Figure 5 - pnasal
figure(5); clf(figure(5)); fig = gcf;
fig.Color = [1 1 1]; fig.Units = 'inches';
fig.Position = [22  -0.7   16   5];

% flow predy vs pnasal predy
subplot(1,3,1)
facealpha = 0.05; % was 0.04
facecolor = [0 0 0]; %was [0 0 0]
scatter(predy_pnasal.*100,predy_flow(BB_).*100,2,'filled','markerfacealpha',facealpha, 'markerfacecolor', facecolor); hold on;
%[F6A_r, ~] =corr(predy_pnasal,predy_flow(BB_)); 
%w = ones(length(predy_pnasal),1);
%[F6A_Rsq,~,~,~,~]=glmfitFast(predy_pnasal,predy_flow(BB_),w,1);
[F6A_Rsq, F6A_r] = UnivariateStats(predy_pnasal,predy_flow(BB_),weights_pnasal); % don't use without weights
text(110, 15, ['R^2 = ', num2str(F6A_Rsq,sigdig)]); 
%text(110, 5, ['r = ', num2str(round(F6A_r,2))]);
xlim([0 150]); ylim([0 150]); 
ax = gca; set(ax,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz);
ax.XTick=[0:25:150]; xticklabels(ax, {'0', '25', '50', '75', '100', '125'});
ax.YTick=[0:25:150]; yticklabels(ax, {'0', '25', '50', '75', '100', '125'});
xlabel('Pnasal Shape Predicted {\itflow:drive} (%)', 'FontSize', FntSz, 'FontName', 'Arial Narrow');
ylabel('Flow Shape Predicted {\itflow:drive} (%)', 'FontSize', FntSz, 'FontName', 'Arial Narrow');
text(-45, 148, 'A', 'FontSize', 20, 'FontWeight', 'Bold');
axis square

% scatter with box overlay
subplot(1,3,2);
scatter(100*predy_pnasal,100*Yval_pnasal,2,'filled','markerfacealpha',0.06, 'markerfacecolor', [0 0 0]); hold on;
dx=0.2; xbins=[0 0.3:dx:0.9 1.5];
clear medianX medianY upperIQRY lowerIQRY upperIQRY2 lowerIQRY2
for i=1:length(xbins)-1
    Ix=predy_pnasal>xbins(i)&predy_pnasal<xbins(i+1);
    medianX(i)=prctile(predy_pnasal(Ix),50);
    medianY(i)=prctile(Yval_pnasal(Ix),50);
    upperIQRY(i)=prctile(Yval_pnasal(Ix),75);
    lowerIQRY(i)=prctile(Yval_pnasal(Ix),25);
    upperIQRY2(i)=prctile(Yval_pnasal(Ix),90);
    lowerIQRY2(i)=prctile(Yval_pnasal(Ix),10);
end
plotbarwitherrorsmedian(100*medianY,100*medianX,100*upperIQRY2,100*lowerIQRY2,zeros(length(medianY),3),NaN,100*0.005,100*0.01);
plotbarwitherrorsmedian(100*medianY,100*medianX,100*upperIQRY,100*lowerIQRY,zeros(length(medianY),3),NaN,100*dx/2,100*0.01);
xlim([0 150]); 
ax = gca;set(ax,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz);
ax.XTick=[0:25:150]; xticklabels(ax, {'0', '25', '50', '75', '100', '125'});
ax.YTick=[0:25:150]; yticklabels(ax, {'0', '25', '50', '75', '100', '125'});
xlabel('Pnasal Shape Predicted {\itflow:drive} (%)', 'FontSize', FntSz, 'FontName', 'Arial Narrow');
ylabel('Gold Standard {\itflow:drive} (%)', 'FontSize', FntSz, 'FontName', 'Arial Narrow');
%[F6B_r_1,~]=corr(predy_pnasal,Yval_pnasal);  % don't use corr if we have weights
%[F6B_Rsq_2,~,~,~,~]=glmfitFast(predy_pnasal,Yval_pnasal,weights_pnasal,1);
[F6B_Rsq, F6B_r] = UnivariateStats(predy_pnasal,Yval_pnasal,weights_pnasal);
text(110, 15, ['R^2 = ', num2str(F6B_Rsq,sigdig)]); 
%text(110, 5, ['r = ', num2str(round(F6B_r,2))]); 
text(-45, 148, 'B', 'FontSize', 20, 'FontWeight', 'Bold'); % Add label to plot space
axis square

subplot(1,3,3);
[F6C_r, ~] = plotregressionwithSEM(PredY_pnasal_avg.*100, Gtest_pnasal_avg.*100);
%[F6C_r, ~] =corr(PredY_pnasal_avg,Gtest_pnasal_avg); % don't need corr if have plotreg
w = ones(length(Gtest_pnasal_avg),1);
[F6C_Rsq,~,~,~,~]=glmfitFast(PredY_pnasal_avg, Gtest_pnasal_avg,w,1);
%[F6C_Rsq_2, F6C_r_2] = UnivariateStats(PredY_pnasal_avg,Gtest_pnasal_avg); % don't use without weights
xlim([-5 110]); ylim([-5 110]);
ax = gca; set(ax,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz);
ax.XTick=[0:25:125];
ax.YTick=[0:25:125];
xlabel({'Pnasal Shape Predicted {\itflow:drive}','median (%)'}, 'FontSize', FntSz, 'FontName', 'Arial Narrow');
ylabel({'Gold Standard {\itflow:drive}','median (%)'}, 'FontSize', FntSz, 'FontName', 'Arial Narrow'); 
axis square
text(80, 7, ['R^2 = ', num2str(F6C_Rsq,sigdig)]); 
%text(80, 0, ['r = ', num2str(round(F6C_r,2))]); 
text(-45, 108, 'C', 'FontSize', 20, 'FontWeight', 'Bold'); % Add label to plot space

str = ['..\Figures\Figure_5',settings.experimentnumber];

switch savefigas
    case 'saveasTIFF'; print(fig, str, '-dtiff', '-r1000');
    case 'saveasPNG'; saveas(fig, str, 'png');
    case 'saveasFIG'; savefig(str);
end
if closefigs; close(fig); end

SummaryStats_Fig5 = [F6A_r, F6A_Rsq, F6B_r, F6B_Rsq, F6C_r, F6C_Rsq];

% flow vs pnasal
[F6A_MAE] = CalcMAE(predy_pnasal,predy_flow(BB_)); 

% per breath
[F6B_MAE] = CalcMAE(predy_pnasal,Yval_pnasal); % all breaths
% 100*mean(abs(PredY_pnasal-Gtest_pnasal)) % sleep only breaths, including Ap.

% per patient
% 100*mean(abs(PredY_flow_avg-Gtest_flow_avg)) 
[F6C_MAE] = CalcMAE(PredY_flow_avg, Gtest_flow_avg);

%%  Figure E2 - flow VE  Vs pnasal VE
pnasal_veveup = PtData_pnasal.VE(BB_) ./ PtData_pnasal.Veup(BB_);
flow_veveup = PtData_flow.VE(BB_) ./ PtData_flow.Veup(BB_);
pnasal_veveup_t = pnasal_veveup.^1.5;

squareplot = 1;
facealpha = 0.05;
facecolor = [0.1 0.1 0.1];

figure(82); clf(figure(82)); fig = gcf;
fig.Color = [1 1 1]; fig.Units = 'inches';
if squareplot
    fig.Position = [22  1   8   8];
else
    fig.Position = [0.5  0.5   12  4];
end

%0.5
if squareplot; subplot(2,2,1); else; subplot(1,4,1); end
pnasal_veveup_05 = pnasal_veveup_t.^0.5;
scatter(pnasal_veveup_05, flow_veveup,2,'filled','markerfacealpha',facealpha, 'markerfacecolor', facecolor); hold on;
[SmMed_p, SmMed_f, SmBias_050, ~, LgMed_p, LgMed_f, LgBias_050, ~] = ...
    getBiasForPnasalVsFlow(pnasal_veveup_05, flow_veveup);
plot(SmMed_p, SmMed_f, 'rs'); plot(LgMed_p, LgMed_f, 'rs');
[r_1, p_1] =corr(pnasal_veveup_05, flow_veveup);
[rsq, r_dm] = UnivariateStats(pnasal_veveup_05, flow_veveup, weights(BB_));
mae_1 = CalcMAE(pnasal_veveup_05, flow_veveup);
xlim([-0.1 2.5]); ylim([-0.1 2.5]);
ax = gca; set(ax,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz);
ax.XTick=[0:0.5:2.0]; xticklabels(ax, {'0', '50', '100', '150', '200'});
ax.YTick=[0:0.5:2.0]; yticklabels(ax, {'0', '50', '100', '150', '200'});
xlabel(['Nasal pressure ventilation (% eupnea)']); 
ylabel(['Pneumotach ventilation (% eupnea)']);
text(1.35, 0.4, ['R^2 = ', num2str(rsq,sigdig)]);
text(1.35, 0.2, ['Bias (small) = ', num2str(SmBias_050,sigdig)]);
text(1.35, 0.0, ['Bias (large) = ', num2str(LgBias_050,sigdig)]);
%hndl_ls = lsline;  hndl_ls.LineStyle = '-'; hndl_ls.Color = [0 0 0];
hndl_ref = refline(1,0); hndl_ref.LineStyle = '--'; hndl_ref.Color = [0.25 0.25 0.25];
xlim([-0.1 2.5]); ylim([-0.1 2.5]);
%title(['0.5']);
axis square

%0.67
if squareplot; subplot(2,2,2); else; subplot(1,4,2); end
scatter(pnasal_veveup, flow_veveup,2,'filled','markerfacealpha',facealpha, 'markerfacecolor', facecolor); hold on;
[SmMed_p, SmMed_f, SmBias_067, ~, LgMed_p, LgMed_f, LgBias_067, ~] = ...
    getBiasForPnasalVsFlow(pnasal_veveup, flow_veveup);
plot(SmMed_p, SmMed_f, 'rs'); plot(LgMed_p, LgMed_f, 'rs');
[r_1, p_1] =corr(pnasal_veveup, flow_veveup);
[rsq, r_dm] = UnivariateStats(pnasal_veveup, flow_veveup, weights(BB_));
mae_1 = CalcMAE(pnasal_veveup, flow_veveup);
xlim([-0.1 2.5]); ylim([-0.1 2.5]);
ax = gca; set(ax,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz);
ax.XTick=[0:0.5:2.0]; xticklabels(ax, {'0', '50', '100', '150', '200'});
ax.YTick=[0:0.5:2.0]; yticklabels(ax, {'0', '50', '100', '150', '200'});
xlabel(['Nasal pressure ventilation (% eupnea)']); 
ylabel(['Pneumotach ventilation (% eupnea)']);
text(1.35, 0.4, ['R^2 = ', num2str(rsq,sigdig)]);
text(1.35, 0.2, ['Bias (small) = ', num2str(SmBias_067,sigdig)]);
text(1.35, 0.0, ['Bias (large) = ', num2str(LgBias_067,sigdig)]);
clear hndl_ls;
%hndl_ls = lsline;  hndl_ls.LineStyle = '-'; hndl_ls.Color = [0 0 0];
hndl_ref = refline(1,0); hndl_ref.LineStyle = '--'; hndl_ref.Color = [0.25 0.25 0.25];
%title(['0.67']);
axis square

%0.75
if squareplot; subplot(2,2,3); else; subplot(1,4,3); end
pnasal_veveup_075 = pnasal_veveup_t.^0.75;
scatter(pnasal_veveup_075, flow_veveup,2,'filled','markerfacealpha',facealpha, 'markerfacecolor', facecolor); hold on;
[SmMed_p, SmMed_f, SmBias_075, ~, LgMed_p, LgMed_f, LgBias_075, ~] = ...
    getBiasForPnasalVsFlow(pnasal_veveup_075, flow_veveup);
plot(SmMed_p, SmMed_f, 'rs'); plot(LgMed_p, LgMed_f, 'rs');
[r_1, p_1] =corr(pnasal_veveup_075, flow_veveup);
[rsq, r_dm] = UnivariateStats(pnasal_veveup_075, flow_veveup, weights(BB_));
mae_1 = CalcMAE(pnasal_veveup_075, flow_veveup);
xlim([-0.1 2.5]); ylim([-0.1 2.5]);
ax = gca; set(ax,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz);
ax.XTick=[0:0.5:2.0]; xticklabels(ax, {'0', '50', '100', '150', '200'});
ax.YTick=[0:0.5:2.0]; yticklabels(ax, {'0', '50', '100', '150', '200'});
xlabel(['Nasal pressure ventilation (% eupnea)']); 
ylabel(['Pneumotach ventilation (% eupnea)']);
text(1.35, 0.4, ['R^2 = ', num2str(rsq,sigdig)]);
text(1.35, 0.2, ['Bias (small) = ', num2str(SmBias_075,sigdig)]);
text(1.35, 0.0, ['Bias (large) = ', num2str(LgBias_075,sigdig)]);
%hndl_ls = lsline;  hndl_ls.LineStyle = '-'; hndl_ls.Color = [0 0 0];
hndl_ref = refline(1,0); hndl_ref.LineStyle = '--'; hndl_ref.Color = [0.25 0.25 0.25];
%title(['0.75']);
axis square

%1
if squareplot; subplot(2,2,4); else; subplot(1,4,4); end
pnasal_veveup_1 = pnasal_veveup_t.^1;
scatter(pnasal_veveup_1, flow_veveup,2,'filled','markerfacealpha',facealpha, 'markerfacecolor', facecolor); hold on;
[SmMed_p, SmMed_f, SmBias_100, ~, LgMed_p, LgMed_f, LgBias_100, ~] = ...
    getBiasForPnasalVsFlow(pnasal_veveup_1, flow_veveup);
plot(SmMed_p, SmMed_f, 'rs'); plot(LgMed_p, LgMed_f, 'rs');
[r_1, p_1] =corr(pnasal_veveup_1, flow_veveup);
[rsq, r_dm] = UnivariateStats(pnasal_veveup_1, flow_veveup, weights(BB_));
mae_1 = CalcMAE(pnasal_veveup_1, flow_veveup);
xlim([-0.1 2.5]); ylim([-0.1 2.5]);
ax = gca; set(ax,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz);
ax.XTick=[0:0.5:2.0]; xticklabels(ax, {'0', '50', '100', '150', '200'});
ax.YTick=[0:0.5:2.0]; yticklabels(ax, {'0', '50', '100', '150', '200'});
xlabel(['Nasal pressure ventilation (% eupnea)']); 
ylabel(['Pneumotach ventilation (% eupnea)']);
text(1.35, 0.4, ['R^2 = ', num2str(rsq,sigdig)]);
text(1.35, 0.2, ['Bias (small) = ', num2str(SmBias_100,sigdig)]);
text(1.35, 0.0, ['Bias (large) = ', num2str(LgBias_100,sigdig)]);
%hndl_ls = lsline;  hndl_ls.LineStyle = '-'; hndl_ls.Color = [0 0 0];
hndl_ref = refline(1,0); hndl_ref.LineStyle = '--'; hndl_ref.Color = [0.25 0.25 0.25];
%title([1']);
axis square

% Add labels A B to plot space
subplot(2,2,1); hold on;
text(-0.75, 2.6, 'A', 'FontSize', 20, 'FontWeight', 'Bold');
subplot(2,2,2); hold on;
text(-0.75, 2.6, 'B', 'FontSize', 20, 'FontWeight', 'Bold');
subplot(2,2,3); hold on;
text(-0.75, 2.6, 'C', 'FontSize', 20, 'FontWeight', 'Bold');
subplot(2,2,4); hold on;
text(-0.75, 2.6, 'D', 'FontSize', 20, 'FontWeight', 'Bold');

str = ['..\Figures\', 'Figure_E2_wBias'];
switch savefigas     
    case 'saveasTIFF'; print(fig, str, '-dtiff', '-r1000');     
    case 'saveasPNG'; saveas(fig, str, 'png');      
    case 'saveasFIG'; savefig(str);  
end
if closefigs; close(fig); end

BiasResults = [SmBias_050, SmBias_067, SmBias_075, SmBias_100;
LgBias_050, LgBias_067, LgBias_075, LgBias_100];

BiasResults_Tbl = array2table(round(BiasResults,2), ...
    'VariableNames', {'half', 'twothirds' 'threequart' 'one'}, ...
    'RowName', {'smallBB', 'largeBB'});

%% Histogram of ftrs selected during training 
ftr_array = NaN(54, ftrnum);
for subj=1:54 % set pt num or all, no data for PT=1
    if ~ismember(subj, settings.Flow_list) %PT_list)
        continue
    end
    ftr_array(subj,:) = labels_Step_Subj_flow{subj,ftrnum}; % only one list of labels, same for both flow and pnasal
end
ftr_array = ftr_array(settings.Flow_list,:);
ftr_array_linear = ftr_array(:);
ctrs = 1:1:size(Amatrix2_flow,2);  % ToDo: make sure max ctrs covers all features
[counts, ~] = hist(ftr_array_linear, ctrs);
[~,I_uw] = sort(counts, 'descend');

figure(33); clf(figure(33)); fig = gcf;
fig.Color = [1 1 1]; fig.Units = 'inches';
fig.Position = [1   4   12    4.5];
bar(counts(I_uw), 'facealpha',0.2, 'edgealpha', 0.2); hold on;
xlim([0 PlotUptoNFtrs]);
ax = gca; ax.TickDir = 'out'; ax.FontSize=9; ax.XTick=[1:1:60]; xtickangle(ax,90);
LabelsOrdered_uw = Labels(I_uw);
lbls_uw = LabelsOrdered_uw(1:ftrnum);
lbls = regexprep(lbls_uw,'[_,:{}]','');
xticklabels(ax, lbls(1:ftrnum));
ylabel('Frequency'); box off
title(['Histogram of Features (simple count, ', num2str(ftrnum),' ftr model)']);
str = ['..\Figures\CorrelationTesting\', 'HistogramOfLinRegPnasalFeatures_SimpleCount',settings.experimentnumber]; % 

switch savefigas     
    case 'saveasTIFF'; print(fig, str, '-dtiff', '-r1000');     
    case 'saveasPNG'; saveas(fig, str, 'png');      
    case 'saveasFIG'; savefig(str);  
end
if closefigs; close(fig); end


%% Most important features using DM method
% Scores are based on 1/(x+1), x=rank per loop, higher scores are best.
% Features that regularly show up early (low rank e.g. First), get the biggest scores.
% Features that are occasionally very good, but sometimes late, are still considered.
% feature 3 is first 37 times, second 2 times. appears 39/41 times.
% while 15 never appears in top three, but is present 41 times.
% (side note, the variable 'If' is not my idea, this creation is SS...)
%nnz(ftr_array(:,1)==3)
%nnz(ftr_array_linear==3)
%Labels_Complete(3)

clearvars temp x If score score1 scoredata
score=zeros(1,size(Amatrix2_flow,2));
maxscore=0; %if large score is good
%maxscore = (size(labels_Step_Subj,2))^2; %if small score is good
itersize = size(labels_Step_Subj_flow,1);
try
    for i=1:itersize
        score1 = maxscore + 0*score;
        If = [];
        for j=1:ftrnum %size(labels_Step_Subj,2)
            temp = labels_Step_Subj_flow{i,j};
            x=sum(temp'==If,2);
            temp(x==1)=[];
            if ~isempty(temp)
                If(end+1)=temp;
            end
        end
        for j=1:length(If)
            score1(If(j)) = 1./((j-1)+1);
            %score1(If(j)) = (j-1).^0.5;
            %score1(If(j)) = 1./(j-1).^0.5;
        end
        score=score+score1;
    end
    scoredata = [score;1:length(score)]';
    % scoredata = sortrows(scoredata,'descend'); % descend command only works for table data
    scoredata = sortrows(scoredata,-1); % descending sort of col 1
    I_w=find(isnan(scoredata(:,1)));
    temp=scoredata(I_w,:);
    scoredata(I_w,:)=[];
    scoredata=[scoredata;temp];
    I_w = scoredata(:,2);
catch me
    disp(me.getReport);
end

clearvars If

% Modified histogram of ftrs selected during training
figure(34); clf(figure(34)); fig = gcf;
fig.Color = [1 1 1]; fig.Units = 'inches';
fig.Position = [2.2   2   12    4.5];
bar(scoredata(:,1), 'facealpha',0.2, 'edgealpha', 0.2);
xlim([0 PlotUptoNFtrs]);
ax = gca; ax.TickDir = 'out'; ax.FontSize=9; ax.XTick=[1:1:60]; xtickangle(ax,90);

LabelsOrdered_w = Labels(I_w);
lbls_w = LabelsOrdered_w(1:ftrnum);
lbls = regexprep(lbls_w,'[_,:{}]','');
xticklabels(ax,lbls(1:ftrnum));
ylabel('Weighted Score');
box off
title(['Histogram of Features (weighted count, ', num2str(ftrnum),' ftr model)']);
str = ['..\Figures\CorrelationTesting\', 'HistogramOfLinRegPnasalFeatures_WeightedCount',settings.experimentnumber]; %

switch savefigas     
    case 'saveasTIFF'; print(fig, str, '-dtiff', '-r1000');     
    case 'saveasPNG'; saveas(fig, str, 'png');      
    case 'saveasFIG'; savefig(str);  
end
if closefigs; close(fig); end

%% compare the two histogram ways of "top 50" ftrs
x_over = ismember(lbls_w(1:ftrnum), lbls_uw(1:ftrnum));
nnz(x_over) % this shows the number of features that occur in both lists 
    
%% Find best features through reverse stepwise regression
% starting with unique features from all 25 ftr models
ftr_array_2 = NaN(54, ftrnum);
for subj=1:54 % set pt num or all, no data for PT=1
    if ~ismember(subj, settings.Flow_list) %PT_list)
        continue
    end
   ftr_array_2(subj,:) = labels_Step_Subj_flow{subj,ftrnum};
end
% Start with the list of unique ftrs selected during training at set stop point
ftr_array_2 = ftr_array_2(settings.Flow_list,:);
uniqueftrs = unique(ftr_array_2);

% do backwards elimination
RemoveList = [];
while 1   
    [~,Pvals,~,~,~]=glmfitFast(Amatrix2_flow(:,uniqueftrs),GS,weights,1); 
    Pvals(1)=[];
    
    %remove least important
    [maxp,remi]=max(Pvals);

    disp(['Removing Ftr: ', num2str(uniqueftrs(remi)), ', p= ' num2str(maxp)]);
    RemoveList = [RemoveList; uniqueftrs(remi)];
    
    uniqueftrs(remi)=[];
   
    if isempty(uniqueftrs) || length(uniqueftrs)<1
        break
    end
end
% Remove list is currently in the order of first removed to last removed,
% so we really want those last removed, hence flip
I_r = fliplr(RemoveList');
LabelsOrdered_r = Labels(I_r);
lbls_r = LabelsOrdered_r(1:ftrnum);

%% Find best features through reverse stepwise regression
% starting with all features in all pts
ftr_array_3 = NaN(54, maxDataSize);
for subj=1:54 % set pt num or all, no data for PT=1
    if ~ismember(subj, settings.Flow_list) %PT_list)
        continue
    end
   ftr_array_3(subj,:) = labels_Step_Subj_flow{subj,maxDataSize};
end

% Start with the list of all unique ftrs selected during training from all pts
ftr_array_3 = ftr_array_3(settings.Flow_list,:);
uniqueftrs = unique(ftr_array_3);
%ftr_array_3_linear = ftr_array_3(:);
%searchftr = find(ftr_array_3_linear==126)

% Start with all features
uniqueftrs = [1:1:495];

% chuck out those that do not correlate between flow and pnasal
numofftrs = length(uniqueftrs);
Rsq_vals = NaN(1,numofftrs);
Rsq_vals_flip = NaN(1,numofftrs);
for ft=1:numofftrs
    % coefficient of variation, GOF
    % this is what we did in training.
     % [Rsq,~,~,~,~] = glmfitFast(Amatrix2_train(:,ft), Amatrix2_pnasal(~Isubj,ft), weights_train,1);
    [Rsq,~,~,~,~] = glmfitFast(Amatrix2_flow(:,ft), Amatrix2_pnasal(:,ft), weights,1);
    [Rsq2,~,~,~,~] = glmfitFast(Amatrix2_pnasal(:,ft), Amatrix2_flow(:,ft), weights,1);
    Rsq_vals(ft) = Rsq;
    Rsq_vals_flip(ft) = Rsq2;
end
% find bad correlations
Bad_Rsq_vals_i = Rsq_vals<=0.5;
NumBadFtrsRsq = nnz(Bad_Rsq_vals_i);
str = ['Number of bad ftrs by Rsq: ', num2str(NumBadFtrsRsq)]; disp(str);
% remove bad features
FtrsToExclude = find(Bad_Rsq_vals_i); %FtrsToExcl_withFlip = FtrsToExclude;
uniqueftrs(FtrsToExclude) = [];
if 0
Labels(uniqueftrs)
Labels(FtrsToExclude)
end
RemoveList = [];
labels_Step_Subj_flow_final = [];
beta_array_flow_final=[];  %
Rsq_array_flow_final=[];
Err_array_flow_final=[];
predy_flow_final=[];
predy_pnasal_final=[];
while 1     % do backwards elimination
    
    [Rsq,Pvals,MSE_final,beta_final,predy_final_]=glmfitFast(Amatrix2_flow(:,uniqueftrs),GS,weights,1); 
    Pvals(1)=[]; predy_final_(predy_final_>maxG)=maxG; predy_final_(predy_final_<0)=0;
    
     %% store the output - flow
    labels_Step_Subj_flow_final{1,length(uniqueftrs)} = uniqueftrs;
    beta_array_flow_final{1,length(uniqueftrs)} = beta_final;
    %[Rsq_array_flow_final(1,length(uniqueftrs)),~,~,~,~] = glmfitFast(predy_final_,GS, weights,1);
    Rsq_array_flow_final(1,length(uniqueftrs)) = Rsq;
    Err_array_flow_final(1,length(uniqueftrs)) = nanmean(weights.*abs(predy_final_-GS));
    predy_flow_final(:,length(uniqueftrs)) = predy_final_;

    %% use beta from flow, apply to pnasal, save output
    predy_pnasal_ = [ones(size(Amatrix2_pnasal,1),1) Amatrix2_pnasal(:,uniqueftrs)]*beta_final;
    predy_pnasal_(predy_pnasal_>maxG)=maxG; predy_pnasal_(predy_pnasal_<0) = 0;
    predy_pnasal_final(:,length(uniqueftrs)) = predy_pnasal_;
    
    %remove least important
    [maxp,remi]=max(Pvals);

    disp(['Removing Ftr: ', num2str(uniqueftrs(remi)), ', p= ' num2str(maxp)]);
    RemoveList = [RemoveList; uniqueftrs(remi)];
    
    uniqueftrs(remi)=[];
   
    if isempty(uniqueftrs) || length(uniqueftrs)<1
        break
    end
end
% Remove list is currently in the order of first removed to last removed,
% so we really want those last removed, hence flip
I_r_full = fliplr(RemoveList');
LabelsOrdered_r_full = Labels(I_r_full);
lbls_r_full = LabelsOrdered_r_full(1:ftrnum);

% save after this run
if 0
savefilename = [settings.filename_flow(1:end-4),  '_FinalRun2.mat'];
str=['Saving ' savefilename]; disp(str);
try
save(savefilename); 
catch me
    disp(me.message);
end
end



%% compare the top 25 from the methods
I_uw_topN = (I_uw(1:ftrnum));
I_w_topN = (I_w(1:ftrnum));
I_r_topN = (I_r(1:ftrnum));
I_r_full_topN = (I_r_full(1:ftrnum));

% use set logic for finding the common ftrs 
% note that these tables are in no particular order, just sorted by number
I_c1 = intersect(I_uw_topN', I_w_topN);
I_common = intersect(I_c1, I_r_topN');
LabelsOrdered_c = Labels(I_common);
common_features = table(I_common, LabelsOrdered_c, ...
        'VariableNames', {'Num', 'Ftr'});   
I_c2 = union(I_uw_topN', I_w_topN);
I_common2 = union(I_c1, I_r_topN');
LabelsOrdered_c2 = Labels(I_common2);
common_features2 = table(I_common2, LabelsOrdered_c2, ...
        'VariableNames', {'Num', 'Ftr'});

%% select which top 25 features to use for final model run
NfeaturesOpt=ftrnum;
if 0 % either of the histogram methods
    if 0 % unweighted histogram
        I = I_uw;
    else % weighted histogram
        I = I_w';
    end
    LabelsOrdered = Labels(I(1:NfeaturesOpt));
    TopFtrs = I(1:NfeaturesOpt);
else % else using backward elimination of unique ftrs at threshold
    LabelsOrdered = Labels(I_r(1:NfeaturesOpt));
    TopFtrs = I_r(1:NfeaturesOpt);
end

%% Find Top X unique features based on scores
if 0
    % if using transformed data, need to strip off the suffix
    if settings.TransformTheData; lengthsuf=6; else; lengthsuf=0; end
    MaxUniquefeatures=50;
    LabelsOrdered_NoSuffix=LabelsOrdered;
    for i=1:length(LabelsOrdered_NoSuffix)
        LabelsOrdered_NoSuffix{i}=LabelsOrdered_NoSuffix{i}(1:length(LabelsOrdered_NoSuffix{i})-lengthsuf);
    end
    LabelsOrdered_NoSuffixUnique=unique(LabelsOrdered_NoSuffix,'stable');
    LabelsOrdered_NoSuffixUnique(MaxUniquefeatures+1:end)=[];
    
    LabelsOrdered_NoSuffixUnique_sorted = sortrows(LabelsOrdered_NoSuffixUnique,1);
end

%% Final models using all data, using selected N optimal features (NfeaturesOpt)
% flow
Ftrs = TopFtrs;
Labels_ = [{'Intercept'};LabelsOrdered];
FtrVals_flow = Amatrix2_flow(:,Ftrs);
[Rsq,Pvals,RMSE,betas_flow,FinalMdlPredY_flow]=glmfitFast(FtrVals_flow,GS,weights,1); 
FinalMdlPredY_flow(FinalMdlPredY_flow>maxG)=maxG; % set upper limit
FinalMdlPredY_flow(FinalMdlPredY_flow<0) = 0; % set lower limit
[Rsq_flow, r_flow] = UnivariateStats(FinalMdlPredY_flow, GS, weights);
finalmdl_flow = fitglm(FtrVals_flow,GS,'weights',weights);
finalmdl_summary_table_flow = table([0;Ftrs'], Labels_, betas_flow, finalmdl_flow.Coefficients.SE(1:end), Pvals, ...
        'VariableNames', {'FtrNum', 'FtrName', 'Betas', 'SE', 'Pval'}); 
str=['Final model r=', num2str(r_flow),' ; Rsquared=', num2str(Rsq_flow)]; disp(str)
Summary_FinalMdl_flow = [r_flow,Rsq_flow];

% pnasal
FtrVals_pnasal = Amatrix2_pnasal(:,Ftrs);
[Rsq,Pvals,RMSE,betas_pnasal,FinalMdlPredY_pnasal]=glmfitFast(FtrVals_pnasal(BB_,:),GS(BB_),weights(BB_),1); 
FinalMdlPredY_pnasal(FinalMdlPredY_pnasal>maxG)=maxG; % set upper limit
FinalMdlPredY_pnasal(FinalMdlPredY_pnasal<0) = 0; % set lower limit
[Rsq_pnasal, r_pnasal] = UnivariateStats(FinalMdlPredY_pnasal, GS(BB_), weights(BB_));
finalmdl_pnasal = fitglm(FtrVals_pnasal(BB_,:),GS(BB_),'weights',weights(BB_));
finalmdl_summary_table_pnasal = table([0;Ftrs'], Labels_, betas_pnasal, finalmdl_pnasal.Coefficients.SE(1:end), Pvals, ...
        'VariableNames', {'FtrNum', 'FtrName', 'Betas', 'SE', 'Pval'}); 
    
% ToDo: correct calc of FinalMdlPredY_pnasal, using pnasal features, AND flow betas.
FinalMdlPredY_pnasal_wFlowBeta = [ones(nnz(BB_),1) FtrVals_pnasal(BB_,:)]*betas_flow;
FinalMdlPredY_pnasal_wFlowBeta(FinalMdlPredY_pnasal_wFlowBeta>maxG)=maxG; % set upper limit
FinalMdlPredY_pnasal_wFlowBeta(FinalMdlPredY_pnasal_wFlowBeta<0) = 0; % set lower limit
[Rsq_pnasal_wFlowBeta, r_pnasal_wFlowBeta] = UnivariateStats(FinalMdlPredY_pnasal_wFlowBeta, GS(BB_), weights(BB_));

str=['Final model r=', num2str(r_pnasal),' ; Rsquared=', num2str(Rsq_pnasal)]; disp(str)
str=['Final model r=', num2str(r_pnasal_wFlowBeta),' ; Rsquared=', num2str(Rsq_pnasal_wFlowBeta)]; disp(str)

Summary_FinalMdl_pnasal = [r_pnasal,Rsq_pnasal];

%% Output summary data
maxDataSize
settings.experimentnumber
SummaryStats_FF = [SummaryStats_Fig4, SummaryStats_FigE2, Summary_FinalMdl_flow];
SummaryStats_FP = [SummaryStats_Fig5, SummaryStats_FigE3, Summary_FinalMdl_pnasal];
SummaryStats = [SummaryStats_FF;SummaryStats_FP];

SummaryStatsFinal_25 = [Summary_FinalMdl_flow, Summary_FinalMdl_pnasal];
SummaryStatsFinalT_25 = array2table(SummaryStatsFinal_25, 'VariableNames', {'Flow_r','Flow_Rsq','Pnasal_r','Pnasal_Rsq'})

finalmdl_summary_table_flow
finalmdl_summary_table_pnasal % this should use same betas as flow... fix!

%% reduced N feature model
Ftrs = I_common2;
% step 1 - find top 5 ftrs univariate ftrs from top 25-ish list
RsqVsFD = NaN(length(Ftrs),1);
RsqVsFD_pnasal = NaN(length(Ftrs),1);
FtrVsFtrPnasal = NaN(length(Ftrs),1);
for ft=1:length(Ftrs)   
    % Rsq vs flow:drive
    FtrVal = Amatrix2_flow(:,Ftrs(ft));
    [Rsq,~,~,~]=glmfitFast(FtrVal,GS,weights,1);
    RsqVsFD(ft) = Rsq;
    
    % Rsq in pnasal vs flow:drive
    FtrVal_pnasal = Amatrix2_pnasal(BB_,Ftrs(ft));
    [Rsq,~,~,~]=glmfitFast(FtrVal_pnasal,GS(BB_),weights(BB_),1);
    RsqVsFD_pnasal(ft) = Rsq;
    
    % Rsq, Ftr(pneumo) vs Ftr(pnasal)
    [Rsq,~,~,~]=glmfitFast(FtrVal(BB_),FtrVal_pnasal,weights(BB_),1);
    FtrVsFtrPnasal(ft) = Rsq;
end
[ia, ib] = sort(RsqVsFD,'descend');
SortedFtrs = Ftrs(ib);
Labels_sorted = Labels(SortedFtrs);
Labels_sorted(1:10)
RsqVsFD_sorted = RsqVsFD(ib);
RsqVsFD_pnasal_sorted = RsqVsFD_pnasal(ib);

T2 = table(Ftrs(ib), Labels_sorted, ...
        'VariableNames', {'Num', 'Ftr'});

% figure
if 0
NftrsToPlot = 27;
facealpha = 0.05; % was 0.08
facecolor = [0.2 0.2 0.2]; %was [0.1 0.1 0.1]
for ft=11:NftrsToPlot
    FtrVal = Amatrix2_flow(:,SortedFtrs(ft));
    FtrVal_pnasal = Amatrix2_pnasal(BB_,SortedFtrs(ft));

    figure(101); clf(figure(101)); fig = gcf;
    fig.Color = [1 1 1]; fig.Units = 'inches';
    fig.Position = [-19   5  16   4.5];
    
    subplot(1,3,1);
    set(gca,'box','off','tickdir','out','fontname','arial narrow');
    scatter(FtrVal,GS,2,'filled','markerfacealpha', facealpha, 'markerfacecolor', facecolor); hold on;
    xlabel('Pneumotach feature value'); ylabel('{\itflow:drive}');
    str=['R^2=',num2str(RsqVsFD_sorted(ft))]; title(str); axis square
    
    subplot(1,3,2);
    scatter(FtrVal_pnasal,FtrVal(BB_),2,'filled','markerfacealpha', facealpha, 'markerfacecolor', facecolor); hold on;
    xlabel('Pnasal feature value)'); ylabel('Pneumotach feature value');   
    str=['R^2=',num2str(FtrVsFtrPnasal(ft))]; title(str); axis square
    
    subplot(1,3,3);
    scatter(FtrVal_pnasal,GS(BB_),2,'filled','markerfacealpha', facealpha, 'markerfacecolor', facecolor); hold on;
    xlabel('Pnasal feature value'); ylabel('{\itflow:drive}');
    str=['R^2=',num2str(RsqVsFD_pnasal_sorted(ft))]; title(str); axis square
    
    title_str = ['Ftr ', num2str(ft), ' -- ', Labels_sorted{ft}];
    title_str = regexprep(title_str,'[_,:{\\}]',' ');
    suptitle(title_str);
    save_str = ['..\Figures\Univariate\TopNftrs_', title_str]; %
    saveas(fig, save_str, 'png');
end
end

%% DLM final model
% pseudo-arbitrarily set "good" ftrs
if 0
goodFtrs = [204, 198, 230, 5, 201, 166];
goodFtrs = [204, 230, 5, 201, 166];
goodFtrs = [204, 230, 5, 201];
goodFtrs = [204, 230, 5, 166]; % 166 (TiTtot) doesn't add anything
goodFtrs = [204, 230, 5]; % good three
end
% pick one
%goodFtrs = (T2.Num(1:5))'; % the top five univariate
goodFtrs = I_r(1:5); % the top five ftrs from reverse elimination

% Final flow model using all data, uses selected N optimal features (NfeaturesOpt)
Ftrs = goodFtrs; % Ftrs = TopFtrs; % Ftrs from 25 model above
LabelsOrdered = Labels(Ftrs);
Labels_ = [{'Intercept'};LabelsOrdered];
FtrVals_flow = Amatrix2_flow(:,Ftrs);
[Rsq,Pvals,RMSE,betas_flow,FinalMdlPredY_flow]=glmfitFast(FtrVals_flow,GS,weights,1); 
FinalMdlPredY_flow(FinalMdlPredY_flow>maxG)=maxG; % set upper limit
FinalMdlPredY_flow(FinalMdlPredY_flow<0) = 0; % set lower limit
[Rsq_flow, r_flow] = UnivariateStats(FinalMdlPredY_flow, GS, weights);
finalmdl_flow = fitglm(FtrVals_flow,GS,'weights',weights);
finalmdl_summary_table_flow = table([0;Ftrs'], Labels_, betas_flow, finalmdl_flow.Coefficients.SE(1:end), Pvals, ...
        'VariableNames', {'FtrNum', 'FtrName', 'Betas', 'SE', 'Pval'}); 

% Final pnasal model using all data, uses same selected N optimal features (NfeaturesOpt)
FtrVals_pnasal = Amatrix2_pnasal(:,Ftrs);
[Rsq,Pvals,RMSE,betas_pnasal,FinalMdlPredY_pnasal]=glmfitFast(FtrVals_pnasal(BB_,:),GS(BB_),weights(BB_),1); 
FinalMdlPredY_pnasal(FinalMdlPredY_pnasal>maxG)=maxG; % set upper limit
FinalMdlPredY_pnasal(FinalMdlPredY_pnasal<0) = 0; % set lower limit
[Rsq_pnasal, r_pnasal] = UnivariateStats(FinalMdlPredY_pnasal, GS(BB_), weights(BB_));
finalmdl_pnasal = fitglm(FtrVals_pnasal(BB_,:),GS(BB_),'weights',weights(BB_));
finalmdl_summary_table_pnasal = table([0;Ftrs'], Labels_, betas_pnasal, finalmdl_pnasal.Coefficients.SE(1:end), Pvals, ...
        'VariableNames', {'FtrNum', 'FtrName', 'Betas', 'SE', 'Pval'}); 

Summary_FinalMdl_flow_dlm = [r_flow,Rsq_flow];
Summary_FinalMdl_pnasal_dlm = [r_pnasal,Rsq_pnasal];
SummaryStatsFinal_5 = [Summary_FinalMdl_flow_dlm, Summary_FinalMdl_pnasal_dlm];
SummaryStatsFinalT_5 = array2table(SummaryStatsFinal_5, 'VariableNames', {'Flow_r','Flow_Rsq','Pnasal_r','Pnasal_Rsq'})

%% testing 3d scatterplot
if 0
if length(Ftrs)<4
    criteriaR = GS>=0.7;
    PredR = FinalMdlPredY_flow>=0.7;
    Ilist_mode = Ftrs;
    C = NaN(length(criteriaR),3);
    err = NaN(length(criteriaR),1);
    
    % this is a terribly inefficient method to set colours for criteria
    for k = 1:length(criteriaR)
        if criteriaR(k)==1
            C(k,1:3)=[1,0,0];
        else
            C(k,:)=[0,0,1];
        end
        % also set a flag if we failed to predict this one
        err(k) = (criteriaR(k)~=PredR(k));
    end
    
    % restrict the size down a little, a bit more manageable to view
    subset = 1:50:size(Amatrix2_flow,1);
    C = C(subset,:); err = err(subset,:);
    S = ones(length(C),1)*4;
    X = Amatrix2_flow(subset,Ilist_mode(1));
    Y = Amatrix2_flow(subset,Ilist_mode(2));
    Z = Amatrix2_flow(subset,Ilist_mode(3));
    
    figure(210); clf(figure(210));
    h = scatter3(X, Y, Z, S, C, 'filled'); hold on;
    xlabel([Labels(Ilist_mode(1))]);
    ylabel([Labels(Ilist_mode(2))]);
    zlabel([Labels(Ilist_mode(3))]);
    % these limits will almost certainly need changing for different Ilist
    %xlim([0 1.2]);
    %zlim([0 0.0025]);
    
    % mark the misclassifications
    scatter3(X(logical(err)), Y(logical(err)), Z(logical(err)));
    box on
end   
end

%% Table E3 and E5 - the real final model
NumFtrs = 25; % set as 25 or 5 or 4 or 3
FinalFtrs = labels_Step_Subj_flow_final{NumFtrs}; % note the order here isn't most important.. it's just what's left.
FinalBetas = beta_array_flow_final{NumFtrs};

% re-order the final model betas so that they are in correct elimination order
%I_r_full % is the elimination order
%FinalFtrs % is the order of betas at the set feature count
%[4, 3, 5, 2, 1]

[check,index] = ismember(I_r_full(1:length(FinalFtrs))', FinalFtrs);
if nnz(check) == length(check) % check must be all ones, otherwise there is some error
    % the first beta is the intercept term, keep this aside
    intercept = FinalBetas(1); FinalBetas_ = FinalBetas(2:end);
    FinalBetas_sorted = FinalBetas_(index);
    FinalBetas_sorted = [intercept; FinalBetas_sorted];
    FinalFtrs_sorted = FinalFtrs(index);
else
    disp('failed initial check');
end

LabelsFinal_ = Labels(I_r_full(1:NumFtrs));
LabelsFinal = [{'Intercept'};LabelsFinal_];

FtrVals_flow = Amatrix2_flow(:,FinalFtrs_sorted);
FinalPredy_flow = [ones(nnz(BB),1) FtrVals_flow]*FinalBetas_sorted;
FinalPredy_flow(FinalPredy_flow>maxG)=maxG; % set upper limit
FinalPredy_flow(FinalPredy_flow<0) = 0; % set lower limit
[Rsq_flow, r_flow] = UnivariateStats(FinalPredy_flow, GS, weights);

FtrVals_pnasal = Amatrix2_pnasal(BB_,FinalFtrs_sorted);
FinalPredy_pnasal = [ones(nnz(BB_),1) FtrVals_pnasal]*FinalBetas_sorted;
FinalPredy_pnasal(FinalPredy_pnasal>maxG)=maxG; % set upper limit
FinalPredy_pnasal(FinalPredy_pnasal<0) = 0; % set lower limit
[Rsq_pnasal, r_pnasal] = UnivariateStats(FinalPredy_pnasal, GS(BB_), weights(BB_));

% calc univariate performance (bivariate)
% for each ftr, get r and Rsq for:
%   ftr_flow vs GS flow:drive
%   ftr_flow vs ftr_pnasal
%   bias (median pnasal / median flow)

FinalFtrs_sorted=I_r_full(1:NumFtrs)';

r_FtrFlowVsFD = NaN(length(FinalFtrs_sorted),1);
Rsq_FtrFlowVsFD = NaN(length(FinalFtrs_sorted),1);
r_FtrFlowVsFtrPneumo = NaN(length(FinalFtrs_sorted),1);
Rsq_FtrFlowVsFtrPneumo = NaN(length(FinalFtrs_sorted),1);
bias =  NaN(length(FinalFtrs_sorted),1);

for ft=1:length(FinalFtrs_sorted)   
    % FtrFlow vs flow:drive
    FtrVal_pneumo = Amatrix2_flow(:,FinalFtrs_sorted(ft)); 
    [Rsq_FtrFlowVsFD(ft),~,~,~]=glmfitFast(FtrVal_pneumo,GS,weights,1); 
    R_ = weightedcorrs([FtrVal_pneumo,GS],weights);
    r_FtrFlowVsFD(ft) = R_(1,2);
    
    % FtrFlow vs FtrPnasal
    FtrVal_pnasal = Amatrix2_pnasal(BB_,FinalFtrs_sorted(ft)); 
    %[Rsq_FtrFlowVsFtrPneumo(ft),~,~,~]=glmfitFast(FtrVal_pnasal,FtrVal_pneumo(BB_),weights(BB_),1);
    [Rsq_FtrFlowVsFtrPneumo(ft),~,~,~]=glmfitFast(FtrVal_pneumo(BB_),FtrVal_pnasal,weights_pnasal,1);
    %[Rsquared_w(ft), R_w(ft)] = UnivariateStats(FtrVal_pneumo(BB_),FtrVal_pnasal,weights(BB_));
    %[Rsquared_uw(ft), R_uw(ft)] = UnivariateStats(FtrVal_pneumo(BB_),FtrVal_pnasal);
    R_ = weightedcorrs([FtrVal_pnasal,FtrVal_pneumo(BB_)],weights(BB_));
    r_FtrFlowVsFtrPneumo(ft) = R_(1,2);
    
    % also calculate bias, as (median value Pneumotach) / (median value Pnasal)
    bias(ft) = median(FtrVal_pnasal) ./ median(FtrVal_pneumo);
end

% do this only to get stats (SEM) on betas in the final model
[~,~,stats] = glmfit(Amatrix2_flow(:,FinalFtrs_sorted),GS,'normal','weights',weights);
% the betas here, and the betas doing the reverse elimination should be the same
if isequal(round(stats.beta,4), round(FinalBetas_sorted,4)) % if not equal, something is wrong
    SEM = stats.se;
    Pvals = stats.p;
    Feature = ones(length(FinalBetas_sorted),1);
    FinalModel = [FinalBetas_sorted, SEM, [NaN;r_FtrFlowVsFD],[NaN;Rsq_FtrFlowVsFD], [NaN;r_FtrFlowVsFtrPneumo], [NaN;Rsq_FtrFlowVsFtrPneumo], [NaN;bias]];
    FinalModelTable = array2table(...
        [Feature, FinalBetas_sorted, SEM, Pvals, ...
        [NaN;r_FtrFlowVsFD], [NaN;Rsq_FtrFlowVsFD], ...
        [NaN;r_FtrFlowVsFtrPneumo],...
        [NaN;Rsq_FtrFlowVsFtrPneumo], ...
        [NaN;bias]], ...
        'VariableNames', {'Feature', 'Betas','SEM','Pvals', ...
        'RFtrFlowVsFD','RsqFtrFlowVsFD',...
        'RFtrFlowVsFtrPnasal','RsqFtrFlowVsFtrPnasal',...
        'bias'});
    FinalModelTable.Feature = LabelsFinal;
else
    disp('failed to match betas');
end


%% Final 25 feature Model
if 0
FinalFtrs_sorted
LabelsFinal
FinalBetas_sorted
end
save('FinalModel_FtrsAndBetas.mat', 'FinalFtrs_sorted','LabelsFinal','FinalBetas_sorted');
save('FinalModel_Table.mat', 'FinalModelTable');

%% Figure E4 and E5 - what does the final model look like? 
% and how well does it perform?

FlowModel = 0; % Set one for flow, zero for pnasal

if FlowModel % Flow
    %FinalMdlPredY = FinalMdlPredY_flow;
    %FinalMdlPredY = predy_flow_final(:,ftrnum);
    FinalMdlPredY = FinalPredy_flow; 
    GS_ = GS;
    predy = predy_flow;
    weights_ = weights;
    Gtest_avg = Gtest_flow_avg;
    plotstr = 'Figure_E4_Flow';
    FlowModel = 1; % switch
else % Pnasal
    %FinalMdlPredY = FinalMdlPredY_pnasal;
    %FinalMdlPredY = predy_pnasal_final(BB_,ftrnum);
    FinalMdlPredY = FinalPredy_pnasal; 
    GS_ = GS(BB_);
    predy = predy_pnasal(BB_);
    weights_ = weights(BB_);
    Gtest_avg = Gtest_pnasal_avg;
    plotstr = 'Figure_E5_Pnasal';
end

% get the cross validated performance (unweighted, values presented in text are weighted, and better))
if FlowModel
    CrossValPerf = RsqL1O_flow(NumFtrs)
    %CrossValPerf = RsqL1O_flow(25)
else
    CrossValPerf = RsqL1O_pnasal(NumFtrs)
    %CrossValPerf = RsqL1O_pnasal(25)
end

figure(102); clf(figure(102)); fig = gcf;
fig.Color = [1 1 1]; fig.Units = 'inches';
fig.Position = [-19.5 -0.5  10  10];
%fig.Position = [0.5 0.5  10  10];

% full model vs reduced ftr model
subplot(2,2,1);
set(gca,'box','off','tickdir','out','fontname','arial narrow');
scatter(FinalMdlPredY.*100,predy.*100,2,'filled','markerfacealpha', facealpha, 'markerfacecolor', facecolor); hold on;
hndl_ref = refline(1,0); hndl_ref.LineStyle = '--'; hndl_ref.Color = [0.25 0.25 0.25];
ax = gca; set(ax,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz);
ylim([-5 150]); xlim([-5 150]); ax.YTick=[0:25:125]; ax.XTick=[0:25:125];
xlabel('Reduced feature model {\itflow:drive}', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
ylabel('Full feature model {\itflow:drive}', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
[Rsq, r] = UnivariateStats(FinalMdlPredY, predy, weights_);
axis square; text(90, 20, ['R^2=',num2str(Rsq,sigdig)]); %text(90, 10, ['r=',num2str(r,2)]);

% confusion matrix
subplot(2,2,2);
customcmap = GetCustomColorMap('gray'); % 'SS'
g = NaN*GS_;
ghat = NaN*GS_;
classes = 5;
switch classes
    case 5
        classcutoffs = [0.9 0.7 0.5 0.3];
        labeltext={'Normal','Mild','Moderate','Severe','V.Severe'};
    case 4
        classcutoffs = [0.9 0.5 0.3];
        labeltext={'Normal','Mild','Moderate','Severe'};
    case 3
        classcutoffs = [0.9 0.5];
        labeltext={'Normal','Intermediate','Flow Limited'};
    otherwise
        classcutoffs = [0.7]; %normal mild moderate severe v.severe / normal borderline mild moderate severe
        labeltext={'Normal','Flow Limited'};
end
Nclasses=length(classcutoffs)+1;
g(GS_>classcutoffs(1))=1;
ghat(FinalMdlPredY>classcutoffs(1))=1;
for i=2:Nclasses
    g(GS_<classcutoffs(i-1))=i;
    ghat(FinalMdlPredY<classcutoffs(i-1))=i;
end
[C,order] = confusionmat(g,ghat);
%sumactual=sum(C')';
sumactual=sum(C,2);
sumestimated=sum(C);
%C_Factual=C./sum(C')'
C_Factual=C./sum(C,2)*100; %rows are actual, cols are estimated
C_Festimated=C./sum(C)*100;
C_Total = (C./sum(C(:)))*100;
if 1
    AccA_C_Factual = C_Factual.*diag(ones(1,length(C)));
    AccB_C_Factual = C_Factual.*(diag(ones(1,length(C)-1),1)+diag(ones(1,length(C)-1),-1));
    AccAactual = mean(sum(AccA_C_Factual,2));
    AccBactual = mean(sum(AccA_C_Factual,2)) + mean(sum(AccB_C_Factual,2));
    AccA_C_Festimated = C_Festimated.*diag(ones(1,length(C)));
    AccB_C_Festimated = C_Festimated.*(diag(ones(1,length(C)-1),1)+diag(ones(1,length(C)-1),-1));
    AccAestimated = mean(sum(AccA_C_Festimated));
    AccBestimated = mean(sum(AccA_C_Festimated)) + mean(sum(AccB_C_Festimated));
    ACCs = [AccAactual AccBactual; AccAestimated AccBestimated]
    % first row is actual, second row is estimate - use second row
    %second col gives accuracy if accepting next-category error as correct
end
x = order'; y = order';
if 1; C1 = C_Festimated; else; C1 = C_Factual; end
%C1 = C_Total;
xSplit = diff(x)/2;                 % Find edge points
ySplit = diff(y)/2;
xEdges = [x(1)-xSplit(1) x(2:end)-xSplit x(end)+xSplit(end)];
yEdges = [y(1)-ySplit(1) y(2:end)-ySplit y(end)+ySplit(end)];
[XGrid, YGrid] = meshgrid(xEdges,yEdges);
YGrid = flipud(YGrid);              % To match expected behavior
XGrid = fliplr(XGrid);
C2 = [[C1 zeros(size(C1,1),1)] ; zeros(1,size(C1,2)+1)];% Last row/col ignored
pcolor(XGrid,YGrid,(1-(1-C2/100).^2)*100)
hold on                             % Plot original data points
[X,Y] = meshgrid(x,y);
set(gcf,'colormap',customcmap);
C1 = flipud(C1);
C1 = fliplr(C1);
C1 = C1';
for i=1:size(C1,1)
    for j=1:size(C1,2)
        if C1(i,j)<((max(max(C1))-min(min(C1)))/2+min(min(C1)))
            textcolor=[1 1 1];
        else
            textcolor=[0 0 0];
        end
        text(x(i),y(j),num2str(round(C1(i,j))),'color',textcolor,'horizontalalignment','center','fontname','arial narrow')
    end
end
ax = gca; set(ax,'fontname','arial narrow','FontSize', TickFntSz);
yticks(y); yticklabels(gca,fliplr(labeltext));
xticks(x); xticklabels(gca,fliplr(labeltext));
xlabel('Reduced feature model classification', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
ylabel('Gold standard classification', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
axis square; %title([num2str(AccAestimated,3),'    ',num2str(AccBestimated,3)]);

% GS vs reduced ftr model
% add box plot over
subplot(2,2,3); hold on; set(gca,'box','off','tickdir','out','fontname','arial narrow');
scatter(FinalMdlPredY.*100,GS_.*100,2,'filled','markerfacealpha', facealpha, 'markerfacecolor', facecolor); hold on;
dx=0.2; xbins=[0 0.3:dx:0.9 1.5];
clear medianX medianY upperIQRY lowerIQRY upperIQRY2 lowerIQRY2
for i=1:length(xbins)-1
    Ix=FinalMdlPredY>xbins(i)&FinalMdlPredY<xbins(i+1);
    medianX(i)=prctile(FinalMdlPredY(Ix),50);
    medianY(i)=prctile(GS_(Ix),50);
    upperIQRY(i)=prctile(GS_(Ix),75);
    lowerIQRY(i)=prctile(GS_(Ix),25);
    upperIQRY2(i)=prctile(GS_(Ix),90);
    lowerIQRY2(i)=prctile(GS_(Ix),10);
end
hold on;
plotbarwitherrorsmedian(100*medianY,100*medianX,100*upperIQRY2,100*lowerIQRY2,zeros(length(medianY),3),NaN,100*0.005,100*0.01);
plotbarwitherrorsmedian(100*medianY,100*medianX,100*upperIQRY,100*lowerIQRY,zeros(length(medianY),3),NaN,100*dx/2,100*0.01);
ax = gca; set(ax,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz);
ylim([-5 150]); xlim([-5 150]); ax.YTick=[0:25:125]; ax.XTick=[0:25:125];
xlabel('Reduced feature model {\itflow:drive}', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
ylabel('Gold Standard {\itflow:drive}', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
Rsq = glmfitFast(FinalMdlPredY,GS_,weights_,1);
%text(100, 20, ['R^2=',num2str(Rsq, sigdig)]); axis square;

% patient median - sleep, but not incl apnea and low flow breaths
subplot(2,2,4);
FinalMdlPtAvg = NaN(54,1); 
if FlowModel
    for subj=1:54 % subj=11
        if ~ismember(subj, settings.Flow_list) %PT_list)
            continue
        end
        Isubj=(PtData_flow.PT==subj)&(PtData_flow.Hypnog<4)&(PtData_flow.Ar==0); % find the sleep only breaths that belong to this pt
        PredY_pt = FinalMdlPredY(Isubj);
        % add apneoa breaths as 0.1 to PredY_pt
        indAP = find(RemovedBB_flow.RemovedBB_Apnoea.Pt == subj);
        indLF = find(RemovedBB_flow.RemovedBB_LowFlow.Pt == subj);
        numApBB = RemovedBB_flow.RemovedBB_Apnoea{indAP,2};
        if isempty(numApBB); numApBB=0; end
        numLFBB = RemovedBB_flow.RemovedBB_LowFlow{indLF,2};
        if isempty(numLFBB); numLFBB=0; end
        numRemBB = numApBB + numLFBB;
        PredY_pt = [PredY_pt;ones(numRemBB,1)*0.09];
        FinalMdlPtAvg(subj) = median(PredY_pt);
    end
    FinalMdlPtAvg = FinalMdlPtAvg(settings.Flow_list);
else % pnasal model
     for subj=1:54 % subj=11
        if ~ismember(subj, settings.Pnasal_list) %PT_list)
            continue
        end
        Isubj=(PtData_pnasal.PT==subj)&(PtData_pnasal.Hypnog<4)&(PtData_pnasal.Ar==0); % find the sleep only breaths that belong to this pt
        PredY_pt = FinalMdlPredY(Isubj(BB_));
        % add apneoa breaths as 0.1 to PredY_pt
        indAP = find(RemovedBB_pnasal.RemovedBB_Apnoea.Pt == subj);
        indLF = find(RemovedBB_pnasal.RemovedBB_LowFlow.Pt == subj);
        numApBB = RemovedBB_pnasal.RemovedBB_Apnoea{indAP,2};
        if isempty(numApBB); numApBB=0; end
        numLFBB = RemovedBB_pnasal.RemovedBB_LowFlow{indLF,2};
        if isempty(numLFBB); numLFBB=0; end
        numRemBB = numApBB + numLFBB;
        PredY_pt = [PredY_pt;ones(numRemBB,1)*0.09];
        FinalMdlPtAvg(subj) = median(PredY_pt);
     end
     FinalMdlPtAvg = FinalMdlPtAvg(settings.Pnasal_list);
end
[Rvalue,Pvalue,Slope,Intercept,Xmodel_sem]=plotregressionwithSEM(FinalMdlPtAvg.*100, Gtest_avg.*100);
ylim([-5 125]); xlim([-5 125]);
ax = gca; set(ax,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz);
ax.YTick=[0:25:125]; ax.XTick=[0:25:125];
ylabel({'Gold Standard {\itflow:drive}','median (%)'}, 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
xlabel({'Reduced feature model {\itflow:drive}', 'median (%)'}, 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow'); axis square
[Rsq,p,MSE,beta,ypred] = glmfitFast(Gtest_avg,FinalMdlPtAvg(:,1),ones(length(Gtest_avg),1),1);
%text(100, 20, ['R^2=',num2str(Rsq, sigdig)]);
text(0,110,[{'Sleep breaths only,','including apnea and','low flow breaths'}]);

% supstr = ['Flow model with reduced features, N=',num2str(length(Ftrs)) ];
% suptitle(supstr);

%savestr = ['..\Figures\', plotstr,' Model with N',num2str(length(Ftrs))];
savestr = ['..\Figures\', plotstr];
% Figure E4 and Figure E5
switch savefigas
    case 'saveasTIFF'; print(fig, savestr, '-dtiff', '-r1000');
    case 'saveasPNG'; saveas(fig, savestr, 'png');
    case 'saveasFIG'; savefig(savestr);
end
if closefigs; close(fig); end

%% Table E1 - univariate for all ftrs, not just final models
% calc univariate performance (bivariate)
% for each ftr, get r and Rsq for:
%   ftr_flow vs GS flow:drive
%   ftr_pnasal vs GS flow:drive
%   ftr_flow vs ftr_pnasal

SuppFtrList = {
'Flattening (N=28)';...
'MIF_PIF';...
'Ali_InspFlat';...
'Ali_ExpFlat';...
'Ali_InspExpFlat';...
'AA_InspFlat9020';...
'AA_InspFlat8020';...
'AA_InspFlat7020';...
'AA_InspFlat6020';...
'AA_InspFlat5020';...
'AA_InspFlat90Ti';...
'AA_InspFlat80Ti';...
'AA_InspFlat70Ti';...
'AA_InspFlat60Ti';...
'AA_InspFlat50Ti';...
'AA_ExpFlat9020';...
'AA_ExpFlat8020';...
'AA_ExpFlat7020';...
'AA_ExpFlat6020';...
'AA_ExpFlat5020';...
'AA_ExpFlat90Te';...
'AA_ExpFlat80Te';...
'AA_ExpFlat70Te';...
'AA_ExpFlat60Te';...
'AA_ExpFlat50Te';...
'MIF50';...
'MEF50';...
'Teschler';...
'MPPW'; ...     % now called 'MostPromPeakW';...
'Scooping (deviation away from normal round contour) (N=8)';...
'AA_NED';...
'AreaUnderPeaksI';...   % was called 'SS_Area';...
'QuadI';...      % was called SinI
'QuadI50';...    % was called SinI50
'QuadE';...      % was called SinE
'InvParabI';...
'EllipseI';...
'HypcosI';...
'Asymmetry (N=9)';...
'AsymIndex';...
'SkewDistInsp';...
'KurtDistInsp';...
'SkewDistExp';...
'KurtDistExp';...
'AsymmetryInsp';...
'KurtDataInsp';...
'AsymmetryExp';...
'KurtDataExp';...
'Timing and volume ratio measures (N=28)';...
'Ti_Ttot';...
'Ti_Te';...
'TTran_i_Ti';...
'TTran_i_Ttot';...
'TTran_e_Te';...
'TTran_e_Ttot';...
'VTi_VTe';...
'VTi_VT';...
'VTe_VT';...
'VPEF_VTe ';... % 'KaplanIEvol';... %
'VPIF_VTi';...   % 'MorgensVPIFVTi';...
'InspVol_03Ti';... %'MorgensV03';...
'PIF_MIF';...
'FTi';...
'RTi';...
'DTi';...
'FTe';...
'RTe';...
'DTe';...
'PIF_PEF';...
'MIF_MEF';...
'MIF50_MEF50';...
'SeriesIEflow';...
'SeriesIEtime';...
'KaplanIEvol';...
'TpeakI_Ti';...
'TpeakE_Te';...
'TpeakI_TpeakE';...
'Fluttering (N=12)';...
'InspFlutPowOrig';...
'ExpFlutPowOrig';...
'InspExpFlutPowOrig';...
'InspExpFlutPowOrig_Sum';...
'InspFlutPow4to7';...
'ExpFlutPow4to7';...
'InspExpFlutPow4to7';...
'InspExpFlutPow4to7_Sum';...
'InspFlutPow8to12';...
'ExpFlutPow8to12';...
'InspExpFlutPow8to12';...
'InspExpFlutPow8to12_Sum';...
};

BivariateStats_r = [];
BivariateStats_Rsq = [];
FeatureName_list = [];
PneumoVsPnasal_Rsq = []; % just for interest
PnasalVsPneumo_Rsq = []; % just for interest
for ft = 1:length(SuppFtrList) % ft=2
    str = char(SuppFtrList(ft));
    % find the string in the full list of feature labels
    temp=find(startsWith(FeatureNames{:,2},str));
    if ~isempty(temp) % if we find something, 
        % expand it out to include transforms
        temp = [temp; temp+165; temp+330];
        
        % find Rsq for each
        Rsq_temp = NaN(length(temp),1);
        for foundftr=1:length(temp)
            [Rsq_temp(foundftr)] = glmfitFast(Amatrix2_flow(:,temp(foundftr)),GS, weights, 1);            
        end
        % get the best performing transform and timing variant
        [Flow_Rsq, ind] = max(Rsq_temp); TheOne = temp(ind);
        % get the corresponding r value
        r_ = weightedcorrs([Amatrix2_flow(:,TheOne), GS], weights); % use external fn
        r_temp = (r_(1,2));
        % get the value for pnasal for same ftr/timing/transform
        Pnasal_Rsq = glmfitFast(Amatrix2_pnasal(BB_,TheOne),GS(BB_), weights_pnasal, 1);  
        Bias = median(Amatrix2_pnasal(BB_,TheOne)) ./ median(Amatrix2_flow(:,TheOne));
        PneumoVsPnasal_Rsq_ = glmfitFast(Amatrix2_flow(BB_,TheOne),Amatrix2_pnasal(BB_,TheOne), weights_pnasal, 1);  
        PnasalVsPneumo_Rsq_ = glmfitFast(Amatrix2_pnasal(BB_,TheOne),Amatrix2_flow(BB_,TheOne), weights_pnasal, 1);
        
        % add to data to keep
        BivariateStats_r = [BivariateStats_r; r_temp];
        BivariateStats_Rsq = [BivariateStats_Rsq; [Flow_Rsq Pnasal_Rsq Bias]];
        PneumoVsPnasal_Rsq = [PneumoVsPnasal_Rsq; PneumoVsPnasal_Rsq_];
        PnasalVsPneumo_Rsq = [PnasalVsPneumo_Rsq; PnasalVsPneumo_Rsq_];
        
        % Num       Timing      Transform   Number
        %1 - 81     Original,   1           80
        %82 - 165   Transition, 1           83
        %166 - 246  Original,   0.5         80
        %247 - 330  Transition, 0,5         83
        %331 - 411  Original,   2           80
        %412 - 495  Transition, 2           83    
        %
        % shouldn't this be 81 original, 84 transition ?
        
        FeatureName_list = [FeatureName_list; Labels(TheOne)];
        
    else % if we don't find it, then it's a heading, just add nan to list
        BivariateStats_r = [BivariateStats_r; [NaN]];
        BivariateStats_Rsq = [BivariateStats_Rsq; [NaN NaN NaN]];
        FeatureName_list = [FeatureName_list; {''}];
        PneumoVsPnasal_Rsq = [PneumoVsPnasal_Rsq; NaN];
        PnasalVsPneumo_Rsq = [PnasalVsPneumo_Rsq; NaN];
    end
end

% compile a table to migrate into excel, for some further processing
Table_E1_ExportData = [BivariateStats_r BivariateStats_Rsq];
Table_E1 = array2table([[1:length(Table_E1_ExportData)]', Table_E1_ExportData],...
    'VariableNames', {'Feature','r','Rsq_flow','Rsq_pnasal','Bias'});
Table_E1.Feature = FeatureName_list;

Table_E1_Test = [BivariateStats_r, BivariateStats_Rsq, PneumoVsPnasal_Rsq, PnasalVsPneumo_Rsq];
Table_E1 = array2table([[1:length(Table_E1_Test)]', Table_E1_Test],...
    'VariableNames', {'Feature','r','Rsq_flow','Rsq_pnasal','Bias','PneumoVsPnasal','PnasalVsPneumo'});
Table_E1.Feature = FeatureName_list;


%% get the AHI data
[AHI_perPT, AHI_perPT_table] = getAHI_postanalysis();
if nnz(settings.Flow_list ~= find(~isnan(AHI_perPT(:,1)))) > 0; keyboard; end % check consistency
AHI_perPT_ = AHI_perPT(~isnan(AHI_perPT(:,1)),1);


%% AHI Scatter plot - Figure 1
% Figure 1, and calculate residual SD
figure(1);clf(figure(1)); fig = gcf;
fig.Color = [1 1 1]; fig.Units = 'inches';
fig.Position = [ 1 1  4.5  4.5];
% gold standard
[Rvalue,Pvalue,Slope,Intercept,Xmodel_sem]=plotregressionwithSEM(AHI_perPT_(:,1), Gtest_flow_avg.*100);
ylim([-5 125]); xlim([-5 105]);
ax = gca; set(ax,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz);
ax.YTick=[0:25:125]; ax.XTick=[0:20:100];
ylabel({'Gold Standard {\itflow:drive}','median (%)'}, 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
xlabel('Apnea-Hypopnea Index (evts/hour)', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow'); axis square
%[Rsq,p,MSE,beta,ypred] = glmfitFast(Gtest_flow_avg,AHI_perPT_(:,1),ones(length(Gtest_flow_avg),1),1);
%[Rsq_,p_,MSE_,beta_,ypred_] = glmfitFast(AHI_perPT_(:,1),Gtest_flow_avg,ones(length(Gtest_flow_avg),1),1);
%determine residual SD
linefit=Slope*AHI_perPT_(:,1)+Intercept;
absdiff = abs(AHI_perPT_(:,1)-linefit);
stdabsdiff = std(absdiff);
meanabsdiff = mean(absdiff);
str = ['..\Figures\Figure_1_ScatterHist\Figure_1',settings.experimentnumber];
switch savefigas
    case 'saveasTIFF'; print(fig, str, '-dtiff', '-r1000');
    case 'saveasPNG'; saveas(fig, str, 'png');
    case 'saveasFIG'; savefig(str);
end
if closefigs; close(fig); end

%% After adjusting for AHI
% Median
[b,dev,stats]=glmfit([AHI_perPT_(:,1) PredY_flow_avg],Gtest_flow_avg);
p_1_adj = stats.p(3);
r_1_adj = stats.coeffcorr(1,3);

%% Histograms
savefigas = ''; %'saveasTIFF';
doHistogramsPlots = 1;
if doHistogramsPlots
    if 1; edges = [0:0.1:1.1]; else edges = xbins; end
    % use data that has Apnea and LowFlow breaths included:
    % So, PredY_flow and Gtest_flow.
    
    %for subj=1:54 % set pt num or all, no data for PT=1,  % subj=2
    for subj=[10 21 22 34 54] 
        if ~ismember(subj, settings.Flow_list) %PT_list)
            continue
        end
        Isubj = PTarray_flow==subj;      % get the pt BB's
        Gtest_pt = Gtest_flow(Isubj);    % actual
        PredY_pt = PredY_flow(Isubj);    % pred
        Pt_Median = median(Gtest_pt);
        
        Gtest_pt(Gtest_pt>1.1) = 1.09;
        
        figure(300+subj); clf(figure(300+subj));fig = gcf;
        fig.Color = [1 1 1]; % set background colour to white
        fig.Units = 'inches';
        
        %LabelFntSz = 24;
        fig.Position = [19   0.5   3    2.5];
        h1 = histogram(Gtest_pt,edges, 'facealpha',1, 'edgealpha', 0, 'normalization', 'probability');%'pdf'
        hold on;
        currentYlim = ylim();
        plot([Pt_Median, Pt_Median], [0, currentYlim(2)],'k-', 'linewidth', 2);
        ax = gca;  set(ax,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz);
        ax.XTick=[0:0.25:1.25];
        xticklabels(ax, {'0','25','50','75','100','125'});
        xlabel('{\itflow:drive} (%)', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
        ylabel('Frequency', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow'); %yticks([]); % not using 'Relative probability'
        xlim([-0.05 1.15]);
        supstr = ['pt ',num2str(subj),', AHI ', num2str(round(AHI_perPT(subj,1))), ', median ' num2str(round(median(Gtest_pt*100)))];
        
        % save
        savestr = regexprep(supstr,'[_,:{}]','');
        str = ['..\Figures\Figure_1_ScatterHist\', savestr];
        switch savefigas
            case 'saveasTIFF'; print(fig, str, '-dtiff', '-r1000');
            case 'saveasPNG'; saveas(fig, str, 'png');
            case 'saveasFIG'; savefig(str);
        end
        if closefigs; close(fig); end
    end
end


%% GS vs Pred, GS vs FlowOnly, and GS vs DriveOnly

flow_only = PtData_flow.VE./PtData_flow.Veup;
drive_only = PtData_flow.DriveEdi./PtData_flow.Veup;

figure(300); clf(figure(300))
subplot(1,3,1)
scatter(predy_flow,Yval_flow, 2,'filled','markerfacealpha',0.06, 'markerfacecolor', [0 0 0]); hold on;
subplot(1,3,2)
scatter(flow_only,Yval_flow, 2,'filled','markerfacealpha',0.06, 'markerfacecolor', [0 0 0]); hold on;
subplot(1,3,3)
scatter(drive_only,Yval_flow, 2,'filled','markerfacealpha',0.06, 'markerfacecolor', [0 0 0]); hold on;

[rsq1,p1sq] = glmfitFast(predy_flow, Yval_flow,  weights,1);
[rsq2,p2sq] = glmfitFast(flow_only,  Yval_flow,  weights,1);
[rsq3,p3sq] = glmfitFast(drive_only, Yval_flow,  weights,1);


%% -----------------------  Sample traces ----------------------------



%% Load subject waveform data for VEVdrive, Flow and Edi plot
addNaNgaps=1;
SDU = 1;
% read spreadsheet
AnalyzeDataSpreadsheet = 'AnalyzeDataSpreadsheet.xlsx';
[num,patients,~] = xlsread(AnalyzeDataSpreadsheet,1,'F3:G56');
ftrnum = 25;
TickFntSz = 12;
LabelFntSz = 18;
FntSz = 18;
sigdig = '%.2f';
closefigs = 0;
snip = 10; %set snip number, snip 10 is for paper figures
%snip = 11; % for ATS 2018 talk
%snip = 12; % for ATS 2018 talk

%for n= [ 2 6 14 17 18 21 28 29 34 44 45 53 ] % full list
%for n= [ 14 29 34 45]   % list of those with at least 2 regions of interest
%for n= [ 14 34 45]      % list of those with at least 3 regions of interest
%for n = [ 14]            % list of those with at least 5 regions of interest
%for n = 5 % central
%for n = 5 % interesting
%for n = [2 5 14 21 45]
% for n = 45
%n = 10; % 10 is interesting too, a bit fl, but a bit central too.
for n =  45%[2 5 14 21 45] 
clearvars Edi Flow StarttimeSpike Data1
studyname = char(patients{n});
studynameNoExt = char(patients{n}(1:end-4));
load(['C:\PSG_Data\FlowDrive\SourceMat 20171123\',studyname], 'Edi', 'Flow', 'StarttimeSpike');
try
    load(['C:\PSG_Data\FlowDrive\Converted\',studynameNoExt,'_XHz.mat']);
catch GetEvtsFail
    disp(GetEvtsFail.getReport);
end
EventsAr = DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'EventsAr')==1));
EventsResp = DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'EventsResp')==1));
% arousal   = 1
% ap-O      = 2
% ap-C      = 3
% hyp-O     = 4
% mixed     = 5
% hyp-C     = 6

%ArColor
%RespOColor
%RespCColor

if ~(exist('StarttimeSpike', 'var') == 1)
    StarttimeSpike = 0;
end
Time = StarttimeSpike:0.008:StarttimeSpike+(length(Flow.values)-1)*0.008;
Isubj = PtData_flow.PT==n;
Data1 = [PtData_flow.BB_time(Isubj) PtData_flow.BB_Ttot(Isubj) predyL1O_array_flow(Isubj,ftrnum) GS(Isubj)];
cols = [NaN 1 2];
if addNaNgaps
    tol2=0.1;
    i=1;
    M=size(Data1,2);
    while i<(size(Data1,1)-1)
        if (Data1(i,cols(2))+Data1(i,cols(3))+tol2)<Data1(i+1,cols(2))
            Data1 = [Data1(1:i,:); NaN*ones(1,M); Data1((i+1):size(Data1,1),:)];
            %keyboard
            i=i+1;
        end
        i=i+1;
    end
end

if Flow.length<=length(Time)
    Flow.values(end:length(Time))=0;
else
    Flow.values(length(Time):end)=[];
end
if Edi.length<=length(Time)
    Edi.values(end:length(Time))=0;
else
    Edi.values(length(Time):end)=[];
end

% Data columns: (1) Time, (2) BBTtot, (3) PredY, (4) Gtest
% Clip Gtest to 1.05
Data1(:,4) = min(Data1(:,4), 1.05);

dsf=5; dt=Flow.interval;
FlowF=Flow.values;
if 1
    filter_HFcutoff_butter0 = 12.5;
    filter_order0 = 1;
    [B_butter0,A_butter0] = butter(filter_order0,[filter_HFcutoff_butter0]/(1/dt/2),'low');
    FlowF = filtfilt(B_butter0,A_butter0,FlowF); %filtfilt, otherwise flow signal is right-shifted
end

%% make the plot for the subject data
figure(200+n); clf(200+n); fig = gcf;
fig.Color = [1 1 1]; fig.Units = 'inches';
fig.Position = [24 2 12 7.5];  % position [ x0 y0 width height]
TickFntSz = 12;
threeplot = 0;

% Arousal
ax(1)=subplot(4,1,1); set(gca,'Position',[0.06 0.94 0.90 0.03]);
% turn continuous time Arousal events into list with duration
EvtsAr = downsample(EventsAr,dsf);
Timeds = downsample(Time,dsf);
I2 = find(diff(EvtsAr)==-1);
I1 = find(diff(EvtsAr)==1);
[I1,I2] = TidyStartEndEventList(I1,I2,length(Timeds));
It1 = Timeds(I1); It2 = Timeds(I2);
currentAxes = gca; % switch to current axes
% then mark these on the plot
for pn = 1:length(It1)
    BBmarker = patch([It1(pn) It2(pn) It2(pn) It1(pn)],...
        [0 0 1 1], [.8 .8 .8],...
        'FaceAlpha',0.3, 'EdgeColor','none', 'Parent',currentAxes);%
end
set(gca,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz);
set(gca,'xtick',[],'box','off');
set(gca,'xcolor',[1 1 1])
if n == 14 || n ==5
    set(gca,'ytick',[],'box','off');
    set(gca,'ycolor',[1 1 1])
else
    set(gca,'yticklabel',[]);
    set(gca,'ytick',[]);
    %ylabel('Ar', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
end
% flow:drive
ax(2)=subplot(4,1,2); set(gca,'Position',[0.06 0.64 0.90 0.27]);
if SDU
    stairs(Data1(:,1),100*Data1(:,4),'k', 'LineWidth', 1); hold('on'); % actual
    %stairs(Data1(:,1),100*Data1(:,3),'r', 'LineWidth', 2); % pred
else
    stairs(Data1(:,1),100*Data1(:,4),'k'); hold('on'); % actual
    stairs(Data1(:,1),100*Data1(:,3),'r'); % pred
end

plot([Data1(1,1) Data1(end,1)],[0 0],'k:');
plot([Data1(1,1) Data1(end,1)],100*[1 1],'k:');
plot([Data1(1,1) Data1(end,1)],50*[1 1],'k:');
set(gca,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz);
set(gca,'xtick',[],'box','off');
set(gca,'xcolor',[1 1 1])
if n == 14 || n ==5
    set(gca,'ytick',[],'box','off');
    set(gca,'ycolor',[1 1 1])
else
    %axc = gca; axc.YTick=[0:50:100]; yticklabels(axc, {'0', '50', '100'});
    axc = gca; axc.YTick=[0:25:100]; yticklabels(axc, {'0','25','50','75','100'});
    ylabel('{\itflow:drive} (%)', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
end
ylim([-5 105]);

% Flow (with resp event shading)
ax(3)=subplot(4,1,3); set(gca,'Position',[0.06 0.34 0.90 0.27]);
if SDU
    plot(downsample(Time,dsf),downsample(FlowF,dsf),'k', 'LineWidth', 1); hold on
else
    plot(downsample(Time,dsf),downsample(FlowF,dsf),'k'); hold on
end
plot(downsample(Time,dsf),zeros(length(downsample(FlowF,dsf)),1),'k:');
set(gca,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz);
set(gca,'xtick',[],'box','off');
set(gca,'xcolor',[1 1 1])
if n == 14 || n ==5
    set(gca,'ytick',[],'box','off');
    set(gca,'ycolor',[1 1 1])
else
    ylabel('Flow (L/s)', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
end
if 0; ylim([-0.8 0.8]); else; ylim([-1.1 1.1]); end
% turn continuous time Resp events into list with duration
EvtsResp = downsample(EventsResp,dsf);
Timeds = downsample(Time,dsf);
I2 = find(diff(EvtsResp)==-4);  % 4 is picking up hyp-O only, could expand to include more evt types
I1 = find(diff(EvtsResp)==4);
% enable this if doing central hypop
if n == 5
    I2 = find(diff(EvtsResp)==-6);  % 6 is hyp-C
    I1 = find(diff(EvtsResp)==6);
end
% enable this if doing obs ap
if n == 21
    I2 = find(diff(EvtsResp)==-2);  % 2 is Ap-O
    I1 = find(diff(EvtsResp)==2);
end
[I1,I2] = TidyStartEndEventList(I1,I2,length(Timeds));
It1 = Timeds(I1);  It2 = Timeds(I2);
currentAxes = gca; % switch to current axes
ylimits = currentAxes.YLim; % call this explicitly, fails otherwise
% then mark these on the plot
for pn = 1:length(It1)
    BBmarker = patch([It1(pn) It2(pn) It2(pn) It1(pn)],...
        [ylimits(1)+0.1 ylimits(1)+0.1 ylimits(2)-0.1 ylimits(2)-0.1], [.8 .8 .8],...
        'FaceAlpha',0.3, 'EdgeColor','none', 'Parent',currentAxes);%
end

% Edi
ax(4)=subplot(4,1,4); set(gca,'Position',[0.06 0.04 0.90 0.27]);
if SDU
    plot(downsample(Time,dsf),downsample(Edi.values,dsf),'k', 'LineWidth', 1); hold on;
else
    plot(downsample(Time,dsf),downsample(Edi.values,dsf),'k'); hold on;
end
set(gca,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz);
set(gca,'xtick',[],'box','off');
set(gca,'xcolor',[1 1 1])
axc = gca; 
if n == 14 || n ==5
    set(gca,'ytick',[],'box','off');
    set(gca,'ycolor',[1 1 1])
else
    axc.YTick=[0:15:60]; yticklabels(axc, {'0','15','30','45','60'});
    ylabel('EMGdi (uV)', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
    ylabel('Diaphragm EMG (\muV)', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
end
ylim([-2 62]);

linkaxes(ax,'x');

%% modify the plot
switch n
    case 30 % for TG, SDU2018
        switch snip
            case 10
                v1 = 7960; v2 = 8080;
                %v1 = 10985; v2 = v1+120;
        end
        useStart = 1;
    case 47 % for TG, SDU2018
        switch snip
            case 10
                v1 = 25330; v2 = 25400;
        end
        useStart = 1;
    case 2
        switch snip
            case 1
                v1 = 27000; v2 = 27350;
            case 10
                %v1 = 27025; v2 = 27145; % original export
                v1 = 27060; v2 = 27180; % shifted left a litte
        end
        
        useStart = 1;
    case 5 % central hypop and SHR
        switch snip
            case 1
                v1 = 21660; v2 = 21700; % central  % export two minute window, clip later
                v1 = 22950; v2 = 23150; % Full
                v1 = 22950; v2 = 22965; % 1 of 5
                v1 = 22980; v2 = 22995; % 2 of 5
                v1 = 23010; v2 = 23025; % 3 of 5
                v1 = 23055; v2 = 23070; % 4 of 5
                v1 = 23075; v2 = 23090; % 5 of 5
                useStart = 1;
            case 10
                v1 = 21655; v2 = 21775; % central  % export two minute window, clip later
                useStart = 1;
        end
    case 6
        v1 = 2900; v2 = 3600;
        useStart = 1;
    case 10
        v1 = 19500; v2 = 19620;
        useStart = 1;
    case 14 % off the scale, big flow and massive edi, still "flow limited"
        % multiple
        S14 = StarttimeSpike; StarttimeSpike = S14;
        if threeplot
            set(ax(2),'ylim', [-1.75 1.75]);
            set(ax(3),'ylim', [-2 60]);
        else
            set(ax(3),'ylim', [-1.75 1.75]);
            set(ax(4),'ylim', [-2 60]);
        end
        switch snip
            case 1
                v1 = 5200; v2 = 6400; % big window
                useStart = 1;
            case 2
                v1 = 85487; v2 = 85658; % subclinical events
                useStart = 0;
            case 3
                v1 = 86027; v2 = 86198; % stable FL
                useStart = 0;
            case 4
                v1 = 86385; v2 = 86556; % stable FL
                useStart = 0;
            case 5
                v1 = 23560;  v2 = 23680; % maybe no good
                useStart = 1;
            case 10
                %v1 = 86060; v2 = 86180; % stable FL
                v1 = 86090; v2 = 86210; % stable FL
                useStart = 0;
        end
    case 17
        v1 = 83450; v2 = 83700;
        useStart = 0;
    case 18
        v1 = 105600; v2 = 106000;
        useStart = 0;
    case 21
        switch snip
            case 1
                %v1 = 83900; %v2 = 85000; % not defined
                v1 = 84450; v2 = 84820;
                useStart = 0;
            case 10
                %v1 = 84670; v2 = 84790; % original
                v1 = 84730; v2 = 84850; % shifted left a little
                useStart = 0;
        end
    case 28
        v1 = 19460; v2 = 19580;
        useStart = 1;
    case 29
        % multiple
        S29 = StarttimeSpike; StarttimeSpike = S29;
        switch snip
            case 1
                v1 = 13400; v2 = 13600;
            case 2
                v1 = 15800; v2 = 16000;
        end
        useStart = 1;
    case 34
        % multiple
        switch snip
            case 1
                v1 = 90410; v2 = 90560;
            case 2
                v1 = 100000;  v2 = 100400;
            case 3
                v1 = 101705; v2 = 102100;
        end
        useStart = 0;
    case 44
        v1 = 96500; v2 = 97600;
        useStart = 0;
    case 45
        % multiple
        switch snip
            case 1
                v1 = 86450; v2 = 87000;
            case 2
                v1 = 91380; v2 = 91600;
            case 3
                v1 = 91450; v2 = 91515;
            case 10
                v1 = 91391; v2 = v1+120;  % 2 minutes
                %v1 = 91405; v2 = v1+120 %v2 = 91485; % for TG, SDU2018
            case 11
                v1 = 91451; v2 = 91511;  % 1 minute
            case 12
                v1 = 91480; v2 = 91490;  % 10 seconds, uses if 0 below, fig 2450
        end
        useStart = 0;
    case 53
        v1 = 103200;  v2 = 104200;
        useStart = 0;
end
if useStart
    xlim([StarttimeSpike+v1 StarttimeSpike+v2])
else
    xlim([v1 v2]);
end
% xlim([-inf inf]);

if 0
    figure(2450); clf(figure(2450)); fig = gcf;
    fig.Color = [1 1 1]; fig.Units = 'inches';
    fig.Position = [-12  1   10    7];
    plot(downsample(Time,dsf),downsample(FlowF,dsf),'k'); hold on
    plot(downsample(Time,dsf),zeros(length(downsample(FlowF,dsf)),1),'k:');
    set(gca,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz);
    set(gca,'xtick',[],'box','off');
    set(gca,'xcolor',[1 1 1])
    ylabel('Flow (L/s)', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
    ylim([-0.61 0.61]);    xlim([v1 v2]);
    str = ['..\Figures\SampleData\SampleData_Pt ', num2str(n), ', snip ', num2str(snip), 'singleBreath' ]; %
    print(fig, str, '-dtiff', '-r1000');
end

%% save the plot
str = ['..\Figures\SDU2018\SampleData_Pt ', num2str(n), ', snip ', num2str(snip),'nopred']; %
savefigas = 'saveasPNG';%'saveasTIFF';
switch savefigas
    case 'saveasTIFF'; print(fig, str, '-dtiff', '-r1000');
    case 'saveasPNG'; saveas(fig, str, 'png');
    case 'saveasFIG'; savefig(str);
end
if closefigs; close(fig); end

end
