%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
close; clear; clc
addpath('C:\Users\uqdmann\Dropbox\PUPbeta_git\');
cd('C:\Users\uqdmann\Dropbox\PUPbeta_git\FlowAndEdiPesQAO'); % code
datadir = 'C:\PSG_Data\FlowDrive\Analyzed\ExclBadR\'; % data
experimentnumber = '_n57'; 
load([datadir, 'PnasalDrive_only25Hz_ReNormalized_FeatureSpace_AutoRef2_Edi_Clean_LinRegWorkspace_TrainPnasal_TestPnasal', experimentnumber, '.mat']);

% probably should just load variables of interest to save time/space/electrons
% 'PtData','Amatrix2','FeatureNames','etc...'

%% Tidy up
clearvars -except ...
'Amatrix2_flow' ...
'Amatrix2_pnasal' ...
'Amatrix_flow' ...
'Amatrix_pnasal' ...
'Bad_Rsq_vals_i' ...
'beta_array_pnasal' ...
'ErrTrain_array_pnasal' ...
'experimentnumber' ...
'FeatureNames' ...
'flow_bb_to_match_pnasal' ...
'GS' ...
'Ind' ...
'Isubj' ...
'Labels' ...
'labels_flow' ...
'labels_Step_Subj_pnasal' ...
'LabelsOrderedOpt' ...
'maxDataSize' ...
'maxG' ...
'numofftrs' ...
'Original_FeatureNames' ...
'predyL1O_array_pnasal' ...
'PtData_flow' ...
'PtData_pnasal' ...
'RemovedFtrs' ...
'RsqTrain_array_pnasal' ...
'settings' ...
'weights' 

if settings.experimentnumber == '_n51' % some extra processing that was not originally done    
    % tidy up table labels
    PtData_Labels = PtData_pnasal.Properties.VariableNames';
    for n = 1:length(PtData_Labels)
        PtData_Labels{n}=regexprep(PtData_Labels{n,1},'_p','');
    end
    PtData_pnasal.Properties.VariableNames = PtData_Labels';
    PtData_flow.Properties.VariableNames = PtData_Labels';  
    load('BB_pnasaltoflow');
end

%% Use All breaths for classification performance plots
SleepOnly = 0;  % 0 for All breaths, 1 for sleep only breaths
IncludePnasal = 1;
[BB, BB_] = getAllorSleepOnlyBB(SleepOnly, IncludePnasal, PtData_flow, PtData_pnasal);
str = ['Flow - Using ',num2str(nnz(BB)),' of ', num2str(nnz(~isnan(PtData_flow.PT)))]; disp(str);
str = ['Pnasal - Using ',num2str(nnz(BB_)),' of ', num2str(nnz(~isnan(PtData_pnasal.PT)))]; disp(str);

% Select data to use
PlotUptoNFtrs = min([50 maxDataSize]); % normally 100

ftrnum=min([25 maxDataSize]); % normally 25
ErrorMaxN = size(predyL1O_array_pnasal,2); % normally 100

% pnasal
predy_pnasal = predyL1O_array_pnasal(BB_,ftrnum);
Yval_pnasal = GS(BB_);
weights_pnasal = weights(BB_);

% flow
% load from flow file
expnum = regexprep(experimentnumber,'[_n]','');
expnum_int = str2num(expnum)-50;
experimentnumber_flow = ['_n0', num2str(expnum_int)]; % use the current experiment number, less 50

datadir = 'C:\PSG_Data\FlowDrive\Analyzed\ExclBadR\'; % data
load([datadir, 'FlowDrive_only25Hz_ReNormalized_FeatureSpace_AutoRef2_Edi_Clean_LinRegWorkspace_TrainFlow_TestFlowAndPnasal', experimentnumber_flow, '.mat'], 'predyL1O_array_flow');
ftrnum_in_flow = min([ftrnum (size(predyL1O_array_flow,2))]); % as close as possible
predy_flow = predyL1O_array_flow(flow_bb_to_match_pnasal,ftrnum_in_flow);

savefigas =  'saveasPNG'; % ''; % 'saveasTIFF'; % options are saveasPNG, saveasFIG, and saveasTIFF
closefigs = 0;

TickFntSz = 12;
LabelFntSz = 18;
FntSz = 18;

%% Process test data. flow and pnasal, for all pts combined
for i=1:size(predyL1O_array_pnasal,2)% normally size(predyL1O_array,2) PlotUptoNFtrs   
    [RsqL1O_pnasal(i), RL1O_pnasal(i)] = ...
        UnivariateStats(predyL1O_array_pnasal(BB_,i),GS(BB_), weights(BB_)); % new method for R and Rsq
    ErrL1O_pnasal(i) = nanmean(weights(BB_).*abs(predyL1O_array_pnasal((BB_),i)-GS(BB_)));
end

%% For the average of all pts combined
figure(1); clf(figure(1));
subplot(1,2,1);
for pt = 1:54
    if ismember(pt, settings.Pnasal_list)
        stairs(RsqTrain_array_pnasal(pt,:), 'r'); hold on;
        stairs(ErrTrain_array_pnasal(pt,:), 'k');
    end
end
xlim([0, PlotUptoNFtrs]); ylim([0,1]);
title('Pnasal - Training'); legend('R squared', 'MAE');

subplot(1,2,2);
stairs(RL1O_pnasal,'b'); hold on; 
stairs(RsqL1O_pnasal,'r'); 
stairs(ErrL1O_pnasal, 'k');
xlim([0, PlotUptoNFtrs]); ylim([0,1]);
title('Pnasal - Testing'); legend('R','R squared','MAE');

%% Figure E3 - flow
% no Figure E3 for Train Pnasal - Test Pnasal
SummaryStats_FigE3 = [NaN, NaN, NaN, NaN];

%% Figure E4 - pnasal
figure(84); clf(figure(84)); fig = gcf;
fig.Color = [1 1 1]; fig.Units = 'inches';
fig.Position = [1  1   16   5];

% performance vs number of features
subplot(1,3,1); % pnasal
stairs([RL1O_pnasal(1:PlotUptoNFtrs);RsqL1O_pnasal(1:PlotUptoNFtrs);ErrL1O_pnasal(1:PlotUptoNFtrs)]'); hold on;
currentYlim=ylim();
plot([ftrnum, ftrnum], [0, currentYlim(2)],'k-', 'linewidth', 1.5);
ax = gca; set(ax,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz);
xlabel('Number of Features', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
ylabel('Performance (Testing)', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow'); 
ylim([0 1]); xlim([0 50]); axis square;

% proportion of total error
subplot(1,3,2);
a_out = mean(ErrTrain_array_pnasal(settings.Pnasal_list,1:ErrorMaxN),1);
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

str = ['..\Figures\CorrelationTesting\Figure_E4',settings.experimentnumber]; 

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

for subj=1:54
    if ismember(subj, settings.Pnasal_list) % add in pnasal RemBB
        Isubj=(PtData_pnasal.PT==subj); % find all the breaths that belong to this pt, this is just for counting
        numBBinTest_pnasal(subj,1) = nnz(Isubj); % the total number of breaths for this pt
        Isubj=(PtData_pnasal.PT==subj)&(PtData_pnasal.Hypnog<4)&(PtData_pnasal.Ar==0); % find the sleep only breaths that belong to this pt
        numBBinTest_pnasal(subj,2) = nnz(Isubj); % the number of sleep breaths for this pt
        Gtest_pt_pnasal = GS(Isubj);
        PredY_pt_pnasal = predyL1O_array_pnasal(Isubj,ftrnum);
        % add apneoa breaths as 0.1 to both Gtest_pt and PredY_pt
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
        Gtest_pt_flow = GS(Isubj);
        PredY_pt_flow = predyL1O_array_flow(Isubj,ftrnum_in_flow);
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

Gtest_flow_avg=Gtest_flow_avg(settings.Flow_list);
PredY_flow_avg=PredY_flow_avg(settings.Flow_list);
Gtest_pnasal_avg=Gtest_pnasal_avg(settings.Pnasal_list);
PredY_pnasal_avg=PredY_pnasal_avg(settings.Pnasal_list);


%% checking and cleaning
if 0
    if length(PredY_pnasal) ~= length(Gtest_pnasal); keyboard; end              % check same length
    if nnz(isnan(PredY_pnasal))>0 || nnz(isnan(GS))>0; keyboard; end            % check for any NaN's
    if nnz(PredY_pnasal<0)>0; Predy(PredY_pnasal<0)=0; end                      % force lower limit
    if nnz(PredY_pnasal>maxG)>0; Predy(PredY_pnasal>maxG)=1.5; end              % force upper limit
end

%% Figure 4 - flow
% no Figure E3 for Train Pnasal - Test Pnasal
SummaryStats_Fig4 = [NaN, NaN, NaN, NaN, NaN, NaN];

%% Figure 6 - pnasal
figure(6); clf(figure(6)); fig = gcf;
fig.Color = [1 1 1]; fig.Units = 'inches';
fig.Position = [1  1   16   5];

% flow predy vs pnasal predy
subplot(1,3,1)
facealpha = 0.05; % was 0.04
facecolor = [0 0 0]; %was [0 0 0]
scatter(predy_pnasal.*100,predy_flow(BB_).*100,2,'filled','markerfacealpha',facealpha, 'markerfacecolor', facecolor); hold on;
%[F6A_r, ~] =corr(predy_pnasal,predy_flow(BB_)); 
%w = ones(length(predy_pnasal),1);
%[F6A_Rsq,~,~,~,~]=glmfitFast(predy_pnasal,predy_flow(BB_),w,1);
[F6A_Rsq, F6A_r] = UnivariateStats(predy_pnasal,predy_flow(BB_),weights_pnasal); % don't use without weights
text(110, 15, ['R^2 = ', num2str(round(F6A_Rsq,2))]); 
text(110, 5, ['r = ', num2str(round(F6A_r,2))]);
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
text(110, 15, ['R^2 = ', num2str(round(F6B_Rsq,2))]); 
text(110, 5, ['r = ', num2str(round(F6B_r,2))]); 
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
text(80, 7, ['R^2 = ', num2str(round(F6C_Rsq,2))]); 
text(80, 0, ['r = ', num2str(round(F6C_r,2))]); 

text(-45, 108, 'C', 'FontSize', 20, 'FontWeight', 'Bold'); % Add label to plot space

str = ['..\Figures\CorrelationTesting\Figure_6',settings.experimentnumber];

switch savefigas
    case 'saveasTIFF'; print(fig, str, '-dtiff', '-r1000');
    case 'saveasPNG'; saveas(fig, str, 'png');
    case 'saveasFIG'; savefig(str);
end
if closefigs; close(fig); end

SummaryStats_Fig6 = [F6A_r, F6A_Rsq, F6B_r, F6B_Rsq, F6C_r, F6C_Rsq];

%%  Figure E1 - flow VE  Vs pnasal VE
pnasal_veveup = PtData_pnasal.VE(BB_) ./ PtData_pnasal.Veup(BB_);
flow_veveup = PtData_flow.VE(BB_) ./ PtData_flow.Veup(BB_);
pnasal_veveup_t = pnasal_veveup.^1.5;

squareplot = 1;
sigdig = '%.2f'; %2;
facealpha = 0.05;
facecolor = [0.1 0.1 0.1];

figure(28); clf(figure(28)); fig = gcf;
fig.Color = [1 1 1]; fig.Units = 'inches';
if squareplot
    fig.Position = [-11  1   8   8];
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

str = ['..\Figures\', 'Figure_E1_wBias'];
% switch savefigas     
%     case 'saveasTIFF'; print(fig, str, '-dtiff', '-r1000');     
%     case 'saveasPNG'; saveas(fig, str, 'png');      
%     case 'saveasFIG'; savefig(str);  
% end
if closefigs; close(fig); end

BiasResults = [SmBias_050, SmBias_067, SmBias_075, SmBias_100;
LgBias_050, LgBias_067, LgBias_075, LgBias_100];

BiasResults_Tbl = array2table(round(BiasResults,2), ...
    'VariableNames', {'half', 'twothirds' 'threequart' 'one'}, ...
    'RowName', {'smallBB', 'largeBB'} ...
  )

%% Histogram of ftrs selected during training 
ftr_array = NaN(54, ftrnum);
for subj=1:54 % set pt num or all, no data for PT=1
    if ~ismember(subj, settings.Pnasal_list) %PT_list)
        continue
    end
    ftr_array(subj,:) = labels_Step_Subj_pnasal{subj,ftrnum}; % only one list of labels, same for both flow and pnasal
end
ftr_array = ftr_array(settings.Pnasal_list,:);
ftr_array_linear = ftr_array(:);
ctrs = 1:1:size(Amatrix2_pnasal,2);  % ToDo: make sure max ctrs covers all features
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
score=zeros(1,size(Amatrix2_pnasal,2));
maxscore=0; %if large score is good
%maxscore = (size(labels_Step_Subj,2))^2; %if small score is good
itersize = size(labels_Step_Subj_pnasal,1);
try
    for i=1:itersize
        score1 = maxscore + 0*score;
        If = [];
        for j=1:ftrnum %size(labels_Step_Subj,2)
            temp = labels_Step_Subj_pnasal{i,j};
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
    
%% Alternatively, find best features through reverse stepwise regression
ftr_array_2 = NaN(54, ftrnum);
for subj=1:54 % set pt num or all, no data for PT=1
    if ~ismember(subj, settings.Pnasal_list) %PT_list)
        continue
    end
   ftr_array_2(subj,:) = labels_Step_Subj_pnasal{subj,ftrnum};
end
% Start with the list of unique ftrs selected during training at set stop point
ftr_array_2 = ftr_array_2(settings.Pnasal_list,:);
uniqueftrs = unique(ftr_array_2);

% do backwards elimination
RemoveList = [];
while 1   
    [~,Pvals,~,~,~]=glmfitFast(Amatrix2_pnasal(:,uniqueftrs),GS,weights,1); 
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
TopElimFtrs = fliplr(RemoveList');


%% run model using only top 25 features
NfeaturesOpt=ftrnum;
if 0 % either of the histogram methods
    if 1 % unweighted histogram
        I = I_uw;
    else % weighted histogram
        I = I_w';
    end
    LabelsOrdered = Labels(I(1:NfeaturesOpt));
    TopFtrs = I(1:NfeaturesOpt);
else % else using backward elimination of unique ftrs at threshold
    LabelsOrdered = Labels(TopElimFtrs(1:NfeaturesOpt));
    TopFtrs = TopElimFtrs(1:NfeaturesOpt);
end

%% Final pnasal model using all data, uses selected N optimal features (NfeaturesOpt)
% No final model for flow
Summary_FinalMdl_flow = [NaN,NaN];

%% Final pnasal model using all data, uses selected N optimal features (NfeaturesOpt)
Ftrs = TopFtrs;
Labels = [{'Intercept'};LabelsOrdered];
FtrVals_pnasal = Amatrix2_pnasal(:,Ftrs);
[Rsq,Pvals,RMSE,betas_pnasal,FinalMdlPredY_pnasal]=glmfitFast(FtrVals_pnasal(BB_,:),GS(BB_),weights(BB_),1); 
FinalMdlPredY_pnasal(FinalMdlPredY_pnasal>maxG)=maxG; % set upper limit
FinalMdlPredY_pnasal(FinalMdlPredY_pnasal<0) = 0; % set lower limit
[Rsq_2, r_2] = UnivariateStats(FinalMdlPredY_pnasal, GS(BB_), weights(BB_));
finalmdl_pnasal = fitglm(FtrVals_pnasal(BB_,:),GS(BB_),'weights',weights(BB_));
finalmdl_summary_table_pnasal = table([0;Ftrs'], Labels, betas_pnasal, finalmdl_pnasal.Coefficients.SE(1:end), Pvals, ...
        'VariableNames', {'FtrNum', 'FtrName', 'Betas', 'SE', 'Pval'}); 
str=['Final model r=', num2str(r_2),' ; Rsquared=', num2str(Rsq_2)]; disp(str)
Summary_FinalMdl_pnasal = [r_2,Rsq_2];

%% Output summary data
maxDataSize
settings.experimentnumber
SummaryStats_FF = [SummaryStats_Fig4, SummaryStats_FigE3, Summary_FinalMdl_flow];
SummaryStats_FP = [SummaryStats_Fig6, SummaryStats_FigE4, Summary_FinalMdl_pnasal];
SummaryStats = [SummaryStats_FF;SummaryStats_FP];


%% 
% could end here for the moment
disp('done');

%% Find univariate performance of 'best' flow features
if 0
rVsFD = NaN(length(Ftrs),1);
rVsPneumo = NaN(length(Ftrs),1);

r2VsFD = NaN(length(Ftrs),1);
r2VsPneumo = NaN(length(Ftrs),1);

bias =  NaN(length(Ftrs),1);

for ft=1:length(Ftrs)   
    % R vs flow:drive
    FtrVal = Amatrix_pnasal_matched(:,Ftrs(ft));
    [Rsq,~,~,~]=glmfitFast(FtrVal,Gtest_matched,weights,1);
    r2VsFD(ft) = Rsq;
    R_ = weightedcorrs([FtrVal,Gtest_matched],weights);
    rVsFD(ft) = R_(1,2);
    
    % R vs equivalent feature measured with Pneumotach
    PneumoFtrVal = Amatrix_flow_matched(:,Ftrs(ft));
    [Rsq,~,~,~]=glmfitFast(FtrVal,PneumoFtrVal,weights,1);
    r2VsPneumo(ft) = Rsq;
    
    R_ = weightedcorrs([FtrVal,PneumoFtrVal],weights);
    rVsPneumo(ft) = R_(1,2);
    
    % also calculate bias, as (median value Pneumotach) / (median value Pnasal)
    bias(ft) = median(PneumoFtrVal) / median(FtrVal);
    
    % figure
    if 0
        facealpha = 0.05; % was 0.08
        facecolor = [0.2 0.2 0.2]; %was [0.1 0.1 0.1]
        
        figure(101); clf(figure(101)); fig = gcf;
        fig.Color = [1 1 1]; fig.Units = 'inches';
       fig.Position = [-19   5  12   4.5];
       
       subplot(1,2,1); 
       set(gca,'box','off','tickdir','out','fontname','arial narrow');
       scatter(FtrVal,Gtest_matched,2,'filled','markerfacealpha', facealpha, 'markerfacecolor', facecolor); hold on;
       xlabel('Feature value (Pneumotach)'); ylabel('flow:drive');
       str=['R=',num2str(rVsFD(ft)), ',   R^2=',num2str(r2VsFD(ft))]
       title(str); axis square
       
       subplot(1,2,2);
       scatter(FtrVal,PneumoFtrVal,2,'filled','markerfacealpha', facealpha, 'markerfacecolor', facecolor); hold on;
       xlabel('Feature value (Pneumotach)'); ylabel('Feature value (Pnasal)');
       str=['R=',num2str(rVsPneumo(ft)), ',   R^2=',num2str(r2VsPneumo(ft))]
       title(str); axis square
       
       title_str = ['Ftr ', num2str(ft), ' -- ', LabelsOrdered{ft}];
       title_str = regexprep(title_str,'[_,:{\\}]',' ');
       suptitle(title_str);
       save_str = ['..\Figures\Univariate\', 'Pnasal_', title_str]; %
       saveas(fig, save_str, 'png');
    end
end

Univar_summary = table(rVsFD, rVsPneumo, bias, ...
        'VariableNames', {'RvsFlowDrive', 'RvsPneumo', 'Bias'});

Univar_summary_ = table([NaN;rVsFD], [NaN;rVsPneumo], [NaN;bias], ...
        'VariableNames', {'RvsFlowDrive', 'RvsPneumo', 'Bias'});

FinalModelDetails = [finalmdl_summary_table,Univar_summary_];

end


