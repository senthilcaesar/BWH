%% External Validation
% load oral appliance data and feature set
% find breaths likely and not likely flow limited
% (by obs-Ap, obs-hyp, and Ar)
% load final model
% run model on OA data to get estimated flow:drive (FD_est)
% Do bee swarm plot and confusion matrix
% get 'pseudo' R value as summary performance metric

close all
clear
clc

mypath = 'C:\Users\uqdmann\Dropbox\QAO\FlowAndEdiPesQAO_CodeForPublication';
addpath(mypath);
cd(mypath);

%% options
TransformTheData = 1;   % set as 1 to do tranforms, or 0 to use unadjusted data
addextratransform = 0; 	% set as 1 to do extra transforms
ReplaceApneaLowFlow = 1; % set as 1 to add apnea and low flow breaths back into predy and gtest
OA_data = 1;

%datadir = 'C:\PSG_Data\FlowDrive\FeatureSpaces\';
%datadir = '..\FeatureSpaces\';  %

%% open file
if OA_data
    datadir = 'C:\PSG_Data\OralAppliance\Analyzed\Publication\';
    filename = [datadir, 'FlowOA_25Hz_ReNormalized_FeatureSpace_NoRef_VEVeup_Clean.mat'];
else
    datadir = 'C:\PSG_Data\FlowDrive_MESA\Analyzed\';
    filename = [datadir, 'MESA_FlowDrive_20180718_FeatureSpace_NoRef_VEVeup_Clean_new.mat'];
end
str=['Loading ' filename]; disp(str);

try
    load(filename);
catch me
    disp(me.getReport);
end

%clearvars -except PtData FeatureNames Amatrix %% tidy up, if opening one of the old files

%% Just in case, look at NaN again
allnanrows = sum(isnan(Amatrix),2)==size(Amatrix,2);
if nnz(allnanrows)~=0
    str = ['Removing ', num2str(sum(allnanrows)), ' breaths that contain NaN''s']; disp(str);
    Amatrix(allnanrows,:)=[];
    PtData(allnanrows,:)=[];
    Fnan=sum(isnan(Amatrix)|isinf(Amatrix))/size(Amatrix,1);
    if nnz(Fnan)~=0
        disp('NaN''s or non-finite data remains - further investigation required');
        keyboard;
    end
else
    str = ['Zero NaN-breaths were found']; disp(str);
end

%% use only a given set of features, determined from previous run
if 0
    load('FtrsToUse','LabelsOrderedOpt');
    FtrsToInclude = LabelsOrderedOpt;
    Ind = [];
    for i=1:length(FtrsToInclude)
        temp=find(startsWith(FeatureNames.Name,FtrsToInclude(i)));
        if ~isempty(temp)
            Ind(end+1)=temp;
        end
    end
    FeatureNames = FeatureNames(Ind,:);
    Amatrix = Amatrix(:,Ind);
end

%% set the data matrix to use. either (0) unadjusted or (1) tranformed
if TransformTheData
    [Amatrix2, Labels] = DoDataMatTransform(Amatrix, FeatureNames, addextratransform);
else
    Amatrix2 = [Amatrix];
    Labels = FeatureNames.Name;
end

%% Prepare Oral Appliance data
Gtest_ = ones(size(PtData,1),1)*2; % set a vector of two's,
% that we will mark FL as zero, nonFL as one, and remove the twos

% NonFL breaths are:
%  Arousal, by at least 2 breaths (so using A1, which is threshold, more like about 4 breaths)
%  VE > Veup*0.9, VE is good
%  Etype == 0, not a clinically scored event
if OA_data % normal criteria
    Gtest_(PtData.A1==1 & PtData.VE>(PtData.Veup*0.9) & PtData.Etype==0)=1; 
else % criteria for MESA data
    Gtest_(PtData.Hypnog==4 & PtData.VE>(PtData.Veup*0.9) & PtData.Etype==0)=1; 
end

% FL breaths are:
%  NotAr
%  VE < Veup*0.7,
%  Etype 2 = ApO, 4 = HypO. Could extend to include 5 = M.
% use the clinically scored events (hypop/obstructive) as true labels
if OA_data
    Gtest_(PtData.NotAr==1 & PtData.VE<(PtData.Veup*0.7) & (PtData.Etype==2 | PtData.Etype==4))=0;
else
    Gtest_(PtData.VE<(PtData.Veup*0.5) & (PtData.Etype==2 | PtData.Etype==4))=0;
end

exclude = Gtest_==2;
str = [num2str(length(exclude)), ' breaths available for analysis']; disp(str)
str = [num2str(nnz(Gtest_==1)), ' breaths were labelled as NonFL']; disp(str)
str = [num2str(nnz(Gtest_==0)), ' breaths were labelled as FL']; disp(str)
str = [num2str(nnz(exclude)), ' breaths were unlabelled and removed']; disp(str)

% What is the breakdown of the Gtest_==2 breaths, i.e. neither FL or NonFL
VEVEup = PtData.VE./PtData.Veup;
str = [num2str(nnz(VEVEup>0.7 & VEVEup<0.9)),...
    ' of the removed breaths had 70% > VEVeup > 90%']; disp(str)
Ar1 = PtData.A1==0; % find arousal state
str = [num2str(nnz(VEVEup(Gtest_==2)>0.9&Ar1(Gtest_==2))),...
    ' of the removed breaths had VEVeup >90%, but were not far enough away from sleep']; disp(str)

% Ar    - Arousal breath
% NotAr - Ar=1 is always NotAr=0, Ar=0 is NotAr=1 but only after delay of two breaths
% A2    - Arousal breath, at least 2 breaths away from sleep
% A1    - Arousal breath, at least 'threshold' breaths away from sleep (this average to about 4 breaths)

% remove breaths that are neither FL, nor NonFL
Gtest_(exclude)=[];
PtData(exclude,:)=[];
Amatrix(exclude,:)=[];
Amatrix2(exclude,:)=[];

% optionally add apnea and low flow breaths in as 0.1 breaths
if ReplaceApneaLowFlow
    numRemBB_Ap = sum(RemovedBB_Apnoea.NumBB);
    numRemBB_LF = sum(RemovedBB_LowFlow.NumBB);
    numRemBB = numRemBB_Ap + numRemBB_LF;
else
    numRemBB = 0;
end

%% add apnea and low flow breaths back in
Gtest_ = [Gtest_;zeros(numRemBB,1)];

%% Weights
clear Ndata
if 0
    maxG=1.5;
    dx=0.2;
    xbins=[0 0.3:dx:0.9 maxG];
    %xbins=[0 0.1:0.05:1.1 1.5];
else
    xbins=[-Inf 0.5 +Inf];  % should this be 0.7, NO - data is zeros and ones, so makes no difference
end

for i=1:length(xbins)-1
    Ix=Gtest_>xbins(i)&Gtest_<=xbins(i+1);
    Ndata(i)=sum(Ix);
end
%Ndata = Ndata.^0.5;
weightsbins = 1./(Ndata);
%weightsbins = [2 1 1 1 0.5];
weightsbins = weightsbins/mean(weightsbins);
%weightsbins = weightsbins/weightsbins(end);

weights = NaN*Gtest_;
for i=1:length(xbins)-1
    Ix=Gtest_>=xbins(i)&Gtest_<=xbins(i+1);
    weights(Ix)=weightsbins(i);
end
weights = weights/nanmean(weights); % nnz(~isfinite(weights))

useweights=1;
if ~useweights %overwrite, use no weights
    weights_ = ones(length(weights),1);
end

%% Load final model - just placeholder data for time being
if 0
    % five feature model
    ModelFtrs = [231 230 201 200 198]; % Labels(231)
    ModelBetas = [1.3 -2.32 -4.17 -0.24 -0.68 -0.3];
else
    % 25 feature model
    Final25FtrModel = load('FinalModel_FtrsAndBetas.mat'); %Final25FtrModel = load('FinalModel_Table.mat'); 
    ModelFtrs = Final25FtrModel.FinalFtrs_sorted';
    ModelBetas = Final25FtrModel.FinalBetas_sorted';
end
maxG = 1.5;

%% run final model on OA features
predy = [ones(size(Amatrix2,1),1) Amatrix2(:,ModelFtrs)]*ModelBetas';
predy(predy>maxG)=maxG; % set upper limit
predy(predy<0) = 0; % set lower limit

%% add apnea and low flow breaths back in
predy = [predy;ones(numRemBB,1)*0.09];

%% do bee swarm plot
TickFntSz = 12;
LabelFntSz = 16; 
FntSz = 16; 
savefigas = '';
closefigs = 0;
Yval = Gtest_;

figure(86); clf(figure(86)); hold on;
fig = gcf; fig.Color = [1 1 1]; fig.Units = 'inches';
fig.Position = [29.5  2.5   9   4.5];
subplot(1,2,1); % probability plot (bee swarm plot)
if 1
    scatter(100*predy,100*Yval+10*randn(length(Yval),1),4,'filled','markerfacealpha',0.2, 'markerfacecolor', [0.1 0.1 0.1]);
    hold on;
    dx=0.2; xbins=[0:0.1:1];
    clear meanX meanY upperY lowerY
    for i=1:length(xbins)-1
        Ix=predy>xbins(i)&predy<xbins(i+1);
        meanX(i)=mean(predy(Ix));
        meanY(i)=nansum(Yval(Ix).*weights(Ix))/nansum(weights(Ix));
        upperY(i)=1*nanstd(Yval(Ix))+meanY(i);
        lowerY(i)=meanY(i)-1*nanstd(Yval(Ix));
    end
    ax = gca; set(ax,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz);axis square 
    ax.XTick=[0:25:150]; xticklabels(ax, {'0', '25', '50', '75', '100', '125'});
    yticks([0,100]);
    yticklabels({'Obstructed','Normal'}); ytickangle(ax,90);
    xlim([0 150]); axis square 
    %xlabel('Probability of being Normal or Obstructed using Flow Shape (%)');
    xlabel('Estimated {\itflow:drive} (%)', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');  
    ylabel('Gold Standard Classification', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
else % alternate probability plot, plotSpread
    Ndatapoints=Inf;
    if Ndatapoints>length(predy), Ndatapoints=length(predy); end
    data=100*predy(1:Ndatapoints);
    catIdx=Yval(1:Ndatapoints);
    data={data(catIdx~=1),data(catIdx==1)};
    %figure(5); clf(5);
    largeNfactor=2*(Ndatapoints/300)^0.5; %magic equation
    if largeNfactor<1,largeNfactor=1; end
    [~,data1,temp]=plotSpread(data,'xNames',{'F.limited','Normal'},'categoryIdx',catIdx,...
        'categoryMarkers',{'.','.'},'categoryColors',{'k','k'},'binWidth',0.1/largeNfactor,...
        'xyOri','flipped','magicNumber',1-(0.1/largeNfactor));
    ax = gca; set(ax,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz);axis square 
    xlabel('Pr Normal using flow shape LogReg, %')
    scatter(data1(:,2),100*(data1(:,1)-1),1,'filled','markerfacealpha',0.1)
    set(gca,'tickdir','out');
end

currentYlim=ylim();
plot([70, 70], [currentYlim(1), currentYlim(2)],'r-', 'linewidth', 1.5);

subplot(1,2,2); %confusion matrix
g = NaN*Yval;
ghat = NaN*predy;
classcutoffs = [0.70]; %normal mild moderate severe v.severe / normal borderline mild moderate severe
Nclasses=length(classcutoffs)+1;
% DMconfusionMat, 1 = NFL, 0 = FL
g_ = Yval; % true/label/class % nnz(g_)
ghat_ = NaN*predy; % prediciton
ghat_(predy>classcutoffs(1))=1;
ghat_(predy<=classcutoffs(1))=0;
trueNFL = nnz(g_==1 & ghat_ ==1);
falseNFL = nnz(g_==0 & ghat_ ==1);
trueFL = nnz(g_==0 & ghat_ ==0);
falseFL = nnz(g_==1 & ghat_ ==0);
dlmXtab = [trueNFL; falseNFL; trueFL; falseFL]
t = (dlmXtab./sum(dlmXtab))*100;
C_Total = [t(1), t(4); t(2), t(3)];
g(Yval>classcutoffs(1))=1;
ghat(predy>classcutoffs(1))=1;
for i=2:Nclasses
    g(Yval<classcutoffs(i-1))=i;
    ghat(predy<classcutoffs(i-1))=i;
end
[C,order] = confusionmat(g,ghat);
%sumactual=sum(C')';
sumactual=sum(C,2);
sumestimated=sum(C);
%C_Factual=C./sum(C')'
C_Factual=C./sum(C,2)*100 %rows are actual, cols are estimated
C_Festimated=C./sum(C)*100
C_Total2 = (C./sum(C(:)))*100; % line from lin reg
AccA_C_Factual = C_Factual.*diag(ones(1,length(C)));
AccB_C_Factual = C_Factual.*(diag(ones(1,length(C)-1),1)+diag(ones(1,length(C)-1),-1));
AccAactual = mean(sum(AccA_C_Factual,2));
AccBactual = mean(sum(AccA_C_Factual,2)) + mean(sum(AccB_C_Factual,2));
AccA_C_Festimated = C_Festimated.*diag(ones(1,length(C)));
AccB_C_Festimated = C_Festimated.*(diag(ones(1,length(C)-1),1)+diag(ones(1,length(C)-1),-1));
AccAestimated = mean(sum(AccA_C_Festimated));
AccBestimated = mean(sum(AccB_C_Festimated)) + mean(sum(AccB_C_Festimated));
ACCs = [AccAactual AccBactual; AccAestimated AccBestimated]
%second col gives accuracy if accepting next-category error as correct

x = order'; y = order';
if 1; C1 = C_Festimated; else; C1 = C_Factual; end
%C1 = C_Total; % using the total proportion of breaths in each division
%C = [0 1 2 3 ; 1 2 3 4 ; 2 3 4 5];
xSplit = diff(x)/2;                 % Find edge points
ySplit = diff(y)/2;
xEdges = [x(1)-xSplit(1) x(2:end)-xSplit x(end)+xSplit(end)];
yEdges = [y(1)-ySplit(1) y(2:end)-ySplit y(end)+ySplit(end)];
[XGrid, YGrid] = meshgrid(xEdges,yEdges);
YGrid = flipud(YGrid);              % To match expected behavior
XGrid = fliplr(XGrid);
C2 = [[C1 zeros(size(C1,1),1)] ; zeros(1,size(C1,2)+1)];% Last row/col ignored
customcmap = GetCustomColorMap('gray');
set(gcf,'colormap',customcmap);
pcolor(XGrid,YGrid,C2)
hold on                             % Plot original data points
[X,Y] = meshgrid(x,y);
C1 = flipud(C1);
C1 = fliplr(C1);
C1 = C1';
for i=1:size(C1,1)
    for j=1:size(C1,2)
        if C1(i,j)<18%((max(max(C1))-min(min(C1)))/2+min(min(C1)))
            textcolor=[1 1 1];
        else
            textcolor=[0 0 0];
        end
        text(x(i),y(j),num2str(round(C1(i,j),1)),'color',textcolor,'horizontalalignment','center','fontname','arial narrow')
    end
end
labeltext={'Normal','Obstructed'};
ax = gca; set(ax,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz); axis square
yticks(y); yticklabels(gca,fliplr(labeltext));
ytickangle(ax,90);
xticks(x); xticklabels(gca,fliplr(labeltext));
ylabel('Gold Standard Classification', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
xlabel('Estimated Classification', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow'); 

str = ['..\Figures\', 'Figure_E6_ExtVal_OA_cutat070'];
%str = ['..\Figures\', 'Figure_E6_ExtVal_OA_cutat070_5Ftrs'];
switch savefigas
    case 'saveasTIFF'; print(fig, str, '-dtiff', '-r1000');
    case 'saveasPNG'; saveas(fig, str, 'png'); 
    case 'saveasFIG'; savefig(str); 
end
if closefigs; close(fig); end

