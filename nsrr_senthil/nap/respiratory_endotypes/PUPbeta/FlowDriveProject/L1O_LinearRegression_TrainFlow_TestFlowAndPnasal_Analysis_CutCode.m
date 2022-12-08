

%% Find Top X unique features based on scores
% if using transformed data, need to strip off the suffix
if settings.TransformTheData; lengthsuf=4; else; lengthsuf=0; end
MaxUniquefeatures=50;
LabelsOrdered_NoSuffix=LabelsOrdered;
for i=1:length(LabelsOrdered_NoSuffix)
    LabelsOrdered_NoSuffix{i}=LabelsOrdered_NoSuffix{i}(1:length(LabelsOrdered_NoSuffix{i})-lengthsuf);
end
LabelsOrdered_NoSuffixUnique=unique(LabelsOrdered_NoSuffix,'stable');
LabelsOrdered_NoSuffixUnique(MaxUniquefeatures+1:end)=[];

%%






%% error

% error per breath, nasal pred to flow pred
100*mean(abs(predy_flow-predy_pnasal))
[rho, pval] = corr(predy_flow, predy_pnasal);

% error per breath nasal pred to gold std
100*mean(abs(predy_flow-Yval))
[rho, pval] = corr(predy_flow, Yval); 

% error per pt median, nasal pred to gold std
100*mean(abs(PredY_avg-Gtest_avg))


%% get the AHI data
[AHI_perPT, AHI_perPT_table] = getAHI_postanalysis();
AHI_perPT_ = AHI_perPT(settings.Pnasal_list,1);
AHI_perPT_FlowPts = AHI_perPT(~isnan(AHI_perPT(:,1)),1);


%% histogram of pred vs actual VEVdrive, and do per pt processing in loop
% get per pt median data
% get per pt <threshold data

if 1; edges = [0:0.1:1.1]; else edges = settings.xbins; end

figure(22); clf(figure(22)); fig = gcf;
fig.Color = [1 1 1]; fig.Units = 'inches';
fig.Position = [0.5   0.5   19    5];

subplot(1,3,1);
histogram(GS,edges, 'facealpha',0.2, 'edgealpha', 0.2, 'normalization','probability'); hold on;
currentYlim = ylim();
plot([median(GS), median(GS)], [0, currentYlim(2)],'k-', 'linewidth', 2);
titlestr = ['Actual (median = ', num2str(round(median(GS),2)), ')']; 
title(titlestr); xlabel('{\itflow:drive} (%)'); ylabel('Relative probability'); %yticks([]);
xlim([-0.05 1.15]); ax = gca; ax.XTick=[0:0.2:1]; xticklabels(ax, {'0', '20', '40', '60', '80', '100'});
box off

subplot(1,3,2);
histogram(PredY_pnasal,edges, 'facealpha',0.2, 'edgealpha', 0.2, 'normalization','probability'); hold on;
currentYlim = ylim();
plot([median(PredY_pnasal), median(PredY_pnasal)], [0, currentYlim(2)],'k-', 'linewidth', 2);
titlestr = ['Pnasal Predicted (median = ', num2str(round(median(PredY_pnasal),2)), ')']; 
title(titlestr);xlabel('{\itflow:drive} (%)'); ylabel('Relative probability'); 
xlim([-0.05 1.15]); ax = gca; ax.XTick=[0:0.2:1]; xticklabels(ax, {'0', '20', '40', '60', '80', '100'});
box off

subplot(1,3,3);
histogram(PredY_flow,edges, 'facealpha',0.2, 'edgealpha', 0.2, 'normalization','probability'); hold on;
currentYlim = ylim();
plot([median(PredY_flow), median(PredY_flow)], [0, currentYlim(2)],'k-', 'linewidth', 2);
titlestr = ['Flow Predicted (median = ', num2str(round(median(PredY_flow),2)), ')']; 
title(titlestr);xlabel('{\itflow:drive} (%)'); ylabel('Relative probability'); 
xlim([-0.05 1.15]); ax = gca; ax.XTick=[0:0.2:1]; xticklabels(ax, {'0', '20', '40', '60', '80', '100'});
box off

supstr=['Histograms of Actual and Pnasal Predicted {\itflow:drive} (All pts, incl Ap)'];
suptitle(supstr);
str = ['..\Figures\', supstr];
% switch savefigas     
%     case 'saveasTIFF'; print(fig, str, '-dtiff', '-r1000');     
%     case 'saveasPNG'; saveas(fig, str, 'png');      
%     case 'saveasFIG'; savefig(str);  
% end
if closefigs; close(fig); end

%% histograms for individuals
Gtest_avg = NaN(54,1); PredY_avg = NaN(54,1);
thres1 = 0.5; thres2 = 0.7;
Gtest_thres1 = NaN(54,1); Gtest_thres2 = NaN(54,1);
PredY_thres1 = NaN(54,1); PredY_thres2 = NaN(54,1);

doHistogramsPlots = 0;

for subj=1:54 % set pt num or all, no data for PT=1
    if ~ismember(subj, settings.Pnasal_list) %PT_list)
        continue
    end
    
    Isubj = PTarray==subj;      % get the pt BB's
    Gtest_pt_flow = GS(Isubj);    % actual
    PredY_pt_pnasal = PredY_pnasal(Isubj);    % pred
    
    % medians
    Gtest_avg(subj) = median(Gtest_pt_flow);
    PredY_avg(subj) = median(PredY_pt_pnasal);
    
    
    % below threshold
    Gtest_thres1(subj) = nnz(Gtest_pt_flow<thres1);
    Gtest_thres2(subj) = nnz(Gtest_pt_flow<thres2);
    PredY_thres1(subj) = nnz(PredY_pt_pnasal<thres1);
    PredY_thres2(subj) = nnz(PredY_pt_pnasal<thres2);
    
    if doHistogramsPlots
        figure(300+subj); clf(figure(300+subj));fig = gcf;
        fig.Color = [1 1 1]; % set background colour to white
        fig.Units = 'inches';
        fig.Position = [0.5   0.5   12    5];
        
        subplot(1,2,1);
        h1 = histogram(Gtest_pt_flow,edges, 'facealpha',0.2, 'edgealpha', 0.2, 'normalization','pdf');%'probability');
        hold on;
        currentYlim = ylim();
        %plot([median(Gtest_pt), median(Gtest_pt)], [0, max(h1.Values)],'k-', 'linewidth', 2);
        plot([median(Gtest_pt_flow), median(Gtest_pt_flow)], [0, currentYlim(2)],'k-', 'linewidth', 2);
        titlestr = ['Actual (median = ', num2str(round(median(Gtest_pt_flow),2)), ')'];
        title(titlestr); xlabel('{\itflow:drive} (%)'); ylabel('Relative probability'); %yticks([]);
        xlim([-0.05 1.15]); ax = gca; ax.XTick=[0:0.2:1]; xticklabels(ax, {'0', '20', '40', '60', '80', '100'});
        box off

        subplot(1,2,2);
        h2 = histogram(PredY_pt_pnasal,edges, 'facealpha',0.2, 'edgealpha', 0.2, 'normalization','pdf');%'probability');
        hold on;
        currentYlim = ylim();
        %plot([median(PredY_pt), median(PredY_pt)], [0, max(h2.Values)],'k-', 'linewidth', 2);
        plot([median(PredY_pt_pnasal), median(PredY_pt_pnasal)], [0, currentYlim(2)],'k-', 'linewidth', 2);
        titlestr = ['Predicted (median = ', num2str(round(median(PredY_pt_pnasal),2)), ')'];
        title(titlestr); xlabel('{\itflow:drive} (%)'); ylabel('Relative probability');
        xlim([-0.05 1.15]); ax = gca; ax.XTick=[0:0.2:1]; xticklabels(ax, {'0', '20', '40', '60', '80', '100'});
        box off
        
        supstr = ['Histograms of Actual and Pnasal Predicted {\itflow:drive}, pt ',num2str(subj),', AHI ', num2str(round(AHI_perPT(subj,1))), ', withApBB'];
        suptitle(supstr);
        
        % save
        str = ['..\Figures\', supstr];
        switch savefigas
            case 'saveasTIFF'; print(fig, str, '-dtiff', '-r1000');
            case 'saveasPNG'; saveas(fig, str, 'png');
            case 'saveasFIG'; savefig(str);
        end
        if closefigs; close(fig); end
    end
end

% tidy up
exclude = isnan(Gtest_avg); %a=find(exclude==0)
Gtest_avg(exclude) = []; PredY_avg(exclude) = [];
Gtest_thres1(exclude) = []; Gtest_thres2(exclude) = [];
PredY_thres1(exclude) = []; PredY_thres2(exclude) = [];
Gtest_thres1_pFL = Gtest_thres1./numBBinTest(:,5); % Proportion FL
Gtest_thres2_pFL = Gtest_thres2./numBBinTest(:,5);
PredY_thres1_pFL = PredY_thres1./numBBinTest(:,5);
PredY_thres2_pFL = PredY_thres2./numBBinTest(:,5);




%% Novel metrics, Median VE:Vdrive during sleep and Time with severe obstruction during sleep
% aka, Proportion of breaths FL
figure(27); clf(figure(27)); fig = gcf;
fig.Color = [1 1 1]; fig.Units = 'inches';
fig.Position = [-14   1   12    4.5];
%fig.Position = [0.5   0.5   12    4.5];

subplot(1,2,1);
[r_1, p_1] = plotregressionwithSEM(PredY_avg.*100, Gtest_avg.*100);
xlim([-5 110]); ylim([-5 110]);
ax = gca; ax.FontSize = 12; FntSz = 18;
ax.XTick=[0:25:125];
ax.YTick=[0:25:125];
xlabel('Pnasal Shape Predicted {\itflow:drive} (%)', 'FontSize', FntSz, 'FontName', 'Arial Narrow');
ylabel('Gold Standard {\itflow:drive} (%)', 'FontSize', FntSz, 'FontName', 'Arial Narrow'); 
axis square
str=['r = ', num2str(round(r_1,2))]; text(80, 7, str);
% titlestr = ['Patient Median {\itflow:drive}']; title(titlestr);

subplot(1, 2, 2);
[r_3, p_3] = plotregressionwithSEM(PredY_thres2_pFL.*100, Gtest_thres2_pFL.*100);
xlim([-5 110]); ylim([-5 110]); 
ax = gca; ax.FontSize = 12; FntSz = 18;
ax.XTick=[0:25:125]; 
ax.YTick=[0:25:125];
xlabel('Pnasal Shape Predicted % Flow Limited', 'FontSize', FntSz, 'FontName', 'Arial Narrow');
ylabel('Gold Standard % Flow Limited', 'FontSize', FntSz, 'FontName', 'Arial Narrow'); 
axis square
str=['r = ', num2str(round(r_3,2))]; text(80, 7, str);
% titlestr=['Proportion with Mod\Sev FL']; title(titlestr);

% subplot(1, 3, 3);
% [r_2, p_2] = plotregressionwithSEM(PredY_thres1_pFL.*100, Gtest_thres1_pFL.*100);
% xlim([-5 125]); xlabel('Pnasal Predicted %FL');
% ylim([-5 125]); ylabel('Gold Standard %FL');
% ax = gca; ax.XTick=[0:20:120];
% axis square
% str=['r = ', num2str(round(r_2,2))]; text(80, 7, str);
% titlestr=['Proportion with Severe FL']; title(titlestr);

% Add labels A B to plot space
subplot(1,2,1); hold on;
text(-45, 108, 'A', 'FontSize', 20, 'FontWeight', 'Bold');
subplot(1,2,2); hold on;
text(-45, 108, 'B', 'FontSize', 20, 'FontWeight', 'Bold');

str = ['..\Figures\', 'Figure_NovelMetrics'];
switch savefigas     
    case 'saveasTIFF'; print(fig, str, '-dtiff', '-r1000');     
    case 'saveasPNG'; saveas(fig, str, 'png');      
    case 'saveasFIG'; savefig(str);  
end
if closefigs; close(fig); end
end

%% After adjusting for AHI

% Median
[b,dev,stats]=glmfit([AHI_perPT_(:,1) PredY_avg],Gtest_avg);
p_1_adj = stats.p(3);
r_1_adj = stats.coeffcorr(1,3);

% mod\sev
[b,dev,stats]=glmfit([AHI_perPT_(:,1) PredY_thres2],Gtest_thres2);
p_2_adj = stats.p(3);
r_2_adj = stats.coeffcorr(1,3);

% severe
[b,dev,stats]=glmfit([AHI_perPT_(:,1) PredY_thres1],Gtest_thres1);
p_3_adj = stats.p(3);
r_3_adj = stats.coeffcorr(1,3);







%%







if 0

%%
if 0
    %% Pred VEVdrive vs Actual VEVdrive, for Flow(all), Flow(matched) and Pnasal
    % plot is technically "Actual (x-axis) Vs Pred (y-axis)"
    figure(20); clf(figure(20));
    fig = gcf;
    fig.Color = [1 1 1];
    fig.Units = 'inches';
    fig.Position = [  -12.2    8.5   12    4.5];
    
    subplot(1,3,1);
    scatter(predyL1O_array_pnasal_matched(:,ftrnum), Gtest_matched, 2, 'filled','markerfacealpha',0.2);
    hold on; lsline; title('Pnasal'); ylabel('Actual VEVdrive'); xlabel('Pred VEVdrive');
    xlim([-0.05 1.55]); ylim([-0.05 1.55]); axis square;
    subplot(1,3,2);
    scatter(predyL1O_array_flow_matched(:,ftrnum), Gtest_matched, 2, 'filled','markerfacealpha',0.2);
    hold on; lsline; title('Flow (matched)'); xlabel('Pred VEVdrive');
    xlim([-0.05 1.55]); ylim([-0.05 1.55]); axis square;
    subplot(1,3,3);
    scatter(predyL1O_array(:,ftrnum), Gtest_All, 2, 'filled','markerfacealpha',0.2);
    hold on; lsline; title('Flow (All)'); xlabel('Pred VEVdrive');
    xlim([-0.05 1.55]); ylim([-0.05 1.55]); axis square;
    suptitle(['Actual Vs Pred VEVdrive at ', num2str(ftrnum), ' ftrs']);
    %savestr = [datastr,' ', supstr];
end

%% GS VEVdrive vs Pnasal VEVdrive (not pred VE:Vdrive estimate, just using pnasal to get VE)
figure(16); clf(figure(16)); fig = gcf;
%[r_t, p_t] = plotregressionwithSEM(PtData_matched.g_Edi_Adj, PtData_pnasal_matched.g_Edi_Adj);
scatter(PtData_flow_matched.g_Edi_Adj, PtData_pnasal_matched.g_Edi_Adj,2,'filled','markerfacealpha',0.2); hold on;
[r,p] = corr(PtData_flow_matched.g_Edi_Adj, PtData_pnasal_matched.g_Edi_Adj);
scatter(PT_avg_gtest, PT_avg_predy,20,'filled');
xlabel('Flow VEVdrive (Gold standard)');ylabel('Pnasal VEVdrive');
xlim([-0.05 1.55]);ylim([-0.05 1.55]);
titlestr=['Flow VEVdrive (Gold standard) Vs Pnasal VEVdrive']; title(titlestr);
str = ['..\Figures\', titlestr];
%saveas(fig, str, 'png'); %savefig(str);


%% Time below threshold - WITH/WITHOUT adding Apnoea and LowFlow breaths back in
%clinscorerange = unique(PtData.Etype) % clinical scoring
% clinical scoring codes
% 2 obstructive ap
% 3 central ap
% 4 hypopnoea ob
% 5 mixed

% This has the updated removed counts
% updatecounts = 0;
% if updatecounts
%     myfilename = ['C:\PSG_Data\FlowDrive\FeatureSpaces\FlowDrive_only25Hz_FeatureSpace_AutoRef2_Edi_Clean.mat'];
%     load(myfilename, 'RemovedBB_Apnoea', 'RemovedBB_LowFlow');
% end

ptstats = [];
sleeponly = 1;
addbackAP = 1;
thres = 0.5;
thres2 = 0.7;
numBBinTest = NaN(54,4);

% testing
if 0
    rn = randperm(10000);
    rn = rn/10000;
    idx = randperm(10000);
    actualBelowThres = rn(idx<5000)<thres;
    numActualBBbelowthreshold = nnz(actualBelowThres)
    
    figure(1); clf(figure(1));
    stairs([1:1:height(PtData_matched)], PtData_matched.Hypnog, 'b'); hold on;
    stairs([1:1:height(PtData_matched)], PtData_pnasal_matched.Hypnog, 'r');
    t = find(abs(PtData_matched.Hypnog-PtData_pnasal_matched.Hypnog)>0);
    plot(t, ones(length(t),1),'bx');
    nnz(PtData_matched.Hypnog==4)
end

for subj=1:54
    if ~ismember(subj, settings.Pnasal_list)
        continue
    end
    if sleeponly % make Isubj
        Isubj=(PtData_flow_matched.PT==subj)&(PtData_flow_matched.Hypnog<4)&(PtData_flow_matched.Ar==0);
        numBBinTest(subj,1) = nnz(Isubj);
        %Isubj=(PtData_matched.PT==subj)&(PtData_matched.Hypnog<4)&(PtData_matched.NotAr==1);
        %numBBinTest(subj,2) = nnz(Isubj);
    else
        Isubj=(PtData_flow_matched.PT==subj);
        numBBinTest(subj,3) = nnz(Isubj);
        % Isubj=(PtData_matched.PT==subj)&(PtData_matched.Hypnog<4);
        % numBBinTest(subj,4) = nnz(Isubj);
    end
    
    clinscore = PtData_flow_matched.Etype(Isubj);
    numclinscored = nnz(clinscore==2|clinscore==4); % obstructive and hypopnoea
    predBelowThres = predyL1O_array(Isubj,ftrnum)<thres;
    numPredBBbelowthreshold = nnz(predBelowThres);
    actualBelowThres = Gtest_matched(Isubj)<thres;
    numActualBBbelowthreshold = nnz(actualBelowThres);
    predBelowThres2 = predyL1O_array(Isubj,ftrnum)<thres2;
    numPredBBbelowthreshold2 = nnz(predBelowThres2);
    actualBelowThres2 = Gtest_All(Isubj)<thres2;
    numActualBBbelowthreshold2 = nnz(actualBelowThres2);
    numBBtotal = size(predyL1O_array(Isubj),1);
    if addbackAP
        withAPstr = [' (incl Ap-O breaths)'];
        indAP = find(RemovedBB_Apnoea.Pt == subj);
        indLF = find(RemovedBB_LowFlow.Pt == subj);
        if ~isempty(indAP)
            numclinscored = numclinscored+sum(RemovedBB_Apnoea{indAP,[4 6 7]});
            numPredBBbelowthreshold = numPredBBbelowthreshold+RemovedBB_Apnoea{indAP,2};
            numActualBBbelowthreshold = numActualBBbelowthreshold+RemovedBB_Apnoea{indAP,2};
            numPredBBbelowthreshold2 = numPredBBbelowthreshold2+RemovedBB_Apnoea{indAP,2};
            numActualBBbelowthreshold2 = numActualBBbelowthreshold2+RemovedBB_Apnoea{indAP,2};
            numBBtotal = numBBtotal+RemovedBB_Apnoea{indAP,2};
        end
        if ~isempty(indLF)
            numclinscored = numclinscored+sum(RemovedBB_LowFlow{indLF,[4 6 7]});
            numPredBBbelowthreshold = numPredBBbelowthreshold+RemovedBB_LowFlow{indLF,2};
            numActualBBbelowthreshold = numActualBBbelowthreshold+RemovedBB_LowFlow{indLF,2};
            numPredBBbelowthreshold2 = numPredBBbelowthreshold2+RemovedBB_LowFlow{indLF,2};
            numActualBBbelowthreshold2 = numActualBBbelowthreshold2+RemovedBB_LowFlow{indLF,2};
            numBBtotal = numBBtotal+RemovedBB_LowFlow{indLF,2};
        end
    else
        withAPstr = [''];
    end
    Clinpercent = 100*(numclinscored/numBBtotal);
    Predpercent = 100*(numPredBBbelowthreshold/numBBtotal);
    Actualpercent = 100*(numActualBBbelowthreshold/numBBtotal);
    Predpercent2 = 100*(numPredBBbelowthreshold2/numBBtotal);
    Actualpercent2 = 100*(numActualBBbelowthreshold2/numBBtotal);
    ptstats = [ptstats; [subj,numBBtotal,numclinscored,...
        numPredBBbelowthreshold,numActualBBbelowthreshold,...
        numPredBBbelowthreshold2,numActualBBbelowthreshold2,...
        Clinpercent,...
        Predpercent,Actualpercent,...
        Predpercent2,Actualpercent2]];
end
ptsummary_withAp = array2table(ptstats, 'VariableNames', {'PT','TotalBB','Clin_FL',...
    'Pred_FL','Actual_FL',...
    'Pred_FL2','Actual_FL2',...
    'Clin_percent',...
    'Pred_percent','Actual_percent',...
    'Pred_percent2','Actual_percent2'});
sum(ptsummary_withAp.TotalBB)
numBBinTest = numBBinTest(settings.Pnasal_list,:);

%% average Gtest vs average PredY for each pt (medians...)
figure(24); clf(figure(24)); fig = gcf;
fig.Color = [1 1 1];
fig.Units = 'inches';
fig.Position = [ -12.2 -2.5  12.0000    4.5000];
settings.Pnasal_list = unique(PtData_flow_matched.PT); PT_listWpnasal = settings.Pnasal_list;
subplot(1, 2, 1);
Gtest_avg = NaN(54,1);
PredY_avg = NaN(54,1);
for i=ftrnum % 1:size(predyL1O_array,2) say at 20 ftrs
    for subj=1:54 %length(PT_list)
        if ~ismember(subj, settings.Pnasal_list)
            continue
        end
        %Isubj = PtData.PT==subj; % have not used weights(Isubj) ??
        if 0
            Isubj=(PtData_matched.PT==subj)&(PtData_matched.Hypnog<4)&(PtData_matched.Ar==0); % sleep only
        else
            Isubj=(PtData_flow_matched.PT==subj)&(PtData_flow_matched.NotAr==1);
        end
        Gtest_avg(subj) = median(Gtest_matched(Isubj));
        PredY_avg(subj) = median(predyL1O_array(Isubj,i));
        
        % make the data that includes the excluded Ap breaths
        indx = find(ptsummary_withAp.PT==subj);
        numPredFL = ptsummary_withAp.Pred_FL(indx);
        numActualFL = ptsummary_withAp.Actual_FL(indx);
        
        Gtest_avg_wAp(subj) = median([Gtest_matched(Isubj);ones(numActualFL,1)*0.1]);
        PredY_avg_wAp(subj) = median([predyL1O_array(Isubj,i);ones(numPredFL,1)*0.1]);
        
        scatter(predyL1O_array((Isubj),i),Gtest_matched(Isubj),2,'filled','markerfacealpha',0.5); hold on;
    end
end

% do once only.
Gtest_avg_wAp = Gtest_avg_wAp(settings.Pnasal_list);
PredY_avg_wAp = PredY_avg_wAp(settings.Pnasal_list);

scatter(PredY_avg, Gtest_avg, 50, 'filled','markerfacealpha',1); hold on;
% lsline;
xlim([0 1.5]); xlabel('PredY VE:VDrive');
ylim([0 1.5]); ylabel('Actual VE:VDrive');
title('(Sleep only)Indiv breaths and patient medians');
axis square

subplot(1, 2, 2);
[r_t, p_t] = plotregressionwithSEM(PredY_avg, Gtest_avg);
%scatter(PredY_avg, Gtest_avg, 50, 'filled'); hold on;
%hndl_ls = lsline;  hndl_ls.LineStyle = '-'; hndl_ls.Color = [0 0 1];
%hndl_ref = refline(1,0); hndl_ref.LineStyle = '--'; hndl_ref.Color = [0.5 0.5 0.5];
xlim([0 1.5]); xlabel('PredY VE:VDrive');
ylim([0 1.5]); ylabel('Actual VE:Vdrive');
axis square
[r, ~] = corr(PredY_avg(~isnan(PredY_avg)), Gtest_avg(~isnan(PredY_avg)));
titlestr = ['(Sleep only) Patient medians, r=', num2str(r)]; title(titlestr);

suptitle(['Actual VEVdrive Vs Pred VEVdrive']);
str = ['..\Figures\', 'MedianGtestVsMedianPredY_AutoRef2_SleepOnly_PnasalTnT'];
%saveas(fig, str, 'png'); %savefig(str);


%% median VEVdrive vs AHI
AHI_perPT_ = AHI_perPT(~isnan(AHI_perPT(:,1)),1);

excludelist = find(isnan(PredY_avg(:,1)));
AHI_perPT_ = AHI_perPT(:,1);
AHI_perPT_(excludelist) = [];
PredY_avg_ = PredY_avg;
PredY_avg_(excludelist) = [];

figure(25);clf(figure(25)); fig = gcf;
fig.Color = [1 1 1];
fig.Units = 'inches';
fig.Position = [ -12.2 8  12    4.5];

subplot(1,2,1)
%[r_t, p_t]=plotregressionwithSEM(PredY_avg, AHI_perPT(:,1));
scatter(PredY_avg, AHI_perPT(:,1), 50, 'filled'); hold on;
xlim([-0.05 1.55]); ylim([-5 105]); xlabel('Pred VE:VDrive');
ylabel('AHI'); title('AHI Vs Pred VEVdrive')
axis square

subplot(1,2,2)
scatter(Gtest_avg, AHI_perPT(:,1), 50, 'filled'); hold on;
xlim([-0.05 1.55]); ylim([-5 105]); xlabel('Actual VE:VDrive');
ylabel('AHI'); title('AHI Vs Actual VEVdrive')
axis square
suptitle(['AHI Vs Average VEVdrive']);

str = ['..\Figures\', 'AHI Vs Average VEVdrive_SleepOnly_PnasalTnT'];
%saveas(fig, str, 'png'); %savefig(str);

%% Proportion FL VS AHI
figure(26);clf(figure(26)); fig = gcf;
fig.Color = [1 1 1];
fig.Units = 'inches';
fig.Position = [1   2.5   12    4.5];

subplot(1,2,1)
%[r_t, p_t] = plotregressionwithSEM(ptsummary_withAp_n.Pred_percent, AHI_perPT(:,1));
scatter(ptsummary_withAp.Pred_percent2, AHI_perPT_(:,1),50,'filled','markerfacealpha',1);
% hndl_ls = lsline;  hndl_ls.LineStyle = '-'; hndl_ls.Color = [0 0 1];
% hndl_ref = refline(1,0); hndl_ref.LineStyle = '--'; hndl_ref.Color = [0.5 0.5 0.5];
xlim([-5 105]); xlabel('%FL (NumFL BB/TotalBB)');
ylim([-5 105]); ylabel('AHI');
title(['AHI Vs Pred %FL']);

subplot(1,2,2)
scatter(ptsummary_withAp.Actual_percent2, AHI_perPT_(:,1),50,'filled','markerfacealpha',1);
% hndl_ls = lsline;  hndl_ls.LineStyle = '-'; hndl_ls.Color = [0 0 1];
% hndl_ref = refline(1,0); hndl_ref.LineStyle = '--'; hndl_ref.Color = [0.5 0.5 0.5];
xlim([-5 105]); xlabel('%FL (NumFL BB/TotalBB)');
ylim([-5 105]); ylabel('AHI');
title(['AHI Vs Actual %FL']);

suptitlestr=['AHI Vs Proportion of sleep breaths classified as FL', withAPstr];
suptitle(suptitlestr);
str = ['..\Figures\', 'AHI Vs FL class_SleepOnly_PnasalTnT'];
%saveas(fig, str, 'png'); %savefig(str);

%% Novel metrics, Median VE:Vdrive during sleep and Time with severe obstruction during sleep
% aka, Proportion of breaths FL
figure(27); clf(figure(27)); fig = gcf;
fig.Color = [1 1 1];
fig.Units = 'inches';
%fig.Position = [-12.2   -3   12    4.5];

subplot(1,3,1);
if 1
    [r_1, p_1] = plotregressionwithSEM(PredY_avg_wAp, Gtest_avg_wAp);
else % without Ap BB
    [r_1, p_1] = plotregressionwithSEM(PredY_avg_, Gtest_avg_);
end
%hndl_ls = lsline;  hndl_ls.LineStyle = '-'; hndl_ls.Color = [0 0 1];
%hndl_ref = refline(1,0); hndl_ref.LineStyle = '--'; hndl_ref.Color = [0.5 0.5 0.5];
xlabel('Pnasal Predicted {\itflow:drive}');
ylabel('Gold Standard {\itflow:drive}');
axis square
titlestr = ['Patient Median {\itflow:drive}, r =', num2str(round(r_1,2))]; title(titlestr);

subplot(1, 3, 2);
if 1
    [r_2, p_2] = plotregressionwithSEM(ptsummary_withAp.Pred_percent, ptsummary_withAp.Actual_percent);
else % remove the one outlier
    [r_2, p_2] = plotregressionwithSEM(ptsummary_withAp.Pred_percent([1:3,5:17]), ptsummary_withAp.Actual_percent([1:3,5:17]));
end
%scatter(ptsummary_withAp.Pred_percent, ptsummary_withAp.Actual_percent,50,'filled','markerfacealpha',1);
%hndl_ls = lsline;  hndl_ls.LineStyle = '-'; hndl_ls.Color = [0 0 1];
%hndl_ref = refline(1,0); hndl_ref.LineStyle = '--'; hndl_ref.Color = [0.5 0.5 0.5];
xlim([-5 105]); xlabel('Pnasal Predicted %FL');
ylim([-5 105]); ylabel('Gold Standard %FL');
axis square
titlestr=['Proportion of FL breaths (0.5), r = ', num2str(round(r_2,2))]; title(titlestr);

subplot(1, 3, 3);
[r_3, p_3] = plotregressionwithSEM(ptsummary_withAp.Pred_percent2, ptsummary_withAp.Actual_percent2);
%scatter(ptsummary_withAp.Pred_percent, ptsummary_withAp.Actual_percent,50,'filled','markerfacealpha',1);
%hndl_ls = lsline;  hndl_ls.LineStyle = '-'; hndl_ls.Color = [0 0 1];
%hndl_ref = refline(1,0); hndl_ref.LineStyle = '--'; hndl_ref.Color = [0.5 0.5 0.5];
xlim([-5 105]); xlabel('Pnasal Predicted %FL');
ylim([-5 105]); ylabel('Gold Standard %FL');
axis square
titlestr=['Proportion of FL breaths (0.7), r = ', num2str(round(r_3,2))]; title(titlestr);

%suptitle(['Novel metrics (during sleep)']);

str = ['..\Figures\', 'Figure_6_']; %NovelMetrics_AutoRef2_SleepOnly_TnTPnasal'];
%saveas(fig, str, 'png'); %savefig(str);

if 0
    addpath('C:\Users\uqdmann\My Documents\MATLAB\random\BlandAltman');
    BlandAltman(PredY_avg, Gtest_avg);
end

%% AHI vs Novel Metrics
% grab from plot 25 and 26
figure(32);clf(figure(32)); fig = gcf;
fig.Color = [1 1 1];
fig.Units = 'inches';
%fig.Position = [ -12.2 8  12    4.5];
fig.Position = [ 2 3  12    4.5];

subplot(1,3,1)
[r_1, p_1]=plotregressionwithSEM(Gtest_avg_wAp', AHI_perPT_(:,1));
xlim([-0.05 1.55]); ylim([-5 105]); xlabel('Gold Standard {\itflow:drive}');
ylabel('AHI'); title(['Median {\itflow:drive} Vs AHI, r = -', num2str(round(r_1,2))]);
axis square

subplot(1,3,2)
%[r_2, p_2]=plotregressionwithSEM(ptsummary_withAp.Pred_percent, AHI_perPT_(:,1));
r_2 = corr(ptsummary_withAp.Actual_percent, AHI_perPT_(:,1));
scatter(ptsummary_withAp.Actual_percent, AHI_perPT_(:,1),50,'filled','markerfacealpha',1);
%hndl_ls = lsline;  hndl_ls.LineStyle = '-'; hndl_ls.Color = [0 0 1];
%hndl_ref = refline(1,0); hndl_ref.LineStyle = '--'; hndl_ref.Color = [0.5 0.5 0.5];
xlim([-5 105]); xlabel('%FL (NumFL BB/TotalBB)');
ylim([-5 105]); ylabel('AHI');
title(['AHI Vs Pred %FL']);
axis square

subplot(1,3,3)
[r_3, p_3]=plotregressionwithSEM(ptsummary_withAp.Actual_percent2, AHI_perPT_(:,1));
%r_2 = corr(ptsummary_withAp.Actual_percent, AHI_perPT_(:,1));
%scatter(ptsummary_withAp.Actual_percent, AHI_perPT_(:,1),50,'filled','markerfacealpha',1);
%hndl_ls = lsline;  hndl_ls.LineStyle = '-'; hndl_ls.Color = [0 0 1];
%hndl_ref = refline(1,0); hndl_ref.LineStyle = '--'; hndl_ref.Color = [0.5 0.5 0.5];
xlim([-5 105]); xlabel('Gold Standard %');
ylim([-5 105]); ylabel('AHI');
title(['Proportion of FL breaths (0.7) Vs AHI, r = ', num2str(round(r_3,2))]);
axis square

%suptitlestr=['AHI Vs Novel Metrics'];
%suptitle(suptitlestr);
str = ['..\Figures\', 'AHIvsNovelMetrics_AutoRef2_SleepOnly_PnasalTnT'];
%saveas(fig, str, 'png'); %savefig(str);

%% Proportion of breaths classified as FL (<threshold) Pred Vs Actual VEVdrive
if 0
    figure(28); clf(figure(28));
    plot([ones(length(ptsummary_withAp.Pred_percent),1), 2*ones(length(ptsummary_withAp.Actual_percent),1)], ...
        [ptsummary_withAp.Pred_percent, ptsummary_withAp.Actual_percent], 'bo'); hold on;
    plot([1 2], [ptsummary_withAp.Pred_percent, ptsummary_withAp.Actual_percent], 'color',[0.5 0.5 0.5], 'linewidth',0.5);
    plot([0.95 2.05], [nanmean(ptsummary_withAp.Pred_percent), nanmean(ptsummary_withAp.Actual_percent)], '--or', 'linewidth', 2.5);
    plot([0.95 2.05], [nanmedian(ptsummary_withAp.Pred_percent), nanmedian(ptsummary_withAp.Actual_percent)], '--om', 'linewidth', 2.5);
    xlim([0.75 2.25]); ylim([-5 105]);
    set(gca,'xtick',[1 2]);
    set(gca,'xticklabel',{'Pred', 'Actual'});
    ylabel('% of breaths FL'); xlabel('Method');
    titlestr = ['Proportion of sleep breaths classified as Flow Limited']; suptitle(titlestr);
    title('(Excluded Ap and LowFlow breaths not restored in this plot)');
    fig = gcf;
    fig.Color = [1 1 1]; % set background colour to white
    fig.Units = 'inches';
    %fig.Position = [-12 1 4.5 4.5];
    str = ['..\Figures\', titlestr, '_ManRef2_PredVsActual'];
    %saveas(fig, str, 'png'); %savefig(str);
    
    figure(29); clf(figure(29));
    plot([0.5*ones(length(ptsummary_withAp.Pred_percent),1)], ...
        [(ptsummary_withAp.Pred_percent - ptsummary_withAp.Actual_percent)], 'bo'); hold on;
    boxplot([(ptsummary_withAp.Pred_percent - ptsummary_withAp.Actual_percent)],...
        [(ones(length(ptsummary_withAp.Pred_percent),1))])
    
    xlim([0 1.5]);
    set(gca,'xtick',[0.5 1]);
    set(gca,'xticklabel',{'Indiv delta', 'summary box'});
    ylabel('Predicted % - Actual %');
    suptitle('Difference between Pred and Actual');
    fig = gcf;
    fig.Color = [1 1 1]; % set background colour to white
    fig.Units = 'inches';
    %fig.Position = [-12 1 4.5 4.5];
    str = ['..\Figures\', titlestr, '_ManRef2_PredVsActual_boxplot'];
    %saveas(fig, str, 'png'); %savefig(str);
    
    figure(30); clf(figure(30));
    scatter(ptsummary_withAp.Pred_FL, ptsummary_withAp.Actual_FL,50,'filled','markerfacealpha',1);
    hndl_ls = lsline;  hndl_ls.LineStyle = '-'; hndl_ls.Color = [0 0 1];
    hndl_ref = refline(1,0); hndl_ref.LineStyle = '--'; hndl_ref.Color = [0.5 0.5 0.5];
    xlim([0 4000]); xlabel('PredY, Number of FL breaths');
    ylim([0 4000]); ylabel('Actual, Number of FL breaths');
    suptitlestr=['Number of sleep breaths classified as FL'];
    suptitle(suptitlestr);
    [r,~] = corr(ptsummary_withAp.PredNumFL, ptsummary_withAp.ActualNumFL);
    titlestr=['FL if VE:Vdrive < ',num2str(thres), ', r= ', num2str(r)];
    title(titlestr);
    savestr = ['..\Figures\', suptitlestr, '_ManRef2_Scatter'];
    fig = gcf;
    %saveas(fig, savestr, 'png'); %savefig(str);
end


%% Median VEVdrive VS Proportion FL
Gtest_avg_ = Gtest_avg(~isnan(Gtest_avg));
PredY_avg_ = PredY_avg(~isnan(PredY_avg));
figure(31); clf(figure(31)); fig = gcf;
fig.Color = [1 1 1];
fig.Units = 'inches';
fig.Position = [-12.2   -3   12    4.5];

subplot(1,2,1)
%[r_t, p_t] = plotregressionwithSEM(PredY_avg_, ptsummary_withAp.Pred_percent);
scatter(PredY_avg_, ptsummary_withAp.Pred_percent,50,'filled','markerfacealpha',1);
hndl_ls = lsline;  hndl_ls.LineStyle = '-'; hndl_ls.Color = [0 0 1];
%hndl_ref = refline(1,0); hndl_ref.LineStyle = '--'; hndl_ref.Color = [0.5 0.5 0.5];
xlim([-0.05 1.55]); xlabel('VEVdrive');
ylim([-5 105]); ylabel('%FL (NumFL BB/TotalBB)');
titlestr=['Pred %FL Vs Pred VEVdrive']; title(titlestr);

subplot(1,2,2)
scatter(Gtest_avg_, ptsummary_withAp.Actual_percent,50,'filled','markerfacealpha',1);
hndl_ls = lsline;  hndl_ls.LineStyle = '-'; hndl_ls.Color = [0 0 1];
%hndl_ref = refline(1,0); hndl_ref.LineStyle = '--'; hndl_ref.Color = [0.5 0.5 0.5];
xlim([-0.05 1.55]); xlabel('VEVdrive');
ylim([-5 105]); ylabel('%FL (NumFL BB/TotalBB)');
titlestr=['Actual %FL Vs Actual VEVdrive']; title(titlestr);

suptitlestr=['%FL Vs Average VEVdrive'];
suptitle(suptitlestr);


%%
if 0
    figure(33); clf(figure(33));
    plot([ones(length(ptsummary_withAp.Pred_percent),1), 2*ones(length(ptsummary_withAp.Clin_percent),1)], ...
        [ptsummary_withAp.Pred_percent, ptsummary_withAp.Clin_percent], 'bo'); hold on;
    plot([1 2], [ptsummary_withAp.Pred_percent, ptsummary_withAp.Clin_percent], 'color',[0.5 0.5 0.5], 'linewidth',0.5);
    plot([0.95 2.05], [nanmean(ptsummary_withAp.Pred_percent), nanmean(ptsummary_withAp.Clin_percent)], '--or', 'linewidth', 2.5);
    plot([0.95 2.05], [nanmedian(ptsummary_withAp.Pred_percent), nanmedian(ptsummary_withAp.Clin_percent)], '--om', 'linewidth', 2.5);
    xlim([0.75 2.25]); ylim([-5 105]);
    set(gca,'xtick',[1 2]);
    set(gca,'xticklabel',{'Pred FL', 'Clinical'});
    ylabel('% of breaths FL'); xlabel('Method');
    titlestr = ['Proportion of sleep breaths classified as Flow Limited']; suptitle(titlestr);
    title('(Excluded Ap and LowFlow breaths are restored in this plot)');
    fig = gcf;
    fig.Color = [1 1 1]; % set background colour to white
    fig.Units = 'inches';
    %fig.Position = [-12 1 4.5 4.5];
    str = ['..\Figures\', titlestr, '_ManRef2_withAddins'];
    %saveas(fig, str, 'png'); %savefig(str);
end


%% Load subject waveform data for VEVdrive, Flow and Edi plot
predy_flow(predy_flow<0)=0;
predy_flow(predy_flow>maxG)=1.5;

% read spreadsheet
AnalyzeDataSpreadsheet = 'AnalyzeDataSpreadsheet.xlsx';
[num,patients,~] = xlsread(AnalyzeDataSpreadsheet,1,'F3:G56');
for n=1 %41:54
    if ~ismember(n, settings.Pnasal_list)
        continue
    end
    clearvars Edi Flow StarttimeSpike
    % Plot FL values over time for an example subject (e.g. compare with Spike data). Random check (n=8, 1313) Looks great.
    % n=3; % 1313 is now 14
    % [2;3;4;5;6;8;9;10;11;14;16;17;18;19;20;21;22;23;24;26;27;28;29;30;32;33;34;35;36;38;39;40;41;43;44;45;46;47;50;53;54]
    studyname = char(patients{n});
    %load('J:\PEOPLE\FACULTY\SANDS\Phenotyping from Pes and Flow\SpikeFiles\1313.mat', 'Edi', 'Flow');
    load(['C:\PSG_Data\FlowDrive\SourceMat 20171123\',studyname], 'Edi', 'Flow', 'StarttimeSpike');
    if ~(exist('StarttimeSpike', 'var') == 1)
        StarttimeSpike = 0;
    end
    Time = StarttimeSpike:0.008:StarttimeSpike+(length(Flow.values)-1)*0.008;
    
    PTtime = PtData.BB_time(PtData.PT==n);
    PTtime = PTtime-PTtime(1)+StarttimeSpike;
    Data1 = [PtData.BB_time(PtData.PT==n) PtData.BB_Ttot(PtData.PT==n) predy_flow(PtData.PT==n) Gtest_All(PtData.PT==n)];
    addNaNgaps=1;
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
    % if Flow.length<=length(Time)
    %     Flow.values(end:length(Time))=NaN;
    % else
    %     Flow.values(length(Time):end)=[];
    % end
    if Edi.length<=length(Time)
        Edi.values(end:length(Time))=NaN;
    else
        Edi.values(length(Time):end)=[];
    end
    
    % Data columns: (1) Time, (2) BBTtot, (3) PredY, (4) Gtest.
    
    %% make the plot for the subject data
    figure(200+n); clf(200+n) %ax(1) = subplot(2,1,1);
    ax(1)=subplot(3,1,1);
    stairs(Data1(:,1),100*Data1(:,4),'k'); % actual
    hold('on')
    stairs(Data1(:,1),100*Data1(:,3),'r'); % pred
    plot([Data1(1,1) Data1(end,1)],[0 0],'k:');
    plot([Data1(1,1) Data1(end,1)],100*[1 1],'k:');
    plot([Data1(1,1) Data1(end,1)],50*[1 1],'k:');
    set(gca,'xtick',[],'box','off');
    ylabel('{\itflow:drive} (%)'); ylim([-5 105]);
    
    dsf=5; dt=Flow.interval;
    FlowF=Flow.values;
    if 1
        filter_HFcutoff_butter0 = 12.5;
        filter_order0 = 1;
        [B_butter0,A_butter0] = butter(filter_order0,[filter_HFcutoff_butter0]/(1/dt/2),'low');
        FlowF = filtfilt(B_butter0,A_butter0,FlowF); %filtfilt, otherwise flow signal is right-shifted
    end
    %ax(2)=subplot(3,1,2); plot(Time,Flow.values,'g'); hold('on');
    ax(2)=subplot(3,1,2); plot(downsample(Time,dsf),downsample(FlowF,dsf),'k');
    set(gca,'xtick',[],'box','off');
    ylabel('Flow (L/s)');
    
    ax(3)=subplot(3,1,3); plot(downsample(Time,dsf),downsample(Edi.values,dsf)); hold on;
    box('off');
    ylabel('Edi (uV)');
    
    linkaxes(ax,'x');
    suptitle([studyname]);
end

%% mods to plot
plot(StarttimeSpike+5400, 0, 'r^');
plot(StarttimeSpike+23620, 0, 'r^');

xlim([StarttimeSpike inf])
xlim([StarttimeSpike+5200 StarttimeSpike+6400]) % old indices
xlim([StarttimeSpike+5400 StarttimeSpike+5800]) % old indices
xlim([StarttimeSpike+23560 StarttimeSpike+23680]) % old indices

