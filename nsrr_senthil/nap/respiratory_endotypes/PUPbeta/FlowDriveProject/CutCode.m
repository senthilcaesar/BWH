%% average Gtest vs average PredY for each pt (medians...)
figure(24); clf(figure(24)); fig = gcf;
fig.Color = [1 1 1]; 
fig.Units = 'inches';
%fig.Position = [ -12.2 -2.5  12.0000    4.5000];

subplot(1, 2, 1);
Gtest_avg = NaN(54,1);
PredY_avg = NaN(54,1);
Gtest_avg_wAp = NaN(54,1);
PredY_avg_wAp = NaN(54,1);

for i=ftrnum % 1:size(predyL1O_array,2) say at 20 ftrs
    PT_list = unique(PtData.PT);
    for subj=1:54 %length(PT_list)
        if ~ismember(subj, PT_list)            
            continue
        end
        %Isubj = PtData.PT==subj; % have not used weights(Isubj) ??
        Isubj=(PtData.PT==subj)&(PtData.Hypnog<4)&(PtData.Ar==0); % sleep only
        Gtest_avg(subj) = median(Gtest_All(Isubj));
        PredY_avg(subj) = median(predyL1O_array(Isubj,i));
        
        % make the data that includes the excluded Ap breaths
        indx = find(ptsummary_withAp.PT==subj);
        numPredFL = ptsummary_withAp.Pred_FL(indx);
        numActualFL = ptsummary_withAp.Actual_FL(indx);
        
        Gtest_avg_wAp(subj) = median([Gtest_All(Isubj);ones(numActualFL,1)*0.1]);
        PredY_avg_wAp(subj) = median([predyL1O_array(Isubj,i);ones(numPredFL,1)*0.1]);
       
        scatter(predyL1O_array((Isubj),i),Gtest_All(Isubj),2,'filled','markerfacealpha',0.5); hold on;
    end
end

scatter(PredY_avg, Gtest_avg, 50, 'filled','markerfacealpha',1); hold on;
scatter(PredY_avg_wAp, Gtest_avg_wAp, 50, 'filled','markerfacealpha',1); hold on;
% lsline;
xlim([0 1.5]); xlabel('PredY VE:VDrive'); 
ylim([0 1.5]); ylabel('Actual VE:VDrive');
title('(Sleep only)Indiv breaths and patient medians');
axis square

subplot(1, 2, 2);
[r_t, p_t] = plotregressionwithSEM(PredY_avg_wAp, Gtest_avg_wAp);
%scatter(PredY_avg, Gtest_avg, 50, 'filled'); hold on; 
%hndl_ls = lsline;  hndl_ls.LineStyle = '-'; hndl_ls.Color = [0 0 1];
%hndl_ref = refline(1,0); hndl_ref.LineStyle = '--'; hndl_ref.Color = [0.5 0.5 0.5];
xlim([0 1.5]); xlabel('PredY VE:VDrive');
ylim([0 1.5]); ylabel('Actual VE:Vdrive');
axis square
[r, ~] = corr(PredY_avg_wAp(~isnan(PredY_avg_wAp)), Gtest_avg_wAp(~isnan(PredY_avg_wAp)));
titlestr = ['(Sleep only) Patient medians, r=', num2str(r)]; title(titlestr);

suptitle(['Actual VEVdrive Vs Pred VEVdrive']);
str = ['..\Figures\', 'MedianGtestVsMedianPredY_AutoRef2_SleepOnly'];
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
PT_list = unique(PtData.PT);

sleeponly = 1;
addbackAP = 1;
thres1 = 0.5;
thres2 = 0.7;

numBBinTest = NaN(54,4);

for subj=1:54
    if ~ismember(subj, PT_list)
        continue
    end
    if sleeponly % make Isubj
        Isubj=(PtData.PT==subj)&(PtData.Hypnog<4)&(PtData.Ar==0);
%         numBBinTest(subj,3) = nnz(Isubj);
%         Isubj=(PtData.PT==subj)&(PtData.Hypnog<4)&(PtData.NotAr==1);
%         numBBinTest(subj,4) = nnz(Isubj);
    else
        Isubj=(PtData.PT==subj);
%         numBBinTest(subj,1) = nnz(Isubj);
%         Isubj=(PtData.PT==subj)&(PtData.Hypnog<4);
%         numBBinTest(subj,2) = nnz(Isubj);
    end
    

    clinscore = PtData.Etype(Isubj);
    numclinscored = nnz(clinscore==2|clinscore==4); % obstructive and hypopnoea
    predBelowThres = predyL1O_array(Isubj,ftrnum)<thres1;
    numPredBBbelowthreshold = nnz(predBelowThres);
    actualBelowThres = Gtest_All(Isubj)<thres1;
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

