%%
load('D:\Dropbox (Partners HealthCare)\MAD-OX snore\Analyzed\EventAnalyzed\MADOXBaseline_All.mat')
% load('C:\Users\bwhha\Dropbox (Partners HealthCare)\MAD-OX snore\Analyzed\EventAnalyzed\MADOXBaseline_All.mat')

notalldata = true;
%%
KeepIdx = true(size(BreathDataTableAll,1),1);
if notalldata
    % Evaluate relationship between Pes effort, ventilation and snoring
    KeepIdx = BreathDataTableAll.hypnog_B <= 4 & (BreathDataTableAll.Etype == 4 | BreathDataTableAll.Etype == 0);
    BreathDataTableAll = BreathDataTableAll(KeepIdx,:);
    BreathFLDataTableAll = BreathFLDataTableAll(KeepIdx,:);
    BreathSnoreAll = BreathSnoreAll(KeepIdx,:);
end

%%
% Initialize ands sample variables
PesEffort = BreathDataTableAll.DeltaPes;
VI = BreathDataTableAll.VI;
WinVar = BreathSnoreAll.WinVarMax_i;
SnoredB = 10.*(log10(WinVar./((0.00002).^2)));
quantiles = 5:5:95;
quantiles2 = 5:5:100;

%% Correlations
[r1,p1] = corr(PesEffort, SnoredB,'rows','pairwise','Type','Spearman')
[r2,p2] = corr(VI, SnoredB,'rows','pairwise','Type','Spearman')
[r3,p3] = corr(VI./PesEffort, SnoredB,'rows','pairwise','Type','Spearman')

%% Pes vs Snore lookup without ventilation
if 1
    % Identify cutoffs using prctiles
    PesEffortCut = prctile(PesEffort, quantiles);
    VICut = prctile(VI, quantiles);

    % Seperate input variables into their appropriate prctile
    IDa = [sum(VI > VICut,2), ...
           sum(PesEffort > PesEffortCut,2)];
   IDa(isnan(PesEffort) | isnan(VI)) = nan;

    % Find the mean Snore db for each prctile
    PredSnoredB = nan(length(PesEffortCut), 1);
    SSE = nan(length(PesEffortCut));
    Ncount = 57;
    SampledSnoredB = nan(size(SnoredB));
    for j = 0:length(PesEffortCut)
        % sample snore data in bin 
        SampledSnore = SnoredB(IDa(:,2)==j);

        % check that count > Ncount
        if length(SampledSnore) < Ncount
            PredSnoredB(j+1) = nan;
            SSE(j+1) = nan;
            continue
        end
        SampledSnoredB(IDa(:,2)==j) = SampledSnore;
        % calculate median snore
        PredSnoredB(j+1) = nanmedian(SampledSnore);

        % Compute sum of squared differences relative to predicted
        SSE(j+1) = nansum((SampledSnore - PredSnoredB(j+1)).^2);
    end

    SSEall = nansum(nansum(SSE));
    % SST = nansum((SnoredB - nanmedian(SnoredB)).^2);
    SST = nansum((SampledSnoredB - nanmean(SampledSnoredB)).^2);

    Rsq = 1-(SSEall/SST)
end

%% Assess how well Pes and Ventilation explain snoring
% Identify cutoffs using prctiles
PesEffortCut = prctile(PesEffort, quantiles);
VICut = prctile(VI, quantiles);

% Seperate input variables into their appropriate prctile
IDa = [sum(VI > VICut,2), ...
       sum(PesEffort > PesEffortCut,2)];
IDa(isnan(PesEffort) | isnan(VI)) = nan;

% Find the mean Snore db for each prctile
PredSnoredB = nan(length(PesEffortCut), length(VICut));
SSE = nan(length(PesEffortCut), length(VICut));
RsqPerVI = nan(length(VICut),1);
Ncount = 57;
SampledSnoredB = nan(size(SnoredB));
for i = 0:length(VICut)
    for j = 0:length(PesEffortCut)
        % sample snore data in bin 
        SampledSnore = SnoredB(IDa(:,1)==i & IDa(:,2)==j);
                        
        % check that count > Ncount
        if length(SampledSnore) < Ncount
            PredSnoredB(i+1,j+1) = nan;
            SSE(i+1,j+1) = nan;
            continue
        end
        SampledSnoredB(IDa(:,1)==i & IDa(:,2)==j) = SampledSnore;
        % calculate median snore
        PredSnoredB(i+1,j+1) = nanmedian(SampledSnore);
        
        % Compute sum of squared differences relative to predicted
        SSE(i+1,j+1) = nansum((SampledSnore - PredSnoredB(i+1,j+1)).^2);
    end
    % Fit per ventilation prctile
    SampledSnorePerVI = SampledSnoredB(IDa(:,1)==i);
    SSEperVI = nansum(SSE(i+1,:));
    SSTperVI = nansum((SampledSnorePerVI - nanmean(SampledSnorePerVI)).^2);
    RsqPerVI(i+1) = 1 - SSEperVI/SSTperVI;
end

SSEall = nansum(nansum(SSE));
% SST = nansum((SnoredB - nanmedian(SnoredB)).^2);
SST = nansum((SampledSnoredB - nanmean(SampledSnoredB)).^2);

Rsq = 1-(SSEall/SST)

%% Plot figures
% VIEdges = linspace(0,ceil(max(VI)),bins);
VIEdges = prctile(VI,0:5:100);
discTemp = discretize(VI,VIEdges);
discTemp(isnan(discTemp)) = length(VIEdges)+1;
VIEdgesNaN = [VIEdges,NaN];
VIdisc = VIEdgesNaN(discTemp)';
VIdisc = round(VIdisc,2);

% VdriveNormEdges = linspace(0,ceil(max(VdriveNorm)),bins);
PesEffortEdges = prctile(PesEffort,0:5:100);
discTemp = discretize(PesEffort,PesEffortEdges);
discTemp(isnan(discTemp)) = length(PesEffortEdges)+1;
PesEffortEdgesNaN = [PesEffortEdges,NaN];
PesEffortDisc = PesEffortEdgesNaN(discTemp)';
PesEffortDisc = round(PesEffortDisc,1);
tbl = table(VIdisc,PesEffortDisc,SampledSnoredB);

figure('Position', [1686 546 1672 420])
h = heatmap(tbl,'PesEffortDisc','VIdisc','ColorVariable','SampledSnoredB','ColorMethod','Median');
ylabel('Ventilation (fraction of eupnea)'); xlabel('Esophageal pressure swing (cmH2O)')
h.YDisplayData = flipud(h.YDisplayData);  % equivalent to 'YDir', 'Reverse'
h.CellLabelFormat = '%0.1f';

% Plot snore vs drive per ventilation

N = 20;
C = repmat(linspace(0.75,0.2,N).',1,3);

figure, hold on
for ii = 1:2:length(PesEffortCut)
    SnorePlotVar = PredSnoredB(ii,1:end-1);
    plot(PesEffortCut, SnorePlotVar, '-o', 'Color', C(ii,:),...
        'LineWidth',2,'MarkerFaceColor',C(ii,:), 'MarkerSize', 2)
end
xlabel('Esophageal Pressure (cmH2O)')
ylabel('Snore Power (dB)')
hold off

%% It is more practical to build a model that predicts esophgeal pressure

% Identify cutoffs using prctiles
SnoreCut = prctile(SnoredB, quantiles);
VICut = prctile(VI, quantiles);

% Seperate input variables into their appropriate prctile
IDa = [sum(VI > VICut,2), ...
       sum(SnoredB > SnoreCut,2)];
IDa(isnan(SnoredB) | isnan(VI)) = nan;

% Find the mean Snore db for each prctile
PredEffort = nan(length(SnoreCut), length(VICut));
SSE = nan(length(SnoreCut), length(VICut));
RsqPerVI = nan(length(VICut),1);
PredEffortLong = nan(length(PesEffort),1);
SampledEffortAll = nan(size(PesEffort));
Ncount = 57;
for i = 0:length(VICut)
    for j = 0:length(SnoreCut)
        % sample snore data in bin 
        SampledEffort = PesEffort(IDa(:,1)==i & IDa(:,2)==j);
                        
        % check that count > Ncount
        if length(SampledEffort) < Ncount
            PredEffort(i+1,j+1) = nan;
            SSE(i+1,j+1) = nan;
            continue
        end
        SampledEffortAll(IDa(:,1)==i & IDa(:,2)==j) = SampledEffort;
        
        % calculate median snore
        PredEffort(i+1,j+1) = nanmedian(SampledEffort);
        PredEffortLong(IDa(:,1)==i & IDa(:,2)==j) = PredEffort(i+1,j+1);
        
        % Compute sum of squared differences relative to predicted
        SSE(i+1,j+1) = nansum((SampledEffort - PredEffort(i+1,j+1)).^2);
        
        % sample data in the bin
        SampledEffort = PesEffort(IDa(:,1)==i & IDa(:,2)==j);
    end
    % Fit per ventilation prctile
    SampledEffortPerVI = SampledEffortAll(IDa(:,1)==i);
    SSEperVI = nansum(SSE(i+1,:));
    SSTperVI = nansum((SampledEffortPerVI - nanmean(SampledEffortPerVI)).^2);
    RsqPerVI(i+1) = 1 - SSEperVI/SSTperVI;
end

SSEall = nansum(nansum(SSE));
SST = nansum((SampledEffortAll - nanmean(SampledEffortAll)).^2);

Rsq = 1-(SSEall/SST)


%% Plot figures
% VIEdges = linspace(0,ceil(max(VI)),bins);
VIEdges = prctile(VI,quantiles2);
discTemp = discretize(VI,VIEdges);
discTemp(isnan(discTemp)) = length(VIEdges)+1;
VIEdgesNaN = [VIEdges,NaN];
VIdisc = VIEdgesNaN(discTemp)';
VIdisc = round(VIdisc,2);

% VdriveNormEdges = linspace(0,ceil(max(VdriveNorm)),bins);
SnoreEdges = prctile(SnoredB,quantiles2);
discTemp = discretize(SnoredB,SnoreEdges);
discTemp(isnan(discTemp)) = length(SnoreEdges)+1;
SnoreEdgesNaN = [SnoreEdges,NaN];
SnoreDisc = SnoreEdgesNaN(discTemp)';
SnoreDisc = round(SnoreDisc,1);
tbl = table(VIdisc,SnoreDisc,SampledEffortAll);

figure('Position', [1686 546 1672 420])
h = heatmap(tbl,'SnoreDisc','VIdisc','ColorVariable','SampledEffortAll','ColorMethod','Median');
ylabel('Ventilation (Fraction of Eupnea)'); xlabel('Snore Power (dB)')
h.YDisplayData = flipud(h.YDisplayData);  % equivalent to 'YDir', 'Reverse'
h.CellLabelFormat = '%0.1f';

% Plot snore vs drive per ventilation

N = 20;
C = repmat(linspace(0.75,0.2,N).',1,3);

figure, hold on
for ii = 1:2:length(SnoreCut)
    EffortPlotVar = PredEffort(ii,1:end-1);
    plot(SnoreCut, EffortPlotVar, '-o', 'Color', C(ii,:),...
        'LineWidth',2,'MarkerFaceColor',C(ii,:), 'MarkerSize', 2)
end
xlabel('Snore Power')
ylabel('Esophageal Pressure (Fraction of Eupneic Pressure)')
hold off


%% 
Subject = BreathDataTableAll.Subject;
tbl = table(PesEffort, SnoredB, VIdisc, Subject);
modelspec = 'SnoredB ~ PesEffort + VIdisc';
mdllm = fitlm(tbl,modelspec)
PredPesLM = predict(mdllm);

% [rpc,fig,sstruct] = BlandAltman(PesEffort,PredEffortLME,{'True Pes Effort','Predicted Pes Effort'});

Subject = BreathDataTableAll.Subject;
tbl2 = table(PesEffort, SampledSnoredB, VIdisc, Subject);
modelspec2 = 'SampledSnoredB ~ PesEffort + VIdisc';
mdllm2 = fitlm(tbl2,modelspec2)
PredPesLM2 = predict(mdllm2);

% [rpc,fig,sstruct] = BlandAltman(PesEffort,PredEffortLME,{'True Pes Effort','Predicted Pes Effort'});
%% Generate time series
subject = '1922N0';
load(['D:\Dropbox (Partners HealthCare)\MAD-OX snore\Converted\',subject,'_XHz.mat']);
% load(['C:\Users\bwhha\Dropbox (Partners HealthCare)\MAD-OX snore\Converted\',subject,'_XHz.mat']);

Fs = 125;
SubIdx = strcmp(BreathDataTableAll.Subject, subject);
BreathDataTableSub = BreathDataTableAll(SubIdx,:);
PredEffortSub = PredEffortLong(SubIdx,:);
TruePesEffortSub = PesEffort(SubIdx,:);
VIsub = VI(SubIdx,:);
SnoreSub = SampledSnoredB(SubIdx,:);
Time = DataEventHypnog_Mat(:,strcmp(ChannelsList,'Time'));

if notalldata
    KeepIdxSub = KeepIdx(SubIdx);
end

PredPesEffort = nan(size(DataEventHypnog_Mat,1),1);
TruePesEffort = nan(size(DataEventHypnog_Mat,1),1);
VISub_ts = nan(size(DataEventHypnog_Mat,1),1);
SnoreIntensity = nan(size(DataEventHypnog_Mat,1),1);

for brNum = 1:size(BreathDataTableSub,1)
    % Store Ventilation (relative to VEupnea) in time domain
    startIdx = round(BreathDataTableSub.Time_start(brNum)*Fs) - round(Time(1)*Fs);
    endIdx = round(BreathDataTableSub.Time_end(brNum)*Fs) - round(Time(1)*Fs);
    PredPesEffort(startIdx:endIdx,1) =  PredEffortSub(brNum);   
    TruePesEffort(startIdx:endIdx,1) =  TruePesEffortSub(brNum);
    VISub_ts(startIdx:endIdx,1) = VIsub(brNum);
    SnoreIntensity(startIdx:endIdx,1) = SnoreSub(brNum);
end
DataEventHypnog_Mat(:,end+1:end+2) = [PredPesEffort, TruePesEffort];
ChannelsList(end+1:end+2) = {'PredPes','TruePes'};
DataEventHypnog_Mat(:,end+1) = 10.*(log10((DataEventHypnog_Mat(:,strcmp(ChannelsList, 'WinVar')))./((0.00002).^2)));
ChannelsList(end+1) = {'SnoredB'};

% [rpc,fig,sstruct] = BlandAltman(PredEffortSub,TruePesEffortSub,{'True Pes Effort','Predicted Pes Effort'});

% Ensemble average hypopneas
BreathDataTableSub = NumberApneaHypopneaBreaths(BreathDataTableSub);

% Resample VE_td
Fs_RS = 4;
newTime = Time(1):(1/4):Time(end);
PredPes_RS = interp1(Time, PredPesEffort, newTime, 'previous');
TruePes_RS = interp1(Time, TruePesEffort, newTime, 'previous');
VE_tdRS = interp1(Time, VISub_ts, newTime, 'previous');
SnoreI_RS = interp1(Time, SnoreIntensity, newTime, 'previous');

% Sample data for next set of operations
termEvtIdx = find(BreathDataTableSub.HypopNum == 1)'; % Find terminal breath of hypopnea
brCount = 0;
dTime = 150;

TruePes_AllSub = nan(length(termEvtIdx), dTime*Fs_RS*2+1);
PredPes_AllSub = nan(length(termEvtIdx), dTime*Fs_RS*2+1);
VE_AllSub = nan(length(termEvtIdx), dTime*Fs_RS*2+1);
SnoreI_AllSub = nan(length(termEvtIdx), dTime*Fs_RS*2+1);

disp(['Total number of respiratory events: ' num2str(length(termEvtIdx))])

for evtBreath = termEvtIdx
    brCount = brCount+1;

    % Get event breath indices and check size
    evtBrthIdx = round(BreathDataTableSub.Time_end(evtBreath)*Fs_RS) - round(newTime(1)*Fs_RS);

    % Conditions for setting pre and post 
    if BreathDataTableSub.Time_end(evtBreath) - BreathDataTableSub.Time_start(1) < dTime %start condition
        idxPreEvt = round(BreathDataTableSub.Time_end(evtBreath)*Fs_RS) - ...
            (round(BreathDataTableSub.Time_end(evtBreath)*Fs_RS) - ...
            round(BreathDataTableSub.Time_start(1)*Fs_RS)) -...
            round(newTime(1)*Fs_RS);
        idxPostEvt = round(BreathDataTableSub.Time_end(evtBreath)*Fs_RS) + ...
            round(dTime*Fs_RS) - round(newTime(1)*Fs_RS);
    elseif BreathDataTableSub.Time_end(end) - BreathDataTableSub.Time_start(evtBreath) < dTime %end condition
        idxPreEvt = round(BreathDataTableSub.Time_end(evtBreath)*Fs_RS) - ...
            round(dTime*Fs_RS) - round(newTime(1)*Fs_RS);
        idxPostEvt = round(BreathDataTableSub.Time_end(evtBreath)*Fs_RS) + ...
            (round(BreathDataTableSub.Time_end(end)*Fs_RS) - ...
            round(BreathDataTableSub.Time_start(evtBreath)*Fs_RS)) - ...
            round(newTime(1)*Fs_RS);
    else
        idxPreEvt = round(BreathDataTableSub.Time_end(evtBreath)*Fs_RS) - ...
            round(dTime*Fs_RS) - round(newTime(1)*Fs_RS);
        idxPostEvt = round(BreathDataTableSub.Time_end(evtBreath)*Fs_RS) + ...
            round(dTime*Fs_RS) - round(newTime(1)*Fs_RS); % Time_end of terminal event at 601 (dTime*Fs_RS + 1)
    end

%     if 0
%         figure('Position', [54 545 560 420]), plot(VE_tdRS(idxPreEvt:idxPostEvt)), hold on
%         plot(601, VE_tdRS(idxPreEvt+600), 'o'), hold off
%     end

    % adjust event breath location
    evtBreathMagn = VE_tdRS(evtBrthIdx);
    if evtBreathMagn > 0.9 % bring event breath back to before arousal breath
        newEvtIdxRel = find(VE_tdRS(idxPreEvt:idxPreEvt+dTime*Fs_RS) < 0.95, 1, 'last');
        EvtIdxDiff = newEvtIdxRel - (dTime*Fs_RS + 1);
    else % bring event breath to breath before arousal breath
        tempIdx = find(VE_tdRS(idxPreEvt+dTime*Fs_RS:idxPostEvt) > 0.95,1,'first');
        newEvtIdxRel = tempIdx - 1;
        EvtIdxDiff = newEvtIdxRel;
    end

    %
    newEvtIdx = evtBrthIdx + EvtIdxDiff-1;
    newIdxPreEvt = idxPreEvt + EvtIdxDiff-1;
    newIdxPostEvt = idxPostEvt + EvtIdxDiff-1;

    if newIdxPostEvt > length(VE_tdRS)
        newIdxPostEvt = length(VE_tdRS);
    elseif newIdxPreEvt < 0
        newIdxPreEvt = 1;
    end

    % shift pre/post event index such that event is always at 601
    evtBrthIdx2 = dTime*Fs_RS + 1;
    idxPreEvt2 = evtBrthIdx2 - (newEvtIdx - newIdxPreEvt); 
    idxPostEvt2 = evtBrthIdx2 + (newIdxPostEvt - newEvtIdx);

    % Store Pes
    TruePesTemp = nan(1, length(1:dTime*Fs_RS*2+1));
    TruePesTemp(idxPreEvt2:idxPostEvt2) = TruePes_RS(newIdxPreEvt:newIdxPostEvt);
    TruePes_AllSub(brCount,1:dTime*Fs_RS*2+1) = TruePesTemp;
    
    PredPesTemp = nan(1, length(1:dTime*Fs_RS*2+1));
    PredPesTemp(idxPreEvt2:idxPostEvt2) = PredPes_RS(newIdxPreEvt:newIdxPostEvt);
    PredPes_AllSub(brCount,1:dTime*Fs_RS*2+1) = PredPesTemp;
    
    VEtemp = nan(1, length(1:dTime*Fs_RS*2+1));
    VEtemp(idxPreEvt2:idxPostEvt2) = VE_tdRS(newIdxPreEvt:newIdxPostEvt);
    VE_AllSub(brCount,1:dTime*Fs_RS*2+1) = VEtemp;
    
    SnoreItemp = nan(1, length(1:dTime*Fs_RS*2+1));
    SnoreItemp(idxPreEvt2:idxPostEvt2) = SnoreI_RS(newIdxPreEvt:newIdxPostEvt);
    SnoreI_AllSub(brCount,1:dTime*Fs_RS*2+1) = SnoreItemp;


    if 0
        figure('Position', [630 545 560 420]), plot(VEtemp(idxPreEvt2:idxPostEvt2)), hold on
        plot(601, VEtemp(601), 'o'), hold off
    end
    close all
    
end

% Plot ensemble average
newTimeRS = -dTime:(1/Fs_RS):dTime; %size(VEStruct.(subnames{1}),2)/Fs_RS-dtRS;
VEToPlot = 100.*VE_AllSub;
TruePesToPlot = TruePes_AllSub;
PredPesToPlot = PredPes_AllSub;
SnoreToPlot = SnoreI_AllSub;

ax1 = subplot(5,1,1); hold on
ax2 = subplot(5,1,2); hold on
ax3 = subplot(5,1,3); hold on
ax4 = subplot(5,1,4); hold on
ax5 = subplot(5,1,5); hold on
set(gcf, 'Position', [1112 42 560 924])

for ii = 1:size(VEToPlot,1)
    % Plot VE
    plot(ax1, newTimeRS, VEToPlot(ii,:),'-','Color', [0,0,0,0.05])

    plot(ax2, newTimeRS, SnoreToPlot(ii,:),'-','Color', [0,0,0,0.05])

    plot(ax3, newTimeRS, TruePesToPlot(ii,:),'-','Color', [0,0,0,0.05])
    
    plot(ax4, newTimeRS, PredPesToPlot(ii,:),'-','Color', [0,0,0,0.05])
end

%Plot mean VE
subplot(ax1)
VEToPlot_Avg = nanmean(VEToPlot, 1);
plot(newTimeRS, VEToPlot_Avg, 'k-','LineWidth',2)
% xlabel('Time (s)')
ylabel('Ventilation (% Eupnea)')
% title('Ventila')
ylim([-10, 250])
yticks([0 100 200 300 400])
set(gca, 'FontSize', 10)
hold off

%Plot mean snore
subplot(ax2)
SnoreToPlot_Avg = nanmean(SnoreToPlot, 1);
plot(newTimeRS, SnoreToPlot_Avg, 'k-','LineWidth',2)
% xlabel('Time (s)')
ylabel('Mean Snore I (dB)')
% title('Ventila')
% ylim([-10, 250])
% yticks([0 100 200 300 400])
set(gca, 'FontSize', 10)
hold off

%Plot mean True Pes
subplot(ax3)
TruePes_Avg = nanmean(TruePesToPlot, 1);
plot(newTimeRS, TruePes_Avg, 'k-','LineWidth',2)
% xlabel('Time (s)')
% ylabel('True Pes (cmH_2O)')
ylabel('True Flow Drive')
%     title(subnames{subnum})
ylim([0, max(TruePes_Avg)*1.5])
% yticks([0 100 200 300 400])
set(gca, 'FontSize', 10)
hold off

%Plot mean Pred Pes
subplot(ax4)
PredPes_Avg = nanmean(PredPesToPlot, 1);
plot(newTimeRS, PredPes_Avg, 'k-','LineWidth',2)
% xlabel('Time (s)')
% ylabel('Pred Pes (cmH_2O)')
ylabel('Pred Flow Drive')
%     title(subnames{subnum})
ylim([0, max(PredPes_Avg)*1.5])
% yticks([0 100 200 300 400])
set(gca, 'FontSize', 10)
hold off

%Plot mean True and Pred Pes on single plot
subplot(ax5), hold on
plot(newTimeRS, TruePes_Avg, 'r-','LineWidth',2)
plot(newTimeRS, PredPes_Avg, 'b-.','LineWidth',2)
xlabel('Time (s)')
ylabel('True vs. Pred Pes')
% yticks([0 100 200 300 400])
set(gca, 'FontSize', 10)
hold off
