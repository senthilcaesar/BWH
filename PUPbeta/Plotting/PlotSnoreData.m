function PlotSnoreData(BreathDataTable,BreathSnoreTable, BreathFLDataTable)
bins = 20;

%% Mean Inspiratory Snore 
SnoreVarTemp = 'WinVarMean_i';
%% All breaths (wake and sleep)
VI = BreathDataTable.VI;
Vdrive = BreathDataTable.VdrivePes;
VdriveNorm = BreathDataTable.VdrivePesNorm;
SnoreVar = BreathSnoreTable.(SnoreVarTemp);
SnoreVar_dB = 10.*(log10(SnoreVar./(0.00002).^2));
PlotSnore(VI,Vdrive,VdriveNorm,SnoreVar_dB,bins, 'Mean Insp Snore Power During Sleep and Wake')

%% During sleep only
SleepIdx = BreathDataTable.hypnog_B < 4; 
VI = BreathDataTable.VI(SleepIdx);
Vdrive = BreathDataTable.VdrivePes(SleepIdx);
VdriveNorm = BreathDataTable.VdrivePesNorm(SleepIdx);
SnoreVar = BreathSnoreTable.(SnoreVarTemp)(SleepIdx);
SnoreVar_dB = 10.*(log10(SnoreVar./(0.00002).^2));
PlotSnore(VI, Vdrive, VdriveNorm,SnoreVar_dB,bins,'Mean Insp Snore Power During Sleep Only')

%% During sleep only and arousal breaths removed
SleepIdx = BreathDataTable.hypnog_B < 4 &  BreathDataTable.ARei == 0;
VI = BreathDataTable.VI(SleepIdx);
Vdrive = BreathDataTable.VdrivePes(SleepIdx);
VdriveNorm = BreathDataTable.VdrivePesNorm(SleepIdx);
SnoreVar = BreathSnoreTable.(SnoreVarTemp)(SleepIdx);
SnoreVar_dB = 10.*(log10(SnoreVar./(0.00002).^2));
PlotSnore(VI, Vdrive, VdriveNorm,SnoreVar_dB,bins,'Mean Insp Snore Power During Sleep, No Arousals')


%% Max Inspiratory Snore 
SnoreVarTemp = 'WinVarMax_i';
%% During sleep only and arousal breaths removed
SleepIdx = BreathDataTable.hypnog_B < 4 &  BreathDataTable.ARei == 0;
VI = BreathDataTable.VI(SleepIdx);
Vdrive = BreathDataTable.VdrivePes(SleepIdx);
VdriveNorm = BreathDataTable.VdrivePesNorm(SleepIdx);
SnoreVar = BreathSnoreTable.(SnoreVarTemp)(SleepIdx);
SnoreVar_dB = 10.*(log10(SnoreVar./(0.00002).^2));
PlotSnore(VI, Vdrive, VdriveNorm,SnoreVar_dB,bins,'Max Insp Snore Power During Sleep, No Arousals')

%% Spectral centroid 
SnoreVarTemp = 'centroidMean_i';
%% During sleep only and arousal breaths removed
SleepIdx = BreathDataTable.hypnog_B < 4 &  BreathDataTable.ARei == 0;
VI = BreathDataTable.VI(SleepIdx);
Vdrive = BreathDataTable.VdrivePes(SleepIdx);
VdriveNorm = BreathDataTable.VdrivePesNorm(SleepIdx);
SnoreVar = BreathSnoreTable.(SnoreVarTemp)(SleepIdx);
% SnoreVar_dB = 10.*(log10(SnoreVar./(0.00002).^2));
PlotSnore(VI, Vdrive, VdriveNorm,SnoreVar,bins,'Spectral Centroid Inspiratory Snore Power During Sleep, No Arousals')

%% Power < 750 Hz
SnoreVarTemp = 'pLess750Mean_i';
%% During sleep only and arousal breaths removed
SleepIdx = BreathDataTable.hypnog_B < 4 &  BreathDataTable.ARei == 0;
VI = BreathDataTable.VI(SleepIdx);
Vdrive = BreathDataTable.VdrivePes(SleepIdx);
VdriveNorm = BreathDataTable.VdrivePesNorm(SleepIdx);
SnoreVar = BreathSnoreTable.(SnoreVarTemp)(SleepIdx);
% SnoreVar_dB = 10.*(log10(SnoreVar./(0.00002).^2));
PlotSnore(VI, Vdrive, VdriveNorm,SnoreVar,bins,'Snore Power at <750Hz Snore Power During Sleep, No Arousals')

%% Power > 750 Hz
SnoreVarTemp = 'pGrtr750Mean_i';
%% During sleep only and arousal breaths removed
SleepIdx = BreathDataTable.hypnog_B < 4 &  BreathDataTable.ARei == 0;
VI = BreathDataTable.VI(SleepIdx);
Vdrive = BreathDataTable.VdrivePes(SleepIdx);
VdriveNorm = BreathDataTable.VdrivePesNorm(SleepIdx);
SnoreVar = BreathSnoreTable.(SnoreVarTemp)(SleepIdx);
% SnoreVar_dB = 10.*(log10(SnoreVar./(0.00002).^2));
PlotSnore(VI, Vdrive, VdriveNorm,SnoreVar,bins,'Snore Power at >750Hz Snore Power During Sleep, No Arousals')


%% Mean Expiratory Snore 
SnoreVarTemp = 'WinVarMean_e';
%% During sleep only and arousal breaths removed
SleepIdx = BreathDataTable.hypnog_B < 4 &  BreathDataTable.ARei == 0;
VI = BreathDataTable.VI(SleepIdx);
Vdrive = BreathDataTable.VdrivePes(SleepIdx);
VdriveNorm = BreathDataTable.VdrivePesNorm(SleepIdx);
SnoreVar = BreathSnoreTable.(SnoreVarTemp)(SleepIdx);
SnoreVar_dB = 10.*(log10(SnoreVar./(0.00002).^2));
PlotSnore(VI, Vdrive, VdriveNorm,SnoreVar_dB,bins,'Mean Exp Snore Power During Sleep, No Arousals')

%% Meax Expiratory Snore 
SnoreVarTemp = 'WinVarMax_e';
%% During sleep only and arousal breaths removed
SleepIdx = BreathDataTable.hypnog_B < 4 &  BreathDataTable.ARei == 0;
VI = BreathDataTable.VI(SleepIdx);
Vdrive = BreathDataTable.VdrivePes(SleepIdx);
VdriveNorm = BreathDataTable.VdrivePesNorm(SleepIdx);
SnoreVar = BreathSnoreTable.(SnoreVarTemp)(SleepIdx);
SnoreVar_dB = 10.*(log10(SnoreVar./(0.00002).^2));
PlotSnore(VI, Vdrive, VdriveNorm,SnoreVar_dB,bins,'Max Exp Snore Power During Sleep, No Arousals')

%% Centroid Expiration
SnoreVarTemp = 'centroidMean_e';
%% During sleep only and arousal breaths removed
SleepIdx = BreathDataTable.hypnog_B < 4 &  BreathDataTable.ARei == 0;
VI = BreathDataTable.VI(SleepIdx);
Vdrive = BreathDataTable.VdrivePes(SleepIdx);
VdriveNorm = BreathDataTable.VdrivePesNorm(SleepIdx);
SnoreVar = BreathSnoreTable.(SnoreVarTemp)(SleepIdx);
SnoreVar_dB = 10.*(log10(SnoreVar./(0.00002).^2));
PlotSnore(VI, Vdrive, VdriveNorm,SnoreVar_dB,bins,'Spectral Centroid Expiratory Snore During Sleep, No Arousals')


%% Correlation with flow shape
SnoreVarTemp = 'WinVarMean_e';
% SleepIdx = BreathDataTable.hypnog_B < 4 &  BreathDataTable.ARei == 0;
SleepIdx = BreathDataTable.Etype == 4;
SnoreVar = BreathSnoreTable.(SnoreVarTemp)(SleepIdx);

SnoreVar_dB = 10.*(log10(SnoreVar./(0.00002).^2));
% SnoreVar_dB = SnoreVar;
FlowShape = BreathFLDataTable.AsymmetryExp_O(SleepIdx);

figure
scatter(FlowShape, SnoreVar_dB)

[R,P] = corr(FlowShape, SnoreVar_dB,'Type','Spearman','Rows','pairwise')

for i = 1:size(BreathFLDataTable,2)
    [R,P] = corr(BreathFLDataTable{SleepIdx,i}, SnoreVar_dB,'Type','Spearman', 'Rows','pairwise');
    correlations(i,1:2) = [R,P];
end
[corrsort,corrsorti] = sort(correlations(:,1));

function PlotSnore(VI,Vdrive,VdriveNorm,SnoreVar,bins,TitleName)
% Discretize for heatmap
% VIEdges = linspace(0,ceil(max(VI)),bins);
VIEdges = prctile(VI,5:5:100);
discTemp = discretize(VI,VIEdges);
discTemp(isnan(discTemp)) = length(VIEdges)+1;
VIEdgesNaN = [VIEdges,NaN];
VIdisc = VIEdgesNaN(discTemp)';

% VdriveEdges = linspace(0,ceil(max(Vdrive)),bins);
VdriveEdges = prctile(Vdrive,5:5:100);
discTemp = discretize(Vdrive,VdriveEdges);
discTemp(isnan(discTemp)) = length(VdriveEdges)+1;
VdriveEdgesNaN = [VdriveEdges,NaN];
VdriveDisc = VdriveEdgesNaN(discTemp)';

% VdriveNormEdges = linspace(0,ceil(max(VdriveNorm)),bins);
VdriveNormEdges = prctile(VdriveNorm,5:5:100);
discTemp = discretize(VdriveNorm,VdriveNormEdges);
discTemp(isnan(discTemp)) = length(VdriveNormEdges)+1;
VdriveNormEdgesNaN = [VdriveNormEdges,NaN];
VdriveNormDisc = VdriveNormEdgesNaN(discTemp)';

% make table
tbl = table(VIdisc,VdriveDisc,VdriveNormDisc,SnoreVar);

figure('Position', [1686 546 1672 420])
heatmap(tbl,'VdriveDisc','VIdisc','ColorVariable','SnoreVar','ColorMethod','Median');
title(TitleName)
ylabel('Ventilation (Fraction of Eupnea)'); xlabel('Drive (Intended Ventilation, L/min)')
% set(gca,'view',[90 -90])

figure('Position', [-3 365 1672 598])
heatmap(tbl,'VdriveNormDisc','VIdisc','ColorVariable','SnoreVar','ColorMethod','Median');
title(TitleName)
ylabel('Ventilation (Fraction of Eupnea)'); xlabel('Drive (Intended Ventilation, Frac. Eupnea)')

figure('Position', [1683 40 560 420])
scatter(VI, SnoreVar)
title(TitleName)
ylabel('Ventilation (Frac. Eupnea)'); xlabel('Snore Power (Pa)')

figure('Position', [2250 41 560 420])
scatter(Vdrive, SnoreVar)
title(TitleName)
xlabel('Drive (Intended Ventilation, L/min)'); ylabel('Snore Power (Pa)')

figure('Position',[2799 41 560 420])
scatter(VdriveNorm, SnoreVar)
title(TitleName)
xlabel('Drive (Intended Ventilation, Frac. Eupnea)'); ylabel('Snore Power (Pa)')

figure('Position', [9 46 560 420])
histogram(VI, VIEdges)

figure('Position', [595 67 560 420])
histogram(VdriveNorm, VdriveNormEdges)
set(gca,'view',[90 -90])
