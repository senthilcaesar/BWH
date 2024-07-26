% Settings
hypop = 0; % use hypopnea data
flowlim = 1; removeApnea = 1; % use flow limitation data / remove apnea
%Load data
load('BreathTableBig.mat')
if hypop == 1 && ~exist('BreathDataTableFinal')
    load('HypopTables.mat') %hypopnea breath data
elseif flowlim == 1 && ~exist('BreathDataTableFinal')
    load('FLTables.mat') %flow limited breath data
end

% If flow limited, remove apneas
if removeApnea
    apneaIdx = BreathDataTableFinal.Etype == 2;
    BreathFLDataTableFinal(apneaIdx,:) = [];
    InspArrayFinal(apneaIdx,:) = [];
    BreathDataTableFinal(apneaIdx,:) = [];
end

%load subject info
if ~exist('subjectInfo')
    subjectInfo = readtable('AnalyzeDataSpreadsheet.xlsx','Sheet', 1,'Range','B2:CT34');
    subjectInfo.labels = subjectInfo.PercentReduction >= 70;
end

filepath = 'J:\PEOPLE\FACULTY\SANDS\OralApplianceMM2018\Converted';
fnames = struct2cell(dir(filepath));
Fs = 125;
for subnum = 23:size(subjectInfo,1)
    close all
    if subjectInfo.TOTAL_SUP_AIH(subnum) < 20
        continue
    end
    
    subId = str2num(subjectInfo.MATFilename{subnum}(1:end-8))
    subIdx = BreathDataTableFinal.Subject == subId;
    
    
    %get filename
    filenameIdx = contains(fnames(1,:), num2str(subId));
    filename = fnames{1,filenameIdx};
    
    %load file
    load([filepath, '\', filename])
    
    subIdx = BreathDataTableFinal.Subject == subId;
    
    BreathDataTableFinalSub = BreathDataTableFinal(subIdx,:);
    BreathFLDataTableFinalSub = BreathFLDataTableFinal(subIdx,:);
    
    FLFlow = nan(size(DataEventHypnog_Mat,1),1);
    Feature1 = nan(size(DataEventHypnog_Mat,1),1);
    Feature2 = nan(size(DataEventHypnog_Mat,1),1);
    Feature3 = nan(size(DataEventHypnog_Mat,1),1);
    for ii = 1:size(BreathDataTableFinalSub,1)
        Sample0 = BreathDataTableFinalSub.Time0(ii)*Fs;       
        BBIdx = round(Sample0+BreathDataTableFinalSub.BB_i_start(ii)-1):...
            round(Sample0+BreathDataTableFinalSub.BB_i_end(ii));
        FLFlow(BBIdx) = DataEventHypnog_Mat(BBIdx,2);
        Feature1(BBIdx) = BreathFLDataTableFinalSub.SS_Area_O(ii);
        Feature2(BBIdx) = BreathFLDataTableFinalSub.DV_NED1_O(ii);
        Feature3(BBIdx) = BreathFLDataTableFinalSub.DV_NED2_O(ii);
    end
    

    if flowlim == 1
        subIdx2 = BreathDataTableBig.Subject == subId;
        BreathDataTableBigSub = BreathDataTableBig(subIdx2,:);
        for jj = 1:size(BreathDataTableBigSub,1)
            Sample0 = BreathDataTableBigSub.Time0(jj)*Fs;
            if isnan(Sample0)
                continue
            end
            BBIdx = round(Sample0+BreathDataTableBigSub.BB_i_start(jj)-1):...
                 round(Sample0+BreathDataTableBigSub.BB_i_end(jj));
            FlowDrive(BBIdx) = BreathDataTableBigSub.FlowDrive(jj);
        end
    end
    
    if subjectInfo.labels(subnum) == 1
        response = 'R';
    else
        response = 'NR';
    end
    
    cutline = ones(length(DataEventHypnog_Mat(:,2)),1)*0.7;
    % Plot
%     figure('pos', [-1665 102 1651 851])
    figure('pos', [-1665 507 1651 446])
%     ax1 = subplot(211);
    r = stairs(FlowDrive,'Color', [0.6 0.6 0.6], 'LineWidth', 2);
    hold on
    plot(cutline, '--k')
    plot(DataEventHypnog_Mat(:,2),'b')
    plot(FLFlow, 'r')
    stairs(Feature1,'g', 'LineWidth', 2);
    stairs(Feature2,'m', 'LineWidth', 2);
    stairs(Feature3,'k', 'LineWidth', 2);
    title([num2str(subId), ' ', response])
    ylabel('Flow')
    ylim([-1.5 1.5])
    set(gca,'Box','On','FontSize',12)
    
%     ax2 = subplot(212);
%     stairs(FlowDrive,'b', 'LineWidth', 2)
%     ylabel('Flow')
%     ylim([0 1.25])
%     
%     linkaxes([ax1,ax2], 'x')
end