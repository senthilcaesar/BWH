close
clear
clc

showfigures = 1;

% path to data
LocalDataDir = 'C:\PSG_Data\QAO';
%DataPnasal = [LocalDataDir, '\FS_FD_PtData_wPnasal_OriginalT_noSQRT.mat'];
%DataFlow = [LocalDataDir, '\FS_FD_PtData_wFlow_OriginalT.mat'];

%DataPnasal = [LocalDataDir, '\FS_FD_g_wPnasal_OriginalT_noSQRT.mat'];
DataPnasal = [LocalDataDir, '\FS_FD_g_wPnasal_OriginalT_wSQRT.mat'];
DataFlow = [LocalDataDir, '\FS_FD_g_wFlow_OriginalT.mat'];

load(DataPnasal);
PtDataPnasal = PtData;
clear PtData

load(DataFlow);
PtDataFlow = PtData;
clear PtData

PnasalDat = [];
FlowDat = [];
BBstats = [];

% work on one PT at a time
for pt=25 %[3, 5, 10, 17, 22, 23, 25, 26]
    str=['Processing patient: ', num2str(pt)]; disp(str);
    
    PtFlow = PtDataFlow(PtDataFlow.PT==pt,:);
    PtPnasal = PtDataPnasal(PtDataPnasal.PT==pt,:);
    
    % just using time
    PtFlowTime = PtFlow.BB_time;
    PtPnasalTime = PtPnasal.BB_time;

    %% remove offset (usually present)
    F_1 = PtFlowTime(1); 
    P_1 = PtPnasalTime(1); 
    if P_1 < F_1
        loc = find(PtPnasalTime<PtFlowTime(1), 1, 'last');
        PtPnasalTime(1:loc) = [];
    else
        loc = find(PtFlowTime<PtPnasalTime(1), 1, 'last');
        PtFlowTime(1:loc) = [];
    end  
    
    if showfigures
        figure(1); clf(figure(1)); subplot(1,2,1);
        stairs(PtPnasalTime); hold on;
        stairs(PtFlowTime); axis square;
    end
    
    %% attempt one - a perfect match
    if 0
        % find perfect matches (using pnasal indices)
        matchedBB = intersect(PtPnasalTime, PtFlowTime);
        matchedBB2 = ismember(PtPnasalTime, PtFlowTime);
        wherearethey_inPnasal = find(matchedBB2);
        % find where these are in the flow indices
        matchedBB3 = ismember(PtFlowTime, PtPnasalTime);
        wherearethey_inFlow = find(matchedBB3);
        % perfect matches
        matches = [wherearethey_inPnasal, wherearethey_inFlow];
        % what is the gap
        deltas = diff(matches, 1,2);
    end
    
    %% attempt two - an approx match
    if 0
        % round times to effectively allow a margin or tolerance
        PtPnasalTimeRnd = round(PtPnasalTime, 1);
        PtFlowTimeRnd = round(PtFlowTime, 1);
        % find approximate matches (using pnasal indices)
        matchedBB = intersect(PtPnasalTimeRnd, PtFlowTimeRnd);
        matchedBB2 = ismember(PtPnasalTimeRnd, PtFlowTimeRnd);
        wherearethey_inPnasal = find(matchedBB2);
        % find where these are in the flow indices
        matchedBB3 = ismember(PtFlowTimeRnd, PtPnasalTimeRnd);
        wherearethey_inFlow = find(matchedBB3);
        % approximate matches
        matches = [wherearethey_inPnasal, wherearethey_inFlow];
        % what is the gap
        deltas = diff(matches, 1,2);
    end
    
    %% attempt three - a tolerable match
    % look for breaths within a given tolerance range
    tolerance = 0.5;
    BB_tolerance = @(T) [T-tolerance,T+tolerance];
    
    % step through PnasalTime, looking for tolerable matches in FlowTime
    matchedBB_PtoF = NaN(length(PtPnasalTime),1);
    for bb=1:length(PtPnasalTime)
        [range]=BB_tolerance(PtPnasalTime(bb)); % get range
        [candidates] = find(PtFlowTime>range(1) & PtFlowTime<range(2));
        switch size(candidates,1)
            case 0
                %str=['No match for breath ', num2str(bb)]; disp(str);
                matchedBB_PtoF(bb) = NaN;
            case 1
                matchedBB_PtoF(bb) = candidates;
            otherwise % more than one
                str=['Multiple matches for breath ', num2str(bb)]; disp(str);
        end
    end
    
    if 0 % figure
        
    xlim([0 3E4]);
    xlim([range(1)-30 range(2)+30]);
    
    ax5(1) = subplot(2,1,1);hold on;
    plot(PtFlowTime, 0, 'gx');
    
    plot(range(1), 0, 'r^');
    plot(range(2), 0, 'r^');
    
    ax5(2) = subplot(2,1,2); hold on;
    plot(PtPnasalTime, 0, 'go');
    
    plot(range(1), 0, 'b^');
    plot(range(2), 0, 'b^');
    
    end
    %nnz(isnan(matchedBB_PtoF))
    %abs(length(PtPnasalTime)-length(PtFlowTime))
    
    %% set up indices into both data
    matches=[];
    % step through each row in matchedBB_PtoF, and get corresponding FlowData
    for bb=1:length(matchedBB_PtoF)
        row = matchedBB_PtoF(bb);
        if isnan(row)
            continue
        else
            matches=[matches;[row,bb]]; %  Flow, Pnasal
        end
    end
    
    PtPnasal_sorted = PtPnasal(matches(:,2),:);
    PtFlow_sorted = PtFlow(matches(:,1),:);
    
    if showfigures
        figure(1); subplot(1,2,2);
        stairs(PtPnasal_sorted.BB_time); hold on;
        stairs(PtFlow_sorted.BB_time); axis square;
    end
    
    %% recreate one large table for each
    PnasalDat = [PnasalDat; PtPnasal_sorted];
    FlowDat = [FlowDat; PtFlow_sorted];
    
    
    %% table of breath stats
    BBstats = [BBstats; [pt, length(PtPnasalTime), length(PtFlowTime), size(matches,1)]];
end

BBstats

%% random testing
% it would appear that on occasion, we may be mapping to an adjacent breath
% (this should be checked visually, looking at flow and BB_time marks)
figure(2); clf(figure(2));
scatter(PnasalDat.MIF50, FlowDat.MIF50);
refline(1,0);
xlabel('Pnasal based'); ylabel('Flow based'); title('MIF50');

figure(3); clf(figure(3));
scatter(PnasalDat.g_Edi_Adj, FlowDat.g_Edi_Adj);
refline(1,0);
xlabel('Pnasal based'); ylabel('Flow based'); title('gEdiAdj');

figure(4); clf(figure(4));
scatter(PnasalDat.BB_Ttot, FlowDat.BB_Ttot);
refline(1,0);
xlabel('Pnasal based'); ylabel('Flow based'); title('Ttot');

%% plot both Flow and Pnasal, and both times, for one pt
pt = 25;

% keep the breath times, and breath data table for this pt
tP = PnasalDat{PnasalDat.PT==pt,:};
tF = FlowDat{FlowDat.PT==pt,:};

% load the DataHypnogMat for this pt
LocalDirectory = 'C:\Users\uqdmann\Documents\MATLAB\QAO\FlowDrive';
DataSpreadsheet = [LocalDirectory,'\AnalyzeDataSpreadsheet.xlsx'];
[~,~,raw] = xlsread(DataSpreadsheet,1,'B3:C32');
matfile = [char(raw(pt,2)),'\',char(raw(pt,1)),'.mat'];
load(matfile);

time = DataEventHypnog_Mat(:,1);
Flow = DataEventHypnog_Mat(:,4);
Pnasal_raw = DataEventHypnog_Mat(:,22);
Pnasal = Pnasal_raw-nanmean(Pnasal_raw);

%% determine the good regions (as the not NaN regions)
good_data = ~(any(isnan([Flow Pnasal]),2));

figure(6); clf(figure(6));
plot(time, Flow, 'b'); hold on;
plot(time, Pnasal, 'r');
plot(time, good_data, 'g');

%%
bb = size(tP,1);
%for bb=1:size(tP,1)
    % get the row (breath) data
    tP_ = tP(1:bb,:);
    tF_ = tF(1:bb,:);
    
    % get the start and end times for both
    P_start = tP_(:,16);
%    P_end = tP_(17)+P_start;
    F_start = tF_(:,16);
%    F_end = tF_(17)+F_start;
    
    % set upper and lower limits
%    lwr = min([P_start, P_end, P_start, P_end])-10;
%    upr = max([P_start, P_end, P_start, P_end])+10;
    
    % show it
    figure(5); clf(figure(5));
    ax5(1) = subplot(2,1,1);
    plot(time, Flow, 'b'); hold on;
    plot(time, good_data);
    plot(F_start, 0, 'r^');
%    plot(F_end, 0, 'rv');
    ylabel('Flow');
    ax5(2) = subplot(2,1,2); hold on;
    plot(time, Pnasal, 'r');
    plot(P_start, 0, 'b^');
%    plot(P_end, 0, 'bv');
    ylabel('Pnasal');
    linkaxes(ax5, 'x');
%    xlim([lwr, upr]);
     xlabel('Time');
%end