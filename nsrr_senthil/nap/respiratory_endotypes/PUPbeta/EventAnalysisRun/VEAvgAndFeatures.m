% Get published versions of VE features
function [VEFeaturesAllPos, VEFeaturesSup, VEFeaturesNSup, VEtable] = VEAvgAndFeatures(BreathDataTable,DataEventHypnog_Mat,ChannelsList,Fs)

global settings n

% Identify apnea and hypopnea number (to easily get terminal event)
BreathDataTable_original = BreathDataTable;
BreathDataTable = IsolateApneaAndHypopnea(BreathDataTable_original);

%% Generate time series ventilation signal
Time = DataEventHypnog_Mat(:,1); 
Flow = DataEventHypnog_Mat(:,2); % Isolate flow from matrix

% Initialize matrices to store resamples flow and mean flow volume
VE_td = nan(size(Flow,1),1);

% Loop through all breaths, sample VE, then store in VE array that
% aligns with time array. VE will be read to resample
for brNum = 1:size(BreathDataTable,1)
    % Store Ventilation (relative to VEupnea) in time domain
    startIdx = round(BreathDataTable.Time_start(brNum)*Fs - Time(1)*Fs);
    endIdx = round(BreathDataTable.Time_end(brNum)*Fs - Time(1)*Fs);
    VE_td(startIdx:endIdx,1) =  BreathDataTable.VI(brNum); 
end

% Resample VE_td
Fs_RS = 4;
newTime = Time(1):(1/4):Time(end);
VE_tdRS = interp1(Time, VE_td, newTime, 'previous');

% Sample data for next set of operations
termEvtIdx = find(BreathDataTable.HypopNum == 1 | ...
    BreathDataTable.ApneaNum == 1)'; % Find terminal events
brCount = 0;
dTime = 150;

% Get position info an convert
PosRaw = DataEventHypnog_Mat(:, contains(ChannelsList, 'Pos'));

if ~isempty(PosRaw)
%     [~,~,poscodesdatabase] = xlsread(MasterSpreadsheet,3,'B2:J100');
%     [~,~,protocol] = xlsread(MasterSpreadsheet,1,'AF4:AF999');
    positioncodesout = PositionCodeLookup(settings.poscodesdatabase,'Output');
    positioncodes = PositionCodeLookup(settings.poscodesdatabase,settings.protocol{n});    
    Pos = PositionTranslator(positioncodes,positioncodesout,PosRaw);
else
    Pos = ones(size(DataEventHypnog_Mat,1),1);
    disp('No position data')
end

% Initialize arrays
VE = nan(length(termEvtIdx), dTime*Fs_RS*2+1);
SleepStage = nan(length(termEvtIdx),1);
PosArray = nan(length(termEvtIdx),1);
TtotArray = nan(length(termEvtIdx),1);


%% Ensemble Average VE During Respiratory Events
for evtBreath = termEvtIdx
    brCount = brCount+1;

    % Get event breath indices and check size
    evtBrthIdx = round(BreathDataTable.Time_end(evtBreath)*Fs_RS) - round(newTime(1)*Fs_RS);

    % Conditions for setting pre and post 
    if BreathDataTable.Time_end(evtBreath) - BreathDataTable.Time_start(1) < dTime %start condition
        idxPreEvt = round(BreathDataTable.Time_end(evtBreath)*Fs_RS) - ...
            (round(BreathDataTable.Time_end(evtBreath)*Fs_RS) - ...
            round(BreathDataTable.Time_start(1)*Fs_RS)) -...
            round(newTime(1)*Fs_RS);
        idxPostEvt = round(BreathDataTable.Time_end(evtBreath)*Fs_RS) + ...
            round(dTime*Fs_RS) - round(newTime(1)*Fs_RS);
    elseif BreathDataTable.Time_end(end) - BreathDataTable.Time_start(evtBreath) < dTime %end condition
        idxPreEvt = round(BreathDataTable.Time_end(evtBreath)*Fs_RS) - ...
            round(dTime*Fs_RS) - round(newTime(1)*Fs_RS);
        idxPostEvt = round(BreathDataTable.Time_end(evtBreath)*Fs_RS) + ...
            (round(BreathDataTable.Time_end(end)*Fs_RS) - ...
            round(BreathDataTable.Time_start(evtBreath)*Fs_RS)) - ...
            round(newTime(1)*Fs_RS);
    else
        idxPreEvt = round(BreathDataTable.Time_end(evtBreath)*Fs_RS) - ...
            round(dTime*Fs_RS) - round(newTime(1)*Fs_RS);
        idxPostEvt = round(BreathDataTable.Time_end(evtBreath)*Fs_RS) + ...
            round(dTime*Fs_RS) - round(newTime(1)*Fs_RS); % Time_end of terminal event at 601 (dTime*Fs_RS + 1)
    end

    if 0
        figure('Position', [54 545 560 420]), plot(VE_tdRS(idxPreEvt:idxPostEvt)), hold on
        plot(601, VE_tdRS(idxPreEvt+600), 'o'), hold off
    end

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

    % Store VE
    VEtemp = nan(1, length(1:dTime*Fs_RS*2+1));
    VEtemp(idxPreEvt2:idxPostEvt2) = VE_tdRS(newIdxPreEvt:newIdxPostEvt);
    VE(brCount,1:dTime*Fs_RS*2+1) = VEtemp;

    if 0
        figure('Position', [630 545 560 420]), plot(VEtemp(idxPreEvt2:idxPostEvt2)), hold on
        plot(601, VEtemp(601), 'o'), hold off
    end

    % Store Ttot (estimate based on 40 breaths through event) not ideal
    % but easy to follow and likely reasonable estimate
    if evtBreath-20 <= 0
        brPre = 0; brPost = 20;
    elseif evtBreath+20 >= size(BreathDataTable,1)
        brPost = size(BreathDataTable,1)-evtBreath;
    else
        brPre = 20; brPost = 20;
    end

    TtotArray(brCount,1) = nanmean(BreathDataTable.Time_end(evtBreath-brPre:evtBreath+brPost) - ...
        BreathDataTable.Time_start(evtBreath-brPre:evtBreath+brPost));

    % Find sleep position during event
    PosArray(brCount, 1) = Pos(round(BreathDataTable.Time_end(evtBreath)*Fs)...
        - round(Time(1)*Fs));

    % Find sleep stage during event
    SleepStage(brCount, 1) = BreathDataTable.hypnog_B(evtBreath);
end

VEtable = table(VE,SleepStage,PosArray);

%% Calculate VE features
VEFeaturesAllPos = ComputeVEFeatures(VE,SleepStage,PosArray,TtotArray,'All');
VEFeaturesSup = ComputeVEFeatures(VE,SleepStage,PosArray,TtotArray,'Supine');
VEFeaturesNSup = ComputeVEFeatures(VE,SleepStage,PosArray,TtotArray,'Nonsupine');

%%
% SubTag(regexp(SubTag,'[-]'))=[]; % remove dashes (invalid characters)
% VEStruct.(['a',SubTag]) = VE_AllSub;
% SleepType.(['a',SubTag]).SleepStage = SleepStage;
% SleepType.(['a',SubTag]).PosArray = PosArray;
% SleepType.(['a',SubTag]).Ttot = TtotArray;
end

function BreathDataTableNew = IsolateApneaAndHypopnea(BreathDataTable)
    BreathDataTableNew = BreathDataTable;
    
    %% Isolate Hypopnea Breaths
    %Generate large arrays with only hypopnea breath
    %Add new column indicating hypopnea number
    %Can later filter by hypopnea number
    Etype4 = BreathDataTable.Etype;
    Etype4(Etype4 ~= 4 & Etype4 ~= 6) = 0;
%     Etype4(Etype4 ~= 4) = 0;
    
    HypopShiftDown = circshift(Etype4,1);
    HypopShiftUp = circshift(Etype4,-1);
    HypopDiffDown = Etype4 - HypopShiftDown;
    HypopDiffUp = Etype4 - HypopShiftUp;

    HypopFirstIdx = find(HypopDiffDown == 4 | HypopDiffDown == 6);
    HypopLastIdx = find(HypopDiffUp == 4 | HypopDiffDown == 6);
    
%     HypopFirstIdx = find(HypopDiffDown == 4);
%     HypopLastIdx = find(HypopDiffUp == 4);

    %Generate variable indicating hypopnea number to be stored in table
    HypopNum = zeros(size(Etype4,1),1);
    NumHypop = zeros(size(Etype4,1),1);
    %Loop through hypopnea indices and fill hypopnea number variable
    %Note: 1 is the terminal hypopnea, 2 is second last, etc
    for ii = 1:length(HypopFirstIdx)
        totalHypop = HypopLastIdx(ii) - HypopFirstIdx(ii) + 1;
        HypopNum(HypopFirstIdx(ii):HypopLastIdx(ii),1) = totalHypop:-1:1;
        NumHypop(HypopFirstIdx(ii):HypopLastIdx(ii),1) = totalHypop;
    end

    % Add to table
    BreathDataTableNew.HypopNum = HypopNum;
    BreathDataTableNew.NumHypop = NumHypop;

    %% Isolate Apnea breaths
    Etype2 = BreathDataTable.Etype;
    Etype2(Etype2 ~= 2 & Etype2 ~=3 & Etype2 ~= 5) = 0;
%     Etype2(Etype2 ~= 2 & Etype2 ~= 5) = 0; %exclude central events
%     Etype2(Etype2 ~= 2) = 0; %exclude central events
    
    ApneaShiftDown = circshift(Etype2,1);
    ApneaShiftUp = circshift(Etype2,-1);
    ApneaDiffDown = Etype2 - ApneaShiftDown;
    ApneaDiffUp = Etype2 - ApneaShiftUp;
% 
    ApneaFirstIdx = find(ApneaDiffDown == 2 | ApneaDiffDown == 3 | ApneaDiffDown == 5);
    ApneaLastIdx = find(ApneaDiffUp == 2 | ApneaDiffUp == 3 | ApneaDiffUp == 5);
    
%     ApneaFirstIdx = find(ApneaDiffDown == 2 | ApneaDiffDown == 5); % exclude central
%     ApneaLastIdx = find(ApneaDiffUp == 2 | ApneaDiffUp == 5);

    %Generate variable indicating apnea number to be stored in table
    ApneaNum = zeros(size(Etype2,1),1);
    NumApnea = zeros(size(Etype2,1),1);
    %Loop through apnea indices and fill apnea number variable
    %Note: 1 is the terminal hypopnea, 2 is second last, etc
    for ii = 1:length(ApneaFirstIdx)
        totalApnea = ApneaLastIdx(ii) - ApneaFirstIdx(ii) + 1;
        ApneaNum(ApneaFirstIdx(ii):ApneaLastIdx(ii),1) = totalApnea:-1:1;
        NumApnea(ApneaFirstIdx(ii):ApneaLastIdx(ii),1) = totalApnea;
    end

    % Add to table
    BreathDataTableNew.ApneaNum = ApneaNum;
    BreathDataTableNew.NumApnea = NumApnea;
end

function VEFeatures = ComputeVEFeatures(VE,SleepStage,PosArray,Ttot,PosIn)
    %% Compute Features from VE data
    hypopIdx = 601;
    Fs_RS = 4; % Sampling frequency
    dtRS = 1/Fs_RS;

    % options
    Stage = 'NREM';
    Position = PosIn;
    
    if strcmp(Stage, 'NREM') && strcmp(Position, 'Supine')
        VErows = SleepStage ~=3 & PosArray == 1;
    elseif strcmp(Stage, 'NREM') && strcmp(Position, 'Nonsupine')
        VErows = SleepStage ~=3 & PosArray ~= 1;
    elseif strcmp(Stage, 'NREM') && strcmp(Position, 'All')
        VErows = SleepStage ~=3;
    end
        
    VEtoAnalyze = VE(VErows,:);
    Ttot_Avg = nanmean(Ttot(VErows));
    meanVE = nanmean(VEtoAnalyze, 1);
    medianVE = median(VEtoAnalyze, 1, 'omitnan');
    
    % Number of breaths included in mean
    numBreaths = sum(~isnan(VEtoAnalyze(:,hypopIdx)));
%     numBreaths
    
    %% Initialize variables
    minVE_mean=nan; minVE_2_mean=nan; minVE_3_mean=nan;
    minVE_4_mean=nan; minVE_5_mean=nan; meanVEDist_mean=nan; minVE_med=nan;  minVE_2_med=nan; 
    minVE_3_med=nan; minVE_4_med=nan; minVE_5_med=nan; Dist=nan; VRAavg=nan; LG=nan;
    dt1_mean=nan; dt1_med=nan; dt2_mean=nan; dt2_med=nan; dt5_mean=nan; dt5_med=nan;
    meanVE95to1_mean=nan; meanVE95to1_med=nan; meanVE95to2_mean=nan; meanVE95to2_med=nan;
    meanVE95to5_mean=nan; meanVE95to5_med=nan; slp95to1_mean=nan; slp95to1_med=nan;
    slp95to2_mean=nan; slp95to2_med=nan; slp95to5_mean=nan; slp95to5_med=nan;
    risePostMin_mean=nan; risePostMin_med=nan; dtPostMin_mean=nan; dtPostMin_med=nan;
    slopePostMin_mean=nan; slopePostMin_med=nan;

    if sum(VErows,1) > 5 % only run if sufficient data - otherwise nan ^
    
        %% Mean and minimum ventilation through event 
        % Minimum VE @ or before hypopnea
        [minVE_mean, minVEidx_mean] = min(meanVE(hypopIdx-40:hypopIdx+40));
        [minVE_med, minVEidx_med] = min(medianVE(hypopIdx-40:hypopIdx+40));

        % Lowest VE of 2, 3, 4, and 5 @ or before hypopnea
        minVE_2_mean = mean(mink2(meanVE(hypopIdx-40:hypopIdx+40), 2));
        minVE_2_med = mean(mink2(medianVE(hypopIdx-40:hypopIdx+40), 2));
        minVE_3_mean = mean(mink2(meanVE(hypopIdx-40:hypopIdx+40), 3));
        minVE_3_med = mean(mink2(medianVE(hypopIdx-40:hypopIdx+40), 3));
        minVE_4_mean = mean(mink2(meanVE(hypopIdx-40:hypopIdx+40), 4));
        minVE_4_med = mean(mink2(medianVE(hypopIdx-40:hypopIdx+40), 4));
        minVE_5_mean = mean(mink2(meanVE(hypopIdx-40:hypopIdx+40), 5));
        minVE_5_med = mean(mink2(medianVE(hypopIdx-40:hypopIdx+40), 5));

        % Find when VE begins to dip below eupnea
        EupBreathTemp = find(meanVE(hypopIdx:-1:1) >= 0.95, 1, 'first');
        EupBreathL_mean = hypopIdx - EupBreathTemp+1;

    %     EupBreathTemp = find(medianVE(hypopIdx:-1:1) >= 95, 1, 'first');
    %     EupBreathL_med = hypopIdx - EupBreathTemp+1;

        EupBreathTemp = find(meanVE(hypopIdx:1:end) >= 0.95, 1, 'first');
        EupBreathR_mean = hypopIdx + EupBreathTemp-1;

    %     EupBreathTemp = find(medianVE(hypopIdx:1:end) >= 95, 1, 'first');
    %     EupBreathR_med = hypopIdx + EupBreathTemp-1;

        % Mean VE of all breaths surrounding hypopnea below eupnea
        meanVEDist_mean = 1 - nanmean(meanVE(EupBreathL_mean+1:EupBreathR_mean-1));
    %     meanVEDist_med = 100 - nanmean(medianVE(EupBreathL_mean+1:EupBreathR_med-1));

        % Disturbance
        newTime = 0:dtRS:length(EupBreathL_mean:EupBreathR_mean)*dtRS-dtRS;

        % Compute discrete integral up to the second last point.
        % This is because the transition terminal event to arousal is so rapid. 
        % Therefore, need create an intermediate point that occurs right at 95% VEupnea
        % and add that bit to the disturbance 

        Dist_temp = trapz(newTime(1:end-1), meanVE(EupBreathL_mean:EupBreathR_mean-1));

         % Find time where 95% VEupnea is 
        m = (meanVE(EupBreathR_mean) - meanVE(EupBreathR_mean-1))/...
            (newTime(end) - newTime(end-1));
        t_mid = ((0.95 - meanVE(EupBreathR_mean-1))/m) +  newTime(end-1);
        endBit = 0.95*(t_mid - newTime(end-1));

        %% Dynamic features of reduction in VE during event
        % Find index of the Xth prctile
        Prc1Idx_mean = find(meanVE(hypopIdx-40:hypopIdx+40) <=...
            prctile(meanVE(hypopIdx-40:hypopIdx+40), 1), 1, 'first')...
            + hypopIdx-40+1;
        Prc1Idx_med = find(medianVE(hypopIdx-40:hypopIdx+40) <=...
            prctile(medianVE(hypopIdx-40:hypopIdx+40), 1), 1, 'first')...
            + hypopIdx-40+1;
        Prc2Idx_mean = find(meanVE(hypopIdx-40:hypopIdx+40) <=...
            prctile(meanVE(hypopIdx-40:hypopIdx+40), 2.5), 1, 'first')...
            + hypopIdx-40+1;
        Prc2Idx_med = find(medianVE(hypopIdx-40:hypopIdx+40) <=...
            prctile(medianVE(hypopIdx-40:hypopIdx+40), 2.5), 1, 'first')...
            + hypopIdx-40+1;
        Prc5Idx_mean = find(meanVE(hypopIdx-40:hypopIdx+40) <=...
            prctile(meanVE(hypopIdx-40:hypopIdx+40), 5), 1, 'first')...
            + hypopIdx-40+1;
        Prc5Idx_med = find(medianVE(hypopIdx-40:hypopIdx+40) <=...
            prctile(medianVE(hypopIdx-40:hypopIdx+40), 5), 1, 'first')...
            + hypopIdx-40+1;

        % Time between 95% Veup and Xth prctile (normalized to Ttot so effectively units of breaths)
        dt1_mean = (Prc1Idx_mean - EupBreathL_mean + 1)*dtRS/Ttot_Avg;
        dt1_med = (Prc1Idx_med - EupBreathL_mean + 1)*dtRS/Ttot_Avg;
        dt2_mean = (Prc2Idx_mean - EupBreathL_mean + 1)*dtRS/Ttot_Avg;
        dt2_med = (Prc2Idx_med - EupBreathL_mean + 1)*dtRS/Ttot_Avg;
        dt5_mean = (Prc5Idx_mean - EupBreathL_mean + 1)*dtRS/Ttot_Avg;
        dt5_med = (Prc5Idx_med - EupBreathL_mean + 1)*dtRS/Ttot_Avg;

        % Mean amplitude of data between 95% Veup and Xth prctile
        meanVE95to1_mean = nanmean(meanVE(EupBreathL_mean:Prc1Idx_mean));
        meanVE95to1_med = nanmean(medianVE(EupBreathL_mean:Prc1Idx_med));
        meanVE95to2_mean = nanmean(meanVE(EupBreathL_mean:Prc2Idx_mean));
        meanVE95to2_med = nanmean(medianVE(EupBreathL_mean:Prc2Idx_med));
        meanVE95to5_mean = nanmean(meanVE(EupBreathL_mean:Prc5Idx_mean));
        meanVE95to5_med = nanmean(medianVE(EupBreathL_mean:Prc5Idx_med));

        % Slope of decline in VE during event
        slp95to1_mean = (meanVE(EupBreathL_mean) - meanVE(Prc1Idx_mean))...
            /dt1_mean;
        slp95to1_med = (medianVE(EupBreathL_mean) - meanVE(Prc1Idx_med))...
            /dt1_med;
        slp95to2_mean = (meanVE(EupBreathL_mean) - meanVE(Prc2Idx_mean))...
            /dt2_mean;
        slp95to2_med = (medianVE(EupBreathL_mean) - meanVE(Prc2Idx_med))...
            /dt2_med;
        slp95to5_mean = (meanVE(EupBreathL_mean) - meanVE(Prc5Idx_mean))...
            /dt5_mean;
        slp95to5_med = (medianVE(EupBreathL_mean) - meanVE(Prc5Idx_med))...
            /dt5_med;

        %% Increase in VE between minimum and arousal
        % Find point where VE increases rapidly
        tempVE = meanVE(EupBreathL_mean-1:EupBreathR_mean+1);
        VEshift = circshift(tempVE,1);

        diffVE = tempVE - VEshift;
        meanDiffVE = nanmean(abs(diffVE(1:end-3)));
        stdDiffVE = nanstd(abs(diffVE(1:end-3)));

        IdxAr = EupBreathL_mean-1 + find(diffVE > meanDiffVE+5*stdDiffVE, 1, 'first')-1;

        if ~isempty(IdxAr)
            risePostMin_mean = meanVE(IdxAr-1) - minVE_mean;
            risePostMin_med = medianVE(IdxAr-1) - minVE_med;

            dtPostMin_mean = ((IdxAr-1) - minVEidx_mean)*dtRS/Ttot_Avg;
            dtPostMin_med = ((IdxAr-1) - minVEidx_med)*dtRS/Ttot_Avg;

            slopePostMin_mean = risePostMin_mean/dtPostMin_mean;
            slopePostMin_med = risePostMin_med/dtPostMin_med;
        end

        %% Loop gain and arousal threshold
         % Compute disturbance
        Dist = Dist_temp + endBit;

        % Ventilatory response to arousal 
        [VRAavg, VRAidx] = max(meanVE(hypopIdx:hypopIdx+55));

        % Surrogate loop gain
        LG = VRAavg/Dist;

        % Plot
        if 0
           figure(5)
           plot(meanVE); hold on
           plot(EupBreathL_mean, meanVE(EupBreathL_mean+1), 'ro')
           plot(EupBreathL_mean+t_mid*Fs_RS, 0.95, 'rx')
           plot(minVEidx_mean+hypopIdx-40-1, minVE_mean, 'gx')
           plot(VRAidx+hypopIdx-1, VRAavg, 'm*')
           plot([Prc1Idx_mean Prc2Idx_mean Prc5Idx_mean], ...
               meanVE([Prc1Idx_mean Prc2Idx_mean Prc5Idx_mean]), 'm^')
           plot(IdxAr, meanVE(IdxAr), 'go')
           hold off
           close all
        end

    end
    
    % Create Table
    VEFeatureTable = table(minVE_mean, minVE_2_mean, minVE_3_mean,...
        minVE_4_mean, minVE_5_mean, meanVEDist_mean,minVE_med, minVE_2_med,...
        minVE_3_med, minVE_4_med, minVE_5_med, Dist, VRAavg, LG,...
        dt1_mean, dt1_med, dt2_mean, dt2_med, dt5_mean, dt5_med,...
        meanVE95to1_mean, meanVE95to1_med, meanVE95to2_mean, meanVE95to2_med,...
        meanVE95to5_mean, meanVE95to5_med, slp95to1_mean, slp95to1_med,...
        slp95to2_mean, slp95to2_med, slp95to5_mean, slp95to5_med,...
        risePostMin_mean, risePostMin_med, dtPostMin_mean, dtPostMin_med,...
        slopePostMin_mean, slopePostMin_med);
    
    VEFeatures = table2struct(VEFeatureTable);

end

%% Function to get min K values
function y = mink2(A,k)
A = sort(A);
y = A(1:k);
end