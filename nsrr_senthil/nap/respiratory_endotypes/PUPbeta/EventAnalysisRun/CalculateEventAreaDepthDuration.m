% added by EAS on 2022-02-18
function output = CalculateEventAreaDepthDuration(BreathDataTableFulls, Evts, patientNumber)

    allSettings.startStopCriteria.MANUALLY_SCORED = 1;
    allSettings.startStopCriteria.UNDER_100_VI = 2;
    allSettings.startStopCriteria.UNDER_70_VI = 3;
    EXCLUDE_ZERO_VI_BREATHS = false; % if true, excludes all 0 and NaN ventilation breaths. If false, includes them and and excludes Nan ventilation breaths
    EXCLUDE_LONG_OUTLIER_BREATHS = true; % if true, excludes all breaths with durations over Q3 + 1.5xIQR
    stateCodes.SLEEP_STATE_WAKE = 4;
    stateCodes.SLEEP_STATE_N1 = 2;
    stateCodes.SLEEP_STATE_N2 = 1;
    stateCodes.SLEEP_STATE_N3 = 0;
    stateCodes.SLEEP_STATE_REM = 3;
    stateCodes.SLEEP_STATE_UNKNOWN = 8;
    eventCodes.EVENT_NONE_CODE = 0;
    eventCodes.EVENT_AROUSAL_CODE = 1;
    eventCodes.EVENT_OBS_APNEA_CODE = 2;
    eventCodes.EVENT_CENTRAL_APNEA_CODE = 3;
    eventCodes.EVENT_HYPOPNEA_CODE = 4;
    eventCodes.EVENT_MIXED_CODE = 5;
    BREATH_CONTINUOUS_SEQUENCE_THRESHOLD = 1; % Breaths have to occur within this period (in sec) to be considered as continuous

    breathTable = BreathDataTableFulls;
    VI = breathTable.VI' * 100; % units are % Eupnea 
    breathTimes = breathTable.Ttot';
    breathStartTimes = breathTable.Time_start';
    breathEndTimes = breathTable.Time_end';
    apneaBreaths = breathTable.ApneaB';
    breathEventType = breathTable.Etype';
    breathSleepState = breathTable.hypnog_B';
    breathPosition = breathTable.pos_B';
    breathStableBreathing = breathTable.StableBreathing';
    zeroValues = (VI == 0);
    nanValues = isnan(VI); % find NaN VI values
    if EXCLUDE_ZERO_VI_BREATHS
        % exclude NaN and 0 and apneaB VI breaths
        nanOrZero = bitor(zeroValues, nanValues);
        apneaFlagLogical = logical(apneaBreaths);
        nanOrZeroOrApnea = bitor(nanOrZero, apneaFlagLogical);
        VI(nanOrZeroOrApnea) = [];
        breathStartTimes(nanOrZeroOrApnea) = [];
        breathEndTimes(nanOrZeroOrApnea) = [];
        breathTimes(nanOrZeroOrApnea) = [];
        breathEventType(nanOrZeroOrApnea) = [];
        breathSleepState(nanOrZeroOrApnea) = [];
        breathPosition(nanOrZeroOrApnea) = [];
        breathStableBreathing(nanOrZeroOrApnea) = [];
    else
        % exclude only NaN VI breaths
        VI(nanValues) = [];
        breathStartTimes(nanValues) = [];
        breathEndTimes(nanValues) = [];
        breathTimes(nanValues) = [];
        breathEventType(nanValues) = [];
        breathSleepState(nanValues) = [];
        breathPosition(nanValues) = [];
        breathStableBreathing(nanValues) = [];
    end

    %% breath duration histogram and outlier removal
    medianDuration = median(breathTimes);
    meanDuration = mean(breathTimes);
    standardDevDuration = std(breathTimes);
    breathCount = size(breathTimes,2);
    if EXCLUDE_LONG_OUTLIER_BREATHS    
        longBreathOutlierIndexes = getLongBreathOutliers(breathTimes);
        beforeRemoveOutlierVI = VI;
        VI(longBreathOutlierIndexes) = []; 
        breathStartTimes(longBreathOutlierIndexes) = [];
        breathEndTimes(longBreathOutlierIndexes) = [];
        breathTimes(longBreathOutlierIndexes) = [];    
        breathEventType(longBreathOutlierIndexes) = [];
        breathSleepState(longBreathOutlierIndexes) = [];
        breathPosition(longBreathOutlierIndexes) = [];
        breathStableBreathing(longBreathOutlierIndexes) = [];
    else
        beforeRemoveOutlierVI = [];
        longBreathOutlierIndexes = [];
    end
    
    hypnogram.states = Evts.Hypnogram;
    hypnogram.times = Evts.Hypnogram_t;
    hypnogram.durations = zeros(size(hypnogram.times)) + 30;
    hypnogram.stateCodes = stateCodes;
    events = Evts.Table1;
    eventAndStageIndexes = getBreathEventAndStageIndexes(hypnogram, events, eventCodes, breathStartTimes);
    
    wakeVI = getVIEvent(stateCodes.SLEEP_STATE_WAKE, hypnogram.states, hypnogram.times, hypnogram.durations, VI, breathStartTimes);
    wakeNoEventVI = getVIWakeNoEvent(eventAndStageIndexes, VI);
    wakeNoEventSurroundedWakeVI = getVIWakeNoEventSurroundedWake(eventAndStageIndexes, VI);
    wakeNoEventSurroundedTwoWakeVI = getVIWakeNoEventSurroundedTwoWake(eventAndStageIndexes, VI);
    sleepVI = getVINotEvent(stateCodes.SLEEP_STATE_WAKE, hypnogram.states, hypnogram.times, hypnogram.durations, VI, breathStartTimes);
    sleepNoEventVI = getVISleepNoEvent(eventAndStageIndexes, VI);
    N1VI = getVIEvent(stateCodes.SLEEP_STATE_N1, hypnogram.states, hypnogram.times, hypnogram.durations, VI, breathStartTimes);
    N2VI = getVIEvent(stateCodes.SLEEP_STATE_N2, hypnogram.states, hypnogram.times, hypnogram.durations, VI, breathStartTimes);
    N3VI = getVIEvent(stateCodes.SLEEP_STATE_N3, hypnogram.states, hypnogram.times, hypnogram.durations, VI, breathStartTimes);
    REMVI = getVIEvent(stateCodes.SLEEP_STATE_REM, hypnogram.states, hypnogram.times, hypnogram.durations, VI, breathStartTimes);
    arousalVI = getVIEvent(eventCodes.EVENT_AROUSAL_CODE, events.EventCodes, events.EventStart, events.EventDuration, VI, breathStartTimes);
    arousalSleepVI = getVIArousalSleep(eventAndStageIndexes, VI);
    arousalSleepNoEventVI = getVIArousalSleepNoEvent(eventAndStageIndexes, VI);
    stableBreathingThreeMinVI = getVIBreathsStableBreathing(3, VI, breathStableBreathing);
    [hypopneaVI,hypopneaVIIndexes] = getVIEvent(eventCodes.EVENT_HYPOPNEA_CODE, events.EventCodes, events.EventStart, events.EventDuration, VI, breathStartTimes);
    [obsApneaVI,obsApeaVIIndexes] = getVIEvent(eventCodes.EVENT_OBS_APNEA_CODE, events.EventCodes, events.EventStart, events.EventDuration, VI, breathStartTimes);
    respEventVI = [hypopneaVI obsApneaVI];
    %% generate respiratory event tables for each event (6 different criteria permutations)
    capVI = true; startStopCriteria = allSettings.startStopCriteria.MANUALLY_SCORED;
    hypopneaVITable.CAPPED_VI_START_STOP_MANUALLY_SCORED = getVIAreaDepthDurationEachEvent(allSettings,hypopneaVIIndexes, VI, breathStartTimes, breathEndTimes, BREATH_CONTINUOUS_SEQUENCE_THRESHOLD, patientNumber, 'hypopneas', capVI, startStopCriteria); % for hypopneas
    apneaVITable.CAPPED_VI_START_STOP_MANUALLY_SCORED = getVIAreaDepthDurationEachEvent(allSettings,obsApeaVIIndexes, VI, breathStartTimes, breathEndTimes, BREATH_CONTINUOUS_SEQUENCE_THRESHOLD, patientNumber, 'apneas', capVI, startStopCriteria); % for obstructive apneas
    capVI = true; startStopCriteria = allSettings.startStopCriteria.UNDER_100_VI;
    hypopneaVITable.CAPPED_VI_START_STOP_UNDER_100_VI = getVIAreaDepthDurationEachEvent(allSettings,hypopneaVIIndexes, VI, breathStartTimes, breathEndTimes, BREATH_CONTINUOUS_SEQUENCE_THRESHOLD, patientNumber, 'hypopneas', capVI, startStopCriteria); % for hypopneas
    apneaVITable.CAPPED_VI_START_STOP_UNDER_100_VI = getVIAreaDepthDurationEachEvent(allSettings,obsApeaVIIndexes, VI, breathStartTimes, breathEndTimes, BREATH_CONTINUOUS_SEQUENCE_THRESHOLD, patientNumber, 'apneas', capVI, startStopCriteria); % for obstructive apneas
    capVI = true; startStopCriteria = allSettings.startStopCriteria.UNDER_70_VI;
    hypopneaVITable.CAPPED_VI_START_STOP_UNDER_70_VI = getVIAreaDepthDurationEachEvent(allSettings,hypopneaVIIndexes, VI, breathStartTimes, breathEndTimes, BREATH_CONTINUOUS_SEQUENCE_THRESHOLD, patientNumber, 'hypopneas', capVI, startStopCriteria); % for hypopneas
    apneaVITable.CAPPED_VI_START_STOP_UNDER_70_VI = getVIAreaDepthDurationEachEvent(allSettings,obsApeaVIIndexes, VI, breathStartTimes, breathEndTimes, BREATH_CONTINUOUS_SEQUENCE_THRESHOLD, patientNumber, 'apneas', capVI, startStopCriteria); % for obstructive apneas
    capVI = false; startStopCriteria = allSettings.startStopCriteria.MANUALLY_SCORED;
    hypopneaVITable.UNCAPPED_VI_START_STOP_MANUALLY_SCORED = getVIAreaDepthDurationEachEvent(allSettings,hypopneaVIIndexes, VI, breathStartTimes, breathEndTimes, BREATH_CONTINUOUS_SEQUENCE_THRESHOLD, patientNumber, 'hypopneas', capVI, startStopCriteria); % for hypopneas
    apneaVITable.UNCAPPED_VI_START_STOP_MANUALLY_SCORED = getVIAreaDepthDurationEachEvent(allSettings,obsApeaVIIndexes, VI, breathStartTimes, breathEndTimes, BREATH_CONTINUOUS_SEQUENCE_THRESHOLD, patientNumber, 'apneas', capVI, startStopCriteria); % for obstructive apneas
    capVI = false; startStopCriteria = allSettings.startStopCriteria.UNDER_100_VI;
    hypopneaVITable.UNCAPPED_VI_START_STOP_UNDER_100_VI = getVIAreaDepthDurationEachEvent(allSettings,hypopneaVIIndexes, VI, breathStartTimes, breathEndTimes, BREATH_CONTINUOUS_SEQUENCE_THRESHOLD, patientNumber, 'hypopneas', capVI, startStopCriteria); % for hypopneas
    apneaVITable.UNCAPPED_VI_START_STOP_UNDER_100_VI = getVIAreaDepthDurationEachEvent(allSettings,obsApeaVIIndexes, VI, breathStartTimes, breathEndTimes, BREATH_CONTINUOUS_SEQUENCE_THRESHOLD, patientNumber, 'apneas', capVI, startStopCriteria); % for obstructive apneas
    capVI = false; startStopCriteria = allSettings.startStopCriteria.UNDER_70_VI;
    hypopneaVITable.UNCAPPED_VI_START_STOP_UNDER_70_VI = getVIAreaDepthDurationEachEvent(allSettings,hypopneaVIIndexes, VI, breathStartTimes, breathEndTimes, BREATH_CONTINUOUS_SEQUENCE_THRESHOLD, patientNumber, 'hypopneas', capVI, startStopCriteria); % for hypopneas
    apneaVITable.UNCAPPED_VI_START_STOP_UNDER_70_VI = getVIAreaDepthDurationEachEvent(allSettings,obsApeaVIIndexes, VI, breathStartTimes, breathEndTimes, BREATH_CONTINUOUS_SEQUENCE_THRESHOLD, patientNumber, 'apneas', capVI, startStopCriteria); % for obstructive apneas    
    combinedVITable = concatVITables(hypopneaVITable, apneaVITable);
    tableFilename = ['Resp event VI tables for patient ' num2str(patientNumber) '.mat'];
    %save(tableFilename, 'hypopneaVITable','apneaVITable');
    %% derive these from breath table events and sleep states - not recommended because can't have multiple events overlapping per breath
    BTarousalSleepNoEventVI = VI((breathEventType == eventCodes.EVENT_AROUSAL_CODE) & (breathSleepState ~= stateCodes.SLEEP_STATE_WAKE));
    BTwakeNoEvent = VI((breathEventType == eventCodes.EVENT_NONE_CODE) & (breathSleepState == stateCodes.SLEEP_STATE_WAKE));
    BTsleepNoEventVI = VI((breathEventType == eventCodes.EVENT_NONE_CODE) & (breathSleepState ~= stateCodes.SLEEP_STATE_WAKE));
    BThypopneaVI = VI(breathEventType == eventCodes.EVENT_HYPOPNEA_CODE);

    sleepVICount1 = size(sleepVI,2);
    sleepVICount2 = size(N1VI,2) + size(N2VI,2) + size(N3VI,2) + size(REMVI,2);
    if (sleepVICount1 ~= sleepVICount2)
        disp(['WARNING: Patient ' num2str(patientNumber) ' sleep breath counts do not match - possible other sleep states in data']);
    end
    sleepVI = [N1VI N2VI N3VI REMVI]; % safter way to calculate sleep VI because not influenced by unknown sleep states (code 8)

    hypopneaCount = sum(events.EventCodes == eventCodes.EVENT_HYPOPNEA_CODE);
    obsApneaCount = sum(events.EventCodes == eventCodes.EVENT_OBS_APNEA_CODE);
    centralApneaCount = sum(events.EventCodes == eventCodes.EVENT_CENTRAL_APNEA_CODE);
    mixedApneaCount = sum(events.EventCodes == eventCodes.EVENT_MIXED_CODE);

    medianVI = median(VI);
    meanVI = mean(VI);
    standardDevVI = std(VI);
    breathCountVI = size(VI,2);

    medianSleepVI = median(sleepVI);
    meanSleepVI = mean(sleepVI);
    PUP_AHI = Evts.AHIdata2.AllSleepAllPAHI(1);
    PUP_HB = Evts.HBtotal;
    disp(['Patient ' num2str(patientNumber) ': median sleep VI = ' num2str(medianSleepVI) ', mean sleep VI = ' num2str(meanSleepVI) ', AHI = ' num2str(PUP_AHI) ', total hypopneas = ' num2str(hypopneaCount) ', total obstructive apneas = ' num2str(obsApneaCount) ', total central apneas = ' num2str(centralApneaCount) ', total mixed apneas = ' num2str(mixedApneaCount) ', longest breath = ' num2str(max(breathTimes)), ', outlier breaths removed = ', num2str(size(longBreathOutlierIndexes,2))]);

    %% Metric calculation for table
    N1 = sum(hypnogram.states == hypnogram.stateCodes.SLEEP_STATE_N1);
    N2 = sum(hypnogram.states == hypnogram.stateCodes.SLEEP_STATE_N2);
    N3 = sum(hypnogram.states == hypnogram.stateCodes.SLEEP_STATE_N3);
    REM = sum(hypnogram.states == hypnogram.stateCodes.SLEEP_STATE_REM);
    TST_sec = (N1 + N2 + N3 + REM) * 30;

    %% time in seconds under x % of eupnea
    VITimeUnderT = generateTimeUnderMetrics(VI,breathTimes,[100, 90, 80, 70, 60, 50, 40, 30, 20, 10]);

    %% THIS IS WHERE CAN CHANGE THE SELECTION CRITERIA I.E. WHICH OF 6 CRITERIA PERMUTATIONS ARE USED
    % ON 2022-02-03 PHIL SAID HE PREFERS
    % UNCAPPED_VI_START_STOP_UNDER_100_VI AND UNCAPPED_VI_START_STOP_UNDER_70_VI
    sleepVIT = getMetricsParameter(sleepVI,TST_sec);
    respEventVIT = getMetricsParameter(respEventVI,TST_sec);
    VIEventArea = combinedVITable.UNCAPPED_VI_START_STOP_UNDER_100_VI.EventArea;
    VIEventAreaT = getMetricsParameter(VIEventArea,TST_sec);
    VIEventDepth = combinedVITable.UNCAPPED_VI_START_STOP_UNDER_100_VI.EventDepth;
    VIEventDepthT = getMetricsParameter(VIEventDepth,TST_sec);
    VIEventDuration = combinedVITable.UNCAPPED_VI_START_STOP_UNDER_100_VI.EventDuration;
    VIEventDurationT = getMetricsParameter(VIEventDuration,TST_sec);
    VIEventBreathCount = combinedVITable.UNCAPPED_VI_START_STOP_UNDER_100_VI.EventBreathCount;
    VIEventBreathCountT = getMetricsParameter(VIEventBreathCount,TST_sec);
    output.sleepVIT = sleepVIT;
    output.respEventVIT = respEventVIT;
    output.VIEventAreaT = VIEventAreaT;
    output.VIEventDepthT = VIEventDepthT;
    output.VIEventDurationT = VIEventDurationT;
    output.VIEventBreathCountT = VIEventBreathCountT;
    output.VITimeUnderT = VITimeUnderT;
    output.PUP_AHI = PUP_AHI;
    output.PUP_HB = PUP_HB;
end

function output = generateTimeUnderMetrics(VI, breathTimes, levels)
    for level = levels
        indexes = (VI < level);
        timeUnder = sum(breathTimes(indexes)); % time in seconds under level % eupnea
        eval(sprintf('output.T%sVI = timeUnder;',num2str(level)));
    end
end

function output = getMetricsParameter(data,TST_sec)
    output.median = median(data);
    output.mean = mean(data);
    output.min = min(data);
    output.max = max(data);
    output.stdDev = std(data);
    output.IQR = iqr(data);
    output.sumNormTST = sum(data) / TST_sec;
end

function output = concatVITables(hypopneaVITable, apneaVITable)
    output.CAPPED_VI_START_STOP_MANUALLY_SCORED = table;
    output.CAPPED_VI_START_STOP_UNDER_100_VI  = table;
    output.CAPPED_VI_START_STOP_UNDER_70_VI  = table;
    output.UNCAPPED_VI_START_STOP_MANUALLY_SCORED  = table;
    output.UNCAPPED_VI_START_STOP_UNDER_100_VI  = table;
    output.UNCAPPED_VI_START_STOP_UNDER_70_VI  = table;
    if ~isempty(hypopneaVITable.CAPPED_VI_START_STOP_MANUALLY_SCORED)
        output.CAPPED_VI_START_STOP_MANUALLY_SCORED = [output.CAPPED_VI_START_STOP_MANUALLY_SCORED ; hypopneaVITable.CAPPED_VI_START_STOP_MANUALLY_SCORED];
    end
    if ~isempty(apneaVITable.CAPPED_VI_START_STOP_MANUALLY_SCORED)
        output.CAPPED_VI_START_STOP_MANUALLY_SCORED = [output.CAPPED_VI_START_STOP_MANUALLY_SCORED ; apneaVITable.CAPPED_VI_START_STOP_MANUALLY_SCORED];
    end
    if ~isempty(hypopneaVITable.CAPPED_VI_START_STOP_UNDER_100_VI) 
        output.CAPPED_VI_START_STOP_UNDER_100_VI = [output.CAPPED_VI_START_STOP_UNDER_100_VI ; hypopneaVITable.CAPPED_VI_START_STOP_UNDER_100_VI];
    end
    if ~isempty(apneaVITable.CAPPED_VI_START_STOP_UNDER_100_VI) 
        output.CAPPED_VI_START_STOP_UNDER_100_VI = [output.CAPPED_VI_START_STOP_UNDER_100_VI ; apneaVITable.CAPPED_VI_START_STOP_UNDER_100_VI];
    end
    if ~isempty(hypopneaVITable.CAPPED_VI_START_STOP_UNDER_70_VI) 
        output.CAPPED_VI_START_STOP_UNDER_70_VI = [output.CAPPED_VI_START_STOP_UNDER_70_VI ; hypopneaVITable.CAPPED_VI_START_STOP_UNDER_70_VI];
    end
    if ~isempty(apneaVITable.CAPPED_VI_START_STOP_UNDER_70_VI) 
        output.CAPPED_VI_START_STOP_UNDER_70_VI = [output.CAPPED_VI_START_STOP_UNDER_70_VI ; apneaVITable.CAPPED_VI_START_STOP_UNDER_70_VI];
    end
    if ~isempty(hypopneaVITable.UNCAPPED_VI_START_STOP_MANUALLY_SCORED) 
        output.UNCAPPED_VI_START_STOP_MANUALLY_SCORED = [output.UNCAPPED_VI_START_STOP_MANUALLY_SCORED ; hypopneaVITable.UNCAPPED_VI_START_STOP_MANUALLY_SCORED];
    end
    if ~isempty(apneaVITable.UNCAPPED_VI_START_STOP_MANUALLY_SCORED) 
        output.UNCAPPED_VI_START_STOP_MANUALLY_SCORED = [output.UNCAPPED_VI_START_STOP_MANUALLY_SCORED ; apneaVITable.UNCAPPED_VI_START_STOP_MANUALLY_SCORED];
    end
    if ~isempty(hypopneaVITable.UNCAPPED_VI_START_STOP_UNDER_100_VI) 
        output.UNCAPPED_VI_START_STOP_UNDER_100_VI = [output.UNCAPPED_VI_START_STOP_UNDER_100_VI ; hypopneaVITable.UNCAPPED_VI_START_STOP_UNDER_100_VI];
    end
    if ~isempty(apneaVITable.UNCAPPED_VI_START_STOP_UNDER_100_VI)
        output.UNCAPPED_VI_START_STOP_UNDER_100_VI = [output.UNCAPPED_VI_START_STOP_UNDER_100_VI ; apneaVITable.UNCAPPED_VI_START_STOP_UNDER_100_VI];
    end
    if ~isempty(hypopneaVITable.UNCAPPED_VI_START_STOP_UNDER_70_VI) 
        output.UNCAPPED_VI_START_STOP_UNDER_70_VI = [output.UNCAPPED_VI_START_STOP_UNDER_70_VI ; hypopneaVITable.UNCAPPED_VI_START_STOP_UNDER_70_VI];
    end
    if ~isempty(apneaVITable.UNCAPPED_VI_START_STOP_UNDER_70_VI) 
        output.UNCAPPED_VI_START_STOP_UNDER_70_VI = [output.UNCAPPED_VI_START_STOP_UNDER_70_VI ; apneaVITable.UNCAPPED_VI_START_STOP_UNDER_70_VI];
    end
end

function longBreathOutlierIndexes = getLongBreathOutliers(breathTimes)
    quartiles = quantile(breathTimes,3);
    Q3 = quartiles(3);
    IQR = iqr(breathTimes);
    durationLimit = Q3 + (1.5 * IQR);
    longBreathOutlierIndexes = find(breathTimes > durationLimit);
end

%% get ventilation breaths during a specified event type
function [eventVI, indexes] = getVIEvent(eventCode, events, eventStarts, eventDurations, VI, breathStartTimes)
    indexes = find(events == eventCode);
    starts = eventStarts(indexes);
    ends = starts + eventDurations(indexes);
    indexes = (breathStartTimes > starts) & (breathStartTimes < ends);
    indexes = sum(indexes, 1);
    indexes = find(indexes > 0);
    eventVI = VI(indexes);
end

%% get ventilation breaths in all events other than this event
function eventVI = getVINotEvent(eventCode, events, eventStarts, eventDurations, VI, breathStartTimes)
    indexes = find(events ~= eventCode);
    starts = eventStarts(indexes);
    ends = starts + eventDurations(indexes);
    indexes = (breathStartTimes > starts) & (breathStartTimes < ends);
    indexes = sum(indexes, 1);
    indexes = find(indexes > 0);
    eventVI = VI(indexes);
end

%% get ventilation breaths during arousals during sleep epochs
function arousalSleepVI = getVIArousalSleep(eventAndStageIndexes, VI)
    % indexes during sleep
    sleepIndexes = eventAndStageIndexes.sleepIndexes;
    % indexes during arousal
    arousalIndexes = eventAndStageIndexes.ArousalIndexes;
    % indexes during sleep and arousal
    sleepANDarousalIndexes = intersect(arousalIndexes, sleepIndexes);
    arousalSleepVI = VI(sleepANDarousalIndexes);
end

%% get ventilation breaths during arousals during sleep epochs when no other events occuring
function arousalSleepNoEventVI = getVIArousalSleepNoEvent(eventAndStageIndexes, VI)
    % indexes during sleep
    sleepIndexes = eventAndStageIndexes.sleepIndexes;
    % indexes during arousal
    arousalIndexes = eventAndStageIndexes.ArousalIndexes;
    % indexes during all other events
    eventIndexesExcludingArousal = union(eventAndStageIndexes.ObsApneaIndexes, union(eventAndStageIndexes.CentralApneaIndexes, union(eventAndStageIndexes.HypopneaIndexes, eventAndStageIndexes.MixedEventIndexes)));
    % indexes during arousal but not other events
    arousalANDOtherEvents = intersect(arousalIndexes, eventIndexesExcludingArousal);
    arousalNoOtherEvent = setxor(arousalIndexes, arousalANDOtherEvents);
    % indexes during sleep and not other events (but includes arousals)
    sleepANDOtherEvents = intersect(sleepIndexes, eventIndexesExcludingArousal);
    sleepNoOtherEvent = setxor(sleepIndexes, sleepANDOtherEvents);
    % indexes during arousals and sleep and no other events
    arousalANDSleepNoOtherEventIndexes = intersect(sleepNoOtherEvent, arousalNoOtherEvent);
    arousalSleepNoEventVI = VI(arousalANDSleepNoOtherEventIndexes);
end

function sleepNoEventVI = getVISleepNoEvent(eventAndStageIndexes, VI)
    % indexes during sleep
    sleepIndexes = eventAndStageIndexes.sleepIndexes;
    % indexes during events
    eventIndexes = eventAndStageIndexes.eventIndexes;
    % indexes during sleep but not events
    sleepANDEventIndexes = intersect(eventIndexes, sleepIndexes);
    sleepNoEventIndexes = setxor(sleepANDEventIndexes, sleepIndexes);
    sleepNoEventVI = VI(sleepNoEventIndexes);
end

function wakeNoEventVI = getVIWakeNoEvent(eventAndStageIndexes, VI)
    % indexes during wake
    wakeIndexes = eventAndStageIndexes.WakeIndexes;
    % indexes during events
    eventIndexes = eventAndStageIndexes.eventIndexes;
    % indexes during wake but not events
    wakeANDEventIndexes = intersect(eventIndexes, wakeIndexes);
    wakeNoEventIndexes = setxor(wakeANDEventIndexes, wakeIndexes);
    wakeNoEventVI = VI(wakeNoEventIndexes);
end

function wakeNoEventSurroundedWakeVI = getVIWakeNoEventSurroundedWake(eventAndStageIndexes, VI)
    % indexes during wake
    wakeIndexes = eventAndStageIndexes.WakeIndexes;
    % indexes during events
    eventIndexes = eventAndStageIndexes.eventIndexes;
    % indexes during wake but not events
    wakeANDEventIndexes = intersect(eventIndexes, wakeIndexes);
    wakeNoEventIndexes = setxor(wakeANDEventIndexes, wakeIndexes);
    surroundedWake = zeros(size(wakeNoEventIndexes));
    % ignore the first and last indexes because they can't possibly
    % be surrounded by wake
    for i = 2 : (size(wakeNoEventIndexes,2) - 1)
        % check if this index is surrounded by wake on either side
        if ((wakeNoEventIndexes(i - 1) == (wakeNoEventIndexes(i) - 1)) && ((wakeNoEventIndexes(i + 1) == (wakeNoEventIndexes(i) + 1))))
            surroundedWake(i) = 1;
        end
    end
    surroundedWakeIndexes = wakeNoEventIndexes(find(surroundedWake == 1));
    wakeNoEventSurroundedWakeVI = VI(surroundedWakeIndexes);
end

function wakeNoEventSurroundedTwoWakeVI = getVIWakeNoEventSurroundedTwoWake(eventAndStageIndexes, VI)
    % indexes during wake
    wakeIndexes = eventAndStageIndexes.WakeIndexes;
    % indexes during events
    eventIndexes = eventAndStageIndexes.eventIndexes;
    % indexes during wake but not events
    wakeANDEventIndexes = intersect(eventIndexes, wakeIndexes);
    wakeNoEventIndexes = setxor(wakeANDEventIndexes, wakeIndexes);
    surroundedWake = zeros(size(wakeNoEventIndexes));
    % ignore the first and last indexes because they can't possibly
    % be surrounded by wake
    for i = 3 : (size(wakeNoEventIndexes,2) - 2)
        % check if this index is surrounded by wake on either side
        if ((wakeNoEventIndexes(i - 1) == (wakeNoEventIndexes(i) - 1)) && ((wakeNoEventIndexes(i + 1) == (wakeNoEventIndexes(i) + 1))))
            if ((wakeNoEventIndexes(i - 2) == (wakeNoEventIndexes(i) - 2)) && ((wakeNoEventIndexes(i + 2) == (wakeNoEventIndexes(i) + 2))))
                surroundedWake(i) = 1;
            end
        end
    end
    surroundedWakeIndexes = wakeNoEventIndexes(find(surroundedWake == 1));
    wakeNoEventSurroundedTwoWakeVI = VI(surroundedWakeIndexes);
end

function indexes = getBreathIndicesEventOrStage(data, code, starts, durations, breathStartTimes)
    indexes = find(data == code);
    starts = starts(indexes);
    ends = starts + durations(indexes);
    indexes = (breathStartTimes >= starts) & (breathStartTimes <= ends);
    indexes = sum(indexes, 1);
    indexes = find(indexes > 0);
end

function stableBreathingMinVI = getVIBreathsStableBreathing(threshold, VI, breathStableBreathing)
    stableBreathingMinVI = VI(breathStableBreathing >= threshold);
end

%% get indexes of wake, sleep stages, respiratory events, arousals
function eventAndStageIndexes = getBreathEventAndStageIndexes(hypnogram, events, eventCodes, breathStartTimes)
    eventAndStageIndexes.WakeIndexes = getBreathIndicesEventOrStage(hypnogram.states, hypnogram.stateCodes.SLEEP_STATE_WAKE, hypnogram.times, hypnogram.durations, breathStartTimes);    
    eventAndStageIndexes.N1Indexes = getBreathIndicesEventOrStage(hypnogram.states, hypnogram.stateCodes.SLEEP_STATE_N1, hypnogram.times, hypnogram.durations, breathStartTimes);
    eventAndStageIndexes.N2Indexes = getBreathIndicesEventOrStage(hypnogram.states, hypnogram.stateCodes.SLEEP_STATE_N2, hypnogram.times, hypnogram.durations, breathStartTimes);
    eventAndStageIndexes.N3Indexes = getBreathIndicesEventOrStage(hypnogram.states, hypnogram.stateCodes.SLEEP_STATE_N3, hypnogram.times, hypnogram.durations, breathStartTimes);
    eventAndStageIndexes.REMIndexes = getBreathIndicesEventOrStage(hypnogram.states, hypnogram.stateCodes.SLEEP_STATE_REM, hypnogram.times, hypnogram.durations, breathStartTimes);
    eventAndStageIndexes.ArousalIndexes = getBreathIndicesEventOrStage(events.EventCodes, eventCodes.EVENT_AROUSAL_CODE, events.EventStart, events.EventDuration, breathStartTimes);
    eventAndStageIndexes.ObsApneaIndexes = getBreathIndicesEventOrStage(events.EventCodes, eventCodes.EVENT_OBS_APNEA_CODE, events.EventStart, events.EventDuration, breathStartTimes);
    eventAndStageIndexes.CentralApneaIndexes = getBreathIndicesEventOrStage(events.EventCodes, eventCodes.EVENT_CENTRAL_APNEA_CODE, events.EventStart, events.EventDuration, breathStartTimes);
    eventAndStageIndexes.HypopneaIndexes = getBreathIndicesEventOrStage(events.EventCodes, eventCodes.EVENT_HYPOPNEA_CODE, events.EventStart, events.EventDuration, breathStartTimes);
    eventAndStageIndexes.MixedEventIndexes = getBreathIndicesEventOrStage(events.EventCodes, eventCodes.EVENT_MIXED_CODE, events.EventStart, events.EventDuration, breathStartTimes);
    eventAndStageIndexes.sleepIndexes = union(union(union(eventAndStageIndexes.N1Indexes, eventAndStageIndexes.N2Indexes), eventAndStageIndexes.N3Indexes), eventAndStageIndexes.REMIndexes);
    eventAndStageIndexes.eventIndexes = union(eventAndStageIndexes.ArousalIndexes, union(eventAndStageIndexes.ObsApneaIndexes, union(eventAndStageIndexes.CentralApneaIndexes, union(eventAndStageIndexes.HypopneaIndexes, eventAndStageIndexes.MixedEventIndexes))));
end



%% can handle hypopneas or apneas as input
function eventTable = getVIAreaDepthDurationEachEvent(allSettings,respEventVIIndexes, VI, breathStartTimes, breathEndTimes, BREATH_CONTINUOUS_SEQUENCE_THRESHOLD, patientNumber, eventType, capVI, startStopCriteria)
    eventTable = [];
    eventBreathCount = size(respEventVIIndexes,2); 
    if (eventBreathCount > 0)
        consecutive = [0 diff(respEventVIIndexes)];
        eventStarts = [respEventVIIndexes(1)];
        eventEnds = [];
        for i = 2:eventBreathCount
            if (consecutive(i) ~= 1)
                eventEnds = [eventEnds respEventVIIndexes(i-1)];
                eventStarts = [eventStarts respEventVIIndexes(i)];
            end
        end
        eventEnds = [eventEnds respEventVIIndexes(end)];
        eventLengths = eventEnds - eventStarts + 1;
        usableEvents = 0;
        totalEvents = size(eventLengths,2);
        for eventNumber = 1:totalEvents
            breathIndexesDuringEvent = eventStarts(eventNumber):eventEnds(eventNumber);
            if (breathIndexesDuringEvent(1) < 1) || (breathIndexesDuringEvent(end) > size(VI,2))
                continue; % skip this event because indexes are out of bounds
            end
            % how do I know the breaths before and after
            % the event are continuous (in sequence)? 
            % answer: by ensuring breath n start == breath n-1 end
            startTimes = breathStartTimes(breathIndexesDuringEvent);
            startTimes(1) = [];
            endTimes = breathEndTimes(breathIndexesDuringEvent);
            endTimes(end) = [];
            difference = abs(startTimes - endTimes);
            difference = difference <= BREATH_CONTINUOUS_SEQUENCE_THRESHOLD;
            if ~all(difference)
                continue; % skip this event because breaths aren't continuous (not in sequence)
            end
            if (startStopCriteria ~= allSettings.startStopCriteria.MANUALLY_SCORED)
                if (startStopCriteria == allSettings.startStopCriteria.UNDER_100_VI)
                    thresholdVI = 100;
                elseif (startStopCriteria == allSettings.startStopCriteria.UNDER_70_VI)
                    thresholdVI = 70;
                end
                % first breath is the first under the threshold 
                foundStart = false;
                VIDuringEvent = VI(breathIndexesDuringEvent);
                for index = 1:size(breathIndexesDuringEvent,2)
                    if (VIDuringEvent(index) < thresholdVI)
                        foundStart = true;
                        breathIndexesDuringEvent = breathIndexesDuringEvent(index:end);
                        break;
                    end
                end
                if ~foundStart
                    continue; % no breaths are under the threshold so skip breath
                end
                % last breath is the last under the threshold
                % if up to here it's impossible not to find something
                % because at least 1 breath is under the threshold
                % from the start index search above
                VIDuringEvent = VI(breathIndexesDuringEvent);
                for index = size(breathIndexesDuringEvent,2):-1:1
                    if (VIDuringEvent(index) < thresholdVI)
                        breathIndexesDuringEvent = breathIndexesDuringEvent(1:index);
                        break;
                    end
                end
            end
            VIDuringEvent = VI(breathIndexesDuringEvent); % recalculate this based on new indexes
            VIDuringEventCapped = VIDuringEvent;
            if capVI
                % cap VI to 100% during events
                VIDuringEventCapped(VIDuringEvent > 100) = 100;
            end
            eventBreathCount = size(breathIndexesDuringEvent,2);
            eventDepth = 100 - min(VIDuringEventCapped); % depth from eupnea
            breathDurations = breathEndTimes(breathIndexesDuringEvent) - breathStartTimes(breathIndexesDuringEvent);
            eventDuration = sum(breathDurations);
            eventArea = sum(breathDurations .* (100 - VIDuringEventCapped)); % lost ventilation from eupnea
            eventTable = [eventTable ; [eventArea, eventDepth, eventDuration, eventBreathCount]];
            usableEvents = usableEvents + 1;
            % plot breath ventilation for each usable event
            if 0
                startTimes = breathStartTimes(breathIndexesDuringEvent);
                eventStart = startTimes(1);
                startTimes = startTimes - eventStart;
                endTimes = breathEndTimes(breathIndexesDuringEvent);
                endTimes = endTimes - eventStart;
                totalPlotTimes = [startTimes endTimes(end)];
                totalPlotVI = [VIDuringEvent VIDuringEvent(end)];
                totalPlotVICapped = [VIDuringEventCapped VIDuringEventCapped(end)];
                fig = figure; hold on; fig.WindowState = 'maximized';
                stairs(totalPlotTimes,totalPlotVICapped,'-r'); % to get a line that continues to the end of the last breath
                stairs(startTimes,VIDuringEventCapped,'-ro');
                yline(100,'k--');
                yline(0,'k--');
                ylabel('VI (% eupnea)');
                ylim([-10 110]);
                xlabel('Time (sec)')
                xlim([startTimes(1) endTimes(end)]);
                plotTitle = sprintf('Event area = %.3f s%%, event depth = %.3f %%, event duration = %.3f s, event breath count = %d', eventArea, eventDepth, eventDuration, eventBreathCount);
                title(plotTitle);
                input('Press Enter for next plot');
                close(fig);
            end
        end
        text = sprintf("Patient %d: %d/%d %s were usable for tabulation (%d were unusable)", patientNumber, usableEvents,totalEvents,eventType,totalEvents-usableEvents); disp(text);
    else
        disp("No event breaths were found.");
    end
    if (~isempty(eventTable))
        eventTable = array2table(eventTable);
        eventTable.Properties.VariableNames = {'EventArea' 'EventDepth' 'EventDuration' 'EventBreathCount'};
    else
        eventTable = table;
    end
end

