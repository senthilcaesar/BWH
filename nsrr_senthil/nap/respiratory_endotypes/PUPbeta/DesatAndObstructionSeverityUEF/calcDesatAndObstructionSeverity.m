% EAS added on 2022-01-11
function outputMetrics = calcDesatAndObstructionSeverity(patientID, pupEvents, pupStartTime, pupHypnogram, spo2SampleRate, spo2)
    %% define exclusion criteria
    EXCLUDE_ALL_WAKE_DESATS = 1;
    EXCLUDE_ALL_WAKE_EVENTS = 2;
    EXCLUDE_ALL_EVENTS_BEFORE_SLEEP_ONSET = 3;
    EXCLUDE_ALL_WAKE_EVENTS_EXCEPT_SLEEP_WAKE_TRANSITION = 4;
    %% create table in the same style as STAR event_table.mat
    % Need columns for {'Desaturations'}    {'Apneas'}    {'Hypopneas'}    {'Arousals'}
    % Desaturations and Arousals columns each contain a whatever by 2
    % double (start time and duration
    % Apneas and Hypopneas are the same except they contain 3 column double
    % and the last column is either 0 for hypopnea or 1 for apnea    
    adjustedStartTimes = pupEvents.EventStart - pupStartTime;
    eventDurations = pupEvents.EventDuration;
    eventStartsDurations = [adjustedStartTimes, eventDurations];    
    eventNames = pupEvents.EventName;
    arousalIndexes = find(strcmp(eventNames, 'Arousal (ARO PLM)') | strcmp(eventNames, 'Arousal (ARO RES)') | strcmp(eventNames, 'Arousal (ARO SPONT)') | strcmp(eventNames, 'Arousal (ARO Limb)'));
    hypopneaIndexes = find(strcmp(eventNames, 'Hypopnea'));
    desaturationIndexes = find(strcmp(eventNames, 'SpO2 desaturation'));
    apneaIndexes = find(strcmp(eventNames, 'Central Apnea') | strcmp(eventNames, 'Mixed Apnea') | strcmp(eventNames, 'Obstructive Apnea'));
    Arousals = eventStartsDurations(arousalIndexes,:);
    Arousals = {Arousals};
    Hypopneas = eventStartsDurations(hypopneaIndexes,:);
    Hypopneas = [Hypopneas, Hypopneas(:,1) * 0];
    Hypopneas = {Hypopneas};
    Desaturations = eventStartsDurations(desaturationIndexes,:);
    Desaturations = {Desaturations};
    Apneas = eventStartsDurations(apneaIndexes,:);
    Apneas = [Apneas, Apneas(:,1) * 0 + 1];
    Apneas = {Apneas};
    pupEventTable = table(Desaturations, Apneas, Hypopneas, Arousals);
    %% Convert between ProFusion and Terrill sleep state codes
    TERRILL_WAKE = 4;
    TERRILL_REM = 3;
    TERRILL_N1 = 2;
    TERRILL_N2 = 1;
    TERRILL_N3 = 0;
    TERRILL_UNKNOWN_BAD = 8;
    badCount = sum(pupHypnogram == TERRILL_UNKNOWN_BAD);
    if (badCount > 0)
        text = sprintf('Warning - this patient has %d unknown/bad scored sleep epochs which have been converted to wake epochs in OSAseverity analysis.', badCount); disp(text);
        pupHypnogram(pupHypnogram == TERRILL_UNKNOWN_BAD) = TERRILL_WAKE;
    end
    STAR_WAKE = 0;
    STAR_REM = 4;
    STAR_N1 = 1;
    STAR_N2 = 2;
    STAR_N3 = 3;
    hypnogram = pupHypnogram;
    hypnogram(pupHypnogram == TERRILL_WAKE) = STAR_WAKE;
    hypnogram(pupHypnogram == TERRILL_REM) = STAR_REM;
    hypnogram(pupHypnogram == TERRILL_N1) = STAR_N1;
    hypnogram(pupHypnogram == TERRILL_N2) = STAR_N2;
    hypnogram(pupHypnogram == TERRILL_N3) = STAR_N3;
    nanValues = isnan(spo2); % find NaN values
    spo2(nanValues) = 0; % convert NaN SpO2 values to 0        
    [patient_id, TST_h, TST_min, TSET_h, TSET_min, AHI, ODI, AI, ObsDur, DesDur, DesSev, ObsSev, ave_spo2, min_spo2, t90, ave_depth, med_depth] = DesatAndObstructionSeverity(patientID, pupEventTable, hypnogram, spo2, spo2SampleRate, EXCLUDE_ALL_WAKE_DESATS);
    outputMetrics.EXCLUDE_ALL_WAKE_DESATS = [str2num(patient_id), TST_h, TST_min, TSET_h, TSET_min, AHI, ODI, AI, ObsDur, DesDur, DesSev, ObsSev, ave_spo2, min_spo2, t90, ave_depth, med_depth];
    [patient_id, TST_h, TST_min, TSET_h, TSET_min, AHI, ODI, AI, ObsDur, DesDur, DesSev, ObsSev, ave_spo2, min_spo2, t90, ave_depth, med_depth] = DesatAndObstructionSeverity(patientID, pupEventTable, hypnogram, spo2, spo2SampleRate, EXCLUDE_ALL_WAKE_EVENTS);
    outputMetrics.EXCLUDE_ALL_WAKE_EVENTS = [str2num(patient_id), TST_h, TST_min, TSET_h, TSET_min, AHI, ODI, AI, ObsDur, DesDur, DesSev, ObsSev, ave_spo2, min_spo2, t90, ave_depth, med_depth];
    [patient_id, TST_h, TST_min, TSET_h, TSET_min, AHI, ODI, AI, ObsDur, DesDur, DesSev, ObsSev, ave_spo2, min_spo2, t90, ave_depth, med_depth] = DesatAndObstructionSeverity(patientID, pupEventTable, hypnogram, spo2, spo2SampleRate, EXCLUDE_ALL_EVENTS_BEFORE_SLEEP_ONSET);
    outputMetrics.EXCLUDE_ALL_EVENTS_BEFORE_SLEEP_ONSET = [str2num(patient_id), TST_h, TST_min, TSET_h, TSET_min, AHI, ODI, AI, ObsDur, DesDur, DesSev, ObsSev, ave_spo2, min_spo2, t90, ave_depth, med_depth]; 
    [patient_id, TST_h, TST_min, TSET_h, TSET_min, AHI, ODI, AI, ObsDur, DesDur, DesSev, ObsSev, ave_spo2, min_spo2, t90, ave_depth, med_depth] = DesatAndObstructionSeverity(patientID, pupEventTable, hypnogram, spo2, spo2SampleRate, EXCLUDE_ALL_WAKE_EVENTS_EXCEPT_SLEEP_WAKE_TRANSITION);
    outputMetrics.EXCLUDE_ALL_WAKE_EVENTS_EXCEPT_SLEEP_WAKE_TRANSITION = [str2num(patient_id), TST_h, TST_min, TSET_h, TSET_min, AHI, ODI, AI, ObsDur, DesDur, DesSev, ObsSev, ave_spo2, min_spo2, t90, ave_depth, med_depth];
    outputMetrics.headers = {'patient_id','TST_h','TST_min','TSET_h','TSET_min','AHI','ODI','AI','ObsDur','DesDur','DesSev','ObsSev','ave_spo2','min_spo2','t90','ave_depth','med_depth'};
end % function