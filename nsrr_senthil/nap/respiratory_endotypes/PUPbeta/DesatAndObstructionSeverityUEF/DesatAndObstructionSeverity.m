% EAS added on 2022-01-11. Originally from Samu Kainulainen and modified by EAS
% Used in these publications among others:
% Severe desaturations increase psychomotor vigilance task-based median reaction time and number of lapses in obstructive sleep apnoea patients
% DOI: 10.1183/13993003.01849-2019
% Severity of Desaturations Reflects OSA-Related Daytime Sleepiness Better Than AHI
% DOI: 10.5664/jcsm.7806
function [patient_id TST_h TST_min TSET_h TSET_min AHI ODI AI ObsDur DesDur DesSev ObsSev ave_spo2...
     min_spo2 t90 ave_depth med_depth] = DesatAndObstructionSeverity(patient_id,patient_eventinfo,hypnogram,...
                                    Spo2,Spo2_fs,exclusionCriteria)

%   Inputs:
%       - Row from event table including patient number and necessary event info
%       - SpO2 and sampling frequency
%   Outputs:
%       - Vector comprising all parameters for estimating the severity of OSA
% 
% the linking is done using 60 second delay threshold between resp event 
% and desat -> if no desat occurs in that timeframe, the resp event is not
% considered to relate to any desat. meaning that this part excludes 
% automatically the arousal related apneas and hypopneas or apneas without 
% desat/arousal
% 
% the eventlink matrix construction takes care of resp events in numerical 
% order, so that they do not overlap or mix, and for example no desat can 
% start before resp event etc.

EVENT_TYPE_DESATS = 1;
EVENT_TYPE_RESPIRATORY = 2;
EVENT_TYPE_AROUSAL = 3;

% compute total sleep times from hypnogram
TST_s = sum(hypnogram == 1 | hypnogram == 2 | hypnogram == 3 | hypnogram == 4)*30;
TST_min = TST_s / 60;
TST_h = TST_min / 60;
% verify the patient and print it out
patient_id = patient_id;
% Compute number of arousals and their duration
arousal_info = patient_eventinfo.Arousals{1,1};
if (isempty(arousal_info))
    AI = 0;
    ArDur = 0;
else
    arousal_info = removeEventsBeforeSleepOnset(hypnogram, arousal_info, EVENT_TYPE_AROUSAL, exclusionCriteria); % EAS added
    n_arousals = size(arousal_info,1);
    arousal_duration = sum(arousal_info(:,2));
    % compute arousal parameters
    AI = n_arousals/TST_h;
    ArDur = (arousal_duration/TST_s)*100;
end

% extract respiratory event and desat info
hypopnea_info = patient_eventinfo.Hypopneas{1,1};
if (~isempty(hypopnea_info))
    hypopnea_info_orig = hypopnea_info;
    hypopnea_info = removeEventsBeforeSleepOnset(hypnogram, hypopnea_info, EVENT_TYPE_RESPIRATORY, exclusionCriteria); % EAS added
end
apnea_info = patient_eventinfo.Apneas{1,1};
if (~isempty(apnea_info))
    apnea_info_orig = apnea_info;
    apnea_info = removeEventsBeforeSleepOnset(hypnogram, apnea_info, EVENT_TYPE_RESPIRATORY, exclusionCriteria); % EAS added
end
respiratory_info = [hypopnea_info;apnea_info];

if (isempty(respiratory_info))
    AHI = 0;
    ObsDur = 0;
    ObsSev = 0;
else
    % compute number of resp events from event table
    n_hypopnea = size(hypopnea_info, 1);
    n_apnea = size(apnea_info, 1);
    % compute AHI
    AHI = (n_hypopnea+n_apnea)/TST_h;
    % sort respiratory events based on their starting time
    [~,sortteri] = sort(respiratory_info(:,1));
    respiratory_info = respiratory_info(sortteri,:);
    % compute durations for resp and desats
    respiratory_duration = sum(respiratory_info(:,2));
    % compute obsdur
    ObsDur = (respiratory_duration/TST_s)*100;
end

% desat checkpoint
desat_info = patient_eventinfo.Desaturations{1,1};
if (isempty(desat_info))
    ODI = 0;
    DesDur = 0;
    DesSev = 0;
    ObsSev = 0;
    ave_depth = 0;
    med_depth = 0;
    ave_spo2 = mean(Spo2);
    okt = find(Spo2 > 10);
    min_spo2 = min(Spo2(okt));
    t90 = size(find(Spo2(okt) < 90),1)/Spo2_fs;
else
    if (~isempty(desat_info))
        desat_info_orig = desat_info;
        desat_info = removeEventsBeforeSleepOnset(hypnogram, desat_info, EVENT_TYPE_DESATS, exclusionCriteria); % EAS added
    end
    
    n_desat = size(desat_info,1);
    
    ODI = n_desat/TST_h;
    
    % compute desat durations and DesDur
    desat_duration = sum(desat_info(:,2));
    DesDur = (desat_duration/TST_s)*100;
    
    % time vectors for signals
    N_spo2 = size(Spo2,1);
    t_spo2 = 0:(1/Spo2_fs):(N_spo2-1)/Spo2_fs;
    % conventional SPO2 parameters
    ave_spo2 = mean(Spo2);
    okt = find(Spo2 > 10); % erase artifactant part of SpO2
    min_spo2 = min(Spo2(okt));
    t90 = size(find(Spo2(okt) < 90),1)/Spo2_fs;
    
    
    % detect desaturation events from the signal
    if (isempty(desat_info) == 1) % if no desats, parameters set to 0
        DesSev = 0;
        ave_depth = 0;
        med_depth = 0;
    else
        for ii = 1:size(desat_info,1)
            % indicate individual desaturation starting and ending time
            desat_start = desat_info(ii,1);
            desat_stop = desat_start + desat_info(ii,2);
            % find indeces from the time vector
            start_index = find(ceil(desat_start) == t_spo2);
            stop_index = find(ceil(desat_stop) == t_spo2);
            % extract corresponding part of signal
            ind_desaturation = Spo2(start_index:stop_index);
            if (size(ind_desaturation) < 10 | isempty(ind_desaturation))
                depth(ii,:) = 0;
                DesArea(ii,:) = 0;
            else
            % flip the signal for integration
            flipped_desat = max(ind_desaturation)-ind_desaturation;
            depth(ii,:) = max(flipped_desat);
            % INTEGRATE
            DesArea(ii,:) = trapz(t_spo2(start_index:stop_index),flipped_desat);
            end
        end
        DesSev = (sum(DesArea)/TST_s);
        ave_depth = mean(depth);
        med_depth = median(depth);
    end

    %% compute total sleep event times from hypnogram and remaining events
    eventEpochs = getEventEpochs(size(hypnogram,1), arousal_info, respiratory_info, desat_info)';
    TSET_s = sum(hypnogram == 1 | hypnogram == 2 | hypnogram == 3 | hypnogram == 4 | eventEpochs == 1)*30;
    TSET_min = TSET_s / 60;
    TSET_h = TSET_min / 60;

    %% Computation of Obstruction Severity (ObsSev):
    % go trough events one at the time and link respiratory event to desaturation event.
    %First column in Eventlinks matrix is index of resp event, second is 0=hypopnea/1=apnea,
    % and third index of corresponding desat if threshold criteria is fulfilled
    Eventlinks = zeros(size(respiratory_info,1),3);
    if (isempty(Eventlinks))
        ObsSev = 0;
    else
        % set desaturation delay threshold
        desat_tr = 60;
        
        for jjj = 1:size(respiratory_info,1) % iterate through all respiratory events
            resp_event_start = respiratory_info(jjj,1);
            resp_event_duration = respiratory_info(jjj,2);
            resp_event_stop = resp_event_start + resp_event_duration;
            
            if (jjj+1 > size(respiratory_info,1)) % if the very last last respiratory event
                if (respiratory_info(jjj,3) ~= 0)
                    Eventlinks(jjj,1) = jjj;
                    Eventlinks(jjj,2) = 1;
                    
                    for kk = 1:size(desat_info)
                        if (resp_event_start < desat_info(kk,1) &...
                                resp_event_stop+desat_tr >= desat_info(kk,1))
                            % save the index of desat event // apnea
                            Eventlinks(jjj,3) = kk;
                        end
                    end
                else
                    Eventlinks(jjj,1) = jjj;
                    Eventlinks(jjj,2) = 0;
                    
                    for kk = 1:size(desat_info)
                        if (resp_event_start < desat_info(kk,1) &...
                                resp_event_stop+desat_tr >= desat_info(kk,1))
                            % save the index of desat event // hypopnea
                            Eventlinks(jjj,3) = kk;
                        end
                    end
                end
                
            else % all other respiratory events apart from the last event
                if (respiratory_info(jjj,3) ~= 0) % if event is an apnea
                    Eventlinks(jjj,1) = jjj;
                    Eventlinks(jjj,2) = 1;
                    
                    for kk = 1:size(desat_info) % iterate through each desaturation
                        if (resp_event_start < desat_info(kk,1) & ...
                                resp_event_stop+desat_tr >= desat_info(kk,1)...
                                & respiratory_info(jjj+1,1) > desat_info(kk,1))
                            % save the index of desat event // apnea
                            Eventlinks(jjj,3) = kk;
                        end
                    end
                else % if event is a hypopnea
                    Eventlinks(jjj,1) = jjj;
                    Eventlinks(jjj,2) = 0;
                    for kk = 1:size(desat_info)
                        if (resp_event_start < desat_info(kk,1) &...
                                resp_event_stop+desat_tr >= desat_info(kk,1)...
                                & respiratory_info(jjj+1,1) > desat_info(kk,1)) 
                                % the last check is to prevent a single desaturation 
                                % being linked to multiple respiratory events.
                            % save the index of desat event // hypopnea
                            Eventlinks(jjj,3) = kk;
                        end
                    end
                end
            end
        end
        
        for rho = 1:size(Eventlinks,1)
            if (Eventlinks(rho,3) == 0)
                ObsSev_nominator(rho,:) = 0;
            else
                ObsSev_nominator(rho,:) = respiratory_info(Eventlinks(rho,1),2)*...
                    DesArea(Eventlinks(rho,3));
            end
        end
        ObsSev = sum(ObsSev_nominator)/TST_s;
    end
end

end

% erase events according to given criteria
function output = removeEventsBeforeSleepOnset(hypnogram, events, eventType, exclusionCriteria)
    EXCLUDE_ALL_WAKE_DESATS = 1;
    EXCLUDE_ALL_WAKE_EVENTS = 2;
    EXCLUDE_ALL_EVENTS_BEFORE_SLEEP_ONSET = 3;
    EXCLUDE_ALL_WAKE_EVENTS_EXCEPT_SLEEP_WAKE_TRANSITION = 4;
    EVENT_TYPE_DESATS = 1;
    EVENT_TYPE_RESPIRATORY = 2;
    EVENT_TYPE_AROUSAL = 3;
    if ((exclusionCriteria == EXCLUDE_ALL_WAKE_DESATS) && (eventType ~= EVENT_TYPE_DESATS))
        output = events;
        return;
    end
    excludedEvents = events(:,1) * 0;
    for esa = 1:size(events,1)
        if ((exclusionCriteria == EXCLUDE_ALL_WAKE_DESATS) || (exclusionCriteria == EXCLUDE_ALL_WAKE_EVENTS))
            wake_epochs = find(hypnogram == 0);
            excludedEpochs = wake_epochs;
        elseif (exclusionCriteria == EXCLUDE_ALL_EVENTS_BEFORE_SLEEP_ONSET)
            sleep_epochs = find(hypnogram ~= 0);
            first_sleep_epoch = sleep_epochs(1);
            if (first_sleep_epoch > 1) 
                excludedEpochs = 1 : (first_sleep_epoch - 1);
            else
                output = events;
                return;
            end
        elseif ((exclusionCriteria == EXCLUDE_ALL_WAKE_EVENTS_EXCEPT_SLEEP_WAKE_TRANSITION))
            hypnogramDiffShiftRight = [0; diff(hypnogram)];
            hypnogramDiffShiftLeft = [diff(hypnogram); 0];
            % keep epochs where transition from sleep to wake OR wake to
            % sleep
            sleepWakeTransitions = (hypnogramDiffShiftRight < 0) & (hypnogram == 0);
            wakeSleepTransitions = (hypnogramDiffShiftLeft > 0) & (hypnogram == 0);
            sleep = hypnogram ~= 0;
            excludedEpochs = find(~(sleepWakeTransitions | wakeSleepTransitions | sleep));
        end
        excludedEpochStartTimes = excludedEpochs*30 - 30;
        if (isempty(find( ...
                (events(esa,1) >= excludedEpochStartTimes) &...
                (events(esa,1) < (excludedEpochStartTimes+30))... 
                )))
            excludedEvents(esa) = 0;
        else
            excludedEvents(esa) = find( ...
                (events(esa,1) >= excludedEpochStartTimes) &...
                (events(esa,1) < (excludedEpochStartTimes+30))... 
                );
        end
    end
    includedEvents = find(excludedEvents == 0);
    output = events(includedEvents,:);
end

%% find which epochs contain at least 1 scored event
function eventEpochs = getEventEpochs(epochCount, arousal_info, respiratory_info, desat_info)
    epochStartTimes = (0:(epochCount-1))*30;
    eventEpochs = epochStartTimes * 0;
    for epoch = 1:epochCount
        startTime = epochStartTimes(epoch);
        eventEpochs(epoch) = 0;
        if (eventsInEpoch(arousal_info, startTime) || eventsInEpoch(respiratory_info, startTime) || eventsInEpoch(desat_info, startTime))
            eventEpochs(epoch) = 1;
        end
    end
    % figures for debugging purposes to confirm scored events line up with eventEpochs
    %figure; stairs(eventEpochs,'r'); xlim([1,size(eventEpochs,2)]); ylim([-0.1 1.1]);
    %figure; hold on; stairs(arousal_info(:,1),'r'); stairs(desat_info(:,1),'g'); stairs(respiratory_info(:,1),'b'); ylim([0, epochStartTimes(end) + 30]); view([90 -90]);
end

%% find whether there are any events in a particular epoch
function within = eventsInEpoch(event_info, epochStartTime)
    within = 0;
    checkAllEvents = (event_info(:,1) >= epochStartTime) & (event_info(:,1) < (epochStartTime + 30));
    if (sum(checkAllEvents) > 0)
        within = 1;
    end
end