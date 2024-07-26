function [time_part, Tstart, Tend] = FindTimeForVolume(AUCdata, FlowData, ThresAUC, dt)
%% FindTimeForVolume Determine shortest time to reach volume threshold
% Finds the shortest time, with the maximum average flow, that achieves
% the required volume threshold
%
%Inputs:
% Volume data
% Flow data
% Threshold value
% sample resolution
%
%Outputs:
% real time duration (seconds)
% Start time (index)
% End time (index)


CsumAUC = AUCdata;
% preallocate space for candidates
candidates = zeros(length(AUCdata),4);
for m = 1:(length(CsumAUC))
    % find first value in cumsum to >= threshold
    n = find(CsumAUC>=ThresAUC,1);
    if ~isnan(n)
        % find the length of xdata over this period
        time = (n-m)*dt; %index converted to time
        avgsignal = mean(FlowData(m:n));
        % and store these values
        candidates(m,:) = [m; n; time; avgsignal];
        % then step CsumAUC forward through breath
        CsumAUC = CsumAUC-CsumAUC(m);
    end
end

% remove values with time of zero
% probably not really required, but let's be sure
candidates = candidates(candidates(:,3)~=0,:);

% remove entries where end index < start index
% same as above, shouldn't really be required, but let's be sure
candidates = candidates(candidates(:,2)>candidates(:,1),:);

% then find the biggest mean signal value from the smallest time values
[signal_mean] = max(candidates(candidates(:,3) == min(candidates(:,3)),4));

% and the smallest time value from the biggest mean signal value
[time_part] = min(candidates(candidates(:,4) == max(candidates(:,4)),3));

% catch potential errors
if time_part == dt % time part cannot equal one sample
    time_part = [];
    Tstart = [];
    Tend = [];
else
    % now, these two values intersect and refer to just one point
    % that being the BestFit point
    BestFit = candidates(candidates(:,3) == time_part & candidates(:,4) == signal_mean,:);
    %suppose we find this, and it is unique, that's good
    if size(BestFit) == [1, 4]
        Tstart = BestFit(1);
        Tend = BestFit(2);
    else
        % but, what if  (a) these two values don't intersect, or
        %               (b) produce multiple results
        % pick the one(s) with signal_mean before time_part
        % (looking at the data, this appears to be the best option)
        BestFit = candidates(candidates(:,4) == signal_mean,:);
        if size(BestFit,1) >= 2
            % if we still have multiple rows, just pick the first one
            BestFit = BestFit(1,:);
        end 
        if ~isempty(BestFit)
            Tstart = BestFit(1);
            Tend = BestFit(2);
        else % if all else fails, just return empty array... safety first.
            Tstart = [];
            Tend = [];
        end
    end 
end
end
