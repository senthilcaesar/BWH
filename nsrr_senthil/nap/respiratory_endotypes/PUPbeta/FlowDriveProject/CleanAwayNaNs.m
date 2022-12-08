function PtData = CleanAwayNaNs(PtData)
% just a helper function to clean NaN's while testing features

%% remove breaths that are NaN for ALL features
BBwithALLNans = find(all(isnan(PtData{:,:}),2));
str = ['Removing ', num2str(length(BBwithALLNans)), ' breaths with ALL NaN features']; disp(str);
PtData(BBwithALLNans,:)=[];

%% remove features with high NaN counts (using ~isfinite instead of isnan)
% in current usage, we set this very high because we want to keep all cols
SetThreshold = 0.25;% 0.001; % those with more than a set threshold of NaN
Threshold = round(size(PtData,1)*SetThreshold);
%NaNSpace = isnan(PtData{:,:});
NaNSpace = ~isfinite(PtData{:,:});
ftrsToExclude = find(sum(NaNSpace(:,:))>Threshold);
str = ['Removing ', num2str(length(ftrsToExclude)), ' features with more than ', num2str(SetThreshold*100), '% of NaN data']; disp(str);
PtData(:,ftrsToExclude)=[];
%   %    ~ Ftrs removed
%   3       4
%   1       6
%   0.5     8
%   0.1     10

%% remove any remaining breaths that are NaN for ANY features
BBwithNans = find(any(isnan(PtData{:,:}),2));
str = ['Removing ', num2str(length(BBwithNans)), ' breaths with ANY NaN features']; disp(str);
PtData(BBwithNans,:)=[];

end