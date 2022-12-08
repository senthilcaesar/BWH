% Make variables from BreathDataTable compatible with XHz
% Manually load XHz and EventAnalyzed data

vars = {'VI','Vdr_est'};
XHzStruct = struct();
for mm = 1:length(vars)
    varXHz = interp1(BreathDataTableFulls.Time_end, BreathDataTableFulls.(vars{mm}), DataEventHypnog_Mat(:,1),'previous');
    DataEventHypnog_Mat(:,end+1) = varXHz;
end

ChannelsList = [ChannelsList vars];