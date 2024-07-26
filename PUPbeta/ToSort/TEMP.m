%% Remove aberrant breaths
% Remove breaths before/after 2 and 98th percentile
Time_BB = BreathDataTableBig.Time_end - BreathDataTableBig.Time_start;
Time_Lower = prctile(Time_BB, 2);
Time_Upper = prctile(Time_BB, 98);

BreathDataTableHyp.Time_BB = BreathDataTableHyp.Time_end - ...
                             BreathDataTableHyp.Time_start; 
                         
BreathDataTableFL.Time_BB = BreathDataTableFL.Time_end - ...
                             BreathDataTableFL.Time_start; 
                         
rmIdxHyp = BreathDataTableHyp.Time_BB > Time_Upper;
shortBreath = BreathArrayHyp(rmIdxHyp,:);
figure
for ii = 1:size(shortBreath,1)
    plot(shortBreath(ii,:))
end