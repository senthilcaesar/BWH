function [StTime, BB_Times] = getBBTimes(BreathDataTable, n, win)
%% Get times for breaths in current window
    StTime = table2array(BreathDataTable{n}{win}(1,1));
    BB_start = table2array(BreathDataTable{n}{win}(:,2)) - StTime;
    BB_mid = table2array(BreathDataTable{n}{win}(:,3)) - StTime;
    BB_end = table2array(BreathDataTable{n}{win}(:,4)) - StTime;
    BB_Times = [BB_start, BB_mid, BB_end];
end
%% BreathDataTable contains the following data:
%
% 1 'Time0'
% 2 'Time_start'
% 3 'Time_mid'
% 4 'Time_end'
% 5 'BB_i_start'
% 6 'BB_i_mid'
% 7 'BB_i_end'
% 8 'VI'
% 9 'AR'
% 10 'spo2'
% 11 'pos_B'
% 12 'hypnog_B'
% 13 'Etype'
% 14 'DeltaPes'
% 15 'DeltaPmus'
% 16 'VIpes'
% 17 'DeltaEdi'
% 18 'VIedi'
% 19 'GGpeak'
% 20 'GGtonic'
% 21 'FlowPes_VI'
% 22 'FlowEdi_VI'
