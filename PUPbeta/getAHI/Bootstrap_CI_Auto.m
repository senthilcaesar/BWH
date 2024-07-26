function [AHI, BS_CI_low,BS_CI_high]=Bootstrap_CI_Auto(Evts,TST,alpha_level,iterations)
%% Bootstrap Confidence Interval_Auto
% Code to compute AHI confidence interval using bootstrap
% Adjusted to Raichel's Auto score program
% Dataset: MESA
% Input:
%     Evts: Apnea event file from Raichel's Auto score program - required
%     TST: Total sleep time in hr - required (default: 8)
%     alpha_level: alpha level of confidence interval in % - required (default:5)
%     iterations: bootstrap iteration times - required (default: 3000) 
%
% Output:
%      AHI: Apnea-Hypopnea Index (events/hr)
%      BS_CI_low: bootstrap lower confidence bounds
%      BS_CI_high:bootstrap higher confidence bounds
%
% Corresponding to sleep medicine 2020, Robert J. Thomas, Shuqiang Chen, Uri T. Eden, Michael J. Prerau*
% "Quantifying statistical uncertainty in metrics of sleep disordered breathing"
%
% Shuqiang Chen, Last modified 01/06/2021
%% ********************************************************************
% Set default values
if nargin<2
    TST = 8; 
end

if nargin<3
    alpha_level = 5; % Meaning 95% confidence bounds [2.5% 97.5%]
end

if nargin<4
    iterations = 3000;
end

% Evts file;
Evts_idx = Evts.RespT.InclAHI3a & Evts.RespT.Epochs<4;
Evts_time = Evts.RespT.EventStart(Evts_idx==1);
% double check with  Evts.AHIdata2{1,1}.AllSleepAllPahi(2)
% Evts_idxtemp = Evts.RespT.EventStart(Evts_idx==1& Evts.RespT.Epochs<4);


% Compute IEIs(Substract wake times)
N = length(Evts_time);
AHI = N/(Evts.ArTinfo.TSTmin/60);
IEIs= diff(Evts_time);
wake_time=zeros(length(IEIs),1);
wake_idx=Evts.Hypnogram == 4;

for ii=1:N-1
    estart=Evts.RespT{ii,2};
    eend=Evts.RespT{ii+1,2};
    
    wake_bouts=find(Evts.Hypnogram_t>estart & Evts.Hypnogram_t<eend & wake_idx);
    wake_time(ii)= length(wake_bouts)*30;
end

BS_times=IEIs-wake_time;
BS_times=BS_times/3600;

% Compute the number of iterations guaranteed to go over TST
N_iter= ceil(TST/prctile(BS_times,5))+1;

% Get the boostrap samples
BS_data=BS_times(randi(length(BS_times),iterations,N_iter));
BS_vals=cumsum(BS_data,2)<=(TST);

% Find the number of counts
BS_counts=sum(BS_vals,2);

ptiles=prctile(BS_counts/(TST),[alpha_level/2 100-alpha_level/2]);
BS_CI_low=ptiles(1);
BS_CI_high=ptiles(2);

%     figure
%     histogram(BS_counts/(TST),100);
%     hold on;
%     vline(35,2)
%     
%     vline(BS_CI_low,2,'b');
%     vline(BS_CI_high,2,'b');

end