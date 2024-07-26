function [StableBB_mins, StableBB_table] = getStableBreathingPeriods(SigT,ChannelsList)
% change to SigT
% Finds periods of stable breathing, and returns the duration associated 
% with each stable breathing period.
%
% Inputs: 
%   DataEventHypnogMat (only needs: Time, EventsAr, EventsResp and Epochs)
%   ChannelsList
%
% Outputs:
%   StableBB_mins, a column vector (of length Time) showing duration of stable breathing in minutes 
%   StableBB_table, a table showing summary data for each stable breathing period, as
%       index (start and end), time (start and end), duration (seconds), duration (minutes)
%
%
% First, what is stable breathing?
%   Continuous periods of sleep, without:
%       - Arousals (Arousal breaths) (EventsAr~=0)
%       - Respiratory Events (Hypopnea, Apnea) , mixed central obs, (EventsResp~=0)
%   but what about:
%       - desaturation events
%       - limb movements/PLM

EvtsAr = SigT.EventsAr; % nnz(EvtsAr)
EvtsResp = SigT.EventsResp; % nnz(EvtsResp)
EvtsAll = EvtsAr|EvtsResp; % nnz(EvtsAll)
Hypnog = SigT.Epochs;
Sleep = Hypnog<4 & Hypnog >=0; % nnz(Sleep)
% sleep without events is 'sleep' and 'not event'
SleepWOevents = Sleep&~EvtsAll;  % nnz(SleepWOevents)

if 0 % fig to visualise procedure
    figure(1001); clf(figure(1001)); fig = gcf;
    fig.Color = [ 1 1 1 ]; fig.Units = 'Inches';
    fig.Position = [27 1.5 13 9.5];
    ax1001(1) = subplot(2,1,1)
    plot(Sleep); hold on;
    plot(EvtsAll*0.9);
    ylim([-0.1 1.1]);
    ylabel('Sleep=1, Evts=0.9');
    ax1001(2) = subplot(2,1,2)
    plot(SleepWOevents); hold on;
    ylim([-0.1 1.1]);
    ylabel('Stable sleep');
    linkaxes(ax1001, 'x');
end

%                 % remove periods with noise
%                 warning('off'); % turn off, absent signal will cause warning
%                 respwav=DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'Flow')==1));
%                 if ~exist('Time', 'var') % potentially have Time from before, but check and load incase
%                     Time=DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'Time')==1));
%                 end
%                 noisewav = FlowSignalToNoise(Time,respwav,settings.plotfigure);
%                 warning('on');
%                 Noise = noisewav==3;

StableBB_logical = SleepWOevents;% &~Noise; % nnz(StableBB)

% tidy up all the intermediate variables
clear EvtsAr EvtsResp EvtsAll Hypnog Sleep SleepWOevents noisewav Noise

% find continuous blocks
start_ind_temp=find(diff([0 single(StableBB_logical)'])==1);
end_ind_temp=find(diff([0 single(StableBB_logical)'])==-1);

% force into vertical array
start_ind_temp = start_ind_temp(:);
end_ind_temp = end_ind_temp(:);

% error check here to cover condition if end of event is end of study etc.
[start_ind, end_ind] = TidyStartEndEventList(start_ind_temp, end_ind_temp, length(StableBB_logical));

if 0 % fig to visualise procedure
    figure(1001); clf(figure(1001)); fig = gcf;
    fig.Color = [ 1 1 1 ]; fig.Units = 'Inches';
    fig.Position = [27 1.5 13 3];
    plot(StableBB_logical); hold on;
    plot(start_ind(1), 0, 'rx');
    plot(start_ind(2), 0, 'rx');
    plot(end_ind(1), 0, 'ro');
    plot(end_ind(2), 0, 'ro');
    plot(start_ind(end-1), 0, 'rx');
    plot(start_ind(end), 0, 'rx');
    plot(end_ind(end-1), 0, 'ro');
    plot(end_ind(end), 0, 'ro');
    ylim([-0.1 1.1]);
    ylabel('Stable sleep');
end

% find the dT
Time=SigT.Time;
dT = Time(2)-Time(1);

% calculate time of continuous blocks
StableBB_indx = [start_ind, end_ind];
StableBB_time = [Time(start_ind), Time(end_ind)];
StableBB_duration = (StableBB_indx(:,2)-StableBB_indx(:,1)) .* dT;

% remove any blocks below time threshold
minDuration = 60; % 1 minute
TooShort = find(StableBB_duration<minDuration);
StableBB_indx(TooShort,:) = [];
StableBB_time(TooShort,:) = [];
StableBB_duration(TooShort) = [];

% label continuous blocks with duration in minutes (floor)
StableBB_durationMin = floor(StableBB_duration./60);
StableBB_durationMin(StableBB_durationMin>10) = 10; % cap the max duration at 10 minutes
StableBB_mins = zeros(height(SigT),1);
for i = 1: size(StableBB_indx,1)
    StableBB_mins(StableBB_indx(i,1):StableBB_indx(i,2)) = StableBB_durationMin(i);
end

if 0 % fig to visualise procedure
    figure(1001); clf(figure(1001)); fig = gcf;
    fig.Color = [ 1 1 1 ]; fig.Units = 'Inches';
    fig.Position = [27 1.5 13 3];
    plot(StableBB_logical); hold on;
    plot(StableBB_mins);
    ylim([-1 max(StableBB_mins)+1]);
end

% keep the data for continuous blocks in a handy table format
StableBB_table = table(StableBB_indx, StableBB_time, StableBB_duration, StableBB_durationMin);

% tidy up
%clear start_ind end_ind TimeChan StableBB_mins StableBB_logical StableBBblocks StableBBblocksDuration StableBBblocksDurationMin
