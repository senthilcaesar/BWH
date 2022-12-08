function [BreathSnoreDataTable,PSDbreath,EnvelpBreath,EnvelpBreath_n] = ComputeSnoreFeaturesOld(SigT,ZCrTransform,SnoreInterpStruct,SnoreStruct,BBtime,DBthresh)

% SigT replaces DataEventHypnog_Mat to run in AnalysisOne.m-RMA-5/3/2021
% please make modifications in any main/sub functions as needed

global ChannelsList

%% initialize variables
Time = SigT.Time;
TimeSm = SnoreStruct.SnoreTimeSm;
dt=(Time(end)-Time(1))/(length(Time)-1);
Fs = round(1/dt);

Time_i_start = BBtime(:,1);
Time_i_end = BBtime(:,2);
Time_e_start = BBtime(:,2);
Time_e_end = BBtime(:,3);  

BBs=length(Time_i_start);

BreathSnoreStruct = struct();

% figure, plot(SnoreInterpStruct.Spectrogram.Time, nanmean(SnoreInterpStruct.Spectrogram.lpcenv_n,2),'k',...
%     SnoreStruct.Spectrogram.Time, nanmean(SnoreStruct.Spectrogram.lpcenv_n,2),'r--')

SnoreInterpStruct.SnoreDB = SigT.SnoreDB; % add snoreDB to structure

harmonicIdx = SnoreInterpStruct.HNR > 0.5 & SnoreInterpStruct.HarmonicPower > 0.5; % to be analyzed

% interpolate PSD envelope and power
samplepoints = 1:length(SnoreStruct.Spectrogram.Time);
querypoints = linspace(1, length(SnoreStruct.Spectrogram.Time), length(SnoreInterpStruct.SnoreTime));
SnoreInterpStruct.Spectrogram.Power = interp1(samplepoints', SnoreStruct.Spectrogram.Power, querypoints', 'previous');
SnoreInterpStruct.Spectrogram.Time = interp1(samplepoints', SnoreStruct.Spectrogram.Time, querypoints', 'previous');
SnoreInterpStruct.Spectrogram.lpcenv_n = interp1(samplepoints', SnoreStruct.Spectrogram.lpcenv_n, querypoints', 'previous');
SnoreInterpStruct.Spectrogram.lpcenv = interp1(samplepoints', SnoreStruct.Spectrogram.lpcenv, querypoints', 'previous');

% figure, plot(SnoreInterpStruct.Spectrogram.Time, nanmean(SnoreInterpStruct.Spectrogram.lpcenv_n,2),'k',...
%     SnoreStruct.Spectrogram.Time, nanmean(SnoreStruct.Spectrogram.lpcenv_n,2),'r--')
% 
% BreathDataTable = evalin('base','BreathDataTable');
% BreathSnoreTable = evalin('base','BreathSnoreTable');
%% identify and sample variables to be analyzed
PSDbreath = nan(length(Time_i_start),length(SnoreStruct.Spectrogram.Freq));
EnvelpBreath = nan(length(Time_i_start),length(SnoreStruct.Spectrogram.Freq));
EnvelpBreath_n = nan(length(Time_i_start),length(SnoreStruct.Spectrogram.Freq));
SnoreVarnames = fieldnames(SnoreInterpStruct); % make snore variable names to loop through (exclude Spectrogram below
numrows = size(SnoreInterpStruct.SnoreTime,1);

% BB_Times = evalin('caller','BB_Times')
for i = 1:BBs   
    % get breath indices
    idx_i = find(Time > Time_i_start(i) & Time < Time_i_end(i));
    idx_e = find(Time > Time_e_start(i) & Time < Time_e_end(i));
    
    idx_i_sm = find(TimeSm > Time_i_start(i) & TimeSm < Time_i_end(i));
    idx_e_sm = find(TimeSm > Time_e_start(i) & TimeSm < Time_e_end(i));
    
    if DBthresh > 0
        idx_i = idx_i(SnoreInterpStruct.SnoreDB(idx_i) > DBthresh);
        idx_e = SnoreInterpStruct.SnoreDB(idx_e) > DBthresh;
    end
    
%     idx_i2 = BreathDataTable.Time0(i)*125 + BreathDataTable.BB_i_start(i):...
%         BreathDataTable.Time0(i)*125 + BreathDataTable.BB_i_mid(i);
    
    
    % loop through each variable
    for j = 1:length(SnoreVarnames)
        varname = SnoreVarnames{j};
        variable = SnoreInterpStruct.(varname);
        if length(variable) > numrows || sum(strcmp({'SnoreTime','FilteredSnd', 'Spectrogram'},varname))==1
            continue % if not the same size then its just a var I smuggled into the struct
        end
        
        if strcmp(varname, 'SnoreDB')
            variable = 10.^((variable/10))*(0.00002^2);
        end

        nanadd = nan(1,size(variable,2));
        BreathSnoreStruct.([varname,'Mean_i'])(i,:) = nanmean(variable(idx_i,:),1);
        BreathSnoreStruct.([varname,'Mean_e'])(i,:) = nanmean(variable(idx_e,:),1);
        BreathSnoreStruct.([varname,'Median_i'])(i,:) = nanmedian(variable(idx_i,:),1);
        BreathSnoreStruct.([varname,'Median_e'])(i,:) = nanmedian(variable(idx_e,:),1);
        BreathSnoreStruct.([varname,'Max_i'])(i,:) = max([variable(idx_i,:);nanadd],[],1); % nan bit returns nan if empty
        BreathSnoreStruct.([varname,'Max_e'])(i,:) = max([variable(idx_e,:);nanadd],[],1);
        BreathSnoreStruct.([varname,'Var_i'])(i,:) = nanvar(variable(idx_i,:),1);
        BreathSnoreStruct.([varname,'Var_e'])(i,:) = nanvar(variable(idx_e,:),1);
        
        % Half breath metrics
        harmonicIdx_ = harmonicIdx(idx_i)';
        % half indices
        falsemat = false(size(idx_i));
        idx_i_ = find(idx_i);
        halfidx1 = falsemat;
        halfidx1(idx_i_(1:round(length(idx_i_)/2))) = true;
       
        halfidx2 = falsemat;
        halfidx2(idx_i_(round(length(idx_i_)/2)+1:length(idx_i_))) = true;

        % percent of inspiration that is harmonic
        percHarmonic = sum(harmonicIdx_)/sum(idx_i);
        if isinf(percHarmonic), percHarmonic = nan; end
        BreathSnoreStruct.percHarmonic(i,:) = percHarmonic;
        BreathSnoreStruct.percHarmonic_halfI1(i,:) = sum(harmonicIdx_(halfidx1))/sum(halfidx1);
        BreathSnoreStruct.percHarmonic_halfI2(i,:) = sum(harmonicIdx_(halfidx2))/sum(halfidx2);

        % pitch of first and second half of insp
        BreathSnoreStruct.Pitch_halfI1(i,:) = nanmean(SnoreInterpStruct.f0(halfidx1));
        BreathSnoreStruct.Pitch_halfI2(i,:) = nanmean(SnoreInterpStruct.f0(halfidx2));
    end
    % ensemble average PSD 

    PSDdata = SnoreInterpStruct.Spectrogram.Power(idx_i,:);
    PSDdata2 = 10.^((PSDdata/10))*(0.00002^2);
    PSDmean = nanmean(PSDdata2,1);
    PSDbreath(i,:) = 10*(log10(PSDmean/(0.00002^2)));

    % ensemble average envelope
    envelp = SnoreInterpStruct.Spectrogram.lpcenv_n(idx_i,:);
    envelp2 = 10.^((envelp/10))*(0.00002^2); % convert to sound pressure level
    envelpmean = nanmean(envelp2,1);
    EnvelpBreath_n(i,:)  = 10*(log10(envelpmean/(0.00002^2))); % convert back to DB

    % ensemble average envelope
    envelp = SnoreInterpStruct.Spectrogram.lpcenv(idx_i,:);
    envelp2 = 10.^((envelp/10))*(0.00002^2); % convert to sound pressure level
    envelpmean = nanmean(envelp2,1);
    EnvelpBreath(i,:)  = 10*(log10(envelpmean/(0.00002^2))); % convert back to DB
    
    % Compute multiscale entropy
    ZCrBreath = ZCrTransform(idx_i_sm);
    m=2; r=0.15;
    for tau = 1:20
        [MSEtemp,~,~] = multiscaleSampleEntropy(ZCrBreath, m, r, tau);
        BreathSnoreStruct.(['MSE_',num2str(tau)])(i,:) = MSEtemp;
    end

%     idx = SnoreStruct.SnoreTime > BreathDataTable.Time_start(i) &...
%     SnoreStruct.SnoreTime < BreathDataTable.Time_mid(i);

%     % ensemble average PSD
%     PSDdata = SnoreStruct.Spectrogram.Power(idx,:);
%     PSDbreath2(i,:) = nanmean(PSDdata,1);
% 
%     % ensemble average envelope
%     envelp = SnoreStruct.Spectrogram.lpcenv_n(idx,:);
%     EnvelpBreath2(i,:) = nanmean(envelp,1);
% 
%     figure(12), plot(SnoreStruct.Spectrogram.Freq, EnvelpBreath2(i,:), 'k',...
%         SnoreStruct.Spectrogram.Freq,EnvelpBreath(i,:), 'r--');
end

%% Convert SnoreDB back to DB
fldnms = fieldnames(BreathSnoreStruct);
breathswithsnoreDB = find(contains(fldnms, 'SnoreDB'))';
for nn = breathswithsnoreDB
    variable = fldnms{nn};
    BreathSnoreStruct.(variable) = 10*(log10(BreathSnoreStruct.(variable)/(0.00002^2)));
end

%% combine into one big table
%FT = FeatureNames';
BreathSnoreDataTable = struct2table(BreathSnoreStruct);
% diffsnoretable = (BreathSnoreDataTable{:,1:193} - BreathSnoreTable{:,:})./BreathSnoreTable{:,:};