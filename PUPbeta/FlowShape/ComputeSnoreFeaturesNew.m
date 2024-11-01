function [BreathSnoreDataTable,PSDbreath,EnvelpBreath,EnvelpBreath_n,PwelchSmooth] = ComputeSnoreFeaturesNew(SnoreStruct,BBtime,DBthresh)

global ChannelsList

%% initialize variables
% Time = DataEventHypnog_Mat(:,strcmp(ChannelsList,'Time'));
% Time = SnoreStruct.SnoreTime;
% TimeSm = ZCrTransform.SnoreTimeSm;
% ZCrTransform = ZCrTransform.ZCr;
% dt=(Time(end)-Time(1))/(length(Time)-1);
% Fs = round(1/dt);

Time_i_start = BBtime(:,1);
Time_i_end = BBtime(:,2);
Time_e_start = BBtime(:,2);
Time_e_end = BBtime(:,3);  

BBs=length(Time_i_start);

BreathSnoreStruct = struct();

% figure, plot(SnoreInterpStruct.Spectrogram.Time, nanmean(SnoreInterpStruct.Spectrogram.lpcenv_n,2),'k',...
%     Spectrogram.Time, nanmean(Spectrogram.lpcenv_n,2),'r--')

% SnoreInterpStruct.SnoreDB = DataEventHypnog_Mat(:,strcmp(ChannelsList,'SnoreDB')); % add snoreDB to structure

% interpolate PSD envelope and power
% samplepoints = 1:length(Spectrogram.Time);
% querypoints = linspace(1, length(Spectrogram.Time), length(SnoreInterpStruct.SnoreTime));
% SnoreInterpStruct.Spectrogram.Power = interp1(samplepoints', Spectrogram.Power, querypoints', 'previous');
% SnoreInterpStruct.Spectrogram.Time = interp1(samplepoints', Spectrogram.Time, querypoints', 'previous');
% SnoreInterpStruct.Spectrogram.lpcenv_n = interp1(samplepoints', Spectrogram.lpcenv_n, querypoints', 'previous');
% SnoreInterpStruct.Spectrogram.lpcenv = interp1(samplepoints', Spectrogram.lpcenv, querypoints', 'previous');

% figure, plot(SnoreInterpStruct.Spectrogram.Time, nanmean(SnoreInterpStruct.Spectrogram.lpcenv_n,2),'k',...
%     Spectrogram.Time, nanmean(Spectrogram.lpcenv_n,2),'r--')
% 
% BreathDataTable = evalin('base','BreathDataTable');
% BreathSnoreTable = evalin('base','BreathSnoreTable');
%% identify and sample variables to be analyzed
PSDbreath = nan(length(Time_i_start),length(SnoreStruct.Freq));
EnvelpBreath = nan(length(Time_i_start),length(SnoreStruct.lpcFreq));
EnvelpBreath_n = nan(length(Time_i_start),length(SnoreStruct.lpcFreq));
PwelchSmooth = nan(length(Time_i_start),size(SnoreStruct.lpcCellArray{3,1},2));
SnoreVarnames = SnoreStruct.SnoreCellArray{1,1}.Properties.VariableNames; % make snore variable names to loop through (exclude Spectrogram below
% numrows = size(SnoreStruct.SnoreTime,1);

%% Generate window-level data from which to compute breath-level metrics
% This will be updated in the loop to keep array size down and save on RAM
k=2;
if size(SnoreStruct.SnoreCellArray,2) == 1
    k=1;
    SnoreTable = SnoreStruct.SnoreCellArray{1,k};
    SnoreTableSm =SnoreStruct.SnoreCellArray{2,k};
    SpectP = SnoreStruct.SpectrogramCellArray{1,k};
    lpcenv = SnoreStruct.lpcCellArray{1,k};
    lpcenv_n = SnoreStruct.lpcCellArray{2,k};
    pwsmooth = SnoreStruct.lpcCellArray{3,k};
    Time = SnoreStruct.SnoreCellArray{1,k}.timeWin;
else 
    SnoreTable = vertcat(SnoreStruct.SnoreCellArray{1,k-1},SnoreStruct.SnoreCellArray{1,k});
    SnoreTableSm = vertcat(SnoreStruct.SnoreCellArray{2,k-1},SnoreStruct.SnoreCellArray{2,k});
    SpectP = vertcat(SnoreStruct.SpectrogramCellArray{1,k-1},SnoreStruct.SpectrogramCellArray{1,k});
    lpcenv = vertcat(SnoreStruct.lpcCellArray{1,k-1},SnoreStruct.lpcCellArray{1,k});
    lpcenv_n = vertcat(SnoreStruct.lpcCellArray{2,k-1},SnoreStruct.lpcCellArray{2,k});
    pwsmooth = vertcat(SnoreStruct.lpcCellArray{3,k-1},SnoreStruct.lpcCellArray{3,k});
    Time = SnoreStruct.SnoreCellArray{1,k-1}.timeWin;
end


for i = 1:BBs
    % check if window-level arrays need updating
    % this routine allows me to only analyze two concatenated tables at a
    % time which (I think) saves my RAM. Otherwise I would have a massive
    % snore win table with as many as five hundred thousand of rows
    if Time_e_end(i) > Time(end) && size(SnoreStruct.SnoreCellArray,2) > k
        k = k+1;
        SnoreTable = vertcat(SnoreStruct.SnoreCellArray{1,k-1},SnoreStruct.SnoreCellArray{1,k});
        SnoreTableSm = vertcat(SnoreStruct.SnoreCellArray{2,k-1},SnoreStruct.SnoreCellArray{2,k});
        SpectP = vertcat(SnoreStruct.SpectrogramCellArray{1,k-1},SnoreStruct.SpectrogramCellArray{1,k});
        lpcenv = vertcat(SnoreStruct.lpcCellArray{1,k-1},SnoreStruct.lpcCellArray{1,k});
        lpcenv_n = vertcat(SnoreStruct.lpcCellArray{2,k-1},SnoreStruct.lpcCellArray{2,k});
        pwsmooth = vertcat(SnoreStruct.lpcCellArray{3,k-1},SnoreStruct.lpcCellArray{3,k});
    end
    Time = SnoreTable.timeWin;
    TimeSm = SnoreTableSm.timeWinSm;
    ZCrTransform = SnoreTableSm.ZCrTransform;
    harmonicIdx = SnoreTable.HNR > 0.5 & SnoreTable.HarmonicPower > 0.5; % to be analyzed
    
    % get breath indices
    if 1
        idx_i = find(Time > Time_i_start(i) & Time < Time_i_end(i));
        idx_e = find(Time > Time_e_start(i) & Time < Time_e_end(i));

        idx_i_sm = find(TimeSm > Time_i_start(i) & TimeSm < Time_i_end(i));
        idx_e_sm = find(TimeSm > Time_e_start(i) & TimeSm < Time_e_end(i));
    else
        % get breath indices (first half of breath only)
        idx_i = find(Time > Time_i_start(i) & Time < Time_i_start(i)+(Time_i_end(i)-Time_i_start(i))/2);
        idx_e = find(Time > Time_e_start(i) & Time < Time_e_start(i)+(Time_e_end(i)-Time_e_start(i))/2);

        idx_i_sm = find(TimeSm > Time_i_start(i) & TimeSm < Time_i_start(i)+(Time_i_end(i)-Time_i_start(i))/2);
        idx_e_sm = find(TimeSm > Time_e_start(i) & TimeSm < Time_e_start(i)+(Time_e_end(i)-Time_e_start(i))/2);
    end
    if DBthresh > 0
        idx_i = idx_i(SnoreTable.SnoreDB(idx_i) > DBthresh);
%         idx_e = idx_e(SnoreTable.SnoreDB(idx_e) > DBthresh);
    end
%     if DBthresh > 0 % had to change this because DB not reliable with Nox. Will have to think about this
%         if length(idx_i) > 7
%             idx_i = idx_i(3:end-2);
%         end
%         
%         if length(idx_e) > 7
%             idx_e = idx_e(3:end-2);
%         end
%     end
    
%     idx_i2 = BreathDataTable.Time0(i)*125 + BreathDataTable.BB_i_start(i):...
%         BreathDataTable.Time0(i)*125 + BreathDataTable.BB_i_mid(i);
    
    
    %% loop through each variable and compute breath-level features
    for j = 1:length(SnoreVarnames)
        varname = SnoreVarnames{j};
        variable = SnoreTable.(varname);
        if sum(strcmp({'timeWin'},varname))==1
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
        
        if length(idx_i) == 1
            BreathSnoreStruct.([varname,'Var_i'])(i,:) = nanadd;
        else
            BreathSnoreStruct.([varname,'Var_i'])(i,:) = nanvar(variable(idx_i,:),1);            
        end
        
        if length(idx_i) == 1
            BreathSnoreStruct.([varname,'Var_e'])(i,:) = nanadd;
        else
            BreathSnoreStruct.([varname,'Var_e'])(i,:) = nanvar(variable(idx_e,:),1);  
        end
    end
    
    %% Half breath metrics
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
    BreathSnoreStruct.Pitch_halfI1(i,:) = nanmean(SnoreTable.f0(halfidx1));
    BreathSnoreStruct.Pitch_halfI2(i,:) = nanmean(SnoreTable.f0(halfidx2));
    % ensemble average PSD 
    
    %% Spectral ensembles
    PSDdata = SpectP(idx_i,:);
    PSDdata2 = 10.^((PSDdata/10))*(0.00002^2);
    PSDmean = nanmean(PSDdata2,1);
    PSDbreath(i,:) = 10*(log10(PSDmean/(0.00002^2)));

    % ensemble average envelope
    envelp = lpcenv_n(idx_i,:);
    envelp2 = 10.^((envelp/10))*(0.00002^2); % convert to sound pressure level
    envelpmean = nanmean(envelp2,1);
    EnvelpBreath_n(i,:)  = 10*(log10(envelpmean/(0.00002^2))); % convert back to DB

    % ensemble average envelope
    envelp = lpcenv(idx_i,:);
    envelp2 = 10.^((envelp/10))*(0.00002^2); % convert to sound pressure level
    envelpmean = nanmean(envelp2,1);
    EnvelpBreath(i,:)  = 10*(log10(envelpmean/(0.00002^2))); % convert back to DB
    
    % ensemble average envelope
    pwsmoothtemp = pwsmooth(idx_i,:);
    pwsmooth2 = 10.^((pwsmoothtemp/10))*(0.00002^2); % convert to sound pressure level
    pwsmoothmean = nanmean(pwsmooth2,1);
    PwelchSmooth(i,:)  = 10*(log10(pwsmoothmean/(0.00002^2))); % convert back to DB
    
    % Compute multiscale entropy
    ZCrBreath = ZCrTransform(idx_i_sm);
    m=2; r=0.15;
    for tau = 1:20
        try
            [MSEtemp,~,~] = multiscaleSampleEntropy(ZCrBreath, m, r, tau);
            BreathSnoreStruct.(['MSE_',num2str(tau)])(i,:) = MSEtemp;
        catch
            BreathSnoreStruct.(['MSE_',num2str(tau)])(i,:) = nan;
        end
    end

%     idx = SnoreStruct.SnoreTime > BreathDataTable.Time_start(i) &...
%     SnoreStruct.SnoreTime < BreathDataTable.Time_mid(i);

%     % ensemble average PSD
%     PSDdata = Spectrogram.Power(idx,:);
%     PSDbreath2(i,:) = nanmean(PSDdata,1);
% 
%     % ensemble average envelope
%     envelp = Spectrogram.lpcenv_n(idx,:);
%     EnvelpBreath2(i,:) = nanmean(envelp,1);
% 
%     figure(12), plot(Spectrogram.Freq, EnvelpBreath2(i,:), 'k',...
%         Spectrogram.Freq,EnvelpBreath(i,:), 'r--');
end

%% Convert SnoreDB back to DB
fldnms = fieldnames(BreathSnoreStruct);
breathswithsnoreDB = find(contains(fldnms, 'SnoreDB'))';
for nn = breathswithsnoreDB
    variable = fldnms{nn};
%     variableNew = strrep(variable,'WinVar','SnoreDB');
    BreathSnoreStruct.(variable) = 10*(log10(BreathSnoreStruct.(variable)/(0.00002^2)));
end

%% combine into one big table
%FT = FeatureNames';
BreathSnoreDataTable = struct2table(BreathSnoreStruct);
% diffsnoretable = (BreathSnoreDataTable{:,1:193} - BreathSnoreTable{:,:})./BreathSnoreTable{:,:};