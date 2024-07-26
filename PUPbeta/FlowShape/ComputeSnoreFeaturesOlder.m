function BreathSnoreDataTable = ComputeSnoreFeaturesOld(Time,DataEventHypnog_Mat,SnoreInterpStruct,SnoreStruct,BB_original,DBthresh)

global ChannelsList

%% initialize variables
dt=(Time(end)-Time(1))/(length(Time)-1);
Fs = round(1/dt);

BB_i_start = BB_original(:,1);
BB_i_end = BB_original(:,2);
BB_e_start = BB_original(:,2);
BB_e_end = BB_original(:,3);  

BBs=length(BB_i_start);

BreathSnoreStruct = struct();

%% identify and sample variables to be analyzed
% [FeatureNames]=MakeSnoreFeaturesList(BBs);
SnoreInterpStruct.SnoreDB = DataEventHypnog_Mat(:,strcmp(ChannelsList,'SnoreDB'));
SnoreVarnames = fieldnames(SnoreInterpStruct);
numrows = size(SnoreInterpStruct.SnoreTime,1);
harmonicIdx = SnoreInterpStruct.HNR > 0.5 & SnoreInterpStruct.HarmonicPower > 0.5;

for i = 1:BBs  
    for j = 1:length(SnoreVarnames)
        varname = SnoreVarnames{j};
        variable = SnoreInterpStruct.(varname);
        if length(variable) > numrows || sum(strcmp({'SnoreTime','FilteredSnd'},varname))==1
            continue % if not the same size then its just a var I smuggled into the struct
        end
        
        if strcmp(varname, 'SnoreDB')
            variable = 10.^((variable/10))*(0.00002^2);
        end
        
        idx_i = BB_i_start(i):BB_i_end(i);
        idx_e = BB_e_start(i):BB_e_end(i);
        if DBthresh > 0
            idx_i = SnoreInterpStruct.SnoreDB(idx_i) > DBthresh;
            idx_e = SnoreInterpStruct.SnoreDB(idx_e) > DBthresh;
        end
        nanadd = nan(1,size(variable,2));
        BreathSnoreStruct.([varname,'Mean_i'])(i,:) = nanmean(variable(idx_i,:),1);
        BreathSnoreStruct.([varname,'Mean_e'])(i,:) = nanmean(variable(idx_e,:),1);
        BreathSnoreStruct.([varname,'Median_i'])(i,:) = nanmedian(variable(idx_i,:),1);
        BreathSnoreStruct.([varname,'Median_e'])(i,:) = nanmedian(variable(idx_e,:),1);
        BreathSnoreStruct.([varname,'Max_i'])(i,:) = max([variable(idx_i,:);nanadd],[],1); % nan bit returns nan if empty
        BreathSnoreStruct.([varname,'Max_e'])(i,:) = max([variable(idx_i,:);nanadd],[],1);
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
        BreathSnoreStruct.percHarmonic(i,:) = sum(harmonicIdx_ & idx_i_)/sum(idx_i_);
        BreathSnoreStruct.percHarmonic_halfI1(i,:) = sum(harmonicIdx_ & halfidx1)/sum(halfidx1);
        BreathSnoreStruct.percHarmonic_halfI2(i,:) = sum(harmonicIdx_ & halfidx2)/sum(halfidx2);

        % pitch of first and second half of insp
        BreathSnoreStruct.Pitch_halfI1(i,:) = nanmean(SnoreInterpStruct.f0(halfidx1));
        BreathSnoreStruct.Pitch_halfI2(i,:) = nanmean(SnoreInterpStruct.f0(halfidx2));
        
        
%         % Need to generate a DBvariable for Window by Window data -
%         Decided against this
%         SnoreDBSm = 10*(log10(SnoreStruct.WinVar/(0.00002^2)));
        % ensemble average PSD RESAMPLE THESE SO THEY CAN BE SAME LENGTH AS
        % ABOVE
%         PSDdata = SnoreStruct.Spectrogram.Power(idx_i,:);
%         PSDbreath(i,:) = nanmean(PSDdata,1);
% 
%         % ensemble average envelope
%         envelp = SnoreStruct.Spectrogram.lpcenv(idx_i,:);
%         EnvelpBreath(i,:) = nanmean(envelp,1);
    end
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