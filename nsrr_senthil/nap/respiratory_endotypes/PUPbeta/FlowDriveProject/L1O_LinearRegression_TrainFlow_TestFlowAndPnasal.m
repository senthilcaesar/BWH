% in this code, pnasal is tested with betas trained in flow 

%% start
keyboard % this is just to stop the whole process running accidentally (i.e. pressing F5) 
close all
clear
clc

%% options
settings.useLogitFtrs = 0;       % set as 1 to use ftrs selected from logit, 0 for all ftrs
settings.TransformTheData = 1;   % set as 1 to do tranforms, or 0 to use unadjusted data
settings.addextratransform = 0; 	% set as 1 to do extra transforms
settings.UseGOFtoRemoveFtrs = 1;  % set to one to enable (does nothing if threshold below is set to zero)
settings.RsqThreshold = '1/2';    %   '0' '1/3'      '2/3'   '3/4' 
% the not good univariate predictors step must be done in a different dataset
settings.RemoveUnivariateFtrs = 0; % remove ftrs that are not good univariate predictors of GS
settings.experimentnumber = '_n000';
%settings.datadir = 'C:\PSG_Data\FlowDrive\FeatureSpaces\';
%settings.datadir = 'C:\PSG_Data\FlowDrive\Analyzed\ExclBadR\';
settings.datadir = 'C:\PSG_Data\FlowDrive\Analyzed\Publication\';

%% Turn diary logging on
HowSoonIsNow = datestr(datetime('now'),'yyyymmdd_HHMM');
diaryfilename = ['L1O_LinReg_TrainFlowTestFlowPnasal_', HowSoonIsNow, '.txt'];
diary(diaryfilename);
diary on

%% open clean data files
settings.filename_flow = [settings.datadir, 'FlowDrive_only25Hz_ReNormalized_FeatureSpace_AutoRef2_Edi_Clean.mat']; 
settings.filename_pnasal = [settings.datadir, 'PnasalDrive_only25Hz_ReNormalized_FeatureSpace_AutoRef2_Edi_Clean.mat'];
try
    str=['Loading ' settings.filename_pnasal]; disp(str);
    load(settings.filename_pnasal);%,'Amatrix', 'PtData', 'FeatureNames');
    Amatrix_pnasal = Amatrix;
    PtData_pnasal = PtData;
    FeatureNames_pnasal = FeatureNames;
    RemovedBB_Apnoea_pnasal = RemovedBB_Apnoea;
    RemovedBB_LowFlow_pnasal = RemovedBB_LowFlow;    
    str=['Loading ' settings.filename_flow]; disp(str);
    load(settings.filename_flow);
catch me
    disp(me.getReport);
end

%% Remove NaN breaths
% should already be done, but just in case, look for NaN again
allnanrows = sum(isnan(Amatrix),2)==size(Amatrix,2);
if nnz(allnanrows)~=0
    str = ['Removing ', num2str(sum(allnanrows)), ' breaths that contain NaN''s']; disp(str);
    Amatrix(allnanrows,:)=[];
    PtData(allnanrows,:)=[];
    Fnan=sum(isnan(Amatrix)|isinf(Amatrix))/size(Amatrix,1);
    if nnz(Fnan)~=0
        disp('NaN''s or non-finite data remains - further investigation required');
        keyboard
    end
else
    str = ['Zero NaN-breaths were found in Flow data']; disp(str);
end

allnanrows = sum(isnan(Amatrix_pnasal),2)==size(Amatrix_pnasal,2);
if nnz(allnanrows)~=0
    str = ['Removing ', num2str(sum(allnanrows)), ' breaths that contain NaN''s']; disp(str);
    Amatrix_pnasal(allnanrows,:)=[];
    PtData_pnasal(allnanrows,:)=[];
    Fnan=sum(isnan(Amatrix_pnasal)|isinf(Amatrix_pnasal))/size(Amatrix_pnasal,1);
    if nnz(Fnan)~=0
        disp('NaN''s or non-finite data remains in Pnasal - further investigation required');
        keyboard
    end
else
    str = ['Zero NaN-breaths were found in Pnasal data']; disp(str);
end

%% determine pt lists
settings.Pnasal_list = unique(PtData_pnasal.PT);
settings.Flow_list = unique(PtData.PT);
%% make data tables for regression analysis
%   1. Full length flow data (free of artefact, no matching with pnasal)
%   2. Matched flow and pnasal data (free of artefact)
ChannelsList = {'Flow','Pnasal'};
artdirectory = ['C:\PSG_Data\FlowDrive\SourceMat 20171123'];
% read spreadsheet (options worksheet)
[~,~,raw] = xlsread('C:\Users\uqdmann\Dropbox\PUPbeta_git\FlowAndEdiPesQAO\AnalyzeDataSpreadsheet.xlsx',1,'G3:G56');
Pnasal_summary = [];
PtData_flow_clean = []; Amatrix_flow_clean = [];
PtData_flow_matched = []; Amatrix_flow_matched = [];
PtData_pnasal_matched = []; Amatrix_pnasal_matched = [];
try
for subj=1:54 % loop through all subjects, getting Pnasal and Flow breaths
    
    % get flow first
    if ~ismember(subj, settings.Flow_list)
        str=['No Flow data for Pt ', num2str(subj)]; disp(str); continue
    else
        disp(' '); str=['Processing flow data for Pt ', num2str(subj)]; disp(str);
    end
    
    StarttimeSpike = 0; % reset
    fname = char(raw{subj});
    warning ('off','all');
    load(['C:\PSG_Data\FlowDrive\SourceMat 20171123\' fname], 'StarttimeSpike');
    warning ('on','all');

    Isubj=(PtData.PT==subj);  % Isubj is the logical index into current subject (flow)
    PtData_flow_pt = PtData(Isubj,:);
    PtBB_flow_time = PtData_flow_pt.BB_time-StarttimeSpike;
    Amatrix_flow_pt = Amatrix(Isubj,:);
    str = [num2str(length(PtBB_flow_time)), ' breaths with Flow - before art removal']; disp(str);
    
    % Remove Artefact breaths (in flow first)
    % trim flow, keeping only breaths that are non-artifact
    textfilename=[artdirectory '\' fname(1:end-4) '_' ChannelsList{1} '_art.txt'];
    if exist(textfilename,'file')==2
        displaytext = ['Finding artifact in ' ChannelsList{1}]; disp(displaytext);
        [M] = dlmread(textfilename); col1 = M(:,1); col2 = M(:,2);
        FlowBBtoExclude = false(length(PtBB_flow_time),1);
        for i=1:length(col1)
            lefti=col1(i); if lefti<1, lefti=1; end
            righti=col2(i); if righti>max(PtBB_flow_time), righti=max(PtBB_flow_time); end
            FlowBBtoExclude = FlowBBtoExclude | ((PtBB_flow_time>lefti) & (PtBB_flow_time<righti));
        end
        PtData_flow_pt(FlowBBtoExclude,:) = [];
        Amatrix_flow_pt(FlowBBtoExclude,:) = [];
        PtBB_flow_time(FlowBBtoExclude) = [];
    end
    
    PtFlowTime = PtBB_flow_time;
    if ~isempty(PtFlowTime) % empty for pts without breaths (unlikely, probably an error)
        str = [num2str(length(PtFlowTime)), ' breaths with Flow - after art removal']; disp(str);
    end
    
    % make full length clean flow PtData and Amatrix
    PtData_flow_clean = [PtData_flow_clean; PtData_flow_pt];
    Amatrix_flow_clean = [Amatrix_flow_clean; Amatrix_flow_pt]; 
    
    if ~ismember(subj, settings.Pnasal_list)
        str=['No Pnasal data for Pt ', num2str(subj)]; disp(str); continue
    else
        str=['Processing Pnasal data for Pt ', num2str(subj)]; disp(str);
    end
    
    Isubj_pnasal=(PtData_pnasal.PT==subj);  % Isubj is the logical index into current subject (pnasal)
    PtData_pnasal_pt = PtData_pnasal(Isubj_pnasal,:);
    PtBB_pnasal_time = PtData_pnasal_pt.BB_time-StarttimeSpike;
    Amatrix_pnasal_pt = Amatrix_pnasal(Isubj_pnasal,:);
    str = [num2str(length(PtBB_pnasal_time)), ' breaths with Pnasal - before art removal']; disp(str);
        
    % Remove Artefact breaths (in Pnasal)
    % trim flow and pnasal, keep only breaths that are non-artifact in both
        textfilename=[artdirectory '\' fname(1:end-4) '_' ChannelsList{2} '_art.txt'];
        if exist(textfilename,'file')==2
            displaytext = ['Finding artifact in ' ChannelsList{2}]; disp(displaytext);
            [M] = dlmread(textfilename); col1 = M(:,1); col2 = M(:,2);
            PnasalBBtoExclude = false(length(PtBB_pnasal_time),1);
            for i=1:length(col1)
                lefti=col1(i); if lefti<1, lefti=1; end
                righti=col2(i); if righti>max(PtBB_pnasal_time), righti=max(PtBB_pnasal_time); end
                PnasalBBtoExclude = PnasalBBtoExclude | ((PtBB_pnasal_time>lefti) & (PtBB_pnasal_time<righti));
            end
            PtData_pnasal_pt(PnasalBBtoExclude,:) = [];
            Amatrix_pnasal_pt(PnasalBBtoExclude,:) = [];
            PtBB_pnasal_time(PnasalBBtoExclude) = [];
        end

    % match breaths from pnasal to real flow data
    PtPnasalTime = PtBB_pnasal_time;
    if ~isempty(PtPnasalTime) % empty for pts without pnasal
        str = [num2str(length(PtPnasalTime)), ' breaths with Pnasal - after art removal']; disp(str);
        % step through PnasalTime, looking for tolerable matches in FlowTime
        % look for breaths within a given tolerance range
        tolerance = 0.5; % 1/2 of one second
        BB_tolerance = @(T) [T-tolerance,T+tolerance];
        matchedBB_PtoF = NaN(length(PtPnasalTime),1);
        for bb=1:length(matchedBB_PtoF)
            [range]=BB_tolerance(PtPnasalTime(bb)); % get range
            [candidates] = find(PtFlowTime>range(1) & PtFlowTime<range(2));
            switch size(candidates,1)
                case 0
                    %str=['No match for breath ', num2str(bb)]; disp(str);
                    matchedBB_PtoF(bb) = NaN;
                case 1
                    matchedBB_PtoF(bb) = candidates;
                otherwise % more than one
                    str=['Multiple matches for breath ', num2str(bb)]; disp(str);
            end
        end
        str = [num2str(nnz(isfinite(matchedBB_PtoF))), ' breaths matched on time']; disp(str);
    end
    pnasal_exclusions = isnan(matchedBB_PtoF);
    matchedBB_PtoF_clean = matchedBB_PtoF;
    matchedBB_PtoF_clean(pnasal_exclusions)=[];
    
    totalBB = min(length(PtFlowTime),length(PtPnasalTime));
    percentmatched = (nnz(isfinite(matchedBB_PtoF))/totalBB)*100;
    Pnasal_summary = [Pnasal_summary; [subj, length(PtFlowTime),length(PtPnasalTime), nnz(isfinite(matchedBB_PtoF)), percentmatched]];
    
    % get matched output for current pt
    PtData_flow_pt = PtData_flow_pt(matchedBB_PtoF_clean,:);
    Amatrix_flow_pt = Amatrix_flow_pt(matchedBB_PtoF_clean,:);
    PtData_pnasal_pt(pnasal_exclusions,:) = [];
    Amatrix_pnasal_pt(pnasal_exclusions,:) = [];
    
    % make a matched output for PtData and Amatrix, for flow and pnasal
    PtData_flow_matched = [PtData_flow_matched; PtData_flow_pt];
    PtData_pnasal_matched =  [PtData_pnasal_matched; PtData_pnasal_pt];
    
    Amatrix_flow_matched = [Amatrix_flow_matched; Amatrix_flow_pt];
    Amatrix_pnasal_matched = [Amatrix_pnasal_matched; Amatrix_pnasal_pt];
    
end
catch MakeDataTables
    disp(MakeDataTables.getReport);
    keyboard
end

%% summary stats thus far
if 1
    length(unique(PtData_flow_clean.PT)) % number of pts remaining in clean data
    length(unique(PtData_flow_matched.PT)) % number of pts in matched data
    length(unique(PtData_pnasal_matched.PT)) % number of pts in matched data (same as above)
    length(PtData.PT) - length(PtData_flow_clean.PT) % number of flow breaths removed as artefact
end

%% set up breath ID's and merge matched pnasal data with full clean data
B_i_m = (PtData_flow_matched.PT * 1E8) + (PtData_flow_matched.BB_time*1E3); % breath ID's for matched breaths
B_i_f = (PtData_flow_clean.PT * 1E8) + (PtData_flow_clean.BB_time*1E3); % breath ID's for flow breaths

% ensure uniqueness
if length(unique(B_i_m))~=length(B_i_m); keyboard; end
if length(unique(B_i_f))~=length(B_i_f); keyboard; end

% add unique IDs to mats, and make into tables
p = [array2table(B_i_m) PtData_pnasal_matched];
f = [array2table(B_i_f) PtData_flow_clean];
% merge tables using unique IDs
PtData_combined = outerjoin(f, p, 'Keys',[1], 'MergeKeys',1);

% do it again for Amatrix's
p = [array2table(B_i_m) array2table(Amatrix_pnasal_matched)];
f = [array2table(B_i_f) array2table(Amatrix_flow_clean)];
% merge tables using unique IDs
Amatrix_combined = outerjoin(f, p, 'Keys',1, 'MergeKeys',1); % join flow and pnasal, using 1st col from each as key

%% split back into individual mats for flow and pnasal, but now they are row matched
PtData_flow = PtData_combined(:,[2:18]);
PtData_pnasal = PtData_combined(:,[19:35]);
Amatrix_flow = Amatrix_combined{:,[2:166]};
Amatrix_pnasal = Amatrix_combined{:,[167:331]};

%% tidy up
clearvars -except settings PtData_flow PtData_pnasal Amatrix_flow Amatrix_pnasal FeatureNames

%% Gold standard flow:drive 
% this is the wrong way to do this. 
% should change Amatrix instead
ModelTarget = 'flowdrive';
switch ModelTarget
    case 'flowdrive'
        GS = PtData_flow.g_Edi_Adj_f;
    case 'flowonly'
        GS = PtData_flow.VE_f./PtData_flow.Veup_f;
    case 'driveonly'
        GS = ((PtData_flow.DriveEdi_f./100)./PtData_flow.Veup_f); 
end
% double check no nans or ~finite values in Gtest
if length(isfinite(GS))~=length(GS)
    find(~isfinite(GS));
    keyboard
end
maxG=1.5;
GS(GS>maxG)=maxG;

%% Set the starting features
% option one, determined from logistic regression
% keep the original Amatrix (before selecting specific features, and transforms
% keep the origianl FeatureNames variable before overwriting
% same for pnasal, if included
Original_FeatureNames = FeatureNames;
if settings.useLogitFtrs
    %load('FtrsToUse','LabelsOrderedOpt');
    %load('FtrsToUse_Normalized','LabelsOrderedOpt');
    load('FtrsToUse_ReNormalized_2','LabelsOrderedOpt');
    
    Ind = [];
    for i=1:length(LabelsOrderedOpt)
        temp=find(startsWith(FeatureNames.Name,LabelsOrderedOpt(i)));
        if ~isempty(temp)
            Ind(end+1)=temp;
        end
    end
    FeatureNames = FeatureNames(Ind,:);
    Amatrix_flow = Amatrix_flow(:,Ind);
    Amatrix_pnasal = Amatrix_pnasal(:,Ind);  
end

%% Add transforms
if settings.TransformTheData
    [Amatrix2_flow, Labels] = DoDataMatTransform(Amatrix_flow, FeatureNames, settings.addextratransform);
    [Amatrix2_pnasal, ~] = DoDataMatTransform(Amatrix_pnasal, FeatureNames, settings.addextratransform);
else
    if 1 % this is the normal behaviour
        Amatrix2_flow = Amatrix_flow;
        Amatrix2_pnasal = Amatrix_pnasal;
        Labels = FeatureNames.Name;
    else % this is for testing flow only or drive only
        
        % flow only
        %Amatrix2_flow = PtData_flow.VE_f./PtData_flow.Veup_f;
        %Amatrix2_pnasal = PtData_pnasal.VE_p./PtData_pnasal.Veup_p;
        %Labels = 'Flow';
        
        % drive only
        Amatrix2_flow = ((PtData_flow.DriveEdi_f./100)./PtData_flow.Veup_f); 
        Amatrix2_pnasal = ((PtData_pnasal.DriveEdi_p./100)./PtData_pnasal.Veup_p);         
        Labels = 'Drive';
    end
end

%% Set up Weights
settings.dx=0.2;
settings.xbins=[0 0.3:settings.dx:0.9 maxG];
%xbins=[0 0.1:0.05:1.1 1.5];

for i=1:length(settings.xbins)-1
    Ix=GS>settings.xbins(i)&GS<=settings.xbins(i+1);
    Ndata(i)=sum(Ix);
end
weightsbins = 1./(Ndata);
weightsbins = weightsbins/mean(weightsbins);

weights = NaN*GS;
for i=1:length(settings.xbins)-1
    Ix=GS>=settings.xbins(i)&GS<=settings.xbins(i+1);
    weights(Ix)=weightsbins(i);
end
weights = weights/nanmean(weights); % nnz(~isfinite(weights))

settings.useweights=1;
if ~settings.useweights %overwrite, use no weights
    weights = ones(length(weights),1);
end

%% Leave one subject out loop
% linear regression, leave one out cross validation
% train in flow, test in flow and test in pnasal
RemovedFtrs=[];  % list of features that fail correlation test
RsqTrain_array_flow=[];  %
ErrTrain_array_flow=[];  %
predyL1O_array_flow=[];  %
predyL1O_array_pnasal=[];  %
labels_Step_Subj_flow=[];  %
beta_array_flow=[];  %
t_start_L1O = clock;
if settings.RemoveUnivariateFtrs
    Bad_Rsq_vals_i=NaN(54,50);
else
    Bad_Rsq_vals_i=NaN(54,size(Amatrix2_flow,2));
end
try
    warning ('off','all');
    for subj=1:54 % subj=2 subj=3
        if ~ismember(subj, settings.Flow_list)
            str=['No flow data for Pt ', num2str(subj)]; disp(str); continue
        else
            disp(' '); str=['Hold out Pt ', num2str(subj)]; disp(str);
        end
        
        %% set up the L1O subj
        Isubj=(PtData_flow.PT_f==subj);  % Isubj is the logical index of the L1O patient, can use same indx for all now
        
        % set up training data - flow (always)
        Gtest_train = GS(~Isubj);
        weights_train = weights(~Isubj);
        colofones_train = ones(nnz(~Isubj),1);
        Amatrix2_train = Amatrix2_flow(~Isubj,:);
        
        % set up test data - flow and pnasal
        Gtest_test=GS(Isubj);
        weights_test = weights(Isubj);
        colofones_test = ones(nnz(Isubj),1);
        Amatrix2_flow_test = Amatrix2_flow(Isubj,:);
        Amatrix2_pnasal_test = Amatrix2_pnasal(Isubj,:);
        %PtData_flow_test = PtData_flow(Isubj,:);
        %PtData_pnasal_test = PtData_pnasal(Isubj,:);
        
        Ftr_indx_flow = 1:size(Amatrix2_train,2); % this gets reduced n cell below
        
        %% remove features that are not good univariate predictors of GS
        if settings.RemoveUnivariateFtrs
            % find top 5 ftrs univariate ftrs from top 25-ish list
            RsqVsFD = NaN(length(Ftr_indx_flow),1);
            for ft=1:length(Ftr_indx_flow)
                % Rsq vs flow:drive
                FtrVal = Amatrix2_flow(:,Ftr_indx_flow(ft));
                [Rsq,~,~,~]=glmfitFast(FtrVal,GS,weights,1);
                RsqVsFD(ft) = Rsq;
            end
            [ia, ib] = sort(RsqVsFD,'descend');
            Ftr_indx_flow = Ftr_indx_flow(ib(1:50));
        end
        
        %% remove features that do not correlate between flow and pnasal
        % makes most sense for pnasal patients, but can be applied to Flow pts
        if settings.UseGOFtoRemoveFtrs
            % calc correlations for flow and pnasal features
            numofftrs = length(Ftr_indx_flow);
            Rsq_vals = NaN(1,numofftrs);
            for ft=1:numofftrs
                % coefficient of variation, GOF
                [Rsq,~,~,~,~] = glmfitFast(Amatrix2_train(:,ft), Amatrix2_pnasal(~Isubj,ft), weights_train,1);
                Rsq_vals(ft) = Rsq;
            end
            
            % what does the data look like?
            if 0
                figure(); stairs(Rsq_vals);
            end
            
            % set a threshold for what is a 'bad' correlation
            switch settings.RsqThreshold
                case '0'
                    settings.Rsq_Threshold = -Inf;
                case '1/3'
                    settings.Rsq_Threshold = 1/3;
                case '1/2'
                    settings.Rsq_Threshold = 1/2;
                case '2/3'
                    settings.Rsq_Threshold = 2/3;
                case '3/4'
                    settings.Rsq_Threshold = 3/4;
            end
            
            % find bad correlations
            Bad_Rsq_vals_i(subj,:) = Rsq_vals<settings.Rsq_Threshold;
            NumBadFtrsRsq = nnz(Bad_Rsq_vals_i(subj,:));
            str = ['Number of bad ftrs by Rsq: ', num2str(NumBadFtrsRsq)]; disp(str);
            
            % report number of features rejected
            RemovedFtrs = [RemovedFtrs; ...
                [subj, NumBadFtrsRsq]];
            
            % remove bad features
            FtrsToExclude = find(Bad_Rsq_vals_i(subj,:));
            Ftr_indx_flow(FtrsToExclude) = [];
        end
        
        Labels_flow = Labels(Ftr_indx_flow); % Start with the full list, less removed ftrs
        
        %% do backwards elimination
        while ~isempty(Ftr_indx_flow)
            %% train in flow            
            [~,Pvals_train_flow,~,b_train_flow,predytrain_flow]=glmfitFast(Amatrix2_train(:,Ftr_indx_flow),Gtest_train,weights_train,1); %faster
            Pvals_train_flow(1)=[];
            
            %% save the training data - flow
            predytrain_flow(predytrain_flow>maxG)=maxG; predytrain_flow(predytrain_flow<0)=0;
            [RsqTrain_array_flow(subj,length(Ftr_indx_flow)),~,~,~,~] = glmfitFast(predytrain_flow,Gtest_train, weights_train,1); 
            ErrTrain_array_flow(subj,length(Ftr_indx_flow)) = nanmean(weights_train.*abs(predytrain_flow-Gtest_train));
            
            %% use beta from training, apply to test, save data - flow
            predyL1O_flow = [colofones_test Amatrix2_flow_test(:,Ftr_indx_flow)]*b_train_flow;
            predyL1O_flow(predyL1O_flow>maxG)=maxG;
            predyL1O_flow(predyL1O_flow<0) = 0;
            predyL1O_array_flow(Isubj,length(Ftr_indx_flow)) = predyL1O_flow;
            
            %% use beta from training in flow, apply to test data in pnasal, save data - pnasal
            predyL1O_pnasal = [colofones_test Amatrix2_pnasal_test(:,Ftr_indx_flow)]*b_train_flow;
            predyL1O_pnasal(predyL1O_pnasal>maxG)=maxG;
            predyL1O_pnasal(predyL1O_pnasal<0) = 0;
            predyL1O_array_pnasal(Isubj,length(Ftr_indx_flow)) = predyL1O_pnasal;
            
            %% store the labels and beta array - flow
            labels_Step_Subj_flow{subj,length(Ftr_indx_flow)} = Ftr_indx_flow;
            beta_array_flow{subj,length(Ftr_indx_flow)} = b_train_flow;
            
            [maxp_flow,remi_flow]=max(Pvals_train_flow);
            
            if length(Ftr_indx_flow)==1
                break
            end
            
            disp(['F  - Removing Ftr: ', num2str(Ftr_indx_flow(remi_flow)), ', p= ' num2str(maxp_flow)]);
            
            Ftr_indx_flow(remi_flow)=[];
            Labels_flow(remi_flow)=[];
        end
    end
    warning ('on','all');
catch L1O_loop
    disp(L1O_loop.getReport);
    keyboard
end       

% display processing time
delta_t = etime(clock, t_start_L1O); % delta in seconds
D = duration(0,0,round(delta_t),'Format','hh:mm:ss'); % duration in HH:MM:SS
disp(' '); % add row space for visual clarity in command window
displaytext = ['L1O linear regresion process complete. Total time: ', char(D), ' (hh:mm:ss)'];

% use length of RemovedFtrs to trim data arrays
if settings.RemoveUnivariateFtrs == 1
    if isempty(RemovedFtrs)
        maxDataSize = 50;
    else
        maxDataSize = (50 - max(RemovedFtrs(:,2)));
    end
else
    maxDataSize = size(Amatrix2_flow,2) - max(RemovedFtrs(:,2));
end
predyL1O_array_pnasal = predyL1O_array_pnasal(:,[1:maxDataSize]);
predyL1O_array_flow = predyL1O_array_flow(:,[1:maxDataSize]);
ErrTrain_array_flow = ErrTrain_array_flow(:,[1:maxDataSize]);
RsqTrain_array_flow = RsqTrain_array_flow(:,[1:maxDataSize]);
labels_Step_Subj_flow = labels_Step_Subj_flow(:,[1:maxDataSize]);
beta_array_flow = beta_array_flow(:,[1:maxDataSize]);

% Tidy the PtData tables, by removing the _f and _p from each label
PtData_Labels = PtData_flow.Properties.VariableNames';
for n = 1:length(PtData_Labels)
    PtData_Labels{n}=regexprep(PtData_Labels{n,1},'_f','');
end
PtData_flow.Properties.VariableNames = PtData_Labels';
PtData_pnasal.Properties.VariableNames = PtData_Labels';

HowSoonIsNow = datestr(datetime('now'),'yyyymmdd');
savestr = ['_LinRegWorkspace_TrainFlow_TestFlowAndPnasal_', settings.experimentnumber, '.mat'];
str=['Saving to: ', settings.filename_flow(1:end-4), savestr]; disp(str);
save([settings.filename_flow(1:end-4), savestr]);
disp('done');
diary off

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>> END OF RUN <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------





