%%
% We have previously used the flow-shape predicted flow:drive to estimate the human categorical scoring of flow-limitation certainty.
% However, doing this implies that there is a linear relationship between the two.
% i.e. some flow:drive severity value maps to some flow-limitiaton certainty value (the log odds version)
% But what if there are different shapes that are related to certainty of obstruction vs severity of obstruction?
% The objective of this work is to see if we can build a better model to predict categorical human scoring of flow-limiation certainty.
% We will let it pick new features, from our pool of features.
% Outcomes of interest:
% 1) The actual features used/selected in the model. Are these the same, some overlap, or completely different?
% 2) The performance of the new model.
%    Better - do we want two models.
%    Same - is it likely to be better on new data, and the performance we are getting with FDpred ftr model may be because of training data
%    Worse - ok, that didn't work.
%
% Routine - Leave-one-out (because it's the most conservative way to assess all)
% Approach - SFFS (because I haven't got all year to do this)
% Method - Ordinal Logistic Regression (becuase the data is ordinal, and we want one slope with two intercepts)
% Optimise - Kappa (because there is an imbalance in categories)

% The SFFS method is:
% do one time mdl fit with all, add 'best' one to current list
% calc current kappa with current list
% while progress
%  add process
%   repeat for N of remaining list
%       cycle through current list plus one from remaining list
%       calc new kappa at each step
%   if new kappa better than current kappa, and
%       add improvement delta > add threshold, then
%   set as current, else add progress = 0
%  remove process
%   repeat for N of current list
%       cycle through current list missing one feature
%       calc new kappa at each step
%   if new kappa better than current kappa, and
%       remove improvement delta > remove threshold, then
%   set as current, else remove progress = 0
% if add progress = 0 and remove progress = 0
% then progress = 0, end process.

%% takes about 5 days to run

%% start
close
clear
clc

%% options
settings.useLogitFtrs = 0;       % set as 1 to use ftrs selected from logit, 0 for all ftrs
settings.TransformTheData = 1;   % set as 1 to do tranforms, or 0 to use unadjusted data
settings.addextratransform = 0;   % set as 1 to do extra transforms
settings.UseGOFtoRemoveFtrs = 1;  % set to one to enable (does nothing if threshold below is set to zero)
settings.RsqThreshold = '1/2';    %   '0' '1/3'      '2/3'   '3/4'
settings.RemoveUnivariateFtrs = 0; % remove ftrs that are not good univariate predictors of GS
settings.experimentnumber = '_retrainManual_SFFS';
settings.datadir = '';
settings.savefigs = 0;
settings.mksize = 15; %15 default
settings.mkcol = [0 0 0];
settings.TickFntSz = 12;
settings.LabelFntSz = 18;
settings.TitleFntSz = 20;
settings.sigdig = '%.2f';
settings.doinitialsetup = 1; 

if settings.doinitialsetup
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
    
    % remove pt 1264, which happens to be number 11, at index 9 in flow list
    settings.Flow_list(9)=[];
    % and doesn't occur in pnasal list
    
    
    %% make data tables for regression analysis
    %   1. Full length flow data (free of artefact, no matching with pnasal)
    %   2. Matched flow and pnasal data (free of artefact)
    ChannelsList = {'Flow','Pnasal'};
    artdirectory = ['C:\PSG_Data\FlowDrive\SourceMat 20171123'];
    % read spreadsheet (options worksheet)
    [~,~,raw] = xlsread('AnalyzeDataSpreadsheet_.xlsx',1,'G3:G56');
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
        numpts = length(unique(PtData_flow_clean.PT)) % number of pts remaining in clean data
        numptsmatched = length(unique(PtData_flow_matched.PT)) % number of pts in matched data
        length(unique(PtData_pnasal_matched.PT)) % number of pts in matched data (same as above)
        numartbb= length(PtData.PT) - length(PtData_flow_clean.PT) % number of flow breaths removed as artefact
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
     
    %% load manual scoring in place of Gold Standard
    disp('');
    PtData_rebuild = []; Amatrix_rebuild = [];
    ManScore = [];
    for subj=1:54 % set pt num or all, no data for PT=1, subj=2 
        if ~ismember(subj, settings.Flow_list)
            str=['No flow data for Pt ', num2str(subj)]; disp(str); continue
        end
        
        str=['Pt ', num2str(subj)]; disp(' '); disp(str);
        
        % PtData is all pts, so just get the pt of interest
        Isubj = PtData_combined.PT_f == subj;      % get the subj BB's
        PtData = PtData_combined(Isubj,:); %
        Amatrix_pt = Amatrix_combined(Isubj,:); % same for Amatrix
        
        % value 0 is FL, 2 is intermediate, 4 is NotFL, 5 is shitty data
        % zeros changed to ones, twos are unchanged, fours changed to threes.
        % fives unchanged.
        
        [ManScore_pt, Pt, Amat] = LoadManualscoring(subj,PtData, Amatrix_pt);
        
        % do matching on PtData
        PtData_rebuild = [PtData_rebuild; Pt];
        Amatrix_rebuild = [Amatrix_rebuild; Amat];
        ManScore = [ManScore; ManScore_pt];
        
    end
    
    % re-assign, with only matched values
    HS = ManScore;
    PtData_combined = PtData_rebuild;
    Amatrix_combined = Amatrix_rebuild;
    
    
    %% split back into individual mats for flow and pnasal, but now they are row matched
    PtData_flow = PtData_combined(:,[2:18]);
    PtData_pnasal = PtData_combined(:,[19:35]);
    Amatrix_flow = Amatrix_combined{:,[2:166]};
    Amatrix_pnasal = Amatrix_combined{:,[167:331]};
    
    %% tidy up and save
    clearvars -except settings PtData_flow PtData_pnasal Amatrix_flow Amatrix_pnasal FeatureNames HS
    if 0 
        save('RetrainingWorkspace.mat');
    end
    
    PT=(PtData_pnasal.PT_p);
    PT(isnan(PT))=[];
    unique(PT) % show the pts with pnasal data
else
    load('RetrainingWorkspace.mat'); 
end

%%
settings.verbose = 1; %
settings.findfirstftr = 1; % can turn this off once trained to increase speed, loads previously calc'd data
settings.ThinDataOut = 0;

% double check no nans or ~finite values in Gtest
if nnz(isfinite(HS))~= length(HS)
    GS_nans = find(~isfinite(HS));
    %keyboard
end

% remove rows with nan gold standard
HS(GS_nans,:)=[];
PtData_flow(GS_nans,:)=[];
PtData_pnasal(GS_nans,:)=[];
Amatrix_flow(GS_nans,:)=[];
Amatrix_pnasal(GS_nans,:)=[];

% remove rows was shitty data
SDI = find(HS==5);
HS(SDI,:)=[];
PtData_flow(SDI,:)=[];
PtData_pnasal(SDI,:)=[];
Amatrix_flow(SDI,:)=[];
Amatrix_pnasal(SDI,:)=[];

%% thin data out to make processing faster
% here we are only using every Nth breath
if settings.ThinDataOut
    settings.ThinFactor = 5;
    thin=[1:settings.ThinFactor:length(HS)];
    HS=HS(thin,:);
    PtData_flow=PtData_flow(thin,:);
    PtData_pnasal=PtData_pnasal(thin,:);
    Amatrix_flow=Amatrix_flow(thin,:);
    Amatrix_pnasal=Amatrix_pnasal(thin,:);
end

%% ensure feature limits (flow 1, pnasal 2)
Amatrix_flow = ApplyFeatureLimits(Amatrix_flow, 1);
Amatrix_pnasal = ApplyFeatureLimits(Amatrix_pnasal, 2);

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


%% not entirely categorical, actually is more ordinal
HS_cat = categorical(HS,'Ordinal',true);
cats = categories(HS_cat)
summary(HS_cat)

%% Set up Weights
%settings.dx=0.2;
%settings.xbins=[0 0.3:settings.dx:0.9 maxG]; %original
%settings.xbins=[0 0.3 0.7 1.2]; % with 0.1, 0.5 and 1.0 labels
settings.xbins=[0.5 1.5 2.5 3.5]; % with 1, 2 and 3 labels
Ndata=[];
for i=1:length(settings.xbins)-1
    Ix=HS>settings.xbins(i)&HS<=settings.xbins(i+1);
    Ndata(i)=sum(Ix);
end
weightsbins = 1./(Ndata);
weightsbins = weightsbins/mean(weightsbins);

weights = NaN*HS;
for i=1:length(settings.xbins)-1
    Ix=HS>=settings.xbins(i)&HS<=settings.xbins(i+1);
    weights(Ix)=weightsbins(i);
end
weights = weights/nanmean(weights); % nnz(~isfinite(weights))
%weightsbins_fudge = flipud(unique(weights));

settings.useweights=1;
if ~settings.useweights %overwrite, use no weights
    weights = ones(length(weights),1);
end

%% Leave one subject out loop
% leave one out cross validation
% train in flow, test in flow and test in pnasal
RemovedFtrs=[];  % list of features that fail correlation test
RsqTrain_array_flow=[];  %
RsqTrain_array_flow_pt=[];  %
ErrTrain_array_flow=[];  %
ErrTrain_array_flow_pt=[];  %
predyL1Ocat_array_flow=NaN(size(Amatrix2_flow,1),1);  %
predyL1Ocat_array_pnasal=NaN(size(Amatrix2_flow,1),1);  %
predyL1Oprob_array_flow=NaN(size(Amatrix2_flow,1),3);  %
predyL1Oprob_array_pnasal=NaN(size(Amatrix2_flow,1),3);  %
labels_Step_Subj_flow={};  %
t_start_L1O = clock;
Bad_Rsq_vals_i=NaN(54,size(Amatrix2_flow,2));

% kappa thresholds
addthreshold = 0.005;       % Todo: modify as necessary
removethreshold = 0.005;    % Todo: modify as necessary
       
kappahistory_all={};
training_kappa=NaN(54,1);
test_kappa_flow=NaN(54,1);
test_kappa_pnasal=NaN(54,1);

Restarting = 0; % this is for troubleshooting

progressbar(0,0,0);
progressbar('Pt','Step','Ftrs');
try
    warning ('off','all');
    for subj=1:54 % subj=2 subj=3  
        if ~ismember(subj, settings.Flow_list)
            str=['No flow data for Pt ', num2str(subj)]; disp(str); continue
        else
            disp(' '); str=['Hold out Pt ', num2str(subj)]; disp(str);
        end
        progressbar(subj/54,[],[]);
        
        kappahistory_pt = []; % reset kappa history within pt
        if Restarting
            Restarting = 0; % reset this flag, so we don't end up back here
            disp('loading existing data... (this may take a while)');
            load('EndSFFSforPt4'); % load the work so far
        else
            %% set up the L1O subj
            Isubj=(PtData_flow.PT_f==subj);  % Isubj is the logical index of the L1O patient, can use same indx for all now
            
            % set up training data - flow (always)
            % training data is everything except the current pt
            Gtest_train = HS_cat(~Isubj);
            weights_train = weights(~Isubj);
            colofones_train = ones(nnz(~Isubj),1);
            Amatrix2_train = Amatrix2_flow(~Isubj,:);
            PtData_flow_train = PtData_flow.PT_f(~Isubj,:);
            
            % set up test data - flow and pnasal
            % test data is only the current pt
            Gtest_test=HS_cat(Isubj);
            weights_test = weights(Isubj);
            colofones_test = ones(nnz(Isubj),1);
            Amatrix2_flow_test = Amatrix2_flow(Isubj,:);
            Amatrix2_pnasal_test = Amatrix2_pnasal(Isubj,:);
            if all(isnan(Amatrix2_pnasal_test)); PnasalPT = 0; else PnasalPT = 1; end
            %PtData_flow_test = PtData_flow(Isubj,:);
            %PtData_pnasal_test = PtData_pnasal(Isubj,:);
            
            Ftr_indx_flow = 1:size(Amatrix2_train,2); % this gets reduced n cell below
            
            %% remove features that do not correlate between flow and pnasal
            % makes most sense for pnasal patients, but can be applied to Flow pts
            progressbar('','GOF','');
            progressbar([],1/10,[]);
            if settings.UseGOFtoRemoveFtrs
                % calc correlations for flow and pnasal features
                numofftrs = length(Ftr_indx_flow);
                Rsq_vals = NaN(1,numofftrs);
                for ft=1:numofftrs
                    progressbar([],[],ft/numofftrs);
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
                
                % keep a record of features rejected
                RemovedFtrs = [RemovedFtrs; ...
                    [subj, NumBadFtrsRsq]];
                
                % remove bad features
                FtrsToExclude = find(Bad_Rsq_vals_i(subj,:));
                Ftr_indx_flow(FtrsToExclude) = [];
            end
            
            Labels_flow = Labels(Ftr_indx_flow); % Start with the full list, less removed ftrs
            
            %% sequential forward floating search
            Ftrs = Ftr_indx_flow;
            %if settings.findfirstftr
                % do one loop through (all this is just to get the starting ftr)
                % side note, this may be helpful when finding 'super-features'
                progressbar('','FirstFtr','');
                progressbar([],2/10,[]);
                kappa_array = nan(length(Ftrs),1);
                for fwdcycle = 1:length(kappa_array)
                    progressbar([],[],fwdcycle/length(kappa_array));
                    templist = Ftrs(fwdcycle);
                    if settings.verbose; str=['Testing ftr #', num2str(templist)]; disp(str); end % include indicator of loop
                    [B1,dev,stats]=mnrfit(Amatrix2_train(:,templist),Gtest_train, 'Model', 'ordinal', 'link','logit');
                    [pihat] = mnrval(B1,Amatrix2_train(:,templist), 'Model', 'ordinal');
                    [prob,cat] = max(pihat(:,:),[],2);
                    [tbl_] = crosstab(Gtest_train,cat);
                    try % if it's not 3x3, then we've got a problem
                        tbl=tbl_.*ones(3);
                    catch % this happens when not all cats are in pred or train,
                        % so, work out what we've got, and what we are missing,
                        % and fill it in so that kappa can work at least a bit
                        tbl=FixMissingDataInTbl(tbl_,cat,Gtest_train);
                    end
                    [kappa_array(fwdcycle), ~] = kappaDLMmod(tbl);
                end
            %    save('FirstFtrKappa','kappa_array');
            %else
            %    load('FirstFtrKappa.mat');
            %end % find first feature loop
            
            [currentkappa, newindx] = max(kappa_array);
            currentlist = Ftrs(newindx);
            if settings.verbose; str=['Starting ftr #', num2str(currentlist)]; disp(str); end
            % add this kappa to the record
            kappahistory_pt = [kappahistory_pt, currentkappa];
            
            try
                progress = 1;
                while progress % while we are making progress
                    progressbar([],(length(currentlist)+2)/10,[]);
                    % try adding a single new ftr
                    %   repeat for N of remaining list
                    %  set up the list of available ftrs, (i.e. all those not in the current list)
                    availftrs = Ftrs; indxout = [];
                    for ind = 1:length(currentlist) % there must be a better way...
                        loc = find(Ftrs==currentlist(ind));
                        indxout = [indxout, loc];
                    end
                    availftrs(indxout)=[];
                    
                    kappa_array = nan(length(availftrs),1);
                    progressbar('',['Test adding ftrs, loop ', num2str(length(currentlist))],'');
                    for fwdcycle = 1:length(kappa_array)
                        progressbar([],[],fwdcycle/length(kappa_array));
                        
                        % cycle through current list plus one from remaining list
                        templist = [currentlist, availftrs(fwdcycle)];
                        if settings.verbose; str=['Testing ftrs + ', num2str(templist)]; disp(str); end
                        
                        [B2,dev,stats]=mnrfit(Amatrix2_train(:,templist),Gtest_train, 'Model', 'ordinal', 'link','logit');
                        [pihat] = mnrval(B2,Amatrix2_train(:,templist), 'Model', 'ordinal');
                        [prob,cat] = max(pihat(:,:),[],2);
                        [tbl_, ~, ~, ~] = crosstab(Gtest_train,cat);
                        try % if it's not 3x3, then we've got a problem
                            tbl=tbl_.*ones(3);
                        catch % this happens when not all cats are in pred or train, this happens in pt 12,
                            % so, work out what we've got, and what we are missing,
                            % and fill it in so that kappa can work at least a bit
                            tbl=FixMissingDataInTbl(tbl_,cat,Gtest_train);
                        end
                        [kappa_array(fwdcycle), ~] = kappaDLMmod(tbl); % get cohens kappa valu
                    end
                    % find best new combination
                    [newkappa, newindx] = max(kappa_array);
                    newlist = [currentlist, availftrs(newindx)];
                    % if new kappa better than current kappa, and
                    %   add improvement delta > add threshold, then
                    %  set as current, else add progress = 0
                    adddelta = newkappa - currentkappa;
                    if newkappa > currentkappa && adddelta > addthreshold
                        currentkappa = newkappa;
                        currentlist = newlist;
                        addprogress = 1;
                        if settings.verbose; str=['Adding ftr ', num2str(availftrs(newindx))]; disp(str); end
                    else
                        addprogress = 0;
                    end
                    
                    % try removing a single feature to see if this improves things
                    %  repeat for N of current list
                    % note, there is actually no point doing this if currentlist is only two long,
                    %  because we've already just tried all the individual features
                    if length(currentlist) > 2
                        kappa_array = nan(length(currentlist),1);
                        progressbar('',['Test removing ftrs, loop ', num2str(length(currentlist))],'');
                        for revcycle = 1:length(kappa_array)
                            progressbar([],[],revcycle/length(kappa_array));
                            
                            % cycle through current list missing one feature
                            templist = currentlist; templist(revcycle) = [];
                            if settings.verbose; str=['Testing ftrs - ', num2str(templist)]; disp(str); end
                            
                            [B3]=mnrfit(Amatrix2_train(:,templist),Gtest_train, 'Model', 'ordinal', 'link','logit');
                            [pihat] = mnrval(B3,Amatrix2_train(:,templist), 'Model', 'ordinal');
                            [prob,cat] = max(pihat(:,:),[],2);
                            [tbl_, ~, ~, ~] = crosstab(Gtest_train,cat);
                            try % if it's not 3x3, then we've got a problem
                                tbl=tbl_.*ones(3);
                            catch % this happens when not all cats are in pred or train, this happens in pt 12,
                                % so, work out what we've got, and what we are missing,
                                % and fill it in so that kappa can work at least a bit
                                tbl=FixMissingDataInTbl(tbl_,cat,Gtest_train);
                            end
                            [kappa_array(revcycle), ~] = kappaDLMmod(tbl); % get cohens kappa value
                        end
                        
                        % find best new combination
                        [newkappa, newindx] = max(kappa_array);
                        newlist = currentlist; newlist(newindx)=[];
                        % if new kappa better than current kappa, and
                        %   remove improvement delta > remove threshold, then
                        %  set as current, else remove progress = 0
                        removedelta = newkappa - currentkappa;
                        if newkappa > currentkappa && removedelta > removethreshold
                            currentkappa = newkappa;
                            currentlist = newlist;
                            removeprogress = 1;
                            if settings.verbose; str=['Removing ftr ', num2str(newlist(newindx))]; disp(str); end
                        else
                            removeprogress = 0;
                        end
                    else
                        removeprogress = 0; %
                    end
                    
                    % if add progress = 0 and remove progress = 0
                    if addprogress == 0 && removeprogress == 0
                        progress = 0; % then progress = 0, end process.
                    else
                        kappahistory_pt = [kappahistory_pt, currentkappa];
                    end
                end % while progress within this subj
            catch SFFS_error
                disp(SFFS_error.getReport)
            end
            progressbar('','PtModel','');
            progressbar([],9/10,[]);
            
            % save kappahistory and currentlist
            % the idea is that if it breaks, we can jump back in.
            kappahistory_all{subj}=kappahistory_pt; % keep it all
            labels_Step_Subj_flow{subj} = currentlist;
            
            %Labels(currentlist)
            %save('EndSFFSforPt4');
        end
             
        %% make one final model for this pt
        [B4]=mnrfit(Amatrix2_train(:,currentlist),Gtest_train, 'Model', 'ordinal', 'link','logit');
        [pihat] = mnrval(B4,Amatrix2_train(:,currentlist), 'Model', 'ordinal');
        
        %% Save training performance
        [prob,cat] = max(pihat(:,:),[],2);
        [tbl_, ~, ~, ~] = crosstab(Gtest_train,cat);
        try % if it's not 3x3, then we've got a problem
            tbl=tbl_.*ones(3);
        catch % if it's not complete, fill it in with zeros where appropriate
            tbl=FixMissingDataInTbl(tbl_,cat,Gtest_train);
        end
        training_kappa(subj) = kappaDLMmod(tbl); % get cohens kappa value
        
        %% apply trained model in test data - flow
        [pihat] = mnrval(B4,Amatrix2_flow_test(:,currentlist), 'Model', 'ordinal');
        [prob,cat] = max(pihat(:,:),[],2);
        [tbl_, ~, ~, ~] = crosstab(Gtest_test,cat);
        try % if it's not 3x3, then we've got a problem
            tbl=tbl_.*ones(3);
        catch % if it's not complete, fill it in with zeros where appropriate
            tbl=FixMissingDataInTbl(tbl_,cat,Gtest_test);
        end
        %Save testing data - flow
        test_kappa_flow(subj) = kappaDLMmod(tbl); % get cohens kappa value
        predyL1Ocat_array_flow(Isubj) = cat;
        predyL1Oprob_array_flow(Isubj,:) = pihat;
        
        if PnasalPT
            %% apply trained model in test data - pnasal
            % need to exclude breaths that don't have matching pnasal
            %pnasalbbs = ~(isnan(PtData_pnasal.PT_p(:,1)));
            % not needed here, because selecting based on pt
                        
            [pihat] = mnrval(B4,Amatrix2_pnasal_test(:,currentlist), 'Model', 'ordinal');
            [prob,cat] = max(pihat(:,:),[],2);
            [tbl_, ~, ~, ~] = crosstab(Gtest_test,cat);
            try % if it's not 3x3, then we've got a problem
                tbl=tbl_.*ones(3);
            catch % if it's not complete, fill it in with zeros where appropriate
                tbl=FixMissingDataInTbl(tbl_,cat,Gtest_test);
            end
            %Save testing data - flow
            test_kappa_pnasal(subj) = kappaDLMmod(tbl); % get cohens kappa value
            predyL1Ocat_array_pnasal(Isubj) = cat;
            predyL1Oprob_array_pnasal(Isubj,:) = pihat;      
            
        end     
        % cut text 1 - moved below
        
%         % bar chart
%         [GS_counts, ~] = hist(Gtest_test);
%         [pred_counts, labels] = hist(categorical(cat));
%         
%         GS_proportions = GS_counts./sum(GS_counts);
%        
%         pred_proportions = pred_counts./sum(pred_counts);
%         
%         figure()
%         axis off;
%         bar([1 ,2], [GS_proportions; pred_proportions], 'stacked');
%         xticklabels({'GS','Pred'});
%         yticklabels({});
%         box off
       
        
    end % pt loop 
    warning ('on','all');
catch L1O_loop
    keyboard
    disp(L1O_loop.getReport);
end

% display processing time
delta_t = etime(clock, t_start_L1O); % delta in seconds
D = duration(0,0,round(delta_t),'Format','hh:mm:ss'); % duration in HH:MM:SS
disp(' '); % add row space for visual clarity in command window
displaytext = ['L1O linear regresion process complete. Total time: ', char(D), ' (hh:mm:ss)'];
disp(displaytext);

% Tidy the PtData tables, by removing the _f and _p from each label
PtData_Labels = PtData_flow.Properties.VariableNames';
for n = 1:length(PtData_Labels)
    PtData_Labels{n}=regexprep(PtData_Labels{n,1},'_f','');
end
PtData_flow.Properties.VariableNames = PtData_Labels';
PtData_pnasal.Properties.VariableNames = PtData_Labels';

% HowSoonIsNow = datestr(datetime('now'),'yyyymmdd');
% savestr = ['_LinRegWorkspace_TrainFlow_TestFlowAndPnasal_', settings.experimentnumber, '.mat'];
% str=['Saving to: ', settings.filename_flow(1:end-4), savestr]; disp(str);
% save(['C:\tempdata\',settings.filename_flow(1:end-4), savestr]);

save('EndofSFFSrun_OLR_TrainFlowTestFlowandPnasal_withPihat');

disp('done');

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>> END OF RUN <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% ------------------------------------------------------------------------


if 0
close; clear; clc
tic
disp('Loading a big file... (I''ll need a few moments)');
pt5 = load('EndofSFFSrun_PT5');
load('EndofSFFSrun');
disp('OK, file loaded, carry on');
toc

%% add pt5 back into main data
kappahistory_all{5}=pt5.kappahistory_all{5};
labels_Step_Subj_flow{5}=pt5.labels_Step_Subj_flow{5};
Isubj=(PtData_flow.PT==5);
predyL1Ocat_array_flow(Isubj)=pt5.predyL1Ocat_array_flow(Isubj);
predyL1Ocat_array_pnasal(Isubj)=pt5.predyL1Ocat_array_pnasal(Isubj);
predyL1Oprob_array_flow(Isubj)=pt5.predyL1Oprob_array_flow(Isubj);
predyL1Oprob_array_pnasal(Isubj)=pt5.predyL1Oprob_array_pnasal(Isubj);

training_kappa(5) = pt5.training_kappa(5);
test_kappa_pnasal(5) = pt5.test_kappa_pnasal(5);
test_kappa_flow(5) = pt5.test_kappa_flow(5);
end

%%

close; clear; clc
tic

addpath(genpath('M:\Dropbox\QAO\FLcertainty\RetrainCertaintyFromFeaturesUsingSFFS'));

disp('Loading a big file... (I''ll need a few moments)');
load('EndofSFFSrun_OLR_TrainFlowTestFlowandPnasal.mat');
disp('OK, file loaded, carry on');
toc

% also load AHI data
load('AHI_info');

%% post processing for SFFS
mksize = 15; %15 default
mkcol = [0 0 0];
TickFntSz = 12;
LabelFntSz = 18;
TitleFntSz = 20;
sigdig = '%.2f';

savefigs = 0;
closefigs = 1;

%% how often is stable breathing flow limited

%% how many arousals (%age) have flow limitation

ArousalBB = PtData_flow.Ar;
HS_duringAr = HS_cat(logical(ArousalBB));

histogram(HS_duringAr)

gt = array2table(HS_duringAr,'VariableNames',{'Var1'});
g = summary(gt);
GS_cats = g.Var1.Counts';
GS_props = GS_cats/(sum(GS_cats));
% overall, of breaths marked as Ar
% 5.79% were certain-FL
% 21.14% were possible-FL, and
% 73.07% were normal (non-FL).

%% how many arousals (%age) have flow limitation - per pt
GS_inAr_cats = NaN(54,3);
GS_inAr_props = NaN(54,3);
Pred_inAr_cats = NaN(54,3);
Pred_inAr_props = NaN(54,3);

for subj=1:54 % loop through all subjects
    if ismember(subj, settings.Flow_list)
    Isubj=(PtData_flow.PT==subj) & (PtData_flow.Ar==1);  % Isubj is the logical index into current subject (flow)
    %nnz(Isubj)
    HS_duringAr = HS_cat(Isubj);
    gt = array2table(HS_duringAr,'VariableNames',{'Var1'});
    g = summary(gt);
    GS_inAr_cats(subj,:) = g.Var1.Counts';
    GS_inAr_props(subj,:) = GS_inAr_cats(subj,:)/(sum(GS_inAr_cats(subj,:)));
    
    % and again for model cats
    Pred_duringAr = predyL1Ocat_array_flow(Isubj);
    Predt = array2table(categorical(Pred_duringAr),'VariableNames',{'Var1'});
    Pred = summary(Predt);
    Pred_inAr_cats(subj,:) = Pred.Var1.Counts';
    Pred_inAr_props(subj,:) = Pred_inAr_cats(subj,:)/(sum(Pred_inAr_cats(subj,:)));
    
    end
end

GS_inAr_cats = GS_inAr_cats(settings.Flow_list,:);
GS_inAr_props = GS_inAr_props(settings.Flow_list,:);

Pred_inAr_cats = Pred_inAr_cats(settings.Flow_list,:);
Pred_inAr_props = Pred_inAr_props(settings.Flow_list,:);

%%
figure(11); clf(figure(11));
fig = gcf; fig.Color=[1 1 1];
fig.Units = 'Inches'; fig.Position = [10 3.5 11.5 9.5];
rows=7; cols=6;
for i=1:40
    subj = If(i);
    GS_inArPt = GS_inAr_props(subj,:); 
    subplot(rows, cols, i);
    axis off;
    bar([GS_inArPt]); hold on;
    xticklabels({'CFL','PFL','N'}); %yticklabels({});
    %titlestr = ['FL score in AR']; title(titlestr);
    titlestr = ['AHI = ', num2str(AHIflow(subj),2)]; title(titlestr);
    box off
end
text(8, 0.5, [{'Each plot is all Ar breaths from one Pt.'},...
        {'CFL certain-FL'},...
        {'PFL possible-FL'},...
        {'N   normal'}]);

suptitle('Visual scoring of flow-limitation during Arousals'); 

% saveas(fig, ['Figures\FL_during_Arousals_withAHI'], 'png');

%% and again, but stacked bars, one showing visual, the other showing model
figure(1); clf(figure(1));
fig = gcf; fig.Color=[1 1 1];
fig.Units = 'Inches'; fig.Position = [10 3.5 11.5 9.5];
rows=7; cols=6;
for i=1:40
    subj = If(i);
    GS_inArPt = GS_inAr_props(subj,:); 
    Pred_inArPt = Pred_inAr_props(subj,:); 
    subplot(rows, cols, i);
    axis off;
    
    bar([1, 2], [GS_inArPt; Pred_inArPt], 'stacked'); hold on;
    xticklabels({'HS','Pred'}); yticklabels({});
    titlestr = ['AHI = ', num2str(AHIflow(subj),2)]; title(titlestr);

    box off
end

text(5, 0.5, [{'Each plot is all Ar breaths from one Pt.'},...
         {'Blue is certain-FL'},...
         {'Red is possible-FL'},...
         {'Yellow is normal'}]);

% suptitle('Visual scoring of flow-limitation during Arousals'); 

% saveas(fig, ['Figures\FL_during_Arousals_withAHI_HSandPred'], 'png');


%%
figurenum=100;
titlestr = ['Proportion of normal breaths during arousals _ visual scoring and model output'];
save = 1;
boxplotover = 0;

doBeeswarmplot(GS_inAr_props(:,3), Pred_inAr_props(:,3), figurenum, titlestr, save, boxplotover);

%% get real flow for the two extremes
% OpenInAr, 1237 @ 19847 secs + 45 seconds
% ClosedInAr, 1469 @ 23825secs + 50 seconds

%FName1 = ['E:\PSG_Data\FlowDrive\Complete Studies\DPW1237SSPhenotype\1237.smr']; % Normal in Ar
%FName2 = ['E:\PSG_Data\FlowDrive\Complete Studies\DPW1469SSPhenotype\1469Phenotype.smr']; % FL in Ar

FName1 = 'E:\PSG_Data\FlowDrive\Converted\1237_XHz.mat';  % Normal in Ar
FName2 = 'E:\PSG_Data\FlowDrive\Converted\1469_XHz.mat';  % FL in Ar

F1 = load(FName1);
F2 = load(FName2);

dt = F1.DataEventHypnog_Mat(2,1)-F1.DataEventHypnog_Mat(1,1);
Fs = 1/dt;

flow_1_full = F1.DataEventHypnog_Mat(:,2); clear F1;
flow_2_full = F2.DataEventHypnog_Mat(:,2); clear F2;

F1T1 = single(19847*Fs);
F1T2 = floor(43.65*Fs)+F1T1;

F2T1 = single(23825*Fs);
F2T2 = floor(50.25*Fs)+F2T1;

flow_1 = flow_1_full(F1T1:F1T2);
flow_2 = flow_2_full(F2T1:F2T2);



%%
figure(101); clf(figure(101));
fig = gcf; fig.Color=[1 1 1];
fig.Units = 'Inches'; fig.Position = [28 3 11.5 9.5];

subplot(2,1,1)
plot(smooth(flow_1,3));
box off
ylim([-0.75 0.75])
ax = gca; ax.TickDir = 'out';
ylabel('Flow (L/sec)', 'FontSize', 14, 'FontName', 'Arial Narrow');
ax.XColor = [1 1 1];

subplot(2,1,2)
plot(smooth(flow_2,5));
box off
ylim([-0.75 0.75])
ax = gca; ax.TickDir = 'out';
ylabel('Flow (L/sec)', 'FontSize', 14, 'FontName', 'Arial Narrow');
ax.XColor = [1 1 1];

ScaleBar = floor(Fs*10);
SBarStart = 5000;
line([SBarStart+1:(SBarStart+ScaleBar)], ones(ScaleBar,1)*-0.6, 'LineWidth', 2, 'Color', 'k')
text(SBarStart+(0.5*ScaleBar), -0.65, '10 Seconds', 'HorizontalAlignment', 'center')


% saveas(fig, ['Figures\FlowShapesDuringArousals'], 'png');


%% Plot- vertical stack of proportions of normal breaths, and real flow shapes
figure(103); clf(figure(103));
fig = gcf; fig.Color=[1 1 1];
fig.Units = 'Inches'; fig.Position = [28 3 16 6.3];

subplot(2,4,[1 5])

plotXdata=[GS_inAr_props(:,3); Pred_inAr_props(:,3)];
plotYdata= [zeros(length(GS_inAr_props(:,3)),1); ones(length(Pred_inAr_props(:,3)),1)];
LabelFntSz = 14;
facealpha = 1;
scatter(1+(1*plotYdata+0.025*randn(length(plotYdata),1)),100*plotXdata,10,...
    'filled','markerfacealpha',facealpha, 'markerfacecolor', [0 0 0]);

ax = gca; set(ax,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', LabelFntSz);
xticks([1,2]);
xlim([0.5, 2.5]);
xticklabels({'Visual Scoring','Model Predicted'});
ax.YTick=[0:25:100]; 
ylabel('Proportion of ''Normal'' breaths during Arousals (%)', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');

subplot(2,4,[2:4])
plot(smooth(flow_1,3));
box off
ylim([-0.75 0.75])
ax = gca; ax.TickDir = 'out';
ylabel('Flow (L/sec)', 'FontSize', 14, 'FontName', 'Arial Narrow');
ax.XColor = [1 1 1];

subplot(2,4,[6:8])
plot(smooth(flow_2,5));
box off
ylim([-0.75 0.75])
ax = gca; ax.TickDir = 'out';
ylabel('Flow (L/sec)', 'FontSize', 14, 'FontName', 'Arial Narrow');
ax.XColor = [1 1 1];

ScaleBar = floor(Fs*10);
SBarStart = 5000;
line([SBarStart+1:(SBarStart+ScaleBar)], ones(ScaleBar,1)*-0.6, 'LineWidth', 2, 'Color', 'k')
text(SBarStart+(0.5*ScaleBar), -0.65, '10 Seconds', 'HorizontalAlignment', 'center')


% saveas(fig, ['Figures\ProportionNormalInAr_and_FlowShapesDuringAr'], 'png');





%% convert cell array of labels per step to array/matrix
labelsperstep = [];
numftrsused = NaN(length(labels_Step_Subj_flow),1);

for j=1:length(labels_Step_Subj_flow)
    if isempty(labels_Step_Subj_flow{j})
        continue
    end
    labelsperstep = [labelsperstep, labels_Step_Subj_flow{j}];
    numftrsused(j) = length(labels_Step_Subj_flow{j});
    
end
nansum(numftrsused) % should be same length as labelsperstep
emptyarray = [];
FtrsPerSubj = nan(length(labels_Step_Subj_flow),max(numftrsused)); % apriori knowledge of numftsused
for j=1:length(labels_Step_Subj_flow)
    if isempty(labels_Step_Subj_flow{j})
        emptyarray = [emptyarray, j];
        continue
    end
    steps = length(labels_Step_Subj_flow{j});
    FtrsPerSubj(j,1:steps) = labels_Step_Subj_flow{j};
end
missingpt = intersect(emptyarray, settings.Flow_list) % this identifies the missing pt

%% remove nan rows
FtrsPerSubj = FtrsPerSubj(settings.Flow_list,:);

%% determine the most frequently used features
uniquelabels = unique(labelsperstep);
uniquelabels2 = unique(FtrsPerSubj); uniquelabels2(isnan(uniquelabels2))=[];
isequal(uniquelabels, uniquelabels2')
histdata = [];
for j=1:length(uniquelabels)
    histdata(j) = nansum(labelsperstep(:)==uniquelabels(j));
end

clipatftrs = 23; % was 29;
labelsfreq = [uniquelabels', histdata'];
[a,i]=sort(labelsfreq(:,2),1,'descend');
topFtrsfromSFFS = labelsfreq(i(1:clipatftrs),:); % clip at 10, but if they are all ones, then how to rank?
featuresWeLike = Labels(topFtrsfromSFFS(:,1))

%% Histogram of ftrs selected during training
figure(33); clf(figure(33)); fig = gcf;
fig.Color = [1 1 1]; fig.Units = 'inches';
fig.Position = [2.2   2   12    4.5];
bar(topFtrsfromSFFS(:,2), 'facealpha',0.2, 'edgealpha', 0.2);
xlim([0 clipatftrs+1]);
ax = gca; ax.TickDir = 'out'; ax.FontSize=9; ax.XTick=[1:1:60]; xtickangle(ax,90);

lbls_w = featuresWeLike;
lbls = regexprep(lbls_w,'[_,:{}]','');
set(ax, 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
xticklabels(ax,lbls);
ylabel('Frequency', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
box off
title(['Histogram of Features']);
str = ['..\Figures\CorrelationTesting\', 'HistogramOfLinRegPnasalFeatures_WeightedCount',settings.experimentnumber]; %

if savefigs
    saveas(fig, ['Figures\FeatureSelectedDuringSFFS_OLR'], 'png');
end

%% histogram of number ftrs to stop training (by SFFS method)
if length(numftrsused) ~= length(settings.Flow_list)
    numftrsused = numftrsused(settings.Flow_list);
end

uniquenums = unique(numftrsused);
uniquenums(isnan(uniquenums))=[];
histdatan = [];
for j=1:length(uniquenums)
    histdatan(j) = sum(numftrsused(:)==uniquenums(j));
end
ftrnumfreq = [uniquenums, histdatan'];

figure(34); clf(figure(34)); fig = gcf;
fig.Color = [1 1 1]; fig.Units = 'inches';
fig.Position = [1   4   12    4.5];
subplot(1,2,1)
ax = gca; ax.FontSize=9;
bar(ftrnumfreq(:,2), 'facealpha',0.2, 'edgealpha', 0.2); hold on;
xlim([0 length(uniquenums)+1]);
Ftrsfreq = num2cell(ftrnumfreq(:,1)); xticklabels(ax, Ftrsfreq);
xlabel('Number of features selected during training', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
ylabel('Frequency', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
 ax.TickDir = 'out';  %ax.XTick=[1:1:60]; %xtickangle(ax,90);
box off
%title(['Histogram of Number of Features used in SFFS']);
if savefigs
    str = ['Figures\HistogramOfFtrsUsedInSFFS_OLR']; %
    saveas(fig, str, 'png');
end


%% analysis of kappa during training
KappaPerSubj = nan(length(kappahistory_all),10); % guess max number of kappa 
for j=1:length(kappahistory_all)
    if isempty(kappahistory_all{j})
        continue
    end
    steps = length(kappahistory_all{j});
    KappaPerSubj(j,1:steps) = kappahistory_all{j};
end
% remove nan rows and extra cols
KappaPerSubj = KappaPerSubj(settings.Flow_list,:);
KappaPerSubj = KappaPerSubj(:,1:7);
%nanmedian(KappaPerSubj,1)
%nanstd(KappaPerSubj,1)
figure(34); 
% clf(figure(34)); fig = gcf;
% fig.Color = [1 1 1]; fig.Units = 'inches';
% fig.Position = [1   4   12    4.5];
subplot(1,2,2)
ax = gca; ax.FontSize=9;
boxplot(KappaPerSubj) % shows the average kappa per added feature
ylim([0.38 0.57]); 
ylabel('Cohens Kappa', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow'); 
xlabel('Number of features in training model', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
ax.TickDir = 'out'; %ax.XTick=[1:1:60]; %xtickangle(ax,90);
box off

fig = gcf;
if savefigs
    %saveas(fig, ['Figures\KappaPerFeatureAddedSFFS_OLR'], 'png');
    saveas(fig, ['Figures\ForPaper_Sup_HistogramandKappaPerFeatureAddedSFFS_OLR'], 'png');
end


%% analysis of what features occured at each step of cross-validation
clearvars C1 C2
% was x = [1:1:42]; y = [1:1:29];
x = [1:1:41]; y = [1:1:23];
C1 = zeros(length(y),length(x));
pt = 0;
for j=1:length(labels_Step_Subj_flow)
    if isempty(labels_Step_Subj_flow{j})
        continue
    end
    pt=pt+1;
    t = labels_Step_Subj_flow{j};
    indx = ismember(uniquelabels,t);
    C1(indx,pt) = 1;
end
nnz(C1) % should be same length as labelsperstep(:)
C1 = ~C1; % invert

% set up a grid, unique pts wide, unique labels high
figure(35); clf(figure(35)); fig = gcf;
fig.Color = [1 1 1]; fig.Units = 'inches';
fig.Position = [1  0.5   18    8];
customcmap = GetCustomColorMap('gray'); % 'SS'

xSplit = diff(x)/2;                 % Find edge points
ySplit = diff(y)/2;
xEdges = [x(1)-xSplit(1) x(2:end)-xSplit x(end)+xSplit(end)];
yEdges = [y(1)-ySplit(1) y(2:end)-ySplit y(end)+ySplit(end)];
[XGrid, YGrid] = meshgrid(xEdges,yEdges);
YGrid = flipud(YGrid);              % To match expected behavior
XGrid = fliplr(XGrid);
%C1 = rand(length(y),length(x));
C2 = [[C1 zeros(size(C1,1),1)] ; zeros(1,size(C1,2)+1)];% Last row/col ignored
%pcolor(XGrid,YGrid,(1-(1-C2/100).^2)*100)
pcolor(XGrid,YGrid,C2)
hold on                             % Plot original data points
[X,Y] = meshgrid(x,y);
set(gcf,'colormap',customcmap);

ax = gca; set(ax,'fontname','arial narrow','FontSize', TickFntSz);
yticks(y);
if 0 % the ftr index number
    yyaxis right
    yticks(0:1./(length(yEdges)-1):1);
    yticklabels(gca,flipud(num2str(uniquelabels')));
else % the name
    UniqueLabelStrings = Labels(uniquelabels);
    UniqueLabelStrings = regexprep(UniqueLabelStrings,'[_,:{}]','');
    yticklabels(gca,flipud(UniqueLabelStrings));
end
xticks(x); %xticklabels(gca,fliplr(labeltext));
xlabel('cross-validation step', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
ylabel('feature terms', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');

if savefigs
    saveas(fig, ['Figures\PrepForPaper_featuresduringcrossvalidaiton'], 'png');
end

%%
% feature categories
% 1 - Flattening
% 2 - Scooping
% 3 - Asymmetry
% 4 - Timing and Vol ratio
% 5 - Fluttering
clearvars FtrCats
FtrCats = table(uniquelabels', UniqueLabelStrings, 'VariableNames', {'FtrNum','FtrName'});
FtrCats.Category=NaN(height(FtrCats),1);
% FtrCats.Category=[  4     1     1     2     4     4     4     1     4     1     1     2     2 ...
%      2     2     2     4     5     5     5     5     2     2     4     5     4 ...
%      4     1     1]';
%  
FtrCats.Category=[4 2 4 1 4 1 1 2 2 2 2 2 5 5 5 2 2 4 4 5 5 1 1]';
FtrCats.Frequency = labelsfreq(:,2);
%FtrCatsSorted = sortrows(FtrCats,{'Category'});
FtrCatsSorted = sortrows(FtrCats,{'Category', 'Frequency'}, {'ascend','descend'} );

%FtrCatsSorted= FtrCatsSorted(FtrCatsSorted.Frequency>=11,:);


%% analysis of what features occured at each step of cross-validation, but...
% now we are sorting by ftr category

clearvars C1 C2
x = [1:1:40]; y = [1:1:height(FtrCatsSorted)];
C1 = zeros(length(y),length(x));
pt = 0;
for j=1:length(labels_Step_Subj_flow)
    if isempty(labels_Step_Subj_flow{j})
        continue
    end
    pt=pt+1;
    t = labels_Step_Subj_flow{j};
    %indx = ismember(FtrCatsSortedTrimmed.FtrNum,t);
    indx = ismember(FtrCatsSorted.FtrNum,t);
    C1(indx,pt) = 1;
end
nnz(C1) % should be same length as labelsperstep(:)
C1 = ~C1; % invert

% set up a grid, unique pts wide, unique labels high
figure(35); clf(figure(35)); fig = gcf;
fig.Color = [1 1 1]; fig.Units = 'inches';
fig.Position = [1  0.5   10    8];
%fig.Position = [1  0.5   10    3.5]; % for short list
customcmap = GetCustomColorMap('gray'); % 'SS'

xSplit = diff(x)/2;                 % Find edge points
ySplit = diff(y)/2;
xEdges = [x(1)-xSplit(1) x(2:end)-xSplit x(end)+xSplit(end)];
yEdges = [y(1)-ySplit(1) y(2:end)-ySplit y(end)+ySplit(end)];
[XGrid, YGrid] = meshgrid(xEdges,yEdges);
YGrid = flipud(YGrid);              % To match expected behavior
XGrid = fliplr(XGrid);
%C1 = rand(length(y),length(x));
C2 = [[C1 zeros(size(C1,1),1)] ; zeros(1,size(C1,2)+1)];% Last row/col ignored
%pcolor(XGrid,YGrid,(1-(1-C2/100).^2)*100)
pcolor(XGrid,YGrid,C2)
hold on                             % Plot original data points
[X,Y] = meshgrid(x,y);
set(gcf,'colormap',customcmap);
xlim([0.5 40.5])
ax = gca; set(ax,'fontname','arial narrow','FontSize', TickFntSz);
yticks(y);
if 0 % the ftr index number
    yyaxis right
    yticks([0:1./(height(FtrCatsSorted)):1]+(0.5*(1./(height(FtrCatsSorted)))));
    yticklabels(gca,flipud(num2str(FtrCatsSorted.Category)));   
    ax = gca; ax.YColor = [0 0 0];
    ylabel('Feature Category','Color', [0 0 0]);
else % the name
    UniqueLabelStrings = Labels(FtrCatsSorted.FtrNum);
    UniqueLabelStrings = regexprep(UniqueLabelStrings,'[_,:{}]','');
    yticklabels(gca,flipud(UniqueLabelStrings));
    %newy=yticks();
end
%xticks(x); %xticklabels(gca,fliplr(labeltext));
xticks([]);
xlabel('Cross-validation step', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
ylabel('Feature terms', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');

if savefigs
   % saveas(fig, ['Figures\ForPaper_Sup_Microarray_highfreqftrs'], 'png');
    
  %  saveas(fig, ['Figures\ForPaper_Sup_Microarray_allftrs'], 'png');
end


%% cross-validated performance and within-one performance
flow_tbl = crosstab(HS_cat, predyL1Ocat_array_flow); sum(flow_tbl(:))
k=kappaDLMmod(flow_tbl)

pnasalbbs = ~(isnan(PtData_pnasal.PT(:,1)));
pnasal_tbl = crosstab(HS_cat(pnasalbbs), predyL1Ocat_array_pnasal(pnasalbbs));
k=kappaDLMmod(pnasal_tbl)

disagreement_mat = ~(single(eye(3)));
disagreement_mat_flow = disagreement_mat.*flow_tbl;
sum(disagreement_mat_flow(:))/sum(flow_tbl(:))

disagreement_mat_pnasal = disagreement_mat.*pnasal_tbl;
sum(disagreement_mat_pnasal(:))/sum(pnasal_tbl(:))

within_one_mat = single(eye(3));
within_one_mat(1,2) = 1; within_one_mat(2,1) = 1;
within_one_mat(2,3) = 1; within_one_mat(3,2) = 1;
within_one_flow = flow_tbl.*within_one_mat;
sum(within_one_flow(:))/sum(flow_tbl(:))

within_one_pnasal = pnasal_tbl.*within_one_mat;
sum(within_one_pnasal(:))/sum(pnasal_tbl(:))

%% cross-validated test performance
[muHat,sigmaHat,muCI,sigmaCI]=normfit(test_kappa_flow(~isnan(test_kappa_flow)))

[muHat,sigmaHat,muCI,sigmaCI]=normfit(test_kappa_pnasal(~isnan(test_kappa_pnasal)))


%% training performance
[muHat,sigmaHat,muCI,sigmaCI]=normfit(training_kappa(~isnan(training_kappa)))


%% make a final model, using all 23 ftrs
[B23,~,B23stats]=mnrfit(Amatrix2_flow(:,topFtrsfromSFFS(:,1)),HS_cat, 'Model', 'ordinal', 'link','logit');
[pihatB23_f] = mnrval(B23,Amatrix2_flow(:,topFtrsfromSFFS(:,1)), 'Model', 'ordinal');
[~,cat_f] = max(pihatB23_f(:,:),[],2); pred_sffs_23_flow = categorical(cat_f);
tbl_f_23 = crosstab(HS_cat, pred_sffs_23_flow);
k_f_23=kappaDLMmod(tbl_f_23)
% kappa 0.5684

[B23p]=mnrfit(Amatrix2_pnasal(pnasalbbs,topFtrsfromSFFS(:,1)),HS_cat(pnasalbbs), 'Model', 'ordinal', 'link','logit');

[pihatB23_p] = mnrval(B23,Amatrix2_pnasal(pnasalbbs,topFtrsfromSFFS(:,1)), 'Model', 'ordinal');
[~,cat_p] = max(pihatB23_p(:,:),[],2); pred_sffs_23_pnasal = categorical(cat_p);
tbl_p_23 = crosstab(HS_cat(pnasalbbs), pred_sffs_23_pnasal);
k_p_23=kappaDLMmod(tbl_p_23)
% kappa
% with betas fit from pnasal 0.4743
% with betas fit from flow 0.4602


B23stats.beta
FLcert

%% get FD as cat

% train same features to detect FD

% get kappa, ROC and binary acc

logit = @(p) log(p./(1-p));
logitinverse = @(p) 1./(1+exp(-p)); % also called sigmoid

out = logit(1-(pihatB23_f(:,1)));

figure(1); clf(figure(1));
histogram(out(HS_cat=='1')); hold on
histogram(out(HS_cat=='2'));
histogram(out(HS_cat=='3'));

temp = get(gca, 'xtick');
temp_log = logitinverse(temp);

summary(HS_cat)

%%
labels = (HS_cat=='1')*1;
x = out;

% trim out 2
I = HS_cat=='2';
x(I) = [];
labels(I)=[];
addpath('M:\Dropbox\PUPbeta_git\PUPbeta\Predict');
showfigures = 1; figure(100);
[thresX,AUC,SEM,p,posclass,SensSpec]=ROCAUCSEM(labels,x,showfigures);

out2 = out-thresX;

figure(2); clf(figure(2));
histogram(out2(HS_cat=='1')); hold on
histogram(out2(HS_cat=='2'));
histogram(out2(HS_cat=='3'));

% here 0 is the threshold of optimal 
% probability of CFl as new x-axis
% probability of Normal as another new x-axis

out3 = logit(pihatB23_f(:,3)) + (0.5*(B23(2)-B23(1)));

figure(3); clf(figure(3));
histogram(out3(HS_cat=='1')); hold on
histogram(out3(HS_cat=='2'));
histogram(out3(HS_cat=='3'));

% Probability of flow limitation (certain vs normal), try to ctr the data

% is the offset the same, if same, can split offset ...
Temp1 = logit(1-(pihatB23_f(:,1)))
Temp3 = logit((pihatB23_f(:,3)))
Temp1-Temp3 % yes
out4 = logit((pihatB23_f(:,3))) 

%%
% set up proportions of the FD to use as training categories
FD = PtData_flow.g_Edi_Adj;
FD_cat = NaN(length(HS_cat),1);
catscount = sum([HS_cat=='1', HS_cat=='2',HS_cat=='3'])
catscountcum = cumsum(catscount)
cutoffs = [catscountcum(1)/catscountcum(3),catscountcum(2)/catscountcum(3)]
FDcuts = prctile(FD, cutoffs*100);

FD_cat(FD<FDcuts(1))=1;
FD_cat((FD>=FDcuts(1)) & (FD<FDcuts(2)))=2;
FD_cat(FD>=FDcuts(2))=3;

FD_cat = categorical(FD_cat);

% make a final model, using all 23 ftrs, but trained against FD_cat, not HS_cat
[B23_fd]=mnrfit(Amatrix2_flow(:,topFtrsfromSFFS(:,1)),FD_cat, 'Model', 'ordinal', 'link','logit');
[pihatB23_f_fd] = mnrval(B23_fd,Amatrix2_flow(:,topFtrsfromSFFS(:,1)), 'Model', 'ordinal');
[~,cat_f_fd] = max(pihatB23_f_fd(:,:),[],2); pred_sffs_23_flow_fd = categorical(cat_f_fd);
tbl_f_23_fd = crosstab(FD_cat, pred_sffs_23_flow_fd);
k_f_23_fd=kappaDLMmod(tbl_f_23_fd)

% this should, and does perform worse than HS (in kappa assessment)
%k_f_23_fd = 0.3662
%k_f_23    = 0.5684

%(and also in two-class ROC assessment)
labels = (FD_cat=='1')*1;
x = logit(1-(pihatB23_f_fd(:,1)));

% trim out 2
I = FD_cat=='2';
x(I) = [];
labels(I)=[];
showfigures = 1; figure(101);
[thresX_fd,AUC_fd,SEM,p,posclass,SensSpec_fd]=ROCAUCSEM(labels,x,showfigures);
% AUC_fd = 0.9318
% AUC    = 0.9859

%saveas(figure(101), 'ROC_HSvsFD','png');



%%
% get upper and lower limits on the features we like
theData = Amatrix2_flow(:,topFtrsfromSFFS(:,1));
[minval_flow, ~] = min(theData,[],1);
[maxval_flow, ~] = max(theData,[],1);

theData = Amatrix2_pnasal(:,topFtrsfromSFFS(:,1));
[minval_pnasal, ~] = min(theData,[],1);
[maxval_pnasal, ~] = max(theData,[],1);

FLcert = table();
FLcert.betas = B23;
FLcert.ftrNum = [NaN;NaN;[topFtrsfromSFFS(:,1)]];
FLcert.ftrName = [{'InterceptOne';'InterceptTwo'}; featuresWeLike];
FLcert.lowerFlow = [NaN;NaN;minval_flow'];
FLcert.upperFlow = [NaN;NaN;maxval_flow'];
FLcert.lower = [NaN;NaN;minval_pnasal'];
FLcert.upper = [NaN;NaN;maxval_pnasal'];

save('FLcertModel.mat','FLcert');

%% make a final model, using just the top 5 ftrs
[B5]=mnrfit(Amatrix2_flow(:,topFtrsfromSFFS(1:5,1)),HS_cat, 'Model', 'ordinal', 'link','logit');
[pihatB5] = mnrval(B5,Amatrix2_flow(:,topFtrsfromSFFS(1:5,1)), 'Model', 'ordinal');
[~,cat] = max(pihatB5(:,:),[],2); pred_sffs_5 = categorical(cat);

tbl_f_5 = crosstab(HS_cat, pred_sffs_5);
k_f_5=kappaDLMmod(tbl_f_5)
% kappa 
% 6 features 0.5556
% 5 features 0.5530
% 4 features 0.5423
% 3 features 0.5159
% 2 features 0.4849
% 1 feature 0.4076



%% do proportions, both on cross-validated data and final pred data
% get per pt level data 
doFlow = 0; % don't change here, changes during processing
pts = unique(PtData_flow.PT);
HSPropsPerPt = nan(54,3);
HSPropsPerPtPnasal = nan(54,3);
FinalMdl_kappa_flow = nan(length(pts),1);
FinalMdl_kappa_pnasal = nan(length(pts),1);

options={'Flow','Pnasal'};
ptprocessing=NaN(54,2);
for outerloop=1:length(options)
    str = ['Outerloop = ', num2str(outerloop)]; disp(str);
    DataIn = options{outerloop}; 
    switch DataIn
        case 'Flow'
            PredPropsPerPt_Flow_L1O = nan(length(pts),3);
            PredPropsPerPt_Flow_FinalMdl = nan(length(pts),3);
            FinalMdl_kappa_flow = nan(length(pts),1);
            L1O = predyL1Ocat_array_flow;
            FinalMdl = pred_sffs_23_flow;                
            GS_test=HS_cat;
            PtD = PtData_flow;
            doFlow = 1;
        case 'Pnasal'
            PredPropsPerPt_Pnasal_L1O = nan(length(pts),3);
            PredPropsPerPt_Pnasal_FinalMdl = nan(length(pts),3);
            FinalMdl_kappa_Pnasal = nan(length(pts),1);
            L1O = predyL1Ocat_array_pnasal(pnasalbbs);
            FinalMdl = pred_sffs_23_pnasal;
            GS_test=HS_cat(pnasalbbs);
            PtD = PtData_pnasal(pnasalbbs,:);
            doFlow = 0;
    end
    
    % loop through all pts 
    for pt=1:54
        disp(['Patient ', num2str((pt))]);
        %isubj = PtD.PT==pt;
        isubj = PtD.PT==pt & PtD.Hypnog<4;
        str = ['Pt ', num2str(pt), ' has ', num2str(nnz(isubj)), ' breaths']; disp(str);

        % Do L1O pred first
        label_L1O = L1O(isubj);
        if isempty(label_L1O)
            str = ['No L1O data for pt ', num2str(pt)]; disp(str);
            continue
        end
        if iscategorical(label_L1O)
            lt = array2table(label_L1O,'VariableNames',{'Var1'});
            l = summary(lt);
            pred_cats = l.Var1.Counts';
        else 
            pred_cats = [nnz(label_L1O==1) ...
            nnz(label_L1O==2)...
            nnz(label_L1O==3)];
        end
        pred_props_L1O = pred_cats./sum(pred_cats);
        
        
        % Do FinalMdl pred second
        label_FinalMdl = FinalMdl(isubj);
        if isempty(label_FinalMdl)
            str = ['No final mdl data for pt ', num2str(pt)]; disp(str);
            continue
        end
        if iscategorical(label_FinalMdl)
            lt = array2table(label_FinalMdl,'VariableNames',{'Var1'});
            l = summary(lt);
            pred_cats = l.Var1.Counts';
        else 
            pred_cats = [nnz(label_FinalMdl==1) ...
            nnz(label_FinalMdl==2)...
            nnz(label_FinalMdl==3)];
        end
        pred_props_mdl = pred_cats./sum(pred_cats);
       
        % Do GS third (last)       
        if iscategorical(GS_test)
            gt = array2table(GS_test(isubj),'VariableNames',{'Var1'});
            g = summary(gt);
            GS_cats = g.Var1.Counts';
        else
            GS_testpt = GS_test(isubj);
            GS_cats = [nnz(GS_testpt==1) ...
                nnz(GS_testpt==2)...
                nnz(GS_testpt==3)];
        end
        GS_props = GS_cats./sum(GS_cats);
        
        % store values
        if doFlow
            PredPropsPerPt_Flow_L1O(pt,:) = pred_props_L1O;
            PredPropsPerPt_Flow_FinalMdl(pt,:) = pred_props_mdl;
            HSPropsPerPt(pt,:) = GS_props;
            
            % get kappa's per pt on final mdl as well
            tbl_temp = crosstab(GS_test(isubj), label_FinalMdl);
            FinalMdl_kappa_flow(pt)=kappaDLMmod(tbl_temp);        

        else
            PredPropsPerPt_Pnasal_L1O(pt,:) = pred_props_L1O;   
            PredPropsPerPt_Pnasal_FinalMdl(pt,:) = pred_props_mdl;
            HSPropsPerPtPnasal(pt,:) = GS_props;
            
            % get kappa's per pt on final mdl as well
            tbl_temp = crosstab(GS_test(isubj), label_FinalMdl);
            FinalMdl_kappa_pnasal(pt)=kappaDLMmod(tbl_temp);   
            
        end
    end
end



    
%% plot of training and test kappa through each step in cross-validation
indivKappa = [training_kappa, test_kappa_pnasal, test_kappa_flow, FinalMdl_kappa_flow];
indivKappa = indivKappa(settings.Flow_list,:);

figure(1); clf(figure(1)); fig=gcf; fig.Color=[1 1 1]; 
stairs(indivKappa, 'LineWidth', 2, 'LineStyle', '-.'); box off;
legend('Training','Test-Pnasal','Test-Flow','FinalMDL', 'location','SouthWest');
xlabel('Step'); ylabel('Kappa'); %title('Kappa values through cross-validation');
if savefigs
    saveas(fig, ['Figures\KappaValuesThroughCrossValidation_flow'], 'png');
end
    

    
%% do Pt proportions figure, 2x2 subplots, CFL, PFL ann NFL, HS Vs Pred
sigdig = 2;
options={'Flow','Pnasal'};
for outerloop=1:length(options)
    DataIn = options{outerloop}; %'FDm'; % 'FDm' or 'FDp'
    switch DataIn
        case 'Flow'
            figure(25);clf(figure(25)); fig=gcf;
            fig.Units='inches'; fig.Position=[1 0.15 13.5 10];
            fig.Color=[1 1 1];
            PredData = PredPropsPerPt_Flow_L1O;
            HSProps = HSPropsPerPt;
            
        case 'Pnasal'
            figure(26);clf(figure(26)); fig=gcf;
            fig.Units='inches'; fig.Position=[1 0.15 13.5 10];
            fig.Color=[1 1 1];
            PredData = PredPropsPerPt_Pnasal_L1O;
            HSProps = HSPropsPerPtPnasal;
    end
    
    subplot(2,2,1); % HS_CFL Vs FDp_CFL
    [CFL_Rsq,Pvalue,Slope,Intercept,Xmodel_sem, CFL_r] = ...
        plotregressionwithSEM_dlm(PredData(:,1).*100,HSProps(:,1).*100);
    %[Rsq,p,MSE,beta,ypred,r,Rsq_dm] = glmfitFast_dlm(HSPropsPerPt(:,1),PredData(:,1),ones(length(PredData(:,1))),1);
    hold on
    %scatter(HSPropsPerPt(:,1).*100, ypred.*100, 'marker', 'x');
    text(80, 10, {['r = ', num2str(CFL_r,sigdig)];['R^2 = ', num2str(CFL_Rsq,sigdig)]});
    xlim([-5 100]); ylim([-5 100]);
    ylabel('Human Scoring (%)');
    xlabel('Model Predicted (%)');
    title('Proportion of breaths Certain-Flow-Limited');
    axis square
    
    subplot(2,2,2); % HS_PFL Vs FDp_PFL
    [PFL_Rsq,Pvalue,Slope,Intercept,Xmodel_sem, PFL_r] = ...
    plotregressionwithSEM_dlm(PredData(:,2).*100,HSProps(:,2).*100);
    hold on
    text(50, 10, {['r = ', num2str(PFL_r,sigdig)];['R^2 = ', num2str(PFL_Rsq,sigdig)]});
    xlim([-5 100]); ylim([-5 100]);
    ylabel('Human Scoring - proportion PFL');
    xlabel('Predicted - proportion PFL');
    title('Proportion of breaths Possible-Flow-Limited');
    axis square
    
    subplot(2,2,3); % HS_NFL Vs FDp_NFL
    [NFL_Rsq,Pvalue,Slope,Intercept,Xmodel_sem, NFL_r] = ...
    plotregressionwithSEM_dlm(PredData(:,3).*100,HSProps(:,3).*100);
    hold on
    text(50, 10, {['r = ', num2str(NFL_r,sigdig)];['R^2 = ', num2str(NFL_Rsq,sigdig)]});
    xlim([-5 100]); ylim([-5 100]);
    ylabel('Human Scoring (%)');
    xlabel('Model Predicted (%)');
    title('Proportion of breaths Non-Flow-Limited');
    axis square
    
    subplot(2,2,4); % HS_CFLPFL Vs FDp_CFLPFL
    [CandPFL_Rsq,Pvalue,Slope,Intercept,Xmodel_sem, CandPFL_r] = ...
    plotregressionwithSEM_dlm(sum(PredData(:,1:2),2).*100,sum(HSProps(:,1:2),2).*100);
    hold on
    text(80, 10, {['r = ', num2str(CandPFL_r,sigdig)];['R^2 = ', num2str(CandPFL_Rsq,sigdig)]});
    xlim([-5 100]); ylim([-5 100]);
    ylabel('Human Scoring - proportion CFL and PFL');
    xlabel('Predicted - proportion CFL and PFL');
    title('Proportion of breaths Certain and Possible-Flow-Limited');
    axis square
  
    supstr=(['Proportions of breaths in each category, cross-validated, per patient, using ', DataIn]);
    suptitle(supstr);
    if savefigs
        saveas(fig, ['Figures\PerPtProportions_OLR_FlowShapes_SleepBB',DataIn], 'png');
    end
end



%% Pt proportions figures for paper, 2 subplots, CFL, PFL and NFL, HS Vs Pred
mksize = 15; %15 default
mkcol = [0 0 0];
TickFntSz = 12;
LabelFntSz = 16;
TitleFntSz = 16;
sigdig = '%.2f';

L1Odata = 1;

options={'Flow','Pnasal'};
for outerloop=1:length(options)
    DataIn = options{outerloop}; %'FDm'; % 'FDm' or 'FDp'
    switch DataIn
        case 'Flow'
            figure(25);clf(figure(25)); fig=gcf;
            fig.Units='inches'; fig.Position=[1 0.15 10 5.5];
            fig.Color=[1 1 1];
            if L1Odata % L1O pred
                PredData = PredPropsPerPt_Flow_L1O;
            else % final model pred
                PredData = PredPropsPerPt_Flow_FinalMdl;
            end
            HSProps = HSPropsPerPt;
            
        case 'Pnasal'
            figure(26);clf(figure(26)); fig=gcf;
            fig.Units='inches'; fig.Position=[1 0.15 10 5.5];
            fig.Color=[1 1 1];
            if L1Odata % L1O pred
                PredData = PredPropsPerPt_Pnasal_L1O;
            else % final model pred
                PredData = PredPropsPerPt_Pnasal_FinalMdl;
            end
            HSProps = HSPropsPerPtPnasal;
    end
    
    subplot(1,2,1); % HS_CFL Vs FDp_CFL
    ax1 = gca; set(ax1,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz); hold on;
    [CFL_Rsq,Pvalue,Slope,Intercept,Xmodel_sem, CFL_r] = ...
        plotregressionwithSEM_dlm(PredData(:,1).*100,HSProps(:,1).*100);
    % get mean squared error    
    [Rsq,p,MSE,beta,ypred,r,Rsq_dm] = glmfitFast_dlm(HSProps(:,1).*100,PredData(:,1).*100,ones(length(PredData(:,1)),1),1);
    
    hold on
    %scatter(HSPropsPerPt(:,1).*100, ypred.*100, 'marker', 'x');
    %text(80, 10, {['r = ', num2str(CFL_r,sigdig)];['R^2 = ', num2str(CFL_Rsq,sigdig)]});
    %text(80, 10, {['R^2 = ', num2str(CFL_Rsq,sigdig)]});
    text(70, 10, {['RMSE = ', num2str(MSE.^ 0.5,sigdig)];['R^2 = ', num2str(CFL_Rsq,sigdig)]});
    
    xlim([-5 100]); ylim([-5 100]);
    ax1.XTick=[0:20:100]; ax1.YTick=[0:20:100];
    ylabel('Visual Scoring (%)', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
    xlabel('Flow Shape (%)', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
    title('Frequency of Certain FL',  'FontSize', TitleFntSz, 'FontName', 'Arial Narrow');
    axis square

    if 1
        subplot(1,2,2); % HS_CFLPFL Vs FDp_CFLPFL
        ax2 = gca; set(ax2,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz); hold on;
        [CandPFL_Rsq,Pvalue,Slope,Intercept,Xmodel_sem, CandPFL_r] = ...
            plotregressionwithSEM_dlm(sum(PredData(:,1:2),2).*100,sum(HSProps(:,1:2),2).*100);
        % get mean squared error    
        [Rsq,p,MSE,beta,ypred,r,Rsq_dm] = glmfitFast_dlm(HSProps(:,3).*100,PredData(:,3).*100,ones(length(PredData(:,3)),1),1);
    
        hold on
        %text(80, 10, {['r = ', num2str(CandPFL_r,sigdig)];['R^2 = ', num2str(CandPFL_Rsq,sigdig)]});
        %text(80, 10, {['R^2 = ', num2str(CandPFL_Rsq,sigdig)]});
        text(70, 10, {['RMSE = ', num2str(MSE.^ 0.5,sigdig)];['R^2 = ', num2str(CandPFL_Rsq,sigdig)]});
        xlim([-5 100]); ylim([-5 100]);
        ax2.XTick=[0:20:100];ax2.YTick=[0:20:100];
        ylabel('Visual Scoring (%)', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
        xlabel('Flow Shape (%)', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
        title('Frequency of Possible or Certain FL', 'FontSize', TitleFntSz, 'FontName', 'Arial Narrow');
        axis square
    else
        subplot(1,2,2); % HS_NFL Vs FDp_NFL
        ax2 = gca; set(ax2,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz); hold on;
        [NFL_Rsq,Pvalue,Slope,Intercept,Xmodel_sem, NFL_r] = ...
            plotregressionwithSEM_dlm(PredData(:,3).*100,HSProps(:,3).*100);
        hold on
        %text(80, 10, {['r = ', num2str(NFL_r,sigdig)];['R^2 = ', num2str(NFL_Rsq,sigdig)]});
        text(80, 10, {['R^2 = ', num2str(NFL_Rsq,sigdig)]});
        xlim([-5 100]); ylim([-5 100]);
        ax2.XTick=[0:20:100];ax2.YTick=[0:20:100];
        ylabel('Visual Scoring (%)', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
        xlabel('Flow Shape (%)', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
        title('Frequency of Possible or Certain FL', 'FontSize', TitleFntSz, 'FontName', 'Arial Narrow');
        axis square
    end

%     supstr=(['Proportions of breaths in each category, cross-validated, per patient, using ', DataIn]);
%     suptitle(supstr);
    if savefigs
        if L1Odata
            saveas(fig, ['Figures\ForPaper_PerPtProportions_OLR_FlowShapes_SleepBB_CFLandNFL_L1O',DataIn], 'png');            
        else
            saveas(fig, ['Figures\ForPaper_PerPtProportions_OLR_FlowShapes_SleepBB_CFLandNFL_FinalMDL',DataIn], 'png');
        end
    end
end

%%





%% GetAHI
% use getAHI_postprocessing function
AHI_perPT
AHI_perPT_table

%save('AHI_info.mat', 'AHI_perPT','AHI_perPT_table','AHIdata');

%% AHI processing
% FlowData = [HSPropsPerPt, PredPropsPerPt_Flow_L1O,PredPropsPerPt_Flow_FinalMdl];
% PnasalData = [HSPropsPerPtPnasal, PredPropsPerPt_Pnasal_L1O,PredPropsPerPt_Pnasal_FinalMdl];
% AllData = [AHI_perPT, FlowData, PnasalData];

AHIall = AHI_perPT(:,1);
AHIflow = AHI_perPT(settings.Flow_list,1);
AHIpnasal = AHI_perPT(settings.Pnasal_list,1);

%[muHat,sigmaHat]=normfit(AHIflow)
%[muHat,sigmaHat]=normfit(AHIpnasal)

[Ba, Ia] = sort(AHIall, 'descend');
[Bf, If] = sort(AHIflow, 'descend');
[Bp, Ip] = sort(AHIpnasal, 'descend');

%% stacked bar chart for each patient
PredData = PredPropsPerPt_Flow_L1O;
HSProps = HSPropsPerPt;
PredDataFinalMdl = PredPropsPerPt_Flow_FinalMdl;

d=0;
figure(1); clf(figure(1));
fig = gcf; fig.Color=[1 1 1];
fig.Units = 'Inches'; fig.Position = [10 3.5 11.5 9.5];
rows=7; cols=6;

for i=1:54 
    subj = Ia(i);
    % loop through all subjects
    if ismember(subj, settings.Flow_list)

    d=d+1;
    Isubj=(PtData_flow.PT==subj);  % Isubj is the logical index into current subject (flow)
    %nnz(Isubj)
       
    GS_proportions = HSProps(subj,:); 
    pred_proportions = PredData(subj,:);
    FinalMdl_proportion = PredDataFinalMdl(subj,:);
    
    subplot(rows, cols, d);
    axis off;
    bar([1, 2, 3], [GS_proportions; pred_proportions; FinalMdl_proportion], 'stacked'); hold on;
    xticklabels({'HS','Pred', 'Mdl'}); yticklabels({});
    titlestr = ['AHI = ', num2str(AHIall(subj),2)]; title(titlestr);
    box off
    end
end

if savefigs
    str=['Figures\StackedBar_flow'];
    saveas(fig, str, 'png');
    %print(fig, str, '-dtiff', '-r300');
end

%% stacked bar chart for each patient
PredDataP = PredPropsPerPt_Pnasal_L1O;
HSPropsP = HSPropsPerPtPnasal;
PredDataFinalMdl = PredPropsPerPt_Pnasal_FinalMdl;

d=0;
figure(2); clf(figure(2));
fig = gcf; fig.Color=[1 1 1];
fig.Units = 'Inches'; fig.Position = [9 6.5 11.5 3.8];
rows=3; cols=6;

for i=1:54 
    subj = Ia(i);
    % loop through all subjects
    if ismember(subj, settings.Pnasal_list)

    d=d+1;
    Isubj=(PtData_flow.PT==subj);  % Isubj is the logical index into current subject (flow)
    %nnz(Isubj)
       
    GS_proportions = HSPropsP(subj,:); 
    pred_proportions = PredDataP(subj,:);
    FinalMdl_proportion = PredDataFinalMdl(subj,:);
    
    subplot(rows, cols, d);
    axis off;
    bar([1, 2, 3], [GS_proportions; pred_proportions; FinalMdl_proportion], 'stacked'); hold on;
    xticklabels({'HS','Pred', 'Mdl'}); yticklabels({});
    titlestr = ['AHI = ', num2str(AHIall(subj),2)]; title(titlestr);
    box off
    end
end

if savefigs
    str=['Figures\StackedBar_Pnasal'];
    saveas(fig, str, 'png');
    %print(fig, str, '-dtiff', '-r300');
end

%% XY scatter of AHI vs possible and certain FL

% AHI on X
% Visual scoring on Y 
%   (a) freq of certain, and (b) freq of possible or certain

mksize = 15; 
mkcol = [0 0 0];
TickFntSz = 12;
LabelFntSz = 16;
TitleFntSz = 16;
sigdig = '%.2f';
options={'Flow','Pnasal'};
for outerloop=1:length(options)
    DataIn = options{outerloop}; %'FDm'; % 'FDm' or 'FDp'
    switch DataIn
        case 'Flow'
            figure(25);clf(figure(25)); fig=gcf;
            fig.Units='inches'; fig.Position=[1 0.15 10 5.5];
            fig.Color=[1 1 1];
            PredData = AHIall;
            HSProps = HSPropsPerPt;
            
        case 'Pnasal'
            figure(26);clf(figure(26)); fig=gcf;
            fig.Units='inches'; fig.Position=[1 0.15 10 5.5];
            fig.Color=[1 1 1];
            PredData = AHIall;
            HSProps = HSPropsPerPtPnasal;
    end
    
    subplot(1,2,1); % HS_CFL Vs FDp_CFL
    ax1 = gca; set(ax1,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz); hold on;
    [CFL_Rsq,Pvalue,Slope,Intercept,Xmodel_sem, CFL_r] = ...
        plotregressionwithSEM_dlm(PredData,HSProps(:,1).*100);
    % get mean squared error    
    [Rsq,p,MSE,beta,ypred,r,Rsq_dm] = glmfitFast_dlm(HSProps(:,1).*100,PredData,ones(length(PredData),1),1);
    
    hold on
    %scatter(HSPropsPerPt(:,1).*100, ypred.*100, 'marker', 'x');
    %text(80, 10, {['r = ', num2str(CFL_r,sigdig)];['R^2 = ', num2str(CFL_Rsq,sigdig)]});
    %text(80, 10, {['R^2 = ', num2str(CFL_Rsq,sigdig)]});
    text(70, 10, {['RMSE = ', num2str(MSE.^ 0.5,sigdig)];['R^2 = ', num2str(CFL_Rsq,sigdig)]});
    
    xlim([-5 100]); ylim([-5 100]);
    ax1.XTick=[0:20:100]; ax1.YTick=[0:20:100];
    ylabel('Visual Scoring (%)', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
    xlabel('AHI (evts/hour)', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
    title('Frequency of Certain FL',  'FontSize', TitleFntSz, 'FontName', 'Arial Narrow');
    axis square

   
        subplot(1,2,2); % HS_CFLPFL Vs FDp_CFLPFL
        ax2 = gca; set(ax2,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz); hold on;
        [CandPFL_Rsq,Pvalue,Slope,Intercept,Xmodel_sem, CandPFL_r] = ...
            plotregressionwithSEM_dlm(PredData,sum(HSProps(:,1:2),2).*100);
        % get mean squared error    
        [Rsq,p,MSE,beta,ypred,r,Rsq_dm] = glmfitFast_dlm(HSProps(:,3).*100,PredData,ones(length(PredData),1),1);
    
        hold on
        %text(80, 10, {['r = ', num2str(CandPFL_r,sigdig)];['R^2 = ', num2str(CandPFL_Rsq,sigdig)]});
        %text(80, 10, {['R^2 = ', num2str(CandPFL_Rsq,sigdig)]});
        text(70, 10, {['RMSE = ', num2str(MSE.^ 0.5,sigdig)];['R^2 = ', num2str(CandPFL_Rsq,sigdig)]});
        xlim([-5 100]); ylim([-5 100]);
        ax2.XTick=[0:20:100];ax2.YTick=[0:20:100];
        ylabel('Visual Scoring (%)', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
        xlabel('AHI (evts/hour)', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
        title('Frequency of Possible or Certain FL', 'FontSize', TitleFntSz, 'FontName', 'Arial Narrow');
        axis square
   

    if savefigs

            saveas(fig, ['Figures\ForPaper_PerPtProportions_AHI_SleepBB_CFLandNFL_',DataIn], 'png');            
       
    end
end





%% Ternary plots - Figure 502 - difference between Manual and Model
% je suis plus folle que toi
figure(502); clf(figure(502)); fig = gcf;
ax=gca;
fig.Units='Inches';
fig.Position = [21    2    6    5.5];
axmin = 0;
axmax = 1.0;

h1=ternplot_mod(GS_prop(:,1),GS_prop(:,2),GS_prop(:,3), axmin, axmax, 'ro', 'markersize', 10); hold on;
h2=ternplot_mod(pred_prop(:,1),pred_prop(:,2),pred_prop(:,3), axmin, axmax, 'ks', 'markersize', 10);

if 0 % if not already as a fraction of one, then do this step
    [GS_prop(:,1),GS_prop(:,2),GS_prop(:,3)] = fractions(GS_prop(:,1),GS_prop(:,2),GS_prop(:,3));
    [pred_prop(:,1),pred_prop(:,2),pred_prop(:,3)] = fractions(pred_prop(:,1),pred_prop(:,2),pred_prop(:,3));
end

[GSX, GSY] = terncoords_mod(GS_prop(:,1),GS_prop(:,2),GS_prop(:,3), axmin, axmax);
[predX, predY] = terncoords_mod(pred_prop(:,1),pred_prop(:,2),pred_prop(:,3), axmin, axmax);

% then do line connecting Ref and SHR marker for each patient
% make some gray by adding color option to plot line below
gray = [0.25 0.25 0.25]; %  , 'color', gray
for j = 1:numel(GSX)
    plot([GSX(j) predX(j)],[GSY(j) predY(j)], 'k-', 'LineWidth', 1.25, 'color', gray);
end


handle_label = ternlabel('Certain-FL','Possible-FL','Non-FL');
handle_label(1).FontSize = 14;
handle_label(1).Color = [0,1,0, 0.5]; % Ti green
handle_label(2).FontSize = 14;
handle_label(2).Color = [0,0,1,0.5]; %Te blue
handle_label(3).FontSize = 14;
handle_label(3).Color = [1,0,1, 0.5]; % Ttran magenta

%legend('off');
handle_legend = legend(ax,[h1,h2], {'Manual','Model'});
legend('BoxOff');
handle_legend.FontSize = 12;
if 0
    for j=1:numel(GSX)
        %text(RefX(i)-0.015, RefY(i)-0.005, num2str(i));
        %text(RefX(i)-0.015, RefY(i)-0.005, num2str(AHI(i)));
        %text(RefX(i)-0.02, RefY(i)-0.01, num2str(PtProcessNum2(i)), 'color', 'b');
    end
end

if savefigs
    str=['Figures\Ternary_ManualVsModel'];
    print(fig, str, '-dtiff', '-r300');
end
if closefigs; close(fig); end


%% Ternary plots - Figure 503 - How does AHI look
figure(503)
clf(figure(503));
fig = gcf;
ax=gca;
fig.Units='Inches';
fig.Position = [21    2    6    5.5];
axmin = 0;
axmax = 1.0;

h0=ternplot_mod(GS_prop(:,1),GS_prop(:,2),GS_prop(:,3), axmin, axmax, 'k.', 'markersize', 1); hold on;
%h2=ternplot_mod(pred_prop(:,1),pred_prop(:,2),pred_prop(:,3), axmin, axmax, 'ro', 'markersize', 10);

[GSX, GSY] = terncoords_mod(GS_prop(:,1),GS_prop(:,2),GS_prop(:,3), axmin, axmax);
%[predX, predY] = terncoords_mod(pred_prop(:,1),pred_prop(:,2),pred_prop(:,3), axmin, axmax);

handle_label = ternlabel('Certain-FL','Possible-FL','Non-FL');
handle_label(1).FontSize = 14;
handle_label(1).Color = [0,1,0, 0.5]; % Ti green
handle_label(2).FontSize = 14;
handle_label(2).Color = [0,0,1,0.5]; %Te blue
handle_label(3).FontSize = 14;
handle_label(3).Color = [1,0,1, 0.5]; % Ttran magenta

if 1
    for j=1:numel(GSX)
        mksize = max(1,round(AHI_forPlot(j)/1.5));
        h1 = plot(GSX(j), GSY(j),'r.', 'markersize',mksize);
        if mksize<6
            h10 = plot(GSX(j), GSY(j),'bs', 'markersize',10);
        end
    end
end

legend('off');
handle_legend = legend(ax,[h1,h10], {'Size shows AHI', 'small values'});
legend('BoxOff');
handle_legend.FontSize = 12;

if savefigs
    str=['Figures\Ternary_ManualAndAHI'];
    print(fig, str, '-dtiff', '-r300');
end
if closefigs; close(fig); end


%% Ternary plots - Figure 503 - How does FD look
fd_mod_array = normalise_0to1(FD_array);
figure(504)
clf(figure(504));
fig = gcf;
ax=gca;
fig.Units='Inches';
fig.Position = [21    2    6    5.5];
axmin = 0;
axmax = 1.0;

h0=ternplot_mod(GS_prop(:,1),GS_prop(:,2),GS_prop(:,3), axmin, axmax, 'k.', 'markersize', 1); hold on;
%h2=ternplot_mod(pred_prop(:,1),pred_prop(:,2),pred_prop(:,3), axmin, axmax, 'ro', 'markersize', 10);

[GSX, GSY] = terncoords_mod(GS_prop(:,1),GS_prop(:,2),GS_prop(:,3), axmin, axmax);
%[predX, predY] = terncoords_mod(pred_prop(:,1),pred_prop(:,2),pred_prop(:,3), axmin, axmax);

handle_label = ternlabel('Certain-FL','Possible-FL','Non-FL');
handle_label(1).FontSize = 14;
handle_label(1).Color = [0,1,0, 0.5]; % Ti green
handle_label(2).FontSize = 14;
handle_label(2).Color = [0,0,1,0.5]; %Te blue
handle_label(3).FontSize = 14;
handle_label(3).Color = [1,0,1, 0.5]; % Ttran magenta

if 1
    for j=1:numel(GSX)
        mksize = max(1,round(fd_mod_array(j)*45));
        h1 = plot(GSX(j), GSY(j),'r.', 'markersize',mksize);
        if mksize<6
            h10 = plot(GSX(j), GSY(j),'bs', 'markersize',10);
        end
    end
end

%legend('off');
handle_legend = legend(ax,[h1,h10], {'Size shows FD', 'small values'});
legend('BoxOff');
handle_legend.FontSize = 12;

if savefigs
    str=['Figures\Ternary_ManualAndFD'];
    print(fig, str, '-dtiff', '-r300');
end
if closefigs; close(fig); end


%% check AHI data looks about right
figure();
hist(AHI_perPT(:,1))

%% check FD data looks about right
figure();
hist(FD_array);


%% FlowDrivemeasured vs Human scoring
% index into human score column for certain, possible and non-FL
CFL_i=HS==1; % certain nnz(CFL_i)
PFL_i=HS==2; % possible nnz(PFL_i)
NFL_i=HS==3;   % non-FL nnz(NFL_i)

FDm_CFL = PtData_flow.g_Edi_Adj(CFL_i);  
FDm_PFL = PtData_flow.g_Edi_Adj(PFL_i); 
FDm_NFL = PtData_flow.g_Edi_Adj(NFL_i); 

edges = [0:0.025:1.5];
figure(32); clf(figure(32)); fig=gcf;
fig.Units='inches'; fig.Position=[2.5 3.5 8 5.5];
fig.Color=[1 1 1];
ax = gca;
histogram(FDm_CFL,'DisplayStyle', 'stairs', 'LineWidth', 2); hold on;
histogram(FDm_PFL,'DisplayStyle', 'stairs', 'LineWidth', 2);
histogram(FDm_NFL,'DisplayStyle', 'stairs', 'LineWidth', 2);

% histogram(FDm_CFL, edges); hold on;
% histogram(FDm_PFL, edges);
% histogram(FDm_NFL, edges);

xticks([0 0.3 0.6 0.9 1.2 1.5])
xticklabels({'0', '30', '60', '90', '120', '150'});
xlabel('\it{flow:drive_{measured}} (%)'); %title('CFL');
xlim([-0.00 1.55]);
set(ax,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', settings.TickFntSz);
legend('Certain-flow-limited','Possible-flow-limited','Non-flow-limited','Location','NorthWest');
if 0
    str = ['Figures\HistogramLine_FlowDrivePerHumanCategory'];
    %print(fig, str, '-dtiff', '-r300');
    saveas(fig, str, 'png');
    %savefig(str);
end


%% set up ROC data for uber figure
%DataOut(:,1) % col1 is flow:drive
%DataOut(:,2) % col2 is Manual Score

DataOut(:,1) = PtData_flow.g_Edi_Adj;
DataOut(:,2) = HS;

% index into certain- possible- and non-flow-limited
CFL_i = DataOut(:,2)==1; 
PFL_i = DataOut(:,2)==2; 
NFL_i = DataOut(:,2)==3; 

% CFL Vs NFL
ROC_Data = DataOut;
ROC_Data(PFL_i,:)=[]; % remove the PFL component
IndxToAdj = ROC_Data(:,2)==1; % adjustment, 1 to be 0
ROC_Data(IndxToAdj,2)=0;
IndxToAdj = ROC_Data(:,2)==3; % adjustment, 3 to be 1
ROC_Data(IndxToAdj,2)=1;
figure(100); clf(figure(100));
ax = gca; set(gca,'box','off','tickdir','out','fontname','arial narrow','fontsize',TickFntSz); hold on;
[fpr1,tpr1,thresholds1,auc1,optpt1] = perfcurve(ROC_Data(:,2)', ROC_Data(:,1)',1);
plot(fpr1, tpr1); 
text(optpt1(1),optpt1(2),['\leftarrow AUC = ',num2str(auc1,2)]);

% CFL Vs PFL
ROC_Data = DataOut;
ROC_Data(NFL_i,:)=[]; % remove the NFL component
IndxToAdj = ROC_Data(:,2)==1; % adjustment, 1 to be 0
ROC_Data(IndxToAdj,2)=0;
IndxToAdj = ROC_Data(:,2)==2; % adjustment, 2 to be 1
ROC_Data(IndxToAdj,2)=1;
%plotroc(ROC_Data(:,2)', ROC_Data(:,1)')
[fpr2,tpr2,thresholds2,auc2,optpt2] = perfcurve(ROC_Data(:,2)', ROC_Data(:,1)',1);
plot(fpr2, tpr2);
text(optpt2(1),optpt2(2),['\leftarrow AUC = ',num2str(auc2,2)]);

% PFL Vs NFL
ROC_Data = DataOut;
ROC_Data(CFL_i,:)=[]; % remove the CFL component
IndxToAdj = ROC_Data(:,2)==2; % adjustment, 2 to be 0
ROC_Data(IndxToAdj,2)=0;
IndxToAdj = ROC_Data(:,2)==3; % adjustment, 3 to be 3
ROC_Data(IndxToAdj,2)=1;
%plotroc(ROC_Data(:,2)', ROC_Data(:,1)')
[fpr3,tpr3,thresholds3,auc3,optpt3] = perfcurve(ROC_Data(:,2)', ROC_Data(:,1)',1);
plot(fpr3, tpr3);
text(optpt3(1),optpt3(2),['\leftarrow AUC = ',num2str(auc3,2)]);

axis 'square'
ax.YTick=[0:0.1:1];
ax.XTick=[0:0.1:1];
xlabel('False Positive Rate', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
ylabel('True Positive Rate', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');

text(-0.2, 1.0, 'B', 'FontSize', 20, 'FontName', 'Arial Bold');


%title('Receiver Operating Characteristic Plot');

legend('Certain- Vs Non-Flow-Limited', 'Certain- Vs Possible-Flow-Limited', 'Possible- Vs Non-Flow-Limited', 'Location','southeast');

[h_1,p_1,ci_1,stats_1] = ttest2(DataOut(CFL_i,1),DataOut(NFL_i,1)); % 
[h_2,p_2,ci_2,stats_2] = ttest2(DataOut(CFL_i,1),DataOut(PFL_i,1)); % 
[h_3,p_3,ci_3,stats_3] = ttest2(DataOut(PFL_i,1),DataOut(NFL_i,1)); % 

[auc1; auc2; auc3]




%% uber-figure
% histogram, with box plot below
% x axes bottom - flow:drive
% y axes right - number of breaths fro histogram
% ROC curves
% x axes top - false positve rate
% y axis left - true positive rate


figure(39); clf(figure(39)); fig = gcf;
fig.Units='inches';
fig.Position=[1 1 8.5 9.2];
fig.Color=[1 1 1];
darker = 0.5;

% top subplot, occupying top 2/3rds of space is histogram, with ROC overlay
subplot(2,1,1);
ax1 = gca; set(ax1,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz); hold on;
yyaxis right
histogram(ax1, FDm_PFL,'DisplayStyle', 'stairs', 'LineWidth', 2, 'EdgeColor', [0 1 0]*darker); hold on;
histogram(ax1, FDm_CFL,'DisplayStyle', 'stairs', 'LineWidth', 2, 'EdgeColor', [1 0 0]*darker); hold on;
histogram(ax1, FDm_NFL,'DisplayStyle', 'stairs', 'LineWidth', 2, 'EdgeColor', [0 0 1]*darker);
hold on;
ax1.XTick=[];
ylabel('Number of Breaths');


yyaxis left
hold on
ylabel('True Positive Rate');

ax1.XAxisLocation=('top');    
plot(fpr1, tpr1, 'color', [1 0 1]*darker, 'linewidth', 2, 'linestyle', '-');
plot(fpr2, tpr2, 'color', [1 1 0]*darker, 'linewidth', 2, 'linestyle', '-');
plot(fpr3, tpr3, 'color', [0 1 1]*darker, 'linewidth', 2, 'linestyle', '-');
text(optpt1(1),optpt1(2),['\leftarrow AUC = ',num2str(auc1,2)]);  
text(optpt2(1),optpt2(2),['\leftarrow AUC = ',num2str(auc2,2)]); 
text(optpt3(1),optpt3(2),['\leftarrow AUC = ',num2str(auc3,2)]);
xlabel('False Positive Rate');
ax1.XTick=[0:0.2:1];
ylim([0 1]);
xlim([0 1.55]);

% bottom subplot, last 1/3rd of space is horizontal boxplot
subplot(2,1,2);
ax2 = gca; set(ax2,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz); hold on;
clr = get(ax1,'colororder');
h = boxplot(PtData_flow.g_Edi_Adj,HS,'orientation','horizontal', 'label',{'','',''},'color','rgb');
set(h,'LineWidth',1.5); % or set(h,{'linew'},{2})
set(h(7,:),'MarkerSize',0.1)
xlim([0 1.55]);
box off
% add labels (rotated) to say CFL, PFL and NFL
ax2.YTickLabel={'CFL', 'PFL', 'NFL'};
ax2.XTickLabel={'0', '50','100','150'};
xlabel('\itFlow:drive_{measured} (%)');

ax1.Position = [0.1 0.38 0.8 0.55];
ax2.Position = [0.1 0.07 0.8 0.25];




if 0
    str = ['Figures\ForPaper_UberFigure'];
    %print(fig, str, '-dtiff', '-r300');
    saveas(fig, str, 'png');
    %savefig(str);
end

%% matlab plot colours
 colors = [0    0.4470    0.7410;...
        0.8500    0.3250    0.0980;...
        0.9290    0.6940    0.1250;...
        0.4940    0.1840    0.5560;...
        0.4660    0.6740    0.1880;...
        0.3010    0.7450    0.9330;...
        0.6350    0.0780    0.1840];



%% not-quite-so-uber uber-figure
% histogram, with box plot below
% x axes bottom - flow:drive
% y axes right - number of breaths fro histogram
% ROC AUC data plotted on stat bars above distribution

TickFntSz = 12;
LabelFntSz = 16;
sigdig = '%.2f';

savefigs = 1;

figure(39); clf(figure(39)); fig = gcf;
fig.Units='inches';
fig.Position=[1 0.3 8.5 10];
fig.Color=[1 1 1];
darker = 1;

% top subplot, occupying top 2/3rds of space is histogram, with ROC overlay
subplot(2,1,1);
ax1 = gca; set(ax1,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz); hold on;
h1 = histogram(ax1, FDm_CFL, 'FaceColor', colors(1,:)*darker,'EdgeAlpha', 0); hold on;
h2 = histogram(ax1, FDm_PFL, 'FaceColor', colors(2,:)*darker,'EdgeAlpha', 0); hold on;
h3 = histogram(ax1, FDm_NFL, 'FaceColor', colors(3,:)*darker,'EdgeAlpha', 0); hold on;
ylabel('Number of Breaths', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
xlim([0 1.55]);
%ax1.XTick=[];
xlabel('\itFlow:drive_{measured} (%)', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');

ctrs = [0.42 0.72 0.98];

text(ctrs(1), 1600, 'Certain FL', 'HorizontalAlignment', 'center');
text(ctrs(2), 2050, 'Possible FL', 'HorizontalAlignment', 'center');
text(ctrs(3), 2050, 'Normal', 'HorizontalAlignment', 'center');

%[h_1,p_1,ci_1,stats_1] = ttest2(FDm_CFL, FDm_PFL)
e1 = errorbar(ax1, ctrs(1)+(0.5.*(ctrs(2)-ctrs(1))),2530,0.5.*(ctrs(2)-ctrs(1)),'horizontal'); 
e1.Color = 'black'; e1.LineWidth = 1.5;
text(ctrs(1)+(0.5.*(ctrs(2)-ctrs(1))),2600,['AUC = ',num2str(auc2,2)], 'HorizontalAlignment', 'center');  

%[h_2,p_2,ci_2,stats_2] = ttest2(FDm_PFL, FDm_NFL)
e2 = errorbar(ax1, ctrs(2)+(0.5.*(ctrs(3)-ctrs(2))),2630,0.5.*(ctrs(3)-ctrs(2)),'horizontal');
e2.Color = 'black'; e2.LineWidth = 1.5;
text(ctrs(2)+(0.5.*(ctrs(3)-ctrs(2))),2700,['AUC = ',num2str(auc3,2)], 'HorizontalAlignment', 'center'); 

%[h_3,p_3,ci_3,stats_3] = ttest2(FDm_CFL, FDm_NFL)
e3 = errorbar(ax1, ctrs(1)+(0.5.*(ctrs(3)-ctrs(1))),2780,0.5.*(ctrs(3)-ctrs(1)),'horizontal');
e3.Color = 'black'; e3.LineWidth = 1.5;
text(ctrs(1)+(0.5.*(ctrs(3)-ctrs(1))),2850,['AUC = ',num2str(auc1,2)], 'HorizontalAlignment', 'center');  

% bottom subplot, last 1/3rd of space is horizontal boxplot
subplot(2,1,2);
ax2 = gca; set(ax2,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz); hold on;
clr = get(ax1,'colororder');
clr2=[colors(1,:); colors(2,:); colors(3,:)]*darker;
% h = boxplot(PtData_flow.g_Edi_Adj,HS,'orientation','horizontal', 'label',{'','',''},'color',clr2);
% set(h,'LineWidth',1.5); % or set(h,{'linew'},{2})
% set(h(7,:),'MarkerSize',0.1)
HS_1 = HS==1;
HS_2 = HS==2;
HS_3 = HS==3;

CFL_q = quantile(PtData_flow.g_Edi_Adj(HS_1),[0.25 0.5 0.75]);
PFL_q = quantile(PtData_flow.g_Edi_Adj(HS_2),[0.25 0.5 0.75]);
NFL_q = quantile(PtData_flow.g_Edi_Adj(HS_3),[0.25 0.5 0.75]);

% plot errorbars median
Y = [CFL_q(2); PFL_q(2); NFL_q(2)]; % medians
X = [1;2;3]; % group
ErrorU = [CFL_q(3); PFL_q(3); NFL_q(3)]; % 3rd quartile
ErrorL = [CFL_q(1); PFL_q(1); NFL_q(1)]; % 1st quartile
BarColor=[colors(1,:)*darker;colors(2,:)*darker;colors(3,:)*darker];
baseline=-10;
barwidth=0.8;
linewidth=0.01;
plotbarwitherrorsmedian_dlm(Y,X,ErrorU,ErrorL,BarColor,baseline,barwidth,linewidth)
ylim([0 1.55]);
ax2.YTick=[];
ax2.YAxis.Color=[1 1 1];

box off
% add labels (rotated) to say CFL, PFL and NFL
ax2.XTick=[1:1:3];
ax2.XTickLabel={'Certain FL', 'Possible FL', 'Normal'};

% yyaxis right
% ylabel('\itFlow:drive_{measured} (%)', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
% ylim([0 1.55])
% YTick=[0:0.5:1.5];

%ax2.XTickLabel={'0', '25', '50','75', '100','125','150'};
%xlabel('\itFlow:drive_{measured} (%)', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');

ax1.Position = [0.1 0.39 0.8 0.55];
ax2.Position = [0.1 0.07 0.8 0.25];

camroll(-90)

if savefigs
    str = ['Figures\ForPaper_NotQuiteSoUberFigure'];
    %print(fig, str, '-dtiff', '-r300');
    saveas(fig, str, 'png');
    %savefig(str);
end


