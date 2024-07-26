% in this code, pnasal is tested with betas trained in flow 

%% start
keyboard % this is just to stop the whole process running accidentally (i.e. pressing F5) 
close all
clear
clc

% goto line 800 for post-processing

if 0
    addpath('C:\Users\uqdmann\Dropbox\PUPbeta_git\');
    cd('C:\Users\uqdmann\Dropbox\PUPbeta_git\FlowAndEdiPesQAO');
    load('C:\Users\uqdmann\Dropbox\PUPbeta_git\FeatureSpaces\FlowDrive_only25Hz_FeatureSpace_AutoRef2_Edi_Clean_withPnasalMatched_LinRegWorkspace_TestNTrainPnasal_20180129.mat');
end

%% options
useLogitFtrs = 0;       % set as 1 to use ftrs selected from logit, 0 for all ftrs
SS_fastfit = 1;			% set as 1 to use SS fast glmfit method, 0 for matlab builtin
TransformTheData = 1;   % set as 1 to do tranforms, or 0 to use unadjusted data
addextratransform = 0; 	% set as 1 to do extra transforms
L1O_run = 1;            % set as 1 to run to do the L1O run, or 0 to skip run
ShowFigures = 0;        % set as 1 to show figures, or 0 to not show figures
                     
%datadir = 'C:\PSG_Data\FlowDrive\FeatureSpaces\';
%datadir = '..\FeatureSpaces\';  %
datadir = 'C:\PSG_Data\FlowDrive\Analyzed\ExclBadR\';

useRthreshold = 1; % set to one to enable, then goto line 530 to set R exclusion threshold
experimentnumber = '_exp14';

%%
if L1O_run
    %% Turn diary logging on
    HowSoonIsNow = datestr(datetime('now'),'yyyymmdd_HHMM');
    diaryfilename = ['L1O_LinReg_FD_TnTinPnasal_', HowSoonIsNow, '.txt'];
    diary(diaryfilename);
    diary on
    
    %% open files
    filename = [datadir, 'FlowDrive_only25Hz_ReNormalized_FeatureSpace_AutoRef2_Edi_Clean.mat']; datastr = 'FD Normalized edi 25Hz';
    pnasal_filename = [datadir, 'PnasalDrive_only25Hz_ReNormalized_FeatureSpace_AutoRef2_Edi_Clean.mat']; datastr = 'FlowandPnasalD edi 25Hz';

    try
        str=['Loading ' pnasal_filename]; disp(str);
        load(pnasal_filename,'Amatrix', 'PtData', 'FeatureNames');
        Amatrix_pnasal = Amatrix;
        PtData_pnasal = PtData;
        FeatureNames_pnasal = FeatureNames;
        str=['Loading ' filename]; disp(str);
        load(filename);
    catch me
        disp(me.getReport);
    end
    
    %% Just in case, look at NaN again
    allnanrows = sum(isnan(Amatrix),2)==size(Amatrix,2);
    if nnz(allnanrows)~=0
        str = ['Removing ', num2str(sum(allnanrows)), ' breaths that contain NaN''s']; disp(str);
        Amatrix(allnanrows,:)=[];
        PtData(allnanrows,:)=[];
        Fnan=sum(isnan(Amatrix)|isinf(Amatrix))/size(Amatrix,1);
        if nnz(Fnan)~=0
            disp('NaN''s or non-finite data remains - further investigation required');
            keyboard;
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
            keyboard;
        end
    else
        str = ['Zero NaN-breaths were found in Pnasal data']; disp(str);
    end
    
    %% use only a given set of features, determined from logistic regression
    % keep the original Amatrix (before selecting specific features, and transforms
    % keep the origianl FeatureNames variable before overwriting
    % same for pnasal, if included
    Original_FeatureNames = FeatureNames;
    Original_Amatrix = Amatrix;
    if useLogitFtrs
        %load('FtrsToUse','LabelsOrderedOpt');
        %load('FtrsToUse_Normalized','LabelsOrderedOpt');
        load('FtrsToUse_ReNormalized_2','LabelsOrderedOpt');
        
        FtrsToInclude = LabelsOrderedOpt;
        Ind = [];
        for i=1:length(FtrsToInclude)
            temp=find(startsWith(FeatureNames.Name,FtrsToInclude(i)));
            if ~isempty(temp)
                Ind(end+1)=temp;
            end
        end
        FeatureNames = FeatureNames(Ind,:);
        Amatrix = Amatrix(:,Ind);
        
        % do it again for pnasal
        Original_FeatureNames_pnasal = FeatureNames_pnasal;
        Original_Amatrix_pnasal = Amatrix_pnasal;
        Ind = [];
        for i=1:length(FtrsToInclude)
            temp=find(startsWith(FeatureNames_pnasal.Name,FtrsToInclude(i)));
            if ~isempty(temp)
                Ind(end+1)=temp;
            end
        end
        FeatureNames_pnasal = FeatureNames_pnasal(Ind,:);
        Amatrix_pnasal = Amatrix_pnasal(:,Ind);
    else
        % do a simple check, so that we at least start with the same features
        [~,ia,ib] = intersect(FeatureNames_pnasal.Name, FeatureNames.Name);
        FeatureNames_pnasal = FeatureNames_pnasal(ia,:);
        Amatrix_pnasal = Amatrix_pnasal(:,ia);
        FeatureNames = FeatureNames(ib,:);
        Amatrix = Amatrix(:,ib);
    end
    
    %% setup the data matrix to use.
    if TransformTheData
        [Amatrix2, Labels] = DoDataMatTransform(Amatrix, FeatureNames, addextratransform);
        [Amatrix2_pnasal, ~] = DoDataMatTransform(Amatrix_pnasal, FeatureNames_pnasal, addextratransform);
    else
        Amatrix2 = [Amatrix];
        Labels = FeatureNames.Name;
        Amatrix2_pnasal = Amatrix_pnasal;
    end
    
    %% Preparation (Do this carefully)
    Gtest_All = PtData.g_Edi_Adj;
    Gtest_All_pnasal = PtData_pnasal.g_Edi_Adj; % potentially unused, as we use Flow based Gtest
  
    colofones = ones(length(Gtest_All),1);
   
    
    %% Weights
    % ToDo: check method of weights
    clear Ndata
    maxG=1.5;
    dx=0.2;
    xbins=[0 0.3:dx:0.9 maxG];
    %xbins=[0 0.1:0.05:1.1 1.5];
    
    for i=1:length(xbins)-1
        Ix=Gtest_All>xbins(i)&Gtest_All<=xbins(i+1);
        Ndata(i)=sum(Ix);
    end
    %Ndata = Ndata.^0.5;
    weightsbins = 1./(Ndata);
    %weightsbins = [2 1 1 1 0.5];
    weightsbins = weightsbins/mean(weightsbins);
    %weightsbins = weightsbins/weightsbins(end);
    
    weights = NaN*Gtest_All;
    for i=1:length(xbins)-1
        Ix=Gtest_All>=xbins(i)&Gtest_All<=xbins(i+1);
        weights(Ix)=weightsbins(i);
    end
    weights = weights/nanmean(weights); % nnz(~isfinite(weights))
    
    useweights=1;
    if ~useweights %overwrite, use no weights
        weights = ones(length(weights),1);
    end
    
    %% Leave one subject out loop
    % lin reg leave one out
    % FD, 41 Subjects, 135k breaths, 50 ftrs, 2x transforms, takes 30 min to 1 hr
    
    % use glmfit to do generalized linear model regression
    % use fitglm to create a generalized linear regression model (logistic?)
    clear RvalueTrain Err ErrRms Nfeatures
    clear RsqTrain_array ErrTrain_array ErrRmsTrain_array
    labels_Step_Subj =[];
    alpha=0.05;
    MaxNfeatures=1;
    neverbreak=1;
    t_start_L1O = clock;
    pts_pnasal = unique(PtData.PT);
    %Pnasal_list = [3 5 8 9 18 22 33 34 35 36 38 43 46 47 50 53 54];
    Pnasal_list = unique(PtData_pnasal.PT);
    Labels_Complete = Labels; % save a backup of the complete Label set
    Fwd=0;
    VEVeup_All = [];
    VEVeup_pnasal_array = [];
    VEVeup_flow_array = [];
    Pnasal_summary = [];
    PnasalFlowBB_PTindex=[];
    PtData_matched = [];
    Gtest_matched = [];
    weights_matched = [];
    Amatrix2_pnasal_matched = [];
    Amatrix2_flow_matched = [];
    Amatrix2_p_matched = [];
    PtData_pnasal_matched = [];
    predyL1O_array_pnasal_matched = [];
    predyL1O_array_flow_matched = [];
    ChannelsList = {'Flow','Pnasal'};
    artdirectory = ['C:\PSG_Data\FlowDrive\SourceMat 20171123'];
    % read spreadsheet (options worksheet)
    [~,~,raw] = xlsread('C:\Users\uqdmann\Dropbox\PUPbeta_git\FlowAndEdiPesQAO\AnalyzeDataSpreadsheet.xlsx',1,'G3:G56');
    
    %% go through all patients once to compile the data
    for subj=1:54 
        %% loop through all subjects, getting aligned Pnasal and Flow breaths
        if ~ismember(subj, Pnasal_list)
            str=['No Pnasal for Pt ', num2str(subj)]; disp(str); continue
        end
        
        disp(' '); str=['Getting data for Pt ', num2str(subj)]; disp(str); tic;

        PtHasPnasal = 1;
        
        %% set up the full length Flow data
        Isubj=(PtData.PT==subj);  % Isubj is the logical index of the L1O patient
        
        % set up training data
        Gtest_train = Gtest_All;
        Gtest_train(Gtest_train>maxG)=maxG;
        Gtest_train(Isubj)=[];
        weights_train = weights;
        weights_train(Isubj)=[];
        weights_train=weights_train/nanmean(weights_train);
        colofones_train = colofones;
        colofones_train(Isubj)=[];
        Amatrix2_train = Amatrix2;
        Amatrix2_train(Isubj,:)=[];
        
        %remove all NAN from all rows
        nanx = sum(isnan(Amatrix2_train),2)>0;   % nnz(nanx)
        nany = isnan(Gtest_train);               % nnz(nany)
        if nnz(nanx)>0 || nnz(nany)>0; disp('NaNs found'); end
        Amatrix2_train(nanx|nany,:)=[];
        Gtest_train(nanx|nany)=[];
        weights_train(nanx|nany)=[];
        weights_train=weights_train/nanmean(weights_train);
        colofones_train(nanx|nany)=[];
        
        % set up test data
        Gtest_test = Gtest_All(Isubj);
        Gtest_test(Gtest_test>maxG)=maxG;
        weights_test = weights(Isubj);
        weights_test=weights_test/nanmean(weights_test);
        colofones_test = colofones(Isubj);
        Amatrix2_test = Amatrix2(Isubj,:);
        
        %remove all NAN from all rows
        nanx = sum(isnan(Amatrix2_test),2)>0;   % nnz(nanx)
        nany = isnan(Gtest_test);               % nnz(nany)
        if nnz(nanx)>0 || nnz(nany)>0; disp('NaNs found'); end
        Amatrix2_test(nanx|nany,:)=[];
        Gtest_test(nanx|nany)=[];
        weights_test(nanx|nany)=[];
        weights_test=weights_test/nanmean(weights_test);
        colofones_test(nanx|nany)=[];
        
        
        %% set up the pnasal based data
        % in doing this, we also produce a second set of flow data, but
        % these are matched with the pnasal breaths
        Isubj_pnasal=(PtData_pnasal.PT==subj);
        PnasalFlowBB_PTindex_pt=[];
        
        predyL1O_array_pnasal_pt = [];
        predyL1O_array_flow_pt = [];
        VEVeup_pnasal_array_pt =[];
        VEVeup_flow_array_pt =[];
        Gtest_test_flow = Gtest_test;
        weights_test_flow = weights_test;
        
        StarttimeSpike = 0; % reset
        fname = char(raw{subj});
        load(['C:\PSG_Data\FlowDrive\SourceMat 20171123\' fname], 'StarttimeSpike');
        
        PtData_flow_pt = PtData(Isubj,:);
        PtBB_flow_time = PtData_flow_pt.BB_time-StarttimeSpike;
        Amatrix2_flow_pt = Amatrix2(Isubj,:);
        str = [num2str(length(PtBB_flow_time)), ' breaths with Flow - before art removal']; disp(str);
        
        PtData_pnasal_pt = PtData_pnasal(Isubj_pnasal,:);
        PtBB_pnasal_time = PtData_pnasal_pt.BB_time-StarttimeSpike;
        Amatrix2_pnasal_pt = Amatrix2_pnasal(Isubj_pnasal,:);
        str = [num2str(length(PtBB_pnasal_time)), ' breaths with Pnasal - before art removal']; disp(str);
               
        % trim both flow and pnasal, keep only breaths that are non-artifact in both
        for j=1:length(ChannelsList)
            textfilename=[artdirectory '\' fname(1:end-4) '_' ChannelsList{j} '_art.txt'];
            if exist(textfilename,'file')==2
                displaytext = ['Finding artifact in ' ChannelsList{j}]; disp(displaytext);
                [M] = dlmread(textfilename); col1 = M(:,1); col2 = M(:,2);
                FlowBBtoExclude = false(length(PtBB_flow_time),1);
                PnasalBBtoExclude = false(length(PtBB_pnasal_time),1);
                for i=1:length(col1)
                    lefti=col1(i); if lefti<1, lefti=1; end
                    righti=col2(i); if righti>max(PtBB_flow_time), righti=max(PtBB_flow_time); end
                    FlowBBtoExclude = FlowBBtoExclude | ((PtBB_flow_time>lefti) & (PtBB_flow_time<righti));
                    righti=col2(i); if righti>max(PtBB_pnasal_time), righti=max(PtBB_pnasal_time); end
                    PnasalBBtoExclude = PnasalBBtoExclude | ((PtBB_pnasal_time>lefti) & (PtBB_pnasal_time<righti));
                end
                PtData_flow_pt(FlowBBtoExclude,:) = [];
                Amatrix2_flow_pt(FlowBBtoExclude,:) = [];
                PtBB_flow_time(FlowBBtoExclude) = [];
                Gtest_test_flow(FlowBBtoExclude) = [];
                weights_test_flow(FlowBBtoExclude) = [];
                PtData_pnasal_pt(PnasalBBtoExclude,:) = [];
                Amatrix2_pnasal_pt(PnasalBBtoExclude,:) = [];
                PtBB_pnasal_time(PnasalBBtoExclude) = [];
            end
        end
                 
        % match breaths from pnasal to real flow data
        PtFlowTime = PtBB_flow_time;
        PtPnasalTime = PtBB_pnasal_time;
        if ~isempty(PtPnasalTime) % empty for pts without pnasal
            str = [num2str(length(PtFlowTime)), ' breaths with Flow - after art removal']; disp(str);
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
        
        PtData_flow_pt = PtData_flow_pt(matchedBB_PtoF_clean,:);
        Amatrix2_flow_pt = Amatrix2_flow_pt(matchedBB_PtoF_clean,:);
        Gtest_test_flow = Gtest_test_flow(matchedBB_PtoF_clean);
        weights_test_flow = weights_test_flow(matchedBB_PtoF_clean);
        
        PtData_pnasal_pt(pnasal_exclusions,:) = [];
        Amatrix2_pnasal_pt(pnasal_exclusions,:) = [];
        
        % grow the data to keep
        Gtest_matched = [Gtest_matched; Gtest_test_flow];
        weights_matched = [weights_matched; weights_test_flow];
        
        PtData_matched = [PtData_matched; PtData_flow_pt];
        PtData_pnasal_matched =  [PtData_pnasal_matched; PtData_pnasal_pt];
        
        Amatrix2_flow_matched = [Amatrix2_flow_matched; Amatrix2_flow_pt];
        Amatrix2_pnasal_matched = [Amatrix2_pnasal_matched; Amatrix2_pnasal_pt];
        
        
        %                 figure(3);
        %                 plot(PtData_flow_pt.BB_time-StarttimeSpike, -0.44, 'r^');
        %                 plot(PtData_pnasal_pt.BB_time-StarttimeSpike, -0.54, 'b^');
        %                 grid on
        %                 xlabel('time'); ylabel('Flow (red) and Pnasal (blue)');
        %                 title('First row shows detected BB, second row show matched BB');
        
        %                 % do a VEVeup plot between measures to check alignment
        %                 VEVeup_flow = PtData_flow_pt.VE./PtData_flow_pt.Veup;
        %                 VEVeup_pnasal = PtData_pnasal_pt.VE./PtData_pnasal_pt.Veup;
        %                 figure(); scatter(VEVeup_pnasal, VEVeup_flow, 3);
    end % end of gathering data
    disp(' ');disp(' ');
    
    %% redo Weights
    % ToDo: check method of weights
    clear Ndata
    maxG=1.5;
    dx=0.2;
    xbins=[0 0.3:dx:0.9 maxG];
    %xbins=[0 0.1:0.05:1.1 1.5];
    
    for i=1:length(xbins)-1
        Ix=Gtest_matched>xbins(i)&Gtest_matched<=xbins(i+1);
        Ndata(i)=sum(Ix);
    end
    %Ndata = Ndata.^0.5;
    weightsbins = 1./(Ndata);
    %weightsbins = [2 1 1 1 0.5];
    weightsbins = weightsbins/mean(weightsbins);
    %weightsbins = weightsbins/weightsbins(end);
    
    weights_m = NaN*Gtest_matched;
    for i=1:length(xbins)-1
        Ix=Gtest_matched>=xbins(i)&Gtest_matched<=xbins(i+1);
        weights_m(Ix)=weightsbins(i);
    end
    weights_m = weights_m/nanmean(weights_m); % nnz(~isfinite(weights))
    
    useweights=1;
    if ~useweights %overwrite, use no weights
        weights_m = ones(length(weights_m),1);
    end
    
    %% redo
    NfeatureSteps=100;
    colofones = ones(length(Gtest_matched),1);
    predyL1O = NaN*Gtest_matched;
    predyL1O_array = NaN*ones(length(Gtest_matched),NfeatureSteps);
    predyL1O_array_flow = NaN*ones(length(Gtest_matched),NfeatureSteps);
    
    %% do backwards elimination
    try 
        % new for feature correlations
        RejectedFtrs=[];  % list of features that fail correlation test
        n=1; % counter
        RforAll=[]; % 3d array, each row is feature, each col is tx and each sheet is pt
        
        % then go through all patients again, this time training and testing in pnasal
        for subj=1:54 % subj=3
            if ~ismember(subj, Pnasal_list)
                disp(['Skipping pt ', num2str(subj), ': no pnasal']);
                continue
            end
            
            % Gtest_matched
            % weights_matched
            % PtData_matched
            % Amatrix2_pnasal_matched
            
            Labels_flow = Labels_Complete; % Start with full list
            Isubj=(PtData_matched.PT==subj);  % Isubj is the logical index of the L1O patient
            
            % set up training data
            Gtest_train = Gtest_matched;
            Gtest_train(Gtest_train>maxG)=maxG;
            Gtest_train(Isubj)=[];
            useOldWeights = 0;
            if useOldWeights % these are dragged from the old weights
                weights = weights_matched;
            else % these are the re-calculated weights
                weights = weights_m;
            end
            weights_train = weights;
            weights_train(Isubj)=[];
            weights_train=weights_train/nanmean(weights_train);
            
            colofones_train = colofones;
            colofones_train(Isubj)=[];
            
            Amatrix2_train = Amatrix2_pnasal_matched;
            Amatrix2_train(Isubj,:)=[];
            
            Amatrix2_train_flow = Amatrix2_flow_matched;
            Amatrix2_train_flow(Isubj,:)=[];

            %remove all NAN from all rows
            nanx = sum(isnan(Amatrix2_train),2)>0;   % nnz(nanx)
            nany = isnan(Gtest_train);               % nnz(nany)
            if nnz(nanx)>0 || nnz(nany)>0; disp('NaNs found'); end
            Amatrix2_train(nanx|nany,:)=[];
            Gtest_train(nanx|nany)=[];
            weights_train(nanx|nany)=[];
            weights_train=weights_train/nanmean(weights_train);
            colofones_train(nanx|nany)=[];
            
            % set up test data
            Gtest_test = Gtest_matched(Isubj);
            Gtest_test(Gtest_test>maxG)=maxG;
            weights_test = weights(Isubj);
            weights_test=weights_test/nanmean(weights_test);
            colofones_test = colofones(Isubj);
            Amatrix2_test = Amatrix2_pnasal_matched(Isubj,:);
            Amatrix2_test_flow = Amatrix2_flow_matched(Isubj,:);
            
            %remove all NAN from all rows
            nanx = sum(isnan(Amatrix2_test),2)>0;   % nnz(nanx)
            nany = isnan(Gtest_test);               % nnz(nany)
            if nnz(nanx)>0 || nnz(nany)>0; disp('NaNs found'); end
            Amatrix2_test(nanx|nany,:)=[];
            Gtest_test(nanx|nany)=[];
            weights_test(nanx|nany)=[];
            weights_test=weights_test/nanmean(weights_test);
            colofones_test(nanx|nany)=[];
            
            warning('off');
            maxp=Inf;
            
            Ftr_indx = 1:size(Amatrix2_train,2);
            Ftr_indx_flow = 1:size(Amatrix2_train_flow,2);
            
            %% remove features that do not correlate between flow and pnasal
            % calc correlations for flow and pnasal features
            numofftrs = length(Ftr_indx);
            r_vals = NaN(1,numofftrs); Rsq_vals = NaN(1,numofftrs);
            p_vals = NaN(1,numofftrs); MSE_vals = NaN(1,numofftrs);
            for ft=1:numofftrs
                % correlation coefficient
                r_val = weightedcorrs([Amatrix2_train(:,ft), Amatrix2_train_flow(:,ft)], weights_train);
                r_vals(ft) = r_val(1,2);
                % coefficient of variation
                [Rsq,p,MSE,~,~] = glmfitFast(Amatrix2_train(:,ft), Amatrix2_train_flow(:,ft), weights_train,0);
                %[Rsq_alt,p_alt,MSE_alt,~,~] = glmfitFast(Amatrix2_train_flow(:,ft), Amatrix2_train(:,ft), weights_train,0);
                Rsq_vals(ft) = Rsq;  p_vals(ft) = p(2);  MSE_vals(ft) = MSE;
            end
            MSE_vals(MSE_vals>1)=1; % cap maximum to one
            % what does the data look like?
            if 0 
                figure(); plot(Rsq_vals); hold on; plot(Rsq_vals);
                figure(); plot(p_vals);
                figure(); plot(MSE_vals);
            end
            % set a threshold for what is a 'bad' correlation
            if useRthreshold
                if 0 
                    r_Threshold = 3/4;  
                    Rsq_Threshold = 2/3; 
                else
                    r_Threshold = 0.866;  
                    Rsq_Threshold = 3/4; 
                end
                %r_Threshold = 0.8; % overwrite it
            else % keep all features, no exclusion.
                r_Threshold = -Inf;
                Rsq_Threshold = -Inf; 
            end
            % find bad correlations,
            Bad_r_vals_i = r_vals<r_Threshold; 
            Bad_Rsq_vals_i = Rsq_vals<Rsq_Threshold;
            NumBadFtrsR = nnz(Bad_r_vals_i);
            NumBadFtrsRsq = nnz(Bad_Rsq_vals_i);
            str = ['Patient ', num2str(subj)]; disp(str);
            str = ['Number of bad ftrs by r: ', num2str(NumBadFtrsR)]; disp(str);
            str = ['Number of bad ftrs by Rsq: ', num2str(NumBadFtrsRsq)]; disp(str);

            NumBadFtrsRandRsq = nnz(Bad_r_vals_i & Bad_Rsq_vals_i); % tells how many overlap exactly. 
            % if the number that overlap exactly is = to the lowest number, then do xor to find the extra ones
            NumBadFtrsRxorRsq = nnz(xor(Bad_r_vals_i,Bad_Rsq_vals_i));
            %ExtrasToExclude = find(xor(Bad_r_vals_i,Bad_Rsq_vals_i));
                   
            % report number of features rejected
            RejectedFtrs = [RejectedFtrs; ...
                [subj, NumBadFtrsR, NumBadFtrsRsq, NumBadFtrsRandRsq, NumBadFtrsRxorRsq]];
   
            % remove bad features
            %FtrsToExclude = find(Bad_Rsq_vals_i);
            FtrsToExclude = find(Bad_r_vals_i);
            Ftr_indx(FtrsToExclude) = [];
            Ftr_indx_flow(FtrsToExclude) = [];

            % what transform gives best Rsq flow-pnasal
            % this is the first half, getting all the data
            if 1
            % 495 features, made up of 165 features, in three transform variants
            % 150 features, made up of 50 features, in three transform variants
            totalftrs = length(Rsq_vals); % number of unique features
            numTxs = 3;  % number of transforms
            increment = totalftrs/numTxs; 
            numPts = 17; 
            
            for ftr=1:increment
                R=[];
                for tx_variant=1:3
                    % change 'R_vals' and 'Rsq_vals' as req'd
                    R = [R,r_vals(ftr+(increment*tx_variant-increment))]; 
                end
                RforAll(ftr,:,n) = R; 
            end
            n = n+1;
            end
            
        % temporary end, to get Rsq etc, 
        % if using this 'end', comment out  "end of L1O loop"  below @ ~line 760  
        %end    
        
        % this is the second half of the process for determining the 
        % "best" transform, averaged across all patients
        % this only makes sense to run if the 'end' above is active.
        if 0
        % find average of Rsq for "best" transform
        RsqforAll = RforAll; % make a copy to work with
        meanRsqforAll = mean(RforAll,3); % mean along third dimension
        maxTxPos = [];
        for ft=1:increment
           [~,loc]=max(meanRsqforAll(ft,:));
           maxTxPos = [maxTxPos, loc];
        end
        maxTxPos = maxTxPos';
        end

            %% do backwards elimination
            while length(Ftr_indx_flow)>=MaxNfeatures || maxp>alpha
                %% train in flow
                try
                    [~,Pvals_train_flow,~,b_train_flow,~]=glmfitFast(Amatrix2_train_flow(:,Ftr_indx_flow),Gtest_train,weights_train,1); %faster
                    Pvals_train_flow(1)=[];
                catch F2
                    disp(F2.getReport);
                end
                
                if length(Ftr_indx_flow)<=NfeatureSteps %start saving results
                    
                    %% save the training data - flow
                    predytrain_flow = [colofones_train Amatrix2_train_flow(:,Ftr_indx_flow)]*b_train_flow;
                    predytrain_flow(predytrain_flow>maxG)=maxG;
                    predytrain_flow(predytrain_flow<0)=0; %will be overwritten
                    
                    if 1 % calc and store training error
                        [RsqTrain_array_flow(subj,length(Ftr_indx_flow)), RTrain_array_flow(subj,length(Ftr_indx_flow))] = ...
                            UnivariateStats(predytrain_flow,Gtest_train, weights_train); % new method for R and Rsq 
                        %RsqTrain_array(subj,length(Ftr_indx)) = 1-nansum(weights_train.*(Gtest_train-predytrain).^2)/nansum(weights_train.*(Gtest_train-nanmean(Gtest_train)).^2);
                        ErrTrain_array(subj,length(Ftr_indx_flow)) = nanmean(weights_train.*abs(predytrain_flow-Gtest_train));
                        ErrRmsTrain_array(subj,length(Ftr_indx_flow)) = nanmean((weights_train.*(predytrain_flow-Gtest_train)).^2).^0.5; 
                    end

                    %% use beta from training, apply to test, save data - flow
                    predyL1O__flow = [colofones_test Amatrix2_test_flow(:,Ftr_indx_flow)]*b_train_flow;
                    predyL1O__flow(predyL1O__flow>maxG)=maxG;
                    predyL1O__flow(predyL1O__flow<0) = 0;
                    predyL1O_array_flow(Isubj,length(Ftr_indx_flow)) = predyL1O__flow; 
                     
                    %% store the labels and beta array - flow
                    labels_Step_Subj_flow{subj,length(Ftr_indx_flow)} = Ftr_indx_flow;
                    beta_array_flow{subj,length(Ftr_indx_flow)} = b_train_flow;
                   
                    %% use beta from training in flow, apply to test data in pnasal, save data - pnasal
                    predyL1O_ = [colofones_test Amatrix2_test(:,Ftr_indx_flow)]*b_train_flow;
                    predyL1O_(predyL1O_>maxG)=maxG;
                    predyL1O_(predyL1O_<0) = 0;
                    predyL1O_array(Isubj,length(Ftr_indx_flow)) = predyL1O_;

                end
                
                [maxp_flow,remi_flow]=max(Pvals_train_flow);

                if length(Ftr_indx_flow)==1 || (~neverbreak && maxp<alpha && (length(Ftr_indx_flow)-1)<=MaxNfeatures)
                    break
                end
                
                disp(['F  - Removing Ftr: ', num2str(Ftr_indx_flow(remi_flow)), ', p= ' num2str(maxp_flow)]);
                
                Ftr_indx_flow(remi_flow)=[];
                Labels_flow(remi_flow)=[];
            end

            toc
       
       % comment this 'end' out if just doing feature correlation testing
       end  % end of L1O loop,
        
    catch me
        disp(me.getReport);
    end
    
    % display processing time
    delta_t = etime(clock, t_start_L1O); % delta in seconds
    D = duration(0,0,round(delta_t),'Format','hh:mm:ss'); % duration in HH:MM:SS
    disp(' '); % add row space for visual clarity in command window
    
    displaytext = ['L1O linear regresion process complete. Total time: ', char(D), ' (hh:mm:ss)'];
    HowSoonIsNow = datestr(datetime('now'),'yyyymmdd');
    
    savestr = ['_LinRegWorkspace_TestNTrainPnasal_', HowSoonIsNow, experimentnumber, '.mat'];
    
    disp(displaytext);
    
    str=['Saving to: ', filename(1:end-4), savestr]; disp(str);
    save([filename(1:end-4), savestr]);
    
    diary off
    
    keyboard
    
end % end of L1O run

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>> END OF RUN <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------


%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
close all
clear
clc
addpath('C:\Users\uqdmann\Dropbox\PUPbeta_git\');
cd('C:\Users\uqdmann\Dropbox\PUPbeta_git\FlowAndEdiPesQAO');
%datadir = 'C:\PSG_Data\FlowDrive\FeatureSpaces\';
datadir = 'C:\PSG_Data\FlowDrive\Analyzed\';
%load('C:\Users\uqdmann\Dropbox\PUPbeta_git\FeatureSpaces\FlowDrive_only25Hz_FeatureSpace_AutoRef2_Edi_Clean_withPnasalMatched_LinRegWorkspace_TestNTrainPnasal_20180129.mat');
%load('C:\PSG_Data\FlowDrive\FeatureSpaces\FlowDrive_only25Hz_Normalized_FeatureSpace_AutoRef2_Edi_Clean_withPnasalMatched_LinRegWorkspace_TestNTrainPnasal_20180227.mat');
%load('C:\PSG_Data\FlowDrive\FeatureSpaces\FlowDrive_only25Hz_Normalized_FeatureSpace_AutoRef2_Edi_Clean_withPnasalMatched_LinRegWorkspace_TestNTrainPnasal_20180302.mat');
%load([datadir, 'FlowDrive_only25Hz_ReNormalized_FeatureSpace_AutoRef2_Edi_Clean_LinRegWorkspace_TestNTrainPnasal_20180430_ReNormalizedFtrSet.mat']);
load([datadir, 'FlowDrive_only25Hz_ReNormalized_FeatureSpace_AutoRef2_Edi_Clean_LinRegWorkspace_TestNTrainPnasal_20180530.mat']);

if 1
    close; clear; clc
    addpath('C:\Users\uqdmann\Dropbox\PUPbeta_git\');
    cd('C:\Users\uqdmann\Dropbox\PUPbeta_git\FlowAndEdiPesQAO');
    datadir = 'C:\PSG_Data\FlowDrive\Analyzed\ExclBadR\';
   
    experimentnumber = '_exp11';
    load([datadir, 'FlowDrive_only25Hz_ReNormalized_FeatureSpace_AutoRef2_Edi_Clean_LinRegWorkspace_TestNTrainPnasal', experimentnumber, '.mat']);
    %exp = '_exp10'; % rewrite for saving, process for exp 6-10 uses data from exp1-5
    
    % should just load variables of interest
    % 'PtData','Amatrix2','FeatureNames','etc...'
    
    pnasal = 0 ; % Flow/Pnasal switch
end

%% Use All breaths for classification performance plots
SleepOnly = 0;  % 0 for All breaths, 1 for sleep only breaths
IncludePnasal = 1;
[BB, BB_] = getAllorSleepOnlyBB(SleepOnly, IncludePnasal, PtData, PtData_matched);
str = ['Flow - Using ',num2str(nnz(BB)),' of ', num2str(nnz(PtData.PT))]; disp(str);
str = ['Pnasal - Using ',num2str(nnz(BB_)),' of ', num2str(nnz(PtData_matched.PT))]; disp(str);

% Select data to use
PlotUptoNFtrs = 40; % normally 100
ftrnum=25; % 25
pts_flow = unique(PtData.PT); % the patients with flow data
pts_pnasal = unique(PtData_matched.PT); % the patients with pnasal data

predy_flow = predyL1O_array_flow(BB_,ftrnum);
predy_pnasal = predyL1O_array(BB_,ftrnum);
Yval = Gtest_matched(BB_);

savefigas = '';% 'saveasTIFF'; % options are saveasPNG, saveasFIG, and saveasTIFF
closefigs = 0;

TickFntSz = 12;
LabelFntSz = 18;
FntSz = 18;

%% Process test data. (flow and pnasal)
% for all pts combined
clear RsqL1O_ RsqL1O_dlm R_dlm R_dlm2 ErrL1O_ ErrL1Orms_
for i=1:PlotUptoNFtrs % normally size(predyL1O_array,2)
    %RsqL1O_(i) = 1-nansum(weights(BB_).*(Gtest_matched(BB_)-predyL1O_array(BB_,i)).^2)/nansum(weights(BB_).*(Gtest_matched(BB_)-nanmean(Gtest_matched(BB_))).^2);
    [RsqL1O_pnasal(i), RL1O_pnasal(i)] = ...
       UnivariateStats(predyL1O_array(BB_,i),Gtest_matched(BB_), weights(BB_)); % new method for R and Rsq     
    [RsqL1O_flow(i), RL1O_flow(i)] = ...
       UnivariateStats(predyL1O_array_flow(BB_,i),Gtest_matched(BB_), weights(BB_)); % new method for R and Rsq 
    ErrL1O_pnasal(i) = nanmean(weights(BB_).*abs(predyL1O_array((BB_),i)-Gtest_matched(BB_)));
    ErrL1Orms_pnasal(i) = nanmean((weights(BB_).*(predyL1O_array((BB_),i)-Gtest_matched(BB_))).^2).^0.5;
    ErrL1O_flow(i) = nanmean(weights(BB_).*abs(predyL1O_array_flow((BB_),i)-Gtest_matched(BB_)));
    ErrL1Orms_flow(i) = nanmean((weights(BB_).*(predyL1O_array_flow((BB_),i)-Gtest_matched(BB_))).^2).^0.5;
end



%% For the average of all pts combined, compare old R and Rsq with new R and Rsq
figure(1); clf(figure(1));
subplot(1,2,1);
stairs(RL1O_pnasal.^0.5,'g'); % is wrong, only testing
hold on; 
stairs(RL1O_pnasal,'b'); 
stairs(RsqL1O_pnasal,'r'); 
stairs(ErrL1O_pnasal, 'k');
xlim([0, PlotUptoNFtrs]); ylim([0,1]);
title('Pnasal'); legend('wRong','R','R squared','MAE');

subplot(1,2,2);
stairs(RsqL1O_flow.^0.5,'g'); % is wrong, only testing
hold on;
stairs(RL1O_flow,'b'); 
stairs(RsqL1O_flow,'r');
stairs(ErrL1O_flow, 'k');
xlim([0, PlotUptoNFtrs]); ylim([0,1]);
title('Flow'); legend('wRong','R','R squared','MAE');
suptitle('average all patients during testing');

%% per pt training data
figure(2); clf(figure(2));
subplot(1,2,1);
% for pt = 1:54
%     if ismember(pt, pts_pnasal)
%         plot(RsqTrain_array(pt,:).^0.5, 'g'); % is wrong, only testing
%         hold on;
%         plot(RTrain_array(pt,:), 'b');
%         plot(RsqTrain_array(pt,:), 'r');
%         %plot(ErrRmsTrain_array(pt,:), 'g');
%         plot(ErrTrain_array(pt,:), 'k');
%         
%         hold on;
%     end
% end
% xlim([0, PlotUptoNFtrs]); ylim([0,1]);
% title('Pnasal'); legend('wRong','R','R squared','MAE');

subplot(1,2,2);
for pt = 1:54
    if ismember(pt, pts_pnasal)
        plot(RsqTrain_array_flow(pt,:).^0.5, 'g'); % is wrong, only testing
        hold on;
        plot(RTrain_array_flow(pt,:), 'b'); 
        plot(RsqTrain_array_flow(pt,:), 'r');
         plot(ErrTrain_array(pt,:), 'k');
    end
end
xlim([0, PlotUptoNFtrs]); ylim([0,1]);
title('Flow'); legend('wRong','R','R squared');
suptitle('per patient during training');


%% Figure E4
% set max number of features for training error
% temp set at 49, which is the min num of ftrs that any pt started with
ErrorMaxN = 49; % PlotUptoNFtrs; % normally 100

figure(84); clf(figure(84)); fig = gcf;
fig.Color = [1 1 1]; fig.Units = 'inches';
fig.Position = [-10  -1   9   8];

% performance vs number of features
subplot(2,2,1);
stairs([RL1O_flow(1:PlotUptoNFtrs);RsqL1O_flow(1:PlotUptoNFtrs);ErrL1O_flow(1:PlotUptoNFtrs)]'); hold on;
legend('R', 'Rsq');
currentYlim=ylim();
plot([ftrnum, ftrnum], [0, currentYlim(2)],'k-', 'linewidth', 1.5);
ax = gca; set(ax,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz);
xlabel('Number of Features', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
ylabel('Performance (Testing)', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow'); 
ylim([0 1]); xlim([0 50]); axis square; title('Flow');

subplot(2,2,3);
stairs([RL1O_pnasal(1:PlotUptoNFtrs);RsqL1O_pnasal(1:PlotUptoNFtrs);ErrL1O_pnasal(1:PlotUptoNFtrs)]'); hold on;
legend('R', 'Rsq');
currentYlim=ylim();
plot([ftrnum, ftrnum], [0, currentYlim(2)],'k-', 'linewidth', 1.5);
ax = gca; set(ax,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz);
xlabel('Number of Features', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
ylabel('Performance (Testing)', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow'); 
ylim([0 1]); xlim([0 50]); axis square; title('Pnasal');

% confusion mat
subplot(2,2,2);
customcmap = GetCustomColorMap('gray'); % 'SS'
g = NaN*Yval;
ghat = NaN*predy_flow;
classes = 5;
switch classes
    case 5
        classcutoffs = [0.9 0.7 0.5 0.3];
        labeltext={'Normal','Mild','Moderate','Severe','V.Severe'};
    case 4
        classcutoffs = [0.9 0.5 0.3];
        labeltext={'Normal','Mild','Moderate','Severe'};
    case 3
        classcutoffs = [0.9 0.5];
        labeltext={'Normal','Intermediate','Flow Limited'};
    otherwise
        classcutoffs = [0.7]; %normal mild moderate severe v.severe / normal borderline mild moderate severe
        labeltext={'Normal','Flow Limited'};
end
Nclasses=length(classcutoffs)+1;
g(Yval>classcutoffs(1))=1;
ghat(predy_flow>classcutoffs(1))=1;
for i=2:Nclasses
    g(Yval<classcutoffs(i-1))=i;
    ghat(predy_flow<classcutoffs(i-1))=i;
end
[C,order] = confusionmat(g,ghat);
%sumactual=sum(C')';
sumactual=sum(C,2);
sumestimated=sum(C);
%C_Factual=C./sum(C')'
C_Factual=C./sum(C,2)*100; %rows are actual, cols are estimated
C_Festimated=C./sum(C)*100;
C_Total = (C./sum(C(:)))*100;
 
if 1
    AccA_C_Factual = C_Factual.*diag(ones(1,length(C)));
    AccB_C_Factual = C_Factual.*(diag(ones(1,length(C)-1),1)+diag(ones(1,length(C)-1),-1));
    AccAactual = mean(sum(AccA_C_Factual,2));
    AccBactual = mean(sum(AccA_C_Factual,2)) + mean(sum(AccB_C_Factual,2));
    AccA_C_Festimated = C_Festimated.*diag(ones(1,length(C)));
    AccB_C_Festimated = C_Festimated.*(diag(ones(1,length(C)-1),1)+diag(ones(1,length(C)-1),-1));
    AccAestimated = mean(sum(AccA_C_Festimated));
    AccBestimated = mean(sum(AccA_C_Festimated)) + mean(sum(AccB_C_Festimated));
    ACCs = [AccAactual AccBactual; AccAestimated AccBestimated]
    % first row is actual, second row is estimate - use second row
    %second col gives accuracy if accepting next-category error as correct
end

x = order'; y = order';
if 1; C1 = C_Festimated; else; C1 = C_Factual; end
%C1 = C_Total;
xSplit = diff(x)/2;                 % Find edge points
ySplit = diff(y)/2;
xEdges = [x(1)-xSplit(1) x(2:end)-xSplit x(end)+xSplit(end)];
yEdges = [y(1)-ySplit(1) y(2:end)-ySplit y(end)+ySplit(end)];
[XGrid, YGrid] = meshgrid(xEdges,yEdges);
YGrid = flipud(YGrid);              % To match expected behavior
XGrid = fliplr(XGrid);
C2 = [[C1 zeros(size(C1,1),1)] ; zeros(1,size(C1,2)+1)];% Last row/col ignored
pcolor(XGrid,YGrid,(1-(1-C2/100).^2)*100)
hold on                             % Plot original data points
[X,Y] = meshgrid(x,y);
colormap(gcf,'hot')
%colorbar();
set(gcf,'colormap',customcmap);
C1 = flipud(C1);
C1 = fliplr(C1);
C1 = C1';
for i=1:size(C1,1)
    for j=1:size(C1,2)
        if C1(i,j)<((max(max(C1))-min(min(C1)))/2+min(min(C1)))
            textcolor=[1 1 1];
        else
            textcolor=[0 0 0];
        end
        text(x(i),y(j),num2str(round(C1(i,j))),'color',textcolor,'horizontalalignment','center','fontname','arial narrow')
    end
end
ax = gca; set(ax,'fontname','arial narrow','FontSize', TickFntSz);
yticks(y); yticklabels(gca,fliplr(labeltext));
xticks(x); xticklabels(gca,fliplr(labeltext));
xlabel('Flow Shape Classification', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
ylabel('Gold Standard Classification', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
axis square; title('Flow');

% confusion mat
subplot(2,2,4);
customcmap = GetCustomColorMap('gray'); % 'SS'
g = NaN*Yval;
ghat = NaN*predy_pnasal;
classes = 5;
switch classes
    case 5
        classcutoffs = [0.9 0.7 0.5 0.3];
        labeltext={'Normal','Mild','Moderate','Severe','V.Severe'};
    case 4
        classcutoffs = [0.9 0.5 0.3];
        labeltext={'Normal','Mild','Moderate','Severe'};
    case 3
        classcutoffs = [0.9 0.5];
        labeltext={'Normal','Intermediate','Flow Limited'};
    otherwise
        classcutoffs = [0.7]; %normal mild moderate severe v.severe / normal borderline mild moderate severe
        labeltext={'Normal','Flow Limited'};
end
Nclasses=length(classcutoffs)+1;
g(Yval>classcutoffs(1))=1;
ghat(predy_pnasal>classcutoffs(1))=1;
for i=2:Nclasses
    g(Yval<classcutoffs(i-1))=i;
    ghat(predy_pnasal<classcutoffs(i-1))=i;
end
[C,order] = confusionmat(g,ghat);
%sumactual=sum(C')';
sumactual=sum(C,2);
sumestimated=sum(C);
%C_Factual=C./sum(C')'
C_Factual=C./sum(C,2)*100; %rows are actual, cols are estimated
C_Festimated=C./sum(C)*100;
C_Total = (C./sum(C(:)))*100;
 
if 1
    AccA_C_Factual = C_Factual.*diag(ones(1,length(C)));
    AccB_C_Factual = C_Factual.*(diag(ones(1,length(C)-1),1)+diag(ones(1,length(C)-1),-1));
    AccAactual = mean(sum(AccA_C_Factual,2));
    AccBactual = mean(sum(AccA_C_Factual,2)) + mean(sum(AccB_C_Factual,2));
    AccA_C_Festimated = C_Festimated.*diag(ones(1,length(C)));
    AccB_C_Festimated = C_Festimated.*(diag(ones(1,length(C)-1),1)+diag(ones(1,length(C)-1),-1));
    AccAestimated = mean(sum(AccA_C_Festimated));
    AccBestimated = mean(sum(AccA_C_Festimated)) + mean(sum(AccB_C_Festimated));
    ACCs = [AccAactual AccBactual; AccAestimated AccBestimated]
    % first row is actual, second row is estimate - use second row
    %second col gives accuracy if accepting next-category error as correct
end

x = order'; y = order';
if 1; C1 = C_Festimated; else; C1 = C_Factual; end
%C1 = C_Total;
xSplit = diff(x)/2;                 % Find edge points
ySplit = diff(y)/2;
xEdges = [x(1)-xSplit(1) x(2:end)-xSplit x(end)+xSplit(end)];
yEdges = [y(1)-ySplit(1) y(2:end)-ySplit y(end)+ySplit(end)];
[XGrid, YGrid] = meshgrid(xEdges,yEdges);
YGrid = flipud(YGrid);              % To match expected behavior
XGrid = fliplr(XGrid);
C2 = [[C1 zeros(size(C1,1),1)] ; zeros(1,size(C1,2)+1)];% Last row/col ignored
pcolor(XGrid,YGrid,(1-(1-C2/100).^2)*100)
hold on                             % Plot original data points
[X,Y] = meshgrid(x,y);
colormap(gcf,'hot')
%colorbar();
set(gcf,'colormap',customcmap);
C1 = flipud(C1);
C1 = fliplr(C1);
C1 = C1';
for i=1:size(C1,1)
    for j=1:size(C1,2)
        if C1(i,j)<((max(max(C1))-min(min(C1)))/2+min(min(C1)))
            textcolor=[1 1 1];
        else
            textcolor=[0 0 0];
        end
        text(x(i),y(j),num2str(round(C1(i,j))),'color',textcolor,'horizontalalignment','center','fontname','arial narrow')
    end
end
ax = gca; set(ax,'fontname','arial narrow','FontSize', TickFntSz);
yticks(y); yticklabels(gca,fliplr(labeltext));
xticks(x); xticklabels(gca,fliplr(labeltext));
xlabel('Pnasal Shape Classification', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
ylabel('Gold Standard Classification', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
axis square; title('Pnasal');

str = ['..\Figures\CorrelationTesting\Figure_E4_F',experimentnumber]; 

savefigas = 'saveasPNG';
switch savefigas     
    case 'saveasTIFF'; print(fig, str, '-dtiff', '-r1000');     
    case 'saveasPNG'; saveas(fig, str, 'png');      
    case 'saveasFIG'; savefig(str);  
end 
if closefigs; close(fig); end
savefigas = '';

%% back-up tables
Old_ApBB = RemovedBB_Apnoea;
Old_LFBB = RemovedBB_LowFlow;

%% Add apnea and low flow breaths back in for everything hereafter
load([datadir, 'PnasalDrive_only25Hz_ReNormalized_FeatureSpace_AutoRef2_Edi_Clean.mat'],'RemovedBB_Apnoea');
load([datadir, 'PnasalDrive_only25Hz_ReNormalized_FeatureSpace_AutoRef2_Edi_Clean.mat'],'RemovedBB_LowFlow');
pts_pnasal = unique(PtData_matched.PT);
numBBinTest = NaN(54,6);
PredY_pnasal = []; PredY_flow = []; Gtest = []; PTarray=[];
for subj=1:54
    if ~ismember(subj, pts_pnasal)
        continue
    end
    Isubj=(PtData_matched.PT==subj); % find all the breaths that belong to this pt, this is just for counting
    numBBinTest(subj,1) = nnz(Isubj); % the number of breaths for this pt
    
    Isubj=(PtData_matched.PT==subj)&(PtData_matched.Hypnog<4)&(PtData_matched.Ar==0); % find the sleep only breaths that belong to this pt
    numBBinTest(subj,2) = nnz(Isubj); % the number of sleep breaths for this pt
    
    Gtest_pt = Gtest_matched(Isubj);
    PredY_pt_pnasal = predyL1O_array(Isubj,ftrnum); % this is pnasal
    PredY_pt_flow = predyL1O_array_flow(Isubj,ftrnum);
    
    % add apneoa breaths as 0.1 to both Gtest_pt and PredY_pt
    indAP = find(RemovedBB_Apnoea.Pt == subj);
    indLF = find(RemovedBB_LowFlow.Pt == subj);
    numApBB = RemovedBB_Apnoea{indAP,2};
    if isempty(numApBB); numApBB=0; end
    numBBinTest(subj,3) = numApBB;
    numLFBB = RemovedBB_LowFlow{indLF,2};
    if isempty(numLFBB); numLFBB=0; end
    numBBinTest(subj,4) = numLFBB;
    numRemBB = numApBB + numLFBB;
    
    PredY_pt_pnasal = [PredY_pt_pnasal;ones(numRemBB,1)*0.09];
    PredY_pt_flow = [PredY_pt_flow;ones(numRemBB,1)*0.09];
    numBBinTest(subj,5) = length(PredY_pt_pnasal);
    Gtest_pt = [Gtest_pt;ones(numRemBB,1)*0.09];
    numBBinTest(subj,6) = length(Gtest_pt);
    
    PTarray = [PTarray; ones(length(PredY_pt_pnasal),1)*subj]; % same as PtData.PT but includes apnea breaths
    PredY_pnasal = [PredY_pnasal; PredY_pt_pnasal];
    PredY_flow = [PredY_flow; PredY_pt_flow];
    Gtest = [Gtest; Gtest_pt]; 
end
%
% tidy up table, remove nan rows
% rows in numBBinTest are each patient
% cols in numBBinTest are:
% Total BB, Sleep BB, Ap BB, Low flow BB, Total BB to use (with AP)
numBBinTest = numBBinTest(pts_pnasal,:);
if length(PredY_pnasal) ~= length(Gtest); keyboard; end                % check same length
if nnz(isnan(PredY_pnasal))>0 || nnz(isnan(Gtest))>0; keyboard; end    % check for any NaN's
if nnz(PredY_pnasal<0)>0; Predy(PredY_pnasal<0)=0; end                      % force lower limit
if nnz(PredY_pnasal>maxG)>0; Predy(PredY_pnasal>maxG)=1.5; end                % force upper limit

%% get the AHI data
[AHI_perPT, AHI_perPT_table] = getAHI_postanalysis();
AHI_perPT_ = AHI_perPT(pts_pnasal,1);

%% histogram of pred vs actual VEVdrive, and do per pt processing in loop
% get per pt median data
% get per pt <threshold data

if 1; edges = [0:0.1:1.1]; else edges = xbins; end

figure(22); clf(figure(22)); fig = gcf;
fig.Color = [1 1 1]; fig.Units = 'inches';
fig.Position = [0.5   0.5   19    5];

subplot(1,3,1);
histogram(Gtest,edges, 'facealpha',0.2, 'edgealpha', 0.2, 'normalization','probability'); hold on;
currentYlim = ylim();
plot([median(Gtest), median(Gtest)], [0, currentYlim(2)],'k-', 'linewidth', 2);
titlestr = ['Actual (median = ', num2str(round(median(Gtest),2)), ')']; 
title(titlestr); xlabel('{\itflow:drive} (%)'); ylabel('Relative probability'); %yticks([]);
xlim([-0.05 1.15]); ax = gca; ax.XTick=[0:0.2:1]; xticklabels(ax, {'0', '20', '40', '60', '80', '100'});
box off

subplot(1,3,2);
histogram(PredY_pnasal,edges, 'facealpha',0.2, 'edgealpha', 0.2, 'normalization','probability'); hold on;
currentYlim = ylim();
plot([median(PredY_pnasal), median(PredY_pnasal)], [0, currentYlim(2)],'k-', 'linewidth', 2);
titlestr = ['Pnasal Predicted (median = ', num2str(round(median(PredY_pnasal),2)), ')']; 
title(titlestr);xlabel('{\itflow:drive} (%)'); ylabel('Relative probability'); 
xlim([-0.05 1.15]); ax = gca; ax.XTick=[0:0.2:1]; xticklabels(ax, {'0', '20', '40', '60', '80', '100'});
box off

subplot(1,3,3);
histogram(PredY_flow,edges, 'facealpha',0.2, 'edgealpha', 0.2, 'normalization','probability'); hold on;
currentYlim = ylim();
plot([median(PredY_flow), median(PredY_flow)], [0, currentYlim(2)],'k-', 'linewidth', 2);
titlestr = ['Flow Predicted (median = ', num2str(round(median(PredY_flow),2)), ')']; 
title(titlestr);xlabel('{\itflow:drive} (%)'); ylabel('Relative probability'); 
xlim([-0.05 1.15]); ax = gca; ax.XTick=[0:0.2:1]; xticklabels(ax, {'0', '20', '40', '60', '80', '100'});
box off

supstr=['Histograms of Actual and Pnasal Predicted {\itflow:drive} (All pts, incl Ap)'];
suptitle(supstr);
str = ['..\Figures\', supstr];
% switch savefigas     
%     case 'saveasTIFF'; print(fig, str, '-dtiff', '-r1000');     
%     case 'saveasPNG'; saveas(fig, str, 'png');      
%     case 'saveasFIG'; savefig(str);  
% end
if closefigs; close(fig); end

%% histograms for individuals
Gtest_avg = NaN(54,1); PredY_avg = NaN(54,1);
thres1 = 0.5; thres2 = 0.7;
Gtest_thres1 = NaN(54,1); Gtest_thres2 = NaN(54,1);
PredY_thres1 = NaN(54,1); PredY_thres2 = NaN(54,1);

doHistogramsPlots = 0;

for subj=1:54 % set pt num or all, no data for PT=1
    if ~ismember(subj, pts_pnasal) %PT_list)
        continue
    end
    
    Isubj = PTarray==subj;      % get the pt BB's
    Gtest_pt = Gtest(Isubj);    % actual
    PredY_pt_pnasal = PredY_pnasal(Isubj);    % pred
    
    % medians
    Gtest_avg(subj) = median(Gtest_pt);
    PredY_avg(subj) = median(PredY_pt_pnasal);
    
    % below threshold
    Gtest_thres1(subj) = nnz(Gtest_pt<thres1);
    Gtest_thres2(subj) = nnz(Gtest_pt<thres2);
    PredY_thres1(subj) = nnz(PredY_pt_pnasal<thres1);
    PredY_thres2(subj) = nnz(PredY_pt_pnasal<thres2);
    
    if doHistogramsPlots
        figure(300+subj); clf(figure(300+subj));fig = gcf;
        fig.Color = [1 1 1]; % set background colour to white
        fig.Units = 'inches';
        fig.Position = [0.5   0.5   12    5];
        
        subplot(1,2,1);
        h1 = histogram(Gtest_pt,edges, 'facealpha',0.2, 'edgealpha', 0.2, 'normalization','pdf');%'probability');
        hold on;
        currentYlim = ylim();
        %plot([median(Gtest_pt), median(Gtest_pt)], [0, max(h1.Values)],'k-', 'linewidth', 2);
        plot([median(Gtest_pt), median(Gtest_pt)], [0, currentYlim(2)],'k-', 'linewidth', 2);
        titlestr = ['Actual (median = ', num2str(round(median(Gtest_pt),2)), ')'];
        title(titlestr); xlabel('{\itflow:drive} (%)'); ylabel('Relative probability'); %yticks([]);
        xlim([-0.05 1.15]); ax = gca; ax.XTick=[0:0.2:1]; xticklabels(ax, {'0', '20', '40', '60', '80', '100'});
        box off

        subplot(1,2,2);
        h2 = histogram(PredY_pt_pnasal,edges, 'facealpha',0.2, 'edgealpha', 0.2, 'normalization','pdf');%'probability');
        hold on;
        currentYlim = ylim();
        %plot([median(PredY_pt), median(PredY_pt)], [0, max(h2.Values)],'k-', 'linewidth', 2);
        plot([median(PredY_pt_pnasal), median(PredY_pt_pnasal)], [0, currentYlim(2)],'k-', 'linewidth', 2);
        titlestr = ['Predicted (median = ', num2str(round(median(PredY_pt_pnasal),2)), ')'];
        title(titlestr); xlabel('{\itflow:drive} (%)'); ylabel('Relative probability');
        xlim([-0.05 1.15]); ax = gca; ax.XTick=[0:0.2:1]; xticklabels(ax, {'0', '20', '40', '60', '80', '100'});
        box off
        
        supstr = ['Histograms of Actual and Pnasal Predicted {\itflow:drive}, pt ',num2str(subj),', AHI ', num2str(round(AHI_perPT(subj,1))), ', withApBB'];
        suptitle(supstr);
        
        % save
        str = ['..\Figures\', supstr];
        switch savefigas
            case 'saveasTIFF'; print(fig, str, '-dtiff', '-r1000');
            case 'saveasPNG'; saveas(fig, str, 'png');
            case 'saveasFIG'; savefig(str);
        end
        if closefigs; close(fig); end
    end
end

% tidy up
exclude = isnan(Gtest_avg); %a=find(exclude==0)
Gtest_avg(exclude) = []; PredY_avg(exclude) = [];
Gtest_thres1(exclude) = []; Gtest_thres2(exclude) = [];
PredY_thres1(exclude) = []; PredY_thres2(exclude) = [];
Gtest_thres1_pFL = Gtest_thres1./numBBinTest(:,5); % Proportion FL
Gtest_thres2_pFL = Gtest_thres2./numBBinTest(:,5);
PredY_thres1_pFL = PredY_thres1./numBBinTest(:,5);
PredY_thres2_pFL = PredY_thres2./numBBinTest(:,5);


%% Figure_6 in manuscript
figure(21); clf(figure(21)); fig = gcf;
fig.Color = [1 1 1]; fig.Units = 'inches';
fig.Position = [1  1   16   5];

% flow predy vs pnasal predy
subplot(1,3,1)
facealpha = 0.05; % was 0.04
facecolor = [0 0 0]; %was [0 0 0]
scatter(predy_pnasal.*100,predy_flow.*100,2,'filled','markerfacealpha',facealpha, 'markerfacecolor', facecolor); hold on;
[r_1_6a, p_1] =corr(predy_pnasal,predy_flow); 
%[r_2, r_dm] = UnivariateStats(predy,predy_flow); % don't use without weights
w = ones(length(predy_flow),1);
[Rsq_6a,p,MSE,beta,ypred]=glmfitFast(predy_pnasal,predy_flow,w,1);
text(110, 15, ['r = ', num2str(round(r_1_6a,2))]);
str=['R^2 = ', num2str(round(Rsq_6a,2))]; text(110, 5, str);
xlim([0 150]); ylim([0 150]); 
    ax = gca; set(ax,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz);
    ax.XTick=[0:25:150]; xticklabels(ax, {'0', '25', '50', '75', '100', '125'});
    ax.YTick=[0:25:150]; yticklabels(ax, {'0', '25', '50', '75', '100', '125'});
    xlabel('Pnasal Shape Predicted {\itflow:drive} (%)', 'FontSize', FntSz, 'FontName', 'Arial Narrow');
    ylabel('Flow Shape Predicted {\itflow:drive} (%)', 'FontSize', FntSz, 'FontName', 'Arial Narrow'); 
text(-45, 148, 'A', 'FontSize', 20, 'FontWeight', 'Bold');
axis square

% scatter with box overlay
subplot(1,3,2);
scatter(100*predy_pnasal,100*Yval,2,'filled','markerfacealpha',0.06, 'markerfacecolor', [0 0 0]); hold on;
dx=0.2; xbins=[0 0.3:dx:0.9 1.5];
clear medianX medianY upperIQRY lowerIQRY upperIQRY2 lowerIQRY2
for i=1:length(xbins)-1
    Ix=predy_pnasal>xbins(i)&predy_pnasal<xbins(i+1);
    medianX(i)=prctile(predy_pnasal(Ix),50);
    medianY(i)=prctile(Yval(Ix),50);
    upperIQRY(i)=prctile(Yval(Ix),75);
    lowerIQRY(i)=prctile(Yval(Ix),25);
    upperIQRY2(i)=prctile(Yval(Ix),90);
    lowerIQRY2(i)=prctile(Yval(Ix),10);
end
plotbarwitherrorsmedian(100*medianY,100*medianX,100*upperIQRY2,100*lowerIQRY2,zeros(length(medianY),3),NaN,100*0.005,100*0.01);
plotbarwitherrorsmedian(100*medianY,100*medianX,100*upperIQRY,100*lowerIQRY,zeros(length(medianY),3),NaN,100*dx/2,100*0.01);
xlim([0 150]); 
ax = gca;set(ax,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz);
ax.XTick=[0:25:150]; xticklabels(ax, {'0', '25', '50', '75', '100', '125'});
ax.YTick=[0:25:150]; yticklabels(ax, {'0', '25', '50', '75', '100', '125'});
xlabel('Pnasal Shape Predicted {\itflow:drive} (%)', 'FontSize', FntSz, 'FontName', 'Arial Narrow');
ylabel('Gold Standard {\itflow:drive} (%)', 'FontSize', FntSz, 'FontName', 'Arial Narrow');
%[r_1, p_1] =corr(predy,Yval);  % don't use corr if have weights
[r_2, r_dm_6b] = UnivariateStats(predy_pnasal,Yval,weights);
[Rsq_6b,p,MSE,beta,ypred]=glmfitFast(predy_pnasal,Yval,weights,1);
str=['r = ', num2str(round(r_dm_6b,2))]; text(110, 15, str);
str=['R^2 = ', num2str(round(Rsq_6b,2))]; text(110, 5, str);
text(-45, 148, 'B', 'FontSize', 20, 'FontWeight', 'Bold'); % Add label to plot space
axis square

subplot(1,3,3);
[r_1_6c, p_1] = plotregressionwithSEM(PredY_avg.*100, Gtest_avg.*100);
% [r_c, p_c] =corr(PredY_avg,Gtest_avg); % don't need corr if have plotreg
% [r_2, r_dm] = UnivariateStats(PredY_avg,Gtest_avg); % don't use without weights
w = ones(length(Gtest_avg),1);
[Rsq_6c,p,MSE,beta,ypred]=glmfitFast(PredY_avg, Gtest_avg,w,1);
xlim([-5 110]); ylim([-5 110]);
ax = gca; set(ax,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz);
ax.XTick=[0:25:125];
ax.YTick=[0:25:125];
xlabel({'Pnasal Shape Predicted {\itflow:drive}','median (%)'}, 'FontSize', FntSz, 'FontName', 'Arial Narrow');
ylabel({'Gold Standard {\itflow:drive}','median (%)'}, 'FontSize', FntSz, 'FontName', 'Arial Narrow'); 
axis square
str=['r = ', num2str(round(r_1_6c,2))]; text(80, 7, str);
str=['R^2 = ', num2str(round(Rsq_6c,2))]; text(80, 0, str);
text(-45, 108, 'C', 'FontSize', 20, 'FontWeight', 'Bold'); % Add label to plot space

str = ['..\Figures\CorrelationTesting\Figure_6_F',experimentnumber];

savefigas = 'saveasPNG';
switch savefigas
    case 'saveasTIFF'; print(fig, str, '-dtiff', '-r1000');
    case 'saveasPNG'; saveas(fig, str, 'png');
    case 'saveasFIG'; savefig(str);
end
if closefigs; close(fig); end
savefigas = '';

%% error

% error per breath, nasal pred to flow pred
100*mean(abs(predy_flow-predy_pnasal))
[rho, pval] = corr(predy_flow, predy_pnasal);

% error per breath nasal pred to gold std
100*mean(abs(predy_flow-Yval))
[rho, pval] = corr(predy_flow, Yval); 

% error per pt median, nasal pred to gold std
100*mean(abs(PredY_avg-Gtest_avg))

if 0
%% Novel metrics, Median VE:Vdrive during sleep and Time with severe obstruction during sleep
% aka, Proportion of breaths FL
figure(27); clf(figure(27)); fig = gcf;
fig.Color = [1 1 1]; fig.Units = 'inches';
fig.Position = [-14   1   12    4.5];
%fig.Position = [0.5   0.5   12    4.5];

subplot(1,2,1);
[r_1, p_1] = plotregressionwithSEM(PredY_avg.*100, Gtest_avg.*100);
xlim([-5 110]); ylim([-5 110]);
ax = gca; ax.FontSize = 12; FntSz = 18;
ax.XTick=[0:25:125];
ax.YTick=[0:25:125];
xlabel('Pnasal Shape Predicted {\itflow:drive} (%)', 'FontSize', FntSz, 'FontName', 'Arial Narrow');
ylabel('Gold Standard {\itflow:drive} (%)', 'FontSize', FntSz, 'FontName', 'Arial Narrow'); 
axis square
str=['r = ', num2str(round(r_1,2))]; text(80, 7, str);
% titlestr = ['Patient Median {\itflow:drive}']; title(titlestr);

subplot(1, 2, 2);
[r_3, p_3] = plotregressionwithSEM(PredY_thres2_pFL.*100, Gtest_thres2_pFL.*100);
xlim([-5 110]); ylim([-5 110]); 
ax = gca; ax.FontSize = 12; FntSz = 18;
ax.XTick=[0:25:125]; 
ax.YTick=[0:25:125];
xlabel('Pnasal Shape Predicted % Flow Limited', 'FontSize', FntSz, 'FontName', 'Arial Narrow');
ylabel('Gold Standard % Flow Limited', 'FontSize', FntSz, 'FontName', 'Arial Narrow'); 
axis square
str=['r = ', num2str(round(r_3,2))]; text(80, 7, str);
% titlestr=['Proportion with Mod\Sev FL']; title(titlestr);

% subplot(1, 3, 3);
% [r_2, p_2] = plotregressionwithSEM(PredY_thres1_pFL.*100, Gtest_thres1_pFL.*100);
% xlim([-5 125]); xlabel('Pnasal Predicted %FL');
% ylim([-5 125]); ylabel('Gold Standard %FL');
% ax = gca; ax.XTick=[0:20:120];
% axis square
% str=['r = ', num2str(round(r_2,2))]; text(80, 7, str);
% titlestr=['Proportion with Severe FL']; title(titlestr);

% Add labels A B to plot space
subplot(1,2,1); hold on;
text(-45, 108, 'A', 'FontSize', 20, 'FontWeight', 'Bold');
subplot(1,2,2); hold on;
text(-45, 108, 'B', 'FontSize', 20, 'FontWeight', 'Bold');

str = ['..\Figures\', 'Figure_NovelMetrics'];
switch savefigas     
    case 'saveasTIFF'; print(fig, str, '-dtiff', '-r1000');     
    case 'saveasPNG'; saveas(fig, str, 'png');      
    case 'saveasFIG'; savefig(str);  
end
if closefigs; close(fig); end
end

%% After adjusting for AHI

% Median
[b,dev,stats]=glmfit([AHI_perPT_(:,1) PredY_avg],Gtest_avg);
p_1_adj = stats.p(3);
r_1_adj = stats.coeffcorr(1,3);

% mod\sev
[b,dev,stats]=glmfit([AHI_perPT_(:,1) PredY_thres2],Gtest_thres2);
p_2_adj = stats.p(3);
r_2_adj = stats.coeffcorr(1,3);

% severe
[b,dev,stats]=glmfit([AHI_perPT_(:,1) PredY_thres1],Gtest_thres1);
p_3_adj = stats.p(3);
r_3_adj = stats.coeffcorr(1,3);

if 0
%%  Figure E1 - flow VE  Vs pnasal VE
pnasal_veveup = PtData_pnasal_matched.VE ./ PtData_pnasal_matched.Veup;
flow_veveup = PtData_matched.VE ./ PtData_matched.Veup;
pnasal_veveup_t = pnasal_veveup.^1.5;

squareplot = 1;
sigdig = 2;
facealpha = 0.05;
facecolor = [0.1 0.1 0.1];

figure(28); clf(figure(28)); fig = gcf;
fig.Color = [1 1 1]; fig.Units = 'inches';
if squareplot
    fig.Position = [-11  1   8   8];
else
    fig.Position = [0.5  0.5   12  4];
end

%0.5
if squareplot; subplot(2,2,1); else; subplot(1,4,1); end
pnasal_veveup_05 = pnasal_veveup_t.^0.5;
scatter(pnasal_veveup_05, flow_veveup,2,'filled','markerfacealpha',facealpha, 'markerfacecolor', facecolor); hold on;
[r_1, p_1] =corr(pnasal_veveup_05, flow_veveup);
[rsq, r_dm] = UnivariateStats(pnasal_veveup_05, flow_veveup, weights);
mae_1 = CalcMAE(pnasal_veveup_05, flow_veveup);
xlim([-0.1 2.5]); ylim([-0.1 2.5]);
xlabel(['Pnasal VE/Veup']); ylabel(['Flow VE/Veup']);
text(1.7, 0.4, ['MAE = ', num2str(round(mae_1,sigdig))]);
text(1.7, 0.2, ['r = ', num2str(round(r_dm,sigdig))]);
text(1.7, 0.0, ['R^2 = ', num2str(round(rsq,sigdig))]);
%hndl_ls = lsline;  hndl_ls.LineStyle = '-'; hndl_ls.Color = [0 0 0];
hndl_ref = refline(1,0); hndl_ref.LineStyle = '--'; hndl_ref.Color = [0.25 0.25 0.25];
xlim([-0.1 2.5]); ylim([-0.1 2.5]);
%title(['0.5']);
axis square

%0.67
if squareplot; subplot(2,2,2); else; subplot(1,4,2); end
scatter(pnasal_veveup, flow_veveup,2,'filled','markerfacealpha',facealpha, 'markerfacecolor', facecolor); hold on;
[r_1, p_1] =corr(pnasal_veveup, flow_veveup);
[rsq, r_dm] = UnivariateStats(pnasal_veveup, flow_veveup, weights);
mae_1 = CalcMAE(pnasal_veveup, flow_veveup);
xlim([-0.1 2.5]); ylim([-0.1 2.5]);
xlabel(['Pnasal VE/Veup']); ylabel(['Flow VE/Veup']);
text(1.7, 0.4, ['MAE = ', num2str(round(mae_1,sigdig))]);
text(1.7, 0.2, ['r = ', num2str(round(r_dm,sigdig))]);
text(1.7, 0.0, ['R^2 = ', num2str(round(rsq,sigdig))]);
clear hndl_ls;
%hndl_ls = lsline;  hndl_ls.LineStyle = '-'; hndl_ls.Color = [0 0 0];
hndl_ref = refline(1,0); hndl_ref.LineStyle = '--'; hndl_ref.Color = [0.25 0.25 0.25];
%title(['0.67']);
axis square

%0.75
if squareplot; subplot(2,2,3); else; subplot(1,4,3); end
pnasal_veveup_075 = pnasal_veveup_t.^0.75;
scatter(pnasal_veveup_075, flow_veveup,2,'filled','markerfacealpha',facealpha, 'markerfacecolor', facecolor); hold on;
[r_1, p_1] =corr(pnasal_veveup_075, flow_veveup);
[rsq, r_dm] = UnivariateStats(pnasal_veveup_075, flow_veveup, weights);
mae_1 = CalcMAE(pnasal_veveup_075, flow_veveup);
xlim([-0.1 2.5]); ylim([-0.1 2.5]);
xlabel(['Pnasal VE/Veup']); ylabel(['Flow VE/Veup']);
text(1.7, 0.4, ['MAE = ', num2str(round(mae_1,sigdig))]);
text(1.7, 0.2, ['r = ', num2str(round(r_dm,sigdig))]);
text(1.7, 0.0, ['R^2 = ', num2str(round(rsq,sigdig))]);
%hndl_ls = lsline;  hndl_ls.LineStyle = '-'; hndl_ls.Color = [0 0 0];
hndl_ref = refline(1,0); hndl_ref.LineStyle = '--'; hndl_ref.Color = [0.25 0.25 0.25];
%title(['0.75']);
axis square

%1
if squareplot; subplot(2,2,4); else; subplot(1,4,4); end
pnasal_veveup_1 = pnasal_veveup_t.^1;
scatter(pnasal_veveup_1, flow_veveup,2,'filled','markerfacealpha',facealpha, 'markerfacecolor', facecolor); hold on;
[r_1, p_1] =corr(pnasal_veveup_1, flow_veveup);
[rsq, r_dm] = UnivariateStats(pnasal_veveup_1, flow_veveup, weights);
mae_1 = CalcMAE(pnasal_veveup_1, flow_veveup);
xlim([-0.1 2.5]); ylim([-0.1 2.5]);
xlabel(['Pnasal VE/Veup']); ylabel(['Flow VE/Veup']);
text(1.7, 0.4, ['MAE = ', num2str(round(mae_1,sigdig))]);
text(1.7, 0.2, ['r = ', num2str(round(r_dm,sigdig))]);
text(1.7, 0.0, ['R^2 = ', num2str(round(rsq,sigdig))]);
%hndl_ls = lsline;  hndl_ls.LineStyle = '-'; hndl_ls.Color = [0 0 0];
hndl_ref = refline(1,0); hndl_ref.LineStyle = '--'; hndl_ref.Color = [0.25 0.25 0.25];
%title([1']);
axis square

% Add labels A B to plot space
subplot(2,2,1); hold on;
% plot(1,1,'rs');
% currentYlim=ylim();
% currentXlim=xlim();
% plot([1,1], [0, currentYlim(2)],'r:', 'linewidth', 1);
% plot([0, currentXlim(2)],[1,1] ,'r:', 'linewidth', 1);
text(-0.8, 2.4, 'A', 'FontSize', 20, 'FontWeight', 'Bold');
subplot(2,2,2); hold on;
text(-0.8, 2.4, 'B', 'FontSize', 20, 'FontWeight', 'Bold');
subplot(2,2,3); hold on;
text(-0.8, 2.4, 'C', 'FontSize', 20, 'FontWeight', 'Bold');
subplot(2,2,4); hold on;
text(-0.8, 2.4, 'D', 'FontSize', 20, 'FontWeight', 'Bold');

%suptitle(['Pnasal VE/Veup Vs Flow VE/Veup']);
str = ['..\Figures\', 'Figure_E1_wExcl'];
switch savefigas     
    case 'saveasTIFF'; print(fig, str, '-dtiff', '-r1000');     
    case 'saveasPNG'; saveas(fig, str, 'png');      
    case 'saveasFIG'; savefig(str);  
end
if closefigs; close(fig); end
end


%% Histogram of ftrs selected during training 
ftr_array = NaN(54, ftrnum);
for subj=1:54 % set pt num or all, no data for PT=1
    if ~ismember(subj, pts_pnasal) %PT_list)
        continue
    end

    ftr_array(subj,:) = labels_Step_Subj_flow{subj,ftrnum}; % only one list of labels, same for both flow and pnasal

end
ftr_array = ftr_array(pts_pnasal,:);
uniqueftrs = unique(ftr_array);
ftr_array_linear = ftr_array(:);
ctrs = 1:1:size(Amatrix2,2);  % ToDo: make sure max ctrs covers all features
[counts, ~] = hist(ftr_array_linear, ctrs);
[~,I_uw] = sort(counts, 'descend');

figure(33); clf(figure(33)); fig = gcf;
fig.Color = [1 1 1]; fig.Units = 'inches';
fig.Position = [1   4   12    4.5];
bar(counts(I_uw), 'facealpha',0.2, 'edgealpha', 0.2); hold on;

xlim([0 40]);
ax = gca;
ax.TickDir = 'out';
ax.FontSize=9;
ax.XTick=[1:1:60];
xtickangle(ax,90);

LabelsOrdered_uw = Labels_Complete(I_uw);
lbls_uw = LabelsOrdered_uw(1:100);
lbls = regexprep(lbls_uw,'[_,:{}]','');
xticklabels(ax, lbls(1:25));

ylabel('Frequency'); 
box off
title(['Histogram of Features (simple count, 25 ftr model)']);
str = ['..\Figures\CorrelationTesting\', 'HistogramOfLinRegPnasalFeatures_SimpleCount',experimentnumber]; % 

savefigas = 'saveasPNG';
switch savefigas     
    case 'saveasTIFF'; print(fig, str, '-dtiff', '-r1000');     
    case 'saveasPNG'; saveas(fig, str, 'png');      
    case 'saveasFIG'; savefig(str);  
end
if closefigs; close(fig); end
savefigas = '';
% e.g. feature 3 is first 37 times, second 2 times. appears 39/41 times.
% while 15 never appears in top three, but is present 41 times.
%nnz(ftr_array(:,1)==3)
%nnz(ftr_array_linear==3)
%Labels_Complete(3)


%% Most important features using DM method
% Scores are based on 1/(x+1), x=rank per loop, higher scores are best.
% Features that regularly show up early (low rank e.g. First), get the biggest scores.
% Features that are occasionally very good, but sometimes late, are still considered.
% feature 3 is first 37 times, second 2 times. appears 39/41 times.
% while 15 never appears in top three, but is present 41 times.
% (side note, the variable 'If' is not my idea, this creation is SS...)
%nnz(ftr_array(:,1)==3)
%nnz(ftr_array_linear==3)
%Labels_Complete(3)

clearvars temp x If score score1 scoredata
score=zeros(1,size(Amatrix2,2));
maxscore=0; %if large score is good
%maxscore = (size(labels_Step_Subj,2))^2; %if small score is good

itersize = size(labels_Step_Subj_flow,1);

try
    for i=1:itersize
        score1 = maxscore + 0*score;
        If = [];
        for j=1:ftrnum %size(labels_Step_Subj,2)

            temp = labels_Step_Subj_flow{i,j};

            x=sum(temp'==If,2);
            temp(x==1)=[];
            if ~isempty(temp)
                If(end+1)=temp;
            end
        end
        for j=1:length(If)
            score1(If(j)) = 1./((j-1)+1);
            %score1(If(j)) = (j-1).^0.5;
            %score1(If(j)) = 1./(j-1).^0.5;
        end
        score=score+score1;
    end
    
    scoredata = [score;1:length(score)]';
    % scoredata = sortrows(scoredata,'descend'); % descend cmd only works for tble data
    scoredata = sortrows(scoredata,-1); % descending sort of col 1
    I_w=find(isnan(scoredata(:,1)));
    temp=scoredata(I_w,:);
    scoredata(I_w,:)=[];
    scoredata=[scoredata;temp];
    I_w = scoredata(:,2);
catch me
    disp(me.getReport);
end

% Modified histogram of ftrs selected during training
figure(34); clf(figure(34)); fig = gcf;
fig.Color = [1 1 1]; fig.Units = 'inches';
fig.Position = [2.2   2   12    4.5];
bar(scoredata(:,1), 'facealpha',0.2, 'edgealpha', 0.2);
xlim([0 40]);
ax = gca;
ax.TickDir = 'out';
ax.FontSize=9;
ax.XTick=[1:1:60];
xtickangle(ax,90);

LabelsOrdered_w = Labels_Complete(I_w);
lbls_w = LabelsOrdered_w(1:100);
lbls = regexprep(lbls_w,'[_,:{}]','');
xticklabels(ax,lbls(1:50));
ylabel('Weighted Score');
box off
title(['Histogram of Features (weighted count, 25 ftr model)']);
str = ['..\Figures\CorrelationTesting\', 'HistogramOfLinRegPnasalFeatures_WeightedCount',experimentnumber]; %
savefigas = 'saveasPNG';
switch savefigas     
    case 'saveasTIFF'; print(fig, str, '-dtiff', '-r1000');     
    case 'saveasPNG'; saveas(fig, str, 'png');      
    case 'saveasFIG'; savefig(str);  
end
if closefigs; close(fig); end
savefigas = '';

%% compare the two histogram ways of "top 50" ftrs
x_over = ismember(lbls_w(1:ftrnum), lbls_uw(1:ftrnum));
nnz(x_over) % this shows the number of features that occur in both lists
  
    
%% Alternatively, find best features through reverse stepwise regression
ftrnum_ = ftrnum; %6; 

ftr_array_2 = NaN(54, ftrnum_);
for subj=1:54 % set pt num or all, no data for PT=1
    if ~ismember(subj, pts_pnasal) %PT_list)
        continue
    end

   ftr_array_2(subj,:) = labels_Step_Subj_flow{subj,ftrnum_};

end
% Start with the list of unique ftrs selected during training at set stop point
ftr_array_2 = ftr_array_2(pts_pnasal,:);
uniqueftrs = unique(ftr_array_2);

% do backwards elimination
RemoveList = [];
while 1   
    [~,Pvals,~,~,~]=glmfitFast(Amatrix2_pnasal_matched(:,uniqueftrs),Gtest_matched,weights,1); 
    Pvals(1)=[];
    
    %remove least important
    [maxp,remi]=max(Pvals);

    disp(['Removing Ftr: ', num2str(uniqueftrs(remi)), ', p= ' num2str(maxp)]);
    RemoveList = [RemoveList; uniqueftrs(remi)];
    
    uniqueftrs(remi)=[];
   
    if isempty(uniqueftrs) || length(uniqueftrs)<1
        break
    end
end

% Remove list is currently in the order of first removed to last removed,
% so we really want those last removed, hence flip
TopElimFtrs = fliplr(RemoveList');


%% run model using only top 25 features
NfeaturesOpt=25;
if 0 % either of the histogram methods
    if 1 % unweighted histogram
        I = I_uw;
    else % weighted histogram
        I = I_w';
    end
    LabelsOrdered = Labels_Complete(I(1:NfeaturesOpt));
    TopFtrs = I(1:NfeaturesOpt);
else % else using backward elimination of unique ftrs at threshold
    LabelsOrdered = Labels_Complete(TopElimFtrs(1:NfeaturesOpt));
    TopFtrs = TopElimFtrs(1:NfeaturesOpt);
end


%% Final model using all data, uses selected N optimal features (NfeaturesOpt)
Ftrs = TopFtrs;
Labels = [{'Intercept'};LabelsOrdered];

FtrVals = Amatrix2_pnasal_matched(:,Ftrs);

tic
[Rsq,Pvals,RMSE,betas,FinalMdlPredY]=glmfitFast(FtrVals,Gtest_matched,weights,1); %faster
R_ = weightedcorrs([FinalMdlPredY,Gtest_matched],weights); R = R_(1,2);
mdl_pred_1 = [ones(length(weights),1) FtrVals]*betas;
mdl_pred_1(mdl_pred_1>maxG)=maxG; % set upper limit
mdl_pred_1(mdl_pred_1<0) = 0; % set lower limit
%Pvals(1)=[]; % remove intercept term
%betas(1)=[]; % remove intercept term
toc

tic
finalmdl = fitglm(FtrVals,Gtest_matched,'weights',weights);
mdl_pred_2 = predict(finalmdl, FtrVals);
mdl_pred_2(mdl_pred_2>maxG)=maxG; % set upper limit
mdl_pred_2(mdl_pred_2<0) = 0; % set lower limit
toc

%[finalmdl.Coefficients.Estimate, betas]
%[finalmdl.Coefficients.pValue, Pvals]

% test mdl_prediction against Gtest_
if 0
figure(4); clf(figure(4)); fig = gcf;
fig.Color = [1 1 1]; fig.Units = 'inches';
fig.Position = [0.5   5   12    4.5];
subplot(1,2,1);
scatter(mdl_pred_1,Gtest_matched,2,'filled','markerfacealpha',0.4); hold on
lsline();
ylabel('Gtest'); xlabel('Mdl Pred 1');
subplot(1,2,2);
scatter(mdl_pred_2,Gtest_matched,2,'filled','markerfacealpha',0.4); hold on
lsline()
ylabel('Gtest'); xlabel('Mdl Pred 2');
end

%finalmdl.Coefficients.tStat(2:end),
finalmdl_summary_mat = [[0;Ftrs'], betas, finalmdl.Coefficients.SE(1:end), Pvals];
finalmdl_summary_table = table([0;Ftrs'], Labels, betas, finalmdl.Coefficients.SE(1:end), Pvals, ...
        'VariableNames', {'FtrNum', 'FtrName', 'Betas', 'SE', 'Pval'});
  
str=['Final model r=', num2str(R),' ; Rsquared=', num2str(Rsq)]; disp(str)

summarystats = [r_1_6a, Rsq_6a, r_dm_6b, Rsq_6b, r_1_6c, Rsq_6c, RL1O_pnasal(ftrnum), RsqL1O_pnasal(ftrnum),  AccAestimated, AccBestimated, R, Rsq];
    
% Find univariate performance of 'best' features
rVsFD = NaN(length(Ftrs),1);
rVsPneumo = NaN(length(Ftrs),1);

r2VsFD = NaN(length(Ftrs),1);
r2VsPneumo = NaN(length(Ftrs),1);

bias =  NaN(length(Ftrs),1);

facealpha = 0.05; % was 0.08
facecolor = [0.2 0.2 0.2]; %was [0.1 0.1 0.1]

for ft=1:length(Ftrs)   
    % R vs flow:drive
    FtrVal = Amatrix2_pnasal_matched(:,Ftrs(ft));
    [Rsq,~,~,~]=glmfitFast(FtrVal,Gtest_matched,weights,1);
    r2VsFD(ft) = Rsq;
    R_ = weightedcorrs([FtrVal,Gtest_matched],weights);
    rVsFD(ft) = R_(1,2);
    
    % R vs equivalent feature measured with Pneumotach
    PneumoFtrVal = Amatrix2_flow_matched(:,Ftrs(ft));
    [Rsq,~,~,~]=glmfitFast(FtrVal,PneumoFtrVal,weights,1);
    r2VsPneumo(ft) = Rsq;
    
    R_ = weightedcorrs([FtrVal,PneumoFtrVal],weights);
    rVsPneumo(ft) = R_(1,2);
    
    % also calculate bias, as (median value Pneumotach) / (median value Pnasal)
    bias(ft) = median(PneumoFtrVal) / median(FtrVal);
    
    % figure
    if 0
       figure(101); clf(figure(101)); fig = gcf;
       fig.Color = [1 1 1]; fig.Units = 'inches';
       fig.Position = [-19   5  12   4.5];
       
       subplot(1,2,1); 
       set(gca,'box','off','tickdir','out','fontname','arial narrow');
       scatter(FtrVal,Gtest_matched,2,'filled','markerfacealpha', facealpha, 'markerfacecolor', facecolor); hold on;
       xlabel('Feature value (Pneumotach)'); ylabel('flow:drive');
       str=['R=',num2str(rVsFD(ft)), ',   R^2=',num2str(r2VsFD(ft))]
       title(str); axis square
       
       subplot(1,2,2);
       scatter(FtrVal,PneumoFtrVal,2,'filled','markerfacealpha', facealpha, 'markerfacecolor', facecolor); hold on;
       xlabel('Feature value (Pneumotach)'); ylabel('Feature value (Pnasal)');
       str=['R=',num2str(rVsPneumo(ft)), ',   R^2=',num2str(r2VsPneumo(ft))]
       title(str); axis square
       
       title_str = ['Ftr ', num2str(ft), ' -- ', LabelsOrdered{ft}];
       title_str = regexprep(title_str,'[_,:{\\}]',' ');
       suptitle(title_str);
       save_str = ['..\Figures\Univariate\', 'Pnasal_', title_str]; %
       saveas(fig, save_str, 'png');
    end
end

Univar_summary = table(rVsFD, rVsPneumo, bias, ...
        'VariableNames', {'RvsFlowDrive', 'RvsPneumo', 'Bias'});

Univar_summary_ = table([NaN;rVsFD], [NaN;rVsPneumo], [NaN;bias], ...
        'VariableNames', {'RvsFlowDrive', 'RvsPneumo', 'Bias'});

FinalModelDetails = [finalmdl_summary_table,Univar_summary_];


%%
experimentnumber % double check the experiment number
summarystats % this is set up to copy into excel "FeaturesList"

close all

% end here







%% Find Top X unique features based on scores
% if using transformed data, need to strip off the suffix
if TransformTheData; lengthsuf=4; else; lengthsuf=0; end
MaxUniquefeatures=50;
LabelsOrdered_NoSuffix=LabelsOrdered;
for i=1:length(LabelsOrdered_NoSuffix)
    LabelsOrdered_NoSuffix{i}=LabelsOrdered_NoSuffix{i}(1:length(LabelsOrdered_NoSuffix{i})-lengthsuf);
end
LabelsOrdered_NoSuffixUnique=unique(LabelsOrdered_NoSuffix,'stable');
LabelsOrdered_NoSuffixUnique(MaxUniquefeatures+1:end)=[];

%%










%%









%%
if 0
    %% Pred VEVdrive vs Actual VEVdrive, for Flow(all), Flow(matched) and Pnasal
    % plot is technically "Actual (x-axis) Vs Pred (y-axis)"
    figure(20); clf(figure(20));
    fig = gcf;
    fig.Color = [1 1 1];
    fig.Units = 'inches';
    fig.Position = [  -12.2    8.5   12    4.5];
    
    subplot(1,3,1);
    scatter(predyL1O_array_pnasal_matched(:,ftrnum), Gtest_matched, 2, 'filled','markerfacealpha',0.2);
    hold on; lsline; title('Pnasal'); ylabel('Actual VEVdrive'); xlabel('Pred VEVdrive');
    xlim([-0.05 1.55]); ylim([-0.05 1.55]); axis square;
    subplot(1,3,2);
    scatter(predyL1O_array_flow_matched(:,ftrnum), Gtest_matched, 2, 'filled','markerfacealpha',0.2);
    hold on; lsline; title('Flow (matched)'); xlabel('Pred VEVdrive');
    xlim([-0.05 1.55]); ylim([-0.05 1.55]); axis square;
    subplot(1,3,3);
    scatter(predyL1O_array(:,ftrnum), Gtest_All, 2, 'filled','markerfacealpha',0.2);
    hold on; lsline; title('Flow (All)'); xlabel('Pred VEVdrive');
    xlim([-0.05 1.55]); ylim([-0.05 1.55]); axis square;
    suptitle(['Actual Vs Pred VEVdrive at ', num2str(ftrnum), ' ftrs']);
    %savestr = [datastr,' ', supstr];
end

%% GS VEVdrive vs Pnasal VEVdrive (not pred VE:Vdrive estimate, just using pnasal to get VE)
figure(16); clf(figure(16)); fig = gcf;
%[r_t, p_t] = plotregressionwithSEM(PtData_matched.g_Edi_Adj, PtData_pnasal_matched.g_Edi_Adj);
scatter(PtData_matched.g_Edi_Adj, PtData_pnasal_matched.g_Edi_Adj,2,'filled','markerfacealpha',0.2); hold on;
[r,p] = corr(PtData_matched.g_Edi_Adj, PtData_pnasal_matched.g_Edi_Adj);
scatter(PT_avg_gtest, PT_avg_predy,20,'filled');
xlabel('Flow VEVdrive (Gold standard)');ylabel('Pnasal VEVdrive');
xlim([-0.05 1.55]);ylim([-0.05 1.55]);
titlestr=['Flow VEVdrive (Gold standard) Vs Pnasal VEVdrive']; title(titlestr);
str = ['..\Figures\', titlestr];
%saveas(fig, str, 'png'); %savefig(str);


%% Time below threshold - WITH/WITHOUT adding Apnoea and LowFlow breaths back in
%clinscorerange = unique(PtData.Etype) % clinical scoring
% clinical scoring codes
% 2 obstructive ap
% 3 central ap
% 4 hypopnoea ob
% 5 mixed

% This has the updated removed counts
% updatecounts = 0;
% if updatecounts
%     myfilename = ['C:\PSG_Data\FlowDrive\FeatureSpaces\FlowDrive_only25Hz_FeatureSpace_AutoRef2_Edi_Clean.mat'];
%     load(myfilename, 'RemovedBB_Apnoea', 'RemovedBB_LowFlow');
% end

ptstats = [];
pts_pnasal = unique(PtData_matched.PT);
sleeponly = 1;
addbackAP = 1;
thres = 0.5;
thres2 = 0.7;
numBBinTest = NaN(54,4);

% testing
if 0
    rn = randperm(10000);
    rn = rn/10000;
    idx = randperm(10000);
    actualBelowThres = rn(idx<5000)<thres;
    numActualBBbelowthreshold = nnz(actualBelowThres)
    
    figure(1); clf(figure(1));
    stairs([1:1:height(PtData_matched)], PtData_matched.Hypnog, 'b'); hold on;
    stairs([1:1:height(PtData_matched)], PtData_pnasal_matched.Hypnog, 'r');
    t = find(abs(PtData_matched.Hypnog-PtData_pnasal_matched.Hypnog)>0);
    plot(t, ones(length(t),1),'bx');
    nnz(PtData_matched.Hypnog==4)
end

for subj=1:54
    if ~ismember(subj, pts_pnasal)
        continue
    end
    if sleeponly % make Isubj
        Isubj=(PtData_matched.PT==subj)&(PtData_matched.Hypnog<4)&(PtData_matched.Ar==0);
        numBBinTest(subj,1) = nnz(Isubj);
        %Isubj=(PtData_matched.PT==subj)&(PtData_matched.Hypnog<4)&(PtData_matched.NotAr==1);
        %numBBinTest(subj,2) = nnz(Isubj);
    else
        Isubj=(PtData_matched.PT==subj);
        numBBinTest(subj,3) = nnz(Isubj);
        % Isubj=(PtData_matched.PT==subj)&(PtData_matched.Hypnog<4);
        % numBBinTest(subj,4) = nnz(Isubj);
    end
    
    clinscore = PtData_matched.Etype(Isubj);
    numclinscored = nnz(clinscore==2|clinscore==4); % obstructive and hypopnoea
    predBelowThres = predyL1O_array(Isubj,ftrnum)<thres;
    numPredBBbelowthreshold = nnz(predBelowThres);
    actualBelowThres = Gtest_matched(Isubj)<thres;
    numActualBBbelowthreshold = nnz(actualBelowThres);
    predBelowThres2 = predyL1O_array(Isubj,ftrnum)<thres2;
    numPredBBbelowthreshold2 = nnz(predBelowThres2);
    actualBelowThres2 = Gtest_All(Isubj)<thres2;
    numActualBBbelowthreshold2 = nnz(actualBelowThres2);
    numBBtotal = size(predyL1O_array(Isubj),1);
    if addbackAP
        withAPstr = [' (incl Ap-O breaths)'];
        indAP = find(RemovedBB_Apnoea.Pt == subj);
        indLF = find(RemovedBB_LowFlow.Pt == subj);
        if ~isempty(indAP)
            numclinscored = numclinscored+sum(RemovedBB_Apnoea{indAP,[4 6 7]});
            numPredBBbelowthreshold = numPredBBbelowthreshold+RemovedBB_Apnoea{indAP,2};
            numActualBBbelowthreshold = numActualBBbelowthreshold+RemovedBB_Apnoea{indAP,2};
            numPredBBbelowthreshold2 = numPredBBbelowthreshold2+RemovedBB_Apnoea{indAP,2};
            numActualBBbelowthreshold2 = numActualBBbelowthreshold2+RemovedBB_Apnoea{indAP,2};
            numBBtotal = numBBtotal+RemovedBB_Apnoea{indAP,2};
        end
        if ~isempty(indLF)
            numclinscored = numclinscored+sum(RemovedBB_LowFlow{indLF,[4 6 7]});
            numPredBBbelowthreshold = numPredBBbelowthreshold+RemovedBB_LowFlow{indLF,2};
            numActualBBbelowthreshold = numActualBBbelowthreshold+RemovedBB_LowFlow{indLF,2};
            numPredBBbelowthreshold2 = numPredBBbelowthreshold2+RemovedBB_LowFlow{indLF,2};
            numActualBBbelowthreshold2 = numActualBBbelowthreshold2+RemovedBB_LowFlow{indLF,2};
            numBBtotal = numBBtotal+RemovedBB_LowFlow{indLF,2};
        end
    else
        withAPstr = [''];
    end
    Clinpercent = 100*(numclinscored/numBBtotal);
    Predpercent = 100*(numPredBBbelowthreshold/numBBtotal);
    Actualpercent = 100*(numActualBBbelowthreshold/numBBtotal);
    Predpercent2 = 100*(numPredBBbelowthreshold2/numBBtotal);
    Actualpercent2 = 100*(numActualBBbelowthreshold2/numBBtotal);
    ptstats = [ptstats; [subj,numBBtotal,numclinscored,...
        numPredBBbelowthreshold,numActualBBbelowthreshold,...
        numPredBBbelowthreshold2,numActualBBbelowthreshold2,...
        Clinpercent,...
        Predpercent,Actualpercent,...
        Predpercent2,Actualpercent2]];
end
ptsummary_withAp = array2table(ptstats, 'VariableNames', {'PT','TotalBB','Clin_FL',...
    'Pred_FL','Actual_FL',...
    'Pred_FL2','Actual_FL2',...
    'Clin_percent',...
    'Pred_percent','Actual_percent',...
    'Pred_percent2','Actual_percent2'});
sum(ptsummary_withAp.TotalBB)
numBBinTest = numBBinTest(pts_pnasal,:);

%% average Gtest vs average PredY for each pt (medians...)
figure(24); clf(figure(24)); fig = gcf;
fig.Color = [1 1 1];
fig.Units = 'inches';
fig.Position = [ -12.2 -2.5  12.0000    4.5000];
pts_pnasal = unique(PtData_matched.PT); PT_listWpnasal = pts_pnasal;
subplot(1, 2, 1);
Gtest_avg = NaN(54,1);
PredY_avg = NaN(54,1);
for i=ftrnum % 1:size(predyL1O_array,2) say at 20 ftrs
    for subj=1:54 %length(PT_list)
        if ~ismember(subj, pts_pnasal)
            continue
        end
        %Isubj = PtData.PT==subj; % have not used weights(Isubj) ??
        if 0
            Isubj=(PtData_matched.PT==subj)&(PtData_matched.Hypnog<4)&(PtData_matched.Ar==0); % sleep only
        else
            Isubj=(PtData_matched.PT==subj)&(PtData_matched.NotAr==1);
        end
        Gtest_avg(subj) = median(Gtest_matched(Isubj));
        PredY_avg(subj) = median(predyL1O_array(Isubj,i));
        
        % make the data that includes the excluded Ap breaths
        indx = find(ptsummary_withAp.PT==subj);
        numPredFL = ptsummary_withAp.Pred_FL(indx);
        numActualFL = ptsummary_withAp.Actual_FL(indx);
        
        Gtest_avg_wAp(subj) = median([Gtest_matched(Isubj);ones(numActualFL,1)*0.1]);
        PredY_avg_wAp(subj) = median([predyL1O_array(Isubj,i);ones(numPredFL,1)*0.1]);
        
        scatter(predyL1O_array((Isubj),i),Gtest_matched(Isubj),2,'filled','markerfacealpha',0.5); hold on;
    end
end

% do once only.
Gtest_avg_wAp = Gtest_avg_wAp(pts_pnasal);
PredY_avg_wAp = PredY_avg_wAp(pts_pnasal);

scatter(PredY_avg, Gtest_avg, 50, 'filled','markerfacealpha',1); hold on;
% lsline;
xlim([0 1.5]); xlabel('PredY VE:VDrive');
ylim([0 1.5]); ylabel('Actual VE:VDrive');
title('(Sleep only)Indiv breaths and patient medians');
axis square

subplot(1, 2, 2);
[r_t, p_t] = plotregressionwithSEM(PredY_avg, Gtest_avg);
%scatter(PredY_avg, Gtest_avg, 50, 'filled'); hold on;
%hndl_ls = lsline;  hndl_ls.LineStyle = '-'; hndl_ls.Color = [0 0 1];
%hndl_ref = refline(1,0); hndl_ref.LineStyle = '--'; hndl_ref.Color = [0.5 0.5 0.5];
xlim([0 1.5]); xlabel('PredY VE:VDrive');
ylim([0 1.5]); ylabel('Actual VE:Vdrive');
axis square
[r, ~] = corr(PredY_avg(~isnan(PredY_avg)), Gtest_avg(~isnan(PredY_avg)));
titlestr = ['(Sleep only) Patient medians, r=', num2str(r)]; title(titlestr);

suptitle(['Actual VEVdrive Vs Pred VEVdrive']);
str = ['..\Figures\', 'MedianGtestVsMedianPredY_AutoRef2_SleepOnly_PnasalTnT'];
%saveas(fig, str, 'png'); %savefig(str);


%% median VEVdrive vs AHI
[AHI_perPT, AHI_perPT_table] = getAHI_postanalysis();
AHI_perPT_ = AHI_perPT(~isnan(AHI_perPT(:,1)),1);

excludelist = find(isnan(PredY_avg(:,1)));
AHI_perPT_ = AHI_perPT(:,1);
AHI_perPT_(excludelist) = [];
PredY_avg_ = PredY_avg;
PredY_avg_(excludelist) = [];

figure(25);clf(figure(25)); fig = gcf;
fig.Color = [1 1 1];
fig.Units = 'inches';
fig.Position = [ -12.2 8  12    4.5];

subplot(1,2,1)
%[r_t, p_t]=plotregressionwithSEM(PredY_avg, AHI_perPT(:,1));
scatter(PredY_avg, AHI_perPT(:,1), 50, 'filled'); hold on;
xlim([-0.05 1.55]); ylim([-5 105]); xlabel('Pred VE:VDrive');
ylabel('AHI'); title('AHI Vs Pred VEVdrive')
axis square

subplot(1,2,2)
scatter(Gtest_avg, AHI_perPT(:,1), 50, 'filled'); hold on;
xlim([-0.05 1.55]); ylim([-5 105]); xlabel('Actual VE:VDrive');
ylabel('AHI'); title('AHI Vs Actual VEVdrive')
axis square
suptitle(['AHI Vs Average VEVdrive']);

str = ['..\Figures\', 'AHI Vs Average VEVdrive_SleepOnly_PnasalTnT'];
%saveas(fig, str, 'png'); %savefig(str);

%% Proportion FL VS AHI
figure(26);clf(figure(26)); fig = gcf;
fig.Color = [1 1 1];
fig.Units = 'inches';
fig.Position = [1   2.5   12    4.5];

subplot(1,2,1)
%[r_t, p_t] = plotregressionwithSEM(ptsummary_withAp_n.Pred_percent, AHI_perPT(:,1));
scatter(ptsummary_withAp.Pred_percent2, AHI_perPT_(:,1),50,'filled','markerfacealpha',1);
% hndl_ls = lsline;  hndl_ls.LineStyle = '-'; hndl_ls.Color = [0 0 1];
% hndl_ref = refline(1,0); hndl_ref.LineStyle = '--'; hndl_ref.Color = [0.5 0.5 0.5];
xlim([-5 105]); xlabel('%FL (NumFL BB/TotalBB)');
ylim([-5 105]); ylabel('AHI');
title(['AHI Vs Pred %FL']);

subplot(1,2,2)
scatter(ptsummary_withAp.Actual_percent2, AHI_perPT_(:,1),50,'filled','markerfacealpha',1);
% hndl_ls = lsline;  hndl_ls.LineStyle = '-'; hndl_ls.Color = [0 0 1];
% hndl_ref = refline(1,0); hndl_ref.LineStyle = '--'; hndl_ref.Color = [0.5 0.5 0.5];
xlim([-5 105]); xlabel('%FL (NumFL BB/TotalBB)');
ylim([-5 105]); ylabel('AHI');
title(['AHI Vs Actual %FL']);

suptitlestr=['AHI Vs Proportion of sleep breaths classified as FL', withAPstr];
suptitle(suptitlestr);
str = ['..\Figures\', 'AHI Vs FL class_SleepOnly_PnasalTnT'];
%saveas(fig, str, 'png'); %savefig(str);

%% Novel metrics, Median VE:Vdrive during sleep and Time with severe obstruction during sleep
% aka, Proportion of breaths FL
figure(27); clf(figure(27)); fig = gcf;
fig.Color = [1 1 1];
fig.Units = 'inches';
%fig.Position = [-12.2   -3   12    4.5];

subplot(1,3,1);
if 1
    [r_1, p_1] = plotregressionwithSEM(PredY_avg_wAp, Gtest_avg_wAp);
else % without Ap BB
    [r_1, p_1] = plotregressionwithSEM(PredY_avg_, Gtest_avg_);
end
%hndl_ls = lsline;  hndl_ls.LineStyle = '-'; hndl_ls.Color = [0 0 1];
%hndl_ref = refline(1,0); hndl_ref.LineStyle = '--'; hndl_ref.Color = [0.5 0.5 0.5];
xlabel('Pnasal Predicted {\itflow:drive}');
ylabel('Gold Standard {\itflow:drive}');
axis square
titlestr = ['Patient Median {\itflow:drive}, r =', num2str(round(r_1,2))]; title(titlestr);

subplot(1, 3, 2);
if 1
    [r_2, p_2] = plotregressionwithSEM(ptsummary_withAp.Pred_percent, ptsummary_withAp.Actual_percent);
else % remove the one outlier
    [r_2, p_2] = plotregressionwithSEM(ptsummary_withAp.Pred_percent([1:3,5:17]), ptsummary_withAp.Actual_percent([1:3,5:17]));
end
%scatter(ptsummary_withAp.Pred_percent, ptsummary_withAp.Actual_percent,50,'filled','markerfacealpha',1);
%hndl_ls = lsline;  hndl_ls.LineStyle = '-'; hndl_ls.Color = [0 0 1];
%hndl_ref = refline(1,0); hndl_ref.LineStyle = '--'; hndl_ref.Color = [0.5 0.5 0.5];
xlim([-5 105]); xlabel('Pnasal Predicted %FL');
ylim([-5 105]); ylabel('Gold Standard %FL');
axis square
titlestr=['Proportion of FL breaths (0.5), r = ', num2str(round(r_2,2))]; title(titlestr);

subplot(1, 3, 3);
[r_3, p_3] = plotregressionwithSEM(ptsummary_withAp.Pred_percent2, ptsummary_withAp.Actual_percent2);
%scatter(ptsummary_withAp.Pred_percent, ptsummary_withAp.Actual_percent,50,'filled','markerfacealpha',1);
%hndl_ls = lsline;  hndl_ls.LineStyle = '-'; hndl_ls.Color = [0 0 1];
%hndl_ref = refline(1,0); hndl_ref.LineStyle = '--'; hndl_ref.Color = [0.5 0.5 0.5];
xlim([-5 105]); xlabel('Pnasal Predicted %FL');
ylim([-5 105]); ylabel('Gold Standard %FL');
axis square
titlestr=['Proportion of FL breaths (0.7), r = ', num2str(round(r_3,2))]; title(titlestr);

%suptitle(['Novel metrics (during sleep)']);

str = ['..\Figures\', 'Figure_6_']; %NovelMetrics_AutoRef2_SleepOnly_TnTPnasal'];
%saveas(fig, str, 'png'); %savefig(str);

if 0
    addpath('C:\Users\uqdmann\My Documents\MATLAB\random\BlandAltman');
    BlandAltman(PredY_avg, Gtest_avg);
end

%% AHI vs Novel Metrics
% grab from plot 25 and 26
figure(32);clf(figure(32)); fig = gcf;
fig.Color = [1 1 1];
fig.Units = 'inches';
%fig.Position = [ -12.2 8  12    4.5];
fig.Position = [ 2 3  12    4.5];

subplot(1,3,1)
[r_1, p_1]=plotregressionwithSEM(Gtest_avg_wAp', AHI_perPT_(:,1));
xlim([-0.05 1.55]); ylim([-5 105]); xlabel('Gold Standard {\itflow:drive}');
ylabel('AHI'); title(['Median {\itflow:drive} Vs AHI, r = -', num2str(round(r_1,2))]);
axis square

subplot(1,3,2)
%[r_2, p_2]=plotregressionwithSEM(ptsummary_withAp.Pred_percent, AHI_perPT_(:,1));
r_2 = corr(ptsummary_withAp.Actual_percent, AHI_perPT_(:,1));
scatter(ptsummary_withAp.Actual_percent, AHI_perPT_(:,1),50,'filled','markerfacealpha',1);
%hndl_ls = lsline;  hndl_ls.LineStyle = '-'; hndl_ls.Color = [0 0 1];
%hndl_ref = refline(1,0); hndl_ref.LineStyle = '--'; hndl_ref.Color = [0.5 0.5 0.5];
xlim([-5 105]); xlabel('%FL (NumFL BB/TotalBB)');
ylim([-5 105]); ylabel('AHI');
title(['AHI Vs Pred %FL']);
axis square

subplot(1,3,3)
[r_3, p_3]=plotregressionwithSEM(ptsummary_withAp.Actual_percent2, AHI_perPT_(:,1));
%r_2 = corr(ptsummary_withAp.Actual_percent, AHI_perPT_(:,1));
%scatter(ptsummary_withAp.Actual_percent, AHI_perPT_(:,1),50,'filled','markerfacealpha',1);
%hndl_ls = lsline;  hndl_ls.LineStyle = '-'; hndl_ls.Color = [0 0 1];
%hndl_ref = refline(1,0); hndl_ref.LineStyle = '--'; hndl_ref.Color = [0.5 0.5 0.5];
xlim([-5 105]); xlabel('Gold Standard %');
ylim([-5 105]); ylabel('AHI');
title(['Proportion of FL breaths (0.7) Vs AHI, r = ', num2str(round(r_3,2))]);
axis square

%suptitlestr=['AHI Vs Novel Metrics'];
%suptitle(suptitlestr);
str = ['..\Figures\', 'AHIvsNovelMetrics_AutoRef2_SleepOnly_PnasalTnT'];
%saveas(fig, str, 'png'); %savefig(str);

%% Proportion of breaths classified as FL (<threshold) Pred Vs Actual VEVdrive
if 0
    figure(28); clf(figure(28));
    plot([ones(length(ptsummary_withAp.Pred_percent),1), 2*ones(length(ptsummary_withAp.Actual_percent),1)], ...
        [ptsummary_withAp.Pred_percent, ptsummary_withAp.Actual_percent], 'bo'); hold on;
    plot([1 2], [ptsummary_withAp.Pred_percent, ptsummary_withAp.Actual_percent], 'color',[0.5 0.5 0.5], 'linewidth',0.5);
    plot([0.95 2.05], [nanmean(ptsummary_withAp.Pred_percent), nanmean(ptsummary_withAp.Actual_percent)], '--or', 'linewidth', 2.5);
    plot([0.95 2.05], [nanmedian(ptsummary_withAp.Pred_percent), nanmedian(ptsummary_withAp.Actual_percent)], '--om', 'linewidth', 2.5);
    xlim([0.75 2.25]); ylim([-5 105]);
    set(gca,'xtick',[1 2]);
    set(gca,'xticklabel',{'Pred', 'Actual'});
    ylabel('% of breaths FL'); xlabel('Method');
    titlestr = ['Proportion of sleep breaths classified as Flow Limited']; suptitle(titlestr);
    title('(Excluded Ap and LowFlow breaths not restored in this plot)');
    fig = gcf;
    fig.Color = [1 1 1]; % set background colour to white
    fig.Units = 'inches';
    %fig.Position = [-12 1 4.5 4.5];
    str = ['..\Figures\', titlestr, '_ManRef2_PredVsActual'];
    %saveas(fig, str, 'png'); %savefig(str);
    
    figure(29); clf(figure(29));
    plot([0.5*ones(length(ptsummary_withAp.Pred_percent),1)], ...
        [(ptsummary_withAp.Pred_percent - ptsummary_withAp.Actual_percent)], 'bo'); hold on;
    boxplot([(ptsummary_withAp.Pred_percent - ptsummary_withAp.Actual_percent)],...
        [(ones(length(ptsummary_withAp.Pred_percent),1))])
    
    xlim([0 1.5]);
    set(gca,'xtick',[0.5 1]);
    set(gca,'xticklabel',{'Indiv delta', 'summary box'});
    ylabel('Predicted % - Actual %');
    suptitle('Difference between Pred and Actual');
    fig = gcf;
    fig.Color = [1 1 1]; % set background colour to white
    fig.Units = 'inches';
    %fig.Position = [-12 1 4.5 4.5];
    str = ['..\Figures\', titlestr, '_ManRef2_PredVsActual_boxplot'];
    %saveas(fig, str, 'png'); %savefig(str);
    
    figure(30); clf(figure(30));
    scatter(ptsummary_withAp.Pred_FL, ptsummary_withAp.Actual_FL,50,'filled','markerfacealpha',1);
    hndl_ls = lsline;  hndl_ls.LineStyle = '-'; hndl_ls.Color = [0 0 1];
    hndl_ref = refline(1,0); hndl_ref.LineStyle = '--'; hndl_ref.Color = [0.5 0.5 0.5];
    xlim([0 4000]); xlabel('PredY, Number of FL breaths');
    ylim([0 4000]); ylabel('Actual, Number of FL breaths');
    suptitlestr=['Number of sleep breaths classified as FL'];
    suptitle(suptitlestr);
    [r,~] = corr(ptsummary_withAp.PredNumFL, ptsummary_withAp.ActualNumFL);
    titlestr=['FL if VE:Vdrive < ',num2str(thres), ', r= ', num2str(r)];
    title(titlestr);
    savestr = ['..\Figures\', suptitlestr, '_ManRef2_Scatter'];
    fig = gcf;
    %saveas(fig, savestr, 'png'); %savefig(str);
end


%% Median VEVdrive VS Proportion FL
Gtest_avg_ = Gtest_avg(~isnan(Gtest_avg));
PredY_avg_ = PredY_avg(~isnan(PredY_avg));
figure(31); clf(figure(31)); fig = gcf;
fig.Color = [1 1 1];
fig.Units = 'inches';
fig.Position = [-12.2   -3   12    4.5];

subplot(1,2,1)
%[r_t, p_t] = plotregressionwithSEM(PredY_avg_, ptsummary_withAp.Pred_percent);
scatter(PredY_avg_, ptsummary_withAp.Pred_percent,50,'filled','markerfacealpha',1);
hndl_ls = lsline;  hndl_ls.LineStyle = '-'; hndl_ls.Color = [0 0 1];
%hndl_ref = refline(1,0); hndl_ref.LineStyle = '--'; hndl_ref.Color = [0.5 0.5 0.5];
xlim([-0.05 1.55]); xlabel('VEVdrive');
ylim([-5 105]); ylabel('%FL (NumFL BB/TotalBB)');
titlestr=['Pred %FL Vs Pred VEVdrive']; title(titlestr);

subplot(1,2,2)
scatter(Gtest_avg_, ptsummary_withAp.Actual_percent,50,'filled','markerfacealpha',1);
hndl_ls = lsline;  hndl_ls.LineStyle = '-'; hndl_ls.Color = [0 0 1];
%hndl_ref = refline(1,0); hndl_ref.LineStyle = '--'; hndl_ref.Color = [0.5 0.5 0.5];
xlim([-0.05 1.55]); xlabel('VEVdrive');
ylim([-5 105]); ylabel('%FL (NumFL BB/TotalBB)');
titlestr=['Actual %FL Vs Actual VEVdrive']; title(titlestr);

suptitlestr=['%FL Vs Average VEVdrive'];
suptitle(suptitlestr);


%%
if 0
    figure(33); clf(figure(33));
    plot([ones(length(ptsummary_withAp.Pred_percent),1), 2*ones(length(ptsummary_withAp.Clin_percent),1)], ...
        [ptsummary_withAp.Pred_percent, ptsummary_withAp.Clin_percent], 'bo'); hold on;
    plot([1 2], [ptsummary_withAp.Pred_percent, ptsummary_withAp.Clin_percent], 'color',[0.5 0.5 0.5], 'linewidth',0.5);
    plot([0.95 2.05], [nanmean(ptsummary_withAp.Pred_percent), nanmean(ptsummary_withAp.Clin_percent)], '--or', 'linewidth', 2.5);
    plot([0.95 2.05], [nanmedian(ptsummary_withAp.Pred_percent), nanmedian(ptsummary_withAp.Clin_percent)], '--om', 'linewidth', 2.5);
    xlim([0.75 2.25]); ylim([-5 105]);
    set(gca,'xtick',[1 2]);
    set(gca,'xticklabel',{'Pred FL', 'Clinical'});
    ylabel('% of breaths FL'); xlabel('Method');
    titlestr = ['Proportion of sleep breaths classified as Flow Limited']; suptitle(titlestr);
    title('(Excluded Ap and LowFlow breaths are restored in this plot)');
    fig = gcf;
    fig.Color = [1 1 1]; % set background colour to white
    fig.Units = 'inches';
    %fig.Position = [-12 1 4.5 4.5];
    str = ['..\Figures\', titlestr, '_ManRef2_withAddins'];
    %saveas(fig, str, 'png'); %savefig(str);
end


%% Load subject waveform data for VEVdrive, Flow and Edi plot
predy_flow(predy_flow<0)=0;
predy_flow(predy_flow>maxG)=1.5;

% read spreadsheet
AnalyzeDataSpreadsheet = 'AnalyzeDataSpreadsheet.xlsx';
[num,patients,~] = xlsread(AnalyzeDataSpreadsheet,1,'F3:G56');
for n=1 %41:54
    if ~ismember(n, pts_pnasal)
        continue
    end
    clearvars Edi Flow StarttimeSpike
    % Plot FL values over time for an example subject (e.g. compare with Spike data). Random check (n=8, 1313) Looks great.
    % n=3; % 1313 is now 14
    % [2;3;4;5;6;8;9;10;11;14;16;17;18;19;20;21;22;23;24;26;27;28;29;30;32;33;34;35;36;38;39;40;41;43;44;45;46;47;50;53;54]
    studyname = char(patients{n});
    %load('J:\PEOPLE\FACULTY\SANDS\Phenotyping from Pes and Flow\SpikeFiles\1313.mat', 'Edi', 'Flow');
    load(['C:\PSG_Data\FlowDrive\SourceMat 20171123\',studyname], 'Edi', 'Flow', 'StarttimeSpike');
    if ~(exist('StarttimeSpike', 'var') == 1)
        StarttimeSpike = 0;
    end
    Time = StarttimeSpike:0.008:StarttimeSpike+(length(Flow.values)-1)*0.008;
    
    PTtime = PtData.BB_time(PtData.PT==n);
    PTtime = PTtime-PTtime(1)+StarttimeSpike;
    Data1 = [PtData.BB_time(PtData.PT==n) PtData.BB_Ttot(PtData.PT==n) predy_flow(PtData.PT==n) Gtest_All(PtData.PT==n)];
    addNaNgaps=1;
    cols = [NaN 1 2];
    if addNaNgaps
        tol2=0.1;
        i=1;
        M=size(Data1,2);
        while i<(size(Data1,1)-1)
            if (Data1(i,cols(2))+Data1(i,cols(3))+tol2)<Data1(i+1,cols(2))
                Data1 = [Data1(1:i,:); NaN*ones(1,M); Data1((i+1):size(Data1,1),:)];
                %keyboard
                i=i+1;
            end
            i=i+1;
        end
    end
    % if Flow.length<=length(Time)
    %     Flow.values(end:length(Time))=NaN;
    % else
    %     Flow.values(length(Time):end)=[];
    % end
    if Edi.length<=length(Time)
        Edi.values(end:length(Time))=NaN;
    else
        Edi.values(length(Time):end)=[];
    end
    
    % Data columns: (1) Time, (2) BBTtot, (3) PredY, (4) Gtest.
    
    %% make the plot for the subject data
    figure(200+n); clf(200+n) %ax(1) = subplot(2,1,1);
    ax(1)=subplot(3,1,1);
    stairs(Data1(:,1),100*Data1(:,4),'k'); % actual
    hold('on')
    stairs(Data1(:,1),100*Data1(:,3),'r'); % pred
    plot([Data1(1,1) Data1(end,1)],[0 0],'k:');
    plot([Data1(1,1) Data1(end,1)],100*[1 1],'k:');
    plot([Data1(1,1) Data1(end,1)],50*[1 1],'k:');
    set(gca,'xtick',[],'box','off');
    ylabel('{\itflow:drive} (%)'); ylim([-5 105]);
    
    dsf=5; dt=Flow.interval;
    FlowF=Flow.values;
    if 1
        filter_HFcutoff_butter0 = 12.5;
        filter_order0 = 1;
        [B_butter0,A_butter0] = butter(filter_order0,[filter_HFcutoff_butter0]/(1/dt/2),'low');
        FlowF = filtfilt(B_butter0,A_butter0,FlowF); %filtfilt, otherwise flow signal is right-shifted
    end
    %ax(2)=subplot(3,1,2); plot(Time,Flow.values,'g'); hold('on');
    ax(2)=subplot(3,1,2); plot(downsample(Time,dsf),downsample(FlowF,dsf),'k');
    set(gca,'xtick',[],'box','off');
    ylabel('Flow (L/s)');
    
    ax(3)=subplot(3,1,3); plot(downsample(Time,dsf),downsample(Edi.values,dsf)); hold on;
    box('off');
    ylabel('Edi (uV)');
    
    linkaxes(ax,'x');
    suptitle([studyname]);
end

%% mods to plot
plot(StarttimeSpike+5400, 0, 'r^');
plot(StarttimeSpike+23620, 0, 'r^');

xlim([StarttimeSpike inf])
xlim([StarttimeSpike+5200 StarttimeSpike+6400]) % old indices
xlim([StarttimeSpike+5400 StarttimeSpike+5800]) % old indices
xlim([StarttimeSpike+23560 StarttimeSpike+23680]) % old indices

%%



