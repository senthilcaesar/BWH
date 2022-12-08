%% start
keyboard % this is just to stop the whole process running accidentally (i.e. pressing F5)
close all
clear
clc
% goto line 780 for post-processing
if 1
    addpath('C:\Users\uqdmann\Dropbox\PUPbeta_git\');
    cd('C:\Users\uqdmann\Dropbox\PUPbeta_git\FlowAndEdiPesQAO');
    %load('C:\Users\uqdmann\Dropbox\PUPbeta_git\FeatureSpaces\FlowDrive_only25Hz_FeatureSpace_AutoRef2_Edi_Clean_withPnasalMatched_LinRegWorkspace_20180122.mat')
    % with Pnasal matched is using flow beta in pnasal data
    % starts at ~ line 760
end

%% options
Use25HzDS = 1;          % set as 1 for 25Hz data, or 0 to use 125Hz data
EdiDrive = 1;           % set as 1 for Edi Drive data, 0 for Pes drive
useLogitFtrs = 1;       % set as 1 to use ftrs selected from logit, 0 for all ftrs
SS_fastfit = 1;			% set as 1 to use SS fast glmfit method, 0 for matlab builtin
TransformTheData = 1;   % set as 1 to do tranforms, or 0 to use unadjusted data
addextratransform = 0; 	% set as 1 to do extra transforms
ShowFigures = 0;        % set as 1 to show figures, or 0 to not show figures
IncludePnasal = 1;      % set as 1 to include Pnasal, or 0 to ignore pnasal
% pnasal testing in this code is using the flow betas

%datadir = 'C:\PSG_Data\FlowDrive\FeatureSpaces\';
%datadir = '..\FeatureSpaces\';  %
datadir = 'C:\PSG_Data\FlowDrive\Analyzed\';


%% Turn diary logging on
HowSoonIsNow = datestr(datetime('now'),'yyyymmdd_HHMM');
diaryfilename = ['L1O_LinReg_FD_', HowSoonIsNow, '.txt'];
diary(diaryfilename);
diary on

%% open file
if Use25HzDS
    if EdiDrive
        % filename = [datadir, 'FlowDrive_only25Hz_ReNormalized_FeatureSpace_AutoRef2_Edi_Clean.mat']; datastr = 'FD Normalized edi 25Hz';
        %filename = [datadir, 'FlowDrive_only25Hz_FeatureSpace_AutoRef2_Edi_Clean.mat']; datastr = 'FD edi 25Hz';
        filename = [datadir, 'FlowDrive_only25Hz_ReNormalized_FeatureSpace_AutoRef2_Edi_Clean.mat']; datastr = 'FlowandPnasalD edi 25Hz';
        if IncludePnasal % for Pnasal
            %pnasal_filename = [datadir, 'PnasalDrive_only25Hz_exp067_FeatureSpace_AutoRef2_Edi_Clean.mat']; datastr = 'FlowandPnasalD edi 25Hz';
            pnasal_filename = [datadir, 'PnasalDrive_only25Hz_ReNormalized_FeatureSpace_AutoRef2_Edi_Clean.mat']; datastr = 'FlowandPnasalD edi 25Hz';
        end
    else; filename = [datadir, 'FlowDrive_only25Hz_FeatureSpace_Pes_Clean.mat']; datastr = 'FD pes 25Hz';
    end
else
    if EdiDrive; filename = [datadir, 'FlowDrive_only125Hz_FeatureSpace_Edi_Clean.mat']; datastr = 'FD edi 125Hz';
    else; filename = [datadir, 'FlowDrive_only125Hz_FeatureSpace_Pes_Clean.mat']; datastr = 'FD pes 125Hz';
    end
end

try
    if IncludePnasal
        str=['Loading ' pnasal_filename]; disp(str);
        load(pnasal_filename,'Amatrix', 'PtData', 'FeatureNames');
        Amatrix_pnasal = Amatrix;
        PtData_pnasal = PtData;
        FeatureNames_pnasal = FeatureNames;
    end
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

if IncludePnasal
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
end

%% use only a given set of features, determined from logistic regression
% keep the original Amatrix (before selecting specific features, and transforms
% keep the origianl FeatureNames variable before overwriting
% same for pnasal, if included
Original_FeatureNames = FeatureNames;
Original_Amatrix = Amatrix;
if useLogitFtrs
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
    
    if IncludePnasal
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
    end
else
    if IncludePnasal
        % do a simple check, so that we at least start with the same features
        [~,ia,ib] = intersect(FeatureNames_pnasal.Name, FeatureNames.Name);
        FeatureNames_pnasal = FeatureNames_pnasal(ia,:);
        Amatrix_pnasal = Amatrix_pnasal(:,ia);
        FeatureNames = FeatureNames(ib,:);
        Amatrix = Amatrix(:,ib);
    end
end

%% Remove features in "FtrsToExclude" list, using string search
if 0 % turned off, AA_PeaksRatio has been removed earlier
    % FtrsToExclude = {...
    %     'InspFlutPow_Vpeak2_21_O',...
    %     'ExpFlutPow_VpeakE2_22_O',...
    %     'InspFlutPow_Vpeak2_20_T',...
    %     'ExpFlutPow_VpeakE2_21_T',...
    %     'InspExpFlutPow_22_T',...
    %     'AA_PeaksRatio_25_T',...
    %     'MedianFlowChange_26_T',...
    %     'STDFlowChange_27_T',...
    %     'AA_PeaksRatio_24_O',...
    %     'MinFlowChange_25_O',...
    %     'MedianFlowChange_26_O'...
    %     }
    FtrsToExclude = {...
        'AA_PeaksRatio',...
        };
    Ind = [];
    for i=1:length(FtrsToExclude) %needs checking
        temp=find(startsWith(FeatureNames.Name,FtrsToExclude(i)));
        if ~isempty(temp)
            Ind(end+1)=temp;
        end
    end
    FeatureNames(Ind,:)=[];
    Amatrix(:,Ind)=[];
end

%% setup the data matrix to use, and do the transforms
if TransformTheData
    [Amatrix2, Labels] = DoDataMatTransform(Amatrix, FeatureNames, addextratransform);
    % do it again for pnasal
    if IncludePnasal
        [Amatrix2_pnasal, ~] = DoDataMatTransform(Amatrix_pnasal, FeatureNames_pnasal, addextratransform);
    end
else
    Amatrix2 = [Amatrix];
    Labels = FeatureNames.Name;
    if IncludePnasal
        Amatrix2_pnasal = Amatrix_pnasal;
    end
end

%% test feature duration scaling
% note this only works with full feature set.
if 0
    close; clear all; clc;
    datadir = 'C:\PSG_Data\FlowDrive\FeatureSpaces\';
    load([datadir,'FlowDrive_only25Hz_Normalized_FeatureSpace_AutoRef2_Edi_CleanR.mat']);
    
    unique(PtData.PT)
    % Amatrix(:,1:79) should be equivalent to Amatrix(:,159:end)
    %Amatrix = BreathFLDataTable{:,:};
    mat1 = Amatrix(:,1:88);     % standard duration
    mat2 = Amatrix(:,177:end);  % double duration
    mat3 = Amatrix(:,89:176);   % ttran timing
    
    % correlations
    corrVals = NaN(88,1);
    for ftr=1:88
        corrVals(ftr) = corr(mat1(:,ftr), mat2(:,ftr));
    end
    oddities = find(corrVals<=0.90);
    oddlist = [num2cell((oddities)) , num2cell(corrVals(oddities)), FeatureNames.Name(oddities)];
    oddTable = cell2table(oddlist, 'VariableNames', {'Ftr', 'Correlation', 'FtrName'}); oddTable
    
    % absolute deltas, and percent change
    mat_diff = (mat1 - mat2) ./ abs(mat1);
    percentchange = abs(nanmean(mat_diff,1));
    odds2 = find(percentchange>0.05);
    
    oddlist2 = [num2cell((odds2))', num2cell(corrVals(odds2)), num2cell(percentchange(odds2))', FeatureNames.Name(odds2)];
    oddTable2 = cell2table(oddlist2, 'VariableNames', {'Ftr', 'Correlation', 'Change', 'FtrName'}); oddTable2
    
    ftr = 50;
    figure(1); clf(figure(1));
    scatter(mat1(:,ftr), mat2(:,ftr), 3,'filled','markerfacealpha',0.3);
    xlabel('normal duration'); ylabel('double duration');
    
end

%% Preparation (Do this carefully)
Gtest_All = PtData.g_Edi_Adj;
if IncludePnasal
    Gtest_All_pnasal = PtData_pnasal.g_Edi_Adj; % potentially unused, as we use Flow based Gtest
end

%NfeatureSteps=100;
colofones = ones(length(Gtest_All),1);
%predyL1O = NaN*Gtest_All;
%predyL1O_array = NaN*ones(length(Gtest_All),NfeatureSteps);

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
if ~useweights %overwrite, i.e. do not use weights, just equal for all classes
    weights = ones(length(weights),1);
end

%% Leave one subject out loop
% lin reg leave one out
% FD, 41 Subjects, 135k breaths, 50 ftrs, 2x transforms, takes 30 min to 1 hr

% use glmfit to do generalized linear model regression
% use fitglm to create a generalized linear regression model (logistic?)
try
    clear RvalueTrain Err ErrRms Nfeatures
    clear RsqTrain_array ErrTrain_array ErrRmsTrain_array
    labels_Step_Subj =[];
    alpha=0.05;
    MaxNfeatures=1;
    neverbreak=1;
    t_start_L1O = clock;
    PT_list = unique(PtData.PT);
    if IncludePnasal
        %Pnasal_list = [3 5 8 9 18 22 33 34 35 36 38 43 46 47 50 53 54];
        %Pnasal_list = [3 8 9 18 22 35 36 38 43 46 47 50 53 54];
        Pnasal_list = unique(PtData_pnasal.PT); % dynamic, better
    end
    Labels_Complete = Labels; % save a backup of the complete Label set
    Fwd=0;
    VEVeup_All = [];
    if IncludePnasal
        VEVeup_pnasal_array = [];
        VEVeup_flow_array = [];
        Pnasal_summary = [];
        PnasalFlowBB_PTindex=[];
        PtData_matched = [];
        PtData_pnasal_matched = [];
        Gtest_matched = [];
        weights_matched = [];
        predyL1O_array_pnasal_matched = [];
        predyL1O_array_flow_matched = [];
        Amatrix2_flow_matched = [];
        Amatrix2_pnasal_matched = [];
        ChannelsList = {'Flow','Pnasal'};
        artdirectory = ['C:\PSG_Data\FlowDrive\SourceMat 20171123'];
        % read spreadsheet (options worksheet)
        [~,~,raw] = xlsread('C:\Users\uqdmann\Dropbox\PUPbeta_git\FlowAndEdiPesQAO\AnalyzeDataSpreadsheet.xlsx',1,'G3:G56');
    end
    
    for subj=1:54 %subj=2 subj=3
        %% loop through all subjects, leaving one out each time
        if ~ismember(subj, PT_list)
            str=['No data for Pt ', num2str(subj)]; disp(str); continue
        end
        
        disp(' '); str=['Performing analysis, witholding Pt ', num2str(subj)]; disp(str); tic;
        
        PtHasPnasal = 0;
        if IncludePnasal
            if ismember(subj, Pnasal_list)
                str=['Also analysing pnasal data for this pt']; disp(str);
                PtHasPnasal = 1;
            end
        end
        
        %% set up the full length Flow data
        Labels = Labels_Complete; % Start with full list
        Isubj=(PtData.PT==subj);  % Isubj is the logical index of the L1O patient
        
        % set up training data, i.e. all but the left out pt
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
        
        % set up test data, i.e. the left out pt
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
        % in doing this, we also produce a second set of flow data,
        % but these are only those matched with the pnasal breaths
        if PtHasPnasal
            Isubj_pnasal=(PtData_pnasal.PT==subj);
            PnasalFlowBB_PTindex_pt=[];
            
            predyL1O_array_pnasal_pt = [];
            predyL1O_array_flow_pt = [];
            VEVeup_pnasal_array_pt =[];
            VEVeup_flow_array_pt =[];
            Gtest_test_flow = Gtest_test;
            weights_test_flow = weights_test;
            
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
                    if nnz(FlowBBtoExclude)~=0 || nnz(PnasalBBtoExclude)~=0
                        str=['Removing breaths (Flow, Pnasal): ', num2str(nnz(FlowBBtoExclude)),', ',num2str(nnz(PnasalBBtoExclude))]; disp(str);
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
        end
        
        
        
        %% redo weights for different length data
        % all breaths (unchanged number of breaths)
        % pnasal matched breaths (reduced number of breaths)
        
        
        %% redo
        NfeatureSteps=100;
        colofones = ones(length(Gtest_All),1);
        predyL1O = NaN*Gtest_All;
        predyL1O_array = NaN*ones(length(Gtest_All),NfeatureSteps);
        predyL1O_array_pnasal = NaN*ones(length(Gtest_matched),NfeatureSteps);
        
        %% do backwards elimination
        warning('off');
        maxp=Inf;
        
        %% remove features that do not correlate between flow and pnasal
        Ftr_indx = 1:size(Amatrix2_train,2);

        %% do backwards elimination
        while length(Ftr_indx)>=MaxNfeatures || maxp>alpha
            
            if SS_fastfit % nnz(~isfinite(Amatrix2_))
                [~,Pvals_train,~,b_train,Ypred_train]=glmfitFast(Amatrix2_train(:,Ftr_indx),Gtest_train,weights_train,1); %faster
                Pvals_train(1)=[];
            else % matlab built-in
                [b_train,dev,stats]=glmfit(Amatrix2_train(:,Ftr_indx),Gtest_train,'normal','weights',weights_train);
                Pvals_train = stats.p(2:end);
            end
            
            if length(Ftr_indx)<=NfeatureSteps %start saving results
                %% save the training data
                predytrain = [colofones_train Amatrix2_train(:,Ftr_indx)]*b_train;
                predytrain(predytrain>maxG)=maxG;
                predytrain(predytrain<0)=0; %will be overwritten
               
                    % add in R at this step.
                    [RsqTrain_array(subj,length(Ftr_indx)), RTrain_array(subj,length(Ftr_indx))] = ...
                        UnivariateStats(predytrain,Gtest_train, weights_train); % new method for R and Rsq
                    %RsqTrain_array(subj,length(Ftr_indx)) = 1-nansum(weights_train.*(Gtest_train-predytrain).^2)/nansum(weights_train.*(Gtest_train-nanmean(Gtest_train)).^2);
                    ErrTrain_array(subj,length(Ftr_indx)) = nanmean(weights_train.*abs(predytrain-Gtest_train));
                    ErrRmsTrain_array(subj,length(Ftr_indx)) = nanmean((weights_train.*(predytrain-Gtest_train)).^2).^0.5;
           
                
                %% use beta from training, apply to test, save data
                predyL1O_ = [colofones_test Amatrix2_test(:,Ftr_indx)]*b_train;
                predyL1O_(predyL1O_>maxG)=maxG;
                predyL1O_(predyL1O_<0) = 0;
                predyL1O_array(Isubj,length(Ftr_indx)) = predyL1O_;
                
                if 0 % do at end
                    RsqL10_array(subj,length(Ftr_indx)) = 1-...
                        nansum(weights_test.*(Gtest_test-predyL1O_array(Isubj,length(Ftr_indx))).^2)...
                        /nansum(weights_test.*(Gtest_test-nanmean(Gtest_test)).^2);
                    
                    if (RsqL10_array(subj,length(Ftr_indx)) < 0) % preserve sign
                        RL10_array(subj,length(Ftr_indx)) = -1*(abs(RsqL10_array(subj,length(Ftr_indx))).^0.5);
                    else
                        RL10_array(subj,length(Ftr_indx)) = (RsqL10_array(subj,length(Ftr_indx))).^0.5;
                    end
                    ErrL1O_array(subj,length(Ftr_indx)) =  nanmean(weights_test.*abs(predyL1O_-Gtest_test));
                    ErrRmsL1O_array(subj,length(Ftr_indx)) = nanmean((weights_test.*(predyL1O_-Gtest_test)).^2).^0.5;
                end
                
                %% store the labels and beta array
                labels_Step_Subj{subj,length(Ftr_indx)} = Ftr_indx;
                beta_array{subj,length(Ftr_indx)} = b_train;
                
                %% how do the flow betas perform in Pnasal features?
                if PtHasPnasal
                    % some intitial set up
                    PnasalFlowBB_PTindex_pt = subj*ones(length(weights_test_flow),1);
                    
                    predyL1O_pnasal = [ones(length(weights_test_flow),1), Amatrix2_pnasal_pt(:, Ftr_indx)]*b_train;
                    predyL1O_pnasal(predyL1O_pnasal>maxG)=maxG;
                    predyL1O_pnasal(predyL1O_pnasal<0) = 0;
                    
                    predyL1O_flow = [ones(length(weights_test_flow),1), Amatrix2_flow_pt(:, Ftr_indx)]*b_train;
                    predyL1O_flow(predyL1O_flow>maxG)=maxG;
                    predyL1O_flow(predyL1O_flow<0) = 0;
                    
                    if 0
                        % test for (near) equivalence
                        figure(1); clf(figure(1));
                        subplot(1,3,1)
                        scatter(predyL1O_pnasal, Gtest_test_flow, 3);
                        hold on; lsline; title('Pnasal');
                        subplot(1,3,2)
                        scatter(predyL1O_flow, Gtest_test_flow, 3);
                        hold on; lsline; title('Flow (matched)');
                        subplot(1,3,3)
                        scatter(predyL1O_, Gtest_test, 3);
                        hold on; lsline; title('Flow (all)');
                    end
                    
                    predyL1O_array_pnasal_pt = [predyL1O_array_pnasal_pt, predyL1O_pnasal];
                    predyL1O_array_flow_pt = [predyL1O_array_flow_pt, predyL1O_flow];
                    
                    VEVeup_Flow = PtData_flow_pt.VE./PtData_flow_pt.Veup;
                    VEVeup_pnasal = PtData_pnasal_pt.VE./PtData_pnasal_pt.Veup;
                    
                    VEVeup_flow_array_pt = [VEVeup_flow_array_pt, VEVeup_Flow];
                    VEVeup_pnasal_array_pt = [VEVeup_pnasal_array_pt, VEVeup_pnasal];
                    
                    if 0 % VE/Veup figure
                        figure(90); clf(figure(90));
                        scatter(VEVeup_Flow, VEVeup_pnasal, 3,'filled','markerfacealpha',0.5);hold on;
                        hndl_ls = lsline;  hndl_ls.LineStyle = '-'; hndl_ls.Color = [0 0 1];
                        hndl_ref = refline(1,0); hndl_ref.LineStyle = '--'; hndl_ref.Color = [0.5 0.5 0.5];
                        
                        xlabel('Flow VEVeup');ylabel('Pnasal VEVeup'); axis square;
                        [r,~] = corr(VEVeup_Flow, VEVeup_pnasal);
                        str=['VE/Veup for Flow vs Pnasal, r = ', num2str(r)]; suptitle(str);
                        
                        VEVeup_All = [VEVeup_All; [VEVeup_Flow, VEVeup_pnasal]];
                        
                        fig = gcf;  fig.Color = [1 1 1];  fig.Units = 'inches';
                        fig.Position = [-12 1 4.5 4.5];
                        str = ['..\Figures\', 'Testing_ManRef2_VEVeupFlowVsPnasal'];
                        % saveas(fig, str, 'png'); %savefig(str);
                    end
                    
                    if 0&&length(Ftr_indx)==50 % figures
                        
                        % flow_pred to Gtest_flow (but only pnasal matched breaths)
                        figure(100+subj); clf(figure(100+subj));
                        subplot(1,3,1);
                        scatter(Gtest_flow_clean,predyL1O_clean, 2,'filled','markerfacealpha',0.3);
                        
                        xlabel('Actual VEVDrive');ylabel('Flow Predicted VEVdrive');
                        str=['r = ', num2str(rho_1)]; title(str);
                        
                        % match pnasal_pred to Gtest_flow
                        subplot(1,3,2);
                        scatter(Gtest_flow_clean,pnasal_pred_clean, 2,'filled','markerfacealpha',0.3);
                        
                        xlabel('Actual VEVDrive');ylabel('Pnasal Predicted VEVdrive');
                        str=['r = ', num2str(rho_2)]; title(str);
                        
                        % scatter of flow predy vs pnasal predy
                        subplot(1,3,3);
                        scatter(predyL1O_clean, pnasal_pred_clean, 2,'filled','markerfacealpha',0.3);
                        xlabel('Flow Predicted VEVdrive');ylabel('Pnasal Predicted VEVdrive');
                        
                        str=['r = ', num2str(rho_3)]; title(str);
                        
                        str_plot=['IndivPerf Pt', num2str(subj), ' at ', num2str(length(Ftr_indx)), ' features']; suptitle(str_plot);
                        
                        fig = gcf;
                        fig.Color = [1 1 1]; % set background colour to white
                        fig.Units = 'inches';
                        fig.Position = [-12 1 14 4.5];
                        
                        str = ['..\Figures\', str_plot ];
                        saveas(fig, str, 'png'); %savefig(str);
                        
                        if 0
                            
                            % stair plot of predictive values vs real VEVdrive
                            figure(12); clf(figure(12));
                            PtPnasalTime_clean = PtPnasalTime;
                            PtPnasalTime_clean(exclude) = [];
                            
                            predytrain_ = predytrain(Isubj);
                            stairs(PtPnasalTime_clean, Gtest_flow_clean); hold on;
                            stairs(PtPnasalTime_clean, predyL1O_clean);
                            %stairs(PtPnasalTime_clean, predytrain_(matches_clean));
                            stairs(PtPnasalTime_clean, pnasal_pred_clean);
                            %legend('label','flow (L1O)','flow (train)','pnasal');
                            legend('Actual','flow pred','pnasal pred');
                            xlabel('Time');ylabel('Actual and Predicted VEVdrive');
                            
                        end
                    end
                end
                
                %% dlm testing
                if 0 %RsqL10_array(subj,length(Ftr_indx))<0
                    keyboard;
                    
                    %alternate Rsq calculation 1
                    %weights_test_s = weights_test.^0.5;
                    SSres_w = sum(weights_test.*(Gtest_test - predyL1O_array(Isubj,length(Ftr_indx))).^2);
                    SStot_w = sum(weights_test.*(Gtest_test - mean(Gtest_test)).^2);
                    SSreg_w = sum(weights_test.*(predyL1O_array(Isubj,length(Ftr_indx))-mean(Gtest_test)).^2);
                    Rsquared_1_w = 1 - (SSres_w/SStot_w)
                    Rsquared_2_w = (SSreg_w/SStot_w)
                    SStot2_w = SSres_w+SSreg_w;
                    
                    %alternate Rsq calculation 2
                    [mdl] = fitlm(Gtest_test, predyL1O_array(Isubj,length(Ftr_indx)), 'weights', weights_test);
                    Rsq_bycorr = rho^2;
                    
                    %alternate Rsq calculation 3
                    [Rsq_test,Pvals_test,MSE_test,Ypred_test]=glmfitFaster_useknownbeta(Amatrix2_test(:,Ftr_indx),Gtest_test,weights_test,1, b); %fasterer
                    
                    figure(1); clf(figure(1));
                    scatter(Ypred_test, Gtest_test);
                    xlabel('PredY'); ylabel('Actual');
                    lsline;
                    
                    figure(2); clf(figure(2));
                    scatter(Ypred_test, predyL1O_array(Isubj,length(Ftr_indx)));
                    xlabel('PredY from fit'); ylabel('PredY from above');
                    lsline;
                end
            end
            %remove least important
            if 0&&length(Ftr_indx)>200 %speed up at higher levels
                temp=[Pvals_train';1:length(Pvals_train)]';
                temp=sortrows(temp);
                remi = temp(end-20+1:end,2);
                maxp = temp(end-20+1,1);
            else
                [maxp,remi]=max(Pvals_train);
            end
            if length(Ftr_indx)==1 || (~neverbreak && maxp<alpha && (length(Ftr_indx)-1)<=MaxNfeatures)
                break
            end
            disp(['Removing Ftr: ', num2str(Ftr_indx(remi)), ', p= ' num2str(maxp)]);
            Ftr_indx(remi)=[];
            Labels(remi)=[];
        end
        
        % grow a vector of real Gtest values that match with pnasal breaths
        % also grow a matrix of pnasal predy values, with each row
        % being a pnasal breath, and each col being the length of
        % features (as per predyL1O_array for flow breaths)
        if PtHasPnasal
            PtData_matched = [PtData_matched; PtData_pnasal_pt];
            PtData_pnasal_matched =  [PtData_pnasal_matched; PtData_pnasal_pt];
            
            PnasalFlowBB_PTindex = [PnasalFlowBB_PTindex; PnasalFlowBB_PTindex_pt];
            Gtest_matched = [Gtest_matched; Gtest_test_flow];
            weights_matched = [weights_matched; weights_test_flow];
            
            predyL1O_array_flow_pt = fliplr(predyL1O_array_flow_pt); % flip them so they are same direction of other feature arrays
            predyL1O_array_pnasal_pt = fliplr(predyL1O_array_pnasal_pt);
            predyL1O_array_pnasal_matched = [predyL1O_array_pnasal_matched;  predyL1O_array_pnasal_pt];
            predyL1O_array_flow_matched =  [predyL1O_array_flow_matched;  predyL1O_array_flow_pt];
            
            VEVeup_flow_array_pt = fliplr(VEVeup_flow_array_pt);
            VEVeup_pnasal_array_pt = fliplr(VEVeup_pnasal_array_pt);
            VEVeup_flow_array = [VEVeup_flow_array; VEVeup_flow_array_pt];
            VEVeup_pnasal_array = [VEVeup_pnasal_array; VEVeup_pnasal_array_pt];
            
            Amatrix2_flow_matched = [Amatrix2_flow_matched; Amatrix2_flow_pt];
            Amatrix2_pnasal_matched = [Amatrix2_pnasal_matched; Amatrix2_pnasal_pt];
        end
        
        if 0
            Bvals = b(2:end); % excluding intercept term
            Bvals_ = abs(Bvals./std(Amatrix2_(:,If))');
            temp = std(Amatrix2_(:,If));
            if 1
                Data=[Pvals,(1:length(Pvals))'];
                Data2 = sortrows(Data);
                If2 = Data2(:,2);
            else
                Data=[Bvals_,(1:length(Bvals_))'];
                Data2 = sortrows(Data,'descend');
                If2 = Data2(:,2);
            end
            Labels2 = Labels(If2);
            P = Pvals(If2);
            B_ = Bvals_(If2);
            
            predytrain = [colofones_ Amatrix2_(:,If)]*b;
            predytrain(predytrain>maxG)=maxG;
            predytrain(predytrain<0)=0;
            predyL1O(Isubj) = [colofones(Isubj) Amatrix2(Isubj,If)]*b;
            
            Rsq = 1-nansum(w.*(Gtest-predytrain).^2)/nansum(w.*(Gtest-nanmean(Gtest)).^2);
            RvalueTrain(subj) = Rsq^0.5;
            Err(subj) = nanmean(abs(predyL1O(Isubj)-Gtest_(Isubj)));
            ErrRms(subj) = nanmean((predyL1O(Isubj)-Gtest_(Isubj)).^2).^0.5;
            Nfeatures(subj) = length(Labels2);
        end
        toc
    end
catch me
    disp(me.getReport);
end

% display processing time
delta_t = etime(clock, t_start_L1O); % delta in seconds
D = duration(0,0,round(delta_t),'Format','hh:mm:ss'); % duration in HH:MM:SS
disp(' '); % add row space for visual clarity in command window

displaytext = ['L1O linear regresion process complete. Total time: ', char(D), ' (hh:mm:ss)'];
HowSoonIsNow = datestr(datetime('now'),'yyyymmdd');
savestr = ['_LinRegWorkspace_InclPnasal', HowSoonIsNow, '.mat'];

disp(displaytext);

str=['Saving to: ', filename(1:end-4), savestr]; disp(str);
save([filename(1:end-4), savestr]);
%save([filename(1:end-4), '_withPnasalMatched', savestr]);

diary off

keyboard

% end of L1O run



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
%  _                     _  ______      _
% | |                   | | |  _  \    | |
% | |     ___   __ _  __| | | | | |__ _| |_ __ _
% | |    / _ \ / _` |/ _` | | | | / _` | __/ _` |
% | |___| (_) | (_| | (_| | | |/ / (_| | || (_| |
% \_____/\___/ \__,_|\__,_| |___/ \__,_|\__\__,_|
%
%

close all
clear
clc
addpath('C:\Users\uqdmann\Dropbox\PUPbeta_git\');
cd('C:\Users\uqdmann\Dropbox\PUPbeta_git\FlowAndEdiPesQAO');
%datadir = 'C:\PSG_Data\FlowDrive\FeatureSpaces\';
%load('C:\PSG_Data\FlowDrive\FeatureSpaces\FlowDrive_only25Hz_Normalized_FeatureSpace_AutoRef2_Edi_Clean_LinReg_LinRegWorkspace_20180302.mat');
%load('C:\PSG_Data\FlowDrive\Analyzed\FlowDrive_only25Hz_ReNormalized_FeatureSpace_AutoRef2_Edi_Clean_LinRegWorkspace_20180427.mat');

load('C:\PSG_Data\FlowDrive\Analyzed\FlowDrive_only25Hz_ReNormalized_FeatureSpace_AutoRef2_Edi_Clean_LinRegWorkspace_InclPnasal_20180531.mat');
IncludePnasal = 0; % it was included in the analysis, but turned off for current processing
% IncludePnasal = 1;


%% Use All breaths for classification performance plots
SleepOnly = 0;
if ~IncludePnasal
    PtData_matched = [];
end
[BB, BB_] = getAllorSleepOnlyBB(SleepOnly, IncludePnasal, PtData, PtData_matched);
str = ['Flow - Using ',num2str(nnz(BB)),' of ', num2str(nnz(PtData.PT))]; disp(str);

if IncludePnasal
    str = ['Pnasal - Using ',num2str(nnz(BB_)),' of ', num2str(nnz(PtData_matched.PT))]; disp(str);
end

%% Process test data - Flow
% for all pts combined
for i=1:size(predyL1O_array,2)
    [RsqL1O_(i), RL1O_(i)] = ...
        UnivariateStats(predyL1O_array(BB,i),Gtest_All(BB), weights(BB)); % new method for R and Rsq
    %RsqL1O_(i) = 1-nansum(weights(BB).*(Gtest_All(BB)-predyL1O_array(BB,i)).^2)/nansum(weights(BB).*(Gtest_All(BB)-nanmean(Gtest_All(BB))).^2);
    ErrL1O_(i) = nanmean(weights(BB).*abs(predyL1O_array((BB),i)-Gtest_All(BB)));
    ErrL1Orms_(i) = nanmean((weights(BB).*(predyL1O_array((BB),i)-Gtest_All(BB))).^2).^0.5;
end
%RL1O = RsqL1O_.^0.5; % incorrect! SQRT of R-Squared is not R

% % for each pt
% for i=1:size(predyL1O_array,2)
%     PT_list = unique(PtData.PT);
%     for subj=1:54 %length(PT_list)
%         if ~ismember(subj, PT_list); continue; end
%         Isubj = PtData.PT==subj & BB;
%         [RsqL1O_perPT(subj,i), RL1O_perPT(subj,i)] = ...
%             UnivariateStats(predyL1O_array((Isubj),i),Gtest_All(Isubj), weights(Isubj)); % new method for R and Rsq
%
%         ErrL1O_perPT(subj,i) = nanmean(weights(Isubj).*abs(predyL1O_array((Isubj),i)-Gtest_All(Isubj)));
%         ErrL1Orms_perPT(subj,i) = nanmean((weights(Isubj).*(predyL1O_array((Isubj),i)-Gtest_All(Isubj))).^2).^0.5;
%
% %         RsqL10_perPT(subj,i) = 1-nansum(weights(Isubj).*(Gtest_All(Isubj)-predyL1O_array(Isubj,i)).^2)/nansum(weights(Isubj).*(Gtest_All(Isubj)-nanmean(Gtest_All(Isubj))).^2);
% %         % incorrect! SQRT of R-Squared is not R
% %         if (RsqL10_perPT(subj,i) < 0) % preserve sign
% %             RL10_perPT(subj,i) = -1*(abs(RsqL10_perPT(subj,i)).^0.5);
% %         else
% %             RL10_perPT(subj,i) = (RsqL10_perPT(subj,i)).^0.5;
% %         end
%     end
% end

%% sampler of R and Rsq values
% indiv pt Rsq values less critical
if 0
    RsqL1O_sampler_individual = RsqL10_perPT(:,[3, 10, 20, 30, 40, 50]);
    RsqL1O_sampler_combined = RsqL1O_([3, 10, 20, 30, 40, 50]);
    RsqL1O_sampler_individual_trim= [PT_list, RsqL1O_sampler_individual(PT_list,:)];
    RsqL1O_indivPt = array2table(RsqL1O_sampler_individual_trim, 'VariableNames', {'PT','Ftrs3','Ftrs10','Ftrs20','Ftrs30','Ftrs40','Ftrs50'});
    RsqL1O_combinedPt = array2table(RsqL1O_sampler_combined, 'VariableNames', {'Ftrs3','Ftrs10','Ftrs20','Ftrs30','Ftrs40','Ftrs50'});
    % incorrect! SQRT of R-Squared is not R
    RL1O_combinedPt = array2table(RsqL1O_sampler_combined.^0.5, 'VariableNames', {'Ftrs3','Ftrs10','Ftrs20','Ftrs30','Ftrs40','Ftrs50'});
end

%% Process test data - Flow, only breaths matching with pnasal breaths
if IncludePnasal
    % for all pts combined
    for i=1:size(predyL1O_array_flow_matched,2)
        RsqL1O_flow(i) = 1-nansum(weights_matched(BB_).*(Gtest_matched(BB_)-predyL1O_array_flow_matched((BB_),i)).^2)/nansum(weights_matched(BB_).*(Gtest_matched(BB_)-nanmean(Gtest_matched(BB_))).^2);
        ErrL1O_flow(i) = nanmean(weights_matched(BB_).*abs(predyL1O_array_flow_matched((BB_),i)-Gtest_matched(BB_)));
        ErrL1Orms_flow(i) = nanmean((weights_matched(BB_).*(predyL1O_array_flow_matched((BB_),i)-Gtest_matched(BB_))).^2).^0.5;
    end
    
    RL1O_flow = RsqL1O_flow.^0.5; % incorrect! SQRT of R-Squared is not R
    
    % for each pt
    for i=1:size(predyL1O_array_flow_matched,2)
        PT_list = unique(PnasalFlowBB_PTindex);
        for subj=1:54 %length(PT_list)
            if ~ismember(subj, PT_list); continue; end
            Isubj = PnasalFlowBB_PTindex==subj  & BB_;
            RsqL10_perPT_flow(subj,i) = 1-nansum(weights_matched(Isubj).*(Gtest_matched(Isubj)-predyL1O_array_flow_matched(Isubj,i)).^2)/nansum(weights_matched(Isubj).*(Gtest_matched(Isubj)-nanmean(Gtest_matched(Isubj))).^2);
            ErrL1O_perPT_flow(subj,i) = nanmean(weights_matched(Isubj).*abs(predyL1O_array_flow_matched((Isubj),i)-Gtest_matched(Isubj)));
            ErrL1Orms_perPT_flow(subj,i) = nanmean((weights_matched(Isubj).*(predyL1O_array_flow_matched((Isubj),i)-Gtest_matched(Isubj))).^2).^0.5;
            % incorrect! SQRT of R-Squared is not R
            if (RsqL10_perPT_flow(subj,i) < 0) % preserve sign
                RL10_perPT_flow(subj,i) = -1*(abs(RsqL10_perPT_flow(subj,i)).^0.5);
            else
                RL10_perPT_flow(subj,i) = (RsqL10_perPT_flow(subj,i)).^0.5;
            end
        end
    end
end

%% Process test data - Pnasal
if IncludePnasal
    % for all pts combined
    for i=1:size(predyL1O_array_pnasal_matched,2)
        RsqL1O_pnasal(i) = 1-nansum(weights_matched(BB_).*(Gtest_matched(BB_)-predyL1O_array_pnasal_matched((BB_),i)).^2)/nansum(weights_matched(BB_).*(Gtest_matched(BB_)-nanmean(Gtest_matched(BB_))).^2);
        ErrL1O_pnasal(i) = nanmean(weights_matched(BB_).*abs(predyL1O_array_pnasal_matched((BB_),i)-Gtest_matched(BB_)));
        ErrL1Orms_pnasal(i) = nanmean((weights_matched(BB_).*(predyL1O_array_pnasal_matched((BB_),i)-Gtest_matched(BB_))).^2).^0.5;
    end
    RL1O_pnasal = RsqL1O_pnasal.^0.5; % incorrect! SQRT of R-Squared is not R
    
    % % for each pt
    % for i=1:size(predyL1O_array_pnasal_matched,2)
    %     PT_list = unique(PnasalFlowBB_PTindex);
    %     for subj=1:54 %length(PT_list)
    %         if ~ismember(subj, PT_list); continue; end
    %         Isubj = PnasalFlowBB_PTindex==subj & BB_;
    %         RsqL10_perPT_pnasal(subj,i) = 1-nansum(weights_matched(Isubj).*(Gtest_matched(Isubj)-predyL1O_array_pnasal_matched(Isubj,i)).^2)/nansum(weights_matched(Isubj).*(Gtest_matched(Isubj)-nanmean(Gtest_matched(Isubj))).^2);
    %         ErrL1O_perPT_pnasal(subj,i) = nanmean(weights_matched(Isubj).*abs(predyL1O_array_pnasal_matched((Isubj),i)-Gtest_matched(Isubj)));
    %         ErrL1Orms_perPT_pnasal(subj,i) = nanmean((weights_matched(Isubj).*(predyL1O_array_pnasal_matched((Isubj),i)-Gtest_matched(Isubj))).^2).^0.5;
    %
    %         % incorrect! SQRT of R-Squared is not R
    %         if (RsqL10_perPT_pnasal(subj,i) < 0) % preserve sign
    %             RL10_perPT_pnasal(subj,i) = -1*(abs(RsqL10_perPT_pnasal(subj,i)).^0.5);
    %         else
    %             RL10_perPT_pnasal(subj,i) = (RsqL10_perPT_pnasal(subj,i)).^0.5;
    %         end
    %     end
    % end
    % end
end

%% Select data to use in main Figures
PlotUptoNFtrs = 50;
ftrnum=25;

DataIn = 'Flow'; % Flow, FlowP Pnasal
switch DataIn
    case 'Flow'
        R = RL1O_;
        Rsq = RsqL1O_;
        Err = ErrL1O_;
        ErrRms = ErrL1Orms_;
        
        predy = predyL1O_array(BB,ftrnum);
        Yval = Gtest_All(BB);
        
        flow_only = PtData.VE(BB);
        drive_only = PtData.DriveEdi(BB);
        
        supstr = ['Flow (all)'];
    case 'FlowP'
        R = RL1O_flow;
        Rsq = RsqL1O_flow;
        Err = ErrL1O_flow;
        ErrRms = ErrL1Orms_flow;
        predy=predyL1O_array_flow_matched(BB_,ftrnum);
        Yval=Gtest_matched(BB_);
        supstr = ['Flow (matched)'];
    case 'Pnasal'
        R = RL1O_pnasal;
        Rsq = RsqL1O_pnasal;
        Err = ErrL1O_pnasal;
        ErrRms = ErrL1Orms_pnasal;
        predy = predyL1O_array_pnasal_matched(BB_,ftrnum);
        Yval=Gtest_matched(BB_);
        supstr = ['Pnasal'];
end
savefigas = '';%'saveasTIFF'; % options are saveasPNG, saveasFIG, and saveasTIFF
closefigs = 0;

TickFntSz = 12;
LabelFntSz = 18;
FntSz = 18;

%% GS vs Pred, GS vs FlowOnly, and GS vs DriveOnly
% R values (unweighted)
[r1,p1] = corr(predy, Yval);
[r2,p2] = corr(predy, flow_only);
[r3,p3] = corr(predy, drive_only);

% R values (unweighted)
[r1_] = weightedcorrs([predy, Yval], weights); r1_ = r1_(1,2);
[r2_] = weightedcorrs([predy, flow_only], weights); r2_ = r2_(1,2);
[r3_] = weightedcorrs([predy, drive_only], weights); r3_ = r3_(1,2);

% R-Squared values
[rsq1,p1sq] = glmfitFast(predy, Yval, weights,0);
[rsq2,p2sq] = glmfitFast(predy, flow_only, weights,0);
[rsq3,p3sq] = glmfitFast(predy, drive_only, weights,0);

%%
% ______ _
% |  ___(_)
% | |_   _  __ _ _   _ _ __ ___  ___
% |  _| | |/ _` | | | | '__/ _ \/ __|
% | |   | | (_| | |_| | | |  __/\__ \
% \_|   |_|\__, |\__,_|_|  \___||___/
%           __/ |
%          |___/
%

%% show training data error - one plot (two panels)
% Comment next line for runs with corect R
% RTrain_array = RsqTrain_array.^0.5; % incorrect! SQRT of R-Squared is not R
figure(4); clf(figure(4));
fig = gcf; fig.Color = [1 1 1];
fig.Units = 'inches';
fig.Position = [-13   -0.5   12   4.5];
subplot(1,2,1);
pts=unique(PtData.PT);
for pt = 1:54
    if ismember(pt, pts)
        plot(RTrain_array(pt,:), 'k'); hold on;
        plot(RsqTrain_array(pt,:), 'k'); hold on;
        plot(ErrTrain_array(pt,:), 'k'); hold on;
        plot(ErrRmsTrain_array(pt,:), 'k'); hold on;
    end
end
ax = gca; set(ax,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz);
xlabel('Number of Features', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
ylabel('Performance (Training)', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
ylim([0 1]); xlim([0 50]); axis square

% Classification performance
subplot(1,2,2);
plot([R(1:PlotUptoNFtrs);Rsq(1:PlotUptoNFtrs);Err(1:PlotUptoNFtrs);ErrRms(1:PlotUptoNFtrs)]'); hold on;
currentYlim=ylim();
plot([ftrnum, ftrnum], [0, currentYlim(2)],'k-', 'linewidth', 2);
ax = gca; set(ax,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz);
xlabel('Number of Features', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
ylabel('Performance (Testing)', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
ylim([0 1]);axis square
legend('r', 'Rsq','MAE', 'RMSE', 'location','east');
% title('Performance Vs Number of Features (Flow Testing)');
str = ['..\Figures\', 'Figure_E5'];  % was E4

% Add labels A B to plot space
subplot(1,2,1); hold on;
text(-15, 0.99, 'A', 'FontSize', 20, 'FontWeight', 'Bold');
subplot(1,2,2); hold on;
text(-15, 0.99, 'B', 'FontSize', 20, 'FontWeight', 'Bold');

switch savefigas
    case 'saveasTIFF'; print(fig, str, '-dtiff', '-r1000');
    case 'saveasPNG'; saveas(fig, str, 'png');
    case 'saveasFIG'; savefig(str);
end
if closefigs; close(fig); end


%% Pred VEVdrive vs Actual VEVdrive, for Flow(all), Flow(matched) and Pnasal
if IncludePnasal
    % plot is technically "Actual Vs Pred"
    figure(19); clf(figure(19)); fig = gcf;
    fig.Color = [1 1 1];
    fig.Units = 'inches';
    fig.Position = [  -12.2    8.5   12    4.5];
    fig.Position = [  0.5    0.5   12    4.5];
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
    if closefigs; close(fig); end
end


%% Use betas trained in pneumotach flow, tested in pnasal
% need to set/select pnasal above (about line 900)
if 0
    figure(20); clf(figure(20)); fig = gcf;
    fig.Color = [1 1 1]; fig.Units = 'inches';
    fig.Position = [0.5   0.5   19    5];
    facealpha = 0.05; % was 0.08
    facecolor = [0.1 0.1 0.1]; %was [0.5 0.5 0.5]
    
    subplot(1,3,1); set(gca,'box','off','tickdir','out','fontname','arial narrow');
    plot([R(1:PlotUptoNFtrs);Err(1:PlotUptoNFtrs);ErrRms(1:PlotUptoNFtrs)]'); hold on;
    currentYlim=ylim();
    plot([ftrnum, ftrnum], [0, currentYlim(2)],'k-', 'linewidth', 2);
    ax = gca; set(ax,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz);
    xlabel('Number of Features', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
    ylabel('Performance (Testing)', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
    ylim([0 1]);axis square
    legend('r', 'MAE', 'RMSE', 'location','east');
    
    
    % scatter with box overlay
    subplot(1,3,2); set(gca,'box','off','tickdir','out','fontname','arial narrow');
    scatter(100*predy,100*Yval,2,'filled','markerfacealpha', facealpha, 'markerfacecolor', facecolor); hold on;
    dx=0.2; xbins=[0 0.3:dx:0.9 1.5]; lsline();
    clear medianX medianY upperIQRY lowerIQRY upperIQRY2 lowerIQRY2
    for i=1:length(xbins)-1
        Ix=predy>xbins(i)&predy<xbins(i+1);
        medianX(i)=prctile(predy(Ix),50);
        medianY(i)=prctile(Yval(Ix),50);
        upperIQRY(i)=prctile(Yval(Ix),75);
        lowerIQRY(i)=prctile(Yval(Ix),25);
        upperIQRY2(i)=prctile(Yval(Ix),90);
        lowerIQRY2(i)=prctile(Yval(Ix),10);
    end
    plotbarwitherrorsmedian(100*medianY,100*medianX,100*upperIQRY2,100*lowerIQRY2,zeros(length(medianY),3),NaN,100*0.005,100*0.01);
    plotbarwitherrorsmedian(100*medianY,100*medianX,100*upperIQRY,100*lowerIQRY,zeros(length(medianY),3),NaN,100*dx/2,100*0.01);
    xlim([0 150]);
    
    ax = gca; set(ax,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz);
    ax.XTick=[0:25:150]; xticklabels(ax, {'0', '25', '50', '75', '100', '125'});
    ax.YTick=[0:25:150]; yticklabels(ax, {'0', '25', '50', '75', '100', '125'});
    if 1
        xlabel('Flow Shape Predicted {\itflow:drive} (%)', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
        ylabel('Gold Standard {\itflow:drive} (%)', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
    else
        xlabel({'{\itflow:drive} (%)','Flow Predicted'}, 'FontSize', FntSz, 'FontName', 'Arial Narrow');
        ylabel({'{\itflow:drive} (%)','Gold Standard'}, 'FontSize', FntSz, 'FontName', 'Arial Narrow');
    end
    
    axis square
    str=['r = ', num2str(round(R(ftrnum),2))];  text(110, 15, str);
    %title(['Flow Predicted Flow:Drive Vs Gold Standard {\itflow:drive}']);
    % Add labels A B to plot space
    %text(-45, 148, 'A', 'FontSize', 20, 'FontWeight', 'Bold');
    
    
    % confusion mat
    subplot(1,3,3);
    %customcmap = GetCustomColorMap('SS');
    customcmap = GetCustomColorMap('gray');
    g = NaN*Yval;
    ghat = NaN*predy;
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
    ghat(predy>classcutoffs(1))=1;
    for i=2:Nclasses
        g(Yval<classcutoffs(i-1))=i;
        ghat(predy<classcutoffs(i-1))=i;
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
    %xlabel('Flow Predicted {\itflow:drive} (%)', 'FontSize', FntSz, 'FontName', 'Arial Narrow');
    %ylabel('Gold Standard {\itflow:drive} (%)', 'FontSize', FntSz, 'FontName', 'Arial Narrow');
    xlabel('Flow Shape Classification', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
    ylabel('Gold Standard Classification', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
    
    %title(['Predictive Performance']);
    axis square
    
    % Add labels B to plot space
    %text(-1.5, 5.4, 'B', 'FontSize', 20, 'FontWeight', 'Bold');
    
    suptitle(['Using betas from training in pneumotach, then testing in pnasal. (scatter and confusion plot are at 25 ftrs)']);
    
    str = ['..\Figures\Figure_TrainFlow_TestPnasal', ]; %
    switch savefigas
        case 'saveasTIFF'; print(fig, str, '-dtiff', '-r1000');
        case 'saveasPNG'; saveas(fig, str, 'png');
        case 'saveasFIG'; savefig(str);
    end
    if closefigs; close(fig); end
    
end

%% Figure E3
figure(83); clf(figure(83)); fig = gcf;
fig.Color = [1 1 1]; fig.Units = 'inches';
fig.Position = [1  1   16   5];

% performance vs number of features
subplot(1,3,1);
plot([R(1:PlotUptoNFtrs);Rsq(1:PlotUptoNFtrs);Rsq(1:PlotUptoNFtrs).^0.5;Err(1:PlotUptoNFtrs)]'); hold on;
currentYlim=ylim();
plot([ftrnum, ftrnum], [0, currentYlim(2)],'k-', 'linewidth', 1.5);
ax = gca; set(ax,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz);
xlabel('Number of Features', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
ylabel('Performance (Testing)', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
ylim([0 1]); xlim([0 50]); axis square
if 0 % old incorrect R
    text(40, 0.88, {'correlation', 'coefficient'}, 'fontname','arial narrow','FontSize',10, 'FontWeight', 'normal');
    text(40, 0.28, {'mean', 'absolute', 'error'}, 'fontname','arial narrow','FontSize',10, 'FontWeight', 'normal');
else
    text(30, 0.8, {'R-ong - correlation coefficient'}, 'fontname','arial narrow','FontSize',10, 'FontWeight', 'normal');
    text(30, 0.7, {'R-ite - correlation coefficient'}, 'fontname','arial narrow','FontSize',10, 'FontWeight', 'normal');
    text(30, 0.62, {'R-squared - coefficient of determination'}, 'fontname','arial narrow','FontSize',10, 'FontWeight', 'normal');
    text(30, 0.25, {'mean absolute error'}, 'fontname','arial narrow','FontSize',10, 'FontWeight', 'normal');
end

R(ftrnum)
Rsq(ftrnum)

% proportion of absolute error
subplot(1,3,2);
a_out = mean(ErrTrain_array(pts,:),1);
a_diff = diff(a_out);
err100 = a_out(100);
p_of100 = 100.*(err100 ./ a_out(1:100));
plot(100-p_of100, 'k-'); hold on;
currentYlim=ylim();
plot([ftrnum, ftrnum], [0, currentYlim(2)/2],'k-', 'linewidth', 1.5);
ax = gca; set(ax,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz);
xlabel('Number of Features', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
ylabel('Proportion of absolute error (%)', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
ylim([0 30]); xlim([0 100]); axis square

% confusion mat
subplot(1,3,3);
customcmap = GetCustomColorMap('gray'); % 'SS'
g = NaN*Yval;
ghat = NaN*predy;
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
ghat(predy>classcutoffs(1))=1;
for i=2:Nclasses
    g(Yval<classcutoffs(i-1))=i;
    ghat(predy<classcutoffs(i-1))=i;
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
axis square

str = ['..\Figures\Figure_E3_newR', ];
switch savefigas
    case 'saveasTIFF'; print(fig, str, '-dtiff', '-r1000');
    case 'saveasPNG'; saveas(fig, str, 'png');
    case 'saveasFIG'; savefig(str);
end
if closefigs; close(fig); end

%% get individual patient performance in classification
% i.e. exactly correct, and correct within one severity class
classcutoffs = [0.9 0.7 0.5 0.3];
Nclasses=length(classcutoffs)+1;
ftrs = 25;
ExactClass_perPT=[];
WithinOneClass_perPT=[];
MAE_perPT=[];
for subj=1:54 %length(PT_list) % subj=2
    if ~ismember(subj, PT_list); continue; end
    Isubj = PtData.PT==subj & BB;
    Yval = Gtest_All(Isubj);
    predy = predyL1O_array(Isubj,ftrs);
    
    %figure(40); clf(figure(40));
    %scatter(Yval, predy);
    
    g = NaN*Yval;
    ghat = NaN*predy;
    g(Yval>classcutoffs(1))=1;
    ghat(predy>classcutoffs(1))=1;
    for i=2:Nclasses
        g(Yval<classcutoffs(i-1))=i;
        ghat(predy<classcutoffs(i-1))=i;
    end
    [C,order] = confusionmat(g,ghat);
    %sumactual=sum(C')';
    sumactual=sum(C,2);
    sumestimated=sum(C);
    %C_Factual=C./sum(C')'
    C_Factual=C./sum(C,2)*100; %rows are actual, cols are estimated
    C_Festimated=C./sum(C)*100;
    C_Total = (C./sum(C(:)))*100;
    
    AccA_C_Factual = C_Factual.*diag(ones(1,length(C)));
    AccB_C_Factual = C_Factual.*(diag(ones(1,length(C)-1),1)+diag(ones(1,length(C)-1),-1));
    AccAactual = mean(sum(AccA_C_Factual,2));
    AccBactual = mean(sum(AccA_C_Factual,2)) + mean(sum(AccB_C_Factual,2));
    AccA_C_Festimated = C_Festimated.*diag(ones(1,length(C)));
    AccB_C_Festimated = C_Festimated.*(diag(ones(1,length(C)-1),1)+diag(ones(1,length(C)-1),-1));
    AccAestimated = mean(sum(AccA_C_Festimated));
    AccBestimated = mean(sum(AccA_C_Festimated)) + mean(sum(AccB_C_Festimated));
    ACCs = [AccAactual AccBactual; AccAestimated AccBestimated];
    %second col gives accuracy if accepting next-category error as correct
    
    ExactClass_perPT(subj) = AccAestimated;
    WithinOneClass_perPT(subj) = AccBestimated;
    MAE_perPT(subj) = ErrL1O_perPT(subj,i);
    
end

PerPtPerformance = table([1:41]',...
    round(ExactClass_perPT(PT_list),1)',...
    round(WithinOneClass_perPT(PT_list),1)',...
    round(MAE_perPT(PT_list),2)', ...
    'VariableNames', {'PT','Exact','WithinOne','MAE'})


%% Add apnea and low flow breaths back in for everything hereafter
if 1 % need to do this if the individual patient performance is run above
    predy = predyL1O_array(BB,ftrnum);
    Yval = Gtest_All(BB);
end

numBBinTest = NaN(54,6);
PredY = []; Gtest = []; PTarray=[];
RealFlow = 1; % one for pneumotach, zero for pnasal --> better yet, use dedicated function for Pnasal
if RealFlow
    PT_list = unique(PtData.PT);
    totalNumPts = 54;
else
    %PT_list =[3   8 9 18 22       35 36 38 43 46 47 50 53 54]; % that's only 14
    PT_list = [3 5 8 9 18 22 33 34 35 36 38 43 46 47 50 53 54]; % this is full list, N=17 with Pnasal
    totalNumPts = 54;
end
%%
for subj=1:totalNumPts
    if ~ismember(subj, PT_list)
        continue
    end
    
    if RealFlow
        Isubj=(PtData.PT==subj); % find all the breaths that belong to this pt, this is just for counting
        numBBinTest(subj,1) = nnz(Isubj); % the number of breaths for this pt
        
        if 1
            Isubj=(PtData.PT==subj)&(PtData.Hypnog<4)&(PtData.Ar==0); % find the sleep only breaths that belong to this pt
        else
            Isubj=(PtData.PT==subj)&(PtData.NotAr==1);
        end
        numBBinTest(subj,2) = nnz(Isubj); % the number of sleep breaths for this pt
        %PredY_pt = predyL1O_array(Isubj,ftrnum);
        PredY_pt = predy(Isubj);
    else
        Isubj=(PtData_matched.PT==subj); % find all the breaths that belong to this pt, this is just for counting
        numBBinTest(subj,1) = nnz(Isubj); % the number of breaths for this pt
        if 1
            Isubj=(PtData_matched.PT==subj)&(PtData_matched.Hypnog<4)&(PtData_matched.Ar==0); % find the sleep only breaths that belong to this pt
        else
            Isubj=(PtData_matched.PT==subj)&(PtData_matched.NotAr==1);
        end
        numBBinTest(subj,2) = nnz(Isubj); % the number of sleep breaths for this pt
        PredY_pt = predyL1O_array_pnasal_matched(Isubj,ftrnum);
    end
    
    Gtest_pt = Gtest_All(Isubj);
    
    % add apnea breaths as 0.1 to both Gtest_pt and PredY_pt
    indAP = find(RemovedBB_Apnoea.Pt == subj);
    indLF = find(RemovedBB_LowFlow.Pt == subj);
    numApBB = RemovedBB_Apnoea{indAP,2};
    if isempty(numApBB); numApBB=0; end
    numBBinTest(subj,3) = numApBB;
    numLFBB = RemovedBB_LowFlow{indLF,2};
    if isempty(numLFBB); numLFBB=0; end
    numBBinTest(subj,4) = numLFBB;
    numRemBB = numApBB + numLFBB;
    
    PredY_pt = [PredY_pt;ones(numRemBB,1)*0.09];
    numBBinTest(subj,5) = length(PredY_pt);
    Gtest_pt = [Gtest_pt;ones(numRemBB,1)*0.09];
    numBBinTest(subj,6) = length(Gtest_pt);
    
    
    PTarray = [PTarray; ones(length(PredY_pt),1)*subj]; % same as PtData.PT but includes apnea breaths
    PredY = [PredY; PredY_pt];
    Gtest = [Gtest; Gtest_pt];
    
end
%
% tidy up table, remove nan rows
% rows in numBBinTest are each patient
% cols in numBBinTest are:
% Total BB, Sleep BB, Ap BB, Low flow BB, Total BB to use (with AP)
numBBinTest = numBBinTest(PT_list,:);
if length(PredY) ~= length(Gtest); keyboard; end                % check same length
if nnz(isnan(PredY))>0 || nnz(isnan(Gtest))>0; keyboard; end    % check for any NaN's
if (nnz(PredY<0)>0); Predy(PredY<0)=0; end                      % force lower limit
if nnz(PredY>maxG)>0; Predy(PredY>maxG)=1.5; end                % force upper limit

%% get the AHI data
[AHI_perPT, AHI_perPT_table] = getAHI_postanalysis();
if nnz(PT_list ~= find(~isnan(AHI_perPT(:,1)))) > 0; keyboard; end % check consistency
AHI_perPT_ = AHI_perPT(~isnan(AHI_perPT(:,1)),1);

%% histogram of pred vs actual VEVdrive, and do per pt processing in loop
% get per pt median data
% get per pt <threshold data

if 1; edges = [0:0.1:1.1]; else edges = xbins; end

figure(22); clf(figure(22)); fig = gcf;
fig.Color = [1 1 1]; fig.Units = 'inches';
fig.Position = [0.5   0.5   12    5];

subplot(1,2,1);
histogram(Gtest,edges, 'facealpha',0.2, 'edgealpha', 0.2, 'normalization','probability'); hold on;
currentYlim = ylim();
plot([median(Gtest), median(Gtest)], [0, currentYlim(2)],'k-', 'linewidth', 2);
titlestr = ['Actual (median = ', num2str(round(median(Gtest),2)), ')'];
title(titlestr);
ax = gca;  set(ax,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz);
ax.XTick=[0:0.25:1.50]; xticklabels(ax, {'0', '25', '50', '75', '100', '125'});
%yticks([]);
xlabel('Gold Standard {\itflow:drive} (%)', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
ylabel('Relative probability', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
xlim([-0.05 1.15]);
box off

subplot(1,2,2);
histogram(PredY,edges, 'facealpha',0.2, 'edgealpha', 0.2, 'normalization','probability'); hold on;
currentYlim = ylim();
plot([median(PredY), median(PredY)], [0, currentYlim(2)],'k-', 'linewidth', 2);
titlestr = ['Predicted (median = ', num2str(round(median(PredY),2)), ')'];
title(titlestr);
ax = gca;  set(ax,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz);
ax.XTick=[0:0.25:1.50]; xticklabels(ax, {'0', '25', '50', '75', '100', '125'});
%yticks([]);
xlabel('Flow Shape Predicted {\itflow:drive} (%)', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
ylabel('Relative probability', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
xlim([-0.05 1.15]);
box off

supstr=['Histograms of Actual and Flow Shape Predicted {\itflow:drive} (All pts, incl Ap)'];
%suptitle(supstr);
savestr = regexprep(supstr,'[_,:{\\}]','');
str = ['..\Figures\', savestr];
switch savefigas
    case 'saveasTIFF'; print(fig, str, '-dtiff', '-r1000');
    case 'saveasPNG'; saveas(fig, str, 'png');
    case 'saveasFIG'; savefig(str);
end
if closefigs; close(fig); end

%% histograms for individuals
Gtest_avg = NaN(54,1); PredY_avg = NaN(54,1);
thres1 = 0.5; % severe
thres2 = 0.7; % moderate\severe
Gtest_thres1 = NaN(54,1); Gtest_thres2 = NaN(54,1);
PredY_thres1 = NaN(54,1); PredY_thres2 = NaN(54,1);

doHistogramsPlots = 0;

for subj=1:54 % set pt num or all, no data for PT=1,  % subj=2
    if ~ismember(subj, PT_list) %PT_list)
        continue
    end
    
    Isubj = PTarray==subj;      % get the pt BB's
    Gtest_pt = Gtest(Isubj);    % actual
    PredY_pt = PredY(Isubj);    % pred
    
    % medians
    Gtest_avg(subj) = median(Gtest_pt);
    PredY_avg(subj) = median(PredY_pt);
    
    % below threshold
    Gtest_thres1(subj) = nnz(Gtest_pt<thres1);
    Gtest_thres2(subj) = nnz(Gtest_pt<thres2);
    PredY_thres1(subj) = nnz(PredY_pt<thres1);
    PredY_thres2(subj) = nnz(PredY_pt<thres2);
    
    if doHistogramsPlots
        figure(300+subj); clf(figure(300+subj));fig = gcf;
        fig.Color = [1 1 1]; % set background colour to white
        fig.Units = 'inches';
        
        if 0 % two panel figure, actual and predicted
            fig.Position = [0.5   0.5   12    5];
            
            subplot(1,2,1);
            h1 = histogram(Gtest_pt,edges, 'facealpha',0.2, 'edgealpha', 0.2, 'normalization', 'probability');%'pdf'
            hold on;
            currentYlim = ylim();
            %plot([median(Gtest_pt), median(Gtest_pt)], [0, max(h1.Values)],'k-', 'linewidth', 2);
            plot([median(Gtest_pt), median(Gtest_pt)], [0, currentYlim(2)],'k-', 'linewidth', 2);
            titlestr = ['Actual (median = ', num2str(round(median(Gtest_pt),2)), ')'];
            title(titlestr); xlabel('{\itflow:drive} (%)'); ylabel('Relative probability'); %yticks([]);
            xlim([-0.05 1.15]); ax = gca;
            set(ax,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz);
            ax.XTick=[0:0.2:1]; xticklabels(ax, {'0', '20', '40', '60', '80', '100'});
            box off
            
            subplot(1,2,2);
            h2 = histogram(PredY_pt,edges, 'facealpha',0.2, 'edgealpha', 0.2, 'normalization','probability');% 'pdf'); %
            hold on;
            currentYlim = ylim();
            %plot([median(PredY_pt), median(PredY_pt)], [0, max(h2.Values)],'k-', 'linewidth', 2);
            plot([median(PredY_pt), median(PredY_pt)], [0, currentYlim(2)],'k-', 'linewidth', 2);
            titlestr = ['Predicted (median = ', num2str(round(median(PredY_pt),2)), ')'];
            title(titlestr); xlabel('{\itflow:drive} (%)'); ylabel('Relative probability');
            xlim([-0.05 1.15]); ax = gca;
            set(ax,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz);
            ax.XTick=[0:0.2:1]; xticklabels(ax, {'0', '20', '40', '60', '80', '100'});
            box off
            
            supstr = ['Histograms of Actual and Flow Predicted {\itflow:drive} (pt ',num2str(subj),', AHI ', num2str(round(AHI_perPT(subj,1))), ', withApBB)'];
            suptitle(supstr);
            
            % save
            savestr = regexprep(supstr,'[_,:{}]','');
            str = ['..\Figures\', savestr];
            switch savefigas
                case 'saveasTIFF'; print(fig, str, '-dtiff', '-r1000');
                case 'saveasPNG'; saveas(fig, str, 'png');
                case 'saveasFIG'; savefig(str);
            end
            if closefigs; close(fig); end
            
        else % one panel, actual only
            LabelFntSz = 24;
            fig.Position = [-19   -0.5   6    5];
            h1 = histogram(Gtest_pt,edges, 'facealpha',1, 'edgealpha', 0, 'normalization', 'probability');%'pdf'
            hold on;
            currentYlim = ylim();
            %plot([median(Gtest_pt), median(Gtest_pt)], [0, max(h1.Values)],'k-', 'linewidth', 2);
            plot([median(Gtest_pt), median(Gtest_pt)], [0, currentYlim(2)],'k-', 'linewidth', 2);
            %ax.XTick=[0:0.2:1];
            %xticklabels(ax, {'0', '20', '40', '60', '80', '100'});
            ax = gca;  set(ax,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz);
            ax.XTick=[0:0.25:1];
            xticklabels(ax, {'0','25','50','75','100'});
            xlabel('{\itflow:drive} (%)', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
            ylabel('Frequency', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow'); %yticks([]); % not using 'Relative probability'
            xlim([-0.05 1.15]);
            box off
            
            supstr = ['pt ',num2str(subj),', AHI ', num2str(round(AHI_perPT(subj,1))), ', median ' num2str(round(median(Gtest_pt*100)))];
            %title(supstr);
            
            % save
            savestr = regexprep(supstr,'[_,:{}]','');
            str = ['..\Figures\Histograms\AHIandMedian_', savestr];
            
            switch savefigas
                case 'saveasTIFF'; print(fig, str, '-dtiff', '-r1000');
                case 'saveasPNG'; saveas(fig, str, 'png');
                case 'saveasFIG'; savefig(str);
            end
            if closefigs; close(fig); end
            
        end
    end
end

% tidy up
exclude = isnan(Gtest_avg);
Gtest_avg(exclude) = []; PredY_avg(exclude) = [];
Gtest_thres1(exclude) = []; Gtest_thres2(exclude) = [];
PredY_thres1(exclude) = []; PredY_thres2(exclude) = [];
Gtest_thres1_pFL = Gtest_thres1./numBBinTest(:,5); % Proportion sev FL
Gtest_thres2_pFL = Gtest_thres2./numBBinTest(:,5); % Proportion modsev FL
PredY_thres1_pFL = PredY_thres1./numBBinTest(:,5);
PredY_thres2_pFL = PredY_thres2./numBBinTest(:,5);

%% Figure 4
figure(21); clf(figure(21)); fig = gcf;
fig.Color = [1 1 1]; fig.Units = 'inches';
fig.Position = [-19   5  12   4.5];
facealpha = 0.05; % was 0.08
facecolor = [0.1 0.1 0.1]; %was [0.5 0.5 0.5]
% scatter with box overlay
subplot(1,2,1); set(gca,'box','off','tickdir','out','fontname','arial narrow');
scatter(100*predy,100*Yval,2,'filled','markerfacealpha', facealpha, 'markerfacecolor', facecolor); hold on;
dx=0.2; xbins=[0 0.3:dx:0.9 1.5]; lsline();
clear medianX medianY upperIQRY lowerIQRY upperIQRY2 lowerIQRY2
for i=1:length(xbins)-1
    Ix=predy>xbins(i)&predy<xbins(i+1);
    medianX(i)=prctile(predy(Ix),50);
    medianY(i)=prctile(Yval(Ix),50);
    upperIQRY(i)=prctile(Yval(Ix),75);
    lowerIQRY(i)=prctile(Yval(Ix),25);
    upperIQRY2(i)=prctile(Yval(Ix),90);
    lowerIQRY2(i)=prctile(Yval(Ix),10);
end
plotbarwitherrorsmedian(100*medianY,100*medianX,100*upperIQRY2,100*lowerIQRY2,zeros(length(medianY),3),NaN,100*0.005,100*0.01);
plotbarwitherrorsmedian(100*medianY,100*medianX,100*upperIQRY,100*lowerIQRY,zeros(length(medianY),3),NaN,100*dx/2,100*0.01);
xlim([0 150]);

ax = gca; set(ax,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz);
ax.XTick=[0:25:150]; xticklabels(ax, {'0', '25', '50', '75', '100', '125'});
ax.YTick=[0:25:150]; yticklabels(ax, {'0', '25', '50', '75', '100', '125'});
if 1
    xlabel('Flow Shape Predicted {\itflow:drive} (%)', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
    ylabel('Gold Standard {\itflow:drive} (%)', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
else
    xlabel({'{\itflow:drive} (%)','Flow Predicted'}, 'FontSize', FntSz, 'FontName', 'Arial Narrow');
    ylabel({'{\itflow:drive} (%)','Gold Standard'}, 'FontSize', FntSz, 'FontName', 'Arial Narrow');
end

axis square
str=['r = ', num2str(round(R(ftrnum),2))];  text(110, 15, str);
%title(['Flow Predicted Flow:Drive Vs Gold Standard {\itflow:drive}']);
% Add labels A B to plot space
text(-45, 148, 'A', 'FontSize', 20, 'FontWeight', 'Bold');

subplot(1,2,2);
[r_1, p_1] = plotregressionwithSEM(PredY_avg.*100, Gtest_avg.*100);
xlim([-5 110]); ylim([-5 110]);
ax = gca; set(ax,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz);
ax.XTick=[0:25:125];
ax.YTick=[0:25:125];
xlabel({'Flow Shape Predicted {\itflow:drive}','median (%)'}, 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
ylabel({'Gold Standard {\itflow:drive}','median (%)'}, 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
axis square
str=['r = ', num2str(round(r_1,2))]; text(80, 7, str);
text(-45, 108, 'B', 'FontSize', 20, 'FontWeight', 'Bold');

str = ['..\Figures\Figure_4', ]; %
switch savefigas
    case 'saveasTIFF'; print(fig, str, '-dtiff', '-r1000');
    case 'saveasPNG'; saveas(fig, str, 'png');
    case 'saveasFIG'; savefig(str);
end
if closefigs; close(fig); end

% determine error
if 0
    [Rvalue,Pvalue,Slope,Intercept,Xmodel_sem]=plotregressionwithSEM(predy.*100, Yval.*100);
    linefit=Slope*(predy.*100)+Intercept;
    absdiff = abs((predy.*100)-linefit);
    stdabsdiff = std(absdiff);
    meanabsdiff = mean(absdiff);
    
    [rho, pval ] = corr(predy, Yval);
    
    % weighted mean absolute error
    100*mean(ErrL1O_perPT(PT_list,ftrnum))
    
    % per breath
    100*mean(abs(predy-Yval))
    
    % per patient
    100*mean(abs(PredY_avg-Gtest_avg))
    
end
%% Novel metrics, Median VE:Vdrive during sleep and Time with severe obstruction during sleep
% aka, Proportion of breaths FL
figure(27); clf(figure(27)); fig = gcf;
fig.Color = [1 1 1]; fig.Units = 'inches';
fig.Position = [-14   1   12    4.5];
%fig.Position = [0.5   0.5   12    4.5];
subplot(1,2,1);

[r_1, p_1] = plotregressionwithSEM(PredY_avg.*100, Gtest_avg.*100);
xlim([-5 110]); ylim([-5 110]);
ax = gca; ax.FontSize = 12; LabelFntSz = 18;
ax.XTick=[0:25:125];
ax.YTick=[0:25:125];
if 0 % olf labels
    xlabel('Flow Predicted {\itflow:drive}');
    ylabel('Gold Standard {\itflow:drive}');
else
    xlabel('Flow Shape Predicted {\itflow:drive}, median (%)', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
    ylabel('Gold Standard {\itflow:drive}, median (%)', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
end

axis square
str=['r = ', num2str(round(r_1,2))]; text(80, 7, str);
% titlestr = ['Patient Median {\itflow:drive}']; title(titlestr);

subplot(1, 2, 2);
[r_3, p_3] = plotregressionwithSEM(PredY_thres2_pFL.*100, Gtest_thres2_pFL.*100);
xlim([-5 110]); ylim([-5 110]);
ax = gca; ax.FontSize = 12; LabelFntSz = 18;
ax.XTick=[0:25:125];
ax.YTick=[0:25:125];
if 0 % old labels
    xlabel('Flow Predicted %FL');
    ylabel('Gold Standard %FL');
else
    xlabel('Flow Shape Predicted % Flow Limited', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
    ylabel('Gold Standard % Flow Limited', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
end

axis square
str=['r = ', num2str(round(r_3,2))]; text(80, 7, str);
% titlestr=['Proportion with ModSev FL']; title(titlestr);

% subplot(1, 3, 3);
% [r_2, p_2] = plotregressionwithSEM(PredY_thres1_pFL.*100, Gtest_thres1_pFL.*100);
% xlim([-5 125]); xlabel('Flow Predicted %FL');
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

str = ['..\Figures\', 'Figure_old_4'];
switch savefigas
    case 'saveasTIFF'; print(fig, str, '-dtiff', '-r1000');
    case 'saveasPNG'; saveas(fig, str, 'png');
    case 'saveasFIG'; savefig(str);
end
if closefigs; close(fig); end

%% AHI vs Novel Metrics
figure(32);clf(figure(32)); fig = gcf;
fig.Color = [1 1 1]; fig.Units = 'inches';
% fig.Position = [ -12.2 8  12    4.5];
fig.Position = [ 0.5 0.5  10    10];

% gold standard
subplot(2,2,1)
[r_1, p_1]=plotregressionwithSEM(AHI_perPT_(:,1), Gtest_avg.*100);
ylim([-5 125]); xlim([-5 105]);
ax = gca; ax.FontSize = 12; LabelFntSz = 18;
ax.YTick=[0:25:125];
ax.XTick=[0:20:100];
ylabel('Gold Standard {\itflow:drive}, median (%)', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
xlabel('AHI', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
axis square
%str=['r = ', num2str(round(r_1,2))];  text(5, 15, str);
% title(['AHI vs Median {\itflow:drive}']);

subplot(2,2,2)
[r_3, p_3]=plotregressionwithSEM(AHI_perPT_(:,1), Gtest_thres2_pFL.*100);
ylim([-5 125]); xlim([-5 105]);
ax = gca; ax.FontSize = 12; LabelFntSz = 18;
ax.YTick=[0:25:125];
ax.XTick=[0:20:100];
ylabel('Gold Standard % Flow Limited', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
xlabel('AHI', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
axis square
%str=['r = ', num2str(round(r_3,2))];  text(70, 7, str);
% title(['AHI Vs Proportion with ModSev FL']);

% predicted
subplot(2,2,3)
[r_1, p_1]=plotregressionwithSEM(AHI_perPT_(:,1), PredY_avg.*100);
ylim([-5 125]); xlim([-5 105]);
ax = gca; ax.FontSize = 12; LabelFntSz = 18;
ax.YTick=[0:25:125];
ax.XTick=[0:20:100];
ylabel('Flow Shape Predicted {\itflow:drive}, median (%)', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
xlabel('AHI', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
axis square
%str=['r = ', num2str(round(r_1,2))];  text(5, 15, str);
% title(['AHI vs Median {\itflow:drive}']);

subplot(2,2,4)
[r_3, p_3]=plotregressionwithSEM(AHI_perPT_(:,1), PredY_thres2_pFL.*100);
ylim([-5 125]); xlim([-5 105]);
ax = gca; ax.FontSize = 12; LabelFntSz = 18;
ax.YTick=[0:25:125];
ax.XTick=[0:20:100];
ylabel('Flow Shape Predicted % Flow Limited', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
xlabel('AHI', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
axis square
%str=['r = ', num2str(round(r_3,2))];  text(70, 7, str);

% subplot(1,3,3)
% [r_2, p_2]=plotregressionwithSEM(AHI_perPT_(:,1), PredY_thres1_pFL.*100);
% xlim([-5 105]); ylabel('Gold Standard FL %');
% ylim([-5 105]); xlabel('AHI'); axis square
% str=['r = ', num2str(round(r_2,2))];  text(70, 7, str);
% title(['AHI Vs Proportion with Sev FL']);

% Add labels A B C D to plot space
subplot(2,2,1); hold on;
text(-45, 124, 'A', 'FontSize', 20, 'FontWeight', 'Bold');
subplot(2,2,2); hold on;
text(-45, 124, 'B', 'FontSize', 20, 'FontWeight', 'Bold');
subplot(2,2,3); hold on;
text(-45, 124, 'C', 'FontSize', 20, 'FontWeight', 'Bold');
subplot(2,2,4); hold on;
text(-45, 124, 'D', 'FontSize', 20, 'FontWeight', 'Bold');

str = ['..\Figures\', 'Figure_E9']; % Figure sup
switch savefigas
    case 'saveasTIFF'; print(fig, str, '-dtiff', '-r1000');
    case 'saveasPNG'; saveas(fig, str, 'png');
    case 'saveasFIG'; savefig(str);
end
if closefigs; close(fig); end


%% Figure for ASA 2018 abstract
figure(32);clf(figure(32)); fig = gcf;
fig.Color = [1 1 1]; fig.Units = 'inches';
fig.Position = [ -5 5  4.5  4.5];

% predicted
[Rvalue,Pvalue,Slope,Intercept,Xmodel_sem]=plotregressionwithSEM(AHI_perPT_(:,1), PredY_avg.*100);
ylim([-5 125]); xlim([-5 105]);
ax = gca; set(ax,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz);
ax.YTick=[0:25:125];
ax.XTick=[0:20:100];
ylabel('Estimated {\itflow:drive}, median (%)', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
xlabel('Apnea-Hypopnea Index (events/hr)', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
axis square

%determine residual SD
linefit=Slope*AHI_perPT_(:,1)+Intercept;
absdiff = abs(AHI_perPT_(:,1)-linefit);
stdabsdiff = std(absdiff);
meanabsdiff = mean(absdiff);

str = ['..\Figures\', 'ASA2018_Figure_1_scatter']; %
switch savefigas
    case 'saveasTIFF'; print(fig, str, '-dtiff', '-r1000');
    case 'saveasPNG'; saveas(fig, str, 'png');
    case 'saveasFIG'; savefig(str);
end
if closefigs; close(fig); end


%% Figure for ATS 2018 talk
figure(32);clf(figure(32)); fig = gcf;
fig.Color = [1 1 1]; fig.Units = 'inches';
fig.Position = [ -15 5  16  4.5];

% panel one - breath level
subplot(1,3,1);
facealpha = 0.05; % was 0.08
facecolor = [0.1 0.1 0.1]; %was [0.5 0.5 0.5]
% scatter with box overlay
scatter(100*predy,100*Yval,2,'filled','markerfacealpha', facealpha, 'markerfacecolor', facecolor); hold on;
dx=0.2; xbins=[0 0.3:dx:0.9 1.5]; lsline();
clear medianX medianY upperIQRY lowerIQRY upperIQRY2 lowerIQRY2
for i=1:length(xbins)-1
    Ix=predy>xbins(i)&predy<xbins(i+1);
    medianX(i)=prctile(predy(Ix),50);
    medianY(i)=prctile(Yval(Ix),50);
    upperIQRY(i)=prctile(Yval(Ix),75);
    lowerIQRY(i)=prctile(Yval(Ix),25);
    upperIQRY2(i)=prctile(Yval(Ix),90);
    lowerIQRY2(i)=prctile(Yval(Ix),10);
end
plotbarwitherrorsmedian(100*medianY,100*medianX,100*upperIQRY2,100*lowerIQRY2,zeros(length(medianY),3),NaN,100*0.005,100*0.01);
plotbarwitherrorsmedian(100*medianY,100*medianX,100*upperIQRY,100*lowerIQRY,zeros(length(medianY),3),NaN,100*dx/2,100*0.01);

xlim([-5 150]); ylim([-5 150]);
ax = gca; set(ax,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz);
ax.XTick=[0:25:150]; xticklabels(ax, {'0', '25', '50', '75', '100', '125'});
ax.YTick=[0:25:150]; yticklabels(ax, {'0', '25', '50', '75', '100', '125'});
xlabel('Estimated {\itflow:drive} (%)', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
ylabel('Measured {\itflow:drive} (%)', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');


axis square
str=['r = ', num2str(round(R(ftrnum),2))];  text(110, 15, str);
%title(['Flow Predicted Flow:Drive Vs Gold Standard {\itflow:drive}']);
% Add labels A B to plot space
%text(-45, 148, 'A', 'FontSize', 20, 'FontWeight', 'Bold');

% panel two - patient medians
subplot(1,3,2);
[r_1, p_1] = plotregressionwithSEM(PredY_avg.*100, Gtest_avg.*100);

%xlim([-5 110]); ylim([-5 110]);
xlim([-5 150]); ylim([-5 150]);
ax = gca; set(ax,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz);
ax.XTick=[0:25:150]; xticklabels(ax, {'0', '25', '50', '75', '100', '125'});
ax.YTick=[0:25:150]; yticklabels(ax, {'0', '25', '50', '75', '100', '125'});

%xlabel({'Estimated {\itflow:drive}','median (%)'}, 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
%ylabel({'Measured {\itflow:drive}','median (%)'}, 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
xlabel('Estimated {\itflow:drive},median (%)', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
ylabel('Measured {\itflow:drive}, median (%)', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
axis square
str=['r = ', num2str(round(r_1,2))]; text(80, 7, str);
%text(-45, 108, 'B', 'FontSize', 20, 'FontWeight', 'Bold');

% error in patient medians
% mean(abs(PredY_avg-Gtest_avg))*100

% panel three - AHI scatter
subplot(1,3,3)
[r_1, p_1]=plotregressionwithSEM(AHI_perPT_(:,1), Gtest_avg.*100);
%ylim([-5 125]); xlim([-5 105]);
xlim([-5 100]); ylim([-5 150]);
ax = gca; set(ax,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz);
ax.XTick=[0:20:100]; xticklabels(ax, {'0', '20', '40', '60', '80', '100'});
ax.YTick=[0:25:150]; yticklabels(ax, {'0', '25', '50', '75', '100', '125'});
ylabel('Measured {\itflow:drive}, median (%)', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
xlabel('AHI', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
axis square

str = ['..\Figures\', 'ATS2018_Figure_BreathMedianAHI_scatter']; %
switch savefigas
    case 'saveasTIFF'; print(fig, str, '-dtiff', '-r1000');
    case 'saveasPNG'; saveas(fig, str, 'png');
    case 'saveasFIG'; savefig(str);
end
if closefigs; close(fig); end

%% NewFigure 2  and calculate residual SD
figure(32);clf(figure(32)); fig = gcf;
fig.Color = [1 1 1]; fig.Units = 'inches';
fig.Position = [ -5 5  4.5  4.5];

% gold standard
[Rvalue,Pvalue,Slope,Intercept,Xmodel_sem]=plotregressionwithSEM(AHI_perPT_(:,1), Gtest_avg.*100);
ylim([-5 125]); xlim([-5 105]);
ax = gca; set(ax,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz);
ax.YTick=[0:25:125];
ax.XTick=[0:20:100];
ylabel({'Gold Standard {\itflow:drive}','median (%)'}, 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
xlabel('Apnea-Hypopnea Index (evts/hour)', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
axis square

%determine residual SD
linefit=Slope*AHI_perPT_(:,1)+Intercept;
absdiff = abs(AHI_perPT_(:,1)-linefit);
stdabsdiff = std(absdiff);
meanabsdiff = mean(absdiff);

str = ['..\Figures\', 'NewFigure_2_scatter']; % Figure sup
switch savefigas
    case 'saveasTIFF'; print(fig, str, '-dtiff', '-r1000');
    case 'saveasPNG'; saveas(fig, str, 'png');
    case 'saveasFIG'; savefig(str);
end
if closefigs; close(fig); end

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


%%
%  _____            ______         _
% |_   _|           |  ___|       | |
%   | | ___  _ __   | |_ ___  __ _| |_ _   _ _ __ ___  ___
%   | |/ _ \| '_ \  |  _/ _ \/ _` | __| | | | '__/ _ \/ __|
%   | | (_) | |_) | | ||  __/ (_| | |_| |_| | | |  __/\__ \
%   \_/\___/| .__/  \_| \___|\__,_|\__|\__,_|_|  \___||___/
%           | |
%           |_|


%% Histogram of ftrs selected during training
ftr_array = NaN(54, ftrnum);
for subj=1:54 % set pt num or all, no data for PT=1
    if ~ismember(subj, PT_list) %PT_list)
        continue
    end
    ftr_array(subj,:) = labels_Step_Subj{subj,ftrnum};
end
ftr_array = ftr_array(PT_list,:);
uniqueftrs = unique(ftr_array);
ftr_array_linear = ftr_array(:);
ctrs = 1:1:150;  % ToDo: make sure max ctrs covers all features
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
str = ['..\Figures\', 'HistogramOfLinRegPneumotachFeatures_SimpleCount']; %
switch savefigas
    case 'saveasTIFF'; print(fig, str, '-dtiff', '-r1000');
    case 'saveasPNG'; saveas(fig, str, 'png');
    case 'saveasFIG'; savefig(str);
end
if closefigs; close(fig); end

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
try
    for i=1:size(labels_Step_Subj,1)
        score1 = maxscore + 0*score;
        If = [];
        for j=1:50%size(labels_Step_Subj,2)
            temp = labels_Step_Subj{i,j};
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
    %LabelsOrdered = Labels_Complete(scoredata(:,2));
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
str = ['..\Figures\', 'HistogramOfLinRegPneumotachFeatures_WeightedCount']; %
switch savefigas
    case 'saveasTIFF'; print(fig, str, '-dtiff', '-r1000');
    case 'saveasPNG'; saveas(fig, str, 'png');
    case 'saveasFIG'; savefig(str);
end
if closefigs; close(fig); end


%% compare the two histogram ways of "top 50" ftrs
x_over = ismember(lbls_w(1:ftrnum), lbls_uw(1:ftrnum));
nnz(x_over) % this shows the number of features that occur in both lists


%% Alternatively, find best features through reverse stepwise regression
% Unique ftrs selected during training at set stop point
setftrnum = 25;
ftr_array_2 = NaN(54, setftrnum);
for subj=1:54 % set pt num or all, no data for PT=1
    if ~ismember(subj, PT_list) %PT_list)
        continue
    end
    ftr_array_2(subj,:) = labels_Step_Subj{subj,setftrnum};
end
ftr_array_2 = ftr_array_2(PT_list,:);
uniqueftrs = unique(ftr_array_2);

% do backwards elimination
RemoveList = [];
while 1
    [~,Pvals,~,~,~]=glmfitFast(Amatrix2(:,uniqueftrs),Gtest_All,weights,1);
    Pvals(1)=[];
    
    %remove least important
    [maxp,remi]=max(Pvals);
    
    disp(['Removing Ftr: ', num2str(uniqueftrs(remi)), ', p= ' num2str(maxp)]);
    RemoveList = [RemoveList; uniqueftrs(remi)];
    
    uniqueftrs(remi)=[];
    
    if length(uniqueftrs)<1
        break
    end
end

% Remove list is currently in the order of first removed to last removed,
% so we really want those last removed, hence flip
TopElimFtrs = fliplr(RemoveList');
Top25Ftrs = TopElimFtrs(1:25);

%% run model using only top 25 features

if 1 % either of the histogram methods
    NfeaturesOpt=25;
    if 0 % unweighted histogram
        I = I_uw;
    else % weighted histogram
        I = I_w';
    end
    LabelsOrdered = Labels_Complete(I(1:NfeaturesOpt));
    TopFtrs = I(1:NfeaturesOpt);
else % else using backward elimination of unique ftrs at threshold
    LabelsOrdered = Labels_Complete((Top25Ftrs));
    TopFtrs = Top25Ftrs;
end


%% Final model using all data, uses selected N optimal features (NfeaturesOpt)
Ftrs = TopFtrs;
Labels = [{'Intercept'};LabelsOrdered];
FtrVals = Amatrix2(:,Ftrs);

tic
[Rsq,Pvals,RMSE,betas]=glmfitFast(FtrVals,Gtest_All,weights,1); %faster
mdl_pred_1 = [ones(length(weights),1) FtrVals]*betas;
mdl_pred_1(mdl_pred_1>maxG)=maxG; % set upper limit
mdl_pred_1(mdl_pred_1<0) = 0; % set lower limit
%Pvals(1)=[]; % remove intercept term
%betas(1)=[]; % remove intercept term
toc

tic
finalmdl = fitglm(FtrVals,Gtest_All,'weights',weights);
mdl_pred_2 = predict(finalmdl, FtrVals);
mdl_pred_2(mdl_pred_2>maxG)=maxG; % set upper limit
mdl_pred_2(mdl_pred_2<0) = 0; % set lower limit
toc

% test mdl_prediction against Gtest_
if 0
    figure(4); clf(figure(4)); fig = gcf;
    fig.Color = [1 1 1]; fig.Units = 'inches';
    fig.Position = [0.5   5   12    4.5];
    subplot(1,2,1);
    scatter(Gtest_matched, mdl_pred_1,2,'filled','markerfacealpha',0.4);
    xlabel('Gtest'); ylabel('Mdl Pred 1');
    subplot(1,2,2);
    scatter(Gtest_matched, mdl_pred_2,2,'filled','markerfacealpha',0.4);
    xlabel('Gtest'); ylabel('Mdl Pred 2');
end

%finalmdl.Coefficients.tStat(2:end),
finalmdl_summary = [[0;Ftrs'], betas, finalmdl.Coefficients.SE(1:end), Pvals];
finalmdl_summary_wLabels = table([0;Ftrs'], Labels, betas, finalmdl.Coefficients.SE(1:end), Pvals, ...
    'VariableNames', {'FtrNum', 'FtrName', 'Betas', 'SE', 'Pval'});
finalmdl_summary_wLabels_noP = table(Labels, betas, finalmdl.Coefficients.SE(1:end), ...
    'VariableNames', {'FtrName', 'Betas', 'SE'});


%% Find univariate performance of 'best' features (... Bivariate ...)
rVsFD = NaN(length(Ftrs),1);
rVsPnasal = NaN(length(Ftrs),1);

r2VsFD = NaN(length(Ftrs),1);
r2VsPnasal = NaN(length(Ftrs),1);

bias =  NaN(length(Ftrs),1);

facealpha = 0.05; % was 0.08
facecolor = [0.2 0.2 0.2]; %was [0.1 0.1 0.1]

if ~(exist('TnTPnasal','var')==1) % only need to do once
    TnTPnasal = load('C:\PSG_Data\FlowDrive\Analyzed\FlowDrive_only25Hz_ReNormalized_FeatureSpace_AutoRef2_Edi_Clean_LinRegWorkspace_TestNTrainPnasal_20180530.mat',...
        'Amatrix2_flow_matched', 'Amatrix2_pnasal_matched', 'Labels_Complete');
end

for ft=1:length(Ftrs)
    % R vs flow:drive
    FtrVal = Amatrix2(:,Ftrs(ft));
    [Rsq,~,~,~,~]=glmfitFast(FtrVal,Gtest_All,weights,1);
    r2VsFD(ft) = Rsq;
    R_ = weightedcorrs([FtrVal,Gtest_All],weights);
    rVsFD(ft) = R_(1,2); %[r2VsFD(ft).^0.5, rVsFD(ft),corr(FtrVal,Gtest_All)]
    
    % R vs equivalent feature measured with Pnasal
    FtrVal_matched = TnTPnasal.Amatrix2_flow_matched(:,Ftrs(ft));
    PnasalFtrVal = TnTPnasal.Amatrix2_pnasal_matched(:,Ftrs(ft));
    [Rsq,~,~,~,~]=glmfitFast(FtrVal_matched,PnasalFtrVal,weights_matched,1);
    r2VsPnasal(ft) = Rsq;
    R_ = weightedcorrs([FtrVal_matched,PnasalFtrVal],weights_matched);
    rVsPnasal(ft) = R_(1,2);
    
    % also calculate bias, as (median value Pnasal) / (median value Pneumotach)
    bias(ft) = median(PnasalFtrVal) / median(FtrVal_matched);
    
    % figure
    if 0
        figure(101); clf(figure(101)); fig = gcf;
        fig.Color = [1 1 1]; fig.Units = 'inches';
        fig.Position = [-19   5  12   4.5];
        
        subplot(1,2,1);
        set(gca,'box','off','tickdir','out','fontname','arial narrow');
        scatter(FtrVal,Gtest_All,2,'filled','markerfacealpha', facealpha, 'markerfacecolor', facecolor); hold on;
        xlabel('Feature value (Pneumotach)'); ylabel('flow:drive');
        str=['R=',num2str(rVsFD(ft)), ',   R^2=',num2str(r2VsFD(ft))]
        title(str); axis square
        
        subplot(1,2,2);
        scatter(FtrVal_matched,PnasalFtrVal,2,'filled','markerfacealpha', facealpha, 'markerfacecolor', facecolor); hold on;
        xlabel('Feature value (Pneumotach)'); ylabel('Feature value (Pnasal)');
        str=['R=',num2str(rVsPnasal(ft)), ',   R^2=',num2str(r2VsPnasal(ft))]
        title(str); axis square
        
        title_str = ['Ftr ', num2str(ft), ' -- ', LabelsOrdered{ft}];
        title_str = regexprep(title_str,'[_,:{\\}]',' ');
        suptitle(title_str);
        save_str = ['..\Figures\Univariate\', 'Pneumo_', title_str]; %
        saveas(fig, save_str, 'png');
    end
end

Univar_summary = table(rVsFD, rVsPnasal, bias, ...
    'VariableNames', {'RvsFlowDrive', 'RvsPnasal', 'Bias'})

%% Find Top X unique features based on scores
if 0
    % if using transformed data, need to strip off the suffix
    if TransformTheData; lengthsuf=4; else; lengthsuf=0; end
    MaxUniquefeatures=50;
    LabelsOrdered_NoSuffix=LabelsOrdered;
    for i=1:length(LabelsOrdered_NoSuffix)
        LabelsOrdered_NoSuffix{i}=LabelsOrdered_NoSuffix{i}(1:length(LabelsOrdered_NoSuffix{i})-lengthsuf);
    end
    LabelsOrdered_NoSuffixUnique=unique(LabelsOrdered_NoSuffix,'stable');
    LabelsOrdered_NoSuffixUnique(MaxUniquefeatures+1:end)=[];
end

%% Another option, unique ftrs selected during training at set stop point
setftrnum = 5; % 5, this was the first peak in performance
ftr_array_2 = NaN(54, setftrnum);
for subj=1:54 % set pt num or all, no data for PT=1
    if ~ismember(subj, PT_list) %PT_list)
        continue
    end
    ftr_array_2(subj,:) = labels_Step_Subj{subj,setftrnum};
end
ftr_array_2 = ftr_array_2(PT_list,:);
uniqueftrs = unique(ftr_array_2);

FtrValsUnique = Amatrix2(:,uniqueftrs);
tic
[Rsq,Pvals,RMSE,betas]=glmfitFast(FtrValsUnique,Gtest_All,weights,1); %faster
mdl_pred_1 = [ones(length(weights),1) FtrValsUnique]*betas;
mdl_pred_1(mdl_pred_1>maxG)=maxG; % set upper limit
mdl_pred_1(mdl_pred_1<0) = 0; % set lower limit
Pvals(1)=[]; % remove intercept term
betas(1)=[]; % remove intercept term
toc
Labels_unique = Labels_Complete(uniqueftrs);
figure(41); clf(figure(4)); fig = gcf;
fig.Color = [1 1 1]; fig.Units = 'inches';
fig.Position = [-13   -0.5   12    4.5];
subplot(1,2,1);
scatter(Gtest_All, mdl_pred_1,2,'filled','markerfacealpha',0.4);
xlabel('Gtest'); ylabel('Mdl Pred 1');

subplot(1,2,2);
scatter(predy, mdl_pred_1,2,'filled','markerfacealpha',0.4);
xlabel('top 25 ftrs at step 25'); ylabel('unique ftrs (12) at at step 5');
r = Rsq.^0.5;   % incorrect! SQRT of R-Squared is not R

100*mean(abs(predy - mdl_pred_1))

%% *** univariate analysis ***

savefigas = 'saveasPNG';
Ftr=2;
figure(60); clf(figure(60)); fig = gcf;
fig.Color = [1 1 1]; fig.Units = 'inches';
fig.Position = [1   1  7.5   3];
facealpha = 0.05; % was 0.08
facecolor = [0.1 0.1 0.1]; %was [0.5 0.5 0.5]
subplot(1,2,1); set(gca,'box','off','tickdir','out','fontname','arial narrow');
scatter(Amatrix2(:,Ftr+50), Gtest_All,1,'filled','markerfacealpha', facealpha, 'markerfacecolor', facecolor); hold on;
xlabel('Feature (all breaths)'); ylabel('flow:drive'); axis square
hold on
subplot(1,2,2); set(gca,'box','off','tickdir','out','fontname','arial narrow');
scatter(Amatrix2((BB),Ftr+50), Gtest_All(BB),1,'filled','markerfacealpha', facealpha, 'markerfacecolor', facecolor); hold on;
xlabel('Feature (sleep breaths)'); ylabel('flow:drive'); axis square

FtrName = regexprep(FeatureNames{Ftr,2},'[_,:{}]',' ');
suptitle(FtrName);
str = ['..\Figures\Univariate\AllBreathsVsSleepOnly_', num2str(Ftr), '_', char(FtrName)]; %

switch savefigas
    case 'saveasTIFF'; print(fig, str, '-dtiff', '-r1000');
    case 'saveasPNG'; saveas(fig, str, 'png');
    case 'saveasFIG'; savefig(str);
end
if closefigs; close(fig); end




%% ToDo: continue from here
% get R and bias for each ftr (pneumo vs pnasal)
%

if ~(exist('TnTPnasal','var')==1) % only need to do once
    TnTPnasal = load('C:\PSG_Data\FlowDrive\Analyzed\FlowDrive_only25Hz_ReNormalized_FeatureSpace_AutoRef2_Edi_Clean_LinRegWorkspace_TestNTrainPnasal_20180530.mat',...
        'Amatrix2_flow_matched', 'Amatrix2_pnasal_matched', 'Labels_Complete', 'Gtest_matched');
end

% shown as nine panel plots;
% cols are 1. untransformed, 2. sqrt, 3. squared
% rows are 1. Flow feature, 2. pnasal feature, 3, flow vs pnasal

savefigas = 'saveasPNG';
for Ftr=1:50
    figure(60); clf(figure(60)); fig = gcf;
    fig.Color = [1 1 1]; fig.Units = 'inches';
    fig.Position = [0.1  0.1  10   10];
    facealpha = 0.05; % was 0.08
    facecolor = [0.1 0.1 0.1]; %was [0.5 0.5 0.5]
    
    subplot(3,3,1); set(gca,'box','off','tickdir','out','fontname','arial narrow');
    scatter(Amatrix2(:,Ftr), Gtest_All,1,'filled','markerfacealpha', facealpha, 'markerfacecolor', facecolor); hold on;
    ylabel('flow:drive'); xlabel('Feature value (Pneumotach)');%xlabel(['Pneumo \^1']);
    title('Untransformed');
    
    subplot(3,3,2); set(gca,'box','off','tickdir','out','fontname','arial narrow');
    scatter(Amatrix2(:,Ftr+50), Gtest_All,1,'filled','markerfacealpha', facealpha, 'markerfacecolor', facecolor); hold on;
    xlabel('Feature value (Pneumotach)');
    %ylabel('flow:drive'); %xlabel(['Pneumo \^0.5']);
    title('Sqrt');
    
    subplot(3,3,3); set(gca,'box','off','tickdir','out','fontname','arial narrow');
    scatter(Amatrix2(:,Ftr+100), Gtest_All,1,'filled','markerfacealpha', facealpha, 'markerfacecolor', facecolor); hold on;
    xlabel('Feature value (Pneumotach)');
    %ylabel('flow:drive'); %xlabel(['Pneumo \^2']);
    title('Squared');
    
    subplot(3,3,4); set(gca,'box','off','tickdir','out','fontname','arial narrow');
    scatter(TnTPnasal.Amatrix2_pnasal_matched(:,Ftr), Gtest_matched,1,'filled','markerfacealpha', facealpha, 'markerfacecolor', facecolor); hold on;
    xlabel('Feature value (Pnasal)');
    ylabel('flow:drive'); %xlabel(['Pnasal \^1']);
    
    subplot(3,3,5); set(gca,'box','off','tickdir','out','fontname','arial narrow');
    scatter(TnTPnasal.Amatrix2_pnasal_matched(:,Ftr+50), Gtest_matched,1,'filled','markerfacealpha', facealpha, 'markerfacecolor', facecolor); hold on;
    xlabel('Feature value (Pnasal)');
    %ylabel('flow:drive'); %xlabel(['Pnasal \^0.5']);
    
    subplot(3,3,6); set(gca,'box','off','tickdir','out','fontname','arial narrow');
    scatter(TnTPnasal.Amatrix2_pnasal_matched(:,Ftr+100), Gtest_matched,1,'filled','markerfacealpha', facealpha, 'markerfacecolor', facecolor); hold on;
    xlabel('Feature value (Pnasal)');
    %ylabel('flow:drive'); %xlabel(['Pnasal \^2']);
    
    subplot(3,3,7); set(gca,'box','off','tickdir','out','fontname','arial narrow');
    scatter(TnTPnasal.Amatrix2_flow_matched(:,Ftr), TnTPnasal.Amatrix2_pnasal_matched(:,Ftr),1,'filled','markerfacealpha', facealpha, 'markerfacecolor', facecolor); hold on;
    xlabel('Feature value (Pneumotach)'); ylabel('Feature value (Pnasal)');
    
    subplot(3,3,8); set(gca,'box','off','tickdir','out','fontname','arial narrow');
    scatter(TnTPnasal.Amatrix2_flow_matched(:,Ftr+50), TnTPnasal.Amatrix2_pnasal_matched(:,Ftr+50),1,'filled','markerfacealpha', facealpha, 'markerfacecolor', facecolor); hold on;
    xlabel('Feature value (Pneumotach)'); %ylabel('Feature value (Pnasal)');
    
    subplot(3,3,9); set(gca,'box','off','tickdir','out','fontname','arial narrow');
    scatter(TnTPnasal.Amatrix2_flow_matched(:,Ftr+100), TnTPnasal.Amatrix2_pnasal_matched(:,Ftr+100),1,'filled','markerfacealpha', facealpha, 'markerfacecolor', facecolor); hold on;
    xlabel('Feature value (Pneumotach)'); %ylabel('Feature value (Pnasal)');
    
    
    FtrName = regexprep(FeatureNames{Ftr,2},'[_,:{}]',' ');
    suptitle(FtrName);
    
    str = ['..\Figures\Univariate\Flow_Ftr_', num2str(Ftr), '_', char(FtrName)]; %
    
    switch savefigas
        case 'saveasTIFF'; print(fig, str, '-dtiff', '-r1000');
        case 'saveasPNG'; saveas(fig, str, 'png');
        case 'saveasFIG'; savefig(str);
    end
    if closefigs; close(fig); end
    
end

%%




%%

if 0
    %% test performance of selected features
    Ftrs = TopFtrs(1:2);
    FtrVals = Amatrix2(:,Ftrs);
    Labels = LabelsOrderedOpt;
    Yvariable = Gtest_All;
    
    %% simple linear model
    %Mdl = fitcdiscr(FtrVals(1:500:end,:),Yvariable(1:500:end));%,'Weights',weights, 'DiscrimType','linear');
    Mdl = fitlm(FtrVals(1:5:end,:),Yvariable(1:5:end));%,'Weights',weights, 'DiscrimType','linear');
    PredY = predict(Mdl, FtrVals);
    
    %% Boxplots
    % warning, potentially lots and lots of boxplots
    k=18;
    group = [ones(nnz(Yvariable==0),1)*1;ones(nnz(Yvariable==1),1)*2];
    for outer = 1:11
        figure(504+outer); clf(figure(504+outer));
        for inner = 1:9
            if k==116
                break
            end
            subplot(3,3,inner);
            data = [table2array(PtData(Yvariable==0,k)); table2array(PtData(Yvariable==1,k))];
            boxplot(data,group); hold on;
            refline(0,nanmedian(table2array(PtData(Yvariable==0,k))));
            refline(0,nanmedian(table2array(PtData(Yvariable==1,k))));
            str = [num2str(k)];
            xlabel(str);
            k = k+1;
        end
    end
    
    % handy figure closing
    for outer = 1:11
        close(figure(504+outer));
    end
    
    %% >>>>>>>>>>>>>>>>>>> Univariate analysis <<<<<<<<<<<<<<<<<<<<<<<<<<<
    ShowFigures = 0;
    ftrs=18:1:115;
    AUCData=NaN(size(ftrs,2),7);
    for fts=[ftrs]
        x = PtData{:,fts};
        if any(isinf(x))%isnan(x)|
            continue % this includes inf etc
        end
        Mdl_ = fitcdiscr(x,Yvariable,'Weights',weights);
        Mdl_Pred = predict(Mdl_, x);
        if nnz(Mdl_Pred)==0 || nnz(~Mdl_Pred)==0
            continue
        else
            [thresX,AUC,SEM,p,posclass,sensspec]=ROCAUCSEM(Mdl_Pred,Yvariable,ShowFigures); %
            AUCData(fts-17,1)=thresX;
            AUCData(fts-17,2)=AUC;
            AUCData(fts-17,3)=SEM;
            AUCData(fts-17,4)=p;
            AUCData(fts-17,5)=posclass;
            AUCData(fts-17,6)=sensspec(1);
            AUCData(fts-17,7)=sensspec(2);
        end
    end
    
    AUCData(:,1)=ftrs;
    AUCDataSorted=sortrows(AUCData, 2, 'descend');
    exclude_=isnan(AUCDataSorted(:,2))|AUCDataSorted(:,6)==0|AUCDataSorted(:,7)==0;
    AUCDataSortedTrimmed=AUCDataSorted;
    AUCDataSortedTrimmed(exclude_,:)=[];
    AUCDataSortedTrimmed(:,[4,5])=[];
    
end





%% Load subject waveform data for VEVdrive, Flow and Edi plot
if (nnz(predy<0)>0); predy(predy<0)=0; end
if nnz(predy>maxG)>0; predy(predy>maxG)=1.5; end

% read spreadsheet
AnalyzeDataSpreadsheet = 'AnalyzeDataSpreadsheet.xlsx';
[num,patients,~] = xlsread(AnalyzeDataSpreadsheet,1,'F3:G56');

% savefigs = 1; closefigs = 1;
snip = 10; %set snip number, snip 10 is for paper figures
%snip = 11; % for ATS 2018 talk
%snip = 12; % for ATS 2018 talk

ShowHumanScoring = 0; % 0 adds nothing to plot, 1 adds Toms DOS scoring
% this initial import of DOS scoring is best just for viewing

%for n= [ 2 6 14 17 18 21 28 29 34 44 45 53 ] % full list
%for n= [ 14 29 34 45]   % list of those with at least 2 regions of interest
%for n= [ 14 34 45]      % list of those with at least 3 regions of interest
%for n = [ 14]            % list of those with at least 5 regions of interest
%for n = 5 % central
%for n = 5 % interesting
%for n = [2 5 14 21 45]
% for n = 45
for n = 10 % 10 is interesting too, a bit fl, but a bit central too.
    
    %     if ~ismember(n, PT_list)
    %         continue
    %     end
    
    clearvars Edi Flow StarttimeSpike Data1
    % Plot FL values over time for an example subject (e.g. compare with Spike data). Random check (n=8, 1313) Looks great.
    % n=3; % 1313 is now 14
    % [2;3;4;5;6;8;9;10;11;14;16;17;18;19;20;21;22;23;24;26;27;28;29;30;32;33;34;35;36;38;39;40;41;43;44;45;46;47;50;53;54]
    studyname = char(patients{n});
    studynameNoExt = char(patients{n}(1:end-4));
    %load('J:\PEOPLE\FACULTY\SANDS\Phenotyping from Pes and Flow\SpikeFiles\1313.mat', 'Edi', 'Flow');
    load(['C:\PSG_Data\FlowDrive\SourceMat 20171123\',studyname], 'Edi', 'Flow', 'StarttimeSpike');
    
    try
        load(['C:\PSG_Data\FlowDrive\Converted\',studynameNoExt,'_XHz.mat']);
    catch GetEvtsFail
        disp(GetEvtsFail.getReport);
    end
    EventsAr = DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'EventsAr')==1));
    EventsResp = DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'EventsResp')==1));
    % arousal   = 1
    % ap-O      = 2
    % ap-C      = 3
    % hyp-O     = 4
    % mixed     = 5
    % hyp-C     = 6
    
    %ArColor
    %RespOColor
    %RespCColor
    
    if ~(exist('StarttimeSpike', 'var') == 1)
        StarttimeSpike = 0;
    end
    Time = StarttimeSpike:0.008:StarttimeSpike+(length(Flow.values)-1)*0.008;
    %PTtime = PtData.BB_time(PtData.PT==n);
    %PTtime = PTtime-PTtime(1)+StarttimeSpike;
    
    Isubj = PtData.PT==n;
    %Data1 = [PtData.BB_time(PtData.PT==n) PtData.BB_Ttot(PtData.PT==n) predy(PtData.PT==n) Gtest_All(PtData.PT==n)];
    
    if ShowHumanScoring
        [Gtest_pt, PredY_pt, ManScore_pt] = ImportDOSscoring(n,PtData,Gtest_All,predy);
        Data1 = [PtData.BB_time(Isubj) PtData.BB_Ttot(Isubj) PredY_pt Gtest_pt ManScore_pt];
    else
        Data1 = [PtData.BB_time(Isubj) PtData.BB_Ttot(Isubj) predyL1O_array(Isubj,ftrnum) Gtest_All(Isubj)];
    end
    
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
    
    if Flow.length<=length(Time)
        Flow.values(end:length(Time))=0;
    else
        Flow.values(length(Time):end)=[];
    end
    if Edi.length<=length(Time)
        Edi.values(end:length(Time))=0;
    else
        Edi.values(length(Time):end)=[];
    end
    
    % Data columns: (1) Time, (2) BBTtot, (3) PredY, (4) Gtest. and
    % optionally (5) human scoring
    
    % Clip Gtest to 1.05
    Data1(:,4) = min(Data1(:,4), 1.05);
    
    dsf=5; dt=Flow.interval;
    FlowF=Flow.values;
    if 1
        filter_HFcutoff_butter0 = 12.5;
        filter_order0 = 1;
        [B_butter0,A_butter0] = butter(filter_order0,[filter_HFcutoff_butter0]/(1/dt/2),'low');
        FlowF = filtfilt(B_butter0,A_butter0,FlowF); %filtfilt, otherwise flow signal is right-shifted
    end
    
    %% make the plot for the subject data
    figure(200+n); clf(200+n); fig = gcf;
    fig.Color = [1 1 1]; fig.Units = 'inches';
    fig.Position = [-12.5 1 12 7.5];
    % position [ x0 y0 width height]
    
    TickFntSz = 12;
    threeplot = 0;
    if ~threeplot
        ax(1)=subplot(4,1,1); set(gca,'Position',[0.06 0.94 0.90 0.03]);
        % turn continuous time Arousal events into list with duration
        EvtsAr = downsample(EventsAr,dsf);
        Timeds = downsample(Time,dsf);
        I2 = find(diff(EvtsAr)==-1);
        I1 = find(diff(EvtsAr)==1);
        [I1,I2] = TidyStartEndEventList(I1,I2,length(Timeds));
        It1 = Timeds(I1); It2 = Timeds(I2);
        currentAxes = gca; % switch to current axes
        % then mark these on the plot
        for pn = 1:length(It1)
            BBmarker = patch([It1(pn) It2(pn) It2(pn) It1(pn)],...
                [0 0 1 1], [.8 .8 .8],...
                'FaceAlpha',0.3, 'EdgeColor','none', 'Parent',currentAxes);%
        end
        set(gca,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz);
        set(gca,'xtick',[],'box','off');
        set(gca,'xcolor',[1 1 1])
        set(gca,'yticklabel',[]);
        set(gca,'ytick',[]);
        ylabel('Ar', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
    end
    
    if threeplot
        ax(1)=subplot(3,1,1); set(gca,'Position',[0.06 0.68 0.90 0.28]);
    else
        ax(2)=subplot(4,1,2); set(gca,'Position',[0.06 0.64 0.90 0.27]);
    end
    stairs(Data1(:,1),100*Data1(:,4),'k'); % actual
    hold('on')
    stairs(Data1(:,1),100*Data1(:,3),'r'); % pred
    if ShowHumanScoring
        stairs(Data1(:,1),100*Data1(:,5),'g'); % human score
    end
    plot([Data1(1,1) Data1(end,1)],[0 0],'k:');
    plot([Data1(1,1) Data1(end,1)],100*[1 1],'k:');
    plot([Data1(1,1) Data1(end,1)],50*[1 1],'k:');
    set(gca,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz);
    set(gca,'xtick',[],'box','off');
    set(gca,'xcolor',[1 1 1])
    %axc = gca; axc.YTick=[0:50:100]; yticklabels(axc, {'0', '50', '100'});
    axc = gca; axc.YTick=[0:25:100]; yticklabels(axc, {'0','25','50','75','100'});
    ylabel('{\itflow:drive} (%)', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
    ylim([-5 105]);
    
    if threeplot
        % turn continuous time Arousal events into list with duration
        EvtsAr = downsample(EventsAr,dsf);
        Timeds = downsample(Time,dsf);
        I2 = find(diff(EvtsAr)==-1);
        I1 = find(diff(EvtsAr)==1);
        [I1,I2] = TidyStartEndEventList(I1,I2,length(Timeds));
        It1 = Timeds(I1); It2 = Timeds(I2);
        currentAxes = gca; % switch to current axes
        % then mark these on the plot
        for pn = 1:length(It1)
            BBmarker = patch([It1(pn) It2(pn) It2(pn) It1(pn)],...
                [105 105 115 115], [.8 .8 .8],...
                'FaceAlpha',0.3, 'EdgeColor','none', 'Parent',currentAxes);%
        end
    end
    
    if threeplot
        %ax(2)=subplot(3,1,2); plot(Time,Flow.values,'g'); hold('on');
        ax(2)=subplot(3,1,2); set(gca,'Position',[0.06 0.36 0.90 0.28]);
    else
        ax(3)=subplot(4,1,3); set(gca,'Position',[0.06 0.34 0.90 0.27]);
    end
    plot(downsample(Time,dsf),downsample(FlowF,dsf),'k'); hold on
    plot(downsample(Time,dsf),zeros(length(downsample(FlowF,dsf)),1),'k:');
    set(gca,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz);
    set(gca,'xtick',[],'box','off');
    set(gca,'xcolor',[1 1 1])
    ylabel('Flow (L/s)', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
    if 0; ylim([-0.8 0.8]); else; ylim([-1.1 1.1]); end
    
    
    % turn continuous time Resp events into list with duration
    EvtsResp = downsample(EventsResp,dsf);
    Timeds = downsample(Time,dsf);
    I2 = find(diff(EvtsResp)==-4);  % 4 is picking up hyp-O only, could expand to include more evt types
    I1 = find(diff(EvtsResp)==4);
    % enable this if doing central hypop
    if n == 5
        I2 = find(diff(EvtsResp)==-6);  % 6 is hyp-C
        I1 = find(diff(EvtsResp)==6);
    end
    % enable this if doing obs ap
    if n == 21
        I2 = find(diff(EvtsResp)==-2);  % 2 is Ap-O
        I1 = find(diff(EvtsResp)==2);
    end
    
    [I1,I2] = TidyStartEndEventList(I1,I2,length(Timeds));
    It1 = Timeds(I1);  It2 = Timeds(I2);
    currentAxes = gca; % switch to current axes
    ylimits = currentAxes.YLim; % call this explicitly, fails otherwise
    % then mark these on the plot
    for pn = 1:length(It1)
        BBmarker = patch([It1(pn) It2(pn) It2(pn) It1(pn)],...
            [ylimits(1)+0.1 ylimits(1)+0.1 ylimits(2)-0.1 ylimits(2)-0.1], [.8 .8 .8],...
            'FaceAlpha',0.3, 'EdgeColor','none', 'Parent',currentAxes);%
    end
    
    %plot(downsample(Time,dsf),-downsample(EventsAr,dsf),'r');
    %plot(downsample(Time,dsf),downsample(EventsResp,dsf),'g');
    if threeplot
        ax(3)=subplot(3,1,3); set(gca,'Position',[0.06 0.04 0.90 0.28]);
    else
        ax(4)=subplot(4,1,4); set(gca,'Position',[0.06 0.04 0.90 0.27]);
    end
    plot(downsample(Time,dsf),downsample(Edi.values,dsf),'k'); hold on;
    set(gca,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz);
    set(gca,'xtick',[],'box','off');
    set(gca,'xcolor',[1 1 1])
    axc = gca; axc.YTick=[0:15:60]; yticklabels(axc, {'0','15','30','45','60'});
    ylabel('EMGdi (uV)', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
    ylim([-2 62]);
    
    linkaxes(ax,'x');
    %suptitle([studyname]);
    
    %% mods to plot
    switch n
        case 2
            switch snip
                case 1
                    v1 = 27000; v2 = 27350;
                case 10
                    v1 = 27025; v2 = 27145;
            end
            
            useStart = 1;
        case 5 % central hypop and SHR
            switch snip
                case 1
                    v1 = 21660; v2 = 21700; % central  % export two minute window, clip later
                    v1 = 22950; v2 = 23150; % Full
                    v1 = 22950; v2 = 22965; % 1 of 5
                    v1 = 22980; v2 = 22995; % 2 of 5
                    v1 = 23010; v2 = 23025; % 3 of 5
                    v1 = 23055; v2 = 23070; % 4 of 5
                    v1 = 23075; v2 = 23090; % 5 of 5
                    useStart = 1;
                case 10
                    v1 = 21650; v2 = 21770; % central  % export two minute window, clip later
                    useStart = 1;
            end
        case 6
            v1 = 2900; v2 = 3600;
            useStart = 1;
        case 10
            v1 = 19500; v2 = 19620;
            useStart = 1;
            
        case 14 % off the scale, big flow and massive edi, still "flow limited"
            % multiple
            S14 = StarttimeSpike; StarttimeSpike = S14;
            if threeplot
                set(ax(2),'ylim', [-1.75 1.75]);
                set(ax(3),'ylim', [-2 60]);
            else
                set(ax(3),'ylim', [-1.75 1.75]);
                set(ax(4),'ylim', [-2 60]);
            end
            switch snip
                case 1
                    v1 = 5200; v2 = 6400; % big window
                    useStart = 1;
                case 2
                    v1 = 85487; v2 = 85658; % subclinical events
                    useStart = 0;
                case 3
                    v1 = 86027; v2 = 86198; % stable FL
                    useStart = 0;
                case 4
                    v1 = 86385; v2 = 86556; % stable FL
                    useStart = 0;
                case 5
                    v1 = 23560;  v2 = 23680; % maybe no good
                    useStart = 1;
                case 10
                    v1 = 86060; v2 = 86180; % stable FL
                    useStart = 0;
            end
            
        case 17
            v1 = 83450; v2 = 83700;
            useStart = 0;
        case 18
            v1 = 105600; v2 = 106000;
            useStart = 0;
        case 21
            switch snip
                case 1
                    %v1 = 83900; %v2 = 85000; % not defined
                    v1 = 84450; v2 = 84820;
                    useStart = 0;
                case 10
                    v1 = 84670; v2 = 84790;
                    useStart = 0;
            end
        case 28
            v1 = 19460; v2 = 19580;
            useStart = 1;
        case 29
            % multiple
            S29 = StarttimeSpike; StarttimeSpike = S29;
            switch snip
                case 1
                    v1 = 13400; v2 = 13600;
                case 2
                    v1 = 15800; v2 = 16000;
            end
            useStart = 1;
        case 34
            % multiple
            switch snip
                case 1
                    v1 = 90410; v2 = 90560;
                case 2
                    v1 = 100000;  v2 = 100400;
                case 3
                    v1 = 101705; v2 = 102100;
            end
            useStart = 0;
        case 44
            v1 = 96500; v2 = 97600;
            useStart = 0;
        case 45
            % multiple
            switch snip
                case 1
                    v1 = 86450; v2 = 87000;
                case 2
                    v1 = 91380; v2 = 91600;
                case 3
                    v1 = 91450; v2 = 91515;
                case 10
                    v1 = 91391; v2 = 91511;  % 2 minutes
                case 11
                    v1 = 91451; v2 = 91511;  % 1 minute
                case 12
                    v1 = 91480; v2 = 91490;  % 10 seconds, uses if 0 below, fig 2450
            end
            useStart = 0;
        case 53
            v1 = 103200;  v2 = 104200;
            useStart = 0;
    end
    if useStart
        xlim([StarttimeSpike+v1 StarttimeSpike+v2])
    else
        xlim([v1 v2]);
    end
    % xlim([-inf inf]);
    %%
    if 0
        figure(2450); clf(figure(2450)); fig = gcf;
        fig.Color = [1 1 1]; fig.Units = 'inches';
        fig.Position = [-12  1   10    7];
        plot(downsample(Time,dsf),downsample(FlowF,dsf),'k'); hold on
        plot(downsample(Time,dsf),zeros(length(downsample(FlowF,dsf)),1),'k:');
        set(gca,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz);
        set(gca,'xtick',[],'box','off');
        set(gca,'xcolor',[1 1 1])
        ylabel('Flow (L/s)', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
        ylim([-0.61 0.61]);    xlim([v1 v2]);
        str = ['..\Figures\SampleData\SampleData_Pt ', num2str(n), ', snip ', num2str(snip), 'singleBreath' ]; %
        print(fig, str, '-dtiff', '-r1000');
    end
    
    %%
    str = ['..\Figures\SampleData\SampleData_Pt ', num2str(n), ', snip ', num2str(snip)]; % ,', 3 of 5 '
    %str = ['..\Figures\Figure_2'];
    savefigs = 0;
    if savefigs; saveas(fig, str, 'png'); savefig(str); end %
    
    switch savefigas
        case 'saveasTIFF'; print(fig, str, '-dtiff', '-r1000');
        case 'saveasPNG'; saveas(fig, str, 'png');
        case 'saveasFIG'; savefig(str);
    end
    
    if closefigs; close(fig); end
    
end

disp('done');

%% Add labels A B to plot space
% this only works for selected pt...
subplot(3,1,1); hold on;
text(91400, 99, 'A', 'FontSize', 20, 'FontWeight', 'Bold');
subplot(3,1,2); hold on;
text(91400, 0.6, 'B', 'FontSize', 20, 'FontWeight', 'Bold');
subplot(3,1,3); hold on;
text(91400, 45, 'C', 'FontSize', 20, 'FontWeight', 'Bold');

%%
xlim([StarttimeSpike, inf]);
AHI_perPT(45)


%%
fig = gcf;
str = ['..\Figures\SampleData\SampleData_Pt 45 snip 4 withEvts'];
saveas(fig, str, 'png');


%% re-run MIFL on selected areas
global settings
settings.scalingexponent = 0.666666667; % 2/3
settings.sqrt_scaling = 0;
settings.plotfigure = 0;

[time,Vflow,BB_i_start,BB_i_mid,BB_i_end,BB_t,VI,VE,Ttot_B,leak,...
    IEratio,VT,Vpeak,Vpeakmean,Apnea_B,Vflow_out,VTi,VTe] =...
    VEfromFlow_sqrt_V16(Time,FlowF);

figure(77); clf(figure(77));
plot(time, Vflow); hold on
plot(time(BB_i_start), Vflow(BB_i_start), 'rx');
plot(time(BB_i_mid), Vflow(BB_i_mid), 'rv');
plot(time(BB_i_end), Vflow(BB_i_end), 'r.');
xlim([v1 v2]);

%%
t2 = time(BB_i_start);
loc1 = find(t2>91515,1,'first');
loc2 = find(t2>91600,1,'first');

BB_i_start_ = BB_i_start(loc1:loc2);
BB_i_mid_ = BB_i_mid(loc1:loc2);
BB_i_end_ = BB_i_end(loc1:loc2);
BB_original = [BB_i_start_,BB_i_mid_,BB_i_end_];
BB_Ttrans = [];
TiTrans = ones(length(BB_i_start),1);
TeTrans = ones(length(BB_i_start),1);
TAA_ = zeros(length(time),1);
MIFLsettings = [0, 1];
[BreathDataTable] = MIFL(time, Vflow, BB_original, BB_Ttrans, TiTrans, TeTrans, TAA_, Apnea_B, MIFLsettings);

%%










%%










%% redo All or Sleep only breaths
SleepOnly = 1;
IncludePnasal = 1;
[BB, BB_] = getAllorSleepOnlyBB(SleepOnly, IncludePnasal, PtData, PtData_matched);

nnz(BB)
if IncludePnasal
    nnz(BB_)
end



%% Pnasal vs Flow plots - compare pnasal and flow - not including Ap and LF BB
Pnasal_summary_table = array2table(Pnasal_summary, 'VariableNames', {'PT','FlowBB','PnasalBB','MatchedBB','PercentMatched'});
nnz(BB_)

% scatter of flow predy vs pnasal predy - this pnasal predy uses flow betas
figure(23);clf(figure(23)); fig = gcf;
fig.Color = [1 1 1]; fig.Units = 'inches';
fig.Position = [-12.2   4   4.5    4.5];

%subplot(1, 2, 1);
scatter(predyL1O_array_pnasal_matched(BB_,ftrnum).*100, predyL1O_array_flow_matched(BB_,ftrnum).*100, 2,'filled','markerfacealpha',0.3);
xlabel('Pnasal Predicted {\itflow:drive} (%)');ylabel('Flow Predicted {\itflow:drive} (%)');
str=['Predicted {\itflow:drive} with alternate airflow signals']; title(str);
r = corr(predyL1O_array_pnasal_matched(BB_,ftrnum), predyL1O_array_flow_matched(BB_,ftrnum));
lsline()
% subplot(1, 2, 2);
% scatter(VEVeup_pnasal_array(BB_,ftrnum), VEVeup_flow_array(BB_,ftrnum), 2,'filled','markerfacealpha',0.3);
% xlabel('Pnasal VEVeup');ylabel('Flow VEVeup');
% str=['Flow based VS Pnasal based VEVeup']; title(str);
% xlim([0 4]); ylim([0 4]);

%str_plot=['Pred performance at ', num2str(ftrnum), ' features, VEVeup is all breaths']; suptitle(str_plot);

str = ['..\Figures\Figure_5_1'];
switch savefigas
    case 'saveasTIFF'; print(fig, str, '-dtiff', '-r1000');
    case 'saveasPNG'; saveas(fig, str, 'png');
    case 'saveasFIG'; savefig(str);
end
if closefigs; close(fig); end

%%








%% moving time median
if 0
    Yval = Gtest;
    xdata = [0 medianX];
    ydata = [0 medianY];
    %plot(100*xdata,100*ydata,'r');
    x0 = [1 1];
    upper = [1.1 2];
    lower = [0.9 0.5];
    fun = @(x,xdata)x(1)*(xdata.^x(2));
    [Xmodel]=lsqcurvefit(fun,x0,xdata,ydata,lower,upper);
    xin = 0:0.01:1.5;
    yout = fun(Xmodel,xin);
    plot(100*xin,100*yout,'r');
    predybackup = predy;
    predy = fun(Xmodel,predy);
    
    subplot(2,2,4);
    scatter(100*predy,100*Yval,2,'filled','markerfacealpha',0.4)
    xlim([0 150]);
    hold('on')
    
    %dx=0.1; xbins=[0 0.15:dx:1.05 1.5];
    dx=0.2; xbins=[0 0.3:dx:0.9 1.5];
    title(['Rvalue: ' num2str(Rvalue(end),2)])
    %plot binned data
    clear medianX medianY upperIQRY lowerIQRY upperIQRY2 lowerIQRY2
    for i=1:length(xbins)-1
        Ix=predy>xbins(i)&predy<xbins(i+1);
        medianX(i)=prctile(predy(Ix),50);
        medianY(i)=prctile(Yval(Ix),50);
        upperIQRY(i)=prctile(Yval(Ix),75);
        lowerIQRY(i)=prctile(Yval(Ix),25);
        upperIQRY2(i)=prctile(Yval(Ix),90);
        lowerIQRY2(i)=prctile(Yval(Ix),10);
    end
    %medianX = [0.2:dx:1]; % overwrite
    plotbarwitherrorsmedian(100*medianY,100*medianX,100*upperIQRY2,100*lowerIQRY2,zeros(length(medianY),3),NaN,100*0.005,100*0.01)
    plotbarwitherrorsmedian(100*medianY,100*medianX,100*upperIQRY,100*lowerIQRY,zeros(length(medianY),3),NaN,100*dx/2,100*0.01)
end


%% Patient level summary data
% Without adding Ap BB back in
%unique(PtData.PT)
for i=1:54%max(PtData.PT)
    if ~ismember(i, PT_list)
        str=['No data for Pt ', num2str(i)]; disp(str); continue
    end
    
    I=PtData.PT==i&PtData.NotAr==1;
    str=[num2str(i), ', ', num2str(nnz(I))]; disp(str);
    medianG(i) = nanmedian(Gtest_All(I));
    medianGest(i) = nanmedian(predy(I));
    
    denominator = sum(~isnan(Gtest_All(I)));
    
    Fmildmod(i) = sum(Gtest_All(I)<0.9&Gtest_All(I)>0.5)/denominator;
    Fmildmodest(i) = sum(predy(I)<0.9&Gtest_All(I)>0.5)/denominator;
    
    Fmodsev(i) = sum(Gtest_All(I)<0.7)/denominator;
    Fmodsevest(i) = sum(predy(I)<0.7)/denominator;
    
    Fsev(i) = sum(Gtest_All(I)<0.5)/denominator;
    Fsevest(i) = sum(predy(I)<0.5)/denominator;
end

% With adding ApBB back in
% ToDo: with ap breaths put back in
for i=1:54%max(PtData.PT)
    I=PtData.PT==i&PtData.NotAr==1;
    if nnz(I)==0
        continue
    end
    
    indAP = find(RemovedBB_Apnoea.Pt == i);
    indLF = find(RemovedBB_LowFlow.Pt == i);
    if ~isempty(indAP)
        numclinscored = numclinscored+sum(RemovedBB_Apnoea{indAP,[4 6 7]});
        numPredBBbelowthreshold = numPredBBbelowthreshold+RemovedBB_Apnoea{indAP,2};
        numActualBBbelowthreshold = numActualBBbelowthreshold+RemovedBB_Apnoea{indAP,2};
        numPredBBbelowthreshold2 = numPredBBbelowthreshold2+RemovedBB_Apnoea{indAP,2};
        numActualBBbelowthreshold2 = numActualBBbelowthreshold2+RemovedBB_Apnoea{indAP,2};
        numBBtotal = numBBtotal+RemovedBB_Apnoea{indAP,2};
    end
    if ~isempty(indLF)
        
        numPredBBbelowthreshold = numPredBBbelowthreshold+RemovedBB_LowFlow{indLF,2};
        numActualBBbelowthreshold = numActualBBbelowthreshold+RemovedBB_LowFlow{indLF,2};
        numPredBBbelowthreshold2 = numPredBBbelowthreshold2+RemovedBB_LowFlow{indLF,2};
        numActualBBbelowthreshold2 = numActualBBbelowthreshold2+RemovedBB_LowFlow{indLF,2};
        numBBtotal = numBBtotal+RemovedBB_LowFlow{indLF,2};
    end
    
    
    medianG(i) = nanmedian(Gtest_All(I));
    medianGest(i) = nanmedian(predy(I));
    
    denominator = sum(~isnan(Gtest_All(I)));
    
    Fmildmod(i) = sum(Gtest_All(I)<0.9&Gtest_All(I)>0.5)/denominator;
    Fmildmodest(i) = sum(predy(I)<0.9&Gtest_All(I)>0.5)/denominator;
    
    Fmodsev(i) = sum(Gtest_All(I)<0.7)/denominator;
    Fmodsevest(i) = sum(predy(I)<0.7)/denominator;
    
    Fsev(i) = sum(Gtest_All(I)<0.5)/denominator;
    Fsevest(i) = sum(predy(I)<0.5)/denominator;
end

exclude = isnan(medianG);
medianGest(exclude) = [];
medianG(exclude) = [];
Fmodsevest(exclude) = [];
Fmodsev(exclude) = [];
Fsevest(exclude) = [];
Fsev(exclude) = [];
Fmildmodest(exclude) = [];
Fmildmod(exclude) = [];

corr(medianGest',medianG')
corr(Fmodsevest',Fmodsev')
corr(Fsevest',Fsev')
corr(Fmildmodest',Fmildmod')

figure(5); clf(5);
subplot(1,4,1);
[Rtemp,Ptemp]=plotregressionwithSEM(medianGest,medianG); title(num2str(Rtemp)); xlabel('median');
subplot(1,4,2);
[Rtemp,Ptemp]=plotregressionwithSEM(Fmodsevest,Fmodsev); title(num2str(Rtemp)); xlabel('modsev');
subplot(1,4,3);
[Rtemp,Ptemp]=plotregressionwithSEM(Fsevest,Fsev); title(num2str(Rtemp)); xlabel('sev');
subplot(1,4,4);
[Rtemp,Ptemp]=plotregressionwithSEM(Fmildmodest,Fmildmod); title(num2str(Rtemp)); xlabel('mildmod');

PredY_avg; Gtest_avg;
PredY_thres1_pFL; Gtest_thres1_pFL; % Proportion sev FL
PredY_thres2_pFL; Gtest_thres2_pFL; % Proportion modsev FL

r_temp = corr(medianGest',medianG');


%% Correlate with AHI
AHItotal = AHI_perPT(~isnan(AHI_perPT(:,1)),1);
figure(6); clf(6);
subplot(1,1,1);
[Rtemp,Ptemp]=plotregressionwithSEM(AHItotal',medianGest); title(num2str(Rtemp));

[b,dev,stats]=glmfit([AHItotal medianGest'],medianG');
stats.p(3)
stats.coeffcorr(1,3)
[b,dev,stats]=glmfit([AHItotal Fmodsevest'],Fmodsev');
stats.p(3)
stats.coeffcorr(1,3)
[b,dev,stats]=glmfit([AHItotal Fsevest'],Fsev');
stats.p(3)
stats.coeffcorr(1,3)

%% Repeat, sensitivity for Nfeatures
PT_list = unique(PtData.PT);
for n=30 % 1:100
    clear medianG medianGest Fmildmod Fmildmodest Fmodsev Fmodsevest Fsev Fsevest
    predy = predyL1O_array(:,n);
    
    for i=1:54
        if ~ismember(subj, PT_list)
            continue
        end
        I=PtData.PT==i&PtData.NotAr==1;
        medianG(i) = nanmedian(Gtest_All(I));
        medianGest(i) = nanmedian(predy(I));
        
        denominator = sum(~isnan(Gtest_All(I)));
        
        Fmildmod(i) = sum(Gtest_All(I)<0.9&Gtest_All(I)>0.5)/denominator;
        Fmildmodest(i) = sum(predy(I)<0.9&Gtest_All(I)>0.5)/denominator;
        
        Fmodsev(i) = sum(Gtest_All(I)<0.7)/denominator;
        Fmodsevest(i) = sum(predy(I)<0.7)/denominator;
        
        Fsev(i) = sum(Gtest_All(I)<0.5)/denominator;
        Fsevest(i) = sum(predy(I)<0.5)/denominator;
    end
    
    %corr(Fsevest(PT_list)',Fsev(PT_list)')
    
    Rvals(n,:)=[...
        corr(medianGest',medianG'),...
        corr(Fmodsevest',Fmodsev'),...
        corr(Fsevest',Fsev'),...
        corr(Fmildmodest',Fmildmod'),...
        ];
end




%%

figure(7); clf(figure(7));
subplot(1,2,1);
mfl=100*Data1(:,4);
mfl_categories2_ = 100*[0:0.05:1.51];
for jj=2:length(mfl_categories2_)
    mfl_categories2(n,jj-1)=sum((mfl>=mfl_categories2_(jj-1)&mfl<mfl_categories2_(jj)))/length(mfl);
end

bar((mfl_categories2_(1:end-1)+0.025),mfl_categories2(n,:));
box('off');
xlim([0 150]);
ylim([0 0.14]);

subplot(1,2,2);
mfl=100*Data1(:,3);
mfl_categories2_ = 100*[0:0.05:1.51];
for jj=2:length(mfl_categories2_)
    mfl_categories2(n,jj-1)=sum((mfl>=mfl_categories2_(jj-1)&mfl<mfl_categories2_(jj)))/length(mfl);
end

bar((mfl_categories2_(1:end-1)+0.025),mfl_categories2(n,:));
box('off');
xlim([0 150]);
ylim([0 0.14]);
% a very simplistic way to see the tow trend lines
%plot(PtData.BB_time(PtData.PT==n),smooth(Gtest(PtData.PT==n),99),'c-');
%plot(PtData.BB_time(PtData.PT==n),smooth(Gpredicted(PtData.PT==n),99),'m-');
%end


%% Script to explore differences in coefficient p values and loglikelihood ratio test p values for coefficient inclusion
%Results: P values are similar, but not identical
%Appear to approach identity at large N
Ftrs=[1:20];
rangei=1:8000;
mdlu=fitglm(Amatrix2(rangei,Ftrs),Gtest_(rangei),'Distribution','binomial','weights',weights(rangei));
mdlr=fitglm(Amatrix2(rangei,Ftrs(1:end-1)),Gtest_(rangei),'Distribution','binomial','weights',weights(rangei));%
pValueCoefficient=mdlu.Coefficients.pValue(end) %displayed
uLogL=mdlu.LogLikelihood;
rLogL=mdlr.LogLikelihood;
dof=1;
[h,pValueLogLikelihood,stat,cValue] = lratiotest(uLogL,rLogL,dof);
pValueLogLikelihood %displayed

%% Can we speed this up by not making the covariance matrix?
%Result: Not really-.

rangei=1:40000;
I=sum(isnan([Amatrix2 Gtest_ weights]),2)>0;
Amatrix2_=Amatrix2;
Gtest2_=Gtest_;
weights_=weights;
Amatrix2_(I,:)=[];
Gtest2_(I)=[];
weights_(I)=[];
I2=sum(isnan(Amatrix2_),2)>0;
Ftrs=[1:200];
tic
[b,dev,stats]=glmfit(Amatrix2_(rangei,Ftrs),Gtest2_(rangei),'binomial','weights',weights_(rangei));
toc
Pvals = stats.p(2:end);

optionsX=statset('glmfit');
optionsX.TolX=1e-02;
tic
[b,dev,stats]=glmfitFastLR(Amatrix2_(rangei,Ftrs),Gtest2_(rangei),'binomial','weights',weights_(rangei),'options',optionsX);
toc
Pvals = stats.p(2:end);

