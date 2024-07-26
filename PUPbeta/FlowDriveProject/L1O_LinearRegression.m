%% start
close all
clear
clc

if 0
    addpath('C:\Users\uqdmann\Dropbox\PUPbeta_git\');
    cd('C:\Users\uqdmann\Dropbox\PUPbeta_git\FlowAndEdiPesQAO');
    load('C:\Users\uqdmann\Dropbox\PUPbeta_git\FeatureSpaces\FlowDrive_only25Hz_FeatureSpace_AutoRef2_Edi_Clean_withPnasalMatched_LinRegWorkspace_20180122.mat')
end

%% options
Use25HzDS = 1;          % set as 1 for 25Hz data, or 0 to use 125Hz data
EdiDrive = 1;           % set as 1 for Edi Drive data, 0 for Pes drive
useLogitFtrs = 1;       % set as 1 to use ftrs selected from logit, 0 for all ftrs
SS_fastfit = 1;			% set as 1 to use SS fast glmfit method, 0 for matlab builtin
TransformTheData = 1;   % set as 1 to do tranforms, or 0 to use unadjusted data
addextratransform = 0; 	% set as 1 to do extra transforms
L1O_run = 1;            % set as 1 to run to do the L1O run, or 0 to skip run
ShowFigures = 1;        % set as 1 to show figures, or 0 to not show figures
IncludePnasal = 1;      % set as 1 to include Pnasal, or 0 to ignore pnasal

datadir = 'C:\PSG_Data\FlowDrive\FeatureSpaces\';
%datadir = '..\FeatureSpaces\';  %

%% open file
if Use25HzDS
    if EdiDrive
        filename = [datadir, 'FlowDrive_only25Hz_FeatureSpace_AutoRef2_Edi_Clean.mat']; datastr = 'FD edi 25Hz';
        %filename = [datadir, 'PnasalDrive_only25Hz_exp05_FeatureSpace_AutoRef2_Edi_Clean.mat']; datastr = 'FlowandPnasalD edi 25Hz';
        if IncludePnasal % for Pnasal
            pnasal_filename = [datadir, 'PnasalDrive_only25Hz_exp067_FeatureSpace_AutoRef2_Edi_Clean.mat']; datastr = 'FlowandPnasalD edi 25Hz';
            filename = [datadir, 'PnasalDrive_only25Hz_Normalized_FeatureSpace_Auto_Edi_Clean.mat'];
        end
    else; filename = [datadir, 'FlowDrive_only25Hz_FeatureSpace_Pes_Clean.mat']; datastr = 'FD pes 25Hz'; end
else
    if EdiDrive; filename = [datadir, 'FlowDrive_only125Hz_FeatureSpace_Edi_Clean.mat']; datastr = 'FD edi 125Hz';
    else; filename = [datadir, 'FlowDrive_only125Hz_FeatureSpace_Pes_Clean.mat']; datastr = 'FD pes 125Hz';end
end


%%
if L1O_run
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
    if useLogitFtrs
        load('FtrsToUse','LabelsOrderedOpt');
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
    
    %% setup the data matrix to use.
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
    
    %% Preparation (Do this carefully)
    Gtest_All = PtData.g_Edi_Adj;
    if IncludePnasal
        Gtest_All_pnasal = PtData_pnasal.g_Edi_Adj; % potentially unused, as we use Flow based Gtest
    end
    NfeatureSteps=100;
    colofones = ones(length(Gtest_All),1);
    predyL1O = NaN*Gtest_All;
    predyL1O_array = NaN*ones(length(Gtest_All),NfeatureSteps);
    
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
    try
        clear RvalueTrain Err ErrRms Nfeatures
        clear RsqTrain_array ErrTrain_array ErrRmsTrain_array
        labels_Step_Subj =[];
        alpha=0.05;
        MaxNfeatures=1;
        neverbreak=1;
        t_start_L1O = clock;
        PT_list = unique(PtData.PT);
        Pnasal_list = [3 5 8 9 18 22 33 34 35 36 38 43 46 47 50 53 54];
        Labels_Complete = Labels; % save a backup of the complete Label set
        Fwd=0;
        VEVeup_All = [];
        if IncludePnasal
            VEVeup_pnasal_array = [];
            VEVeup_flow_array = [];
            Pnasal_summary = [];
            PnasalFlowBB_PTindex=[];
            PtData_matched = [];
            Gtest_matched = [];
            weights_matched = [];
            predyL1O_array_pnasal_matched = [];
            predyL1O_array_flow_matched = [];
            ChannelsList = {'Flow','Pnasal'};
            artdirectory = ['C:\PSG_Data\FlowDrive\SourceMat 20171123'];
            % read spreadsheet (options worksheet)
            [~,~,raw] = xlsread('C:\Users\uqdmann\Dropbox\PUPbeta_git\FlowAndEdiPesQAO\AnalyzeDataSpreadsheet.xlsx',1,'G3:G56');
        end

        for subj=1:54 % 
            %% loop through all subjects, leaving one out each time
            if ~ismember(subj, PT_list)    
                str=['No data for Pt ', num2str(subj)]; disp(str); continue
            end
            
            disp(' '); str=['Performing analysis, witholding Pt ', num2str(subj)]; disp(str); tic;
            
            PtHasPnasal = 0;     
            if ismember(subj, Pnasal_list)    
                if IncludePnasal
                    str=['Processing pnasal data']; disp(str);
                    PtHasPnasal = 1;
                end
            end
           
            %% set up the full length Flow data
            Labels = Labels_Complete; % Start with full list
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
                
%                 load(['C:\PSG_Data\FlowDrive\SourceMat 20171123\' fname], 'Flow','Pnasal','time');
%                 plottime = [0:0.008:(Flow.length/(1/Flow.interval))];
%                 plottime(end) = []; dsf = 10;
%                 figure(3);clf(figure(3));
%                 plot(downsample(plottime,dsf), downsample(Flow.values,dsf), 'color', [0.5 0 0]); hold on;
%                 plot(downsample(plottime,dsf), downsample(-1.*Pnasal.values,dsf), 'color', [0 0 0.5]);
%                 %plot(PtBB_flow_time,-0.4, 'r^'); hold on;
%                 %plot(PtBB_pnasal_time, -0.5, 'b^');
                
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
        
%                 figure(3);
%                 plot(PtBB_flow_time, -0.42, 'r^');
%                 plot(PtBB_pnasal_time, -0.52, 'b^');

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
            
            %continue
            
            %% do backwards elimination
            warning('off');
            maxp=Inf;
            Ftr_indx = 1:size(Amatrix2_train,2);
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
                    if 0 % do at end
                    RsqTrain_array(subj,length(Ftr_indx)) = 1-nansum(weights_train.*(Gtest_train-predytrain).^2)/nansum(weights_train.*(Gtest_train-nanmean(Gtest_train)).^2);
                    ErrTrain_array(subj,length(Ftr_indx)) = nanmean(weights_train.*abs(predytrain-Gtest_train));
                    ErrRmsTrain_array(subj,length(Ftr_indx)) = nanmean((weights_train.*(predytrain-Gtest_train)).^2).^0.5;
                    end                    
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
                                       
                    %% how does this perform in Pnasal?
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
                        
                        if 0
                            figure(90); clf(figure(90));
                            scatter(VEVeup_Flow, VEVeup_pnasal, 3,'filled','markerfacealpha',0.5);hold on;
                            hndl_ls = lsline;  hndl_ls.LineStyle = '-'; hndl_ls.Color = [0 0 1];
                            hndl_ref = refline(1,0); hndl_ref.LineStyle = '--'; hndl_ref.Color = [0.5 0.5 0.5];
                            
                            xlabel('Flow VEVeup');ylabel('Pnasal VEVeup'); axis square;
                            [r,~] = corr(VEVeup_Flow, VEVeup_pnasal);
                            str=['VE/Veup for Flow vs Pnasal, r = ', num2str(r)]; suptitle(str);
                            
                            VEVeup_All = [VEVeup_All; [VEVeup_Flow, VEVeup_pnasal]];
                            
                            fig = gcf;
                            fig.Color = [1 1 1]; % set background colour to white
                            fig.Units = 'inches';
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
    savestr = ['_LinRegWorkspace_', HowSoonIsNow, '.mat'];
    
    disp(displaytext);
    
    str=['Saving to: ', filename(1:end-4), savestr]; disp(str);
    save([filename(1:end-4), '_withPnasalMatched', savestr]);

end % end of L1O run

%% Most important features using DM method
% Scores are based on 1/(x+1), x=rank per loop, higher scores are best.
% Features that regularly show up early (low rank e.g. First), get the biggest scores.
% Features that are occasionally very good, but sometimes late, are still considered.
score=zeros(1,size(Amatrix2,2));
maxscore=0; %if large score is good
%maxscore = (size(labels_Step_Subj,2))^2; %if small score is good
try
    for i=1:size(labels_Step_Subj,1)
        score1 = maxscore + 0*score;
        Ftr_indx = [];
        for j=1:size(labels_Step_Subj,2)
            temp = labels_Step_Subj{i,j};
            x=sum(temp'==Ftr_indx,2);
            temp(x==1)=[];
            if ~isempty(temp)
                Ftr_indx(end+1)=temp;
            end
        end
        for j=1:length(Ftr_indx)
            score1(Ftr_indx(j)) = 1./((j-1)+1);
            %score1(If(j)) = (j-1).^0.5;
        end
        score=score+score1;
    end
    
    scoredata = [score;1:length(score)]';
    % scoredata = sortrows(scoredata,'descend'); % descend cmd only works for tble data
    scoredata = sortrows(scoredata,-1); % descending sort of col 1
    I=find(isnan(scoredata(:,1)));
    temp=scoredata(I,:);
    scoredata(I,:)=[];
    scoredata=[scoredata;temp];
    LabelsOrdered = Labels_Complete(scoredata(:,2));
catch me
    disp(me.getReport);
end

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

%% Use All breaths or sleep only breaths
SleepOnly = 1;
if SleepOnly
    BB = (PtData.Hypnog<4)&(PtData.Ar==0); % sleep only
    if IncludePnasal
        BB_ = (PtData_matched.Hypnog<4)&(PtData_matched.Ar==0); % sleep only
    end
else
    BB = true(size(predyL1O_array,1),1);
    if IncludePnasal
        BB_ = true(size(PtData_matched,1),1);
    end
end
nnz(BB)
if IncludePnasal
    nnz(BB_)
end

%% Process test data - Flow
% for all pts combined
for i=1:size(predyL1O_array,2)
    RsqL1O_(i) = 1-nansum(weights(BB).*(Gtest_All(BB)-predyL1O_array(BB,i)).^2)/nansum(weights(BB).*(Gtest_All(BB)-nanmean(Gtest_All(BB))).^2);
    ErrL1O_(i) = nanmean(weights(BB).*abs(predyL1O_array((BB),i)-Gtest_All(BB)));
    ErrL1Orms_(i) = nanmean((weights(BB).*(predyL1O_array((BB),i)-Gtest_All(BB))).^2).^0.5;
end
RL1O = RsqL1O_.^0.5;
  
% for each pt
for i=1:size(predyL1O_array,2)
    PT_list = unique(PtData.PT);
    for subj=1:54 %length(PT_list)
        if ~ismember(subj, PT_list); continue; end
        Isubj = PtData.PT==subj & BB;
        RsqL10_perPT(subj,i) = 1-nansum(weights(Isubj).*(Gtest_All(Isubj)-predyL1O_array(Isubj,i)).^2)/nansum(weights(Isubj).*(Gtest_All(Isubj)-nanmean(Gtest_All(Isubj))).^2);
        ErrL1O_perPT(subj,i) = nanmean(weights(Isubj).*abs(predyL1O_array((Isubj),i)-Gtest_All(Isubj)));
        ErrL1Orms_perPT(subj,i) = nanmean((weights(Isubj).*(predyL1O_array((Isubj),i)-Gtest_All(Isubj))).^2).^0.5;
        if (RsqL10_perPT(subj,i) < 0) % preserve sign
            RL10_perPT(subj,i) = -1*(abs(RsqL10_perPT(subj,i)).^0.5);
        else
            RL10_perPT(subj,i) = (RsqL10_perPT(subj,i)).^0.5;
        end
    end
end

%% sampler of R and Rsq values
% indiv pt Rsq values less critical
if 0
RsqL1O_sampler_individual = RsqL10_perPT(:,[3, 10, 20, 30, 40, 50]);
RsqL1O_sampler_combined = RsqL1O_([3, 10, 20, 30, 40, 50]);

RsqL1O_sampler_individual_trim= [PT_list, RsqL1O_sampler_individual(PT_list,:)];
RsqL1O_indivPt = array2table(RsqL1O_sampler_individual_trim, 'VariableNames', {'PT','Ftrs3','Ftrs10','Ftrs20','Ftrs30','Ftrs40','Ftrs50'});

RsqL1O_combinedPt = array2table(RsqL1O_sampler_combined, 'VariableNames', {'Ftrs3','Ftrs10','Ftrs20','Ftrs30','Ftrs40','Ftrs50'});
RL1O_combinedPt = array2table(RsqL1O_sampler_combined.^0.5, 'VariableNames', {'Ftrs3','Ftrs10','Ftrs20','Ftrs30','Ftrs40','Ftrs50'});
end

%% Process test data - Flow, only breaths matching with pnasal breaths
% for all pts combined
for i=1:size(predyL1O_array_flow_matched,2)
    RsqL1O_flow(i) = 1-nansum(weights_matched(BB_).*(Gtest_matched(BB_)-predyL1O_array_flow_matched((BB_),i)).^2)/nansum(weights_matched(BB_).*(Gtest_matched(BB_)-nanmean(Gtest_matched(BB_))).^2);
    ErrL1O_flow(i) = nanmean(weights_matched(BB_).*abs(predyL1O_array_flow_matched((BB_),i)-Gtest_matched(BB_)));
    ErrL1Orms_flow(i) = nanmean((weights_matched(BB_).*(predyL1O_array_flow_matched((BB_),i)-Gtest_matched(BB_))).^2).^0.5;
end
RL1O_flow = RsqL1O_flow.^0.5;

% for each pt
for i=1:size(predyL1O_array_flow_matched,2)
    PT_list = unique(PnasalFlowBB_PTindex);
    for subj=1:54 %length(PT_list)
        if ~ismember(subj, PT_list); continue; end
        Isubj = PnasalFlowBB_PTindex==subj  & BB_; 
        RsqL10_perPT_flow(subj,i) = 1-nansum(weights_matched(Isubj).*(Gtest_matched(Isubj)-predyL1O_array_flow_matched(Isubj,i)).^2)/nansum(weights_matched(Isubj).*(Gtest_matched(Isubj)-nanmean(Gtest_matched(Isubj))).^2);
        ErrL1O_perPT_flow(subj,i) = nanmean(weights_matched(Isubj).*abs(predyL1O_array_flow_matched((Isubj),i)-Gtest_matched(Isubj)));
        ErrL1Orms_perPT_flow(subj,i) = nanmean((weights_matched(Isubj).*(predyL1O_array_flow_matched((Isubj),i)-Gtest_matched(Isubj))).^2).^0.5;
        if (RsqL10_perPT_flow(subj,i) < 0) % preserve sign
            RL10_perPT_flow(subj,i) = -1*(abs(RsqL10_perPT_flow(subj,i)).^0.5);
        else
            RL10_perPT_flow(subj,i) = (RsqL10_perPT_flow(subj,i)).^0.5;
        end
    end
end

%% Process test data - Pnasal
% for all pts combined
for i=1:size(predyL1O_array_pnasal_matched,2)
    RsqL1O_pnasal(i) = 1-nansum(weights_matched(BB_).*(Gtest_matched(BB_)-predyL1O_array_pnasal_matched((BB_),i)).^2)/nansum(weights_matched(BB_).*(Gtest_matched(BB_)-nanmean(Gtest_matched(BB_))).^2);
    ErrL1O_pnasal(i) = nanmean(weights_matched(BB_).*abs(predyL1O_array_pnasal_matched((BB_),i)-Gtest_matched(BB_)));
    ErrL1Orms_pnasal(i) = nanmean((weights_matched(BB_).*(predyL1O_array_pnasal_matched((BB_),i)-Gtest_matched(BB_))).^2).^0.5;
end
RL1O_pnasal = RsqL1O_pnasal.^0.5;

% for each pt
for i=1:size(predyL1O_array_pnasal_matched,2)
    PT_list = unique(PnasalFlowBB_PTindex);
    for subj=1:54 %length(PT_list)
        if ~ismember(subj, PT_list); continue; end
        Isubj = PnasalFlowBB_PTindex==subj & BB_;
        RsqL10_perPT_pnasal(subj,i) = 1-nansum(weights_matched(Isubj).*(Gtest_matched(Isubj)-predyL1O_array_pnasal_matched(Isubj,i)).^2)/nansum(weights_matched(Isubj).*(Gtest_matched(Isubj)-nanmean(Gtest_matched(Isubj))).^2);
        ErrL1O_perPT_pnasal(subj,i) = nanmean(weights_matched(Isubj).*abs(predyL1O_array_pnasal_matched((Isubj),i)-Gtest_matched(Isubj)));
        ErrL1Orms_perPT_pnasal(subj,i) = nanmean((weights_matched(Isubj).*(predyL1O_array_pnasal_matched((Isubj),i)-Gtest_matched(Isubj))).^2).^0.5;
        if (RsqL10_perPT_pnasal(subj,i) < 0) % preserve sign
            RL10_perPT_pnasal(subj,i) = -1*(abs(RsqL10_perPT_pnasal(subj,i)).^0.5);
        else
            RL10_perPT_pnasal(subj,i) = (RsqL10_perPT_pnasal(subj,i)).^0.5;
        end
    end
end

%% Select data to use in main Figures
PlotUptoNFtrs = 50;
ftrnum=25;

DataIn = 'Flow'; % Flow, FlowP Pnasal
switch DataIn
    case 'Flow'
        R = RL1O;
        Rsq = RsqL1O_;
        Err = ErrL1O_;
        ErrRms = ErrL1Orms_;
        predy = predyL1O_array(BB,ftrnum);
        Yval = Gtest_All(BB);  
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

%% Pred VEVdrive vs Actual VEVdrive, for Flow(all), Flow(matched) and Pnasal
% plot is technically "Actual Vs Pred"
figure(20); clf(figure(20)); fig = gcf;
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

%% Classification performance
figure(21); clf(figure(21));

% performance vs number of features plot
subplot(1,3,1);
plot([R(1:PlotUptoNFtrs);Rsq(1:PlotUptoNFtrs);Err(1:PlotUptoNFtrs);ErrRms(1:PlotUptoNFtrs)]');
xlabel('Number of Ftrs'); ylim([0 1]);
legend('RL1O','RsqL1O', 'ErrL1O', 'ErrL1Orms', 'location','southeast');
title('Performance Vs Number of Ftrs');
axis square

% scatter with box overlay
subplot(1,3,2);
scatter(100*predy,100*Yval,2,'filled','markerfacealpha',0.2); hold on;
dx=0.2; xbins=[0 0.3:dx:0.9 1.5];
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
xlim([0 150]); xlabel('Predicted VEVdrive');ylabel('Actual VEVDrive');
str=['r = ', num2str(R(ftrnum))]; title(str);
axis square

% confusion mat
subplot(1,3,3);
customcmap = GetCustomColorMap('SS');
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
C_Factual=C./sum(C,2)*100 %rows are actual, cols are estimated
C_Festimated=C./sum(C)*100
C_Total = (C./sum(C(:)))*100;

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
yticks(y); yticklabels(gca,fliplr(labeltext));
xticks(x); xticklabels(gca,fliplr(labeltext));
xlabel('Flow Shape Classification');
ylabel('Actual');
title(['PredictivePerformance']);
axis square

% tidy up and save
fig = gcf;
fig.Color = [1 1 1]; 
savestr = [datastr,' AutoRef2, ', supstr];
suptitle(savestr);
% set background colour to white
fig.Units = 'inches';
fig.Position = [-12.2    3   12   4.5];
str = ['..\Figures\', savestr];
saveas(fig, str, 'png'); %savefig(str);

%% histogram of pred vs actual VEVdrive
% these plots in general look a bit ordinary, 
% is it because they do not use breath weighting, i.e. there are many more Non-FL breaths
if 1
    edges = [0:0.1:1.5];
else
    edges = xbins;
end
figure(22); clf(figure(22));
subplot(1,2,1);
histogram(Yval,edges, 'facealpha',0.2, 'edgealpha', 0.2, 'normalization','probability');
title('Actual');
xlabel('VEVdrive'); 
xlim([-0.05 1.55]);

subplot(1,2,2);
histogram(predy,edges);
title('Pred'); 
xlabel('VEVdrive'); 
xlim([-0.05 1.55]);
suptitle(['Histograms of Pred and Actual VEVdrive (All pts)',supstr]);

for subj=1:54 % set pt num or all, no data for PT=1
    if ~ismember(subj, PT_list) %PT_list)
        continue
    end
    %Isubj = PtData.PT==subj;
    Isubj=(PtData.PT==subj)&(PtData.Hypnog<4)&(PtData.Ar==0); % sleep only
    Gtest_pt = Gtest_All(Isubj);
    PredY_pt = predyL1O_array(Isubj,ftrnum);
    
    figure(300+subj); clf(figure(300+subj));fig = gcf;
    fig.Color = [1 1 1]; % set background colour to white
    fig.Units = 'inches';
    fig.Position = [-12.2   9   12    4.5];
    
    subplot(1,2,1);
    h1 = histogram(Gtest_pt,edges, 'facealpha',0.2, 'edgealpha', 0.2, 'normalization','pdf');%'probability');
    hold on;
    currentYlim = ylim();
    %plot([median(Gtest_pt), median(Gtest_pt)], [0, max(h1.Values)],'k-', 'linewidth', 2);
    plot([median(Gtest_pt), median(Gtest_pt)], [0, currentYlim(2)],'k-', 'linewidth', 2);
    title('Actual');
    xlabel('VE:Vdrive'); yticks([]);
    xlim([-0.05 1.55]);
    box off
    
    subplot(1,2,2);
    h2 = histogram(PredY_pt,edges, 'facealpha',0.2, 'edgealpha', 0.2, 'normalization','pdf');%'probability');
    hold on;
    currentYlim = ylim();
    %plot([median(PredY_pt), median(PredY_pt)], [0, max(h2.Values)],'k-', 'linewidth', 2);
    plot([median(PredY_pt), median(PredY_pt)], [0, currentYlim(2)],'k-', 'linewidth', 2);
    title('Pred');
    xlabel('VE:Vdrive'); yticks([]);
    xlim([-0.05 1.55]);
    box off
    supstr = ['Histograms of Actual and Pred VEVdrive (Flow), pt ',num2str(subj),', AHI ', num2str(round(AHI_perPT(subj,1)))];
    suptitle(supstr);
    
    % save
    str = ['..\Figures\', supstr];
    saveas(fig, str, 'png'); %savefig(str);
    
    % close
    close(fig);
end

%% scatter of PredY and Gtest at a set point
for i=ftrnum % the number of features
    PT_list = unique(PtData.PT);
    for subj=1:1%54 % set pt num or all, no data for PT=1
        if ~ismember(subj, Pnasal_list) %PT_list)
            continue
        end
        figure(100+subj); clf(figure(100+subj));
        fig = gcf;
        fig.Color = [1 1 1]; % set background colour to white
        fig.Units = 'inches';
        fig.Position = [-12.2   9   12    4.5];
        
        subplot(1,3,1); Isubj = PnasalFlowBB_PTindex==subj;
        scatter(predyL1O_array_pnasal_matched(Isubj,i), Gtest_matched(Isubj), 2, 'filled','markerfacealpha',1);
        hold on; lsline; title('Pnasal'); ylabel('Actual'); xlabel('PredY'); 
        xlim([-0.05 1.55]); ylim([-0.05 1.55]); axis square;
        
        subplot(1,3,2);
        scatter(predyL1O_array_flow_matched(Isubj,i), Gtest_matched(Isubj), 2, 'filled','markerfacealpha',1);
        hold on; lsline; title('Flow (matched)'); xlabel('PredY'); 
        xlim([-0.05 1.55]); ylim([-0.05 1.55]); axis square;
         
        subplot(1,3,3); Isubj = PtData.PT==subj;
        scatter(predyL1O_array(Isubj,i), Gtest_All(Isubj), 2, 'filled','markerfacealpha',1);
        hold on; lsline; title('Flow (All)'); xlabel('PredY'); 
        xlim([-0.05 1.55]); ylim([-0.05 1.55]); axis square;
    end
end

%% Pnasal vs Flow plots - compare pnasal and flow
Pnasal_summary_table = array2table(Pnasal_summary, 'VariableNames', {'PT','FlowBB','PnasalBB','MatchedBB','PercentMatched'});

% scatter of flow predy vs pnasal predy
figure(23);clf(figure(23));
fig = gcf;
fig.Color = [1 1 1]; 
fig.Units = 'inches';
fig.Position = [-12.2   9   12    4.5];

subplot(1, 2, 1);
%[r_t, p_t] = plotregressionwithSEM(predyL1O_array_pnasal_matched(:,ftrnum), predyL1O_array_flow_matched(:,ftrnum));
scatter(predyL1O_array_pnasal_matched(BB_,ftrnum), predyL1O_array_flow_matched(BB_,ftrnum), 2,'filled','markerfacealpha',0.3);
xlabel('Pnasal Predicted VEVdrive');ylabel('Flow Predicted VEVdrive');
str=['Flow based VS Pnasal based PredY']; title(str);

subplot(1, 2, 2);
scatter(VEVeup_pnasal_array(BB_,ftrnum), VEVeup_flow_array(BB_,ftrnum), 2,'filled','markerfacealpha',0.3);
xlabel('Pnasal VEVeup');ylabel('Flow VEVeup');
str=['Flow based VS Pnasal based VEVeup']; title(str);
xlim([0 4]); ylim([0 4]);

str_plot=['Pred performance at ', num2str(ftrnum), ' features, VEVeup is all breaths']; suptitle(str_plot);

str = ['..\Figures\Flow and Pnasal VEVdrive VEVeup'];
saveas(fig, str, 'png'); %savefig(str);


%% average Gtest vs average PredY for each pt (medians...)
figure(24); clf(figure(24)); fig = gcf;
fig.Color = [1 1 1]; 
fig.Units = 'inches';
fig.Position = [ -12.2 -2.5  12.0000    4.5000];

subplot(1, 2, 1);
Gtest_avg = NaN(54,1);
PredY_avg = NaN(54,1);
for i=ftrnum % 1:size(predyL1O_array,2) say at 20 ftrs
    PT_list = unique(PtData.PT);
    for subj=1:54 %length(PT_list)
        if ~ismember(subj, PT_list)            
            continue
        end
        %Isubj = PtData.PT==subj; % have not used weights(Isubj) ??
        Isubj=(PtData.PT==subj)&(PtData.Hypnog<4)&(PtData.Ar==0); % sleep only
        Gtest_avg(subj) = median(Gtest_All(Isubj));
        PredY_avg(subj) = median(predyL1O_array(Isubj,i));
        scatter(predyL1O_array((Isubj),i),Gtest_All(Isubj),2,'filled','markerfacealpha',0.5); hold on;
    end
end

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
str = ['..\Figures\', 'MedianGtestVsMedianPredY_AutoRef2_SleepOnly'];
saveas(fig, str, 'png'); %savefig(str);


%% median VEVdrive vs AHI
[AHI_perPT, AHI_perPT_table] = getAHI_postanalysis();
AHI_perPT_ = AHI_perPT(~isnan(AHI_perPT(:,1)),1);

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
PT_list = unique(PtData.PT);
sleeponly = 1;
addbackAP = 1;
thres = 0.5;
for subj=1:54
    if ~ismember(subj, PT_list)
        continue
    end
    if sleeponly % make Isubj
        Isubj=(PtData.PT==subj)&(PtData.Hypnog<4)&(PtData.Ar==0);
    else
        Isubj=(PtData.PT==subj);
    end
    clinscore = PtData.Etype(Isubj);
    numclinscored = nnz(clinscore==2|clinscore==4); % obstructive and hypopnoea
    predBelowThres = predyL1O_array(Isubj,ftrnum)<thres;
    numPredBBbelowthreshold = nnz(predBelowThres);
    actualBelowThres = Gtest_All(Isubj)<thres;
    numActualBBbelowthreshold = nnz(actualBelowThres);
    numBBtotal = size(predyL1O_array(Isubj),1);
    if addbackAP
        withAPstr = [' (incl Ap-O breaths)'];
        indAP = find(RemovedBB_Apnoea.Pt == subj);
        indLF = find(RemovedBB_LowFlow.Pt == subj);
        if ~isempty(indAP)
            numclinscored = numclinscored+sum(RemovedBB_Apnoea{indAP,[4 6 7]});
            numPredBBbelowthreshold = numPredBBbelowthreshold+RemovedBB_Apnoea{indAP,2};
            numActualBBbelowthreshold = numActualBBbelowthreshold+RemovedBB_Apnoea{indAP,2};
            numBBtotal = numBBtotal+RemovedBB_Apnoea{indAP,2};
        end
        if ~isempty(indLF)
            numclinscored = numclinscored+sum(RemovedBB_LowFlow{indLF,[4 6 7]});
            numPredBBbelowthreshold = numPredBBbelowthreshold+RemovedBB_LowFlow{indLF,2};
            numActualBBbelowthreshold = numActualBBbelowthreshold+RemovedBB_LowFlow{indLF,2};
            numBBtotal = numBBtotal+RemovedBB_LowFlow{indLF,2};
        end
    else
        withAPstr = [''];
    end
    Clinpercent = 100*(numclinscored/numBBtotal);
    Predpercent = 100*(numPredBBbelowthreshold/numBBtotal);
    Actualpercent = 100*(numActualBBbelowthreshold/numBBtotal);
    ptstats = [ptstats; [subj,numBBtotal,numclinscored,numPredBBbelowthreshold,numActualBBbelowthreshold,Clinpercent,Predpercent,Actualpercent]];
end
ptsummary_withAp = array2table(ptstats, 'VariableNames', {'PT','TotalBB','Clin_FL','Pred_FL','Actual_FL','Clin_percent','Pred_percent','Actual_percent'});


%% Proportion FL VS AHI
figure(26);clf(figure(26)); fig = gcf;
fig.Color = [1 1 1]; 
fig.Units = 'inches';
fig.Position = [-12.2   2.5   12    4.5];

subplot(1,2,1)
%[r_t, p_t] = plotregressionwithSEM(ptsummary_withAp_n.Pred_percent, AHI_perPT(:,1));
scatter(ptsummary_withAp.Pred_percent, AHI_perPT_(:,1),50,'filled','markerfacealpha',1);
hndl_ls = lsline;  hndl_ls.LineStyle = '-'; hndl_ls.Color = [0 0 1];
hndl_ref = refline(1,0); hndl_ref.LineStyle = '--'; hndl_ref.Color = [0.5 0.5 0.5];
xlim([-5 105]); xlabel('%FL (NumFL BB/TotalBB)');
ylim([-5 105]); ylabel('AHI');
title(['AHI Vs Pred %FL']);

subplot(1,2,2)
scatter(ptsummary_withAp.Actual_percent, AHI_perPT_(:,1),50,'filled','markerfacealpha',1); 
hndl_ls = lsline;  hndl_ls.LineStyle = '-'; hndl_ls.Color = [0 0 1];
hndl_ref = refline(1,0); hndl_ref.LineStyle = '--'; hndl_ref.Color = [0.5 0.5 0.5];
xlim([-5 105]); xlabel('%FL (NumFL BB/TotalBB)');
ylim([-5 105]); ylabel('AHI');
title(['AHI Vs Actual %FL']);

suptitlestr=['AHI Vs Proportion of sleep breaths classified as FL', withAPstr];
suptitle(suptitlestr);


%% Novel metrics, Median VE:Vdrive during sleep and Time with severe obstruction during sleep
% aka, Proportion of breaths FL
figure(27); clf(figure(27)); fig = gcf;
fig.Color = [1 1 1]; 
fig.Units = 'inches';
fig.Position = [-12.2   -3   12    4.5];

subplot(1,2,1);
[r_1, p_1] = plotregressionwithSEM(PredY_avg, Gtest_avg);
%scatter(PredY_avg, Gtest_avg, 50, 'filled'); hold on; 
%hndl_ls = lsline;  hndl_ls.LineStyle = '-'; hndl_ls.Color = [0 0 1];
%hndl_ref = refline(1,0); hndl_ref.LineStyle = '--'; hndl_ref.Color = [0.5 0.5 0.5];
xlim([0 1.1]); xlabel('PredY VE:VDrive');
ylim([0 1.1]); ylabel('Actual VE:Vdrive');
axis square
titlestr = ['Patient median VEVdrive, r=', num2str(r_1)]; title(titlestr);

subplot(1, 2, 2);
[r_2, p_2] = plotregressionwithSEM(ptsummary_withAp.Pred_percent, ptsummary_withAp.Actual_percent);
%scatter(ptsummary_withAp.Pred_percent, ptsummary_withAp.Actual_percent,50,'filled','markerfacealpha',1); 
%hndl_ls = lsline;  hndl_ls.LineStyle = '-'; hndl_ls.Color = [0 0 1];
%hndl_ref = refline(1,0); hndl_ref.LineStyle = '--'; hndl_ref.Color = [0.5 0.5 0.5];
xlim([-5 105]); xlabel('PredY %FL (NumFL BB/TotalBB)');
ylim([-5 105]); ylabel('Actual %FL (NumFL BB/TotalBB)');
titlestr=['Proportion FL (VE:Vdrive < ',num2str(thres), '), r= ', num2str(r_2)]; title(titlestr);

suptitle(['Novel metrics (during sleep)']);
str = ['..\Figures\', 'NovelMetrics_AutoRef2_SleepOnly'];
saveas(fig, str, 'png'); %savefig(str);


%% AHI vs Novel Metrics
% grab from plot 25 and 26
figure(32);clf(figure(32)); fig = gcf;
fig.Color = [1 1 1]; 
fig.Units = 'inches';
%fig.Position = [ -12.2 8  12    4.5];
fig.Position = [ 2 3  12    4.5];

subplot(1,2,1)
%[r_1, p_1]=plotregressionwithSEM(PredY_avg, AHI_perPT(:,1));
r_1 = corr(PredY_avg_, AHI_perPT_(:,1));
scatter(PredY_avg_, AHI_perPT_(:,1), 50, 'filled'); hold on; 
xlim([-0.05 1.55]); ylim([-5 105]); xlabel('Pred VE:VDrive');
ylabel('AHI'); title('AHI Vs Pred VEVdrive')
axis square
subplot(1,2,2)
%[r_2, p_2]=plotregressionwithSEM(ptsummary_withAp.Pred_percent, AHI_perPT_(:,1));
r_2 = corr(ptsummary_withAp.Pred_percent, AHI_perPT_(:,1));
scatter(ptsummary_withAp.Pred_percent, AHI_perPT_(:,1),50,'filled','markerfacealpha',1); 
%hndl_ls = lsline;  hndl_ls.LineStyle = '-'; hndl_ls.Color = [0 0 1];
%hndl_ref = refline(1,0); hndl_ref.LineStyle = '--'; hndl_ref.Color = [0.5 0.5 0.5];
xlim([-5 105]); xlabel('%FL (NumFL BB/TotalBB)');
ylim([-5 105]); ylabel('AHI');
title(['AHI Vs Pred %FL']);

suptitlestr=['AHI Vs Novel Metrics'];
suptitle(suptitlestr);
str = ['..\Figures\', 'AHIvsNovelMetrics_AutoRef2_SleepOnly'];
saveas(fig, str, 'png'); %savefig(str);

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
if (nnz(predy<0)>0); predy(predy<0)=0; end
if nnz(predy>maxG)>0; predy(predy>maxG)=1.5; end

% read spreadsheet
AnalyzeDataSpreadsheet = 'AnalyzeDataSpreadsheet.xlsx';
[num,patients,~] = xlsread(AnalyzeDataSpreadsheet,1,'F3:G56');
% list of those with multiple regions of interest [ 14 29 34 ]
% full list [ 2 6 14 17 18 21 28 29 34 44 45 53 ]
for n=[34 ]
    if ~ismember(n, PT_list)
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
Data1 = [PtData.BB_time(PtData.PT==n) PtData.BB_Ttot(PtData.PT==n) predy(PtData.PT==n) Gtest_All(PtData.PT==n)];
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
figure(200+n); clf(200+n); fig = gcf;
fig.Color = [1 1 1]; 
ax(1)=subplot(3,1,1);
stairs(Data1(:,1),100*Data1(:,4),'k'); % actual
hold('on')
stairs(Data1(:,1),100*Data1(:,3),'r'); % pred
plot([Data1(1,1) Data1(end,1)],[0 0],'k:');
plot([Data1(1,1) Data1(end,1)],100*[1 1],'k:');
plot([Data1(1,1) Data1(end,1)],50*[1 1],'k:');
set(gca,'xtick',[],'box','off');
ylabel('V_E:V_{drive} (%wake)'); ylim([-5 105]);

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

% mods to plot
switch n
    case 2
        v1 = 27000;
        v2 = 27350;
        useStart = 1;
    case 6
        v1 = 2900;
        v2 = 3600;
        useStart = 1;
    case 14
        % multiple
        S14 = StarttimeSpike;
        if 0
            v1 = 5200;
            v2 = 6400;
        else
            StarttimeSpike = S14; 
            v1 = 23560;
            v2 = 23680;
        end
        useStart = 1;
    case 17
        v1 = 83450;
        v2 = 83700;
        useStart = 0;
    case 18
        v1 = 105600;
        v2 = 106000;
        useStart = 0;
    case 21
        v1 = 83900;
        v2 = 85000; % not defined
        useStart = 0;
    case 28
        v1 = 19460;
        v2 = 19580;
        useStart = 1;
    case 29
        % multiple
        S29 = StarttimeSpike;
        if 0
            v1 = 13400;
            v2 = 13600;
        else
           StarttimeSpike = S29; 
           v1 = 15800;
           v2 = 16000;
        end
        useStart = 1;
    case 34 
        % multiple
        if 0
            v1 = 90410;
            v2 = 90560;
        elseif 0
            v1 = 100000;
            v2 = 100400;
        else
            v1 = 101705;
            v2 = 102100;
        end
        useStart = 0;
    case 44
        v1 = 96500;
        v2 = 97600;
        useStart = 0;
    case 45
        v1 = 86450;
        v2 = 87000;
        useStart = 0;
    case 53
        v1 = 103200;
        v2 = 104200;
        useStart = 0;
end
if useStart
    xlim([StarttimeSpike+v1 StarttimeSpike+v2]) 
else
    xlim([v1 v2]);
end


str = ['..\Figures\SampleData_Pt ', num2str(n), ', snip 3'];
saveas(fig, str, 'png'); %savefig(str);

close(fig);

end











%%
xlim([StarttimeSpike, inf]);










%% Manually pick optimal number of features
NfeaturesOpt=25;
LabelsOrderedOpt = LabelsOrdered(1:NfeaturesOpt);

%% Final model using all data, uses selected N optimal features (NfeaturesOpt)
%not tested
if 0
    If = 1:size(Amatrix2,2);
    Labels = Labels_;
    MaxNfeatures=NfeaturesOpt;
    while length(If)>=MaxNfeatures
        
        [Rsq,Pvals,RMSE,b]=glmfitFast(Amatrix2(:,If),Gtest_,weights,1); %faster
        Pvals(1)=[];
        
        [maxp,maxpi]=max(Pvals);
        if length(If)==MaxNfeatures
            break
        end
        disp(['Removing Ftr: ', num2str(If(maxpi)), ', p= ' num2str(maxp)]);
        If(maxpi)=[];
        Labels(maxpi)=[];
    end
    Labels = Labels_(If)
    a=sortrows( [Pvals (1:NfeaturesOpt)'])
    Labels(a(:,2))
    direction = b(2:end)>0;
end

if 0
    tic
    finalmdl = fitglm(Amatrix2(:,If(1:50)),Gtest_,'weights',weights);
    toc
    tic
    mdl_prediction = predict(finalmdl, Amatrix2(:,If(1:50)));
    toc
    % test mdl_prediction against Gtest_
    figure(4); clf(figure(4));
    scatter(Gtest_, mdl_prediction,2,'filled','markerfacealpha',0.4);
    xlabel('Gtest'); ylabel('Mdl Pred');
end



%% moving time median
if 0
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


%%
AAsummary = [alpha MaxNfeatures median(Nfeatures) Fwd useweights median(RvalueTrain) Rvalue ...
    ACCs(2,:) median(Err) median(ErrRms) ErrL1O ErrL1Orms];



%% Patient level summary data
%Need to put back breaths that are thrown out because of low VE here (apneas)
for i=1:54%max(PtData.PT)
    I=PtData.PT==i&PtData.NotAr==1; nnz(I)
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


%% Correlate with AHI
figure(6); clf(6);
subplot(1,1,1);
[Rtemp,Ptemp]=plotregressionwithSEM(AHItotal',medianGest); title(num2str(Rtemp));
AHItotal = AHI_perPT(~isnan(AHI_perPT(:,1)),1);

[b,dev,stats]=glmfit([AHItotal medianGest'],medianG')
[b,dev,stats]=glmfit([AHItotal Fsevest'],Fsev')


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
If=[1:20];
rangei=1:8000;
mdlu=fitglm(Amatrix2(rangei,If),Gtest_(rangei),'Distribution','binomial','weights',weights(rangei));
mdlr=fitglm(Amatrix2(rangei,If(1:end-1)),Gtest_(rangei),'Distribution','binomial','weights',weights(rangei));%
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
If=[1:200];
tic
[b,dev,stats]=glmfit(Amatrix2_(rangei,If),Gtest2_(rangei),'binomial','weights',weights_(rangei));
toc
Pvals = stats.p(2:end);

optionsX=statset('glmfit');
optionsX.TolX=1e-02;
tic
[b,dev,stats]=glmfitFastLR(Amatrix2_(rangei,If),Gtest2_(rangei),'binomial','weights',weights_(rangei),'options',optionsX);
toc
Pvals = stats.p(2:end);

