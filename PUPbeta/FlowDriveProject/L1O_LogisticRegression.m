%% start
keyboard
close all
clear
clc
addpath('C:\Users\uqdmann\Dropbox\PUPbeta_git\');
cd('C:\Users\uqdmann\Dropbox\PUPbeta_git\FlowAndEdiPesQAO');

%% options
Use25HzDS = 1;          % set as 1 for 25Hz data, or 0 to use 125Hz data
FlowDriveData = 0;      % set as 1 for flow drive data, or 0 to use OA data
EdiDrive = 1;           % set as 1 for Edi Drive data, 0 for Pes drive
LinearReg = 0;			% set as 1 for linear regression, 0 for logistic regression
SS_fastfit = 1;			% set as 1 to use SS fast glmfit method, 0 for matlab builtin
SS_logitreg = 0;        % set as 1 to run the SS version, or 0 for DM alternate version
% the only difference is an explicit listing of
% ",'link','logit'" in the glmfit methods, this is
% the default in 'binomial', so theory says no change.
TransformTheData = 0;   % set as 1 to do tranforms, or 0 to use unadjusted data
addextratransform = 0; 	% set as 1 to do extra transforms
L1O_run = 1;            % set as 1 to run to do the L1O run, or 0 to go
% straight to processing results (loads saved L1O run)
ShowFigures = 0;        % set as 1 to show figures, or 0 to not show figures

%datadir = 'C:\PSG_Data\FlowDrive\FeatureSpaces\';
%datadir = '..\FeatureSpaces\';  %
datadir = 'C:\PSG_Data\OralAppliance\Analyzed\';

%%
if L1O_run
    %% Turn diary logging on
    HowSoonIsNow = datestr(datetime('now'),'yyyymmddHHMM');
    diaryfilename = ['L1O_LogReg_OA_DLM_', HowSoonIsNow, '.txt'];
    diary(diaryfilename);
    diary on
    
    %% open file
    if FlowDriveData
        if Use25HzDS
            if EdiDrive
                filename = [datadir, 'FlowDrive_only25Hz_Normalized_FeatureSpace_AutoRef2_Edi_Clean.mat']; datastr = 'FD Normalized edi 25Hz';
                %filename = [datadir, 'FlowDrive_only25Hz_FeatureSpace_Edi_Clean.mat']; datastr = 'FD edi 25Hz';
            else; filename = [datadir, 'FlowDrive_only25Hz_FeatureSpace_Pes_Clean.mat']; datastr = 'FD pes 25Hz'; end
        else
            if EdiDrive; filename = [datadir, 'FlowDrive_only125Hz_FeatureSpace_Edi_Clean.mat']; datastr = 'FD edi 125Hz';
            else; filename = [datadir, 'FlowDrive_only125Hz_FeatureSpace_Pes_Clean.mat']; datastr = 'FD pes 125Hz';end
        end
    else
        if Use25HzDS; filename = [datadir, 'FlowOA_25Hz_ReNormalized_FeatureSpace_NoRef_VEVeup_Clean.mat']; datastr = 'OA 25Hz';
        else; filename = [datadir, 'FlowOA_only125Hz_FeatureSpace_VEVeup_Clean.mat']; datastr = 'OA 125Hz'; end
    end
    
    str=['Loading ' filename]; disp(str);

    try
        load(filename);
    catch me
        disp(me.getReport);
    end
    
    %load([datadir, 'SA_FD_16.mat']);
    %load([datadir, 'SSblend_FS_FD_56.mat']);
    %clearvars -except PtData FeatureNames Amatrix %% tidy up, if opening one of the old files
    
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
        str = ['Zero NaN-breaths were found']; disp(str);
    end
    
    %% use only a given set of features, determined from previous run
    if 0
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
    end
    
    %% Remove features in "FtrsToExclude" list, using string search
    if 1
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
    
    %% set the data matrix to use. either (0) unadjusted or (1) tranformed
    if TransformTheData
        [Amatrix2, Labels] = DoDataMatTransform(Amatrix, FeatureNames, addextratransform);
    else
        Amatrix2 = [Amatrix];
        Labels = FeatureNames.Name;
    end
    
    %% Preparation (Do this carefully)
    if FlowDriveData
        %Gtest_ = PtData{:,13}; % continuous
        Gtest_ = PtData.g_Edi_Adj;
        % make a handy table of the data
        Gt = table(Gtest_, 'VariableNames', {'G'});
        Atable = array2table(Amatrix, 'VariableNames', FeatureNames{:,2});
        data_in = [Gt Atable];
        if ~LinearReg % if logisitic reg, convert the continuous into binary
            GtestLR=NaN*Gtest_;
            GtestLR(Gtest_>0.7)=1; % NFL  nnz(GtestLR==1)
            GtestLR(Gtest_<=0.7)=0; % FL  nnz(GtestLR==0)
            GtestLR(isnan(Gtest_))=NaN;
            Gtest_=GtestLR;
        end
    else % OralAppliance data
        Gtest_ = ones(size(PtData,1),1)*2; % set a vector of two's,
        % that we will mark FL as zero, nonFL as one, and remove the twos
        
        % NonFL breaths are:
        %  Arousal, by at least 2 breaths (so using A1, which is threshold, more like about 4 breaths)
        %  VE > Veup*0.9, VE is good
        %  Etype == 0, not a clinically scored event
        Gtest_(PtData.A1==1 & PtData.VE>(PtData.Veup*0.9) & PtData.Etype==0)=1;
        
        % FL breaths are:
        %  NotAr
        %  VE < Veup*0.7,
        %  Etype 2 = ApO, 4 = HypO.Could extend to include 5 = M.
        % use the clinically scored events (hypop/obstructive) as true labels
        Gtest_(PtData.NotAr==1 & PtData.VE<(PtData.Veup*0.7) & (PtData.Etype==2 | PtData.Etype==4))=0;
        exclude = Gtest_==2;
        str = [num2str(length(exclude)), ' breaths available for analysis']; disp(str)
        str = [num2str(nnz(Gtest_==1)), ' breaths were labelled as NonFL, and ', num2str(nnz(Gtest_==0)), ' as FL']; disp(str)
        str = [num2str(nnz(exclude)), ' breaths were unlabelled and removed']; disp(str)
        
        % What is the breakdown of the Gtest_==2 breaths, i.e. neither FL or NonFL
        VEVEup = PtData.VE./PtData.Veup;
        str = [num2str(nnz(VEVEup>0.7 & VEVEup<0.9)),...
            ' of the removed breaths had 70% > VEVeup > 90%']; disp(str)
        Ar1 = PtData.A1==0; % find arousal state
        str = [num2str(nnz(VEVEup(Gtest_==2)>0.9&Ar1(Gtest_==2))),...
            ' of the removed breaths had VEVeup >90%, but were not far enough away from sleep']; disp(str)
        
        % Ar    - Arousal breath
        % NotAr - Ar=1 is always NotAr=0, Ar=0 is NotAr=1 but only after delay of two breaths
        % A2    - Arousal breath, at least 2 breaths away from sleep
        % A1    - Arousal breath, at least 'threshold' breaths away from sleep (this average to about 4 breaths)
        
        % remove breaths that are neither FL, nor NonFL
        Gtest_(exclude)=[];
        PtData(exclude,:)=[];
        Amatrix(exclude,:)=[];
        Amatrix2(exclude,:)=[];
    end
    colofones = ones(length(Gtest_),1);
    predyL1O = NaN*Gtest_;
    NfeatureSteps=100;
    predyL1O_array = NaN*ones(length(Gtest_),NfeatureSteps);
    
    %% Weights
    % ToDo: check method of weights
    
    clear Ndata
    if LinearReg
        maxG=1.5;
        dx=0.2;
        xbins=[0 0.3:dx:0.9 maxG];
        %xbins=[0 0.1:0.05:1.1 1.5];
    else
        xbins=[-Inf 0.5 +Inf];  % should this be 0.7
    end
    
    for i=1:length(xbins)-1
        Ix=Gtest_>xbins(i)&Gtest_<=xbins(i+1);
        Ndata(i)=sum(Ix);
    end
    %Ndata = Ndata.^0.5;
    weightsbins = 1./(Ndata);
    %weightsbins = [2 1 1 1 0.5];
    weightsbins = weightsbins/mean(weightsbins);
    %weightsbins = weightsbins/weightsbins(end);
    
    weights = NaN*Gtest_;
    for i=1:length(xbins)-1
        Ix=Gtest_>=xbins(i)&Gtest_<=xbins(i+1);
        weights(Ix)=weightsbins(i);
    end
    weights = weights/nanmean(weights); % nnz(~isfinite(weights))
    
    useweights=1;
    if ~useweights %overwrite, use no weights
        weights = ones(length(weights),1);
    end
    
    %% Leave one subject out loop (moderately slow)
    % log reg leave one out
    %          30 subjects, 80k breaths, takes ~8 hrs
    % DLM. FD, 41 Subjects, 136k breaths, 153 features, no txforms, takes ~10.5 hrs
    % DLM. OA, 29 Subjects, 16k breaths, 153 features, no txforms, takes ~1 hr
    %      OA data contains ~100k breaths, many removed in classification
    % lin reg leave one out
    %      FD, 41 Subjects, 136k breaths, 50 ftrs, 2x transforms, takes 30 min to 1 hr
    
    % use glmfit to do generalized linear model regression
    % use fitglm to create a generalized linear regression model
    try
        clear RvalueTrain Err ErrRms Nfeatures
        clear RsqTrain_array ErrTrain_array ErrRmsTrain_array
        labels_Step_Subj =[];
        alpha=0.05;
        MaxNfeatures=1;
        neverbreak=1;
        t_start_L1O = clock;
        
        Labels_ = Labels;
        Fwd=0; %backwards is crashing, perhaps duplicate data included
        
        pts = unique(PtData.PT);
        for subj=1:length(pts)      
            %% loop through all subjects, leaving one out each time
%             if ~ismember(subj, pts)    
%                 str=['No data for Pt ', num2str(subj)]; disp(str); continue
%             end
            
            disp(' '); str=['Performing analysis, witholding Pt ', num2str(pts(subj))]; disp(str); tic;                     
            
            Labels = Labels_;
            Isubj=PtData.PT==pts(subj);
            Gtest = Gtest_;
            if LinearReg
                Gtest(Gtest>maxG)=maxG;
            end
            Gtest(Isubj)=[];
            w = weights;
            w(Isubj)=[];
            w=w/nanmean(w);
            colofones_ = colofones;
            colofones_(Isubj)=[];
            
            Amatrix2_ = Amatrix2;
            Amatrix2_(Isubj,:)=[];
            
            %remove all NAN from all rows
            nanx = sum(isnan(Amatrix2_),2)>0;   % nnz(nanx)
            nany = isnan(Gtest);                % nnz(nany)
            Amatrix2_(nanx|nany,:)=[];
            Gtest(nanx|nany)=[];
            w(nanx|nany)=[];
            w=w/nanmean(w);
            colofones_(nanx|nany)=[];
            warning('off')
            maxp=Inf;
            if Fwd==0 %backwards elimination
                If = 1:size(Amatrix2_,2);
                while length(If)>=MaxNfeatures || maxp>alpha
                    if LinearReg
                        if SS_fastfit % nnz(~isfinite(Amatrix2_))
                            [Rsq,Pvals,RMSE,b]=glmfitFast(Amatrix2_(:,If),Gtest,w,1); %faster
                            Pvals(1)=[];
                        else % matlab built-in
                            [b,dev,stats]=glmfit(Amatrix2_(:,If),Gtest,'normal','weights',w);
                            Pvals = stats.p(2:end);
                        end
                    else %logistic reg, no fast option
                        if SS_logitreg
                            [b,dev,stats]=glmfit(Amatrix2_(:,If),Gtest,'binomial','weights',w);
                        else
                            if 1
                                [b,dev,stats]=glmfit(Amatrix2_(:,If),Gtest,'binomial','weights',w ,'link','logit');
                            else % testing LDA
                                LDA_Mdl = fitcdiscr(Amatrix2_(:,If),Gtest,'Weights',w);
                                LDA_out = predict(LDA_Mdl, Amatrix2_(:,If));
                                plotconfusion(Gtest',LDA_out');
                                [b,dev,stats]=glmfit(Amatrix2_(:,If),Gtest,'binomial','weights',w ,'link','probit');
                                [b, bint,r,rint,stats] = regress(Gtest,Amatrix2_(:,If));
                            end
                        end
                        Pvals = stats.p(2:end);
                    end
                    
                    if length(If)<=NfeatureSteps %start saving results
                        if LinearReg
                            predyL1O_array(Isubj,length(If)) = [colofones(Isubj) Amatrix2(Isubj,If)]*b;
                            predytrain = [colofones_ Amatrix2_(:,If)]*b;
                            predytrain(predytrain>maxG)=maxG;
                            predytrain(predytrain<0)=0; %will be overwritten
                        else
                            temp = [colofones(Isubj) Amatrix2(Isubj,If)]*b;
                            predyL1O_array(Isubj,length(If)) = 1./(1+exp(-temp));
                            temp = [colofones_ Amatrix2_(:,If)]*b;
                            predytrain = 1./(1+exp(-temp));
                        end
                        RsqTrain_array(subj,length(If)) = 1-nansum(w.*(Gtest-predytrain).^2)/nansum(w.*(Gtest-nanmean(Gtest)).^2);
                        ErrTrain_array(subj,length(If)) = nanmean(w.*abs(predytrain-Gtest));
                        ErrRmsTrain_array(subj,length(If)) = nanmean((w.*(predytrain-Gtest)).^2).^0.5;
                        labels_Step_Subj{subj,length(If)} = If;
                    end
                    %remove least important
                    if 0&&length(If)>200 %speed up at higher levels
                        temp=[Pvals';1:length(Pvals)]';
                        temp=sortrows(temp);
                        remi = temp(end-20+1:end,2);
                        maxp = temp(end-20+1,1);
                    else
                        [maxp,remi]=max(Pvals);
                        
                    end
                    if length(If)==1 || (~neverbreak && maxp<alpha && (length(If)-1)<=MaxNfeatures)
                        break
                    end
                    disp(['Removing Ftr: ', num2str(If(remi)), ', p= ' num2str(maxp)]);
                    If(remi)=[];
                    Labels(remi)=[];
                end
            else %forwards stepwise inclusion
                if 1 %Matlab stock
                    [b,se,Pvals,inmodel,stats,nextstep,history] = stepwisefit(...
                        [Amatrix2_],Gtest,'penter',alpha,'premove',2*alpha,'MaxIter',MaxNfeatures); %no weighting possible here
                    If = [];
                    for i=1:sum(inmodel)
                        temp = find(history.in(i,:));
                        x=sum(temp'==If,2);
                        temp(x==1)=[];
                        if ~isempty(temp)
                            If(end+1)=temp;
                        end
                    end
                else %robust
                    If = [];
                    notinc = 1:size(Amatrix2_,2);
                    while length(If)<=MaxNfeatures %&& maxp<alpha
                        RMSE=NaN*notinc;
                        P=NaN*notinc;
                        tic
                        for i=1:length(notinc)
                            %[~,Pvals,RMSE(i),~]=glmfitFast(Amatrix2_(:,[If notinc(i)]),Gtest,w);
                            %P(i) = Pvals(end);
                            RMSE(i)=glmfitFastRMSE(Amatrix2_(:,[If notinc(i)]),Gtest,w);
                        end
                        toc
                        [~,mini]=min(RMSE);
                        [~,Pvals,~,~]=glmfitFast(Amatrix2_(:,[If notinc(mini)]),Gtest,w,0);
                        if Pvals(end)>alpha
                            disp('end FSR');
                            break
                        end
                        If = [If notinc(mini)];
                        notinc(mini)=[];
                        Labels_(If)
                        %check for remove
                        [~,Pvals,RMSE,~]=glmfitFast(Amatrix2_(:,If),Gtest,w,0);
                        Pvals(1)=[];
                        while max(Pvals)>2*alpha
                            [~,remi]=max(Pvals);
                            notinc = [notinc If(remi)];
                            If(remi)=[];
                            [~,Pvals,RMSE,~]=glmfitFast(Amatrix2_(:,If),Gtest,w,0);
                            Pvals(1)=[];
                        end
                    end
                end
                [Rsq,Pvals,RMSE,b]=glmfitFast(Amatrix2_(:,If),Gtest,w,0);
                Pvals(1)=[];
                Labels = Labels_(If); % Check this order?
            end
            
            Bvals = b(2:end);
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
            
            if LinearReg
                predytrain = [colofones_ Amatrix2_(:,If)]*b;
                predytrain(predytrain>maxG)=maxG;
                predytrain(predytrain<0)=0;
                predyL1O(Isubj) = [colofones(Isubj) Amatrix2(Isubj,If)]*b;
            else
                temp = [colofones_ Amatrix2_(:,If)]*b;
                predytrain = 1./(1+exp(-temp));
                temp = [colofones(Isubj) Amatrix2(Isubj,If)]*b;
                predyL1O(Isubj) = 1./(1+exp(-temp));
                
            end
            
            Rsq = 1-nansum(w.*(Gtest-predytrain).^2)/nansum(w.*(Gtest-nanmean(Gtest)).^2);
            RvalueTrain(subj) = Rsq^0.5;
            Err(subj) = nanmean(abs(predyL1O(Isubj)-Gtest_(Isubj)));
            ErrRms(subj) = nanmean((predyL1O(Isubj)-Gtest_(Isubj)).^2).^0.5;
            Nfeatures(subj) = length(Labels2);
            toc
            %keyboard
        end
    catch me
        disp(me.getReport);
    end
    
    % display processing time
    delta_t = etime(clock, t_start_L1O); % delta in seconds
    D = duration(0,0,round(delta_t),'Format','hh:mm:ss'); % duration in HH:MM:SS
    disp(' '); % add row space for visual clarity in command window
    if LinearReg
        displaytext = ['L1O linear regresion process complete. Total time: ', char(D), ' (hh:mm:ss)'];
        savestr = ['_LinRegWorkspace.mat'];
    else
        displaytext = ['L1O logistic regresion process complete. Total time: ', char(D), ' (hh:mm:ss)'];
        savestr = ['_LogRegWorkspace.mat'];
    end
    disp(displaytext);
    
    str=['Saving to: ', filename(1:end-4), savestr]; disp(str);
    save([filename(1:end-4), savestr]);
    
    diary off
    
    keyboard
    
end % end of L1O run

%%
% >
% >
% >
% >
%

%% load saved data
%str=['Loading ' [filename(1:end-4), '_LogRegWorkspace.mat']]; disp(str);
close all
clear
clc
%datadir = 'C:\PSG_Data\FlowDrive\FeatureSpaces\';
datadir = 'C:\PSG_Data\OralAppliance\Analyzed\';
load([datadir, 'FlowOA_25Hz_ReNormalized_FeatureSpace_NoRef_VEVeup_Clean_LogRegWorkspace.mat']);

savefigs = 0;
closefigs = 0;

ftrnum=50;

TickFntSz = 12;
LabelFntSz = 18; 
FntSz = 18; 

savefigas = '';





%% Histogram of ftrs selected during training
ftr_array = NaN(length(pts), ftrnum); i = 0;

for subj=1:length(pts) %32 % set pt num or all, no data for PT=1
%     if ~ismember(subj, pts) %PT_list)
%         continue
%     end
    ftr_array(subj,:) = labels_Step_Subj{subj,ftrnum};
end

%ftr_array = ftr_array(pts,:);
ftr_array_linear = ftr_array(:);
ftrsfoundat50 = unique(ftr_array_linear);
ctrs = 1:1:165;
[counts, ~] = hist(ftr_array_linear, ctrs);
[~,I_uw] = sort(counts, 'descend');

figure(33); clf(figure(33)); fig = gcf;
fig.Color = [1 1 1]; fig.Units = 'inches';
fig.Position = [2   5   12    4.5];
bar(counts(I_uw), 'facealpha',0.2, 'edgealpha', 0.2);
xlim([0 50]);
ax = gca;
ax.TickDir = 'out';
ax.FontSize=9;
ax.XTick=[1:1:50];
xtickangle(ax,90);

LabelsOrdered_uw = Labels_(I_uw);
lbls_uw = LabelsOrdered_uw(1:100);
lbls = regexprep(lbls_uw,'[_,:{}]','');
xticklabels(ax, lbls(1:50));

ylabel('Frequency');
box off
title(['Histogram of Features (simple count, 50 ftr model)']);
str = ['..\Figures\', 'HistogramOfFeatures_SimpleCount']; %
if savefigs; saveas(fig, str, 'png'); end %savefig(str);
if closefigs; close(fig); end



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
    LabelsOrdered_w = Labels_(scoredata(:,2));
catch me
    disp(me.getReport);
end

% Modified histogram of ftrs selected during training
figure(34); clf(figure(34)); fig = gcf;
fig.Color = [1 1 1]; fig.Units = 'inches';
fig.Position = [2.2   6   12    4.5];
bar(scoredata(:,1), 'facealpha',0.2, 'edgealpha', 0.2);
xlim([0 50]);
ax = gca;
ax.TickDir = 'out';
ax.FontSize=9;
ax.XTick=[1:1:50];
xtickangle(ax,90);

lbls_w =  LabelsOrdered_w(1:100);
lbls = regexprep(lbls_w,'[_,:{}]','');
xticklabels(ax,lbls(1:50));

ylabel('Weighted Score');
box off
title(['Modified Histogram of Features (weighted count, 50 ftr model)']);
str = ['..\Figures\', 'HistogramOfLogRegFeatures_WeightedCount']; %
if savefigs; saveas(fig, str, 'png'); end %savefig(str);
if closefigs; close(fig); end


%% compare the two ways of "top 50" ftrs
x_over = ismember(lbls_w(1:ftrnum), lbls_uw(1:ftrnum));
nnz(x_over)

%% save most important fetaures
if 0 % don't do this everytime...
LabelsOrderedOpt = LabelsOrdered_w(1:50);
save('FtrsToUse_ReNormalized_2', 'LabelsOrderedOpt');
end

%% Find Top X unique features based on scores
% if using transformed data, need to strip off the suffix
if 0
if TransformTheData; lengthsuf=4; else; lengthsuf=0; end
MaxUniquefeatures=50;
LabelsOrdered_NoSuffix=LabelsOrdered;
for i=1:length(LabelsOrdered_NoSuffix)
    LabelsOrdered_NoSuffix{i}=LabelsOrdered_NoSuffix{i}(1:length(LabelsOrdered_NoSuffix{i})-lengthsuf);
end
LabelsOrdered_NoSuffixUnique=unique(LabelsOrdered_NoSuffix,'stable');
LabelsOrdered_NoSuffixUnique(MaxUniquefeatures+1:end)=[];
end


%% What data to use and plot

for i=1:size(predyL1O_array,2)
    RsqL1O_(i) = 1-nansum(weights.*(Gtest_-predyL1O_array(:,i)).^2)/nansum(weights.*(Gtest_-nanmean(Gtest_)).^2);
    ErrL1O_(i) = nanmean(weights.*abs(predyL1O_array(:,i)-Gtest_));
    ErrL1Orms_(i) = nanmean((weights.*(predyL1O_array(:,i)-Gtest_)).^2).^0.5;
end
RL1O = RsqL1O_.^0.5;

if 0 %training, last loop data
    %Gtest = Gtest;
    %predy = predyL1O;
    predy = predytrain;% predyL1O_fsmr;
    %Yval = Gtest(NotI);
    Yval = Gtest;
    Rvalue = RvalueTrain;
    Err_ = nanmean(abs(predy-Gtest));
    STD_ = nanstd(predy-Gtest);
else %test data
    if 0
        RsqL1O = 1-nansum(weights.*(Gtest_-predyL1O).^2)/nansum(weights.*(Gtest_-nanmean(Gtest_)).^2);
        Rvalue = RsqL1O^0.5;
        ErrL1O = nanmean(abs(predyL1O-Gtest_));
        ErrL1Orms = nanmean((predyL1O-Gtest_).^2).^0.5;
        predy = predyL1O;% predyL1O_fsmr;
        Yval = Gtest_;
    else
        Rvalue = RsqL1O_(ftrnum)^0.5;
        ErrL1O = ErrL1O_(ftrnum);
        ErrL1Orms = ErrL1Orms_(ftrnum);
        predy = predyL1O_array(:,ftrnum);%
        Yval = Gtest_;
    end
end



%% Figure E2
% top left, performance during testing
% top right, proportion of absolute error
% bottom left, probability plot
% bottom right, confusion matrix

figure(82); clf(figure(82)); hold on;
fig = gcf; fig.Color = [1 1 1]; fig.Units = 'inches';
fig.Position = [0.5  0.5   9   9];
PlotUptoNFtrs = 100;
%%
subplot(2,2,1); % testing performance
if 1
plot([RL1O(1:PlotUptoNFtrs);ErrL1O_(1:PlotUptoNFtrs)]'); hold on;
currentYlim=ylim();
plot([ftrnum, ftrnum], [0, currentYlim(2)],'k-', 'linewidth', 1.5);
ax = gca; set(ax,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz);
xlabel('Number of Features', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
ylabel('Performance (Testing)', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow'); 
ylim([0 1]); xlim([0 100]); axis square
text(75, 0.88, {'correlation', 'coefficient'}, 'fontname','arial narrow','FontSize',10, 'FontWeight', 'normal');
text(75, 0.28, {'mean', 'absolute', 'error'}, 'fontname','arial narrow','FontSize',10, 'FontWeight', 'normal');
else
plot([RL1O(1:PlotUptoNFtrs);RsqL1O_(1:PlotUptoNFtrs);ErrL1O_(1:PlotUptoNFtrs);ErrL1Orms_(1:PlotUptoNFtrs)]');
xlabel('Number of Ftrs'); ylim([0 1]);
legend('RL1O','RsqL1O', 'ErrL1O', 'ErrL1Orms', 'location','east');
title('Performance Vs Number of Ftrs (OA Testing)');
fig = gcf;
fig.Color = [1 1 1]; % set background colour to white
fig.Units = 'inches';
fig.Position = [7 1 4 4];
%fig.Position = [-12 1 4.5 4.5];
axis square
str = ['..\Figures\', 'PerformanceVsNumFtrs_TestingAvg', datastr];
end
%%
subplot(2,2,2); % proportion of error
a_out = mean(ErrTrain_array(:,:),1);
err100 = a_out(100);
p_of100 = 100.*(err100 ./ a_out(1:100));
plot(100-p_of100, 'k-'); hold on;
currentYlim=ylim();
plot([ftrnum, ftrnum], [0, currentYlim(2)/2],'k-', 'linewidth', 1.5);
ax = gca; set(ax,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz);
xlabel('Number of Features', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
ylabel('Proportion of absolute error (%)', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow'); 
ylim([0 50]); xlim([0 100]); axis square
%%
subplot(2,2,3); % probability plot
if 1
    scatter(100*predy,100*Yval+10*randn(length(Yval),1),4,'filled','markerfacealpha',0.2, 'markerfacecolor', [0.1 0.1 0.1]);
    hold on;
    dx=0.2; xbins=[0:0.1:1];
    clear meanX meanY upperY lowerY
    for i=1:length(xbins)-1
        Ix=predy>xbins(i)&predy<xbins(i+1);
        meanX(i)=mean(predy(Ix));
        meanY(i)=nansum(Yval(Ix).*weights(Ix))/nansum(weights(Ix));
        upperY(i)=1*nanstd(Yval(Ix))+meanY(i);
        lowerY(i)=meanY(i)-1*nanstd(Yval(Ix));
    end
    ax = gca; set(ax,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz);axis square 
    yticks([0,100]);
    yticklabels({'Obstructed','Normal'}); ytickangle(ax,90);
    xlim([0 100]); axis square 
    %xlabel('Probability of being Normal or Obstructed using Flow Shape (%)');
    xlabel('Flow Shape Probability (%)', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');  
    ylabel('Gold Standard Classification', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
else % alternate probability plot, plotSpread
    Ndatapoints=Inf;
    if Ndatapoints>length(predy), Ndatapoints=length(predy); end
    data=100*predy(1:Ndatapoints);
    catIdx=Yval(1:Ndatapoints);
    data={data(catIdx~=1),data(catIdx==1)};
    %figure(5); clf(5);
    largeNfactor=2*(Ndatapoints/300)^0.5; %magic equation
    if largeNfactor<1,largeNfactor=1; end
    [~,data1,temp]=plotSpread(data,'xNames',{'F.limited','Normal'},'categoryIdx',catIdx,...
        'categoryMarkers',{'.','.'},'categoryColors',{'k','k'},'binWidth',0.1/largeNfactor,...
        'xyOri','flipped','magicNumber',1-(0.1/largeNfactor));
    ax = gca; set(ax,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz);axis square 
    xlabel('Pr Normal using flow shape LogReg, %')
    scatter(data1(:,2),100*(data1(:,1)-1),1,'filled','markerfacealpha',0.1)
    set(gca,'tickdir','out');
end
%%
subplot(2,2,4); %confusion matrix
g = NaN*Yval;
ghat = NaN*predy;
classcutoffs = [0.5]; %normal mild moderate severe v.severe / normal borderline mild moderate severe
Nclasses=length(classcutoffs)+1;
% DMconfusionMat, 1 = NFL, 0 = FL
g_ = Yval; % true/label/class % nnz(g_)
ghat_ = NaN*predy; % prediciton
ghat_(predy>classcutoffs(1))=1;
ghat_(predy<=classcutoffs(1))=0;
trueNFL = nnz(g_==1 & ghat_ ==1);
falseNFL = nnz(g_==0 & ghat_ ==1);
trueFL = nnz(g_==0 & ghat_ ==0);
falseFL = nnz(g_==1 & ghat_ ==0);
dlmXtab = [trueNFL; falseNFL; trueFL; falseFL]
t = (dlmXtab./sum(dlmXtab))*100;
C_Total = [t(1), t(4); t(2), t(3)];
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
C_Total2 = (C./sum(C(:)))*100; % line from lin reg
AccA_C_Factual = C_Factual.*diag(ones(1,length(C)));
AccB_C_Factual = C_Factual.*(diag(ones(1,length(C)-1),1)+diag(ones(1,length(C)-1),-1));
AccAactual = mean(sum(AccA_C_Factual,2));
AccBactual = mean(sum(AccA_C_Factual,2)) + mean(sum(AccB_C_Factual,2));
AccA_C_Festimated = C_Festimated.*diag(ones(1,length(C)));
AccB_C_Festimated = C_Festimated.*(diag(ones(1,length(C)-1),1)+diag(ones(1,length(C)-1),-1));
AccAestimated = mean(sum(AccA_C_Festimated));
AccBestimated = mean(sum(AccB_C_Festimated)) + mean(sum(AccB_C_Festimated));
ACCs = [AccAactual AccBactual; AccAestimated AccBestimated]
%second col gives accuracy if accepting next-category error as correct

x = order'; y = order';
if 1; C1 = C_Festimated; else; C1 = C_Factual; end
%C1 = C_Total; % using the total proportion of breaths in each division
%C = [0 1 2 3 ; 1 2 3 4 ; 2 3 4 5];
xSplit = diff(x)/2;                 % Find edge points
ySplit = diff(y)/2;
xEdges = [x(1)-xSplit(1) x(2:end)-xSplit x(end)+xSplit(end)];
yEdges = [y(1)-ySplit(1) y(2:end)-ySplit y(end)+ySplit(end)];
[XGrid, YGrid] = meshgrid(xEdges,yEdges);
YGrid = flipud(YGrid);              % To match expected behavior
XGrid = fliplr(XGrid);
C2 = [[C1 zeros(size(C1,1),1)] ; zeros(1,size(C1,2)+1)];% Last row/col ignored
customcmap = GetCustomColorMap('gray');
set(gcf,'colormap',customcmap);
pcolor(XGrid,YGrid,C2)
hold on                             % Plot original data points
[X,Y] = meshgrid(x,y);
C1 = flipud(C1);
C1 = fliplr(C1);
C1 = C1';
for i=1:size(C1,1)
    for j=1:size(C1,2)
        if C1(i,j)<18%((max(max(C1))-min(min(C1)))/2+min(min(C1)))
            textcolor=[1 1 1];
        else
            textcolor=[0 0 0];
        end
        text(x(i),y(j),num2str(round(C1(i,j),1)),'color',textcolor,'horizontalalignment','center','fontname','arial narrow')
    end
end
labeltext={'Normal','Obstructed'};
ax = gca; set(ax,'box','off','tickdir','out', 'fontname','arial narrow','FontSize', TickFntSz); axis square
yticks(y); yticklabels(gca,fliplr(labeltext));
ytickangle(ax,90);
xticks(x); xticklabels(gca,fliplr(labeltext));
ylabel('Gold Standard Classification', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow');
xlabel('Flow Shape Classification', 'FontSize', LabelFntSz, 'FontName', 'Arial Narrow'); 

str = ['..\Figures\', 'Figure_E2_50Ftrs'];
switch savefigas
    case 'saveasTIFF'; print(fig, str, '-dtiff', '-r1000');
    case 'saveasPNG'; saveas(fig, str, 'png'); 
    case 'saveasFIG'; savefig(str); 
end
if closefigs; close(fig); end


%% show training and testing performance - one plot
Rtrain_aray = RsqTrain_array.^0.5; 
figure(4); clf(figure(4));
subplot(1,2,1);
fig = gcf;
fig.Color = [1 1 1]; 
fig.Units = 'inches';
fig.Position = [1    1  12   4.5];
pts=unique(PtData.PT);
for pt = 1:length(pts)  
    %subplot(2,2,1)
    plot(Rtrain_aray(pt,:), 'k'); hold on;
    %xlabel('Number of Ftrs'); ylim([0.5 1]); title('R');axis square
    %subplot(2,2,2)
%     plot(RsqTrain_array(pt,:), 'k'); hold on;
    %xlabel('Number of Ftrs'); ylim([0.5 1]); title('Rsq');axis square
    %subplot(2,2,3)
    plot(ErrTrain_array(pt,:), 'k'); hold on;
    %xlabel('Number of Ftrs'); ylim([0 0.4]); title('Err');axis square
    %subplot(2,2,4)
    plot(ErrRmsTrain_array(pt,:), 'k'); hold on;
    %xlabel('Number of Ftrs'); ylim([0 0.4]); title('Errsq');axis square
end
xlabel('Number of Features'); ylabel('Performance (Training)'); ylim([0 1]); xlim([0 100]);
box off;
axis square
%ylim([0 1]); title('Performance Vs NumFtrs (OATraining)');axis square
%suptitle('Performance Vs Number of Ftrs (OA Training)');
% fig = gcf;
% fig.Color = [1 1 1]; % set background colour to white
% fig.Units = 'inches';
% fig.Position = [7 1 4 4];
% str = ['..\Figures\', 'PerformanceVsNumFtrs_OATrainingEach', datastr];
% saveas(fig, str, 'png'); %savefig(str);

subplot(1,2,2);
%plot([RL1O(1:PlotUptoNFtrs);RsqL1O_(1:PlotUptoNFtrs);ErrL1O_(1:PlotUptoNFtrs);ErrL1Orms_(1:PlotUptoNFtrs)]'); hold on;
plot([RL1O(1:PlotUptoNFtrs);ErrL1O_(1:PlotUptoNFtrs);ErrL1Orms_(1:PlotUptoNFtrs)]'); hold on;
currentYlim=ylim();
plot([ftrnum, ftrnum], [0, currentYlim(2)],'k-', 'linewidth', 2);
xlabel('Number of Features'); ylabel('Performance (Testing)'); ylim([0 1]); xlim([0 100]);
box off;
legend('r', 'MAE', 'RMSE', 'location','east');
% title('Performance Vs Number of Features (Pnasal Testing)');
axis square
str = ['..\Figures\', 'Figure_E2'];

% Add labels A B to plot space
subplot(1,2,1); hold on;
text(-20, 0.95, 'A', 'FontSize', 20, 'FontWeight', 'Bold');
subplot(1,2,2); hold on;
text(-20, 0.95, 'B', 'FontSize', 20, 'FontWeight', 'Bold');

if savefigs; saveas(fig, str, 'png'); end %savefig(str);
if closefigs; close(fig); end


%% Analysis and plots
figure(2); clf(figure(2)); subplot(1,2,1); axis square;
fig = gcf;
fig.Color = [1 1 1]; 
fig.Units = 'inches';
fig.Position = [1    1  12   4.5];
if LinearReg
    scatter(100*predy,100*Yval,2,'filled','markerfacealpha',0.4);
    hold on;
    %dx=0.1; xbins=[0 0.15:dx:1.05 1.5];
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
    %medianX = [0.2:dx:1]; % overwrite
    plotbarwitherrorsmedian(100*medianY,100*medianX,100*upperIQRY2,100*lowerIQRY2,zeros(length(medianY),3),NaN,100*0.005,100*0.01);
    plotbarwitherrorsmedian(100*medianY,100*medianX,100*upperIQRY,100*lowerIQRY,zeros(length(medianY),3),NaN,100*dx/2,100*0.01);
    xlim([0 150]); xlabel('Pr using flow shape LinReg, FL');
else
    scatter(100*predy,100*Yval+10*randn(length(Yval),1),4,'filled','markerfacealpha',0.2, 'markerfacecolor', [0.1 0.1 0.1]);
    hold on;
    dx=0.2; xbins=[0:0.1:1];
    clear meanX meanY upperY lowerY
    for i=1:length(xbins)-1
        Ix=predy>xbins(i)&predy<xbins(i+1);
        meanX(i)=mean(predy(Ix));
        meanY(i)=nansum(Yval(Ix).*weights(Ix))/nansum(weights(Ix));
        upperY(i)=1*nanstd(Yval(Ix))+meanY(i);
        lowerY(i)=meanY(i)-1*nanstd(Yval(Ix));
    end
    %plot(100*meanX,100*meanY,'k-');
    %plot(100*meanX,100*upperY,'k--');
    %plot(100*meanX,100*lowerY,'k--');
    yticks([0,100]);
    yticklabels({'Obstructed','Normal'});
    xlim([0 100]);  xlabel('Probability of being Normal or Obstructed using Flow Shape (%)')
end
axis square;
% Add labels A B to plot space
text(-20, 130, 'A', 'FontSize', 20, 'FontWeight', 'Bold');
%title(['Rvalue: ' num2str(Rvalue(end),4)])
str=['R value is ', (num2str(Rvalue(end),4))]; disp(str);
   

% create classes from regression results
g = NaN*Yval;
ghat = NaN*predy;
subplot(1,2,2); hold on; axis square;
if LinearReg
    classcutoffs = [0.9 0.5];
    labeltext={'Normal','Mild-Moderate','Severe'};
    %classcutoffs = [0.9 0.7];
    %classcutoffs = [0.7]; %normal mild moderate severe v.severe / normal borderline mild moderate severe
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
    
    %         AccA_C_Factual = C_Factual.*diag(ones(1,length(C)));
    %         AccB_C_Factual = C_Factual.*(diag(ones(1,length(C)-1),1)+diag(ones(1,length(C)-1),-1));
    %         AccAactual = mean(sum(AccA_C_Factual,2));
    %         AccBactual = mean(sum(AccA_C_Factual,2)) + mean(sum(AccB_C_Factual,2));
    %         AccA_C_Festimated = C_Festimated.*diag(ones(1,length(C)));
    %         AccB_C_Festimated = C_Festimated.*(diag(ones(1,length(C)-1),1)+diag(ones(1,length(C)-1),-1));
    %         AccAestimated = mean(sum(AccA_C_Festimated));
    %         AccBestimated = mean(sum(AccB_C_Festimated)) + mean(sum(AccB_C_Festimated));
    %         ACCs = [AccAactual AccBactual; AccAestimated AccBestimated]
    
    
    x = order'; y = order';
    if 1; C1 = C_Festimated; else; C1 = C_Factual; end
    %C = [0 1 2 3 ; 1 2 3 4 ; 2 3 4 5];
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
    colorbar();
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
    labeltext={'Normal','Mild','Moderate','Severe','V.Severe'};
    
else
    classcutoffs = [0.5]; %normal mild moderate severe v.severe / normal borderline mild moderate severe
    Nclasses=length(classcutoffs)+1;
    
    % DMconfusionMat, 1 = NFL, 0 = FL
    g_ = Yval; % true/label/class % nnz(g_)
    ghat_ = NaN*predy; % prediciton
    ghat_(predy>classcutoffs(1))=1;
    ghat_(predy<=classcutoffs(1))=0;
    trueNFL = nnz(g_==1 & ghat_ ==1);
    falseNFL = nnz(g_==0 & ghat_ ==1);
    trueFL = nnz(g_==0 & ghat_ ==0);
    falseFL = nnz(g_==1 & ghat_ ==0);
    dlmXtab = [trueNFL; falseNFL; trueFL; falseFL]
    t = (dlmXtab./sum(dlmXtab))*100;
    C_Total = [t(1), t(4); t(2), t(3)];
    
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
    C_Total2 = (C./sum(C(:)))*100; % line from lin reg
    AccA_C_Factual = C_Factual.*diag(ones(1,length(C)));
    AccB_C_Factual = C_Factual.*(diag(ones(1,length(C)-1),1)+diag(ones(1,length(C)-1),-1));
    AccAactual = mean(sum(AccA_C_Factual,2));
    AccBactual = mean(sum(AccA_C_Factual,2)) + mean(sum(AccB_C_Factual,2));
    AccA_C_Festimated = C_Festimated.*diag(ones(1,length(C)));
    AccB_C_Festimated = C_Festimated.*(diag(ones(1,length(C)-1),1)+diag(ones(1,length(C)-1),-1));
    AccAestimated = mean(sum(AccA_C_Festimated));
    AccBestimated = mean(sum(AccB_C_Festimated)) + mean(sum(AccB_C_Festimated));
    ACCs = [AccAactual AccBactual; AccAestimated AccBestimated]
    %second col gives accuracy if accepting next-category error as correct
    
    x = order'; y = order';
    if 1; C1 = C_Festimated; else; C1 = C_Factual; end
    %C1 = C_Total; % using the total proportion of breaths in each division
    %C = [0 1 2 3 ; 1 2 3 4 ; 2 3 4 5];
    xSplit = diff(x)/2;                 % Find edge points
    ySplit = diff(y)/2;
    xEdges = [x(1)-xSplit(1) x(2:end)-xSplit x(end)+xSplit(end)];
    yEdges = [y(1)-ySplit(1) y(2:end)-ySplit y(end)+ySplit(end)];
    [XGrid, YGrid] = meshgrid(xEdges,yEdges);
    YGrid = flipud(YGrid);              % To match expected behavior
    XGrid = fliplr(XGrid);
    C2 = [[C1 zeros(size(C1,1),1)] ; zeros(1,size(C1,2)+1)];% Last row/col ignored
    %customcmap = GetCustomColorMap('SS');
    customcmap = GetCustomColorMap('gray');
    set(gcf,'colormap',customcmap);
    pcolor(XGrid,YGrid,C2)
    hold on                             % Plot original data points
    [X,Y] = meshgrid(x,y);
    %h=colorbar();
    %set(gcf,'colormap',customcmap);
    %set(h,'FontName','Arial Narrow','Limits',[0 max(max(C2))]);
    C1 = flipud(C1);
    C1 = fliplr(C1);
    C1 = C1';
    for i=1:size(C1,1)
        for j=1:size(C1,2)
            if C1(i,j)<18%((max(max(C1))-min(min(C1)))/2+min(min(C1)))
                textcolor=[1 1 1];
            else
                textcolor=[0 0 0];
            end
            text(x(i),y(j),num2str(round(C1(i,j),1)),'color',textcolor,'horizontalalignment','center','fontname','arial narrow')
        end
    end
    labeltext={'Normal','Obstructed'};
    
end
axis square;
yticks(y); yticklabels(gca,fliplr(labeltext));
xticks(x); xticklabels(gca,fliplr(labeltext));
xlabel('Flow Shape Classification');
ylabel('Actual Classification');

 

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


%% some plots - check univariate performance
ftr = 112;
figure(1);clf(figure(1));scatter((PtData.VE./PtData.Veup), Amatrix2(:,ftr),2,'filled','markerfacealpha',0.4);
xlabel('VEVeup'); xlim([-0.05 3.5]);


%% Final model using all data, uses selected N optimal features (NfeaturesOpt)
%not tested
if 0
    If = 1:size(Amatrix2,2);
    Labels = Labels_;
    MaxNfeatures=NfeaturesOpt;
    while length(If)>=MaxNfeatures
        if LinearReg
            [Rsq,Pvals,RMSE,b]=glmfitFast(Amatrix2(:,If),Gtest_,weights,1); %faster
            Pvals(1)=[];
        else
            if SS_logitreg
                [b,dev,stats]=glmfit(Amatrix2(:,If),Gtest_,'binomial','weights',weights);
            else
                [b,dev,stats]=glmfit(Amatrix2(:,If),Gtest_,'binomial','weights',weights, 'link','logit');
            end
            Pvals = stats.p(2:end);
        end
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
    figure(3); clf(figure(3));
    scatter(Gtest_, mdl_prediction,2,'filled','markerfacealpha',0.4);
    xlabel('Gtest'); ylabel('Mdl Pred');
end




%%







%%






%%






%%
%AAsummary = [alpha MaxNfeatures median(Nfeatures) Fwd useweights median(RvalueTrain) Rvalue ...
%    ACCs(2,:) median(Err) median(ErrRms) ErrL1O ErrL1Orms];

%% Patient level summary data
%Need to put back breaths that are thrown out because of low VE here (apneas)
for i=1:max(PtData.PT)
    I=PtData.PT==i&PtData.NotAr==1;
    medianG(i) = nanmedian(Gtest_(I));
    medianGest(i) = nanmedian(predy(I));
    
    denominator = sum(~isnan(Gtest_(I)));
    
    Fmildmod(i) = sum(Gtest_(I)<0.9&Gtest_(I)>0.5)/denominator;
    Fmildmodest(i) = sum(predy(I)<0.9&Gtest_(I)>0.5)/denominator;
    
    Fmodsev(i) = sum(Gtest_(I)<0.7)/denominator;
    Fmodsevest(i) = sum(predy(I)<0.7)/denominator;
    
    Fsev(i) = sum(Gtest_(I)<0.5)/denominator;
    Fsevest(i) = sum(predy(I)<0.5)/denominator;
end

corr(medianGest',medianG')
corr(Fmodsevest',Fmodsev')
corr(Fsevest',Fsev')
corr(Fmildmodest',Fmildmod')
figure(9); clf(9);
subplot(1,5,1);
[Rtemp,Ptemp]=plotregressionwithSEM(medianGest,medianG); title(num2str(Rtemp));
subplot(1,5,2);
[Rtemp,Ptemp]=plotregressionwithSEM(Fmodsevest,Fmodsev); title(num2str(Rtemp));
subplot(1,5,3);
[Rtemp,Ptemp]=plotregressionwithSEM(Fsevest,Fsev); title(num2str(Rtemp));
subplot(1,5,4);
[Rtemp,Ptemp]=plotregressionwithSEM(Fmildmodest,Fmildmod); title(num2str(Rtemp));

%% Correlate with AHI
figure(10); clf(10);
subplot(1,1,1);
[Rtemp,Ptemp]=plotregressionwithSEM(AHItotal',medianGest); title(num2str(Rtemp));

[b,dev,stats]=glmfit([AHItotal medianGest'],medianG')
[b,dev,stats]=glmfit([AHItotal Fsevest'],Fsev')

%% Repeat, sensitivity for Nfeatures
for n=1:100
    clear medianG medianGest Fmildmod Fmildmodest Fmodsev Fmodsevest Fsev Fsevest
    predy = predyL1O_array(:,n);
    for i=1:max(PtData.PT)
        I=PtData.PT==i&PtData.NotAr==1;
        medianG(i) = nanmedian(Gtest_(I));
        medianGest(i) = nanmedian(predy(I));
        
        denominator = sum(~isnan(Gtest_(I)));
        
        Fmildmod(i) = sum(Gtest_(I)<0.9&Gtest_(I)>0.5)/denominator;
        Fmildmodest(i) = sum(predy(I)<0.9&Gtest_(I)>0.5)/denominator;
        
        Fmodsev(i) = sum(Gtest_(I)<0.7)/denominator;
        Fmodsevest(i) = sum(predy(I)<0.7)/denominator;
        
        Fsev(i) = sum(Gtest_(I)<0.5)/denominator;
        Fsevest(i) = sum(predy(I)<0.5)/denominator;
    end
    n
    Rvals(n,:)=[...
        corr(medianGest',medianG'),...
        corr(Fmodsevest',Fmodsev'),...
        corr(Fsevest',Fsev'),...
        corr(Fmildmodest',Fmildmod'),...
        ];
end



%% Load subject data for plot
predy(predy<0)=0;

%Plot FL values over time for an example subject (e.g. compare with Spike data). Random check (n=8) Looks great.
n=2; % 1313?
%for n=1:30
%ypredAll = predict(SVMModel_reg,Amatrix2(:,If));
%ypred = predict(SVMModel_reg,Amatrix2(:,If));
Data1 = [PtData.BB_time(PtData.PT==n) PtData.BB_Ttot(PtData.PT==n) predy(PtData.PT==n) Gtest_(PtData.PT==n)];

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

figure(8); clf(8) %ax(1) = subplot(2,1,1);
ax(1)=subplot(3,1,1);
stairs(Data1(:,1),100*Data1(:,4),'k');
hold('on')
stairs(Data1(:,1),100*Data1(:,3),'r');
plot([Data1(1,1) Data1(end,1)],[0 0],'k:')
plot([Data1(1,1) Data1(end,1)],100*[1 1],'k:')

set(gca,'xtick',[],'box','off')


load('J:\PEOPLE\FACULTY\SANDS\Phenotyping from Pes and Flow\SpikeFiles\1313.mat', 'Edi', 'Flow');
Time = 0:0.008:(length(Flow.values)-1)*0.008;

ax(2)=subplot(3,1,2); plot(Time,Flow.values,'g'); hold('on');
dsf=5; dt=Flow.interval;
FlowF=Flow.values;
if 1
    filter_HFcutoff_butter0 = 12.5;
    filter_order0 = 1;
    [B_butter0,A_butter0] = butter(filter_order0,[filter_HFcutoff_butter0]/(1/dt/2),'low');
    FlowF = filtfilt(B_butter0,A_butter0,FlowF); %filtfilt, otherwise flow signal is right-shifted
end
ax(2)=subplot(3,1,2); plot(downsample(Time,dsf),downsample(FlowF,dsf),'k');
set(gca,'xtick',[],'box','off')
ax(3)=subplot(3,1,3); plot(Time,Edi.values);
box('off');
linkaxes(ax,'x');

xlim([16200 16900])

%%

figure()
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

