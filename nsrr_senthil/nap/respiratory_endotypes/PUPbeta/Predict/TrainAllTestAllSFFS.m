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

Restarting = 0; % this is for troubleshooting, and starting at pt 2 or 3 or ...

progressbar(0,0,0);
progressbar('Pt','Step','Ftrs');
try
    warning ('off','all');
    for subj=1%:54 % subj=2 subj=3  
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
            Isubj=1==zeros(length(Amatrix2_flow),1);  % Isubj is the logical index of the L1O patient, can use same indx for all now
            
            % set up training data - flow (always)
            % training data is everything except the current pt
            Gtest_train = HS_cat(~Isubj);
            weights_train = weights(~Isubj);
            colofones_train = ones(nnz(~Isubj),1);
            Amatrix2_train = Amatrix2_flow(~Isubj,:);
            PtData_flow_train = PtData_flow.PT_f(~Isubj,:);
            
            % set up test data - flow and pnasal
            % test data is only the current pt
            
            Isubj=0==zeros(length(Amatrix2_flow),1);  % Isubj is the logical index of the L1O patient, can use same indx for all now
            
            Gtest_test=HS_cat(Isubj);
            weights_test = weights(Isubj);
            colofones_test = ones(nnz(Isubj),1);
            Amatrix2_flow_test = Amatrix2_flow(Isubj,:);
            Amatrix2_pnasal_test = Amatrix2_pnasal(Isubj,:);
            
            Isubj=1==zeros(length(Amatrix2_flow),1);  % Isubj is the logical index of the L1O patient, can use same indx for all now
            
            
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
            
            %%
            settings.verbose=1;
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