function [Ilist,Perf_Train,Perf_Test,SVMModel_Train, SVMModel_Test]=svmforwardselect_perpatient_regression(...
    Amatrix,criteriaR,tempignorevars,Yvariable,maxNfeaturesJ,rangevalidate,...
    forcedvars, weights, prior_Ilist, prior_perf1, prior_perf2)

% Written by Scott Sands 2015-10-28
% Modified by dwayne mann, 2017-06-28

%funcsensspecstr='y+(1-x)';
%breakifdeltacostlessthan=0.001;

warning('off', 'all');

if length(maxNfeaturesJ)==1
    maxNfeatures=maxNfeaturesJ;
    if maxNfeatures<0||maxNfeatures>size(Amatrix,2)
        maxNfeatures=size(Amatrix,2);
    end
end
criteriaR_all = criteriaR(:); % keep full list for later
criteriaR_train=criteriaR_all;
criteriaR_train(rangevalidate)=NaN; % set up the training range
criteriaR_leftout=criteriaR_all;
criteriaR_leftout(~rangevalidate)=NaN; % set up the left-out range
if length(maxNfeaturesJ)>1
    maxNfeatures=ceil(median(maxNfeaturesJ));
end
SVMModel_Train{maxNfeatures} = [];
SVMModel_Test{maxNfeatures} = [];

%% Handle prior learning (if it exists)
if ~isempty(prior_Ilist)
    % if we have prior learning, then let's know about it
    Ilist=prior_Ilist;
    Perf_Train=prior_perf1;
    Perf_Test=prior_perf2;
    % return if we've already learnt all that we can
    complete=maxNfeaturesJ==length(Ilist);
    if complete
        return
    else % we still have a bit to learn, so let's continue processing,
        % picking up from the end of the prior learning
        n = length(Ilist)+1;
    end
else
    Ilist=[];
    Perf_Train=[];
    Perf_Test=[];
    n=1;
end

for jj=n:maxNfeatures %fwd stepwise, with this loop now starting at n
    RHO_=NaN*ones(1,size(Amatrix,2));
    PVAL_=NaN*ones(1,size(Amatrix,2));
    
    str=(['Starting parfor loop for cycle ', num2str(jj), ' of ', num2str(maxNfeatures) ]); disp(str);
    tic
    %% where possible, this should run as a parfor loop
    parfor kk=1:size(Amatrix,2) %run through each variable/column
        warning('off', 'all');
        % some checks, such as we can't have the same ftr appearing twice
        if sum(kk==Ilist)||sum(kk==tempignorevars)||ismember(kk,Ilist)
            RHO_(kk)=NaN;
            PVAL_(kk)=NaN;
            continue
        end

        % skip this feature if it's more than 5% NaN or inf
        if nnz(~isfinite(Amatrix(:,kk)))> size(Amatrix,1)*0.05
            %str=['Skipping feature ', num2str(kk)];disp(str);
            continue
        end
        
        %% make a model and use it to make a prediction
        % if Ilist is empty, our test list is simply the current feature
        % otherwise, add the current feature to our existing list
        Ilisttest = [Ilist kk];

        LinearRModel1 = fitrlinear(Amatrix(:,Ilisttest),criteriaR_train(:),'Learner','svm','Regularization','ridge', 'Weights', weights);           
        PredTsvmLinear = predict(LinearRModel1,Amatrix(:,Ilisttest));            
        Ytest = Yvariable;
        [Ytest, PredTsvmLinear] = RemoveNaNs(Ytest, PredTsvmLinear);
        
        [RHO_(kk),PVAL_(kk)] = corr(Ytest, PredTsvmLinear, 'rows', 'all');
        RHO_(kk) = (RHO_(kk))^2;
        
        if 0
            figure(1); clf(figure(1));
            scatter(Ytest, PredTsvmLinear);
            str = ['(fitrlinear) Feature ', num2str(Ilisttest)];
            title(str);
            
            figure(2); clf(figure(2));
            scatter(Yvariable,Amatrix(:,kk));
            str = ['(Raw) Feature ', num2str(kk)];
            title(str);
            
            LinearRModel_simple = fitlm(Amatrix(:,Ilisttest),criteriaR_train(:),'Weights', weights);
            PredTsvmLinear_simple = predict(LinearRModel_simple,Amatrix(:,Ilisttest));
            
            figure(3); clf(figure(3));
            scatter(Yvariable,PredTsvmLinear_simple);
            str = ['(fitlm)Feature ', num2str(Ilisttest)];
            title(str);
        end
    end
    
    toc
    %% find the next best feature to add to the list
    % using the performance from the training data, find the feature with
    % the minimum cost, as determined by the specified cost function

    % Note: should this should be using the performance in the test data,
    % i.e. the model is trained in the training data, but cross-validated
    % with the test data. the performance in this test data directs the
    % next iteration of training.

    [v, temp3i] = max(RHO_);
          
    if jj<=length(forcedvars)
        temp3i=forcedvars(jj);
    end
    
    %% Break processing if we are no longer improving (optional)

    %         if jj>1&&jj>length(forcedvars)&&(CostList(end-1)-CostList(end))<=breakifdeltacostlessthan %min arbitrary 1% change in AUC
    %             break
    %         end

    %% add best feature to Ilist
    Ilist=[Ilist temp3i];

    %% record training and validation performance

    % make a model using the training data, and the current Ilisttest
    SVMModelA = fitrlinear(Amatrix(:,Ilist),criteriaR_train(:),'Learner','svm','Regularization','ridge', 'Weights', weights);
    SVMModel_Train{jj} = SVMModelA;
    PredTsvmLinearA = predict(SVMModelA,Amatrix(:,Ilist));  
    Ytest = Yvariable;
    [Ytest, PredTsvmLinearA] = RemoveNaNs(Ytest, PredTsvmLinearA);
    [RHO,~] = corr(Ytest, PredTsvmLinearA, 'rows', 'all');
    Perf_Train{jj} = RHO;
    
    % make a model using the test data, and the current Ilisttest
    SVMModelB = fitrlinear(Amatrix(:,Ilist),criteriaR_leftout(:),'Learner','svm','Regularization','ridge', 'Weights', weights);
    SVMModel_Test{jj} = SVMModelB;
    PredTsvmLinearB = predict(SVMModelB,Amatrix(:,Ilist));  
    Ytest = Yvariable;
    [Ytest, PredTsvmLinearB] = RemoveNaNs(Ytest, PredTsvmLinearB); 
    [RHO,~] = corr(Ytest, PredTsvmLinearB, 'rows', 'all');
    Perf_Test{jj} = RHO;
end

if 0
Ilist = [80, 62, 66];
SVMModel_plot = fitrlinear(Amatrix(:,Ilist(1:3)),criteriaR(:),'Learner','svm','Regularization','ridge', 'Weights', weights);
SVMModel_all = fitrlinear(Amatrix(:,Ilist),criteriaR(:),'Learner','svm','Regularization','ridge', 'Weights', weights);
PredTsvmLinAll = predict(SVMModel_all,Amatrix(:,Ilist)); 
Ytest = Yvariable;
Ytest(isnan(PredTsvmLinAll))=[];
PredTsvmLinAll(isnan(PredTsvmLinAll))=[];
[RHO,PVAL] = corr(Ytest, PredTsvmLinAll, 'rows', 'all');
figure(3); clf(figure(3));
scatter(Ytest, PredTsvmLinAll);
svm_3d_matlab_vis_mod(SVMModel_plot,Amatrix(:,Ilist(1:3)),criteriaR,0);
end

end

function [DataA, DataB] = RemoveNaNs(DataA, DataB)
exclude = isnan(DataA);
DataA(exclude)=[];
DataB(exclude)=[];
exclude = isnan(DataB);
DataA(exclude)=[];
DataB(exclude)=[];
end