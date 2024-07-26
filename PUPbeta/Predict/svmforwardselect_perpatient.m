function [Ilist,Perf_Train,Perf_Test,SVMModel_Train, SVMModel_Test]=svmforwardselect_perpatient(...
    Amatrix,criteriaR,tempignorevars,Yvariable,maxNfeaturesJ,rangevalidate,...
    forcedvars,coststr,svpmethodstr,kernelfunction,sigma,weights, ...
    prior_Ilist, prior_perf1, prior_perf2)

% Written by Scott Sands 2015-10-28
% Modified by dwayne mann, 2017-06-28

%funcsensspecstr='y+(1-x)';
%breakifdeltacostlessthan=0.001;

warning off all;
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
performance{size(Amatrix,2),1}=[];

CostList=[];
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
    TP_=NaN*ones(1,size(Amatrix,2));
    FP_=NaN*ones(1,size(Amatrix,2));
    TN_=NaN*ones(1,size(Amatrix,2));
    FN_=NaN*ones(1,size(Amatrix,2));
    PPV_=NaN*ones(1,size(Amatrix,2));
    NPV_=NaN*ones(1,size(Amatrix,2));
    Sens_=NaN*ones(1,size(Amatrix,2));
    Spec_=NaN*ones(1,size(Amatrix,2));
    SensplusSpec_=NaN*ones(1,size(Amatrix,2));
    Acc_=NaN*ones(1,size(Amatrix,2));
    phighvslow=NaN*ones(1,size(Amatrix,2));
    phighvslow_=NaN*ones(1,size(Amatrix,2));
    absdeltamean=NaN*ones(1,size(Amatrix,2));
    absdeltamedian=NaN*ones(1,size(Amatrix,2));

    Sens=NaN*ones(1,size(Amatrix,2));
    Spec=NaN*ones(1,size(Amatrix,2));
    SensplusSpec=NaN*ones(1,size(Amatrix,2));
    Acc=NaN*ones(1,size(Amatrix,2));

    str=(['Starting parfor loop for cycle ', num2str(jj), ' of ', num2str(maxNfeatures) ]); disp(str);
    tic
    %% where possible, this should run as a parfor loop
    parfor kk=1:size(Amatrix,2) %run through each variable/column
        % some checks, such as we can't have the same ftr appearing twice
        if sum(kk==Ilist)||sum(kk==tempignorevars)||ismember(kk,Ilist)
            TP_(kk)=NaN;
            FP_(kk)=NaN;
            TN_(kk)=NaN;
            FN_(kk)=NaN;
            PPV_(kk)=NaN;
            NPV_(kk)=NaN;
            Sens_(kk)=NaN;
            Spec_(kk)=NaN;
            SensplusSpec_(kk)=NaN;
            Acc_(kk)=NaN;
            phighvslow(kk)=NaN;
            phighvslow_(kk)=NaN;
            absdeltamean(kk)=NaN;
            absdeltamedian(kk)=NaN;
            continue
        end

        %% make a model and use it to make a prediction
        % if Ilist is empty, our test list is simply the current feature
        % otherwise, add the current feature to our existing list
        Ilisttest = [Ilist kk];

        % new fitcsvm method, and also includes observation weights.
        % make a model using the training data, and the current Ilisttest
        SVMModel1 = fitcsvm(Amatrix(:,Ilisttest),criteriaR_train(:),'Standardize',true,'Solver',svpmethodstr,...
            'KernelFunction',kernelfunction,'KernelScale',sigma, 'Weights', weights);

        % then make a prediction using the current model, for all the data
        PredTsvm = predict(SVMModel1,Amatrix(:,Ilisttest));

        %% assess how good the prediction was
        % determine the performance in the training data
        TP_(kk) = nansum(criteriaR_train(:).*PredTsvm);
        FP_(kk) = nansum((1-criteriaR_train(:)).*(PredTsvm));
        TN_(kk) = nansum((1-criteriaR_train(:)).*(1-PredTsvm));
        FN_(kk) = nansum(criteriaR_train(:).*(1-PredTsvm));
        PPV_(kk) = TP_(kk)/(TP_(kk)+FP_(kk));
        NPV_(kk) = TN_(kk)/(TN_(kk)+FN_(kk));
        Sens_(kk) = TP_(kk)/(TP_(kk)+FN_(kk));
        Spec_(kk) = TN_(kk)/(TN_(kk)+FP_(kk));
        SensplusSpec_(kk) = Sens_(kk)+ Spec_(kk);
        Acc_(kk) = (TP_(kk)+TN_(kk))/(TP_(kk)+FP_(kk)+TN_(kk)+FN_(kk));
        % note: the above is more or less equivalent to the next line
        % performance_train{kk} = PredictiveValue(criteriaR_train,PredTsvm,Yvariable);

        try
            [~,phighvslow(kk)]=ttest2(Yvariable(PredTsvm==1),Yvariable(PredTsvm==0));
            [phighvslow_(kk)]=ranksum(Yvariable(PredTsvm==1),Yvariable(PredTsvm==0));
            absdeltamean(kk) = abs(nanmean(Yvariable(PredTsvm==1))-nanmean(Yvariable(PredTsvm==0)));
            absdeltamedian(kk) = abs(nanmedian(Yvariable(PredTsvm==1))-nanmedian(Yvariable(PredTsvm==0)));
        catch me
            phighvslow(kk)=NaN;
            phighvslow_(kk)=NaN;
            absdeltamean(kk)=NaN;
            absdeltamedian(kk)=NaN;
        end
      
        % detemine the performance in the left-out data
        performance{kk} = PredictiveValue(criteriaR_leftout,PredTsvm,Yvariable);% PredT and Y should be limited to the left out pt only
        Sens(kk) = performance{kk}.Sens_sem_chance_p(1);
        Spec(kk) = performance{kk}.Spec_sem_chance_p(1);
        SensplusSpec(kk) = Sens(kk)+ Spec(kk);
        Acc(kk) = performance{kk}.Acc_sem_chance_p(1);

    end
    toc
    %% find the next best feature to add to the list
    % using the performance from the training data, find the feature with
    % the minimum cost, as determined by the specified cost function

    % Note: should this should be using the performance in the test data,
    % i.e. the model is trained in the training data, but cross-validated
    % with the test data. the performance in this test data directs the
    % next iteration of training.

    tempcost=eval(coststr); %tempp, 1-tempAUC, -tempmaxsensspec
    if sum(~isnan(tempcost))==0 %because min_index of [NaN NaN] is not NaN...
        break
    end
    [~,temp3i]=min(tempcost);
    if jj<=length(forcedvars)
        temp3i=forcedvars(jj);
    end
    CostList=[CostList tempcost(temp3i)];

    %% Break processing if we are no longer improving (optional)

    %         if jj>1&&jj>length(forcedvars)&&(CostList(end-1)-CostList(end))<=breakifdeltacostlessthan %min arbitrary 1% change in AUC
    %             break
    %         end

    %% add best feature to Ilist
    Ilist=[Ilist temp3i];

    %% record training and validation performance

    % record the performance in the training data
    Perf_Train=[Perf_Train; PPV_(temp3i) NPV_(temp3i) Sens_(temp3i) Spec_(temp3i) Acc_(temp3i) TP_(temp3i) FP_(temp3i) TN_(temp3i) FN_(temp3i) phighvslow(temp3i) phighvslow_(temp3i) absdeltamean(temp3i) absdeltamedian(temp3i)];

    % record the performance in the left-out (validation) data
    Perf_Test=[Perf_Test; ...
        performance{temp3i}.PPV_sem_chance_p(1)...
        performance{temp3i}.NPV_sem_chance_p(1)...
        performance{temp3i}.Sens_sem_chance_p(1)...
        performance{temp3i}.Spec_sem_chance_p(1)...
        performance{temp3i}.Acc_sem_chance_p(1)...
        performance{temp3i}.TP_FP_TN_FN...
        performance{temp3i}.phighvslow_ttest...
        performance{temp3i}.phighvslow_ranksum...
        performance{temp3i}.predY_mean-performance{temp3i}.predN_mean...
        performance{temp3i}.predY_median-performance{temp3i}.predN_median];
    
    % make a model using the training data, and the current Ilisttest
    SVMModelA = fitcsvm(Amatrix(:,Ilist),criteriaR_train(:),'Standardize',true,'Solver',svpmethodstr,...
        'KernelFunction',kernelfunction,'KernelScale',sigma, 'Weights', weights);
    SVMModel_Train{jj} = SVMModelA;
    
    % make a model using the test data, and the current Ilisttest
    SVMModelB = fitcsvm(Amatrix(:,Ilist),criteriaR_leftout(:),'Standardize',true,'Solver',svpmethodstr,...
        'KernelFunction',kernelfunction,'KernelScale',sigma, 'Weights', weights); 
    SVMModel_Test{jj} = SVMModelB;
end

end
