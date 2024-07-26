function [PPVNPVSensSpec,table1,YvariableByPredictedOutcome,PredT,Ilist,PPVNPVSensSpecAll,YvariableByPredictedOutcomeAll,PredTAll]=svmleave1out(Amatrix,criteriaR,tempignorevars,Yvariable,maxNfeaturesJ,breakifdeltacostlessthan,rangevalidate,excludefortraining,forcedvars,coststr,svpmethodstr,kernelfunction,sigma)

%Written by Scott Sands 2015-10-28

%funcsensspecstr='y+(1-x)';
%breakifdeltacostlessthan=0.001;
warning off all;
if length(maxNfeaturesJ)==1
    maxNfeatures=maxNfeaturesJ;
    if maxNfeatures<0||maxNfeatures>size(Amatrix,2)
        maxNfeatures=size(Amatrix,2);
    end
end
criteriaR_ = criteriaR(:);
Exclude=[];
clear AUClistX IlistX B1array
PredT = NaN*(1:length(criteriaR))';
for Exclude=rangevalidate
    Ilist=[]; AUClist=[]; plist=[]; CostList=[];
    criteriaR=criteriaR_;
    criteriaR(Exclude)=NaN; %leave1out
    criteriaR(excludefortraining)=NaN; %leaveoutalways
    if length(maxNfeaturesJ)>1
        maxNfeatures=maxNfeaturesJ(Exclude);
    end
    for jj=1:maxNfeatures %fwd stepwise
        PPV_=NaN*ones(1,size(Amatrix,2));
        NPV_=NaN*ones(1,size(Amatrix,2));
        Sens_=NaN*ones(1,size(Amatrix,2));
        SensplusSpec_=NaN*ones(1,size(Amatrix,2));
        Acc_=NaN*ones(1,size(Amatrix,2));
        phighvslow=NaN*ones(1,size(Amatrix,2));
        phighvslow_=NaN*ones(1,size(Amatrix,2));
        absdeltamean=NaN*ones(1,size(Amatrix,2));
        absdeltamedian=NaN*ones(1,size(Amatrix,2));
        for kk=1:size(Amatrix,2) %run through each variable/column
            if sum(kk==Ilist)||sum(kk==tempignorevars)
                PPV_(kk)=NaN;
                NPV_(kk)=NaN;
                Sens_(kk)=NaN;
                SensplusSpec_(kk)=NaN;
                Acc_(kk)=NaN;
                phighvslow(kk)=NaN;
                phighvslow_(kk)=NaN;
                absdeltamean(kk)=NaN;
                absdeltamedian(kk)=NaN;
                continue
            end
            Ilisttest = [Ilist kk];
            SVMModel = svmtrain(Amatrix(:,Ilisttest),criteriaR(:),'method',svpmethodstr,'kernel_function',kernelfunction,'rbf_sigma',sigma);
            PredTsvm = svmclassify(SVMModel,Amatrix(:,Ilisttest));

            TP_ = nansum(criteriaR(:).*PredTsvm);
            FP_ = nansum((1-criteriaR(:)).*PredTsvm);
            TN_ = nansum((1-criteriaR(:)).*(1-PredTsvm));
            FN_ = nansum(criteriaR(:).*(1-PredTsvm));
            PPV_(kk) = TP_/(TP_+FP_);
            NPV_(kk) = TN_/(TN_+FN_);
            Sens_(kk) = TP_/(TP_+FN_);
            Spec_(kk) = TN_/(TN_+FP_);
            SensplusSpec_(kk) = Sens_(kk)+ Spec_(kk);
            Acc_(kk) = (TP_+TN_)/(TP_+FP_+TN_+FN_);
            [~,phighvslow(kk)]=ttest2(Yvariable(PredTsvm==1),Yvariable(PredTsvm==0));
            try
                [phighvslow_(kk)]=ranksum(Yvariable(PredTsvm==1),Yvariable(PredTsvm==0));
            catch me
                phighvslow_(kk)=NaN;
            end
            absdeltamean(kk) = abs(nanmean(Yvariable(PredTsvm==1))-nanmean(Yvariable(PredTsvm==0)));
            absdeltamedian(kk) = abs(nanmedian(Yvariable(PredTsvm==1))-nanmedian(Yvariable(PredTsvm==0)));
        end
        tempcost=eval(coststr); %tempp, 1-tempAUC, -tempmaxsensspec
        if sum(~isnan(tempcost))==0 %because min_index of [NaN NaN] is not NaN...
            break
        end
        [~,temp3i]=min(tempcost);
        if jj<=length(forcedvars)
            temp3i=forcedvars(jj);
        end
        CostList=[CostList tempcost(temp3i)];
        if jj>1&&jj>length(forcedvars)&&(CostList(end-1)-CostList(end))<=breakifdeltacostlessthan %min arbitrary 1% change in AUC
            break
        end
        Ilist=[Ilist temp3i];
    end
    if ~isempty(Exclude)
        IlistX{Exclude}=Ilist;
    else
        IlistX=Ilist;
    end
    
    if 1
        NN=length(Yvariable);
        table1=zeros(size(Amatrix,2)); %rows are variable number, cols are iterations, values are fraction of selections
        for i=1:length(IlistX)
            if ~isempty(Exclude)
                temp=IlistX{i};
            else
                temp=IlistX(i);
            end
            for j=1:length(temp)
                table1(temp(j),j)=table1(temp(j),j)+1/NN;
            end
        end
    end
    
    % now run model [still in leave1out]
    SVMModel = svmtrain(Amatrix(:,Ilist),criteriaR(:),'method',svpmethodstr,'kernel_function',kernelfunction,'rbf_sigma',sigma);
    PredT(Exclude) = svmclassify(SVMModel,Amatrix(Exclude,Ilist));
end


TP_ = nansum(criteriaR_(rangevalidate).*PredT(rangevalidate));
FP_ = nansum((1-criteriaR_(rangevalidate)).*PredT(rangevalidate));
TN_ = nansum((1-criteriaR_(rangevalidate)).*(1-PredT(rangevalidate)));
FN_ = nansum(criteriaR_(rangevalidate).*(1-PredT(rangevalidate)));
PPV_ = TP_/(TP_+FP_);
NPV_ = TN_/(TN_+FN_);
Sens_ = TP_/(TP_+FN_);
Spec_ = TN_/(TN_+FP_);
Acc_ = (TP_+TN_)/(TP_+FP_+TN_+FN_);
PPVNPVSensSpec = [PPV_ NPV_ Sens_ Spec_ Acc_ TP_ FP_ TN_ FN_];
[~,phighvslow]=ttest2(Yvariable(PredT==1),Yvariable(PredT==0));
[phighvslow_]=ranksum(Yvariable(PredT==1),Yvariable(PredT==0));

YvariableByPredictedOutcome=...
    [nanmean(Yvariable(PredT==1)) ...
    nanstd(Yvariable(PredT==1))/nansum(PredT==1)^0.5 ...
    nanmean(Yvariable(PredT==0)) ...
    nanstd(Yvariable(PredT==0))/nansum(PredT==0)^0.5 ...
    phighvslow nanmedian(Yvariable(PredT==1)) nanmedian(Yvariable(PredT==0)) phighvslow_];

%rerun without leave1out
Ilist=[]; AUClist=[]; plist=[]; CostList=[];
criteriaR=criteriaR_;
criteriaR(excludefortraining)=NaN; %leaveoutalways
    if length(maxNfeaturesJ)>1
        maxNfeatures=ceil(median(maxNfeaturesJ));
    end
    for jj=1:maxNfeatures %fwd stepwise %copy-paste from above if edits are made
        PPV_=NaN*ones(1,size(Amatrix,2));
        NPV_=NaN*ones(1,size(Amatrix,2));
        Sens_=NaN*ones(1,size(Amatrix,2));
        SensplusSpec_=NaN*ones(1,size(Amatrix,2));
        Acc_=NaN*ones(1,size(Amatrix,2));
        phighvslow=NaN*ones(1,size(Amatrix,2));
        phighvslow_=NaN*ones(1,size(Amatrix,2));
        absdeltamean=NaN*ones(1,size(Amatrix,2));
        absdeltamedian=NaN*ones(1,size(Amatrix,2));
        for kk=1:size(Amatrix,2) %run through each variable/column
            if sum(kk==Ilist)||sum(kk==tempignorevars)
                PPV_(kk)=NaN;
                NPV_(kk)=NaN;
                Sens_(kk)=NaN;
                SensplusSpec_(kk)=NaN;
                Acc_(kk)=NaN;
                phighvslow(kk)=NaN;
                phighvslow_(kk)=NaN;
                absdeltamean(kk)=NaN;
                absdeltamedian(kk)=NaN;
                continue
            end
            Ilisttest = [Ilist kk];
            SVMModel = svmtrain(Amatrix(:,Ilisttest),criteriaR(:),'method',svpmethodstr,'kernel_function',kernelfunction,'rbf_sigma',sigma);
            PredTsvm = svmclassify(SVMModel,Amatrix(:,Ilisttest));
            
            TP_ = nansum(criteriaR(:).*PredTsvm);
            FP_ = nansum((1-criteriaR(:)).*PredTsvm);
            TN_ = nansum((1-criteriaR(:)).*(1-PredTsvm));
            FN_ = nansum(criteriaR(:).*(1-PredTsvm));
            PPV_(kk) = TP_/(TP_+FP_);
            NPV_(kk) = TN_/(TN_+FN_);
            Sens_(kk) = TP_/(TP_+FN_);
            Spec_(kk) = TN_/(TN_+FP_);
            SensplusSpec_(kk) = Sens_(kk)+ Spec_(kk);
            Acc_(kk) = (TP_+TN_)/(TP_+FP_+TN_+FN_);
            [~,phighvslow(kk)]=ttest2(Yvariable(PredTsvm==1),Yvariable(PredTsvm==0));
            try
                [phighvslow_(kk)]=ranksum(Yvariable(PredTsvm==1),Yvariable(PredTsvm==0));
            catch me
                phighvslow_(kk)=NaN;
            end
            absdeltamean(kk) = abs(nanmean(Yvariable(PredTsvm==1))-nanmean(Yvariable(PredTsvm==0)));
            absdeltamedian(kk) = abs(nanmedian(Yvariable(PredTsvm==1))-nanmedian(Yvariable(PredTsvm==0)));
        end
        tempcost=eval(coststr); %tempp, 1-tempAUC, -tempmaxsensspec
        if sum(~isnan(tempcost))==0 %because min_index of [NaN NaN] is not NaN...
            break
        end
        [~,temp3i]=min(tempcost);
        if jj<=length(forcedvars)
            temp3i=forcedvars(jj);
        end
        CostList=[CostList tempcost(temp3i)];
        if jj>1&&jj>length(forcedvars)&&(CostList(end-1)-CostList(end))<=breakifdeltacostlessthan %min arbitrary 1% change in AUC
            break
        end
        Ilist=[Ilist temp3i];
    end

    % now run model
    SVMModel = svmtrain(Amatrix(:,Ilist),criteriaR(:),'method',svpmethodstr,'kernel_function',kernelfunction,'rbf_sigma',sigma);
    
    PredTAll = NaN*criteriaR;
    PredTAll(rangevalidate) = svmclassify(SVMModel,Amatrix(rangevalidate,Ilist));
    
    %need to set PredTAll
    
TP_ = nansum(criteriaR_(rangevalidate).*PredTAll(rangevalidate));
FP_ = nansum((1-criteriaR_(rangevalidate)).*PredTAll(rangevalidate));
TN_ = nansum((1-criteriaR_(rangevalidate)).*(1-PredTAll(rangevalidate)));
FN_ = nansum(criteriaR_(rangevalidate).*(1-PredTAll(rangevalidate)));
PPV_ = TP_/(TP_+FP_);
NPV_ = TN_/(TN_+FN_);
Sens_ = TP_/(TP_+FN_);
Spec_ = TN_/(TN_+FP_);
Acc_ = (TP_+TN_)/(TP_+FP_+TN_+FN_);
PPVNPVSensSpecAll = [PPV_ NPV_ Sens_ Spec_ Acc_ TP_ FP_ TN_ FN_];

[~,phighvslow]=ttest2(Yvariable(PredTAll==1),Yvariable(PredTAll==0));
[phighvslow_]=ranksum(Yvariable(PredTAll==1),Yvariable(PredTAll==0));

YvariableByPredictedOutcomeAll=...
    [nanmean(Yvariable(PredTAll==1)) ...
    nanstd(Yvariable(PredTAll==1))/nansum(PredTAll==1)^0.5 ...
    nanmean(Yvariable(PredTAll==0)) ...
    nanstd(Yvariable(PredTAll==0))/nansum(PredTAll==0)^0.5 ...
    phighvslow nanmedian(Yvariable(PredTAll==1)) nanmedian(Yvariable(PredTAll==0)) phighvslow_];
