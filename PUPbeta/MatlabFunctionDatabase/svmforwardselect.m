function [Ilist,PPVNPVSensSpecAll,YvariableByPredictedOutcomeAll,PredTAll,tempXX,tempXX2]=svmforwardselect(Amatrix,criteriaR,tempignorevars,Yvariable,maxNfeaturesJ,breakifdeltacostlessthan,rangevalidate,excludefortraining,forcedvars,coststr,svpmethodstr,kernelfunction,sigma)

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
criteriaR_ = criteriaR(:); %keep for later

%rerun without leave1out
Ilist=[]; CostList=[]; tempXX=[]; tempXX2=[];
criteriaR=criteriaR_;
criteriaR(excludefortraining)=NaN; %leave out for training
    if length(maxNfeaturesJ)>1
        maxNfeatures=ceil(median(maxNfeaturesJ));
    end
    for jj=1:maxNfeatures %fwd stepwise %copy-paste from above if edits are made
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
        tic
        parfor kk=1:size(Amatrix,2) %run through each variable/column
            if sum(kk==Ilist)||sum(kk==tempignorevars)
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
            Ilisttest = [Ilist kk];
            
            SVMModel = svmtrain(Amatrix(:,Ilisttest),criteriaR(:),'method',svpmethodstr,'kernel_function',kernelfunction,'rbf_sigma',sigma);
            PredTsvm = svmclassify(SVMModel,Amatrix(:,Ilisttest)); %all
            PredTsvm2 = svmclassify(SVMModel,Amatrix(rangevalidate,Ilisttest)); %test only
            performance{kk} = PredictiveValue(criteriaR_,PredTsvm,Yvariable);
            TP_(kk) = nansum(criteriaR(:).*PredTsvm);
            FP_(kk) = nansum((1-criteriaR(:)).*(PredTsvm));
            TN_(kk) = nansum((1-criteriaR(:)).*(1-PredTsvm));
            FN_(kk) = nansum(criteriaR(:).*(1-PredTsvm));
            PPV_(kk) = TP_(kk)/(TP_(kk)+FP_(kk));
            NPV_(kk) = TN_(kk)/(TN_(kk)+FN_(kk));
            Sens_(kk) = TP_(kk)/(TP_(kk)+FN_(kk));
            Spec_(kk) = TN_(kk)/(TN_(kk)+FP_(kk));
            SensplusSpec_(kk) = Sens_(kk)+ Spec_(kk);
            Acc_(kk) = (TP_(kk)+TN_(kk))/(TP_(kk)+FP_(kk)+TN_(kk)+FN_(kk));
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
        end
        toc
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
        tempXX=[tempXX; PPV_(temp3i) NPV_(temp3i) Sens_(temp3i) Spec_(temp3i) Acc_(temp3i) TP_(temp3i) FP_(temp3i) TN_(temp3i) FN_(temp3i) phighvslow(temp3i) phighvslow_(temp3i) absdeltamean(temp3i) absdeltamedian(temp3i)];
        tempXX2=[tempXX2; ... 
            performance{temp3i}.PPV_sem_chance_p(1)... 
            performance{temp3i}.NPV_sem_chance_p(1)... 
            performance{temp3i}.Sens_sem_chance_p(1)... 
            performance{temp3i}.Spec_sem_chance_p(1)... 
            performance{temp3i}.Acc_sem_chance_p(1)... 
            performance{temp3i}.TP_FP_TN_FN... 
            performance{temp3i}.phighvslow_ttest... 
            performance{temp3i}.phighvslow_ranksum... 
            performance{temp3i}.predY_mean-performance{temp3i}.predN_mean... 
            performance{temp3i}.predY_median-performance{temp3i}.predN_median]
            
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