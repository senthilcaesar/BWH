function [PPVNPVSensSpec,table1,thresoptX,YvariableByPredictedOutcome,PredT,B,Bstats,Ilist,thresopt,PPVNPVSensSpecAll,YvariableByPredictedOutcomeAll,B1array,AUCdata]=logregleave1out(Amatrix,criteriaR,tempignorevars,funcsensspecstr,Yvariable,maxNfeaturesJ,breakifdeltacostlessthan,plotfigs,rangevalidate,excludefortraining,forcedvars)

%Written by Scott Sands 2015-10-28

%funcsensspecstr='y+(1-x)';
%breakifdeltacostlessthan=0.001;

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
thresoptX = NaN*(1:length(criteriaR))';
for Exclude=rangevalidate
    Ilist=[]; AUClist=[]; plist=[]; CostList=[];
    criteriaR=criteriaR_;
    criteriaR(Exclude)=NaN; %leave1out
    criteriaR(excludefortraining)=NaN; %leaveoutalways
    if length(maxNfeaturesJ)>1
        maxNfeatures=maxNfeaturesJ(Exclude);
    end
    for jj=1:maxNfeatures %fwd stepwise
        tempAUC=NaN*ones(1,size(Amatrix,2));
        tempmaxsensspec=NaN*ones(1,size(Amatrix,2));
        tempp=NaN*ones(1,size(Amatrix,2));
        for kk=1:size(Amatrix,2) %run through each variable/column
            if sum(kk==Ilist)||sum(kk==tempignorevars)
                tempp(kk)=NaN;
                tempAUC(kk)=NaN;
                tempmaxsensspec(kk)=NaN;
                continue
            end
            Ilisttest = [Ilist kk];
            [Btemp,~,temp] = glmfit(Amatrix(:,Ilisttest),criteriaR(:),'binomial'); %,'weights',weights
            tempp(kk)=temp.p(end);
            Amatrixtemp=Amatrix(:,Ilisttest);
            p_success=(glmval(Btemp,Amatrix(:,Ilisttest),'logit'));
            DeltaAHIp_logit = log(p_success./(1-p_success)); %logit
            X=p_success;
            [x,y,~,AUC,~] = perfcurve(criteriaR*1,X,1); %need to find the threshold value that gives the OPTROCPT!
            if AUC<0.5,
                AUC=1-AUC;
                y=1-y;
                x=1-x;
            end
            tempAUC(kk)=AUC; %AUC
            tempmaxsensspec(kk)=max(eval(funcsensspecstr));
        end
        tempcost=1-tempAUC; %tempp, 1-tempAUC, -tempmaxsensspec
        if sum(~isnan(tempcost))==0 %because min_index of [NaN NaN] is not NaN...
            break
        end
        [~,temp3i]=min(tempcost);
        if jj<=length(forcedvars)
            temp3i=forcedvars(jj);
        end
        plist=[plist tempp(temp3i)];
        AUClist=[AUClist tempAUC(temp3i)];
        CostList=[CostList tempcost(temp3i)];
        if jj>1&&jj>length(forcedvars)&&(CostList(end-1)-CostList(end))<=breakifdeltacostlessthan %min arbitrary 1% change in AUC
            break
        end
        Ilist=[Ilist temp3i];
    end
    if ~isempty(Exclude)
        IlistX{Exclude}=Ilist;
        AUClistX{Exclude}=AUClist;
    else
        IlistX=Ilist;
        AUClistX=AUClist;
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
    [B,~,~] = glmfit(Amatrix(:,Ilist),criteriaR(:),'binomial'); %,'weights',weights
    B1array(Exclude)=B(2);
    p_success=(glmval(B,Amatrix(:,Ilist),'logit'));
    DeltaAHIp_logit = log(p_success./(1-p_success)); %logit
    X=p_success;
    [x,y,t,AUC,opt] = perfcurve(criteriaR*1,X,1); %need to find the threshold value that gives the OPTROCPT!
    
    %find optimal threshold from ROC.
    funcsensspec = eval(funcsensspecstr); %y+(1-x)%y+(1-x); % min distance 1./sqrt((1-y).^2+(x).^2) % harmonic mean y.*(1-x)./(y+(1-x))
    [~,I]=max(funcsensspec);
    
    thresopt=mean(t(I:(I+1)));
    thresoptX(Exclude)=thresopt;
   
    
    %t_at_highsens=[t(find(y==1,1)) t(find(y>0.8,1))];
    if plotfigs
        figure(105); set(gcf,'color',[1 1 1]);
        plot(p_success(criteriaR==0),Yvariable(criteriaR==0),'r.','markersize',15); box('off'); hold('on');
        plot(p_success(criteriaR==1),Yvariable(criteriaR==1),'.','color',[0 0.5 0],'markersize',15);
        plot(p_success(Exclude),Yvariable(Exclude),'k.','markersize',30); box('off'); hold('on');
        plot(thresopt*[1 1],[min(Yvariable) max(Yvariable)],'k:');
        ylim([min(Yvariable) max(Yvariable)]);
        xlim([0 1]);
        hold('off');
    end
    %were we right?
    if ~isempty(Exclude)
        PredT(Exclude)=p_success(Exclude)>=thresopt;
        if isnan(p_success(Exclude))
            PredT(Exclude)=NaN;
        end
    else
        PredT=p_success>=thresopt;
    end
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
for jj=1:maxNfeatures %fwd stepwise
    tempAUC=NaN*ones(1,size(Amatrix,2));
    tempmaxsensspec=NaN*ones(1,size(Amatrix,2));
    tempp=NaN*ones(1,size(Amatrix,2));
    for kk=1:size(Amatrix,2)
        if sum(kk==Ilist)||sum(kk==tempignorevars)
            tempp(kk)=NaN;
            tempAUC(kk)=NaN;
            tempmaxsensspec(kk)=NaN;
            continue
        end
        Ilisttest = [Ilist kk];
        [Btemp,~,temp] = glmfit(Amatrix(:,Ilisttest),criteriaR(:),'binomial'); %,'weights',weights
        tempp(kk)=temp.p(end);
        p_success=(glmval(Btemp,Amatrix(:,Ilisttest),'logit'));
        DeltaAHIp_logit = log(p_success./(1-p_success)); %logit
        X=p_success;
        [x,y,~,AUC,~] = perfcurve(criteriaR*1,X,1); %need to find the threshold value that gives the OPTROCPT!
        % sensspecweighted=max(2*y/3+(1-x)/3);
        if AUC<0.5,
            AUC=1-AUC;
            y=1-y;
            x=1-x;
        end
        tempAUC(kk)=AUC; %AUC
        tempmaxsensspec(kk)=max(eval(funcsensspecstr));
    end
    tempcost=1-tempAUC; %tempp, 1-tempAUC, -tempmaxsensspec
    if sum(~isnan(tempcost))==0 %because min_index of [NaN NaN] is not NaN...
        break
    end
    [~,temp3i]=min(tempcost);
    if jj<=length(forcedvars)
            temp3i=forcedvars(jj);
    end
    plist=[plist tempp(temp3i)];
    AUClist=[AUClist tempAUC(temp3i)];
    CostList=[CostList tempcost(temp3i)];
    if jj>1&&jj>length(forcedvars)&&(CostList(end-1)-CostList(end))<=breakifdeltacostlessthan %min arbitrary 1% change in AUC
        break
    end
    Ilist=[Ilist temp3i];
end

% now run model
[B,~,Bstats] = glmfit(Amatrix(:,Ilist),criteriaR(:),'binomial'); %,'weights',weights
p_success=(glmval(B,Amatrix(:,Ilist),'logit'));
DeltaAHIp_logit = log(p_success./(1-p_success)); %logit
X=p_success;

[x,y,t,~,opt] = perfcurve(criteriaR*1,X,1); %need to find the threshold value that gives the OPTROCPT!

%only use the testing data for output AUC
[~,~,~,AUC] = perfcurve(criteriaR(rangevalidate)*1,X(rangevalidate),1); %need to find the threshold value that gives the OPTROCPT!
if AUC<0.5, AUC=1-AUC; end
[~,AUCse,AUCp] = AUCci(AUC,sum(criteriaR(rangevalidate)==1),sum(criteriaR(rangevalidate)==0));
AUCdata = [AUC AUCse AUCp]; 

%find optimal threshold from ROC.
funcsensspec = eval(funcsensspecstr); %y+(1-x)%y+(1-x); % min distance 1./sqrt((1-y).^2+(x).^2) % harmonic mean y.*(1-x)./(y+(1-x))
[~,I]=max(funcsensspec);

thresopt=mean(t(I:(I+1)));
%thresoptX(Exclude)=thresopt;

%t_at_highsens=[t(find(y==1,1)) t(find(y>0.8,1))];
if plotfigs
    figure(105); set(gcf,'color',[1 1 1]);
    plot(p_success(criteriaR==0),Yvariable(criteriaR==0),'r.','markersize',15); box('off'); hold('on');
    plot(p_success(criteriaR==1),Yvariable(criteriaR==1),'.','color',[0 0.5 0],'markersize',15);
    %plot(DeltaAHIp_logit(Exclude),Yvariable(Exclude),'k.','markersize',30); box('off'); hold('on');
    plot(thresopt*[1 1],[min(Yvariable) max(Yvariable)],'k:');
    ylim([min(Yvariable) max(Yvariable)]);
    xlim([0 1]);
    hold('off');
end

%were we right?
    PredTAll = NaN*criteriaR;
    PredTAll(rangevalidate) = p_success(rangevalidate)>=thresopt;
    

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
