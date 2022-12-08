function [PPVNPVSensSpec,thresoptX,YvariableByPredictedOutcome,PredT,thresopt,PPVNPVSensSpecAll,YvariableByPredictedOutcomeAll,AUCdata]=univariateleave1out(X,criteriaR,funcsensspecstr,Yvariable)
%still editing this
if 0
    X = Tn;
    ythres=130;
    Yvariable=SBP;
    criteriaR = Yvariable>ythres;
end

PredT=NaN;
thresoptX=NaN;
pos=NaN;
for i=1:length(X)
    criteriaR_=1*criteriaR;
    criteriaR_(i)=NaN;
            [x,y,t,AUC,~] = perfcurve(criteriaR_*1,X,1); %need to find the threshold value that gives the OPTROCPT!
            if AUC<0.5,
                AUC=1-AUC;
                y=1-y;
                x=1-x;
                pos(i)=0;
            else
                pos(i)=1;
            end
            funcsensspec = eval(funcsensspecstr); %y+(1-x)%y+(1-x); % min distance 1./sqrt((1-y).^2+(x).^2) % harmonic mean y.*(1-x)./(y+(1-x))
            [~,I]=max(funcsensspec);
            thresoptX(i)=mean(t(I:(I+1)));
            PredT(i)=X(i)>thresoptX(i);
end
PredT=PredT(:);
%re-run with all data
[x,y,t,AUC,~] = perfcurve(criteriaR*1,X,1); %need to find the threshold value that gives the OPTROCPT!
if AUC<0.5,
    AUC=1-AUC;
    y=1-y;
    x=1-x;
    pos_all=0;
else
    pos_all=1;
end
if AUC<0.5, AUC=1-AUC; end
[~,AUCse,AUCp] = AUCci(AUC,sum(criteriaR(:)==1),sum(criteriaR(:)==0));
AUCdata = [AUC AUCse AUCp]; 

funcsensspec = eval(funcsensspecstr); %y+(1-x)%y+(1-x); % min distance 1./sqrt((1-y).^2+(x).^2) % harmonic mean y.*(1-x)./(y+(1-x))
[~,I]=max(funcsensspec);
thresopt=mean(t(I:(I+1)));
PredTAll=X>thresopt;
%temp = [PredTAll(:) PredT(:)]; % compare classifications

X_ = [nanmean(X(PredT==0)) nanmean(X(PredT==1))];
Y_ = [nanmean(Yvariable(PredT==0)) nanmean(Yvariable(PredT==1))];
X_all = [nanmean(X(PredTAll==0)) nanmean(X(PredTAll==1))];
Y_all = [nanmean(Yvariable(PredTAll==0)) nanmean(Yvariable(PredTAll==1))];
Errors = [nanstd(Yvariable(PredT==0)) nanstd(Yvariable(PredT==1))]./[nansum(PredT==0) nansum(PredT==1)].^0.5;
[~,phighvslow]=ttest2(Yvariable(PredT==0),Yvariable(PredT==1));   
[~,phighvslow_all]=ttest2(Yvariable(PredTAll==0),Yvariable(PredTAll==1));   

if pos_all==0
    PredTAll=1-PredTAll;
    PredT=1-PredT;
end

TP_ = nansum(criteriaR_(:).*PredT(:));
FP_ = nansum((1-criteriaR_(:)).*PredT(:));
TN_ = nansum((1-criteriaR_(:)).*(1-PredT(:)));
FN_ = nansum(criteriaR_(:).*(1-PredT(:)));
PPV_ = TP_/(TP_+FP_);
NPV_ = TN_/(TN_+FN_);
Sens_ = TP_/(TP_+FN_);
Spec_ = TN_/(TN_+FP_);
Acc_ = (TP_+TN_)/(TP_+FP_+TN_+FN_);
PPVNPVSensSpec = [PPV_ NPV_ Sens_ Spec_ Acc_ TP_ FP_ TN_ FN_];

TP_ = nansum(criteriaR(:).*PredTAll(:));
FP_ = nansum((1-criteriaR(:)).*PredTAll(:));
TN_ = nansum((1-criteriaR(:)).*(1-PredTAll(:)));
FN_ = nansum(criteriaR(:).*(1-PredTAll(:)));
PPV_ = TP_/(TP_+FP_);
NPV_ = TN_/(TN_+FN_);
Sens_ = TP_/(TP_+FN_);
Spec_ = TN_/(TN_+FP_);
Acc_ = (TP_+TN_)/(TP_+FP_+TN_+FN_);
PPVNPVSensSpecAll = [PPV_ NPV_ Sens_ Spec_ Acc_ TP_ FP_ TN_ FN_];

[~,phighvslow]=ttest2(Yvariable(PredT==1),Yvariable(PredT==0));

YvariableByPredictedOutcome=...
    [nanmean(Yvariable(PredT==1)) ...
    nanstd(Yvariable(PredT==1))/nansum(PredT==1)^0.5 ...
    nanmean(Yvariable(PredT==0)) ...
    nanstd(Yvariable(PredT==0))/nansum(PredT==0)^0.5 ...
    phighvslow];

[~,phighvslow]=ttest2(Yvariable(PredTAll==1),Yvariable(PredTAll==0));

YvariableByPredictedOutcomeAll=...
    [nanmean(Yvariable(PredTAll==1)) ...
    nanstd(Yvariable(PredTAll==1))/nansum(PredTAll==1)^0.5 ...
    nanmean(Yvariable(PredTAll==0)) ...
    nanstd(Yvariable(PredTAll==0))/nansum(PredTAll==0)^0.5 ...
    phighvslow];
