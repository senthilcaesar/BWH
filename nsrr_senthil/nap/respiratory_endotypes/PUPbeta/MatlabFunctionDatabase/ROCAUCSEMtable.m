function [T,TOptimal,TOptimalI,AUC,direction]=ROCAUCSEMtable(labels,x,outcomecontinuous)

%labels = x<thres;
I = isnan(labels) | isnan(x);
labels(I)=[];
x(I)=[];
outcomecontinuous(I)=[];

%test to find posclass
[~,~,~,AUC0,~] = perfcurve(labels,x,0);
[~,~,~,AUC1,~] = perfcurve(labels,x,1);
posclass = AUC1>AUC0;   
direction = posclass*2 - 1;
%rerun AUC
[X,Y,thr,AUC.value,OPTTHRES] = perfcurve(labels,x,posclass);
AUC2 = 0.5*sum( (X(2:end)-X(1:end-1)).*(Y(2:end)+Y(1:end-1))); %equation used
%%
%X=1-spec; Y=sens
sortedX = flipud(sort(unique(x(:))));

% x=direction*x;

T=[];
epsilon = 0.0000000001*nanstd(x);
T.ThresU = [sortedX(1)+epsilon;sortedX];
T.ThresL = [sortedX;sortedX(end)-epsilon];
T.Thres = mean([T.ThresU T.ThresL],2);

for i=1:length(T.Thres)
    if direction==1
    PredTrue = x>=T.Thres(i);
    PredFalse = x<T.Thres(i);
    else
    PredTrue = x<=T.Thres(i);
    PredFalse = x>T.Thres(i);
    end
    TruePos = PredTrue & labels==1;
    TrueNeg = PredFalse & labels~=1;
    FalsePos = PredTrue & labels~=1;
    FalseNeg = PredFalse & labels==1;
    T.TP(i,1) = sum(TruePos);
    T.TN(i,1) = sum(TrueNeg);
    T.FP(i,1) = sum(FalsePos);
    T.FN(i,1) = sum(FalseNeg);  
    T.YPredTrueMedian(i,1) = median(outcomecontinuous(PredTrue));
    try
    T.YPredTrueMedianP(i,1) = signrank(outcomecontinuous(PredTrue));
    catch me
    T.YPredTrueMedianP(i,1) = NaN;
    end
    T.YPredTrueMean(i,1) = mean(outcomecontinuous(PredTrue));
    [~,T.YPredTrueMeanP(i,1)] = ttest(outcomecontinuous(PredTrue),zeros(length(outcomecontinuous(PredTrue)),1));
end
T=struct2table(T);
T.Se = T.TP./(T.TP+T.FN);
T.Sp = T.TN./(T.TN+T.FP);
T.X = 1-T.Sp;
T.PPV = T.TP./(T.TP+T.FP);
T.NPV = T.TN./(T.TN+T.FN);
T.SePlusSp = T.Se + T.Sp;
T.FrPredPos = (T.TP + T.FP)./(length(x));

%%
% T.x_ = X;
% T.y_ = Y;


[~,TOptimalI]=max(T.SePlusSp);

TOptimal = T(TOptimalI,:);
%find best Xthreshold and Sens/Spec

%Calculate SEM and p
N1 = sum(labels==1); %must be the positive, i.e. abnormal, group being detected
N2 = sum(labels==0);
Q1 = AUC.value/(2-AUC.value);
Q2 = 2*AUC.value^2/(1+AUC.value);
AUC.SEM = ((AUC.value*(1-AUC.value)+(N1-1)*(Q1-AUC.value^2)+(N2-1)*(Q2-AUC.value^2))/(N1*N2))^0.5;
%from Hanley McNeill 1982 
%accessed free at http://www.med.mcgill.ca/epidemiology/hanley/software/Hanley_McNeil_Radiology_82.pdf

AUC.p = 2*[1-normcdf(abs(AUC.value-0.5)./AUC.SEM,0,1)];