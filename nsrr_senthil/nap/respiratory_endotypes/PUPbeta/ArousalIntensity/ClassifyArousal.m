function [ArousalScaleAvg]=ClassifyArousal(FtrMtx,BaseLineFtrs,TrainMatrix)
%[ArousalScaleAvg]=ClassifyArousal(FtrMtx,BaseLineFtrs,TrainMatrix)
%FtrMtx: wavelet feature matrix for the eeg signal during arousal, 
%BaseLineFtrs: wavelet feature matrix for the eeg signal before arousal (same length)
%TrainMatrix: Training data set with 16 columns and 542 observations. Last
%column reperesents the visual scores. 
%ArousalScaleAvg: Arousal score obtained by averaging the result of 8
%different classifiers


DifferenceRatio=FtrMtx(:,1:23)./BaseLineFtrs(:,1:23);
DifferenceRatio(24)=FtrMtx(:,24);

FeaturesToBeRemoved=[zeros(1,4) ones(1,4) zeros(1,5) ones(1,1) zeros(1,5) ones(1,4) zeros(1,1)];
Difference=DifferenceRatio(:,~FeaturesToBeRemoved);


X=TrainMatrix(:,~FeaturesToBeRemoved);
y=TrainMatrix(:,25);

ArousalScale(1) = knnclassify(Difference, X, y, 4,'cityblock');
ArousalScale(2) = knnclassify(Difference, X, y, 3,'cityblock');
ArousalScale(3) = knnclassify(Difference, X, y, 5,'cityblock');

ArousalScale(4) = classify(Difference,X,y,'linear');
ArousalScale(5) = classify(Difference,X,y,'diaglinear');%
ArousalScale(9) = classify(Difference,X,y,'quadratic');
ArousalScale(6) = classify(Difference,X,y,'diagquadratic');%
ArousalScale(8) = classify(Difference,X,y,'mahalanobis');

if 0 %out in 2018
t = classregtree(X,y,'method','classification');
t2 = prune(t,'level',6);
yfit = eval(t2,Difference);
ArousalScale(7)=str2double(yfit{1});
else %written by SS, AA to check
t_=fitctree(X,y);
t_ = prune(t_,'level',6);
ArousalScale(7)=predict(t_,Difference);
end

W1=(1/9)+(1/9 - 1/18)/8;
W2=1/18;
ArousalScaleAvg = round(ArousalScale*[W1*ones(7,1);W2;W1]);








