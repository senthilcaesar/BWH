function [ArousalScaleAvg]=TestArousalIntensity(FtrMtx,BaseLineFtrs,ClassificationModels)

% % % ClassificationModels.Mdl1 =fitcknn(X,y,'NumNeighbors',4,'Distance','cityblock');
% % % ClassificationModels.Mdl2 =fitcknn(X,y,'NumNeighbors',3,'Distance','cityblock');
% % % ClassificationModels.Mdl3 =fitcknn(X,y,'NumNeighbors',5,'Distance','cityblock');
% % % 
% % % ClassificationModels.Mdl4 = fitcdiscr(X,y,'DiscrimType','linear');
% % % ClassificationModels.Mdl5 = fitcdiscr(X,y,'DiscrimType','diaglinear');
% % % ClassificationModels.Mdl9 = fitcdiscr(X,y,'DiscrimType','quadratic');
% % % ClassificationModels.Mdl6 = fitcdiscr(X,y,'DiscrimType','diagquadratic');
% % % ClassificationModels.Mdl8 = fitcdiscr(X,y,'DiscrimType','pseudoquadratic');


DifferenceRatio=FtrMtx(:,1:23)./BaseLineFtrs(:,1:23);
DifferenceRatio(24)=FtrMtx(:,24);

FeaturesToBeRemoved=[zeros(1,4) ones(1,4) zeros(1,5) ones(1,1) zeros(1,5) ones(1,4) zeros(1,1)];
Difference=DifferenceRatio(:,~FeaturesToBeRemoved);


% ArousalScale(1) = predict(ClassificationModels.Mdl1,Difference);
% ArousalScale(2) = predict(ClassificationModels.Mdl2,Difference);
% ArousalScale(3) = predict(ClassificationModels.Mdl3,Difference);
% ArousalScale(4) = predict(ClassificationModels.Mdl4,Difference);
% ArousalScale(5) = predict(ClassificationModels.Mdl5,Difference);
% ArousalScale(9) = predict(ClassificationModels.Mdl9,Difference);
% ArousalScale(6) = predict(ClassificationModels.Mdl6,Difference);
% ArousalScale(8) = predict(ClassificationModels.Mdl8,Difference);
% ArousalScale(7) = predict(ClassificationModels.Mdl7,Difference);
% 
% W1=(1/9)+(1/9 - 1/18)/8;
% W2=1/18;
% ArousalScaleAvg = ArousalScale*[W1*ones(7,1);W2;W1];


% ArousalScale(1) = predict(ClassificationModels.Mdl2,Difference);
ArousalScale(1) = predict(ClassificationModels.Mdl3,Difference);
% ArousalScale(2) = predict(ClassificationModels.Mdl4,Difference);
ArousalScale(2) = predict(ClassificationModels.Mdl9,Difference);
ArousalScale(3) = predict(ClassificationModels.Mdl7,Difference);

W1=1/3;
% W2=1/18;
ArousalScaleAvg = ArousalScale*[W1*ones(3,1)];



