function [ClassificationModels]=TrainArousalIntensity(TrainMatrix)

FeaturesToBeRemoved=[zeros(1,4) ones(1,4) zeros(1,5) ones(1,1) zeros(1,5) ones(1,4) zeros(1,1)];

X=TrainMatrix(:,~FeaturesToBeRemoved);
y=TrainMatrix(:,25);

ClassificationModels.Mdl1 =fitcknn(X,y,'NumNeighbors',4,'Distance','cityblock');
ClassificationModels.Mdl2 =fitcknn(X,y,'NumNeighbors',3,'Distance','cityblock');
ClassificationModels.Mdl3 =fitcknn(X,y,'NumNeighbors',5,'Distance','cityblock');

ClassificationModels.Mdl4 = fitcdiscr(X,y,'DiscrimType','linear');
ClassificationModels.Mdl5 = fitcdiscr(X,y,'DiscrimType','diaglinear');
ClassificationModels.Mdl9 = fitcdiscr(X,y,'DiscrimType','quadratic');
ClassificationModels.Mdl6 = fitcdiscr(X,y,'DiscrimType','diagquadratic');
ClassificationModels.Mdl8 = fitcdiscr(X,y,'DiscrimType','pseudoquadratic');

ClassificationModels.Mdl7ToBePruned= fitctree(X,y);
ClassificationModels.Mdl7 = prune(ClassificationModels.Mdl7ToBePruned,'level',6);







