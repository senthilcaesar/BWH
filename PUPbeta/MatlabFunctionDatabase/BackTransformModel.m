function T = BackTransformModel(mdl,F)
% Estimates backtransformed betas and 95%CIs by assuming zero values for other model covariates
% Models must use 1) mean-subtracted values if interpreting effect based on means of covariates
% Otherwise zero values should be meaningful i.e. 0 = non-smoking, reference body position etc  
%linkpower=0.5
%F.Link = @(mu) mu.^linkpower; F.Derivative= @(mu) linkpower*mu.^(linkpower-1); F.Inverse = @(mu) mu.^(1/linkpower); %link

[~,~,stats] = fixedEffects(mdl);
    ypredarray = [stats.Estimate stats.Lower stats.Upper];
    ypredarray_ = ypredarray + stats.Estimate(1);
    ypredarray__ = F.Inverse(ypredarray_) - F.Inverse(stats.Estimate(1));
    ypredarray__(1,:) = F.Inverse(ypredarray(1,:)); %overwrite constant with original
    ypred = ypredarray__(:,1);
    ypredCI = ypredarray__(:,2:3);
    Ydir = ypred<0;
    temp = [ypred ypredCI];
    temp(Ydir,:)=temp(Ydir,:)*-1;
    tempCI=temp(:,[2 3]);
    YSEMequiv = (temp(:,1) - min(tempCI')')/1.96;
    T = array2table([ypred YSEMequiv stats.pValue ypredCI]);
    T.Properties.RowNames = stats.Name;
    T.Properties.VariableNames = {'Estimate','SE','pValue','Lower','Upper'};