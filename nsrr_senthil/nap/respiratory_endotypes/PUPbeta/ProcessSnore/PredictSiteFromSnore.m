function [PredSnoreTblMLR,PredSnoreTblLR] = PredictSiteFromSnore(BreathSnoreTable)
global settings
%% Predict site of collapse
load([settings.folder 'Dropbox (Partners HealthCare)\MEEI DISE\Dan\Model workspaces\PredSiteMLR_MdlCoeff.mat'],...
    'B','FtrsInMdl','FtrNamesInMdl','FtrNames')
FtrIdx = ismember(BreathSnoreTable.Properties.VariableNames,FtrNames);
FeatureArray = BreathSnoreTable{:,FtrIdx};
[PrVOTEtemp] = mnrval(B,FeatureArray(:,FtrsInMdl), 'Model', 'nominal');
PrVOTE = [PrVOTEtemp(:,4) PrVOTEtemp(:,2) PrVOTEtemp(:,3) PrVOTEtemp(:,1)];
[~,PredVOTE_i] = max(PrVOTE,[],2);
VOTEopts = repmat([1 2 3 4], size(FeatureArray,1),1); 
PredVOTE = VOTEopts == PredVOTE_i;

PredSnoreMLR = [FeatureArray(:,FtrsInMdl) PrVOTE PredVOTE];
PredSnoreTblMLR = array2table(PredSnoreMLR);
PredSnoreTblMLR.Properties.VariableNames = [FtrNamesInMdl ...
    {'PrV_MLR' 'PrO_MLR' 'PrT_MLR' 'PrE_MLR'} ...
    {'PredV_MLR' 'PredO_MLR' 'PredT_MLR' 'PredE_MLR'}]; % 
clear B FtrNames FtrNamesInMdl FtrsInMdl PredSnoreMLR PredVOTE_i VOTEopts

% LR models
load([settings.folder 'Dropbox (Partners HealthCare)\MEEI DISE\Dan\Model workspaces\PredSiteLR_MdlCoeff.mat'],...
    'mdlFinalV','thresoptV','mdlFinalO','thresoptO','mdlFinalT','thresoptT','mdlFinalE','thresoptE')
PrV_LR = predict(mdlFinalV,BreathSnoreTable);
% PredV_LR = PrV_LR > thresoptV;
PrO_LR = predict(mdlFinalO,BreathSnoreTable);
% PredO_LR = PrO_LR > thresoptO;
PrT_LR = predict(mdlFinalT,BreathSnoreTable);
% PredT_LR = PrT_LR > thresoptT;
PrE_LR = predict(mdlFinalE,BreathSnoreTable);
% PredE_LR = PrE_LR > thresoptE;

% Mutually exclusive predictions
PrVOTE2 = [PrV_LR PrO_LR PrT_LR PrE_LR];
[~,PredVOTE_i2] = max(PrVOTE2,[],2);
VOTEopts = repmat([1 2 3 4], size(PrVOTE2,1),1); 
PredVOTE2 = VOTEopts == PredVOTE_i2;

FtrNamesInMdl = unique([mdlFinalV.CoefficientNames mdlFinalO.CoefficientNames ...
    mdlFinalT.CoefficientNames mdlFinalE.CoefficientNames]);
FtrIdx = ismember(BreathSnoreTable.Properties.VariableNames,...
    FtrNamesInMdl);
FtrTblSubset = BreathSnoreTable(:,FtrIdx);
PredSnoreLR = [FtrTblSubset{:,:} PrV_LR PrO_LR PrT_LR PrE_LR ...
    PredVOTE2];
PredSnoreTblLR = array2table(PredSnoreLR);
PredSnoreTblLR.Properties.VariableNames = [FtrTblSubset.Properties.VariableNames ...
    {'PrV_LR' 'PrO_LR' 'PrT_LR' 'PrE_LR'} ...
    {'PredV_LR' 'PredO_LR' 'PredT_LR' 'PredE_LR'}];

