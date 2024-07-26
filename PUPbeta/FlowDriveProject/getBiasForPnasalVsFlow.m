%% this function returns the bias between Pnasal and Flow ventilation
function [SmMed_p, SmMed_f, SmBias1, SmBias2, LgMed_p, LgMed_f, LgBias1, LgBias2] = ...
    getBiasForPnasalVsFlow(pnasal, flow)
% we calculate the median bias for small breaths (i.e. <70% of flow), large
% breaths (i.e. >130% of flow), and also the average which invcludes both
% small and large breaths.
% calculate bias, as 
% 1. pnasal to flow, (median value Pnasal) / (median value Pneumotach flow)
% 2. flow to pnasal, (median value flow) / (median value Pneumotach Pnasal)

% Find small pnasal breaths (use flow as the measure)
SmallBB_i = flow<0.7;
%NumSmallBB = nnz(SmallBB_i);
SmMed_f = median(flow(SmallBB_i));
SmMed_p = median(pnasal(SmallBB_i));
SmBias1 = SmMed_p ./ SmMed_f;
SmBias2 = SmMed_f ./ SmMed_p;

% Find large pnasal breaths (use flow as the measure)
LargeBB_i = flow>1.3;
%NumLargeBB = nnz(LargeBB_i);
LgMed_f = median(flow(LargeBB_i));
LgMed_p = median(pnasal(LargeBB_i));
LgBias1 = LgMed_p ./ LgMed_f;
LgBias2 = LgMed_f ./ LgMed_p;
