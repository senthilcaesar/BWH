function CoeffTable = CalculateBandORFn(DataTable,mdl,ixvar,varargin)

%% Model All
numSD = 2;
CoeffTable = mdl.Coefficients;
CoeffTable.Properties.VariableNames{1} = 'Beta';
CoeffTable.SE = [];
CoeffTable.B_ci = nan(size(mdl.Coefficients,1),2);
CoeffTable.OR = nan(size(mdl.Coefficients,1),1);
CoeffTable.OR_ci = nan(size(mdl.Coefficients,1),2);
CoeffTable = movevars(CoeffTable, 'OR', 'Before', 'tStat');
CoeffTable = movevars(CoeffTable, 'OR_ci', 'Before', 'tStat');
CoeffTable = movevars(CoeffTable, 'B_ci', 'Before', 'OR');

varnames = mdl.Coefficients.Properties.RowNames;

% get CI for intercept
CoeffTable.B_ci(1,1) = mdl.Coefficients.Estimate(1) - 1.96*mdl.Coefficients.SE(1);
CoeffTable.B_ci(1,2) = mdl.Coefficients.Estimate(1) + 1.96*mdl.Coefficients.SE(1);

if isempty(varargin) % Beta is in units of the outcome variable per numSD change in input variable
    OutcomeSD = 1;
else % Beta is in units of the SD change of outcome variable per numSD change in input variable
    outcome = varargin{1};
    OutcomeSD = nanstd(outcome);
end
    
%%
for idx = 1:size(mdl.Coefficients,1)-1
    
    varIdx = strcmp(DataTable.Properties.VariableNames, varnames{idx+1});
    FeatSD = nanstd(DataTable{:,varIdx})*numSD; % per 2SD - equiv to unit change in dichotomous variable
    if sum(varIdx) == 0 % check if its because var is binary (fitlm adds _1 to varname)
        varIdx = strcmp(DataTable.Properties.VariableNames, varnames{idx+1}(1:end-2));
    end

    VarUniq = unique(DataTable{:,varIdx});
    if length(find(~isnan(VarUniq))) == 2 && nansum(VarUniq) == 1 % i.e variable is dichotomous
        FeatSD = 1;
    end
       
    % Special case for AASM grant, can code this in later (might be tricky though)
    if contains(mdl.Coefficients.Properties.RowNames(idx+1,1),':')
        varIdx = strcmp(DataTable.Properties.VariableNames, ixvar);
        FeatSD = nanstd(DataTable{:,varIdx})*numSD; % per numSD
        
        BetaIdx = find(contains(mdl.Coefficients.Properties.RowNames, ixvar));
        Beta1 = mdl.Coefficients.Estimate(BetaIdx(1));
        Beta2 = mdl.Coefficients.Estimate(BetaIdx(2));
        Beta = (Beta1+Beta2)*(FeatSD/OutcomeSD);
        
        SEnew = sqrt(mdl.Coefficients.SE(BetaIdx(1))^2 +...
            mdl.Coefficients.SE(BetaIdx(2))^2 + 2*mdl.CoefficientCovariance(BetaIdx(1),BetaIdx(2)));
        
        BetaCI(1) = ((Beta1+Beta2) - 1.96*SEnew)*(FeatSD/OutcomeSD); 
        BetaCI(2) = ((Beta1+Beta2) + 1.96*SEnew)*(FeatSD/OutcomeSD); 
%         
%         BetaCI2_lower = Beta2 - 1.96*mdl.Coefficients.SE(BetaIdx(2)); 
%         BetaCI2_upper = Beta2 + 1.96*mdl.Coefficients.SE(BetaIdx(2)); 
%         
%         BetaCI(1) = (BetaCI1_lower + BetaCI2_lower)*(FeatSD/OutcomeSD);
%         BetaCI(2) = (BetaCI1_upper + BetaCI2_upper)*(FeatSD/OutcomeSD);
    else
        Beta = mdl.Coefficients.Estimate(idx+1)*(FeatSD/OutcomeSD);
        BetaCI(1) = (mdl.Coefficients.Estimate(idx+1) - 1.96*mdl.Coefficients.SE(idx+1))*(FeatSD/OutcomeSD);
        BetaCI(2) = (mdl.Coefficients.Estimate(idx+1) + 1.96*mdl.Coefficients.SE(idx+1))*(FeatSD/OutcomeSD);
    end

    OR = exp(Beta);
    ORci= exp(BetaCI);

%     if ORci(1) > ORci(2)
%         ORci = fliplr(ORci);
%     end

    CoeffTable.Beta(idx+1) = Beta;
    CoeffTable.B_ci(idx+1,:) = BetaCI;
    CoeffTable.OR(idx+1) = OR;
    CoeffTable.OR_ci(idx+1,:) = ORci;
end

%%
% if isempty(varargin)
%     for idx = 1:size(mdl.Coefficients,1)-1
%         varIdx = strcmp(DataTable.Properties.VariableNames, varnames{idx+1});
%         FeatSD = nanstd(DataTable{:,varIdx});
%         
%         Beta = mdl.Coefficients.Estimate(idx+1)*FeatSD;
%         BetaCI(1) = (mdl.Coefficients.Estimate(idx+1) - 1.96*mdl.Coefficients.SE(idx+1))*FeatSD;
%         BetaCI(2) = (mdl.Coefficients.Estimate(idx+1) + 1.96*mdl.Coefficients.SE(idx+1))*FeatSD;
% 
%         OR = exp(Beta);
%         ORci= exp(BetaCI);
% 
%         if ORci(1) > ORci(2)
%             ORci = fliplr(ORci);
%         end
% 
%         CoeffTable.Beta(idx+1) = Beta;
%         CoeffTable.B_ci(idx+1,:) = BetaCI;
%         CoeffTable.OR(idx+1) = OR;
%         CoeffTable.OR_ci(idx+1,:) = ORci;
%     end
% else