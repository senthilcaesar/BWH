function [dL_opt, dU_opt] = RejectClassExploreFn(Pr, LabelsActual)

%% Asymmetrical d
% d is a parameter that decides the degree to which an observation should 
% be rejected. 
% d = 0.5 - reject none
% d = 0 = reject all

% dL = [0:0.01:0.5]';
dL = [0.05:0.005:0.5]';
dU = 1-dL;

% Rejection rate is the number of rejections (i.e. observations that end up 
% between d and 1-d), normalized to number of observations (i.e. participants)
RejectRate = nan(length(dL)*length(dU),1);

% Normalized sum of the differences in probability and actual label
ErrorRate = nan(length(dL)*length(dU),1);

% Probability output from LR
% Pr = mdlAll.Fitted.Probability;

% E minus R ratio = (ErrorRate(0) - ErrorRate(R))/RejectionRate
ERRatio = nan(length(dL)*length(dU),1);
% Error_0 = 1-ComputeAcc(logical(round(Pr)), LabelsActual); % classification loss
Error_0 = sum(abs(Pr - LabelsActual))/ length(Pr); % Mean absolute error
% Error_0 = sum((Pr - LabelsActual).^2)/ length(Pr); % Mean absolute error
% Error_0 = -(LabelsActual'*log(Pr) + ... % log loss
%             (1-LabelsActual)'*log(1-Pr)) / length(Pr);

dLstore = nan(length(dL)*length(dU),1);
dUstore = nan(length(dL)*length(dU),1);

% initialize counter
count = 0;

for ii = 1:length(dL)
    for jj = 1:length(dU)
        count = count+1;
        dLval = dL(ii);
        dUval = dU(jj);
        
        reject = Pr > dLval & Pr < dUval;
        responder = Pr > dUval;
        nonresponder = Pr < dLval;
      
        % Error rate
        LabelsPred_Rej = false(size(LabelsActual));
        LabelsPred_Rej(responder) = 1; 
        LabelsPred_Rej(reject) = [];
        NonRejIdx = logical(responder+nonresponder);
        LabelsActual_Rej = LabelsActual(NonRejIdx);
        % Classification error/loss
%         ErrorRate(count) = 1-ComputeAcc(LabelsPred_Rej, LabelsActual_Rej);
        
        % Mean absolute error
        ErrorRate(count) = (sum(abs(LabelsActual_Rej - Pr(NonRejIdx))))...
            / sum(NonRejIdx);     
%         
%         % Mean squared error
%         ErrorRate(count) = (sum((LabelsActual_Rej - Pr(NonRejIdx)).^2))...
%             / sum(NonRejIdx); 
%         
%         % Log loss
%         ErrorRate(count) = -(LabelsActual_Rej'*log(Pr(NonRejIdx)) + ...
%             (1-LabelsActual_Rej)'*log(1-Pr(NonRejIdx))) / sum(NonRejIdx);
%         
        
        % Reject rate
        RejectRate(count) = sum(reject) / length(Pr);

        ERRatio(count) = (Error_0 - ErrorRate(count))/RejectRate(count);
        
        % store dL and dU so I can find them
        dLstore(count) = dLval;
        dUstore(count) = dUval;
    end
end

ERRatio(isinf(ERRatio)) = 0;
ERRatio(isnan(ERRatio)) = 0;

[dopt, dopt_i] = max(ERRatio);
maxIdx = find(ERRatio == dopt);
RejectRate(maxIdx);

dL_opt = dLstore(dopt_i);
dU_opt = dUstore(dopt_i);

% dL_opt = 0.33;
% dU_opt = 0.65;
% reject = Pr > dL_opt & Pr < dU_opt;
% responder = Pr > dU_opt;
% nonresponder = Pr < dL_opt;
% 
% % Compute stats
% LabelsPred_Rej = false(size(LabelsActual));
% LabelsPred_Rej(responder) = 1; 
% LabelsPred_Rej(reject) = [];
% NonRejIdx = logical(responder+nonresponder);
% LabelsActual_Rej = LabelsActual(NonRejIdx);
% PercentReduction_Rej = PercentReduction(NonRejIdx);
% PercentReductionT_Rej = OutcomesTableFinal.PercentReduction(NonRejIdx);
% BaselineAHI_Rej = BaselineTableFinal.TotalAHI_Baseline(NonRejIdx);
% TreatmentAHI_Rej = BaselineTableFinal.TotalAHI_Treatment(NonRejIdx);
% 
% [statsTableAll_Rej, TrueFalseTableAll_Rej] = ComputeStats(LabelsPred_Rej, LabelsActual_Rej,...
%     PercentReduction_Rej, PercentReductionT_Rej, BaselineAHI_Rej, TreatmentAHI_Rej);

figure, plot(ErrorRate), hold on, plot(RejectRate), hold off
figure, scatter(RejectRate, ERRatio, 'o', 'MarkerFaceColor', [0.4 0.4 0.4],...
    'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'none')
xlabel('Rejection Rate')
ylabel('\epsilonR')
set(gca, 'FontSize', 15)


if 0
% Plot
% Pr vs. Percent reduction
figure('Position', [512 328 660 529])
scatter(Pr(LabelsActual), PercentReduction(LabelsActual), 75, 'o', ...
    'MarkerEdgeColor', [0 0 0], ...
    'LineWidth',1.5); hold on;
scatter(Pr(~LabelsActual), PercentReduction(~LabelsActual), 75, 'x', ...
    'MarkerEdgeColor', [0 0 0], 'LineWidth', 2);
plot([dL_opt dL_opt], [-200 100], '--', 'Color', [0.7 0.3 0.3])
plot([dU_opt dU_opt], [-200 100], '--', 'Color', [0.7 0.3 0.3])
xlabel('Probability')
ylabel('Percent Reduction')
set(gca,'FontSize',14)
end