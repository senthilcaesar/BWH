% testidx = index of leaveout subject
% FeatureTableFinal.Labels = Response labels
% OutputProbTr = Probability output from logistic regression from training data

% RejectClassExploreFn takes Probabilities and Labels from training set
% Finds optimal lower bound and upper bound
% Probabilities below lower bound (dL) are "non-responders" or 0 
% Probabilities above upper bound (dU) are "responders" or 1
% In between are rejected

[dL, dU] = RejectClassExploreFn(OutputProbTr, FeatureTableFinal.Labels(~testIdx)); 
if OutputProbTe > dL && OutputProbTe < dU
    RejectIdx(leaveout,1) = true;
elseif OutputProbTe > dU
    LabelsPred(leaveout,1) = true;
elseif OutputProbTe < dL
    LabelsPred(leaveout,1) = false;
end
dLstore(leaveout) = dL;
dUstore(leaveout) = dU;