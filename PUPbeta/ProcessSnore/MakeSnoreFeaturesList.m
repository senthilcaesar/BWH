function [FeatureNames]=MakeSnoreFeaturesList(BBs)
%% Populate the feature names variable
% this is a list of all the attribute curently returned
% this function also creates variables for each item in this list, so if 
% you add to this list, a variable of added name will be created back in
% the caler function, and initialised as nan(BBs,1);


FeatureNames = {...
    'SnoreDBMean_i','SnoreDBMean_e',...
    'WinVarMean_i', 'WinVarMax_i', 'WinVarMean_e', 'WinVarMax_e', ... 
    'p0to150Mean_i', 'p0to150Mean_e', 'p150to300Mean_i', 'p150to300Mean_e', ...
    'p300to450Mean_i', 'p300to350Mean_e', 'p450to600Mean_i', 'p450to600Mean_i', ...
    'p600to750Mean_i', 'p600to750Mean_e', 'p750to900Mean_i', 'p750to900Mean_e', ...
    'p0to1kMean_i', 'p0to1kMean_e', 'p1kto2kMean_i', 'p1kto2kMean_e', ...
    'ptotMean_i', 'ptotMean_e','centroidMean_i', 'centroidMean_e',...
    'formant1Mean_i', 'formant1Mean_e', 'formant2Mean_i', 'formant2Mean_e',...
    'formant3Mean_i', 'formant3Mean_e','f0Mean_i','f0Mean_e',...
    'SnoreDBVar_i','SnoreDBVar_e',...
    'WinVarVar_i', 'WinVarMax_i', 'WinVarVar_e', 'WinVarMax_e', ... 
    'p0to150Var_i', 'p0to150Var_e', 'p150to300Var_i', 'p150to300Var_e', ...
    'p300to450Var_i', 'p300to350Var_e', 'p450to600Var_i', 'p450to600Var_i', ...
    'p600to750Var_i', 'p600to750Var_e', 'p750to900Var_i', 'p750to900Var_e', ...
    'p0to1kVar_i', 'p0to1kVar_e', 'p1kto2kVar_i', 'p1kto2kVar_e', ...
    'ptotVar_i', 'ptotVar_e','centroidVar_i', 'centroidVar_e',...
    'formant1Var_i', 'formant1Var_e', 'formant2Var_i', 'formant2Var_e',...
    'formant3Var_i', 'formant3Var_e','f0Var_i','f0Var_e'...
};
% as a way of preventing possible human errors, the variables based upon
% these are automatically enumerated from the FeaturesNames list
for ft = 1:length(FeatureNames)
    % initialise a variable of name "FeaturesNames(ft)" as nan(BBs,1);
    str = char(FeatureNames(ft));
    assignin('caller', str, nan(BBs,1)); % make them back in the caller Fn
end

end
