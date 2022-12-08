%% Simple compare of two .mats

[F_1 , P_1] = uigetfile('*.mat', 'Select the first .mat file)');
[F_2 , P_2] = uigetfile('*.mat', 'Select the second .mat file)', P_1);

visdiff(fullfile([P_1 F_1]), fullfile([P_2 F_2]))
