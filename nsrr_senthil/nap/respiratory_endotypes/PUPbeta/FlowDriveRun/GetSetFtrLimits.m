                                                         
%% Original model

close
clear
clc

%%
load('FinalModel_Table.mat')   
load('FtrNames165.mat')
Pnasal_lims = load('FeatureLimits_Pnasal.mat');
Pneumo_lims = load('FeatureLimits_Pneumo.mat');

%%
upper_flow = NaN(height(FinalModelTable),1);
lower_flow = NaN(height(FinalModelTable),1);
upper_pnasal = NaN(height(FinalModelTable),1);
lower_pnasal = NaN(height(FinalModelTable),1);

%%
for ii = 2:size(FinalModelTable)
    varname = FinalModelTable.Feature{ii};
    varname(length(varname)-3:end)=[];
    
    % reinstate names for those ftrs whose name changed for publication
    if strcmp(varname, 'QuadE_O')==1; varname='SinE_O'; end
    if strcmp(varname, 'QuadI50_O')==1; varname='SinI50_O'; end
    if strcmp(varname, 'QuadI_T')==1; varname='SinI_T'; end
    if strcmp(varname, 'AreaUnderPeaksI_O')==1; varname='SS_Area_O'; end
    if strcmp(varname, 'InspVol_03Ti_O')==1; varname='MorgensV03_O'; end
    if strcmp(varname, 'VPEF_VTe_O')==1; varname='ATS_VPEF_VTe_O'; end
    if strcmp(varname, 'InspFlutPow4to7_O')==1; varname='InspFlutPowOrig_O'; end
    
    index = find(strcmp(FeatureNames.Name, varname));
    
    upper_flow(ii) = Pneumo_lims.centileupper(index);
    lower_flow(ii) = Pneumo_lims.centilelower(index);
    
    upper_pnasal(ii) = Pnasal_lims.centileupper(index);
    lower_pnasal(ii) = Pnasal_lims.centilelower(index);
end
disp('done');

%%
FinalModelTable.upper = upper_pnasal;
FinalModelTable.lower = lower_pnasal;
FinalModelTable.upperFlow = upper_flow;
FinalModelTable.lowerFlow = lower_flow;

%%
save('FinalModel_Table.mat', 'FinalModelTable');



%% 10 Hz model                                                                                                    

close
clear
clc

%%
load('FinalModel_Table10.mat')   
load('FtrNames165.mat')
Pnasal_lims = load('FeatureLimits_Pnasal.mat');
Pneumo_lims = load('FeatureLimits_Pneumo.mat');

%%
upper_flow = NaN(height(FinalModelTable10),1);
lower_flow = NaN(height(FinalModelTable10),1);
upper_pnasal = NaN(height(FinalModelTable10),1);
lower_pnasal = NaN(height(FinalModelTable10),1);

%%
for ii = 2:size(FinalModelTable10)
    varname = FinalModelTable10.Feature{ii};
    varname(length(varname)-3:end)=[];
    
    % reinstate names for those ftrs whose name changed for publication
    if strcmp(varname, 'QuadE_O')==1; varname='SinE_O'; end
    if strcmp(varname, 'QuadI50_O')==1; varname='SinI50_O'; end
    if strcmp(varname, 'QuadI_T')==1; varname='SinI_T'; end
    if strcmp(varname, 'AreaUnderPeaksI_O')==1; varname='SS_Area_O'; end
    if strcmp(varname, 'InspVol_03Ti_O')==1; varname='MorgensV03_O'; end
    if strcmp(varname, 'VPEF_VTe_O')==1; varname='ATS_VPEF_VTe_O'; end
    if strcmp(varname, 'InspFlutPow4to7_O')==1; varname='InspFlutPowOrig_O'; end
    
    index = find(strcmp(FeatureNames.Name, varname));
    
    upper_flow(ii) = Pneumo_lims.centileupper(index);
    lower_flow(ii) = Pneumo_lims.centilelower(index);
    
    upper_pnasal(ii) = Pnasal_lims.centileupper(index);
    lower_pnasal(ii) = Pnasal_lims.centilelower(index);
end
disp('done');

%%
FinalModelTable10.upper = upper_pnasal;
FinalModelTable10.lower = lower_pnasal;
FinalModelTable10.upperFlow = upper_flow;
FinalModelTable10.lowerFlow = lower_flow;

%%
save('FinalModel_Table10.mat', 'FinalModelTable10');