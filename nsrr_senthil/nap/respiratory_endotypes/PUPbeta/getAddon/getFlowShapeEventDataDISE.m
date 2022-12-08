function FlowShapeT = getFlowShapeEventData(MrangeOverride)
global settings AMasterSpreadsheet
%errors: SnoreDB95thcentile is incorrect, Etype is missing

if 0
settings=ImportSettings(settings,AMasterSpreadsheet);

[num,patients,~] = xlsread(AMasterSpreadsheet,1,'AD4:AH10003');
NaNlist = (isnan(num)); num(NaNlist) = 0; % set missing or NaN as 0
analyzelist = logical(num(:,2));
settings.invertflowlist = logical(num(:,1));
settings.protocol = patients(:,3);
settings.patients = patients;
end

M=size(settings.patients,1);
Mrange=1:M;
if exist('MrangeOverride')
    Mrange=MrangeOverride;
end

dir=[settings.AnalyzedDirectory '\EventAnalyzed\'];
% dir=[settings.AnalyzedDirectory '\'];
for i=Mrange
    success(i)=0;
    try
    clear BreathDataTableFulls
    loadstr = [dir settings.savename '_' num2str(i) '.mat'];
    load(loadstr,'EAinfo','BreathDataTableFulls');
%     temp = mergestructs( EAinfo.FlowShape,EAinfo.EvtFtrs.EventDepth)   ; 
    %temp.Patient_ID = [settings.savename '_' num2str(i) ];
%     T = struct2table(temp);
    FLnames = erase(string(fieldnames(EAinfo.FlowShape)),"nanmean_");
    FLidx = ismember(BreathDataTableFulls.Properties.VariableNames,FLnames);
    criteria = BreathDataTableFulls.FlowDrive < 0.6;
    T=varfun(@nanmean,BreathDataTableFulls(criteria,FLidx));
%     T = array2table(TTemp); clear TTemp;
%     T.Properties.VariableNames = BreathDataTableFulls.Properties.VariableNames(FLidx);
    
    if ~exist('FlowShapeT')
        FlowShapeT=T;
        FlowShapeT{:,:}=NaN; %hack initialization of Table
    end
    FlowShapeT{i,:}=T{:,:}; %may need a join function (join, outerjoin, innerjoin) if tables are different widths
    success(i)=1;    
    catch me
        FlowShapeT{i,:}=NaN;
        success(i)=0;
        disp(i)
        disp(me.message)
    end
end

FlowShapeT{success==0,:}=NaN;





