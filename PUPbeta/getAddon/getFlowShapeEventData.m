function FlowShapeT = getFlowShapeEventData(MrangeOverride)
global settings AMasterSpreadsheet
%errors: SnoreDB95thcentile is incorrect, Etype is missing
%
settings=ImportSettings(settings,AMasterSpreadsheet);

[num,patients,~] = xlsread(AMasterSpreadsheet,1,'AD4:AH10003');
NaNlist = (isnan(num)); num(NaNlist) = 0; % set missing or NaN as 0
analyzelist = logical(num(:,2));
settings.invertflowlist = logical(num(:,1));
settings.protocol = patients(:,3);
settings.patients = patients;

M=size(patients,1);
Mrange=1:M;
if exist('MrangeOverride')
    Mrange=MrangeOverride;
end

dir=[settings.workdir 'Analyzed\'];
%dir=[settings.OutputDataDirectory 'EventAnalyzed\'];

for i=Mrange
    success(i)=0;
    try
    clear EAinfo
    loadstr = [settings.AnalyzedDirectory 'EventAnalyzed\' settings.savename '_' num2str(i) '.mat'];
    load(loadstr,'EAinfo');
    temp = mergestructs( EAinfo.FlowShape,EAinfo.EvtFtrs.EventDepth)   ; 
    %temp.Patient_ID = [settings.savename '_' num2str(i) ];
    T = struct2table(temp);
    
    if ~exist('FlowShapeT')
        FlowShapeT=T;
        FlowShapeT{:,:}=NaN; %hack initialization of Table
    end
    FlowShapeT{i,:}=T{:,:}; %may need a join function (join, outerjoin, innerjoin) if tables are different widths
    success(i)=1;    
    catch me
        success(i)=0;
        disp(i)
        disp(me.message)
    end
end

FlowShapeT{success==0,:}=NaN;





