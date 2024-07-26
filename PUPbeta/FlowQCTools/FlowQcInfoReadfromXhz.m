function FlowQcInfoReadfromXhz(MrangeOverride)
% Run start here first
global AMasterSpreadsheet settings

[~,~,MasterWorksheet] = xlsread(AMasterSpreadsheet,1,'T4:BF10000');
Filenames = MasterWorksheet(:,11:12);

I=find(cell2mat(cellfun(@(x)any(~isnan(x)),Filenames(:,1),'UniformOutput',false)));
Filenames((I(end)+1):end,:)=[]; %remove data below last legit entry

Mrange=1:size(Filenames,1);
if exist('MrangeOverride')
    Mrange = MrangeOverride;
end

FlowQcT=[];
for ii=Mrange
    ii
    try
        Convname=[Filenames{ii,2} Filenames{ii,1} '.mat'];
        load(Convname,'Info')
        LowPassPredict=Info.FlowQ.FlowFilterDetect.LowPassPredict_2ndOrderEquivalent(1);
        HighPassPredict=Info.FlowQ.FlowFilterDetect.HighPassPredict_1stOrderEquivalent(1);
        temp=[table(string(Filenames{ii,1})) struct2table(Info.FlowQ) table(LowPassPredict,HighPassPredict,ii)];
        FlowQcT=[FlowQcT;temp];
    catch
%         N=(NaN(1,width(FlowQcT)-1));
%         N=array2table(N);
%         temp=[table(string(Filenames{ii,1})) N];
%         
%         temp.Properties.VariableNames=FlowQcT.Properties.VariableNames;
%         FlowQcT=[FlowQcT;temp];
    end
end

save([settings.workdir '\Summary\FlowQcT.mat'],'FlowQcT');
