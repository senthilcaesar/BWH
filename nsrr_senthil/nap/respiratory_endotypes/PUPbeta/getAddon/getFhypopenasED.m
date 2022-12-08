function FhypopneasED = getFhypopenasED(MrangeOverride)
% RUN StartHere.m first

% changing to run from converted files-1/15/2021


global settings AMasterSpreadsheet ChannelsList

t_startGetData = clock;

%%

[~,path,~] = xlsread(AMasterSpreadsheet,1,'A26'); % only used incase of analyzed file

[~,~,MasterWorksheet] = xlsread(AMasterSpreadsheet,1,'T4:BF10000');
lastrow = find(1*(~isnan(cell2mat(MasterWorksheet(:,9)))),1,'last');

MasterWorksheet(lastrow+1:end,:)=[];

ConvertList = cell2mat(MasterWorksheet(:,9));

NaNlist = (isnan(ConvertList));
ConvertList(NaNlist) = 0; % set missing or NaN as 0


Filenames = MasterWorksheet(:,2:8);


%%

MaxM = size(ConvertList,1); % right now using all the converted files.

try
    settings.Mrange=1:MaxM;
    M = max(settings.Mrange);
end

%%
if exist('MrangeOverride')
    settings.Mrange = MrangeOverride;
end
%%
clear x T x2
% T.AHIorig = nan(M,1);
% T.FhypopneasOrig = nan(M,1);
% T.Fhypopneas4 = nan(M,1);
% T.AHI4 = nan(M,1);
% HBtotal = nan(M,1);
success = zeros(M,1);

FhypopneasED = nan(MaxM,1);
for i=settings.Mrange
    disp(['Processing: ' num2str(i) '/' num2str(settings.Mrange(end))]);
    try
    filepathEA = [settings.AnalyzedDirectory, 'EventAnalyzed\'];
    EA = load([filepathEA, settings.savename, '_', num2str(i),'.mat']);
    ctridx = 100;
    hypopneas = ones(size(EA.Boxes.VI,1),1);% assume all hypopneas
    for jj = 1:size(EA.Boxes.VI,1)
        evtstart = ctridx-find(EA.Boxes.EventsResp(jj,ctridx:-1:1)==0,1,'first')+2;
        evtend = ctridx + find(EA.Boxes.EventsResp(jj,ctridx:1:end)==0,1,'first')-2;

        if any(EA.Boxes.VI(jj,evtstart:evtend)<0.1)
            hypopneas(jj) = 0; % take it 
        end
    end
    FhypopneasED(i) = sum(hypopneas)/size(EA.Boxes.VI,1);
    end
end