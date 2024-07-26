function T = compilePSGReports(MrangeOverride)
% RUN StartHere.m first
% compiles PSGreport data from converted files

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

MaxM = size(ConvertList,1); % right now using all the converted files.
try
    settings.Mrange=1:MaxM;
end

if exist('MrangeOverride')
    settings.Mrange = MrangeOverride(:)';
end

M = max(settings.Mrange);
%%
clear x T x2
success = zeros(M,1);

for i=settings.Mrange
    disp(['Processing: ' num2str(i) '/' num2str(settings.Mrange(end))]);

    tempfilename=Filenames{i,1};
    tempfilename(find(tempfilename=='.'):end)=[];
    filedirC = [settings.workdir 'Converted\' tempfilename '_XHz.mat'];
    if exist(filedirC)==2
        try
            clear Evts
            load(filedirC,'Evts');
            if ~exist('Evts')
                disp('no Evts var');
                continue
            end            

            T(i,:)=Evts.PSGreportSummaryT;

            success(i)=1;
        catch me
        end
    else
        disp(['File absent: ' num2str(i) '/' num2str(settings.Mrange(end)) ', ' filedirC]);
    end
end

%%
T = array2table(T);
T{success(1:height(T))==0,:}=NaN;

%%
% savefilename = [settings.workdir 'Summary\getData_' datestr(now,'yyyymmdd HHMMSS')];
% savefilename(end-6)='T';
% save(savefilename,'T','-v7.3')
% savefilename = [settings.workdir 'Summary\getData'];
% save(savefilename,'T','-v7.3')
% 
% delta_tGetData = etime(clock, t_startGetData); % delta in seconds
% D = duration(0,0,round(delta_tGetData),'Format','hh:mm:ss'); % duration in HH:MM:SS
% disp(' '); % add row space for visual clarity in command window
% displaytext = ['GetData Calculation Complete. Total time: ', char(D), ' (hh:mm:ss)'];
% disp(displaytext);

