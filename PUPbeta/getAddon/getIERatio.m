function getIERatio(MrangeOverride)

% RUN StartHere.m first

% takes lot of time to complete ~30min/file -- almost half of analysis is repeated here
% to get IE ratio

%%
global settings AMasterSpreadsheet

t_start = clock;
[~,path,~] = xlsread(AMasterSpreadsheet,1,'A26');

[num,~,MasterWorksheet] = xlsread(AMasterSpreadsheet,1,'AD4:AH10003');
lastrow = find(cell2mat(cellfun(@(x)any(~isnan(x)),MasterWorksheet(:,1),'UniformOutput',false)),1,'last');
MasterWorksheet(lastrow+1:end,:)=[];
Filenames = MasterWorksheet(:,1:2);

NaNlist = (isnan(num)); num(NaNlist) = 0; % set missing or NaN as 0

settings = ImportSettings(settings,AMasterSpreadsheet);

MaxM = size(num,1);

try
    M = max(settings.Mrange);
catch
    settings.Mrange=1:MaxM;
    M = max(settings.Mrange);
end

if exist('MrangeOverride')
    Mrange=MrangeOverride;
else
    Mrange=settings.Mrange;
end
% intialize variables to zeros/NaN incase FL computation failed or no
% analysis file present
success = zeros(numel(Mrange),1);

%% overcome parfor vs global variable challenge
settings1=settings;

%% try parallel pool to make it faster
parfor n1=1:numel(Mrange)
    
    n=Mrange(n1);
    try
        AnalyzedDir = [path{:} settings1.savename '_' num2str(n) '.mat'];
        
        if exist(AnalyzedDir,'file')==2
            disp(['Processing: ' num2str(n) ': ' settings1.savename '_' num2str(n) '.mat']);
            W = load(AnalyzedDir);
            BreathDataTable=[];
            BreathDataTable2=[];
            temp2=[];
            BreathDataTable=W.BreathDataTable;
            
            if size(BreathDataTable,1)==1 && size(BreathDataTable,2)==1
                BreathDataTable2=BreathDataTable{1,1};
            else
                BreathDataTable2=BreathDataTable;
            end
            
            %%
            % check if IE ratio is already present or not
            temp2=BreathDataTable2;
            IEratiofound=zeros(size(temp2,2),1);
            for ii=1:length(temp2) % remove any empty cells in BreathDataTable
                if isempty(temp2{1,ii})
                    temp2{1,ii}=NaN;
                end
                try
                    IEratiofound(ii)=sum(ismember(temp2{1,ii}.Properties.VariableNames,'IEratio')); %check if IEratio exist.
                end
            end
            BreathDataTable2=temp2;
            
            % if IE ratio not found, calculate and append to
            % BreathdataTable and save Analysis file
            if sum(IEratiofound)==0
                try
                    disp('IEratio not found in BreathDataTable--Calculating IE ratio for windows');
                    
                    % Load Converted files
                    ConvDir=[Filenames{n,2} Filenames{n,1} '.mat'];
                    
                    W2=load(ConvDir);
                    if isfield(W2,'DataEventHypnog_Mat')
                        SigT=array2table(W2.DataEventHypnog_Mat);
                        SigT.Properties.VariableNames = W2.ChannelsList;
                    else
                        SigT=W2.SigT;
                    end
                    
                    
                    % downsample flow
                    if (settings1.downsampledFs~=0)&&(settings1.downsampledFs~=settings1.Fs)
                        if strcmp(settings1.FlowSignal,'SecondPnasal')
                            Flow = SigT.Pnasal;
                            disp('Using Pnasal signal in this analysis');
                        else
                            Flow = SigT.Flow;
                        end
                        displaytext=['Resampling: Flow data from ' num2str(settings1.Fs) ' to ' num2str(settings1.downsampledFs) 'Hz...(and back)'];
                        disp(displaytext);
                        TotalTime = (length(Flow)-1)*(1/settings1.Fs);
                        Time_ = [0:(1/settings1.Fs):TotalTime]';
                        TimeDS=[0:(1/settings1.downsampledFs):TotalTime]';
                        FlowDS = interp1(Time_,Flow,TimeDS,'linear'); % downsample
                        Flow = interp1(TimeDS,FlowDS,Time_,'linear');  % upsample back to original length.
                        Time_=[];
                        TimeDS=[];
                        FlowDS=[];
                        
                        % ensure correct length
                        mismatch = length(Flow)-length(SigT.Flow);
                        if mismatch>=0  % truncate Flow
                            Flow = Flow(1:end-mismatch);
                        else  % add spacer to Flow
                            Flow = [Flow; zeros(abs(mismatch),1)];
                        end
                        if strcmp(settings1.FlowSignal,'SecondPnasal')
                            SigT.Pnasal = Flow; %write it back into DataEventHypnog_Mat
                        else
                            SigT.Flow = Flow;
                        end
                    end
                    
                    
                    % calculate IE ratio
                    for ii=1:length(BreathDataTable2)%47:54
                        BreathDataTableWin=[];
                        BreathDataTableWin=BreathDataTable2{1,ii};
                        
                        disp(['processing window:' num2str(n) ':' num2str(ii) '/' num2str(length(BreathDataTable2))])
                        
                        try
                            starttime =BreathDataTableWin.Time0(1);
                            endtime=starttime+(settings1.windowlength*60);
                            starti= find(SigT.Time==starttime);
%                             disp(num2str(starti))
                            endi= find(SigT.Time>=endtime,1,'first')-1; % some time issues in windows
%                             disp(num2str(endi))
                            FlowWin=Flow(starti:endi);
                            TimeWin=SigT.Time(starti:endi);
                            
                            [~,~,~,~,~,~,~,~,~,~,...
                                IEratio,~,~,~,~,~,~,~,~,~,~,IEratioEstimated] =...
                                VEfromFlow(TimeWin,FlowWin,settings1);
                            
                            
                            BreathDataTableWin.IEratio = IEratio+0*BreathDataTableWin.BB_i_start;
                            BreathDataTableWin.IEratioEstimated = IEratioEstimated+0*BreathDataTableWin.BB_i_start;
                            BreathDataTable2{1,ii}=BreathDataTableWin;
                            
                        catch me
                        end
                        
                    end
                    
                end
                W.BreathDataTable = BreathDataTable2;
                parsave(AnalyzedDir,W); % regular save not working for parallel processing
                success(n1)=1;
                
            elseif sum(IEratiofound)>0
                disp(['IE ratio already exists. Skipping...' num2str(n)])
                success(n1)=0;
            end
            
        end
    catch
        disp('failed IE Ratio Hack');
        success(n1)=0;
    end
end

% save success variable
fname=[settings.workdir,'Summary' filesep 'successIERatio.mat'];
save(fname,'success','-v7.3');

delta_t = etime(clock, t_start); % delta in seconds
D = duration(0,0,round(delta_t),'Format','hh:mm:ss'); % duration in HH:MM:SS
disp(' '); % add row space for visual clarity in command window
displaytext = ['IE Calculation Complete. Total time: ', char(D), ' (hh:mm:ss)'];
disp(displaytext);
end
