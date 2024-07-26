function getSNRWinAnalyzed(MrangeOverride)
% MrangeOverride=[4:250];
%% RUN START HERE FIRST

global settings AMasterSpreadsheet ChannelsList 

if ~isfield(settings,'ImportedSettingsComplete') %|| settings.ImportedSettingsComplete==0
    settings = ImportSettings(settings,AMasterSpreadsheet);
end

%% LOAD AMASTERSPREADSHEET
[~,path,~] = xlsread(AMasterSpreadsheet,1,'A26');
[~,~,MasterWorksheet] = xlsread(AMasterSpreadsheet,1,'AD4:AH10003');
lastrow = find(cell2mat(cellfun(@(x)any(~isnan(x)),MasterWorksheet(:,1),'UniformOutput',false)),1,'last');
MasterWorksheet(lastrow+1:end,:)=[];
Filenames = MasterWorksheet(:,1:2);

%[~,~,raw] = xlsread(AMasterSpreadsheet,2,'C20:C57');
%settings.savename = char(raw{1});

AnalyzeList=cell2mat(MasterWorksheet(:,5));
MaxM = size(AnalyzeList,1); % using analyzed files.

Mrange=1:MaxM;

if isfield(settings,'Mrange')
    Mrange = settings.Mrange;
end

%%
if exist('MrangeOverride','var')
    Mrange = MrangeOverride;
end

M = max(Mrange);
success = zeros(M,1);

%%
for n=Mrange
    clear BreathDataTable2 BreathDataTable
    AnalyzedDir = [settings.AnalyzedDirectory settings.savename '_' num2str(n) '.mat'];
    ConvDir=[settings.ConvertedDirectory Filenames{n,1} '.mat'];
    
    if exist(AnalyzedDir)==2 && exist(ConvDir)==2
        disp(['Processing: ' num2str(n) ': ' settings.savename '_' num2str(n) '.mat']);
        try
             W = load (AnalyzedDir);
             BreathDataTable=W.BreathDataTable;
             
             W2=load(ConvDir);
             if isfield(W2,'DataEventHypnog_Mat')
                SigT=array2table(W2.DataEventHypnog_Mat);
                 SigT.Properties.VariableNames = W2.ChannelsList;
                 W2 = rmfield(W2, {'DataEventHypnog_Mat','ChannelsList'});
                 W2.SigT = SigT;
                 save(ConvDir,'-struct','W2','-v7.3');
             else
                 SigT=W2.SigT;
            end
            
            Time=SigT.Time;
            %         Flow=DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'Flow')==1));
            
            
            %% downsample flow
            if (settings.downsampledFs~=0)&&(settings.downsampledFs~=settings.Fs)
                
                if strcmp(settings.FlowSignal,'SecondPnasal')
                    Flow = SigT.Pnasal;
                    disp('Using Pnasal signal in this analysis');
                else
                    Flow = SigT.Flow;
                end
                
                
                %displaytext=['Resampling: Flow data from ' num2str(settings.Fs) ' to ' num2str(settings.downsampledFs) 'Hz...(and back)'];
                %disp(displaytext);
                %                 set(handletext,'String',displaytext); drawnow;
                TotalTime = (length(Flow)-1)*(1/settings.Fs);
                Time_ = [0:(1/settings.Fs):TotalTime]';
                TimeDS=[0:(1/settings.downsampledFs):TotalTime]';
                FlowDS = interp1(Time_,Flow,TimeDS,'linear'); % downsample
                Flow = interp1(TimeDS,FlowDS,Time_,'linear');  % upsample back to original length.
                clear Time_ TimeDS FlowDS
                
                % ensure correct length
                mismatch = length(Flow)-length(SigT.Flow);
                if mismatch>=0
                    % truncate Flow
                    Flow = Flow(1:end-mismatch);
                else
                    % add spacer to Flow
                    Flow = [Flow; zeros(abs(mismatch),1)];
                end
                %then write it back into DataEventHypnog_Mat
                if strcmp(settings.FlowSignal,'SecondPnasal')
                    SigT.Pnasal = Flow;
                else
                    SigT.Flow= Flow;
                end
                
            end
            
            %%
            %disp('warning: SS has edited code here to deal with BreathDataTable2=BreathDataTable{1,1} error')
            if size(BreathDataTable,1)==1 && size(BreathDataTable,2)==1
                BreathDataTable2=BreathDataTable{1,1};
            else
                BreathDataTable2=BreathDataTable;
            end
            
            temp2=BreathDataTable2;
            for i=1:length(temp2) % remove any empty cells in BreathDataTable
                if isempty(temp2{1,i})
                    temp2{1,i}=NaN;
                end
                try
                    temp2{1,i}(:,{'SNRwindow','SNRbreath'})=[]; %removes them if they already exist. 
                end
            end
            BreathDataTable2= temp2;
            
            
            for i=1:length(BreathDataTable2)
                clear BreathDataTableWin FlowWin TimeWin
                BreathDataTableWin=BreathDataTable2{1,i};
                if isa(BreathDataTableWin,'double') && isnan(BreathDataTableWin)
                    continue
                end
                try
                    starttime =BreathDataTableWin.Time0(1);
                    endtime=starttime+(settings.windowlength*60);
                    starti= find(Time==starttime);
                    endi= find(Time>=endtime,1,'first')-1; % some time issues in windows
                    FlowWin=Flow(starti:endi);
                    TimeWin=Time(starti:endi);
                    BreathDataTableWin = SignalToNoiseWin(BreathDataTableWin,FlowWin,TimeWin);
                    BreathDataTable2{1,i}=BreathDataTableWin;
                    %             BreathDataTable{1,i}=BreathDataTableWin;
                    % BreathDataTable{1,1}{1,i}=BreathDataTableWin;
                catch me
                    disp(['error win:' num2str(i)])
                    disp(me.message)
                end
                %i
            end
            W.BreathDataTable = BreathDataTable2;
%             save(AnalyzedDir,'BreathDataTable','-append')
            save(AnalyzedDir,'-struct','W','-v7.3');
            success(n)=1;
            disp(['Completed n=' num2str(n)]);
        catch
            disp('failed SNR Hack');
        end   
    end
end