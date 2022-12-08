function getConvertQc(MrangeOverride)
% RUN StartHere.m first

global AMasterSpreadsheet settings

t_startGetData = clock;

%% read excel spreadsheet

[~,txt] = xlsread(AMasterSpreadsheet,1,'AD4');
M = size(txt,1);

ID=1:M;

% override if mrange present
if exist('MrangeOverride')
    ID = MrangeOverride(:);
end

success=zeros(length(ID),1);
ConvDir=settings.ConvertedDirectory;
ConvQcT=table();

ChannelsMasterList={'Flow','Thorax','Abdomen','SpO2','EEG1','EEG3','EEG5',...
    'EEG7','EEG9','EEG11','EKG','LOC','ROC','Position'};


for i=1:length(ID)
    try
        clear W temp
        disp(['processing: ' txt{ID(i),1}])
        try
            W = load([ConvDir txt{ID(i),1}]);
        end
        if exist('W','var')
            
            % 1- what signals are present and their fs
            
            temp= W.SigT.Properties.VariableNames;
            P=ismember(ChannelsMasterList,temp);
            SigsPresent=ChannelsMasterList(P);
            
            Fs_Flow=W.ChannelsFs(2); % fs of flow
            
            
            % 2- Events and hypnogram present or absent
            clear temp
            temp=nanmean(W.Evts.Hypnogram,1);
            if temp==8
                HypnogramInfo="No Staging";
            elseif temp==4
                HypnogramInfo="Fully Wake";
            else
                HypnogramInfo="Stages Present";
            end
            
            
            if isempty(W.Evts.RespT)
                RespEvtsInfo="No Resp Events";
            else
                RespEvtsInfo="Resp Events Present";
            end
            
            
            if isempty(W.Evts.ArT)
                ArousalInfo="No Arousals";
            elseif round(nansum(W.Evts.ArT.EventDuration)/height(W.Evts.ArT))==3
                ArousalInfo="Arousals are 3s duration";
            else
                ArousalInfo="Arousals OK";
            end
            
            clear temp
            temp=sum(isnan(W.Evts.SpO2{1,:}));
            if temp/width(W.Evts.SpO2)==1 % means all columns are nans
                SpO2SummaryInfo="No SpO2 Summary Stats";
            else
                SpO2SummaryInfo="SpO2 Summary Stats Present";
            end
            
            
            
            % 3- flow quality information
            LowPassPredict=W.Info.FlowQ.FlowFilterDetect.LowPassPredict_2ndOrderEquivalent(1);
            HighPassPredict=W.Info.FlowQ.FlowFilterDetect.HighPassPredict_1stOrderEquivalent(1);
            W.Info.FlowQ=rmfield(W.Info.FlowQ,'FlowFilterDetect');
            W.Info.LowPassPredict=LowPassPredict;
            W.Info.HighPassPredict=HighPassPredict;
            FlowQcT=struct2table(W.Info.FlowQ);
            
            %max(W.Info.WSinfo.Acc)
            %max(W.Info.WSinfo.AccPred)
            
        end
        
        
        
        
    end
end