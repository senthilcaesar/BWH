function getNewXHzSignalFromEDF()
global settings

%where is source data
AMasterSpreadsheet = [settings.AMasterdir, 'AMasterSpreadsheet.xlsx'];
% read spreadsheet (Files worksheet)
[~,~,raw] = xlsread(AMasterSpreadsheet,1,'T4:BF10000');
ConvertMatFlag = cell2mat(raw(:,9));
%last row is based on whether "Convert?" (ConvertMatFlag) has numeric data
lastrow = find(1*(~isnan(ConvertMatFlag)),1,'last');
raw(lastrow+1:end,:)=[];
ConvertMatFlag(lastrow+1:end,:)=[];
Filenames = raw(:,2:8);

%where is converted XHz data
[~,~,options] = xlsread(AMasterSpreadsheet,2,'C3:C11');
%F_samp = options{1};
%settings.Fs=F_samp;
ExportDataDirectory = char(options(2));
if ExportDataDirectory(end)~=filesep
    ExportDataDirectory=[ExportDataDirectory filesep];
end
settings.FsLeak = double(options{7});

%what new channel to bring in from source
NewChannel='Pulse';
NewChannelXHzName='Pulse';
InterpMethod='linear';

successfullysaved = zeros(size(Filenames,1),1);

    for i=1:size(Filenames,1)
        try
        EDFfilenamedir = [Filenames{i,4} Filenames{i,1}];
        XHzfilenamedir = [ExportDataDirectory [Filenames{i}(1:end-4) '_XHz' '.mat']];
        %check exists here

        [Label,Transducer,Fs] = EDFChannelLabels(EDFfilenamedir);
        
        Label2 = string(Label);
        k = strfind(Label2,NewChannel);
        k = cellfun(@length,k);
        k = find(k==1);
        ExactLabel = Label(k); %check
        strcmp(NewChannel,Label);
        
        [data,Fsdata,~,~,Labeldata,~,~,~,~] = readedfrev3(EDFfilenamedir,k-1,0,Inf); %the '-1' is because EDF channel numbers start from 0.
                
        Fsdata = Fsdata*(1-settings.FsLeak);
        
        if strcmp(NewChannelXHzName,'Pulse')
            InterpMethod = 'nearest'; %overwrite default above
            data(data<20)=NaN;
            data(data>200)=NaN;
        end
        
        matObj = matfile(XHzfilenamedir);
        varlist = who(matObj);
        load(XHzfilenamedir);
        
        
        %settings.Fs = 1./(DataEventHypnog_Mat(2,1)-DataEventHypnog_Mat(1,1));
        
        StartTime = DataEventHypnog_Mat(1,1);
        %[DataEventHypnog_Mat,ChannelsList,ChannelsFs]=runWSanalysis(DataEventHypnog_Mat,ChannelsList,ChannelsFs);
        TimeTemp=(StartTime:(1/Fsdata):StartTime+(length(data)-1)*(1/Fsdata))'; % This is the time vector associated with the 100Hz Flow data.
        Temp = interp1(TimeTemp,data,DataEventHypnog_Mat(:,1),InterpMethod);
        
        if 0
            figure(1)
            plot(DataEventHypnog_Mat(:,1),Temp)
            hold('on')
            
        end
        
        %add code to check if the signal already exists and update, else add to end:
        INewChannel = size(DataEventHypnog_Mat,2) + 1;
        
        ChannelsFs(INewChannel)=Fsdata; %correction of historical error
        ChannelsList{INewChannel}=NewChannelXHzName; %correction of historical error
        DataEventHypnog_Mat(:,INewChannel) = Temp;
        
        save(XHzfilenamedir,varlist{:},'-v7.3')
            i
        successfullysaved(i)=1;
        
        clear DataEventHypnog_Mat ChannelsFs ChannelsList data Temp
        catch me
            i
            'failure'
        end
    end
end