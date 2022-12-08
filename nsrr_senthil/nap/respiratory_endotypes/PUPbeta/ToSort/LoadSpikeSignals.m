

MultiplyPositionChannelby=5; %if supine~=0.33 use value of "6", if supine=2 use a value of "1".
    
    clear channelnameoptions
    channelnameoptions.Evts={'Evts','New_Evts'};
    channelnameoptions.Evts2={'Evts2','New_Evts2'};
    channelnameoptions.Epochs={'Epochs'};
    channelnameoptions.SaO2={'SaO2','SpO2','Sat','Sao2','Spo2','O2sat','o2sat','Sao2fing'};
    channelnameoptions.Position={'Position','Pos','pos','position'};
    %channelnameoptions.PCO2={'PCO2','CO2_Ana','CO2_anal'};
    %channelnameoptions.pO2={'pO2','O2_Ana','O2_anal'};
    channelnameoptions.EEG={'WakeSleep','EEG','EEG_C3_A2','C3_A2'};
    
    channelnameoptions.EEG1={'EEG1','EEG_C3_A2_clean','EEG_C3_A2','C3_A2'};
    channelnameoptions.EEG2={'EEG2','EEG_C4_A1_clean','EEG_C4_A1','C4_A1'};
    channelnameoptions.EEG3={'EEG3','EEG_F3_A2_clean','EEG_F3_A2','F3_A2'};
    channelnameoptions.EEG4={'EEG4','EEG_O2_A1_clean','EEG_O2_A1','O2_A1'};
    channelnameoptions.EKG={'EKG','ECG'};
    
    channelnameoptions.Thorax={'ThNoxRIP','Thorax','RC','Chest','CHEST','Belt2'};
    channelnameoptions.Abdomen={'AbNoxRIP','Abdomen','ABD','Abdom','ABDM','Belt1'};
    
    channelnameoptions.alphaFlow={'alphaFlow'};
    channelnameoptions.kFlow={'kFlow'};
    
    channelnameoptions.Pes={'Pes_clean','Pes','pes','Pepi','pepi'};
    channelnameoptions.Edi={'Edi'};
    channelnameoptions.FlowEdi={'FlowEdi'};
    channelnameoptions.FlowPes={'FlowPes'};
    channelnameoptions.GGpmax={'GGpmax','GGPmax'};
    %channelnameoptions.WakeSleep={'WakeSleep'};
    
    %channelnameoptions.EKG={'EKG','ECG'};
    if 0 %use pnasal signal
        channelnameoptions.Pnasal={'Pnasal','PNasal','Pmask','PMask'}; %'Vflow'
    else %use pneumotach flow
        channelnameoptions.Pnasal={'Vflow','Flow','PNasal','Pmask','PMask'}; %'Vflow'
    end
    channelnameoptions.CPAP={'CPAP','Pmask','PMask'}; 
    channelnamestemp=fieldnames(channelnameoptions);
    
    for i=1:length(channelnamestemp)
        temp=eval(['channelnameoptions.' char(channelnamestemp(i))]);
        foundamatch=0;
        for nn=1:length(temp)
            %Does it exist?
            for j=1:length(w)
                if strcmp(w(j).name,char(temp(nn)))
                    eval([channelnamestemp{i} '=filehandle.' char(temp(nn))]);
                    foundamatch=1;
                    break
                end
            end
            if foundamatch
                break
            end
        end
    end
    
    Channels={'Pnasal','Thorax','Abdomen','SaO2','EEG','Position','CPAP','Pes','Edi','GGpmax','FlowEdi','FlowPes','EEG1','EEG2','EEG3','EEG4','EKG','alphaFlow','kFlow'};
    optionalchannelsstartat=8;
    dt=1/F_samp;
    
    %% Resample and adjust [move EKG and EEG analysis to before the resampling here]
    for i=1:length(Channels)
        displaytext=['Collecting channel: ' Channels{i}];
        try
            disp(displaytext); set(handletext,'String',displaytext); drawnow;
            eval(['Fs=(1/' Channels{i} '.interval);']);
            if round(Fs)~=F_samp
                displaytext=['Resampling: ' Channels{i} ' from ' num2str(round(Fs)) ' to ' num2str(F_samp) 'Hz'];
                disp(displaytext); set(handletext,'String',displaytext); drawnow;
                eval([Channels{i} '.values = resample(' Channels{i} '.values,F_samp,round(Fs));']); %only works if Fsamp / Fs are integer multiples of each other
                eval([Channels{i} '.interval = 1/F_samp;']); %only works if Fsamp / Fs are integer multiples of each other
                eval([Channels{i} '.length = length(' Channels{i} '.values);']); %only works if Fsamp / Fs are integer multiples of each other
            end
        catch me
            displaytext=['No channel: ' Channels{i} ' ' me.message];
            disp(displaytext); set(handletext,'String',displaytext); drawnow;
            if i<optionalchannelsstartat %make empty channel with values = -1;
                eval([Channels{i} '.values = -1+ 0*Pnasal.values;']); %assumes we at least have a Pnasal channel
                eval([Channels{i} '.interval = 1/F_samp;']); 
                eval([Channels{i} '.length = length(' Channels{i} '.values);']); 
            else
                Channels{i}=[]; %make cell empty, later removed...
            end
        end
    end
    
    Channels(find(cellfun(@isempty,Channels)))=[];
    
    %% Get timing, epoch and event data
    
    %get time data
    timeXHz=(Pnasal.start:(1/F_samp):(Pnasal.start+(length(Pnasal.values)/F_samp)-1/F_samp))';
    N=length(timeXHz);
    
    %reserve memory for the large matrix
    DataEventHypnog_Mat=zeros(N,12);
        DataEventHypnog_Mat(:,1)=timeXHz;
        %clear timeXHz;
        
%% New: fix channel lengths if not = N (not the same as Pnasal) -- incorporate this into Alice and others
    for i=1:length(Channels)
        if isempty(Channels{i})
            continue
        end
        eval(['ChannelN(i)=length(' Channels{i} '.values);']);
        if ChannelN(i)>N
            displaytext=['Length of ' Channels{i} ' channel is being altered from ' num2str(ChannelN(i)) ' to ' num2str(N)];
            disp(displaytext); set(handletext,'String',displaytext); drawnow;
            eval([Channels{i} '.values(N+1:end)=[];']);
            eval([Channels{i} '.length=N;']);
        end
        if ChannelN(i)<N
            displaytext=['Length of ' Channels{i} ' channel is being altered from ' num2str(ChannelN(i)) ' to ' num2str(N)];
            disp(displaytext); set(handletext,'String',displaytext); drawnow;
            eval([Channels{i} '.values((ChannelN(i)+1):N)=' Channels{i} '.values(end);']);
            eval([Channels{i} '.length=N;']);
        end
    end