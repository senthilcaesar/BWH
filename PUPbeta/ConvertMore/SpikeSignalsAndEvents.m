%Script that handles the Spike signals and events 

filehandle = matfile(filenameanddir);
w = whos('-file',filenameanddir);


%% Load Spike variables, handle Spike events

%Finds, and loads variable into workspace as structure field name
channelnamestemp=fieldnames(channelnameoptions);

for i=1:length(channelnamestemp)
    temp=eval(['channelnameoptions.' char(channelnamestemp(i))])
    foundamatch=0;
    for n=1:length(temp)
        %Does it exist?
        for j=1:length(w)
            if strcmp(w(j).name,char(temp(n)))
                eval([channelnamestemp{i} '=filehandle.' char(temp(n))]);
                foundamatch=1;
                break
            end
        end
        if foundamatch
            break
        end
    end
end

% if Flow.interval<0.008
%     dt=0.008;
%     %Flow.values=resample(Flow.values,round(dt/Flow.interval),1);
%     Flow.values=downsample(Flow.values,round(dt/Flow.interval));
%     Flow.interval=0.008;
%     Flow.length=length(Flow.values);
% end

%Setup Time array
dt=0.008;
Fs=1/dt;
StartTime=0;
RecordingDuration = eval([channelnamestemp{Neventchannels+1} '.length *' channelnamestemp{Neventchannels+1} '.interval;']);
Time=(StartTime:dt:(StartTime+RecordingDuration-dt))';

clear foundamatch i j n w RecordingDuration StartTime
%% Make all channels the same length as Flow, resample if needed

clear channels_all
for i=Neventchannels+1:length(channelnamestemp)
    channels_all{i-Neventchannels}=channelnamestemp{i};
end

clear channels
for i=1:length(channels_all)
    if exist(channels_all{i},'var')
        channels{i}=channels_all{i};
    end
end

clear channels_all
%Add code here which resamples if channels are not the same sampling rate
for i=1:length(channels)
    if ~isempty(channels{i})
        dtnew=eval([channels{i} '.interval']);
        if dtnew~=dt
            disp(['resampling: ' channels{i} '.interval = ' num2str(eval([channels{i} '.interval']))]);
            if mod(dt/dtnew,1)~=0 %unfinished coding
                eval([channels{i} '.values = resample(' channels{i} '.values,round(1/dt),round(1/dtnew));' ]);
            else
                %eval([channels{i} '.values = resample(' channels{i} '.values,round(dt/' channels{i} '.interval),1);' ]);
                eval([channels{i} '.values = downsample(' channels{i} '.values,round(dt/dtnew));' ]);
            end
            eval([channels{i} '.interval = dt;']);
        end
    end
end

%Make channels the same length (sometimes these are off by up to 20 samples)
for i=1:length(channels)
    if ~isempty(channels{i})
        while length(eval([channels{i} '.values']))~=length(Time)
            if length(eval([channels{i} '.values']))<length(Time)
                eval([channels{i} '.values(end+1)=' channels{i} '.values(end);']); %add a sample to the end if the channel is too short by 1
            elseif length(eval([channels{i} '.values']))>length(Time)
                eval([channels{i} '.values(end)=[];']); %delete a sample from the end if the channel is too long by 1
            end
        end
    end
end

%Rename channel data: for example, Flow = Flow.values; clear Flow structure variable
channels_dt_original = NaN*ones(1,length(channels));
for i=1:length(channels)
    if ~isempty(channels{i})
        channels_dt_original(i) = eval([channels{i} '.interval;']);
        eval(['temp=' channels{i} '.values;']);
        eval([channels{i} '=temp;']);
        clear temp
    end
end
clear ans default dtnew foundamatch 
%% Position channel fix using text
% Fix Pos data if necessary
if ~exist('Position') 
    Position = 0*Time;
end
Position = Position*MultiplyPositionChannelby;
Position = round(Position);

textfilename=[directory filename(1:end-4) '_pos.txt'];
if exist(textfilename,'file')==2
    'found position file'
    %Position.values_original=Position.values;
    [col1, col2] = textread(textfilename,'%d%d%*[^\n]','delimiter',' ');
    for i=1:(length(col1)-1)
        Position(Time>=col1(i)&Time<col1(i+1))=col2(i);
    end
    Position(Time>=col1(end))=col2(end);
end
clear textfilename

%% Remove artifact in all signals using text files
if 0 %causes trouble!
    for j=1:length(channels)
        textfilename=[directory filename '_' channels{j} '_art.txt'];
        if exist(textfilename,'file')==2
            ['found ' channels{j} ' artifact file and removed artifact']
            %Position.values_original=Position.values;
            [col1, col2] = textread(textfilename,'%d%d%*[^\n]','delimiter',' ');
            for i=1:length(col1)
                %dttemp=eval([channels{j} '.interval']);
                lefti=round((col1(i)-Time(1))/dt+1);
                if lefti<1, lefti=1; end
                righti=round((col2(i)-Time(1))/dt+1);
                if righti>length(eval(channels{j})), righti=length(eval(channels{j})); end
                eval([channels{j} '(lefti:righti)=NaN;']);
            end
            Percent_removed = sum(isnan(eval(channels{j})))/length(eval(channels{j}))*100
            figure(5); plot(Time,eval([channels{j}]));
        end
    end
end

%% CPAP
if exist('Pmask')
minabsPmaskforCPAPoff=1.2;
[CPAP,CPAPoff] = CPAPonoroff(Pmask,Time,minabsPmaskforCPAPoff,1);
end
clear minabsPmaskforCPAPoff
%% Epochs / Staging channel
EpochsXHz=interp1((Epochs.times+15),double(Epochs.codes(:,1)'),Time,'nearest','extrap');

%% Events info
%Evts list is just called Evts
if exist('New_Evts','var'), Evts=New_Evts; end
if exist('New_Evts2','var'), Evts2=New_Evts2; end

backupEvts = Evts;

%add extra row of zeros with time = inf just in case there is an event without an end in the first of two concatenated sets of events.
Evts.codes(size(Evts.codes,1)+1,:)=[0 0 0 0];
Evts.times(size(Evts.times,1)+1,:)=Inf;
Evts.text=[Evts.text;['                                            ']];

% new code to copy into other versions of code:
while 1
    if Evts.codes(1,1)==0
        Evts.codes(1,:)=[];
        Evts.times(1)=[];
        Evts.length=Evts.length-1;
        Evts.text(1,:)=[];
    else
        break
    end
end

if exist('Evts2','var')
    while 1
        if Evts2.codes(1,1)==0
            Evts2.codes(1,:)=[];
            Evts2.times(1)=[];
            Evts2.length=Evts.length-1;
            Evts2.text(1,:)=[];
        else
            break
        end
    end
end
if exist('Evts2','var')&&exist('Evts','var')
    %combine events lists
    Evts.times=[Evts.times;Evts2.times];
    Evts.codes=[Evts.codes;Evts2.codes];
    Evts.text=[Evts.text;Evts2.text];
    Evts.length=Evts.length+Evts2.length;
end

%delete second row of zeros if there are two concurrent rows of zeros
for i=length(Evts.codes):-1:2
    if Evts.codes(i)==0&&Evts.codes(i-1)==0
        Evts.codes(i,:)=[];
        Evts.times(i)=[];
        Evts.text(i,:)=[];
        Evts.length=Evts.length-1;
        disp('delete_row_of_zeros');
    end
end

%add a row of zeros if there is a missing row of zeros. Assume end of event is the start of the next event (or end of file/Inf).
for i=length(Evts.codes):-1:1
    if i==length(Evts.codes)
        if Evts.codes(i)>0
            Evts.codes=[Evts.codes;[0 0 0 0]];
            Evts.times=[Evts.times;Inf];
            Evts.text=[Evts.text;['                                            ']];
            Evts.length=Evts.length+1;
            disp('added_row_of_zeros to the end');
        end
    end
    if i>1
        if Evts.codes(i)>0&&Evts.codes(i-1)>0
            Evts.codes=[Evts.codes(1:(i-1),:);[0 0 0 0];Evts.codes(i:end,:)];
            Evts.times=[Evts.times(1:i);Evts.times(i:end)];
            Evts.text=[Evts.text(1:(i-1),:);['                                            '];Evts.text(i:end,:)];
            Evts.length=Evts.length+1;
            disp('added_row_of_zeros');
        end
    end
end

% % Make list of event types
% clear EventTypeList_text EventTypeList_code tempi
% EventTypeList_text=[];
% for i=1:2:size(Evts.codes,1)
%     codematch=0;
%     for j=1:size(EventTypeList_text,1)
%         if strcmp(char(EventTypeList_text(j,:)),Evts.text(i,:))
%             'found code match'
%             codematch=1;
%             break
%         end
%     end
%     if codematch==0
%         tempi=size(EventTypeList_text,1);
%         EventTypeList_code(tempi+1)=Evts.codes(i);
%         EventTypeList_text(tempi+1,:)=Evts.text(i,:);
%     end
% end
% EventTypeList_code=EventTypeList_code'
% char(EventTypeList_text)
% 
% Evts_backup=Evts;
% 
% for i=1:size(Evts.codes,1)
%     if strncmp(Evts.text(i,:),'AR',2)
%         Evts.codes(i,1)=1;
%     end
% end

%remove zero rows:
Evts.codes(2:2:end,:)=[];
    Evts.codes(:,2:4)=[];
Evts.endtimes=Evts.times(2:2:end);
Evts.starttimes=Evts.times(1:2:end);
Evts.durations=Evts.endtimes-Evts.starttimes;
% Evts.timetoprevarousal=Evts.starttimes-[-Inf;Evts.endtimes(1:end-1)];
% Evts.timetonextarousal=[Evts.starttimes(2:end);Inf]-Evts.starttimes;

% Find sleep state and CPAP level for each event (arousal) %%faster code 2015-04-29
for i=1:length(Evts.starttimes)
    I=round((Evts.starttimes(i)-Time(1))/dt+1);
    Evts.epochs(i,:)=EpochsXHz(I);
    if exist('CPAP')
        Evts.CPAP(i,:)=CPAP(I);
    end
end

% Make arousals events in continuous time
EventsArXHz=zeros(length(Time),1); %Arousals.
for x=1:length(Evts.codes)
    if Evts.codes(x)==1
        lefti=round(Evts.starttimes(x)-Time(1))*Fs+1;
        righti=lefti+round(Evts.durations(x))*Fs;
        if righti==Inf
            righti=length(EventsArXHz);
        end
        %EventsArXHz(Time>=Evts.starttimes(x)&Time<Evts.endtimes(x))=Evts.codes(x); %%slower code because of search
        EventsArXHz(lefti:righti)=Evts.codes(x);
    end
end

% Make arousals events in continuous time
EventsRespXHz=zeros(length(Time),1); %Arousals.
for x=1:length(Evts.codes)
    if Evts.codes(x)>1
        lefti=round(Evts.starttimes(x)-Time(1))*Fs+1;
        righti=lefti+round(Evts.durations(x))*Fs;
        if righti==Inf
            righti=length(EventsArXHz);
        end
        %EventsArXHz(Time>=Evts.starttimes(x)&Time<Evts.endtimes(x))=Evts.codes(x); %%slower code because of search
        EventsArXHz(lefti:righti)=Evts.codes(x);
    end
end

clear lefti righti
clear EventTypeList_text EventTypeList_code Evts_backup Evts_backup1 Evts2 I backupEvts codematch foundamatch i j n temp tempi x

%channels{length(channels)+1}='Time';
channels{length(channels)+1}='EventsArXHz';
channels{length(channels)+1}='EventsRespXHz';
channels{length(channels)+1}='EpochsXHz';

