function [Position]=getPosSandman(Filenames,Time,StartTime,settings)

directory = char(Filenames(1,5));
fname = char(Filenames(1,2));

displaytext=['Importing position from Events List:' fname];

Sandmanfilename=[directory fname];
Sfileid=fopen(Sandmanfilename);
[Epochs] = textscan(Sfileid,'%u %s %s %f %f %f %f','whitespace','\t','headerlines',3);
fclose(Sfileid);
%         SandmandEpochs=Epochs(:,1:4); % first col: epoch #, 2: events list, 3:abs start time, 4:dur

EpochNumber = Epochs{1,1};
EventList = Epochs{1,2};
EventTime = Epochs{1,3};
EventDur=Epochs{1,4};
% tempPos=unique(EventList);

%% create epochs table
clear EventStartTime
for j = 1:length(EpochNumber)
    only_event_time = strsplit(EventTime{j});
    event_time = mod(datenum(only_event_time{1},'HH:MM:SS'),1)*86400;
    
              
    if event_time<43200; event_time=event_time+(86400/2); end % converting to 24hr format for anything less than 12am
    if event_time<StartTime % eg as in 1 am
        event_time=event_time+86400/2;
        EventStartTime(j,1)=event_time;
        %                 EventRelTime(j,1) = (event_time-StartTime);
    else
        EventStartTime(j,1)=event_time;
        %                 EventRelTime(j,1) = (event_time-StartTime);
    end
    
end
EpochNumber=double(EpochNumber);
EpochEventT=table(EpochNumber,EventList,EventStartTime,EventDur);


%% initializing position as supine if not mentioned
Position=ones(length(Time),1); 
% Position=ones(EpochNumber(end)*30,1); % fs=1hz;


%% get the actual position labels present in the event text files
PosCategoriesLabels={};
PosCategoriesCodes=[];

pattern={'Supine','Left','Right','Prone','Bathroom In','Upright'};
%using same codes as Profusion Mesa
% removing supine from the list since default is supine
ampcodes=[1 2 0 3 4 4];

for ii=1:length(pattern)
    StageTF = strcmpi(EventList,pattern{ii});
    if sum(StageTF)>0
        PosCategoriesLabels=[PosCategoriesLabels; unique(EventList(StageTF))];
        PosCategoriesCodes=[PosCategoriesCodes;ampcodes(ii)];
      
    end
end

%% get the event start and duration for position
EpochEventT.PosCodes=ones(size(EpochEventT,1),1);

if ~isempty(PosCategoriesLabels)
for ii=1:size(EpochEventT,1)
    for m=1:length(PosCategoriesLabels)
        if strcmpi(EpochEventT.EventList(ii),PosCategoriesLabels{m})
            EpochEventT.PosCodes(ii,:) = PosCategoriesCodes(m);
            % index will be ii
            StfromT=EpochEventT.EventStartTime(ii);
            IndofTime=find(Time<=StfromT,1,'last');
            PosDur=EpochEventT.EventDur(ii)*settings.Fs;
            Position(IndofTime:IndofTime+PosDur-1)=PosCategoriesCodes(m);
          
        end
    end
end
end

if length(Position)>length(Time)
    Position(length(Time)+1:end)=[];
end
figure(111111); clf(111111); plot(Time,Position)


