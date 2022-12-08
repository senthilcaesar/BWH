function [Evts,StartTime] = EventsAndHypnogram(system,Filenames,StartTime,EndTime,handletext,settings) %% Epochs and events
% global settings 
if ~exist('handletext')
    handletext=[];
end

% Each case must generate the following outputs:
% EventName,EventStart,EventDuration,EventLabelOptions,EventCodeOptions,EventExact1orStartsWith2,Hypnogram

% EventName,EventStart,EventDuration are intuitive; lists is the length of the number of events
% EventLabelOptions contains original name options (text) for each type of event
% EventCodeOptions contains matching codes for each each type of event, per
% the following:
% 1=Arousal|2=OA|3=CA|4=OH|5=MA|6=CH|7=OH2%|8=CH2%|9=OHnodesat|10=CHnodesat|11=sustainedFL|12=RERA

% EventExact1orStartsWith2 is a list of 1's or 2's depending on desired match type
% (e.g. use 2 to find Arousal if the EventName contains 'Arousal (spont)', 'Arousal (resp)', etc

% Hypnogram is a simple column list of 30-s Epoch codes

% starttimediff is optional. describes the time difference in s from first
% hypnogram start time and "StartTime". e.g. if the first Epoch starts 20 seconds after the start of file, then starttimediff=20.

% PositionT (optional) is a table that can be created from position event
% lists. Must contain "Time" (same units as the Time signal, describes start times of each position 'event') and "Codes"
% (must be Default codes i.e. use SystemPos='Default' to translate back to
% text positions: 1=Supine,2=Left,3=Right,4=Prone,5=Unknown,6=Upright

%%
%used only potentially for Fs in BrainRT (still needs work)
StartTime=StartTime;
EndTime=EndTime;

%default start time for scoring is start time of recording e.g. Time(1)
starttimediff=0;

% Empty Events List
EventStart=[];
EventDuration=[];
EventName=[];
EventLabelOptions=[];

%%
switch system
    case 'ProfusionXML_MrOs'
        %% ProfusionXml Hypnogram
        directory = char(Filenames(5));
        fname = char(Filenames(2));
        
        S=xml2struct([directory fname]);
        
        %Xml events
        displaytext=['Get Events data:' fname];
        if ~isempty(handletext), disp(displaytext); set(handletext,'String',displaytext); drawnow; end
        
        if ~isfield(S.PSGAnnotation.ScoredEvents,'ScoredEvent')
            displaytext=['Warning: No scored events'];
            if ~isempty(handletext), disp(displaytext); set(handletext,'String',displaytext); drawnow; end
        end
        
        Nevents = length(S.PSGAnnotation.ScoredEvents.ScoredEvent);
        
        clear EventName EventStart EventDuration EventCodesCategory
        for i=1:Nevents
            EventName{i} = S.PSGAnnotation.ScoredEvents.ScoredEvent{i}.EventConcept.Text;
            EventCodesCategory{i} = S.PSGAnnotation.ScoredEvents.ScoredEvent{i}.EventType.Text;
            EventStart(i) = StartTime + str2num(S.PSGAnnotation.ScoredEvents.ScoredEvent{i}.Start.Text);
            EventDuration(i) = str2num(S.PSGAnnotation.ScoredEvents.ScoredEvent{i}.Duration.Text);
        end %W=0,N1=1,N2=2,N3=3,R=5
        EventCodesListUnique = unique(EventName);
        EventCodesCategoryUnique = unique(EventCodesCategory);
        EventLabelOptions= {...
            'ASDA arousal|Arousal (ASDA)',...
            'ASDA',...
            'Central apnea|Central Apnea',...
            'Mixed apnea|Mixed Apnea',...
            'Obstructive apnea|Obstructive Apnea',...
            'Hypopnea|Hypopnea',...
            'Unsure|Unsure',...
            };
        EventCodeOptions=[1 1 3 5 2 4 4];
        EventExact1orStartsWith2 = [1 2 1 1 1 1 1];
        
        displaytext=['Get Hypnogram data: ' fname];
        if ~isempty(handletext), disp(displaytext); set(handletext,'String',displaytext); drawnow; end
        IsStage=NaN*zeros(Nevents,1);
        for i=1:Nevents
            IsStage(i) = strcmp(EventCodesCategory{i},'Stages|Stages');
        end %W=0,N1=1,N2=2,N3=3,R=5
        NepochChanges = sum(IsStage);
        EpochCodes=NaN*zeros(NepochChanges,1);
        Ilist=find(IsStage==1);
        for i=1:NepochChanges
            EpochCodes(i) = str2num(EventName{Ilist(i)}(end));
        end
        EpochDuration=EventDuration(Ilist)/30;
        Nepochs=sum(EpochDuration);
        Hypnogram=NaN*zeros(Nepochs,1);
        iEpochEnd=cumsum(EpochDuration);
        iEpochStart=[0 iEpochEnd(1:end-1)]+1;
        for i=1:length(EpochDuration)
            Hypnogram(iEpochStart(i):iEpochEnd(i))=EpochCodes(i);
        end
        %Convert to Terrill code
        Hypnogram(Hypnogram==0)=-1; %W=-1,N1=1,N2=2,N3=3,R=5
        Hypnogram(Hypnogram==5)=0;  %W=-1,R=0,N1=1,N2=2,N3=3
        Hypnogram=3-Hypnogram; %.MAT: 0=N3,1=N2,2=N1,3=R,4=W
        
        %end ProfusionXML_MrOs version
        
    case 'ProfusionXML'
        %% ProfusionXml Hypnogram
        directory = char(Filenames(5));
        fname = char(Filenames(2));
        
        %Xml events
        displaytext=['Get Events data:' fname];
        if ~isempty(handletext), disp(displaytext); set(handletext,'String',displaytext); drawnow; end
        
        S=xml2struct([directory fname]);
        
        if ~isfield(S.CMPStudyConfig.ScoredEvents,'ScoredEvent')
            displaytext=['Warning: No scored events'];
            if ~isempty(handletext), disp(displaytext); set(handletext,'String',displaytext); drawnow; end
        end
        
        Nevents = length(S.CMPStudyConfig.ScoredEvents.ScoredEvent);
        
        clear EventName EventStart EventDuration
        for i=1:Nevents
            EventName{i} = S.CMPStudyConfig.ScoredEvents.ScoredEvent{i}.Name.Text;
            EventStart(i) = StartTime + str2num(S.CMPStudyConfig.ScoredEvents.ScoredEvent{i}.Start.Text);
            EventDuration(i) = str2num(S.CMPStudyConfig.ScoredEvents.ScoredEvent{i}.Duration.Text);
        end %W=0,N1=1,N2=2,N3=3,R=5
        EventCodesListUnique = unique(EventName);
        
        
        EventLabelOptions= {...
            'Arousal',...
            'Central Apnea',...
            'Mixed Apnea',...
            'Obstructive Apnea',...
            'Hypopnea',... % unconvincing hypop (for MESA, MrOS)
            'Obstructive Hypopnea',... % numom2b
            'Central Hypopnea',...
            'RERA',...
            'Unsure',... % Special case for large cohort studies with >50% reduction in flow (e.g. MESA)
            };
        EventCodeOptions=[1 3 5 2 4 4 6 12 4];
        EventExact1orStartsWith2 = [2 1 1 1 1 1 1 1 1];
        
        
        displaytext=['Get Hypnogram data: ' fname];
        if ~isempty(handletext), disp(displaytext); set(handletext,'String',displaytext); drawnow; end
        Nepochs = length(S.CMPStudyConfig.SleepStages.SleepStage);
        Hypnogram=NaN*zeros(Nepochs,1);
        for i=1:Nepochs
            Hypnogram(i) = str2num(S.CMPStudyConfig.SleepStages.SleepStage{1,i}.Text);
        end %W=0,N1=1,N2=2,N3=3,R=5
        %Convert to Terrill code
        Hypnogram(Hypnogram==0)=-1; %W=-1,N1=1,N2=2,N3=3,R=5
        Hypnogram(Hypnogram==5)=0;  %W=-1,R=0,N1=1,N2=2,N3=3
        Hypnogram=3-Hypnogram; %.MAT: 0=N3,1=N2,2=N1,3=R,4=W
        Hypnogram(Hypnogram==-6)=NaN; % 'unscored' sleep in Profusion shows as -6
        % setting NaN here, then sets 8 below
        
        %end ProfusionXML
        
    case 'ProfusionXML_APN'
        %% ProfusionXml Hypnogram
        directory = char(Filenames(5));
        fname = char(Filenames(2));
        
        %Xml events
        displaytext=['Get Events data:' fname];
        if ~isempty(handletext), disp(displaytext); set(handletext,'String',displaytext); drawnow; end
        
        S=xml2struct([directory fname]);
        
        %S.CMPPSGSCOREDATA.SCOREDEVENTS.SCOREDEVENT{1, 1}.NAME
        
        if ~isfield(S.CMPPSGSCOREDATA.SCOREDEVENTS,'SCOREDEVENT')
            displaytext=['Warning: No scored events'];
            if ~isempty(handletext), disp(displaytext); set(handletext,'String',displaytext); drawnow; end
        end
        
        Nevents = length(S.CMPPSGSCOREDATA.SCOREDEVENTS.SCOREDEVENT);
        
        clear EventName EventStart EventDuration
        for i=1:Nevents
            EventName{i} = S.CMPPSGSCOREDATA.SCOREDEVENTS.SCOREDEVENT{i}.NAME.Text;
            EventStart(i) = StartTime + str2num(S.CMPPSGSCOREDATA.SCOREDEVENTS.SCOREDEVENT{i}.TIME.Text);
            EventDuration(i) = str2num(S.CMPPSGSCOREDATA.SCOREDEVENTS.SCOREDEVENT{i}.DURATION.Text);
        end %W=0,N1=1,N2=2,N3=3,R=5
        EventCodesListUnique = unique(EventName);
        
        EventLabelOptions= {...
            'Arousal',...
            'Central Apnea',...
            'Mixed Apnea',...
            'Obstructive Apnea',...
            'Hypopnea',... % unconvincing hypop (for MESA, MrOS)
            'Obstructive Hypopnea',... % numom2b
            'Central Hypopnea',...
            'RERA',...
            'Unsure',... % Special case for large cohort studies with >50% reduction in flow (e.g. MESA)
            };
        EventCodeOptions=[1 3 5 2 4 4 6 12 4];
        EventExact1orStartsWith2 = [2 1 1 1 1 1 1 1 1];
        
        displaytext=['Get Hypnogram data: ' fname];
        if ~isempty(handletext), disp(displaytext); set(handletext,'String',displaytext); drawnow; end %CMPPSGSCOREDATA.SCOREDEVENTS
        Nepochs = length(S.CMPPSGSCOREDATA.SLEEPSTAGES.SLEEPSTAGE);
        HypnogramT=nan(Nepochs,1);
        for i=1:Nepochs
            HypnogramT(i) = str2num(S.CMPPSGSCOREDATA.SLEEPSTAGES.SLEEPSTAGE{1,i}.Text);
        end %W=0,N1=1,N2=2,N3=3,R=5
        %Convert to Terrill code (0=N3,1=N2,2=N1,3=R,4=W)
        Hypnogram=nan(Nepochs,1);
        
        Hypnogram(HypnogramT==10)=4;
        Hypnogram(HypnogramT==1)=2;
        Hypnogram(HypnogramT==2)=1;
        Hypnogram(HypnogramT==3)=0;
        Hypnogram(HypnogramT==5)=3;
        
        
    case 'BrainRT'
        %%
        directory = char(Filenames(5));
        fname = char(Filenames(2));
        
        %Xml events
        displaytext=['Get Events data:' fname];
        if ~isempty(handletext), disp(displaytext); set(handletext,'String',displaytext); drawnow; end
        
        tic
        S=xml2struct([directory fname]); %slow
        toc
        
        if ~isfield(S.BrainRTResult.ResultElements.EventsResultElement.Events,'Event')
            displaytext=['Warning: No scored events'];
            if ~isempty(handletext), disp(displaytext); set(handletext,'String',displaytext); drawnow; end
        end
        
        %Keep only these sections of the xml struct, and clear S to save memory
        Sevents = S.BrainRTResult.ResultElements.EventsResultElement.Events.Event;
        Shypevents = S.BrainRTResult.ResultElements.HypnogramResultElement.ManualScoring.Events.Event;
        SLightsout = S.BrainRTResult.ResultElements.HypnogramResultElement.LightsOutEvents.Events.Event;
        
        clear S
        
        %make cell array table of events
        Nevents = length(Sevents);
        
        EventCellArray={''};
        formatIn = 'HH:MM:SS.FFF';
        clear EventName EventStart EventDuration
        for i=1:Nevents
            S1 = Sevents{i};
            EventCellArray{i,1}=S1.Type.Text;
            EventCellArray{i,2}=S1.SubType.Text;
            EventCellArray{i,3}=S1.Start.Text;
            EventCellArray{i,4}=S1.End.Text;
            try
                EventCellArray{i,5}=S1.CH1.Text;
                EventCellArray{i,6}=S1.CH2.Text;
            catch me
            end
            EventCellArray{i,7}=S1.Validated.Text;
            
            EventName{i} = [EventCellArray{i,1} '_' EventCellArray{i,2}];
            %duration
            temp = EventCellArray{i,3}(12:end-3);
            temp = mod(datenum(temp,formatIn),1)*86400;
            if temp<86400/2,temp=temp+86400; end
            EventStart(i)=temp;
            temp2 = EventCellArray{i,4}(12:end-3);
            temp2 = mod(datenum(temp2,formatIn),1)*86400;
            if temp2<86400/2,temp2=temp2+86400; end
            EventDuration(i) = temp2-temp;
            EventCellArray{i,8}=num2str(temp2-temp);
        end
        %
        %Possible event channel options, debug only:
        UniqueEventChannels = unique(EventCellArray(cellfun(@isempty,EventCellArray(:,5))==0,5));
        
        EventCodesListUnique = unique(EventName);
        EventLabelOptions= {...
            %                 '133_6',... %arousal
            %                 '133_2',... %arousal maybe
            %                 '133_3',... %arousal maybe
            '135_1',... %arousal maybe
            '129_3',... %event OA
            '129_2',... %event CA
            '129_4',... %event MA
            '129_6',... %event HypC
            '129_7',... %event HypO
            '129_5',... %event Hyp, assumed obstructive
            'Respiration_Obstructive apnea', ... % OA
            'Respiration_Obstructive hypopnea',... % HypO
            'Respiration_Central apnea', ... %CA
            'Respiration_Hypopnea'... %Hyp, assumed obstructive
            };
        EventCodeOptions=[1 2 3 5 6 4 4 2 4 3 4];
        EventExact1orStartsWith2 = [1 1 1 1 1 1 1 1 1 1 1];
        
        
        %Xml hypnogram (manual)
        %directory = char(Filenames(6));
        %fname = char(Filenames(3));
        
        displaytext=['Get Hypnogram data: ' fname];
        if ~isempty(handletext), disp(displaytext); set(handletext,'String',displaytext); drawnow; end
        
        %Shyp=xml2struct([directory filesep fname]); %slow %only if need second separate xml for sleep staging. not needed now.
        
        Nhypevents = length(Shypevents);
        
        HypCellArray={''};
        clear HypCodesList HypStart HypDuration
        for i=1:Nhypevents
            S1 = Shypevents{i};
            HypCellArray{i,1}=S1.Type.Text;
            HypCellArray{i,2}=S1.SubType.Text;
            HypCellArray{i,3}=S1.Start.Text;
            HypCellArray{i,4}=S1.End.Text;
            try
                HypCellArray{i,5}=S1.CH1.Text;
                HypCellArray{i,6}=S1.CH2.Text;
            catch me
            end
            
            HypCodesList{i} = [HypCellArray{i,1} '_' HypCellArray{i,2}];
            %duration
            temp = HypCellArray{i,3}(12:end-3);
            temp = mod(datenum(temp,formatIn),1)*86400;
            if temp<86400/2,temp=temp+86400; end
            HypStart(i,1)=temp;
            temp2 = HypCellArray{i,4}(12:end-3);
            temp2 = mod(datenum(temp2,formatIn),1)*86400;
            if temp2<86400/2,temp2=temp2+86400; end
            HypDuration(i,1) = temp2-temp;
            HypCellArray{i,7}=num2str(temp2-temp);
        end
        
        clear IsStage
        for i=1:Nhypevents
            IsStage(i,1) = strcmp(HypCellArray{i,1},'128');
        end %W=0,N1=1,N2=2,N3=3,R=5
        
        NepochChanges = sum(IsStage);
        EpochCodes=NaN*zeros(NepochChanges,1);
        Ilist=find(IsStage==1);
        for i=1:NepochChanges
            EpochCodes(i) = str2num(HypCellArray{Ilist(i),2});
        end
        iEpochDuration=(round(HypDuration(Ilist)/30));
        
        % disp('warning: code likely to fail here, Time not passed in and out of Events subroutine yet')
        % check time difference between lgihts on and end of EDF
        EndTimeMod = mod(datenum(datestr(EndTime/86400,'HH.MM.SS'),'HH.MM.SS'),1); % old error: 'HH.SS.MM'
        
        EndTimeLOn = mod(datenum(SLightsout.End.Text(12:end-3),formatIn),1);
        temp = [-1 0 +1];
        [~,i1]=min(abs(EndTimeMod - (EndTimeLOn - temp))); %handles the fact that there could be a difference of 1 day; temp(i1) adds/subtracts a day if needed (normally zero)
        EndTimeLOn = EndTimeLOn - temp(i1);
        
        EndTimeLastEvt = mod(datenum(datestr([EventStart(end)+EventDuration(end)]/86400,'HH.MM.SS'),'HH.MM.SS'),1); % old error: 'HH.SS.MM'
        [~,i1]=min(abs(EndTimeMod - (EndTimeLastEvt - temp))); %handles the fact that there could be a difference of 1 day; temp(i1) adds/subtracts a day if needed (normally zero)
        EndTimeLastEvt = EndTimeLastEvt - temp(i1);
        
        EndTimeLastHyp = mod(datenum(datestr([HypStart(end)+HypDuration(end)]/86400,'HH.MM.SS'),'HH.MM.SS'),1); % old error: 'HH.SS.MM'
        [~,i1]=min(abs(EndTimeMod - (EndTimeLastHyp - temp))); %handles the fact that there could be a difference of 1 day; temp(i1) adds/subtracts a day if needed (normally zero)
        EndTimeLastHyp = EndTimeLastHyp - temp(i1);
        
        ErrorEndTimeSec = 86400*(max([EndTimeLOn,EndTimeLastEvt,EndTimeLastHyp])-EndTimeMod); %e.g. Error will be +7200 s if EndTime is early by 2 hrs. Error will be small negative number if correct.
        if abs(ErrorEndTimeSec)>600
            disp(['Warning: Fixing StartTime/EndTime based on max(LightsOn/LastEvt/LastEpoch) minus endEDF. ErrorEndTimeSec=' num2str(ErrorEndTimeSec)]);
            StartTime = StartTime + ErrorEndTimeSec;
            EndTime = EndTime + ErrorEndTimeSec; %these need updating outside this function, also Time needs editing
        end
        
        iEpochStart=(1+round((HypStart(Ilist)-HypStart(1))/30)); %HypStart(Ilist(1))
        iEpochEnd=iEpochStart+iEpochDuration-1;
        checktemp=[iEpochStart iEpochDuration iEpochEnd EpochCodes];
        starttimediff = HypStart(1)-StartTime;
        
        %         if iEpochStart(1)<1 %this could still cause problems.
        %             %scoring starts before recording
        %             starttimediff = starttimediff + (1-iEpochStart(1))*30;
        %         end
        
        %starttimediffN=round((starttimediff/30))
        displaytext=['Start Time Hypnogram minus Signals: ' num2str(starttimediff) ' s'];
        if ~isempty(handletext), disp(displaytext); set(handletext,'String',displaytext); drawnow; end
        
        Nepochs=iEpochEnd(end);
        HypnogramX=NaN*zeros(Nepochs,1);
        for i=1:length(iEpochDuration)
            li=iEpochStart(i);
            ri=iEpochEnd(i);
            if li<1
                li=1;
            end
            if ri>length(HypnogramX)
                ri=length(HypnogramX);
            end
            HypnogramX(li:ri)=EpochCodes(i);
        end
        
        Hypnogram=NaN*zeros(Nepochs,1);
        %Convert to Terrill code (%.MAT: 0=N3,1=N2,2=N1,3=R,4=W)
        Hypnogram(HypnogramX==2)=4;
        Hypnogram(HypnogramX==201)=3;
        Hypnogram(HypnogramX==301)=2;
        Hypnogram(HypnogramX==302)=1;
        Hypnogram(HypnogramX==303)=0;
        Hypnogram(HypnogramX==304)=0;
        
        %         if 1 %set NaN, others, to 8
        %             Hypnogram(isnan(Hypnogram))=8; %no stage code
        %             Hypnogram(Hypnogram<0)=8; %no stage code
        %         end
        %         Hypnogram_t = StartTime+starttimediff+(0:30:(30*(length(Hypnogram)-1)))';
        %
        %             clear S
        %end BrainRT
        
    case 'Minerva'
        %%
        [~,~,eventstxt] = xlsread([directory filesep char(Filenames(2))],1,'A2:N5001'); %max number of events = 5k
        %delete excess
        temp = cell2mat(eventstxt(:,1));
        temp = find(isnan(temp),1);
        eventstxt(temp:end,:)=[];
        EventName = eventstxt(:,6);
        EventStart = cell2mat(eventstxt(:,2)); %warning all data are nearest second
        EventStart=EventStart*86400;
        EventStart(EventStart<43200)=EventStart(EventStart<43200)+86400;
        EventDuration = cell2mat(eventstxt(:,3));
        EventLabelOptions= {...
            'ASDA Arousal' ,...
            'Obs. Apnea',...
            'Cnt. Apnea',...
            'Cnt. Appea',...
            'Obs. Flow <70%',...
            'Obs. Flow <50%',...
            'Cnt. Flow <50%',...
            'Cnt. Flow <70%',...
            'Sustained FL'};
        EventCodeOptions=[1 2 3 3 4 4 6 6 11];
        EventExact1orStartsWith2 = [1 1 1 1 1 1 1 1 1];
        Hypnogram = xlsread([directory char(Filenames(3))],1,'B:B'); %[W=0, 1 2 3, R=5]
        %Convert to Terrill code
        Hypnogram(Hypnogram==0)=-1; %W=-1,N1=1,N2=2,N3=3,R=5
        Hypnogram(Hypnogram==5)=0;  %W=-1,R=0,N1=1,N2=2,N3=3
        Hypnogram=3-Hypnogram; %.MAT: 0=N3,1=N2,2=N1,3=R,4=W
        
    case 'Alice' %Legacy Alice software, exported to text
        %%
        directory = char(Filenames(6));
        fname=char(Filenames(3));
        displaytext=['Get Hypnogram data:' fname];
        if ~isempty(handletext), disp(displaytext); set(handletext,'String',displaytext); drawnow; end
        [num,~,~]=xlsread([directory fname],1,'A2:A86401');
        Hypnogram=15-downsample(num,30); %Alice: W=11,R=12,N1=13,N2=14,N3=15
        
        directory = char(Filenames(5));
        fname = char(Filenames(2));
        displaytext=['Get Events data:' fname];
        if ~isempty(handletext), disp(displaytext); set(handletext,'String',displaytext); drawnow; end
        Events_fname = [directory fname];
        
        [temp,~,~] = xlsread([Events_fname],'c1:c100000'); % Read the column defining the event start time
        Nlines = length(temp);
        
        [~,~,xlseventdataall] = xlsread([Events_fname],['a2:f' num2str(Nlines+1)]); % Read the event type column
        EventName = xlseventdataall(:,1);
        
        %deal with blank (NaN)
        for i=1:length(EventName)
            if isnan(EventName{i})
                EventName{i}='';
            end
        end
        %
        %             [~,~,EventName] = xlsread([Events_fname],'a:a'); % Read the event type column
        %             Nlines=length(EventName);
        %             EventName(1)=[]; % remove the header line of the data
        %
        EventStartList = xlseventdataall(:,3);
        EventStart = NaN*(ones(length(EventStartList),1));
        for i=1:length(EventStartList)
            EventStart(i) = EventStartList{i};
        end
        %[EventStart,~,~] = xlsread([Events_fname],'c:c'); % Read the column defining the event start time
        EventStart(EventStart<0.5)=EventStart(EventStart<0.5)+1;
        EventStart=EventStart*86400;
        
        EventDurationText = xlseventdataall(:,6);
        %[~,~,EventDurationText] = xlsread([Events_fname],['f2:f' num2str(Nlines)]); % Read the column defining the event duration data
        
        EventDuration=NaN*EventStart;
        for m=1:length(EventDurationText)
            if length(EventDurationText{m})>1
                Event_Duration{m}=sscanf((EventDurationText{m}),'%f'); %finds first number in text column
            else
                Event_Duration{m}=EventDurationText{m};
            end
            if isnumeric(Event_Duration{m})
                EventDuration(m)=Event_Duration{m};
            end
        end
        
        %             Evts.codes(2:2:end,:)=[];
        %             Evts.codes(:,2:4)=[];
        %             Evts.endtimes=Evts.times(2:2:end);
        %             Evts.starttimes=Evts.times(1:2:end);
        %             Evts.durations=Evts.endtimes-Evts.starttimes;
        %             Evts.text(2:2:end,:)=[];
        %             Evts.text(:,7:end)=[];
        %             EventDuration = Evts.durations;
        EventLabelOptions={'Central apnea','Central Apnea','Obstructive apnea','Obstructive Apnea',...
            'Hypopnea','u-arousal','Mixed apnea','Mixed Apnea','µ-arousal','Arousal'};
        EventCodeOptions=[3 3 2 2 4 1 5 5 1 1];
        EventExact1orStartsWith2 = [1 1 1 1 1 1 1 1 1 1];
        %             for i=1:length(EventStart)
        %                 EventName{i} = deblank(Evts.text(i,:));
        %             end
        
    case 'AliceTaiwan'
                %%
        directory = char(Filenames(6));
        fname=char(Filenames(3));
        displaytext=['Get Hypnogram data:' fname];
        if ~isempty(handletext), disp(displaytext); set(handletext,'String',displaytext); drawnow; end
        [num,~,~]=xlsread([directory fname],1,'B2:B86401');
        Hypnogram=15-downsample(num,30); %Alice: W=11,R=12,N1=13,N2=14,N3=15
        
        directory = char(Filenames(5));
        fname = char(Filenames(2));
        displaytext=['Get Events data:' fname];
        if ~isempty(handletext), disp(displaytext); set(handletext,'String',displaytext); drawnow; end
        Events_fname = [directory fname];
        
        [temp,~,~] = xlsread([Events_fname],'a6:a100000'); % Read the column defining the event start time
        Nlines = length(temp);
        
        [~,~,xlseventdataall] = xlsread([Events_fname],['a6:d' num2str(Nlines-1)]); % Read the event type column
        EventName = xlseventdataall(:,3);
        
        %deal with blank (NaN)
        for i=1:length(EventName)
            if isnan(EventName{i})
                EventName{i}='';
            end
        end
        %
        %             [~,~,EventName] = xlsread([Events_fname],'a:a'); % Read the event type column
        %             Nlines=length(EventName);
        %             EventName(1)=[]; % remove the header line of the data
        %
        EventStartList = xlseventdataall(:,1);
        EventStart = NaN*(ones(length(EventStartList),1));
        for i=1:length(EventStartList)
            EventStart(i) = EventStartList{i};
        end
        %[EventStart,~,~] = xlsread([Events_fname],'c:c'); % Read the column defining the event start time
        EventStart(EventStart<0.5)=EventStart(EventStart<0.5)+1;
        EventStart=EventStart*86400;
        
        EventDurationText = xlseventdataall(:,4);
        %[~,~,EventDurationText] = xlsread([Events_fname],['f2:f' num2str(Nlines)]); % Read the column defining the event duration data
        
        EventDuration=NaN*EventStart;
        for m=1:length(EventDurationText)
            if length(EventDurationText{m})>1
                Event_Duration{m}=sscanf((EventDurationText{m}),'%f'); %finds first number in text column
            else
                Event_Duration{m}=EventDurationText{m};
            end
            if isnumeric(Event_Duration{m})
                EventDuration(m)=Event_Duration{m};
            end
        end
        
        %             Evts.codes(2:2:end,:)=[];
        %             Evts.codes(:,2:4)=[];
        %             Evts.endtimes=Evts.times(2:2:end);
        %             Evts.starttimes=Evts.times(1:2:end);
        %             Evts.durations=Evts.endtimes-Evts.starttimes;
        %             Evts.text(2:2:end,:)=[];
        %             Evts.text(:,7:end)=[];
        %             EventDuration = Evts.durations;
        EventLabelOptions={'Central apnea','Central Apnea','Obstructive apnea','Obstructive Apnea',...
            'Hypopnea','u-arousal','Mixed apnea','Mixed Apnea','µ-arousal','Arousal',[char(163) 'g-arousal']};
        EventCodeOptions=[3 3 2 2 4 1 5 5 1 1 1];
        EventExact1orStartsWith2 = [1 1 1 1 1 1 1 1 1 1 1];
        %             for i=1:length(EventStart)
        %                 EventName{i} = deblank(Evts.text(i,:));
        %             end
    case 'AliceG3Taiwan'
                %%
        directory = char(Filenames(6));
        fname=char(Filenames(3));
        displaytext=['Get Hypnogram data:' fname];
        if ~isempty(handletext), disp(displaytext); set(handletext,'String',displaytext); drawnow; end
        [~,txt,~]=xlsread([directory fname],1,'C2:C86401');
%         Hypnogram=15-downsample(num,30); %Alice: W=11,R=12,N1=13,N2=14,N3=15

        StageTextOptions = {'WK','N1','N2','N3','REM','NotScored'};
        StageTextCodes = [4,2,1,0,3,8];
        Hypnogram = ones(length(txt),1)*8;
        for jj = 1:length(StageTextOptions)
            epochidx = strcmp(txt, StageTextOptions{jj});
            Hypnogram(epochidx) = StageTextCodes(jj);
        end
        
        directory = char(Filenames(5));
        fname = char(Filenames(2));
        displaytext=['Get Events data:' fname];
        if ~isempty(handletext), disp(displaytext); set(handletext,'String',displaytext); drawnow; end
        Events_fname = [directory fname];
        
        [temp,~,~] = xlsread([Events_fname],'c6:c100000'); % Read the column defining the event start time
        Nlines = length(temp);
        
        [~,~,xlseventdataall] = xlsread([Events_fname],['a2:e' num2str(Nlines-1)]); % Read the event type column
        EventName = xlseventdataall(:,1);
        
        %deal with blank (NaN)
        for i=1:length(EventName)
            if isnan(EventName{i})
                EventName{i}='';
            end
        end
        %
        %             [~,~,EventName] = xlsread([Events_fname],'a:a'); % Read the event type column
        %             Nlines=length(EventName);
        %             EventName(1)=[]; % remove the header line of the data
        %
        EventStartList = xlseventdataall(:,3);
        EventStart = NaN*(ones(length(EventStartList),1));
        for i=1:length(EventStartList)
            EventStart(i) = EventStartList{i};
        end
        %[EventStart,~,~] = xlsread([Events_fname],'c:c'); % Read the column defining the event start time
        EventStart(EventStart<0.5)=EventStart(EventStart<0.5)+1;
        EventStart=EventStart*86400;
        
        EventDurationText = xlseventdataall(:,5);
        %[~,~,EventDurationText] = xlsread([Events_fname],['f2:f' num2str(Nlines)]); % Read the column defining the event duration data
        
        EventDuration=NaN*EventStart;
        for m=1:length(EventDurationText)
            if length(EventDurationText{m})>1
                Event_Duration{m}=sscanf((EventDurationText{m}),'%f'); %finds first number in text column
            else
                Event_Duration{m}=EventDurationText{m};
            end
            if isnumeric(Event_Duration{m})
                EventDuration(m)=Event_Duration{m};
            end
        end
        
        %             Evts.codes(2:2:end,:)=[];
        %             Evts.codes(:,2:4)=[];
        %             Evts.endtimes=Evts.times(2:2:end);
        %             Evts.starttimes=Evts.times(1:2:end);
        %             Evts.durations=Evts.endtimes-Evts.starttimes;
        %             Evts.text(2:2:end,:)=[];
        %             Evts.text(:,7:end)=[];
        %             EventDuration = Evts.durations;
        EventLabelOptions={'Central apnea','Central Apnea','Obstructive apnea','Obstructive Apnea',...
            'Hypopnea','u-arousal','Mixed apnea','Mixed Apnea','µ-arousal','Arousal'};
        EventCodeOptions=[3 3 2 2 4 1 5 5 1 1];
        EventExact1orStartsWith2 = [1 1 1 1 1 1 1 1 1 1];
        %             for i=1:length(EventStart)
        %                 EventName{i} = deblank(Evts.text(i,:));
        %             end
                
    case 'AliceG3'  %check
        %%
        directory = char(Filenames(6));
        fname=char(Filenames(3));
        
        S=xml2struct([directory '\' fname]);
        %Xml Hypnogram
        displaytext=['Get Hypnogram data: ' fname];
        if ~isempty(handletext), disp(displaytext); set(handletext,'String',displaytext); drawnow; end
        
        
        StageTextOptions = {'Wake','NonREM1','NonREM2','NonREM3','REM','NotScored'};
        StageTextCodes = [4,2,1,0,3,8];
        try
            if 1 %user staging
                SScoringData = S.PatientStudy.ScoringData.StagingData.UserStaging;
            else %machine staging
                SScoringData = S.PatientStudy.ScoringData.StagingData.MachineStaging;
            end
            
            Ntransitions = length(SScoringData.NeuroAdultAASMStaging.Stage);
            clear Transitions
            for i=1:Ntransitions
                Transitions{i}.Text = SScoringData.NeuroAdultAASMStaging.Stage{i}.Attributes.Type;
                Transitions{i}.Time = str2num(SScoringData.NeuroAdultAASMStaging.Stage{i}.Attributes.Start)/30 +1;
                Transitions{i}.Code = -1; %Default is unknown
            end
            for i=1:Ntransitions
                for m=1:length(StageTextOptions)
                    if strcmp(StageTextOptions{m},Transitions{i}.Text)==1
                        %[StageTextOptions{m} ' ' num2str(StageTextCodes(m))]
                        Transitions{i}.Code = StageTextCodes(m);
                    end
                end
            end
        catch me
            disp('Error: No staging found in rml file');
            Ntransitions=0;
        end
        RecordingDuration = str2num(S.PatientStudy.Acquisition.Sessions.Session.Duration.Text);
        Nepochs = ceil(RecordingDuration/30);
        
        clear hypnogram
        Hypnogram=NaN*zeros(Nepochs,1);
        for i=1:Ntransitions
            indexi = Transitions{i}.Time;
            Hypnogram(indexi:end) = Transitions{i}.Code;
        end
        
        %hold('off'); stairs(hypnogram_t,hypnogram); hold('on'); plot(Pnasal_Time,EpochsXHz,'r:');
        
        %Xml/rml events [tested and working]
        displaytext=['Get Events data:' fname];
        if ~isempty(handletext), disp(displaytext); set(handletext,'String',displaytext); drawnow; end
        
        %TempEvents = S.PatientStudy.ScoringData.Events.Event;
        %EventTypesOfInterest = {'ObstructiveApnea','Hypopnea','CentralApnea','MixedApnea','Arousal','RelativeDesaturation'};
        
        Nevents = length(S.PatientStudy.ScoringData.Events.Event);
        clear EventStart EventName EventDuration
        for i=1:Nevents
            EventName{i,1} = S.PatientStudy.ScoringData.Events.Event{i}.Attributes.Type;
            EventStart(i,1) = StartTime + str2num(S.PatientStudy.ScoringData.Events.Event{i}.Attributes.Start);
            EventDuration(i,1) = str2num(S.PatientStudy.ScoringData.Events.Event{i}.Attributes.Duration);
        end %W=0,N1=1,N2=2,N3=3,R=5
        
        %Add Events data in additional colums of the matrix:
        %**************************************************************************
        EventLabelOptions={'CentralApnea','ObstructiveApnea','MixedApnea','Hypopnea','ObstructiveHypopnea','Arousal','CentralHypopnea','RERA'};
        EventCodeOptions=[3 2 5 4 4 1 6 12];
        EventExact1orStartsWith2 = [1 1 1 1 1 2 1 1];
        
        if isfield(settings,'RERAmeansOH') && settings.RERAmeansOH %will be a new option below to swap things to OH by choice
            disp('Scored "RERAs" taken to be obstructive hypopneas [USCD clinic]')
            EventLabelOptions=[EventLabelOptions,{'RERA'}];
            EventCodeOptions=[EventCodeOptions 4];
            EventExact1orStartsWith2 = [EventExact1orStartsWith2 1];
        end
        
    case 'AliceG3Machine'  %code takes machine scoring, overwrites with user staging when present
        %%
        directory = char(Filenames(6));
        fname=char(Filenames(3));
        
        S=xml2struct([directory '\' fname]);
        %Xml Hypnogram
        displaytext=['Get Hypnogram data: ' fname];
        if ~isempty(handletext), disp(displaytext); set(handletext,'String',displaytext); drawnow; end
        
        
        StageTextOptions = {'Wake','NonREM1','NonREM2','NonREM3','REM','NotScored'};
        StageTextCodes = [4,2,1,0,3,8];
        
        
        %machine staging:
        SScoringData = S.PatientStudy.ScoringData.StagingData.MachineStaging;
        Ntransitions = length(SScoringData.NeuroAdultAASMStaging.Stage);
        clear Transitions
        for i=1:Ntransitions
            Transitions{i}.Text = SScoringData.NeuroAdultAASMStaging.Stage{i}.Attributes.Type;
            Transitions{i}.Time = str2num(SScoringData.NeuroAdultAASMStaging.Stage{i}.Attributes.Start)/30 +1;
            Transitions{i}.Code = -1; %Default is unknown
        end
        for i=1:Ntransitions
            for m=1:length(StageTextOptions)
                if strcmp(StageTextOptions{m},Transitions{i}.Text)==1
                    %[StageTextOptions{m} ' ' num2str(StageTextCodes(m))]
                    Transitions{i}.Code = StageTextCodes(m);
                end
            end
        end
        RecordingDuration = str2num(S.PatientStudy.Acquisition.Sessions.Session.Duration.Text);
        Nepochs = ceil(RecordingDuration/30);
        clear hypnogram
        Hypnogram=NaN*zeros(Nepochs,1);
        for i=1:Ntransitions
            indexi = Transitions{i}.Time;
            Hypnogram(indexi:end) = Transitions{i}.Code;
        end
        
        try
            if 1 %find user staging and overwrite when epochs are not unscored
                %user staging:
                SScoringDataUser = S.PatientStudy.ScoringData.StagingData.UserStaging;
                Ntransitions = length(SScoringDataUser.NeuroAdultAASMStaging.Stage);
                clear TransitionsUser
                for i=1:Ntransitions
                    TransitionsUser{i}.Text = SScoringDataUser.NeuroAdultAASMStaging.Stage{i}.Attributes.Type;
                    TransitionsUser{i}.Time = str2num(SScoringDataUser.NeuroAdultAASMStaging.Stage{i}.Attributes.Start)/30 +1;
                    TransitionsUser{i}.Code = -1; %Default is unknown
                end
                for i=1:Ntransitions
                    for m=1:length(StageTextOptions)
                        if strcmp(StageTextOptions{m},TransitionsUser{i}.Text)==1
                            %[StageTextOptions{m} ' ' num2str(StageTextCodes(m))]
                            TransitionsUser{i}.Code = StageTextCodes(m);
                        end
                    end
                end
                clear hypnogramUser
                HypnogramUser=NaN*zeros(Nepochs,1);
                for i=1:Ntransitions
                    indexi = TransitionsUser{i}.Time;
                    HypnogramUser(indexi:end) = TransitionsUser{i}.Code;
                end
                HypnogramMachine=Hypnogram;
                I = HypnogramUser~=8;
                Hypnogram(I) = HypnogramUser(I); %overwrite
                HypnogramCompare = table(Hypnogram,HypnogramMachine,HypnogramUser); %for view only
            end
        catch me
            disp('No User Staging Found, Machine Staging is not edited');
        end
        
        %hold('off'); stairs(hypnogram_t,hypnogram); hold('on'); plot(Pnasal_Time,EpochsXHz,'r:');
        
        %Xml/rml events [tested and working]
        displaytext=['Get Events data:' fname];
        if ~isempty(handletext), disp(displaytext); set(handletext,'String',displaytext); drawnow; end
        
        TempEvents = S.PatientStudy.ScoringData.Events.Event;
        EventTypesOfInterest = {'ObstructiveApnea','Hypopnea','CentralHypopnea', 'ObstructiveHypopnea', 'CentralApnea','MixedApnea','Arousal','RelativeDesaturation'};
        
        Nevents = length(S.PatientStudy.ScoringData.Events.Event);
        clear EventStart EventName EventDuration
        for i=1:Nevents
            EventName{i,1} = S.PatientStudy.ScoringData.Events.Event{i}.Attributes.Type;
            EventStart(i,1) = StartTime + str2num(S.PatientStudy.ScoringData.Events.Event{i}.Attributes.Start);
            EventDuration(i,1) = str2num(S.PatientStudy.ScoringData.Events.Event{i}.Attributes.Duration);
        end %W=0,N1=1,N2=2,N3=3,R=5
        
        %Add Events data in additional colums of the matrix:
        %**************************************************************************
        EventLabelOptions={'CentralApnea','ObstructiveApnea','MixedApnea','Hypopnea','ObstructiveHypopnea', 'CentralHypopnea', 'Arousal'};
        EventCodeOptions=[3 2 5 4 4 6 1];
        EventExact1orStartsWith2 = [1 1 1 1 1 1 2];
        
        if isfield(settings,'RERAmeansOH') && settings.RERAmeansOH
            disp('Scored "RERAs" taken to be obstructive hypopneas [USCD clinic]')
            EventLabelOptions=[EventLabelOptions,{'RERA'}];
            EventCodeOptions=[EventCodeOptions 4];
            EventExact1orStartsWith2 = [EventExact1orStartsWith2 1];
        end
        
        if 0 %this currently breaks, and is unused [6/7/2021]
            PositionT = table(EventName,EventStart,EventDuration);
            PositionT.Properties.VariableNames = {'EventName','Time','EventDuration'};
            clear Position IsPosition
            for i=1:Nevents
                if isfield(S.PatientStudy.ScoringData.Events.Event{i},'BodyPosition')
                    Position{i,1} = S.PatientStudy.ScoringData.Events.Event{i}.BodyPosition.Text;
                    IsPosition(i,1)=1;
                else
                    Position{i,1} = [];
                    IsPosition(i,1)=0;
                end
            end
            PositionT.Position=Position;
            PositionT = PositionT(IsPosition==1,:);
            PositionT.EventEnd = PositionT.Time+ PositionT.EventDuration;
            
            
            PosCodes = {'Supine','Left','Right','Prone','Unknown','Upright'}; %not detected: ,'Prone','Unknown','Upright' but leaving these here as placeholders
            NewPosCodes = {'Supine','Left','Right','Prone','Unknown','Upright'}; %from database, hardcoded
            PosCode = [1 2 3 4 5 6]; %based on REMlogic and PUP standard
            
            
            %LTM siesta alice original signal data were supine=0,right=2,left=3
            I = sum(PositionT.Position == string(PosCodes),2)>0; %this might remove positions that are not yet coded
            if sum(I==0)>0
                disp(strjoin(['warning: need to code any uncoded positions:' unique(string(PositionT.Position(:)))']));
            end
            PositionT = PositionT(I,:);
            
            PositionT.PositionOrig = PositionT.Position;
            PositionT.Codes = nan(height(PositionT),1);
            for j=1:length(PosCodes)
                for i=1:height(PositionT)
                    if PositionT.Position(i)==string(PosCodes(j))
                        PositionT.Position{i}=NewPosCodes{j};
                        PositionT.Codes(i)=PosCode(j);
                    end
                end
            end
            
            %Fill position gaps
            DefaultRow = PositionT(1,:);
            DefaultRow.EventName={'BodyPositionDefault'};
            DefaultRow.Time=NaN;
            DefaultRow.EventEnd=NaN;
            DefaultRow.EventDuration=NaN;
            DefaultRow.Position={'Left'};
            DefaultRow.PositionOrig={'Left'};
            DefaultRow.Codes=2;
            
            %Make NewPosT, using interwoven lateral rows
            Irange = 2*[(1:height(PositionT))'];
            NewPosT = repmat(DefaultRow,1+2*height(PositionT),1);
            NewPosT(Irange,:)   = PositionT;
            NewPosT.Time(1)=StartTime;
            NewPosT.EventEnd(end)=EndTime;
            temp = NewPosT.EventEnd(1:end-1); %end event times
            temp2 = NewPosT.Time(2:end); %start next event times (should be equal to end event times, copy these across via cross match when empty)
            temp(isnan(temp))=temp2(isnan(temp));
            NewPosT.EventEnd(1:end-1)=temp;
            temp2(isnan(temp2))=temp(isnan(temp2));
            NewPosT.Time(2:end)=temp2;
            
            NewPosT.EventDuration=NewPosT.EventEnd-NewPosT.Time;
            NewPosT(NewPosT.EventDuration<=0,:)=[]; %remove 0-sec events
            
            PositionT=NewPosT; %overwrite
            clear NewPosT;
        end
        
        %%
        %         clear BPtext BPstart
        %         temp = S.PatientStudy.BodyPositionState.BodyPositionItem;
        %         for i=1:length(temp)
        %             BPtext{i,1}=temp{i}.Attributes.Position;
        %             BPstart(i,1)=str2double(temp{i}.Attributes.Start);
        %         end
        %         PosTtemp = table(BPtext,BPstart);
        %
        
        
    case 'Spike'
        %%
        
        filehandle = matfile([Filenames{5} Filenames{2}]);
        
        
        w = whos(filehandle);
        
        %load events
        echannelnameoptions.Evts={'Evts','New_Evts'};
        echannelnameoptions.Evts2={'Evts2','New_Evts2'};
        echannelnameoptions.Evts3={'Evts3','New_Evts3'};
        echannelnameoptions.Epochs={'Epochs'};
        
        channelnamestemp=fieldnames(echannelnameoptions);
        
        for i=1:length(channelnamestemp)
            temp=eval(['echannelnameoptions.' char(channelnamestemp(i))]);
            foundamatch=0;
            for nn=1:length(temp)
                %Does it exist?
                for j=1:length(w)
                    if strcmp(w(j).name,char(temp(nn)))
                        eval([channelnamestemp{i} '=filehandle.' char(temp(nn)),';']);
                        foundamatch=1;
                        break
                    end
                end
                if foundamatch
                    break
                end
            end
        end
        
        %process events
        %add extra row of zeros with time = inf just in case there is an event without an end in the first of two concatenated sets of events.
        Evts.codes(size(Evts.codes,1)+1,:)=[0 0 0 0];
        Evts.times(size(Evts.times,1)+1,:)=Inf;
        blanktext = repmat(' ',1,size(Evts.text,2));
        Evts.text=[Evts.text;[blanktext]];
        
        
        %delete all start rows without codes at start of list
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
        %repeat for Evts2 separately
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
        %repeat for Evts3 separately
        if exist('Evts3','var')
            while 1
                if Evts3.codes(1,1)==0
                    Evts3.codes(1,:)=[];
                    Evts3.times(1)=[];
                    Evts3.length=Evts.length-1;
                    Evts3.text(1,:)=[];
                else
                    break
                end
            end
        end
        
        for i=1:size(Evts.codes,1)
            Evts.source(i,1) = 1;
        end
        Evts.text(:,49:end)=[]; %SS added 1/25 to handle larger event strings appearing
        if exist('Evts2','var')
            for i=1:size(Evts2.codes,1)
                Evts2.source(i,1) = 2;
            end
            Evts2.text(:,49:end)=[]; %SS added 1/25 to handle larger event strings appearing
        end
        if exist('Evts3','var')
            for i=1:size(Evts3.codes,1)
                Evts3.source(i,1) = 3;
            end
            Evts3.text(:,49:end)=[]; %SS added 1/25 to handle larger event strings appearing
        end
        
        
        if exist('Evts2','var')&&exist('Evts','var')
            %combine events lists
            Evts.times=[Evts.times;Evts2.times];
            Evts.codes=[Evts.codes;Evts2.codes];
            temp = min([size(Evts.text,2) size(Evts2.text,2)]); %SS added 1/25
            Evts.text=[Evts.text(:,1:temp);Evts2.text(:,1:temp)]; %SS updated 1/25
            Evts.length=Evts.length+Evts2.length;
            Evts.source=[Evts.source;Evts2.source];
        end
        
        if exist('Evts3','var')&&exist('Evts','var')
            %combine events lists
            Evts.times=[Evts.times;Evts3.times];
            Evts.codes=[Evts.codes;Evts3.codes];
            temp = min([size(Evts.text,2) size(Evts3.text,2)]);  %SS added 1/25
            Evts.text=[Evts.text(:,1:temp);Evts3.text(:,1:temp)]; %SS updated 1/25
            %Evts.text=[Evts.text;Evts3.text];
            Evts.length=Evts.length+Evts3.length;
            Evts.source=[Evts.source;Evts3.source];
        end
        
        
        %delete second row of zeros if there are two concurrent rows of zeros
        for i=length(Evts.codes):-1:2
            if Evts.codes(i)==0&&Evts.codes(i-1)==0
                Evts.codes(i,:)=[];
                Evts.times(i)=[];
                Evts.source(i)=[];
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
                    Evts.text=[Evts.text;[blanktext]];
                    Evts.length=Evts.length+1;
                    Evts.source(i)=Evts.source(i-1);
                    disp('added_row_of_zeros to the end');
                end
            end
            if i>1
                if Evts.codes(i)>0&&Evts.codes(i-1)>0
                    Evts.codes=[Evts.codes(1:(i-1),:);[0 0 0 0];Evts.codes(i:end,:)];
                    Evts.times=[Evts.times(1:i);Evts.times(i:end)];
                    Evts.text=[Evts.text(1:(i-1),:);[blanktext];Evts.text(i:end,:)];
                    Evts.length=Evts.length+1;
                    Evts.source=[Evts.source(1:i);Evts.source(i:end)];
                    disp('added_row_of_zeros');
                end
            end
        end
        
        %remove zero rows:
        Evts.codes(2:2:end,:)=[];
        Evts.codes(:,2:4)=[];
        Evts.endtimes=Evts.times(2:2:end);
        Evts.starttimes=Evts.times(1:2:end);
        Evts.durations=Evts.endtimes-Evts.starttimes;
        Evts.text(2:2:end,:)=[];
        Evts.text(:,7:end)=[];
        Evts.source=Evts.source(1:2:end);
        
        % check for overlapping events and remove
        %Intention is to prioritize the events list that is shorter i.e. represents corrections
        %This corrections Evts list should be called Evts3, and if not there will be a warning given.
        overlappingeventsoption=1;
        if overlappingeventsoption
            disp('overlappingeventsoption');
            evtssources = [sum(Evts.source==1) sum(Evts.source==2) sum(Evts.source==3)];
            evtsrespsources = [sum(Evts.source==1&Evts.codes~=1) sum(Evts.source==2&Evts.codes~=1) sum(Evts.source==3&Evts.codes~=1)];
            evtsarsources = [sum(Evts.source==1&Evts.codes==1) sum(Evts.source==2&Evts.codes==1) sum(Evts.source==3&Evts.codes==1)];
            
            evtsrespsources2 = evtsrespsources;
            evtsrespsources2(evtsrespsources2==0)=NaN;
            evtsrespsources2(evtsarsources>0)=NaN;
            
            if sum(~isnan(evtsrespsources2))>1
                if 0
                    [~,prioritizesource] = min(evtsrespsources2);
                else
                    prioritizesource = find(~isnan(evtsrespsources2),1,'last');
                end
                % if prioritizesource~=3
                %     disp('warning, not prioritizing Evts3');
                % end
                
                lengthevtscodes=length(Evts.codes);
                for i=length(Evts.codes):-1:1
                    for j=length(Evts.codes):-1:1
                        if i>length(Evts.codes)||j>length(Evts.codes)||i==j||Evts.codes(i)==1||Evts.codes(j)==1
                            continue
                        end
                        a=Evts.starttimes(i);
                        b=Evts.endtimes(i);
                        A=Evts.starttimes(j);
                        B=Evts.endtimes(j);
                        if ~(A>a&&A<b || B>a&&B<b || a>A&&a<B || b>A&&b<B)
                            %no overlap
                            continue
                        end
                        if Evts.source(j)~=prioritizesource
                            disp(['removing event: detected overlap i=' num2str(i) '[' num2str(Evts.source(i)), '], j=' num2str(j) '[' num2str(Evts.source(j)) '], removing j ' num2str(Evts.codes(j)) ', keeping i ' num2str(Evts.codes(i))]);
                            Evts.source(j)=[];
                            Evts.codes(j)=[];
                            Evts.endtimes(j)=[];
                            Evts.starttimes(j)=[];
                            Evts.durations(j)=[];
                            Evts.text(j,:)=[];
                        else
                            disp(['removing event: i=' num2str(i) '[' num2str(Evts.source(i)), '], j=' num2str(j) '[' num2str(Evts.source(j)) '], removing i: ' num2str(Evts.codes(i)) ', keeping j ' num2str(Evts.codes(j))]);
                            Evts.source(i)=[];
                            Evts.codes(i)=[];
                            Evts.endtimes(i)=[];
                            Evts.starttimes(i)=[];
                            Evts.durations(i)=[];
                            Evts.text(i,:)=[];
                        end
                    end
                end
            end
        end
        %  EventStart=Evts.starttimes + StartTime;
        EventDuration = Evts.durations;
        EventStart=Evts.starttimes; % Made by DV on 6/30/21
        
        EventLabelOptions= {...
            'AR' ,...
            'ApO',...
            'ApC',...
            'HypO',...
            'M',...
            'HypC',...
            'HypO2%',...
            'HypC2%',...
            'HypOx',...
            'HypCx',...
            'Ap-C',... %older labels
            'H',... %older labels
            'Hyp-C',... %older labels
            'ap-O',... %older labels
            };
        
        EventCodeOptions=[1 2 3 4 5 6 7 8 9 10 3 4 6 2];
        EventExact1orStartsWith2 = [1 1 1 1 1 1 1 1 1 1 1 1 1 1];
        for i=1:length(EventStart)
            EventName{i} = deblank(Evts.text(i,:));
        end
        
        StartTimeEpochs = Epochs.times(1);
        starttimediff = StartTimeEpochs;    %added by SS, 7/19 to handle start scoring not at start file.
        if StartTimeEpochs>0
            disp(['StartEpochTime=' num2str(starttimediff)]);
        end
        Hypnogram=double(Epochs.codes(:,1));
        
        Hypnogram(Hypnogram==0)=-1; %W=-1,N1=1,N2=2,N3=3,R=5
        Hypnogram(Hypnogram==5)=0;  %W=-1,R=0,N1=1,N2=2,N3=3
        Hypnogram=3-Hypnogram; %.MAT: 0=N3,1=N2,2=N1,3=R,4=W
        Hypnogram(Hypnogram==-3)=NaN;
        
    case 'SpikeDise'
        %%
        
        filehandle = matfile([Filenames{5} Filenames{2}]);
        w = whos(filehandle);
        
        %load events
        echannelnameoptions.Evts={'VOTE'};
        echannelnameoptions.Evts2={'MouthOpenClosed'};
        echannelnameoptions.Evts3={'SideSleep'};
        %         echannelnameoptions.Epochs={'Epochs'};
        
        channelnamestemp=fieldnames(echannelnameoptions);
        
        for i=1:length(channelnamestemp)
            temp=echannelnameoptions.(channelnamestemp{i});
            foundamatch=0;
            for nn=1:length(temp)
                %Does it exist?
                for j=1:length(w)
                    if strcmp(w(j).name,char(temp(nn)))
                        %                         channelnamestemp{i}=filehandle.(temp{nn});
                        eval([channelnamestemp{i} '=filehandle.' char(temp(nn)),';']);
                        foundamatch=1;
                        break
                    end
                end
                if foundamatch
                    break
                end
            end
        end
        
        % added this in for early days of DISE before labelling occured
        % DELETE AFTER LABELLING
        if ~exist('Evts', 'var')
            Evts=NaN;
            Hypnogram=NaN;
        else
            %process events
            %add extra row of zeros with time = inf just in case there is an event without an end in the first of two concatenated sets of events.
            Evts.codes(size(Evts.codes,1)+1,:)=[0 0 0 0];
            Evts.times(size(Evts.times,1)+1,:)=Inf;
            blanktext = repmat(' ',1,size(Evts.text,2));
            Evts.text=[Evts.text;[blanktext]];
            
            %delete all start rows without codes at start of list
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
            %repeat for Evts2 separately
            if exist('Evts2','var') && ~isempty(Evts2.codes)
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
            %repeat for Evts3 separately
            if exist('Evts3','var') && ~isempty(Evts3.codes)
                while 1
                    if Evts3.codes(1,1)==0
                        Evts3.codes(1,:)=[];
                        Evts3.times(1)=[];
                        Evts3.length=Evts.length-1;
                        Evts3.text(1,:)=[];
                    else
                        break
                    end
                end
            end
            
            for i=1:size(Evts.codes,1)
                Evts.source(i,1) = 1;
            end
            if exist('Evts2','var')
                for i=1:size(Evts2.codes,1)
                    Evts2.source(i,1) = 2;
                end
            end
            if exist('Evts3','var')
                for i=1:size(Evts3.codes,1)
                    Evts3.source(i,1) = 3;
                end
            end
            
            
            if exist('Evts2','var')&&exist('Evts','var')&&isfield(Evts2,'source')
                %combine events lists
                Evts.times=[Evts.times;Evts2.times];
                Evts.codes=[Evts.codes;Evts2.codes];
                Evts.text=[Evts.text;Evts2.text];
                Evts.length=Evts.length+Evts2.length;
                Evts.source=[Evts.source;Evts2.source];
            end
            
            if exist('Evts3','var')&&exist('Evts','var')&&isfield(Evts3,'source')
                %combine events lists
                Evts.times=[Evts.times;Evts3.times];
                Evts.codes=[Evts.codes;Evts3.codes];
                Evts.text=[Evts.text;Evts3.text];
                Evts.length=Evts.length+Evts3.length;
                Evts.source=[Evts.source;Evts3.source];
            end
            
            
            %delete second row of zeros if there are two concurrent rows of zeros
            for i=size(Evts.codes,1):-1:2
                if Evts.codes(i)==0&&Evts.codes(i-1)==0
                    Evts.codes(i,:)=[];
                    Evts.times(i)=[];
                    Evts.source(i)=[];
                    Evts.text(i,:)=[];
                    Evts.length=Evts.length-1;
                    disp('delete_row_of_zeros');
                end
            end
            
            
            %add a row of zeros if there is a missing row of zeros. Assume end of event is the start of the next event (or end of file/Inf).
            for i=size(Evts.codes,1):-1:1
                if i==length(Evts.codes)
                    if Evts.codes(i)>0
                        Evts.codes=[Evts.codes;[0 0 0 0]];
                        Evts.times=[Evts.times;Inf];
                        Evts.text=[Evts.text;[blanktext]];
                        Evts.length=Evts.length+1;
                        Evts.source(i)=Evts.source(i-1);
                        disp('added_row_of_zeros to the end');
                    end
                end
                if i>1
                    if Evts.codes(i)>0&&Evts.codes(i-1)>0
                        Evts.codes=[Evts.codes(1:(i-1),:);[0 0 0 0];Evts.codes(i:end,:)];
                        Evts.times=[Evts.times(1:i);Evts.times(i:end)];
                        Evts.text=[Evts.text(1:(i-1),:);[blanktext];Evts.text(i:end,:)];
                        Evts.length=Evts.length+1;
                        Evts.source=[Evts.source(1:i);Evts.source(i:end)];
                        disp('added_row_of_zeros');
                    end
                end
            end
            
            %remove zero rows:
            Evts.codes(2:2:end,:)=[];
            Evts.codes(:,2:4)=[];
            Evts.endtimes=Evts.times(2:2:end);
            Evts.starttimes=Evts.times(1:2:end);
            Evts.durations=Evts.endtimes-Evts.starttimes;
            Evts.text(2:2:end,:)=[];
            % Evts.text(:,7:end)=[];
            %             Evts.text = cellstr(Evts.text);
            Evts.source=Evts.source(1:2:end);
            
            % Overlapping events do not apply for SpikeDISE since Evts,
            % Evts2, and Evts3 are all recording different things
            % check for overlapping events and remove
            %Intention is to prioritize the events list that is shorter i.e. represents corrections
            %This corrections Evts list should be called Evts3, and if not there will be a warning given.
            overlappingeventsoption=0;
            if overlappingeventsoption
                disp('overlappingeventsoption');
                evtssources = [sum(Evts.source==1) sum(Evts.source==2) sum(Evts.source==3)];
                evtsrespsources = [sum(Evts.source==1&Evts.codes~=1) sum(Evts.source==2&Evts.codes~=1) sum(Evts.source==3&Evts.codes~=1)];
                evtsarsources = [sum(Evts.source==1&Evts.codes==1) sum(Evts.source==2&Evts.codes==1) sum(Evts.source==3&Evts.codes==1)];
                
                evtsrespsources2 = evtsrespsources;
                evtsrespsources2(evtsrespsources2==0)=NaN;
                evtsrespsources2(evtsarsources>0)=NaN;
                
                if sum(~isnan(evtsrespsources2))>1
                    if 0
                        [~,prioritizesource] = min(evtsrespsources2);
                    else
                        prioritizesource = find(~isnan(evtsrespsources2),1,'last');
                    end
                    % if prioritizesource~=3
                    %     disp('warning, not prioritizing Evts3');
                    % end
                    
                    lengthevtscodes=length(Evts.codes);
                    for i=length(Evts.codes):-1:1
                        for j=length(Evts.codes):-1:1
                            if i>length(Evts.codes)||j>length(Evts.codes)||i==j||Evts.codes(i)==1||Evts.codes(j)==1
                                continue
                            end
                            a=Evts.starttimes(i);
                            b=Evts.endtimes(i);
                            A=Evts.starttimes(j);
                            B=Evts.endtimes(j);
                            if ~(A>a&&A<b || B>a&&B<b || a>A&&a<B || b>A&&b<B)
                                %no overlap
                                continue
                            end
                            if Evts.source(j)~=prioritizesource
                                disp(['removing event: detected overlap i=' num2str(i) '[' num2str(Evts.source(i)), '], j=' num2str(j) '[' num2str(Evts.source(j)) '], removing j ' num2str(Evts.codes(j)) ', keeping i ' num2str(Evts.codes(i))]);
                                Evts.source(j)=[];
                                Evts.codes(j)=[];
                                Evts.endtimes(j)=[];
                                Evts.starttimes(j)=[];
                                Evts.durations(j)=[];
                                Evts.text(j,:)=[];
                            else
                                disp(['removing event: i=' num2str(i) '[' num2str(Evts.source(i)), '], j=' num2str(j) '[' num2str(Evts.source(j)) '], removing i: ' num2str(Evts.codes(i)) ', keeping j ' num2str(Evts.codes(j))]);
                                Evts.source(i)=[];
                                Evts.codes(i)=[];
                                Evts.endtimes(i)=[];
                                Evts.starttimes(i)=[];
                                Evts.durations(i)=[];
                                Evts.text(i,:)=[];
                            end
                        end
                    end
                end
            end
            %             EventStart=Evts.starttimes + StartTime;
            EventStart=Evts.starttimes; % Made by DV on 6/30/21
            EventDuration = Evts.durations;
            
            %% Codify VOTE classification
            % I decided to do this in post processing
            
            
            %%
            
            %             EventLabelOptions= {...
            %                 'AR' ,...
            %                 'ApO',...
            %                 'ApC',...
            %                 'HypO',...
            %                 'M',...
            %                 'HypC',...
            %                 'HypO2%',...
            %                 'HypC2%',...
            %                 'HypOx',...
            %                 'HypCx',...
            %                 'Ap-C',... %older labels
            %                 'H',... %older labels
            %                 'Hyp-C',... %older labels
            %                 'ap-O',... %older labels
            %                 };
            
            %             EventCodeOptions=[1 2 3 4 5 6 7 8 9 10 3 4 6 2];
            %             EventExact1orStartsWith2 = [1 1 1 1 1 1 1 1 1 1 1 1 1 1];
            %             for i=1:length(EventStart)
            %                 EventName{i} = deblank(Evts.text(i,:));
            %             end
            EventName = cellstr(Evts.text);
            
            Hypnogram = NaN; % no sleep stages in DISE
            %             Hypnogram=double(Epochs.codes(:,1));
            %
            %             Hypnogram(Hypnogram==0)=-1; %W=-1,N1=1,N2=2,N3=3,R=5
            %             Hypnogram(Hypnogram==5)=0;  %W=-1,R=0,N1=1,N2=2,N3=3
            %             Hypnogram=3-Hypnogram; %.MAT: 0=N3,1=N2,2=N1,3=R,4=W
            %             Hypnogram(Hypnogram==-3)=NaN;
        end
        
    case {'NoxMat'}
        %%
        echannelnameoptions.AllEvents={'AllEvents'};
        echannelnameoptions.HypnogramInfo={'HypnogramInfo'};
        
        channelnamestemp=fieldnames(echannelnameoptions);
        
        filehandle = matfile([Filenames{5} Filenames{2}]);
        w = whos(filehandle);
        
        for i=1:length(channelnamestemp)
            temp=eval(['echannelnameoptions.' char(channelnamestemp(i))]);
            foundamatch=0;
            for nn=1:length(temp)
                %Does it exist?
                for j=1:length(w)
                    if strcmp(w(j).name,char(temp(nn)))
                        eval([channelnamestemp{i} '=filehandle.' char(temp(nn)),';']);
                        foundamatch=1;
                        break
                    end
                end
                if foundamatch
                    break
                end
            end
        end
        
        if ~exist('AllEvents')
            disp('Error: Nox Source has no events file [AllEvents]');
        end
        if ~exist('HypnogramInfo')
            disp('Error: Nox Source has no events file [HypnogramInfo]');
        end
        
        %fix tags
        clear tempnewlabelsNoSpace
        for i=1:size(AllEvents.Tags,1)
            temp = AllEvents.Tags(i,:);
            itemp = find(temp==' ');
            temp(itemp)=[];
            tempnewlabelsNoSpace{i,1}=temp;
        end
        AllEvents.Tags = tempnewlabelsNoSpace;
        
        EventName = AllEvents.Tags;
        EventStart = AllEvents.StartEndDuration(1,:); %warning all data are nearest second
        EventStart=EventStart(:);
        %             EventStart=EventStart*86400;
        %             EventStart(EventStart<43200)=EventStart(EventStart<43200)+86400;
        EventDuration = AllEvents.StartEndDuration(3,:);
        EventDuration=EventDuration(:);
        EventLabelOptions= {...
            'arousal' ,...
            'arousal-spontaneous', ...
            'arousal-limbmovement', ...
            'arousal-plm', ...
            'arousal-respiratory', ...
            'apnea-obstructive',...
            'apnea-central',...
            'hypopnea',...
            'hypopnea-obstructive',...
            'hypopnea-central',...
            'apnea-mixed',...
            'rera',...
            };
        EventCodeOptions=[1 1 1 1 1 2 3 4 4 6 5 12];
        EventExact1orStartsWith2 = [2 1 1 1 1 1 1 1 1 1 1 1];
        
        clear tempnewlabelsNoSpace
        for i=1:size(HypnogramInfo.Tags,1)
            temp = HypnogramInfo.Tags(i,:);
            itemp = find(temp==' ');
            temp(itemp)=[];
            tempnewlabelsNoSpace{i,1}=temp;
        end
        HypnogramInfo.Tags = tempnewlabelsNoSpace;
        
        %
        %         %Arousals.StartEndDuration = Arousals.StartEndDuration';
        %         Events.StartEndDuration = Events.StartEndDuration';
        HypnogramInfo.StartEndDuration = HypnogramInfo.StartEndDuration';
        
        
        
        clear tempnewlabelsNoSpace
        
        % EventStart=[Arousals.StartEndDuration(:,1) ; Events.StartEndDuration(:,1)];
        % EventDuration=[Arousals.StartEndDuration(:,3) ; Events.StartEndDuration(:,3)];
        
        
        % EventName = [Arousals.Tags ; Events.Tags];
        
        
        iEpochDuration=(round(HypnogramInfo.StartEndDuration(:,3)/30));
        iEpochStart=(1+round((HypnogramInfo.StartEndDuration(:,1)-StartTime)/30)); %HypStart(Ilist(1))
        iEpochEnd=iEpochStart+iEpochDuration-1;
        EpochCodes = HypnogramInfo.StartEndDuration(:,5);
        %checktemp=[iEpochStart iEpochDuration iEpochEnd EpochCodes];
        starttimediff = HypnogramInfo.StartEndDuration(1,1)-StartTime - 30*(iEpochStart(1)-1); %first epochs filled with NaNs; starttimediff is just residual after rounding
        %starttimediffN=round((starttimediff/30))
        displaytext=['Start Time Hypnogram minus Signals: ' num2str(starttimediff) ' s'];
        if ~isempty(handletext), disp(displaytext); set(handletext,'String',displaytext); drawnow; end
        
        Nepochs=ceil((EndTime-StartTime)/30); %make this based on Time(end)-StartTime...sum(iEpochDuration)
        Hypnogram=NaN*zeros(Nepochs,1);
        for i=1:length(iEpochDuration)
            %                 if EpochCodes(i)==2
            %                     continue
            %                 end
            try
                Hypnogram(iEpochStart(i):iEpochEnd(i))=EpochCodes(i);
            catch me
                'warning: hypnogram starts before signals';
            end
        end
        
     case {'NoxMatSAS'} %code is untested
        %%
        echannelnameoptions.AllEvents={'AllEvents'};
        echannelnameoptions.HypnogramInfo={'HypnogramInfo'};
        
        channelnamestemp=fieldnames(echannelnameoptions);
        
        filehandle = matfile([Filenames{5} Filenames{2}]);
        w = whos(filehandle);
        
        for i=1:length(channelnamestemp)
            temp=eval(['echannelnameoptions.' char(channelnamestemp(i))]);
            foundamatch=0;
            for nn=1:length(temp)
                %Does it exist?
                for j=1:length(w)
                    if strcmp(w(j).name,char(temp(nn)))
                        eval([channelnamestemp{i} '=filehandle.' char(temp(nn)),';']);
                        foundamatch=1;
                        break
                    end
                end
                if foundamatch
                    break
                end
            end
        end
        
        if ~exist('AllEvents')
            disp('Error: Nox Source has no events file [AllEvents]');
        end
        if ~exist('HypnogramInfo')
            disp('Error: Nox Source has no events file [HypnogramInfo]');
        end
        
        %fix tags
        clear tempnewlabelsNoSpace
        for i=1:size(AllEvents.Tags,1)
            temp = AllEvents.Tags(i,:);
            itemp = find(temp==' ');
            temp(itemp)=[];
            tempnewlabelsNoSpace{i,1}=temp;
        end
        AllEvents.Tags = tempnewlabelsNoSpace;
        
        EventName = AllEvents.Tags;
        EventStart = AllEvents.StartEndDuration(1,:); %warning all data are nearest second
        EventStart=EventStart(:);
        %             EventStart=EventStart*86400;
        %             EventStart(EventStart<43200)=EventStart(EventStart<43200)+86400;
        EventDuration = AllEvents.StartEndDuration(3,:);
        EventDuration=EventDuration(:);
        EventLabelOptions= {...
            'arousal' ,...
            'arousal-spontaneous', ...
            'arousal-limbmovement', ...
            'arousal-plm', ...
            'arousal-respiratory', ...
            'apnea-obstructive',...
            'apnea-central',...
            'hypopnea',...
            'hypopnea-obstructive',...
            'hypopnea-central',...
            'apnea-mixed',...
            'rera',...
            };
        EventCodeOptions=[1 1 1 1 1 2 3 4 4 6 5 12];
        EventExact1orStartsWith2 = [2 1 1 1 1 1 1 1 1 1 1 1];
        
        clear tempnewlabelsNoSpace
        for i=1:size(HypnogramInfo.Tags,1)
            temp = HypnogramInfo.Tags(i,:);
            itemp = find(temp==' ');
            temp(itemp)=[];
            tempnewlabelsNoSpace{i,1}=temp;
        end
        HypnogramInfo.Tags = tempnewlabelsNoSpace;
        
        %
        %         %Arousals.StartEndDuration = Arousals.StartEndDuration';
        %         Events.StartEndDuration = Events.StartEndDuration';
        HypnogramInfo.StartEndDuration = HypnogramInfo.StartEndDuration';
        
        
        
        clear tempnewlabelsNoSpace
        
        % EventStart=[Arousals.StartEndDuration(:,1) ; Events.StartEndDuration(:,1)];
        % EventDuration=[Arousals.StartEndDuration(:,3) ; Events.StartEndDuration(:,3)];
        
        
        % EventName = [Arousals.Tags ; Events.Tags];
        
        
        iEpochDuration=(round(HypnogramInfo.StartEndDuration(:,3)/30));
        iEpochStart=(1+round((HypnogramInfo.StartEndDuration(:,1)-StartTime)/30)); %HypStart(Ilist(1))
        iEpochEnd=iEpochStart+iEpochDuration-1;
        EpochCodes = HypnogramInfo.StartEndDuration(:,5);
        %checktemp=[iEpochStart iEpochDuration iEpochEnd EpochCodes];
        starttimediff = HypnogramInfo.StartEndDuration(1,1)-StartTime - 30*(iEpochStart(1)-1); %first epochs filled with NaNs; starttimediff is just residual after rounding
        %starttimediffN=round((starttimediff/30))
        displaytext=['Start Time Hypnogram minus Signals: ' num2str(starttimediff) ' s'];
        if ~isempty(handletext), disp(displaytext); set(handletext,'String',displaytext); drawnow; end
        
        Nepochs=ceil((EndTime-StartTime)/30); %make this based on Time(end)-StartTime...sum(iEpochDuration)
        Hypnogram=NaN*zeros(Nepochs,1);
        for i=1:length(iEpochDuration)
            %                 if EpochCodes(i)==2
            %                     continue
            %                 end
            try
                Hypnogram(iEpochStart(i):iEpochEnd(i))=EpochCodes(i);
            catch me
                'warning: hypnogram starts before signals';
            end
        end
        
     case {'NoxMatT3'} %code is untested
        %%
        echannelnameoptions.AllEvents={'AllEvents'};
        echannelnameoptions.HypnogramInfo={'HypnogramInfo'};
        
        channelnamestemp=fieldnames(echannelnameoptions);
        
        filehandle = matfile([Filenames{5} Filenames{2}]);
        w = whos(filehandle);
        
        for i=1:length(channelnamestemp)
            temp=eval(['echannelnameoptions.' char(channelnamestemp(i))]);
            foundamatch=0;
            for nn=1:length(temp)
                %Does it exist?
                for j=1:length(w)
                    if strcmp(w(j).name,char(temp(nn)))
                        eval([channelnamestemp{i} '=filehandle.' char(temp(nn)),';']);
                        foundamatch=1;
                        break
                    end
                end
                if foundamatch
                    break
                end
            end
        end
        
        if ~exist('AllEvents')
            disp('Error: Nox Source has no events file [AllEvents]');
        end
        if ~exist('HypnogramInfo')
            disp('Error: Nox Source has no events file [HypnogramInfo]');
        end
        
        %fix tags
        clear tempnewlabelsNoSpace
        for i=1:size(AllEvents.Tags,1)
            temp = AllEvents.Tags(i,:);
            itemp = find(temp==' ');
            temp(itemp)=[];
            tempnewlabelsNoSpace{i,1}=temp;
        end
        AllEvents.Tags = tempnewlabelsNoSpace;
        
        EventName = AllEvents.Tags;
        EventStart = AllEvents.StartEndDuration(1,:); %warning all data are nearest second
        EventStart=EventStart(:);
        %             EventStart=EventStart*86400;
        %             EventStart(EventStart<43200)=EventStart(EventStart<43200)+86400;
        EventDuration = AllEvents.StartEndDuration(3,:);
        EventDuration=EventDuration(:);
        EventLabelOptions= {...
            'arousal' ,...
            'arousal-spontaneous', ...
            'arousal-limbmovement', ...
            'arousal-plm', ...
            'arousal-respiratory', ...
            'apnea-obstructive',...
            'apnea-central',...
            'hypopnea',...
            'hypopnea-obstructive',...
            'hypopnea-central',...
            'apnea-mixed',...
            'rera',...
            };
        EventCodeOptions=[1 1 1 1 1 2 3 4 4 6 5 12];
        EventExact1orStartsWith2 = [2 1 1 1 1 1 1 1 1 1 1 1];
        
        clear tempnewlabelsNoSpace
        for i=1:size(HypnogramInfo.Tags,1)
            temp = HypnogramInfo.Tags(i,:);
            itemp = find(temp==' ');
            temp(itemp)=[];
            tempnewlabelsNoSpace{i,1}=temp;
        end
        HypnogramInfo.Tags = tempnewlabelsNoSpace;
        
        %
        %         %Arousals.StartEndDuration = Arousals.StartEndDuration';
        %         Events.StartEndDuration = Events.StartEndDuration';
        HypnogramInfo.StartEndDuration = HypnogramInfo.StartEndDuration';
        
        
        
        clear tempnewlabelsNoSpace
        
        % EventStart=[Arousals.StartEndDuration(:,1) ; Events.StartEndDuration(:,1)];
        % EventDuration=[Arousals.StartEndDuration(:,3) ; Events.StartEndDuration(:,3)];
        
        
        % EventName = [Arousals.Tags ; Events.Tags];
        
        
        iEpochDuration=(round(HypnogramInfo.StartEndDuration(:,3)/30));
        iEpochStart=(1+round((HypnogramInfo.StartEndDuration(:,1)-StartTime)/30)); %HypStart(Ilist(1))
        iEpochEnd=iEpochStart+iEpochDuration-1;
        EpochCodes = HypnogramInfo.StartEndDuration(:,5);
        %checktemp=[iEpochStart iEpochDuration iEpochEnd EpochCodes];
        starttimediff = HypnogramInfo.StartEndDuration(1,1)-StartTime - 30*(iEpochStart(1)-1); %first epochs filled with NaNs; starttimediff is just residual after rounding
        %starttimediffN=round((starttimediff/30))
        displaytext=['Start Time Hypnogram minus Signals: ' num2str(starttimediff) ' s'];
        if ~isempty(handletext), disp(displaytext); set(handletext,'String',displaytext); drawnow; end
        
        Nepochs=ceil((EndTime-StartTime)/30); %make this based on Time(end)-StartTime...sum(iEpochDuration)
        Hypnogram=NaN*zeros(Nepochs,1);
        for i=1:length(iEpochDuration)
            %                 if EpochCodes(i)==2
            %                     continue
            %                 end
            try
                Hypnogram(iEpochStart(i):iEpochEnd(i))=EpochCodes(i);
            catch me
                'warning: hypnogram starts before signals';
            end
        end    
        
        
    case  'Nox'
     % there might be unscored epochs- gaps could exist b/w 2 scored
     % epochs..need to be cautious..there wont be labels
     % there could be partial epochs that started b4r record start and
     % could be deleted in the events list
        %%
        opts = detectImportOptions([char(Filenames(4)) filesep char(Filenames(2))],'NumHeaderLines',2);
        opts.VariableNamesRange = 'A1';
        opts.DataRange = 'A3';
        NoxEvtT = readtable([char(Filenames(4)) filesep char(Filenames(2))],opts,'ReadVariableNames',true);
        
        EventName =  NoxEvtT.Event;
%         EvtNameUnique=unique(EventName);

        EventStart = mod(datenum(NoxEvtT.StartTime),1)*86400;
        EventStart(EventStart<43200) = EventStart(EventStart<43200) + 86400; %
        
        EventDuration = NoxEvtT.Duration;
        NoxEvtT.EventStart=EventStart;
        
         EventLabelOptions= {...
            'Arousal' ,...
            'LM Arousal', ...
            'Respiratory Arousal', ...
            'A. Obstructive',...
            'A. Central',...
            'Hypopnea',...
            'H. Obstructive',...
            'Rera',...
            };
        EventCodeOptions=[1 1 1 2 3 4 4 12];
        EventExact1orStartsWith2 = [1 1 1 1 1 1 1 1];
        
        
        StageLabelOptions = {'Wake','REM','N1','N2','N3'}';
        StageCodeOptions = [4 3 2 1 -1] ; %can not be zeros here
        
         clear match1 match2
        EventName = NoxEvtT.Event;
        for m=1:length(StageLabelOptions)% each column is an event type
            match1(:,m)=strcmpi(string(EventName'),string(StageLabelOptions(m)));
            match2(:,m)=startsWith(string(EventName'),string(StageLabelOptions(m)));
        end
        match = 1*( match1>0 | match2>0);
        temp = match.*StageCodeOptions;
        temp(temp==0)=NaN;
        temp(temp==-1)=0;
        NoxEvtT.EpochCodes=min(temp')'; %just in case two criteria are detected, will leave the lowest code value in list
        
        HypnogramInfo=NoxEvtT(~isnan(NoxEvtT.EpochCodes),:);
        
        Nepochs=ceil((EndTime-StartTime)/30);
        
        
         HypnogramInfo.iEpochDuration=(round(HypnogramInfo.Duration/30)); %number of epochs in current row
        HypnogramInfo.iEpochStart=(1+round((HypnogramInfo.EventStart-StartTime)/30)); %HypStart(Ilist(1))
        HypnogramInfo.iEpochEnd=HypnogramInfo.iEpochStart+HypnogramInfo.iEpochDuration-1;
  
        starttimediff = HypnogramInfo.EventStart(1,1)-StartTime - 30*(HypnogramInfo.iEpochStart(1)-1); %first epochs filled with NaNs; starttimediff is just residual after rounding
        endtimediff = -(HypnogramInfo.EventStart(end,1)-EndTime)./30 ;
        %starttimediffN=round((starttimediff/30))
        displaytext=['Start Time Hypnogram minus Signals: ' num2str(starttimediff) ' s'];
        if ~isempty(handletext), disp(displaytext); set(handletext,'String',displaytext); drawnow; end
        
      
        Hypnogram=NaN*zeros(Nepochs,1);
        for i=1:length(HypnogramInfo.iEpochDuration)
            %                 if EpochCodes(i)==2
            %                     continue 
            %                 end
            try
                Hypnogram(HypnogramInfo.iEpochStart(i):HypnogramInfo.iEpochEnd(i))=HypnogramInfo.EpochCodes(i);
            catch me
                'warning: hypnogram starts before signals';
            end
        end   
   
    case {'RemLogic','RemLogicXML'}  %code is untested
        %%
        
        directory = char(Filenames(5));
        fname = char(Filenames(2));
        
        %Xml events
        displaytext=['Get XML Events data: ' fname];
        if ~isempty(handletext), disp(displaytext); set(handletext,'String',displaytext); drawnow; end
        
        %importing
        S=xml2struct([directory fname]);
        
        HeartBeatVersion=0;
        if isfield(S,'PSGAnnotation')
            S.EventExport = S.PSGAnnotation;
            S = rmfield(S,'PSGAnnotation');
            if isfield(S.EventExport,'ScoredEvents')
                S.EventExport.Events = S.EventExport.ScoredEvents;
                S.EventExport = rmfield(S.EventExport,'ScoredEvents');
            end
            if isfield(S.EventExport.Events,'ScoredEvent')
                S.EventExport.Events.Event = S.EventExport.Events.ScoredEvent;
                S.EventExport.Events = rmfield(S.EventExport.Events,'ScoredEvent');
            end
            HeartBeatVersion=1;
        end
        
        %S.EventExport.Events.Event{1,i}.Type.Text
        
        if ~isfield(S,'EventExport') || ~isfield(S.EventExport.Events,'Event') %not checked
            displaytext=['Warning: No scored events'];
            if ~isempty(handletext), disp(displaytext); set(handletext,'String',displaytext); drawnow; end
        end
        
        Nevents = length(S.EventExport.Events.Event);
        
        EventName={[]};
        EventStart=[];
        EventDuration=[];
        
        clear EventName EventStart EventStop EventDuration
        
        if HeartBeatVersion==0
            formatIn = 'HH:MM:SS.FFF';
            for i=1:Nevents

                    EventName{i,1} = S.EventExport.Events.Event{1,i}.Type.Text;

                    temp = S.EventExport.Events.Event{1,i}.StartTime.Text(12:end);

                temp = mod(datenum(temp,formatIn),1)*86400;
                
                test = [floor(temp/86400*24) floor(mod(temp/86400*24,1)*60) mod(mod(temp/86400*24,1)*60,1)*60];
                if temp<86400/2,temp=temp+86400; end
                EventStart(i,1) = temp;
                %x  y x2
                temp = S.EventExport.Events.Event{1,i}.StopTime.Text(12:end);
                temp = mod(datenum(temp,formatIn),1)*86400;
                if temp<86400/2,temp=temp+86400; end
                EventStop(i,1) = temp;
            end
        else
            %formatIn = 'HH:MM:SS.FFF';
            for i=1:Nevents

                 EventName{i,1} = S.EventExport.Events.Event{1,i}.EventConcept.Text;
  
                 temp = S.EventExport.Events.Event{1,i}.Start.Text; %is time since start recording
               
                temp = StartTime + str2num(temp);
                
                %test = [floor(temp/86400*24) floor(mod(temp/86400*24,1)*60) mod(mod(temp/86400*24,1)*60,1)*60];
                if temp<86400/2,temp=temp+86400; end
                EventStart(i,1) = temp;
                %x  y x2
                
                temp = str2num(S.EventExport.Events.Event{1,i}.Duration.Text);                
                EventStop(i,1) = EventStart(i,1) + temp;
            end
        end
        EventDuration = EventStop - EventStart;
                    
        SleepEventFlag = startsWith(string(EventName),'SLEEP-'); %or should be "SLEEP-" ? as per other version
        %Note HeartBeatVersion has no sleep, so the sleep table is going to
        %be empty
        
        EventTable = table(EventName,EventStart,EventStop,EventDuration,SleepEventFlag);
        SleepTable = EventTable(SleepEventFlag==1,:);
        
        %W=0,N1=1,N2=2,N3=3,R=5
        EventCodesListUnique = unique(EventName);
        EventLabelOptions= {...
            'AROUSAL',...
            'APNEA-CENTRAL',...
            'APNEA-MIXED',...
            'APNEA-OBSTRUCTIVE',...
            'HYPOPNEA',...
            'HYPOPNEA-OBSTRUCTIVE',...
            'HYPOPNEA-CENTRAL',...
            'Obstructive apnea|APNEA-OBSTRUCTIVE',... 
            'Central apnea|APNEA-CENTRAL',...
            'Hypopnea|HYPOPNEA',...
            'Mixed apnea|APNEA-MIXED'...
            };
        EventCodeOptions=[1 3 5 2 4 4 6 2 3 4 5];
        EventExact1orStartsWith2 = [2 1 1 1 1 1 1 1 1 1 1];
        
        displaytext=['Get Hypnogram data: ' fname];
        if ~isempty(handletext), disp(displaytext); set(handletext,'String',displaytext); drawnow; end
        
        if HeartBeatVersion==0
        SleepCategoriesLabels= {...
            'SLEEP-S0',...
            'SLEEP-REM',...
            'SLEEP-S1',...
            'SLEEP-S2',...
            'SLEEP-S3',...
            };
        SleepCodes = [4 3 2 1 0]; %Terrill code (%.MAT: 0=N3,1=N2,2=N1,3=R,4=W)
        
        %get codes
        temp = string(SleepTable{:,1})==SleepCategoriesLabels;
        temp = temp.*SleepCodes;
        temp = max(temp')';
        SleepTable.codes = temp;
        clear temp;
        
        I = SleepTable.EventStart<StartTime;
        SleepTable(I,:)=[];
        
        NepochChanges = size(SleepTable,1);
        iEpochDuration=(round(SleepTable.EventDuration/30));
        iEpochStart=(1+round((SleepTable.EventStart-StartTime)/30)); %HypStart(Ilist(1))
        iEpochEnd=iEpochStart+iEpochDuration-1;
        checktemp=[iEpochStart iEpochDuration iEpochEnd SleepTable.codes];
        
        %first epochs filled with NaNs; starttimediff is just residual after rounding
        
        starttimediff = SleepTable.EventStart(1)-StartTime - 30*(iEpochStart(1)-1); %first epochs filled with NaNs; starttimediff is just residual after rounding
        
        %starttimediffN=round((starttimediff/30))
        displaytext=['Start Time Hypnogram minus Signals: ' num2str(starttimediff) ' s'];
        if ~isempty(handletext), disp(displaytext); set(handletext,'String',displaytext); drawnow; end
        
        Nepochs=ceil((EndTime-StartTime)/30); %make this based on Time(end)-StartTime...sum(iEpochDuration)
        Hypnogram=NaN*zeros(Nepochs,1);
        for i=1:length(iEpochDuration)
            Hypnogram(iEpochStart(i):iEpochEnd(i))=SleepTable.codes(i);
        end
        else
            
            Nepochs=ceil((EndTime-StartTime)/30); %make this based on Time(end)-StartTime...sum(iEpochDuration)
        Hypnogram=4*ones(Nepochs,1);
        ii = find(EventTable.EventName=="Beginning of analysis period|ANALYSIS-START")
        iii = find(EventTable.EventName=="End of analysis period|ANALYSIS-STOP")
        StartAnalysis = EventTable.EventStart(ii);
        EndAnalysis = EventTable.EventStart(iii);
        StartAnalysisEpoch = round((StartAnalysis - StartTime)/30 + 1)
        EndAnalysisEpoch = round((EndAnalysis - StartTime)/30 + 1)
        Hypnogram(StartAnalysisEpoch:EndAnalysisEpoch)=2; %Stage 1 sleep for all analysis period
            
        end
        %end REMlogicXML
        
    case 'RemLogicTextFr' %code is untested
        %%
        directory = char(Filenames(5));
        fname = char(Filenames(2));
        
        clear eventstemp
        fid = fopen([directory '\' fname]);
        i=1;
        while 1
            eventstemp{i,:} = fgetl(fid);  %read in the data
            if eventstemp{i,:}==-1
                eventstemp(i,:) = [];
                break
            end
            i=i+1;
        end
        fclose(fid);   %close the file
        
        HeaderTextList = {'Evénement','Durée[s]'};
        I=[];
        for ii=1:length(HeaderTextList)
            I=[I,contains(eventstemp,HeaderTextList{ii})];
        end
        I=sum(I,2)/length(HeaderTextList);
        NPreHeaderlines=find(I==1)-1;
        HeaderTextExact = eventstemp{I==1};
        
        % Import Table (now that we know where the data starts)
        EventTable_temp = readtable([directory fname],'Delimiter','\t','ReadVariableNames',1,'HeaderLines',NPreHeaderlines);
        
        %             clear timetemp
        %             for i=1:size(EventTable_temp,1)
        %                 eventtemp=EventTable_temp.Heure_hh_mm_ss_(i};
        %                 %Time is in sec since the previous day's midnight.
        %                 tempstr=eventtemp;%eventtemp(I(1)+1:I(2)-1);
        %                 I2=find(tempstr==':');
        %                 timetemp(i,1)=86400*(tempstr(I2(2)+4)=='A')+43200*(tempstr(I2(2)+4)=='P')+3600*mod(str2num(tempstr(1:I2(1)-1)),12)+60*str2num(tempstr(I2(1)+1:I2(1)+2)) + 1*str2num(tempstr(I2(2)+1:I2(2)+2)) + str2num(tempstr(end-3:end));
        %             end
        %             timetemp(timetemp<43200)=timetemp(timetemp<43200)+86400;
        %
        %             clear timechangei
        %             timechangei=find([0; diff(timetemp)]<0);
        %
        %             while ~isempty(timechangei)
        %                 timetemp(timechangei(1):end)=timetemp(timechangei(1):end)+86400;
        %                 timechangei=find([0; diff(timetemp)]<0);
        %             end
        
        % This code is replacing commented section above
        StartTimeTemp = seconds(EventTable_temp.Heure_hh_mm_ss_); % seconds from midnight
        [~,MorningStartIdx] = min(StartTimeTemp); % min will be time at or immediately after midnight
        StartT = [StartTimeTemp(1:MorningStartIdx-1); StartTimeTemp(MorningStartIdx:end) + 86400]; %add a full day to morning time difference
        
        EventDuration = EventTable_temp.Dur_e_s_;
        EventStart = StartT;
        EventStop=EventStart+EventDuration;
        EventName=EventTable_temp.Ev_nement;
        
        EventTable = table(EventName,EventStart,EventStop,EventDuration);
        SleepTable = EventTable;
        
        %W=0,N1=1,N2=2,N3=3,R=5
        EventCodesListUnique = unique(EventName);
        EventLabelOptions= {...
            'AROUSAL',...
            'APNEA-CENTRAL',...
            'APNEA-MIXED',...
            'APNEA-OBSTRUCTIVE',...
            'HYPOPNEA',...
            };
        EventCodeOptions=[1 3 5 2 4];
        EventExact1orStartsWith2 = [2 1 1 1 2];
        
        displaytext=['Get Hypnogram data: ' fname];
        if ~isempty(handletext), disp(displaytext); set(handletext,'String',displaytext); drawnow; end
        
        SleepCategoriesLabels= {...
            'SLEEP-S0',...
            'SLEEP-REM',...
            'SLEEP-S1',...
            'SLEEP-S2',...
            'SLEEP-S3',...
            };
        SleepCodes = [4 3 2 1 0]; %Terrill code (%.MAT: 0=N3,1=N2,2=N1,3=R,4=W)
        
        %get codes
        I = string(SleepTable.EventName)==SleepCategoriesLabels;
        I2 = sum(I,2)>0;
        temp = I.*SleepCodes;
        temp = max(temp')';
        SleepTable.codes = temp;
        
        SleepTable = SleepTable(I2,:);
        temp = SleepTable.EventStart(2:end) - SleepTable.EventStart(1:end-1);
        templast = EndTime-SleepTable.EventStart(end);
        SleepTable.EventDuration(1:end-1) =temp;
        SleepTable.EventDuration(end)=templast;
        SleepTable.EventStop=SleepTable.EventStart+SleepTable.EventDuration;
        clear temp;
        
        NepochChanges = size(SleepTable,1);
        iEpochDuration=(round(SleepTable.EventDuration/30));
        iEpochStart=(1+round((SleepTable.EventStart-StartTime)/30)); %HypStart(Ilist(1))
        iEpochEnd=iEpochStart+iEpochDuration-1;
        checktemp=[iEpochStart iEpochDuration iEpochEnd SleepTable.codes];
        
        %first epochs filled with NaNs; starttimediff is just residual after rounding
        starttimediff = SleepTable.EventStart(1)-StartTime - 30*(iEpochStart(1)-1); %first epochs filled with NaNs; starttimediff is just residual after rounding
        
        %starttimediffN=round((starttimediff/30))
        displaytext=['Start Time Hypnogram minus Signals: ' num2str(starttimediff) ' s'];
        if ~isempty(handletext), disp(displaytext); set(handletext,'String',displaytext); drawnow; end
        
        Nepochs=ceil((EndTime-StartTime)/30); %make this based on EndTime-StartTime...sum(iEpochDuration)
        Hypnogram=NaN*zeros(Nepochs,1);
        for i=1:length(iEpochDuration)
            Hypnogram(iEpochStart(i):iEpochEnd(i))=SleepTable.codes(i);
        end
        
    case 'RemLogicText' %code is under development (i.e. incomplete)
        %%
        directory = char(Filenames(5));
        fname = char(Filenames(2));
        
        clear eventstemp
        fid = fopen([directory '\' fname]);
        i=1;
        while 1
            eventstemp{i,:} = fgetl(fid);  %read in the data
            if eventstemp{i,:}==-1
                eventstemp(i,:) = [];
                break
            end
            i=i+1;
        end
        fclose(fid);   %close the file
        
        HeaderTextList = {'Event','Duration[s]'};
        I=[];
        for ii=1:length(HeaderTextList)
            I=[I,contains(eventstemp,HeaderTextList{ii})];
        end
        I=sum(I,2)/length(HeaderTextList);
        NPreHeaderlines=find(I==1)-1;
        HeaderTextExact = eventstemp{I==1};
        
        % Import Table (now that we know where the data starts)
        EventTable_temp = readtable([directory fname],'Delimiter','\t','ReadVariableNames',1,'HeaderLines',NPreHeaderlines);
        
        StartTimeTemp = seconds(timeofday(EventTable_temp.Time_hh_mm_ss_)); % seconds from midnight
        [~,MorningStartIdx] = min(StartTimeTemp); % min will be time at or immediately after midnight
        StartT = [StartTimeTemp(1:MorningStartIdx-1); StartTimeTemp(MorningStartIdx:end) + 86400]; %add a full day to morning time difference
        %StartT = EventTable_temp.Time_hh_mm_ss_;
        
        
        EventDuration = EventTable_temp.Duration_s_;
        EventStart = StartT;
        EventStop=EventStart+seconds(EventDuration);
        EventName=EventTable_temp.Event;
        
        EventTable = table(EventName,EventStart,EventStop,EventDuration);
        SleepTable = EventTable;
        
        %W=0,N1=1,N2=2,N3=3,R=5
        EventCodesListUnique = unique(EventName);
        EventLabelOptions= {...
            'AROUSAL',...
            'APNEA-CENTRAL',...
            'APNEA-MIXED',...
            'APNEA-OBSTRUCTIVE',...
            'HYPOPNEA',...
            };
        EventCodeOptions=[1 3 5 2 4];
        EventExact1orStartsWith2 = [2 1 1 1 2];
        
        displaytext=['Get Hypnogram data: ' fname];
        if ~isempty(handletext), disp(displaytext); set(handletext,'String',displaytext); drawnow; end
        
        SleepCategoriesLabels= {...
            'SLEEP-S0',...
            'SLEEP-REM',...
            'SLEEP-S1',...
            'SLEEP-S2',...
            'SLEEP-S3',...
            };
        SleepCodes = [4 3 2 1 0]; %Terrill code (%.MAT: 0=N3,1=N2,2=N1,3=R,4=W)
        
        %get codes
        I = string(SleepTable.EventName)==SleepCategoriesLabels;
        I2 = sum(I,2)>0;
        temp = I.*SleepCodes;
        temp = max(temp')';
        temp(~I2) = [];
        SleepTable = SleepTable(I2,:);
        SleepTable.codes = temp;
        
        %         % convert EndTime
        %         [yr, mn, dy, hh, mm, ss] = datevec(EndTime/86400);
        %         EndTime_ = datetime(yr, mn, dy, hh, mm, ss);
        %
        temp = SleepTable.EventStart(2:end) - SleepTable.EventStart(1:end-1);
        %templast = timeofday(EndTime_)-timeofday(SleepTable.EventStart(end));
        templast = (EndTime)-(SleepTable.EventStart(end));
        %SleepTable.EventDuration(1:end-1) = seconds(temp);
        SleepTable.EventDuration(1:end-1) = temp;
        %SleepTable.EventDuration(end)=seconds(templast);
        SleepTable.EventDuration(end)=templast;
        SleepTable.EventStop=SleepTable.EventStart+SleepTable.EventDuration;
        clear temp;
        
        %         % StartTime is passed in at top of fn. need to turn into real time, not int
        %         [yr, mn, dy, hh, mm, ss] = datevec(StartTime/86400);
        %         StartTime_ = timeofday(datetime(yr, mn, dy, hh, mm, ss));
        
        NepochChanges = size(SleepTable,1);
        iEpochDuration=(round(SleepTable.EventDuration/30));
        iEpochStart=(1+round((SleepTable.EventStart-StartTime)/30)); %HypStart(Ilist(1))
        iEpochEnd=iEpochStart+iEpochDuration-1;
        checktemp=[iEpochStart iEpochDuration iEpochEnd SleepTable.codes];
        
        %first epochs filled with NaNs; starttimediff is just residual after rounding
        %
        % DLM says:  this doesn't always work properly.
        % The first scored epoch of sleep doesn't start at an epoch interval.
        % e.g. there may be 48 seconds of the study before
        % the first (30 seconds duration) epoch is recorded.
        % sure, we can fill the 30 seconds prior to epoch start with NaN,
        % but what about the preceding 18 seconds?
        % we also don't want the epochs to start at the study start,
        % because then they will all be offset by 12 seconds...
        
        starttimediff = SleepTable.EventStart(1)-StartTime - 30*(iEpochStart(1)-1); %first epochs filled with NaNs; starttimediff is just residual after rounding
        
        %starttimediffN=round((starttimediff/30))
        displaytext=['Start Time Hypnogram minus Signals: ' num2str(starttimediff) ' s'];
        if ~isempty(handletext), disp(displaytext); set(handletext,'String',displaytext); drawnow; end
        
        Nepochs=ceil((EndTime-StartTime)/30); %make this based on EndTime-StartTime...sum(iEpochDuration)
        Hypnogram=NaN*zeros(Nepochs,1);
        for i=1:length(iEpochDuration)
            Hypnogram(iEpochStart(i):iEpochEnd(i))=SleepTable.codes(i);
        end
        
        
case 'MignotAutostage' %code is under development (i.e. incomplete)
        %%
        directory = char(Filenames(6));
        fname = char(Filenames(3));
        
        clear eventstemp
        fid = fopen([directory '\' fname]);
        i=1;
        while 1
            eventstemp{i,:} = fgetl(fid);  %read in the data
            if eventstemp{i,:}==-1
                eventstemp(i,:) = [];
                break
            end
            i=i+1;
        end
        fclose(fid);   %close the file
        
        temp = str2double(eventstemp);
        Hypnogram = 8+0*temp;
        Hypnogram(temp==0)=4; %W
        Hypnogram(temp==5)=3; %REM
        Hypnogram(temp==1)=2; %N1
        Hypnogram(temp==2)=1; %N2
        Hypnogram(temp==3)=0; %N3
        Hypnogram(temp==4)=0; %N3
        Hypnogram_t = [StartTime:15:(StartTime + (length(Hypnogram)-1)*15)]';
        
%         Nepochs=ceil((EndTime-StartTime)/15); %make this based on EndTime-StartTime...sum(iEpochDuration)
%         Hypnogram=NaN*zeros(Nepochs,1);
%         for i=1:length(iEpochDuration)
%             Hypnogram(iEpochStart(i):iEpochEnd(i))=SleepTable.codes(i);
%         end
%         
%         
%         
%         HeaderTextList = {'Event','Duration[s]'};
%         I=[];
%         for ii=1:length(HeaderTextList)
%             I=[I,contains(eventstemp,HeaderTextList{ii})];
%         end
%         I=sum(I,2)/length(HeaderTextList);
%         NPreHeaderlines=find(I==1)-1;
%         HeaderTextExact = eventstemp{I==1};
%         
%         % Import Table (now that we know where the data starts)
%         EventTable_temp = readtable([directory fname],'Delimiter','\t','ReadVariableNames',1,'HeaderLines',NPreHeaderlines);
%         
%         StartTimeTemp = seconds(timeofday(EventTable_temp.Time_hh_mm_ss_)); % seconds from midnight
%         [~,MorningStartIdx] = min(StartTimeTemp); % min will be time at or immediately after midnight
%         StartT = [StartTimeTemp(1:MorningStartIdx-1); StartTimeTemp(MorningStartIdx:end) + 86400]; %add a full day to morning time difference
%         %StartT = EventTable_temp.Time_hh_mm_ss_;
%         
%         
%         EventDuration = EventTable_temp.Duration_s_;
%         EventStart = StartT;
%         EventStop=EventStart+seconds(EventDuration);
%         EventName=EventTable_temp.Event;
%         
%         EventTable = table(EventName,EventStart,EventStop,EventDuration);
%         SleepTable = EventTable;
%         
%         %W=0,N1=1,N2=2,N3=3,R=5
%         EventCodesListUnique = unique(EventName);
%         EventLabelOptions= {...
%             'AROUSAL',...
%             'APNEA-CENTRAL',...
%             'APNEA-MIXED',...
%             'APNEA-OBSTRUCTIVE',...
%             'HYPOPNEA',...
%             };
%         EventCodeOptions=[1 3 5 2 4];
%         EventExact1orStartsWith2 = [2 1 1 1 2];
%         
%         displaytext=['Get Hypnogram data: ' fname];
%         if ~isempty(handletext), disp(displaytext); set(handletext,'String',displaytext); drawnow; end
%         
%         SleepCategoriesLabels= {...
%             'SLEEP-S0',...
%             'SLEEP-REM',...
%             'SLEEP-S1',...
%             'SLEEP-S2',...
%             'SLEEP-S3',...
%             };
%         SleepCodes = [4 3 2 1 0]; %Terrill code (%.MAT: 0=N3,1=N2,2=N1,3=R,4=W)
%         
%         %get codes
%         I = string(SleepTable.EventName)==SleepCategoriesLabels;
%         I2 = sum(I,2)>0;
%         temp = I.*SleepCodes;
%         temp = max(temp')';
%         temp(~I2) = [];
%         SleepTable = SleepTable(I2,:);
%         SleepTable.codes = temp;
%         
%         %         % convert EndTime
%         %         [yr, mn, dy, hh, mm, ss] = datevec(EndTime/86400);
%         %         EndTime_ = datetime(yr, mn, dy, hh, mm, ss);
%         %
%         temp = SleepTable.EventStart(2:end) - SleepTable.EventStart(1:end-1);
%         %templast = timeofday(EndTime_)-timeofday(SleepTable.EventStart(end));
%         templast = (EndTime)-(SleepTable.EventStart(end));
%         %SleepTable.EventDuration(1:end-1) = seconds(temp);
%         SleepTable.EventDuration(1:end-1) = temp;
%         %SleepTable.EventDuration(end)=seconds(templast);
%         SleepTable.EventDuration(end)=templast;
%         SleepTable.EventStop=SleepTable.EventStart+SleepTable.EventDuration;
%         clear temp;
%         
%         %         % StartTime is passed in at top of fn. need to turn into real time, not int
%         %         [yr, mn, dy, hh, mm, ss] = datevec(StartTime/86400);
%         %         StartTime_ = timeofday(datetime(yr, mn, dy, hh, mm, ss));
%         
%         NepochChanges = size(SleepTable,1);
%         iEpochDuration=(round(SleepTable.EventDuration/30));
%         iEpochStart=(1+round((SleepTable.EventStart-StartTime)/30)); %HypStart(Ilist(1))
%         iEpochEnd=iEpochStart+iEpochDuration-1;
%         checktemp=[iEpochStart iEpochDuration iEpochEnd SleepTable.codes];
%         
%         %first epochs filled with NaNs; starttimediff is just residual after rounding
%         %
%         % DLM says:  this doesn't always work properly.
%         % The first scored epoch of sleep doesn't start at an epoch interval.
%         % e.g. there may be 48 seconds of the study before
%         % the first (30 seconds duration) epoch is recorded.
%         % sure, we can fill the 30 seconds prior to epoch start with NaN,
%         % but what about the preceding 18 seconds?
%         % we also don't want the epochs to start at the study start,
%         % because then they will all be offset by 12 seconds...
%         
%         starttimediff = SleepTable.EventStart(1)-StartTime - 30*(iEpochStart(1)-1); %first epochs filled with NaNs; starttimediff is just residual after rounding
%         
%         %starttimediffN=round((starttimediff/30))
%         displaytext=['Start Time Hypnogram minus Signals: ' num2str(starttimediff) ' s'];
%         if ~isempty(handletext), disp(displaytext); set(handletext,'String',displaytext); drawnow; end
%         
%         Nepochs=ceil((EndTime-StartTime)/30); %make this based on EndTime-StartTime...sum(iEpochDuration)
%         Hypnogram=NaN*zeros(Nepochs,1);
%         for i=1:length(iEpochDuration)
%             Hypnogram(iEpochStart(i):iEpochEnd(i))=SleepTable.codes(i);
%         end
        
    case 'Deltamed' %code is untested
        %%
        directory = char(Filenames(5));
        fname = char(Filenames(2));
        
        clear eventstemp
        fid = fopen([directory '\' fname]);
        i=1;
        while 1
            eventstemp{i,:} = fgetl(fid);  %read in the data
            if eventstemp{i,:}==-1
                eventstemp(i,:) = [];
                break
            end
            i=i+1;
        end
        fclose(fid);   %close the file
        
        HeaderTextList = {'Durée','Annotation'};
        I=[];
        for ii=1:length(HeaderTextList)
            I=[I,contains(eventstemp,HeaderTextList{ii})];
        end
        I=sum(I,2)/length(HeaderTextList);
        NPreHeaderlines=find(I==1)-1;
        %         HeaderTextExact = eventstemp{I==1};
        
        % Use regex to format the data
        eventstempNew = eventstemp;
        for rownum = NPreHeaderlines+1:length(eventstemp)
            temprow = eventstempNew{rownum,1};
            % insert duration into spots with missing data
            newtemprow = regexprep(temprow,'\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s','     00:00:00  ');
            
            % Make 5 space gaps into 6 space gaps
            newtemprow = regexprep(newtemprow,'([0-9])(\s\s\s\s\s)([0-9])','$1      $3');
            
            % Make 4 space gaps into 6 space gaps
            newtemprow = regexprep(newtemprow,'([0-9])(\s\s\s\s)([0-9])','$1      $3');
            
            % Make 4 space gaps into 6 space gaps
            newtemprow = regexprep(newtemprow,'([0-9])(\s\s\s)([0-9])','$1      $3');
            
            %                 % Remove + signs
            %                 newtemprow = regexprep(newtemprow,' \+','');
            %
            % Specific for one of the datasets
            newtemprow = regexprep(newtemprow,'(JMB1)(\s\s)(JMB2)','$1$3');
            
            % Replace ALL single spaced gaps -> replaces two commented
            % out secions above
            % finds spaces with NO space behind and NO space
            % ahead
            newtemprow = regexprep(newtemprow,'(?<!\s)(\s)(?!\s)','');
            
            % Replace h, m with colons
            newtemprow = regexprep(newtemprow,'([0-9])([hm])([0-9])','$1:$3');
            
            % Remove s
            newtemprow = regexprep(newtemprow,'([0-9])([s])','$1');
            
            % add back to events temp
            eventstempNew{rownum,1} = newtemprow;
        end
        
        % save as text
        fid2 = fopen([directory,fname(1:end-4),'_new.txt'],'w');
        for row = 1:length(eventstempNew)
            fprintf(fid2,'%s\n',eventstempNew{row,:});
        end
        fclose(fid2);
        
        % Import Table (now that we know where the data starts)
        EventTable_temp = readtable([directory,fname(1:end-4),'_new.txt'],'Delimiter',' ',...
            'ReadVariableNames',1,'HeaderLines',NPreHeaderlines,...
            'MultipleDelimsAsOne',true);
        % This code is replacing commented section above
        StartTimeTemp = seconds(duration(EventTable_temp.Heurer_elle)); % seconds from midnight
        [~,MorningStartIdx] = min(StartTimeTemp); % min will be time at or immediately after midnight
        StartT = [StartTimeTemp(1:MorningStartIdx-1); StartTimeTemp(MorningStartIdx:end) + 86400]; %add a full day to morning time difference
        
        EventDuration = seconds(duration(EventTable_temp.Dur_e));
        EventStart = StartT;
        EventStop=EventStart+EventDuration;
        EventName=EventTable_temp.Annotation;
        
        EventTable = table(EventName,EventStart,EventStop,EventDuration);
        SleepTable = EventTable;
        
        %W=0,N1=1,N2=2,N3=3,R=5
        EventCodesListUnique = unique(EventName);
        EventLabelOptions= {...
            'Arousal',...
            'ApnéeCentrale',...
            'ApnéeMixte',...
            'ApnéeObstructive',...
            'Hypopnée',...
            'hypopnéeCentrale',...
            };
        EventCodeOptions=[1 3 5 2 4 4];
        EventExact1orStartsWith2 = [2 2 2 2 2 2];
        
        displaytext=['Get Hypnogram data: ' fname];
        if ~isempty(handletext), disp(displaytext); set(handletext,'String',displaytext); drawnow; end
        
        % Import sleep stage table - alway starts at row 6
        SleepStageTable_temp = readtable([directory,fname(1:end-11),'.txt'],'Delimiter','\t',...
            'ReadVariableNames',0,'HeaderLines',5);
        
        % Get start time of sleep stages seconds
        StartTimeTemp = seconds(SleepStageTable_temp.Var1);
        [~,MorningStartIdx] = min(StartTimeTemp); % min will be time at or immediately after midnight
        StartT = [StartTimeTemp(1:MorningStartIdx-1); StartTimeTemp(MorningStartIdx:end) + 86400]; %add a full day to morning time difference
        StageStartTemp = StartT;
        
        SleepStageList = SleepStageTable_temp.Var2;
        StageDurationTemp = nan(size(StageStartTemp));
        SleepStageTable = table(SleepStageList,StageStartTemp,StageDurationTemp);
        SleepStageTable.Properties.VariableNames{'StageDurationTemp'} = 'EventDuration';
        SleepStageTable.Properties.VariableNames{'StageStartTemp'} = 'EventStart';
        
        SleepCategoriesLabels= {...
            'Veille',...
            'S. Paradoxal',...
            'Stade 1',...
            'Stade 2',...
            'Stade 3'...
            };
        SleepCodes = [4 3 2 1 0]; %Terrill code (%.MAT: 0=N3,1=N2,2=N1,3=R,4=W)
        
        %get codes
        I = string(SleepStageTable.SleepStageList)==SleepCategoriesLabels;
        I2 = sum(I,2)>0;
        temp = I.*SleepCodes;
        temp = max(temp')';
        SleepStageTable.codes = temp;
        
        SleepStageTable = SleepStageTable(I2,:);
        temp = SleepStageTable.EventStart(2:end) - SleepStageTable.EventStart(1:end-1);
        templast = EndTime-SleepStageTable.EventStart(end);
        SleepStageTable.EventDuration(1:end-1) =temp;
        SleepStageTable.EventDuration(end)=templast;
        SleepStageTable.EventStop=SleepStageTable.EventStart+SleepStageTable.EventDuration;
        clear temp;
        
        NepochChanges = size(SleepStageTable,1);
        iEpochDuration=(round(SleepStageTable.EventDuration/30));
        iEpochStart=(1+round((SleepStageTable.EventStart-StartTime)/30)); %HypStart(Ilist(1))
        iEpochEnd=iEpochStart+iEpochDuration-1;
        checktemp=[iEpochStart iEpochDuration iEpochEnd SleepStageTable.codes];
        
        %first epochs filled with NaNs; starttimediff is just residual after rounding
        starttimediff = SleepStageTable.EventStart(1)-StartTime - 30*(iEpochStart(1)-1); %first epochs filled with NaNs; starttimediff is just residual after rounding
        
        %starttimediffN=round((starttimediff/30))
        displaytext=['Start Time Hypnogram minus Signals: ' num2str(starttimediff) ' s'];
        if ~isempty(handletext), disp(displaytext); set(handletext,'String',displaytext); drawnow; end
        
        Nepochs=ceil((EndTime-StartTime)/30); %make this based on EndTime-StartTime...sum(iEpochDuration)
        Hypnogram=NaN*zeros(Nepochs,1);
        for i=1:length(iEpochDuration)
            Hypnogram(iEpochStart(i):iEpochEnd(i))=SleepStageTable.codes(i);
        end
        
    case 'SomnoStarEpochsXls'
        %%
        %hypnogram only here:
        directory = char(Filenames(6));
        fname = char(Filenames(3));
        
        displaytext=['Get Hypnogram data:' fname];
        if ~isempty(handletext), disp(displaytext); set(handletext,'String',displaytext); drawnow; end
        [~,~,raw]=xlsread([directory fname],1,'A2:B2000');
        EpochNumber = raw(:,1);
        EpochNumber = cell2mat(EpochNumber);
        Nepochs_ = size(EpochNumber,1) - find(isnan(flipud(EpochNumber)),1,'last')
        raw((Nepochs_+1):end,:)=[];
        
        EpochNumber((Nepochs_+1):end,:)=[];
        
        Compare  = (1:Nepochs_)';
        CheckNoMissingLines = sum(diff((abs(Compare - EpochNumber)))~=0)
        SleepCodeslist = raw(:,2);
        SleepCategoriesLabels= {...
            'U',...
            'W',...
            'N1',...
            'N2',...
            'N3',...
            'R',...
            };
        SleepCategoriesCodes = [8 4 2 1 -1 3]; %Terrill code (%.MAT: 0=N3,1=N2,2=N1,3=R,4=W)
        
        
        clear match
        for m=1:length(SleepCategoriesLabels)% each column is an event type
            match(:,m)=string(SleepCategoriesLabels(m))==string(SleepCodeslist');
        end
        temp = match.*SleepCategoriesCodes;
        temp(temp==0)=NaN;
        temp(temp==-1)=0;
        SleepCodes=min(temp')'; %just in case two criteria are detected, will leave the lowest code value in list
        
        Hypnogram = nan(EpochNumber(end),1);
        
        %first epochs filled with NaNs; starttimediff is just residual after rounding
        starttimediff = 0; %assumed
        
        Nepochs=ceil((EndTime-StartTime)/30); %make this based on EndTime-StartTime...sum(iEpochDuration)
        Hypnogram=NaN*zeros(Nepochs,1);
        for i=1:length(EpochNumber)
            Hypnogram(EpochNumber(i))=SleepCodes(i);
        end
        
    case 'GrassTwinTxt' % this code was written in March 2020, undergoing testing (DLM)
        %%
        directory = char(Filenames(5));
        fname = char(Filenames(2));
        
        displaytext=['Importing events from GrassTwin text file:' fname];
        if ~isempty(handletext), disp(displaytext); set(handletext,'String',displaytext); drawnow; end
        % using textscan here, because dlmream only works for numerical data
        str = fileread([directory, fname]);
        C = textscan(str,'%d16%s%s','Delimiter',',');
        
        if 0 % manually load, only for debugging
            str = fileread('G:\Shared drives\CCHMC\Source\LG_0001.txt');
            % C = textscan(str,'%d16%s%*[^\n]','Delimiter',','); % to read the first two chuncks, and then discard the rest of the line (testing purposes)
            C = textscan(str,'%d16%s%s','Delimiter',',');
        end
        
        % convert the cell arrays to something more friendly
        EpochNumber = C{1,1}; EventTime = C{1,2}; EventDetails = C{1,3};
        
        % turn EventTime into mat times, and call it EventStart
        EventStartT = NaN(size(EventTime,1),1);
        formatIn = 'HH:MM:SS.FFF';
        for i=1:size(EventTime,1)
            temp = mod(datenum(EventTime{i},formatIn),1)*86400;
            if temp<86400/2,temp=temp+86400; end
            EventStartT(i,1) = temp;
        end
        
        % then split EventDetails into Epochs and Events: (keep EventTime with each)
        % Epochs 'Stage - Type'
        % Events (everything else)
        % we don't know the size req'd for each list apriori, so this grows during the for loop (inneficient, but easy)
        EpochList=[];  EventList=[];
        for i=1:size(EventDetails,1)
            % if it startswith 'Stage - ' then add this to epoch list
            if startsWith(string(EventDetails(i)),'Stage - ')
                EpochList = [EpochList; [EpochNumber(i), EventDetails(i), EventStartT(i)]];
            else % otherwise, add it to event list
                EventList = [EventList; [EventStartT(i), EventDetails(i)]];
            end
        end
        
        clearvars C EpochNumber EventTime EventDetails EventStartT formatIn i temp
        % Convert EpochList to be the Hypnogram
        % first, assign Hypnogram number to stage text
        StageTextOptions = {'Stage - Wake','Stage - REM','Stage - N1','Stage - N2','Stage - N3','Stage - No Stage','Stage - 4'};
        StageTextCodes = [4,3,2,1,0,8,0]; %Terrill code: W=4,R=3,N1=2,N2=1,N3=0, unknown=NaN]
        % Note that CCH-MC infant studies score
        % "Stage - 4" for QS, and "Stage - REM" for AS
        % The current approach (c/o Brad Edwards) is set
        % "Stage - 4" as "Stage - N3", and leave REM as is.
        HypCodes = [];
        for i=1:size(EpochList,1)
            for m=1:length(StageTextOptions)
                if strcmp(EpochList{i,2},StageTextOptions{m})
                    HypCodes(i,:) = StageTextCodes(m);
                end
            end
        end
        
        % Not needed here so far because start epochs is the same as start
        % time for the EDF:
        FirstEpochStartTime = EpochList{1,3}
        StartTime
        if 0 %use as needed
            starttimediff = FirstEpochStartTime - StartTime;
        end
        % second, expand/interp to full length (i.e. Nepochs long)
        % i.e. stretch the HypCodes between KnownEpochs
        KnownEpochs = cell2mat(EpochList(:,1));
        HypnogEpochs = 1:1:EpochList{end,1};
        Hypnogram = interp1(single(KnownEpochs),HypCodes,single(HypnogEpochs'),'previous','extrap'); % ToDo: confirm this works for last sample
        
        % then process the Events data
        % Event start is just taken from full list.
        % Event duration and code need teasing out from string, based on type
        % Ar 'Arousal - Dur: xx.x sec. - Type'
        % Resp 'Respiratory Event - Dur: xx.x sec. - Type - Desat xx.x %
        % and assign the event category code as we go
        EventStart = cell2mat(EventList(:,1));
        EventDuration = NaN(size(EventList,1),1);
        clearvars EventName
        for i=1:size(EventList,1)
            str = string(EventList{i,2});
            if startsWith(str,'Arousal - ') || startsWith(str,'* Arousal - ') % multiple options for how this string appears
                % read the first three pieces, and discard the rest of line
                C = textscan(str,'%s%s%s%*[^\n]','Delimiter','-');
                str2 = string(C{1,2}); % get duration
                t1 = regexp(str2,'\d+\.?\d+','match'); % tease out the actual numbers from the string
                EventDuration(i) = str2double(t1); % set event duration
                EventCodesStrTemp = string(C{1,3}); % get event details,
                EventCodesStr = deblank(EventCodesStrTemp); % remove trailing whitespace
                EventName(i) = strjoin(['Arousal -' EventCodesStr]); % prepend 'Arousal - ' and add to event list
                
            elseif startsWith(str,'Respiratory Event - ')
                % read the first three pieces, and discard the rest of line
                C = textscan(str,'%s%s%s%*[^\n]','Delimiter','-');
                str2 = string(C{1,2}); % get duration
                t1 = regexp(str2,'\d+\.?\d+','match'); % tease out the actual numbers from the string
                EventDuration(i) = str2double(t1); % set event duration
                EventCodesStr = string(C{1,3}); % get event details
                EventName(i) = deblank(EventCodesStr); % remove trailing whitespace and add to event list
                
            else
                EventName(i) = "";
                EventDuration(i) = NaN;
                % something different
                % keyboard >>dbquit
            end
        end
        
        EventLabelOptions= {...
            'Arousal',...
            'Central Apnea',...
            'Central Hypopnea',...
            'Mixed Apnea',...
            'Obstructive Apnea',...
            'Obstructive Hypopnea',...
            };
        EventCodeOptions=[1 3 3 5 2 4];
        EventExact1orStartsWith2 = [2 1 1 1 1 1];
        EventName = cellstr(EventName');

        % end of GrassTwinTxt events import
        
        %     case 'SandmanPlusLabview' %%
    case 'ZMachineSynergy' % started by SS 2/16/22
        %%
        directory = char(Filenames(5));
        fname = char(Filenames(2));
        
        displaytext=['Importing events from ZMachine Synergy text file:' fname];
        if ~isempty(handletext), disp(displaytext); set(handletext,'String',displaytext); drawnow; end
        % using textscan here, because dlmream only works for numerical data
        str = fileread([directory, fname]);
        Table1 = readtable([directory, fname]);
        unique(Table1.Annotation)
        
        Table1.HypCodes = nan(height(Table1),1);
        
        %clearvars C EpochNumber EventTime EventDetails EventStartT formatIn i temp
        % Convert EpochList to be the Hypnogram
        % first, assign Hypnogram number to stage text
        StageTextOptions = {'Sleep Stage W','Sleep Stage R','Sleep Stage N1','Sleep Stage N2','Sleep Stage N3','Sleep Stage ?'};
        StageTextCodes = [4,3,2,1,0,8,0]; %Terrill code: W=4,R=3,N1=2,N2=1,N3=0, unknown=NaN]
        % Note that CCH-MC infant studies score
        % "Stage - 4" for QS, and "Stage - REM" for AS
        % The current approach (c/o Brad Edwards) is set
        % "Stage - 4" as "Stage - N3", and leave REM as is.
        
        for i=1:height(Table1)
            for m=1:length(StageTextOptions)
                if strcmp(Table1.Annotation{i},StageTextOptions{m})
                    Table1.HypCodes(i,:) = StageTextCodes(m);
                end
            end
        end
        
        EpochsTable = Table1(~isnan(Table1.HypCodes),:);
        
        % Not needed here so far because start epochs is the same as start
        % time for the EDF:
        FirstEpochStartTime = EpochsTable{1,3}
        StartTime
        if 0 %use as needed
            starttimediff = FirstEpochStartTime - StartTime;
        end
        % second, expand/interp to full length (i.e. Nepochs long)
        % i.e. stretch the HypCodes between KnownEpochs
        KnownEpochs = cell2mat(EpochList(:,1));
        HypnogEpochs = 1:1:EpochList{end,1};
        Hypnogram = interp1(single(KnownEpochs),HypCodes,single(HypnogEpochs'),'previous','extrap'); % ToDo: confirm this works for last sample
        
        % then process the Events data
        % Event start is just taken from full list.
        % Event duration and code need teasing out from string, based on type
        % Ar 'Arousal - Dur: xx.x sec. - Type'
        % Resp 'Respiratory Event - Dur: xx.x sec. - Type - Desat xx.x %
        % and assign the event category code as we go
        EventStart = cell2mat(EventList(:,1));
        EventDuration = NaN(size(EventList,1),1);
        clearvars EventName
        for i=1:size(EventList,1)
            str = string(EventList{i,2});
            if startsWith(str,'Arousal - ') || startsWith(str,'* Arousal - ') % multiple options for how this string appears
                % read the first three pieces, and discard the rest of line
                C = textscan(str,'%s%s%s%*[^\n]','Delimiter','-');
                str2 = string(C{1,2}); % get duration
                t1 = regexp(str2,'\d+\.?\d+','match'); % tease out the actual numbers from the string
                EventDuration(i) = str2double(t1); % set event duration
                EventCodesStrTemp = string(C{1,3}); % get event details,
                EventCodesStr = deblank(EventCodesStrTemp); % remove trailing whitespace
                EventName(i) = strjoin(['Arousal -' EventCodesStr]); % prepend 'Arousal - ' and add to event list
                
            elseif startsWith(str,'Respiratory Event - ')
                % read the first three pieces, and discard the rest of line
                C = textscan(str,'%s%s%s%*[^\n]','Delimiter','-');
                str2 = string(C{1,2}); % get duration
                t1 = regexp(str2,'\d+\.?\d+','match'); % tease out the actual numbers from the string
                EventDuration(i) = str2double(t1); % set event duration
                EventCodesStr = string(C{1,3}); % get event details
                EventName(i) = deblank(EventCodesStr); % remove trailing whitespace and add to event list
                
            else
                EventName(i) = "";
                EventDuration(i) = NaN;
                % something different
                % keyboard >>dbquit
            end
        end
        
        EventLabelOptions= {...
            'Arousal',...
            'Central Apnea',...
            'Central Hypopnea',...
            'Mixed Apnea',...
            'Obstructive Apnea',...
            'Obstructive Hypopnea',...
            };
        EventCodeOptions=[1 3 3 5 2 4];
        EventExact1orStartsWith2 = [2 1 1 1 1 1];
        EventName = cellstr(EventName');

        % end of GrassTwinTxt events import
        
        %     case 'SandmanPlusLabview' %%
    case 'Sandman' %%
        %%
        directory = char(Filenames(5));
        fname = char(Filenames(1,2));
        
        displaytext=['Importing events from Sandman text file:' fname];
        if ~isempty(handletext), disp(displaytext); set(handletext,'String',displaytext); drawnow; end
        
        Sandmanfilename=[directory fname];
        
        %[EpochsTemp] = textscan(Sfileid,'%s');
        
        clear eventstemp
        fid = fopen([directory '\' fname]);
        i=1;
        while 1
            eventstemp{i,:} = fgetl(fid);  %read in the data
            if eventstemp{i,:}==-1
                eventstemp(i,:) = [];
                break
            end
            i=i+1;
        end
        fclose(fid);   %close the file
        
        HeaderTextList = {'Epoch','Event'};
        I=[];
        for ii=1:length(HeaderTextList)
            I=[I,contains(eventstemp,HeaderTextList{ii})];
        end
        I=sum(I,2)/length(HeaderTextList);
        NPreHeaderlines=find(I==1)-1;
        HeaderTextExact = eventstemp{I==1};
        
        % Import Table (now that we know where the data starts)
        EpochEventT = readtable([directory fname],'Delimiter','\t','ReadVariableNames',1,'HeaderLines',NPreHeaderlines);
        %
        %
        %
        %
        %Sfileid=fopen(Sandmanfilename);
        %[EpochEventT] = textscan(Sfileid,'%u %s %s %f %f %f %f','whitespace','\t','headerlines',NPreHeaderlines);
        %fclose(Sfileid);
        %         SandmandEpochs=Epochs(:,1:4); % first col: epoch #, 2: events list, 3:abs start time, 4:dur
        
        %EpochNumber = Epochs.Epoch;
        %EventList = Epochs.Event;
        %EventTime = Epochs.StartTime;
        %EventDur=Epochs.Duration;
        
        clear EventStartTime
        for j = 1:length(EpochEventT.Epoch)
            % only_event_time = strsplit(EventTime(j));
            event_time = mod(datenum(EpochEventT.StartTime(j)),1)*86400;
            if event_time<43200; event_time=event_time+(86400/2); end % converting to 24hr format for anything less than 12am
            if event_time<StartTime % eg as in 1 am
                event_time=event_time+86400/2;
                EventStartTime(j,1)=event_time;
            else
                EventStartTime(j,1)=event_time;
            end
            
        end
        %EpochNumber=double(EpochNumber);
        %EpochEventT=table(EpochNumber,EventList,EventStartTime,EventDur);
        EpochEventT.EventStartTime=EventStartTime;
        
        StageLabelOptions = {'Wake','REM','NREM 1','NREM 2','NREM 3'}';
        StageCodeOptions = [4 3 2 1 -1] ; %can not be zeros here
        
        clear match1 match2
        EventName = EpochEventT.Event;
        for m=1:length(StageLabelOptions)% each column is an event type
            match1(:,m)=strcmpi(string(EventName'),string(StageLabelOptions(m)));
            match2(:,m)=startsWith(string(EventName'),string(StageLabelOptions(m)));
        end
        match = 1*( match1>0 | match2>0);
        temp = match.*StageCodeOptions;
        temp(temp==0)=NaN;
        temp(temp==-1)=0;
        EpochEventT.EpochCodes=min(temp')'; %just in case two criteria are detected, will leave the lowest code value in list
        
        HypnogramT = EpochEventT(~isnan(EpochEventT.EpochCodes),:);
        HypnogramT.DurationDiff = [diff(HypnogramT.EventStartTime);NaN];
        if ~isnan(HypnogramT.Duration(end))
            HypnogramT.DurationDiff(end) = HypnogramT.Duration(end);
        else
            HypnogramT.DurationDiff(end) = 30;
        end
        HypnogramT.EventEndTime = HypnogramT.EventStartTime+HypnogramT.DurationDiff;
        
        Hypnogram_t = (HypnogramT.EventStartTime(1)-30:30:HypnogramT.EventEndTime(end))';
        tempH = [8;HypnogramT.EpochCodes;8];
        tempT = [HypnogramT.EventStartTime(1)-30;HypnogramT.EventStartTime;HypnogramT.EventEndTime(end)];
        Hypnogram = interp1(tempT,tempH,Hypnogram_t+15,'previous','extrap');
        Hypnogram(1)=[];
        Hypnogram_t(1)=[];
        
        %         figure(99);
        %         plot(tempT,tempH,'.'); hold on
        %         plot(Hypnogram_t,Hypnogram,'.-'); hold on
        % EventName,EventStart,EventDuration,EventLabelOptions,EventCodeOptions,EventExact1orStartsWith2,Hypnogram
        starttimediff = Hypnogram_t(1)-StartTime;
        
        
        % respiratory events  arousals
        EventLabelOptions= {...
            'Arousal',...
            'Central Apnea',...
            'Mixed Apnea',...
            'Obstructive Apnea',...
            'Obstructive Hypopnea',...
            'RERA',...
            };
        EventCodeOptions=[1 3 5 2 4 12];
        EventExact1orStartsWith2 = [2 1 2 2 2 1];
        EventList=EpochEventT.Event;
        Idx=zeros(length(EventList),1);
        for ii=1:length(EventLabelOptions)
            StageTF = startsWith(EventList,EventLabelOptions{ii},'IgnoreCase',true);
            Idx=Idx|StageTF;
        end
        EpochEventT.RespEvent=Idx; % just to double check
        EventName=EpochEventT.Event(EpochEventT.RespEvent==1);
        EventStart=EpochEventT.EventStartTime(EpochEventT.RespEvent==1);
        EventDuration=EpochEventT.Duration(EpochEventT.RespEvent==1);
        
        PositionT = EpochEventT(:,{'Event','EventStartTime','Duration'}); %is actually the whole event table, some rows have position event start times
        PositionT.Properties.VariableNames(1:2) = {'Position','Time'};
        PosCodes = {'Supine','Left','Right','Prone','Front','Upright'};
        I = sum(PositionT.Position == string(PosCodes),2)>0;
        PositionT = PositionT(I,:);
        
        NewPosCodes = {'Supine','Left','Right','Prone','Unknown','Upright'}; %from database, hardcoded
        PosCode = [1 2 3 4 5 6]; %from database, hardcoded
        PositionT.PositionOrig = PositionT.Position;
        PositionT.Codes = nan(height(PositionT),1);
        for j=1:length(PosCodes)
            for i=1:height(PositionT)
                if PositionT.Position(i)==string(PosCodes(j))
                    PositionT.Position{i}=NewPosCodes{j};
                    PositionT.Codes(i)=PosCode(j);
                end
            end
        end
        PositionT.DurationDiff = [diff(PositionT.Time);NaN];
        
        % end of Sandman events import
    case 'NKxdf' %AZversion
        %% ProfusionXml Hypnogram
        directory = char(Filenames(5));
        fname = char(Filenames(2));
        
        %Xml events
        displaytext=['Get Events data:' fname];
        if ~isempty(handletext), disp(displaytext); set(handletext,'String',displaytext); drawnow; end
        
        S=xml2struct([directory fname]);
        
        if isfield(S,'html') && isfield(S.html,'body')
            S=S.html.body;
        end
        if isfield(S,'xdf_colon_openxdf')
            fieldname = 'xdf_colon_OpenXDF';
            S = setfield(S,fieldname,getfield(S,lower('xdf_colon_OpenXDF')));
            %add removal of old field here
        end
        
        %Position List to Table
        NoteLog = S.xdf_colon_OpenXDF.xdf_colon_NoteLog.xdf_colon_Note ; %.xdf_colon_Note{1, 21}.xdf_colon_NoteText
        NotelogText=[];
        NoteLogTimes=[];
        for i=1:length(NoteLog)
            NotelogText{i,1}=NoteLog{i}.xdf_colon_NoteText.Text;
            
            timestr = NoteLog{i}.xdf_colon_Time.Text(12:end);
            event_time = round(mod(datenum(timestr(1:5),'HH:MM'),1)*86400) + str2num(timestr(7:end)); %hopeless datenum code! hacks make it work for decimals.
            if event_time<43200; event_time = event_time + 86400; end % add a day for data between midnight and midday
            NoteLogTimes(i,1)=event_time;
        end
        PositionT = table(NotelogText,NoteLogTimes); %is actually the whole event table, some rows have position event start times
        PositionT.Properties.VariableNames = {'Position','Time'};
        PosCodes = {'Supine','Left Side','Right Side','Prone','Unknown','Upright'};
        I = sum(PositionT.Position == string(PosCodes),2)>0;
        PositionT = PositionT(I,:);
        
        NewPosCodes = {'Supine','Left','Right','Prone','Unknown','Upright'}; %from database, hardcoded
        PosCode = [1 2 3 4 5 6]; %from database, hardcoded
        
        PositionT.PositionOrig = PositionT.Position;
        PositionT.Codes = nan(height(PositionT),1);
        for j=1:length(PosCodes)
            for i=1:height(PositionT)
                if PositionT.Position(i)==string(PosCodes(j))
                    PositionT.Position{i}=NewPosCodes{j};
                    PositionT.Codes(i)=PosCode(j);
                end
            end
        end
        
        %Events
        temp = S.xdf_colon_OpenXDF.xdf_colon_ScoringResults.xdf_colon_Scorers.xdf_colon_Scorer;
        Nscorers = length(temp);
        disp(['No. of scorers: ' num2str(Nscorers)])
        for i=1:Nscorers
            if length(temp{i}.xdf_colon_LastName.Text)>0
                ScorerNames{i} = [temp{i}.xdf_colon_FirstName.Text ' ' temp{i}.xdf_colon_LastName.Text];
            else
                ScorerNames{i} = [temp{i}.xdf_colon_FirstName.Text];
            end
        end
        
        %CandidateScorer = find((ScorerNames~="Polysmith")==1);
        disp(['Scorer Names: ' strjoin(ScorerNames,', ')])
        %ScorerJ = CandidateScorer(end);
        
        ScorerInclList = {};
        ScorerExclList = {'Polysmith'};
        if isfield(settings,'ScorerInclList')
            ScorerInclList = settings.ScorerInclList;
        end
        if isfield(settings,'ScorerExclList')
            ScorerExclList = settings.ScorerExclList;
        end
        
        CandidateScorer = find(sum(ScorerNames==string(ScorerExclList'),1)~=1);
        if isempty(ScorerInclList)
            PriorityScorer=[];
        else
        PriorityScorer = sum((ScorerNames==string(ScorerInclList')) .* repmat(1:length(ScorerNames),length(ScorerInclList),1),2)';
            %highest on list is priortiy
        end
        if ~isempty(PriorityScorer)
            ScorerJ = PriorityScorer(1);
            disp(['Note--Found Priority Scorer: ' ScorerNames{CandidateScorer(end)}])
        else
            if length(CandidateScorer)==1
            ScorerJ = CandidateScorer;
            disp(['Warning--No Priority Scorer, Using Only Available Scorer: ' ScorerNames{CandidateScorer(end)}]);
            else
            ScorerJ = CandidateScorer(end);
            disp(['Warning--No Priority Scorer, Using Last Scorer in List: ' ScorerNames{CandidateScorer(end)}]);    
            end
        end
        disp(['Chosen Scorer: ' ScorerNames{ScorerJ}])
        
        ScoringStruct = S.xdf_colon_OpenXDF.xdf_colon_ScoringResults.xdf_colon_Scorers.xdf_colon_Scorer{ScorerJ};
        
        %initialize empty
        EventStart=[];
        EventManvAuto=[];
        EventName=[];
        EventDuration=[];
        EventClasses = {'Apnea','Hypopnea','RERA','Microarousal','CustomEvent'};
        NTIorXDF = {'xdf','xdf','xdf','xdf','nti'};
        for k=1:length(EventClasses)
            try
                Stemp = getfield(ScoringStruct,[NTIorXDF{k} '_colon_' EventClasses{k} 's']);
                Stemp2 = getfield(Stemp,[NTIorXDF{k} '_colon_' EventClasses{k}]); %will crash here if no events in category
                Neventstemp = length(Stemp2);
                    if Neventstemp == 1 %if only a single event (within class e.g. apnea), it wasn't being brought in as a cell, fixed here
                        temp = Stemp2;
                        clear Stemp2;
                        Stemp2{1} = temp;
                    end
            catch me
                continue
            end
            for j=1:Neventstemp
                %starttime
                timestr = Stemp2{j}.xdf_colon_Time.Text(12:end);
                event_time = round(mod(datenum(timestr(1:5),'HH:MM'),1)*86400) + str2num(timestr(7:end)); %hopeless datenum code! hacks make it work for decimals.
                if event_time<43200; event_time = event_time + 86400; end % add a day for data between midnight and midday
                EventStart(end+1,1)=event_time;
                %timecheck = TimeinsecToDatetime(event_time);
                EventManvAuto(end+1,1) = strcmp(Stemp2{j}.nti_colon_Manual.Text,'true');
                EventDuration(end+1,1) = str2num(Stemp2{j}.xdf_colon_Duration.Text);
                try
                    classstr = [' ' Stemp2{j}.xdf_colon_Class.Text];
                catch
                    classstr = '';
                end
                EventName{end+1,1} = [EventClasses{k} classstr];
            end
        end
        
        %for arousals of zero seconds:
        EventDuration(EventDuration<3)=3;
        
        PrelimEventTable = table(EventName,EventStart,EventDuration,EventManvAuto);
        
        Nevents = length(EventName); %needed?
        
        %
        
        %         if ~isfield(S.CMPStudyConfig.ScoredEvents,'ScoredEvent')
        %             displaytext=['Warning: No scored events'];
        %             if ~isempty(handletext), disp(displaytext); set(handletext,'String',displaytext); drawnow; end
        %         end
        
        EventCodesListUnique = unique(EventName);
        
        % normal procedures
        EventLabelOptions= {...
            'Microarousal',...
            'Apnea central',...
            'Apnea mixed',...
            'Apnea obstructive',...
            'Hypopnea obstructive',... %
            'Hypopnea central',... %
            'CustomEvent',... %assumed to be the hypopnea 4% category based on a single example study.
            'RERA ' % has space in it
            };
        EventCodeOptions=[1 3 5 2 4 6 4 12];
        EventExact1orStartsWith2 = [2 1 1 1 1 1 1 1];
        
        
        
        
        displaytext=['Get Hypnogram data: ' fname];
        SStages = ScoringStruct.xdf_colon_SleepStages.xdf_colon_SleepStage;
        
        if ~isempty(handletext), disp(displaytext); set(handletext,'String',displaytext); drawnow; end
        
        %Code assumes that there is one entry per 30-s epoch, not
        %necessarily in order, starting at beginning of the EDF file:
        HypnogramLabel=[];
        HypnogramEpochNumber=[];
        HypnogramInfStage=[];
        HypnogramStgStrength=[];
        for i=1:length(SStages)
            temp = SStages{i};
            HypnogramLabel{end+1,1} = temp.xdf_colon_Stage.Text;
            HypnogramEpochNumber{end+1,1} = str2num(temp.xdf_colon_EpochNumber.Text);
            HypnogramInfStage{end+1,1} = temp.nti_colon_InfStage.Text;
            HypnogramStgStrength{end+1,1} = temp.nti_colon_StgStrength.Text;
        end %W=0,N1=1,N2=2,N3=3,R=5
        
        HypnogramTable = table(HypnogramLabel,HypnogramEpochNumber,HypnogramInfStage,HypnogramStgStrength);
        HypnogramTable = sortrows(HypnogramTable,'HypnogramEpochNumber');
        %Convert to Terrill code %W=-1,N1=1,N2=2,N3=3,R=5
        Hypnogram = nan(height(HypnogramTable),1);%
        Hypnogram(strcmp(HypnogramTable.HypnogramLabel,'W'))=4;
        Hypnogram(strcmp(HypnogramTable.HypnogramLabel,'1'))=2;
        Hypnogram(strcmp(HypnogramTable.HypnogramLabel,'2'))=1;
        Hypnogram(strcmp(HypnogramTable.HypnogramLabel,'3'))=0;
        Hypnogram(strcmp(HypnogramTable.HypnogramLabel,'4'))=0;
        Hypnogram(strcmp(HypnogramTable.HypnogramLabel,'R'))=3;
        %%
    case 'NKxdfB'
        %% NKxdf modified to handle all lowercase xml/struct fieldnames
        directory = char(Filenames(5));
        fname = char(Filenames(2));
        
        %Xml events
        displaytext=['Get Events data:' fname];
        if ~isempty(handletext), disp(displaytext); set(handletext,'String',displaytext); drawnow; end
        
        S=xml2struct([directory fname]);
        
        if isfield(S,'html') && isfield(S.html,'body')
            S=S.html.body;
        end
        
        %Position List to Table
        NoteLog = S.xdf_colon_openxdf.xdf_colon_notelog.xdf_colon_note ; %.xdf_colon_Note{1, 21}.xdf_colon_NoteText
        NotelogText=[];
        NoteLogTimes=[];
        for i=1:length(NoteLog)
            NotelogText{i,1}=strtrim(NoteLog{i}.xdf_colon_notetext.Text);
           
            timestr = strtrim(NoteLog{i}.xdf_colon_time.Text);
            timestr = timestr(12:end);
            event_time = round(mod(datenum(timestr(1:5),'HH:MM'),1)*86400) + str2num(timestr(7:end)); %hopeless datenum code! hacks make it work for decimals.
            if event_time<43200; event_time = event_time + 86400; end % add a day for data between midnight and midday
            NoteLogTimes(i,1)=event_time;
        end
        PositionT = table(NotelogText,NoteLogTimes); %is actually the whole event table, some rows have position event start times
        PositionT.Properties.VariableNames = {'Position','Time'};
        PosCodes = {'Supine','Left Side','Right Side','Prone','Unknown','Upright'};
        I = sum(PositionT.Position == string(PosCodes),2)>0;
        PositionT = PositionT(I,:);
        
        NewPosCodes = {'Supine','Left','Right','Prone','Unknown','Upright'}; %from database, hardcoded
        PosCode = [1 2 3 4 5 6]; %from database, hardcoded
        
        PositionT.PositionOrig = PositionT.Position;
        PositionT.Codes = nan(height(PositionT),1);
        for j=1:length(PosCodes)
            for i=1:height(PositionT)
                if PositionT.Position(i)==string(PosCodes(j))
                    PositionT.Position{i}=NewPosCodes{j};
                    PositionT.Codes(i)=PosCode(j);
                end
            end
        end
        
        %Events
        temp = S.xdf_colon_openxdf.xdf_colon_scoringresults.xdf_colon_scorers.xdf_colon_scorer;
        Nscorers = length(temp);
        for i=1:Nscorers
            ScorerNames{i} = strtrim(temp{i}.xdf_colon_firstname.Text);
        end
        ScorerJ = Nscorers;
        
        CandidateScorer = find((ScorerNames~="Polysmith")==1);
        
        clear ScorerTable
        Nvalid = length(CandidateScorer);
        %run loop through all scorers who are not called "Polysmith";
        %choose the scorer with the most event markings
        %in our brief experience, invalid scorers have had no events at all.
        %may need another strategy if there are often two scorers with reasonable event lists.
        for xx=1:(Nvalid+1) 
        
            ScorerJ = CandidateScorer(xx);
        %disp(['Warning: Scorer assumed as the appropriate scorer: ' ScorerNames{ScorerJ}])
        clear temp
        
        ScoringStruct = S.xdf_colon_openxdf.xdf_colon_scoringresults.xdf_colon_scorers.xdf_colon_scorer{ScorerJ};
        
        %initialize empty
        EventStart=[];
        EventManvAuto=[];
        EventName=[];
        EventDuration=[];
        EventClasses = {'apnea','hypopnea','rera','microarousal','customevent'};
        NTIorXDF = {'xdf','xdf','xdf','xdf','nti'};
        for k=1:length(EventClasses)
            try
                Stemp = getfield(ScoringStruct,[NTIorXDF{k} '_colon_' EventClasses{k} 's']);
                Stemp2 = getfield(Stemp,[NTIorXDF{k} '_colon_' EventClasses{k}]); %will crash here if no events in category
                Neventstemp = length(Stemp2);
            catch me
                continue
            end
            for j=1:Neventstemp
                %starttime
                timestr = strtrim(Stemp2{j}.xdf_colon_time.Text);
                 timestr = timestr(12:end);
                event_time = round(mod(datenum(timestr(1:5),'HH:MM'),1)*86400) + str2num(timestr(7:end)); %hopeless datenum code! hacks make it work for decimals.
                if event_time<43200; event_time = event_time + 86400; end % add a day for data between midnight and midday
                EventStart(end+1,1)=event_time;
                %timecheck = TimeinsecToDatetime(event_time);
                EventManvAuto(end+1,1) = strcmp(Stemp2{j}.nti_colon_manual.Text,'true');
                EventDuration(end+1,1) = str2num(Stemp2{j}.xdf_colon_duration.Text);
                try
                    classstr = [' ' strtrim(Stemp2{j}.xdf_colon_class.Text)];
                catch
                    classstr = '';
                end
                EventName{end+1,1} = [EventClasses{k} classstr];
            end
        end
        
        %for arousals of zero seconds:
        EventDuration(EventDuration<3)=3;
        
        PrelimEventTable = table(EventName,EventStart,EventDuration,EventManvAuto);
        
        Nevents = length(EventName); %needed?
        
        ScorerTable.ScorerName(xx) = ScorerNames(CandidateScorer(xx));
        ScorerTable.Nevents(xx) = Nevents;
        %
        
        %         if ~isfield(S.CMPStudyConfig.ScoredEvents,'ScoredEvent')
        %             displaytext=['Warning: No scored events'];
        %             if ~isempty(handletext), disp(displaytext); set(handletext,'String',displaytext); drawnow; end
        %         end
        
        EventCodesListUnique = unique(EventName);
        
        % normal procedures
        EventLabelOptions= {...
            'microarousal',...
            'apnea central',...
            'apnea mixed',...
            'apnea obstructive',...
            'hypopnea obstructive',... %
            'hypopnea central',... %
            'customevent',... %assumed to be the hypopnea 4% category based on a single example study.
            'rera ' % has space in it
            };
        EventCodeOptions=[1 3 5 2 4 6 4 12];
        EventExact1orStartsWith2 = [2 1 1 1 1 1 1 1];
        
       
        
        displaytext=['Get Hypnogram data: ' fname];
        SStages = ScoringStruct.xdf_colon_sleepstages.xdf_colon_sleepstage;
        
        if ~isempty(handletext), disp(displaytext); set(handletext,'String',displaytext); drawnow; end
        
        %Code assumes that there is one entry per 30-s epoch, not
        %necessarily in order, starting at beginning of the EDF file:
        HypnogramLabel=[];
        HypnogramEpochNumber=[];
        HypnogramInfStage=[];
        HypnogramStgStrength=[];
        for i=1:length(SStages)
            temp = SStages{i};
            HypnogramLabel{end+1,1} = strtrim(temp.xdf_colon_stage.Text);
            HypnogramEpochNumber{end+1,1} = str2double(strtrim(temp.xdf_colon_epochnumber.Text));
            HypnogramInfStage{end+1,1} = strtrim(temp.nti_colon_infstage.Text);
            HypnogramStgStrength{end+1,1} = strtrim(temp.nti_colon_stgstrength.Text);
        end %W=0,N1=1,N2=2,N3=3,R=5
        
        HypnogramTable = table(HypnogramLabel,HypnogramEpochNumber,HypnogramInfStage,HypnogramStgStrength);
        HypnogramTable = sortrows(HypnogramTable,'HypnogramEpochNumber');
        %Convert to Terrill code %W=-1,N1=1,N2=2,N3=3,R=5
        Hypnogram = nan(height(HypnogramTable),1);%
        Hypnogram(strcmp(HypnogramTable.HypnogramLabel,'W'))=4;
        Hypnogram(strcmp(HypnogramTable.HypnogramLabel,'1'))=2;
        Hypnogram(strcmp(HypnogramTable.HypnogramLabel,'2'))=1;
        Hypnogram(strcmp(HypnogramTable.HypnogramLabel,'3'))=0;
        Hypnogram(strcmp(HypnogramTable.HypnogramLabel,'4'))=0;
        Hypnogram(strcmp(HypnogramTable.HypnogramLabel,'R'))=3;
        
        ScorerTable.Nepochs(xx) = sum(~isnan(Hypnogram));
        
            if xx==Nvalid
                disp('Scorer Info--')
                disp(ScorerTable)
                [~,ii]=max(ScorerTable.Nevents);
                CandidateScorer(Nvalid+1)=CandidateScorer(ii);
            end
            
        end
        
        %%
    case 'Danalyzer'
        %Appears to be for staging only at present
        
        directory = char(Filenames(5));
        fname = char(Filenames(2));
        
        %Xml events
        displaytext=['Get Events data:' fname];
        if ~isempty(handletext), disp(displaytext); set(handletext,'String',displaytext); drawnow; end
        
        importscoring=load([directory fname],'stageData');
        
        
        Nevents = 0;
        
        EventName=[];
        EventStart=[];
        EventDuration=[];
        
        % include RERAs as hypopneas
        EventLabelOptions= {...
            'Arousal',...
            'Central Apnea',...
            'Mixed Apnea',...
            'Obstructive Apnea',...
            'Hypopnea',... % unconvincing hypop (for MESA, MrOS)
            'RERA',...
            'Unsure',... % Special case for large cohort studies with >50% reduction in flow (e.g. MESA)
            };
        EventCodeOptions=[1 3 5 2 4 4 4];
        EventExact1orStartsWith2 = [2 1 1 1 1 1 1];
        
        displaytext=['Get Hypnogram data: ' fname];
        HypnogramTemp = importscoring.stageData.stages;
        
        Nepochs = length(HypnogramTemp);
        Hypnogram=NaN*zeros(Nepochs,1);
        
        Hypnogram(HypnogramTemp==0)=4; %W=0,N1=1,N2=2,N3=3,R=5,?=7
        Hypnogram(HypnogramTemp==1)=2; %W=0,N1=1,N2=2,N3=3,R=5,?=7
        Hypnogram(HypnogramTemp==2)=1; %W=0,N1=1,N2=2,N3=3,R=5,?=7
        Hypnogram(HypnogramTemp==3)=0; %W=0,N1=1,N2=2,N3=3,R=5,?=7
        Hypnogram(HypnogramTemp==5)=3; %W=0,N1=1,N2=2,N3=3,R=5,?=7
        Hypnogram(HypnogramTemp==7)=8; %W=0,N1=1,N2=2,N3=3,R=5,?=7
        
        % setting NaN here, then sets 8 below
        
        %end ProfusionXML
        
        %end ProfusionXML
    case 'Michele'
        %%
        directory = char(Filenames(5));
        fname = char(Filenames(2));
        [~,~,eventstxt] = xlsread([directory filesep char(Filenames(2))],1,'A2:D5001'); %max number of events = 5k
        %delete excess
        temp = cell2mat(eventstxt(:,1));
        temp = find(isnan(temp),1);
        eventstxt(temp:end,:)=[];
        EventName = eventstxt(:,2);
        EventStart = cell2mat(eventstxt(:,3)); %warning all data are nearest second
        EventStart=EventStart+StartTime;
        EventDuration = cell2mat(eventstxt(:,4));
        EventLabelOptions= {...
            'A' ,... %1
            'A/L' ,...
            'A/R' ,...
            'O' ,... %4
            'O/A' ,...
            'O/A/D' ,...
            'O/D',...
            'C' ,... %8
            'C/A' ,...
            'C/A/D' ,...
            'C/D',...
            'M' ,... %12
            'M/A' ,...
            'M/A/D' ,...
            'M/D',...
            'H' ,... %16
            'H/A' ,...
            'H/A/D' ,...
            'H/D',...
            'HO' ,... %20
            'HO/A' ,...
            'HO/A/D' ,...
            'HO/D',...
            'RERA/A' ,... %24
            'RERA/A/D'};
        EventCodeOptions=[1 1 1 2 2 2 2 3 3 3 3 5 5 5 5 4 4 4 4 4 4 4 4 12 12];
        EventExact1orStartsWith2 = repmat(1,1,length(EventCodeOptions));
        
        HTable = table(EventName,EventStart,EventDuration);
        HCodes = {'WAKE','Wake','NREM1','Stage 1','NREM2','Stage 2','NREM3','Stage 3','REM'}; % nb.'Unknown Sleep Stage' is faulty (erroneous durations and overlap with actual stages)
        HCodesNum = [4 4 2 2 1 1 0 0 3];
        Iepoch = sum(EventName == string(HCodes),2)>0;
        HTable2 = HTable(Iepoch,:);
        
        HTable2.Codes = repmat(8,height(HTable2),1);
        
        for i=1:length(HCodes)
            HTable2.Codes(HTable2.EventName==string(HCodes(i)))=HCodesNum(i);
        end
        
        
        starttimehyp = (HTable2.EventStart(1) - StartTime);
        starttimediff = starttimehyp - floor(starttimehyp/30)*30;
        
        
        HTable2.EpochNumber = round((HTable2.EventStart - HTable2.EventStart(1))/30) +1 ;
        
        HTable3=table();
        HTable3.EpochNumber = [1:max(HTable2.EpochNumber)]';
        HTable3.Codes = nan(height(HTable3),1);
        HTable3.EventStart = (HTable3.EpochNumber-1)*30 + HTable2.EventStart(1);
        
        HTable3.Codes(HTable2.EpochNumber)=HTable2.Codes;
        HTable3.Codes(isnan(HTable3.Codes))=8;
        if 0 %debug plot
            figure(111); clf(111);
            plot(HTable3.EventStart,HTable3.Codes,'.-');
            hold on
            plot(HTable2.EventStart,HTable2.Codes,'.-');
        end
        %Hypnogram_t=((StartTime+15+starttimediff):30:EndTime)';
        Hypnogram = HTable3.Codes;
        
        %Hypnogram = interp1(HTable2.EventStart+15,HTable2.Codes,Hypnogram_t,'nearest');
        %Hypnogram = [repmat(8,temp,1) ; HTable2.Codes];
        
    case 'ApneaLink'
        
        directory = char(Filenames(5));
        fname = char(Filenames(2));
        
        [Label,Transducer,Fs] = EDFChannelLabels([directory fname]);
        
        
        EventLabelOptions = {'Unclassified apn','Obstructive apne','Hypopnea        ','Central apnea   ','Mixed apnea     '}';
        EventCodeOptions = [2 2 4 3 5];
        
        [Label,Transducer,Fs] = EDFChannelLabels([directory fname]);
        ChannelNumbers = sum((Label == string(EventLabelOptions)) .* repmat(1:length(Label),length(EventLabelOptions),1),2);
        I = ChannelNumbers==0;
        EventLabelOptions(I)=[];
        ChannelNumbers(I)=[];
        
        EventTable=[];
        for i=1:length(EventLabelOptions)
            [SignalTemp,ChannelsFs(i),~,~,LabelCheck{i},~,~,~,~,~] = readedfrev3([directory fname],ChannelNumbers(i)-1,0,Inf);
            I1 = find(diff(SignalTemp)>0);
            I2 = find(diff(SignalTemp)<0);
            [I1,I2] = TidyStartEndEventList(I1,I2,length(SignalTemp));
            if ~isempty(I1)
                EventStartTemp = (I1-1)/ChannelsFs(i) + StartTime;
                EventEndTemp = (I2-1)/ChannelsFs(i) + StartTime;
                EventCodeTemp = EventCodeOptions(i)*ones(length(I1),1);
                EventCodesListTemp = repmat(EventLabelOptions(i),length(I1),1);
                OutTable = table(EventStartTemp,EventEndTemp,EventCodesListTemp,EventCodeTemp);
                EventTable = [EventTable;OutTable];
            end
        end
        EventTable.EventDurationTemp = EventTable.EventEndTemp - EventTable.EventStartTemp;
        
        EventName=EventTable.EventCodesListTemp;
        EventStart=EventTable.EventStartTemp;
        EventDuration=EventTable.EventDurationTemp;
        
        EventExact1orStartsWith2 = EventCodeOptions*0 + 1;
        
        displaytext=['Get Hypnogram data: ' fname];
        SleepSignals = {'Start of evaluat','End of evaluatio'}';
        
        temp = ((Label == string(SleepSignals)) .* repmat(1:length(Label),length(SleepSignals),1))';
        temp(temp==0)=NaN;
        ChannelNumbers = min(temp)';
        I = ChannelNumbers==0 | isnan(ChannelNumbers);
        SleepSignals(I)=[];
        ChannelNumbers(I)=[];
        
        clear ChannelsFs LabelCheck
        i=1
        [SignalTemp,ChannelsFs(i),~,~,LabelCheck{i},~,~,~,~,~] = readedfrev3([directory fname],ChannelNumbers(i)-1,0,Inf);
        I1 = find(diff(SignalTemp)>0);
        I2 = find(diff(SignalTemp)<0);
        [I1,I2] = TidyStartEndEventList(I1,I2,length(SignalTemp));
        EndOfStartT = (I2-1)/ChannelsFs(i);
        
        i=2
        [SignalTemp,ChannelsFs(i),~,~,LabelCheck{i},~,~,~,~,~] = readedfrev3([directory fname],ChannelNumbers(i)-1,0,Inf);
        I1 = find(diff(SignalTemp)>0);
        I2 = find(diff(SignalTemp)<0);
        [I1,I2] = TidyStartEndEventList(I1,I2,length(SignalTemp));
        StartofEndT = (I1-1)/ChannelsFs(i);
        
        Nepochs = ceil((EndTime-StartTime)/30)
        Hypnogram = 4 + zeros(Nepochs,1);
        HypnogramT = 15 + 30*(0:(Nepochs-1))';
        
        Hypnogram(HypnogramT>EndOfStartT & HypnogramT<StartofEndT)=2; %code for stage 1 sleep
        
        
         case 'NSRR' %%
        %%
        %%
        directory = char(Filenames{1,5});
        fname = char(Filenames{1,2});
        
        displaytext=['Importing events from:' fname];disp(displaytext);
        if ~isempty(handletext), disp(displaytext); set(handletext,'String',displaytext); drawnow; end
        
        clear eventstemp
        fid = fopen([directory filesep fname]);
        i=1;
        while 1
            eventstemp{i,:} = fgetl(fid);  %read in the data
            if eventstemp{i,:}==-1
                eventstemp(i,:) = [];
                break
            end
            i=i+1;
        end
        fclose(fid);   %close the file
        
        HeaderTextList = {'class','instance'};
        I=[];
        for ii=1:length(HeaderTextList)
            I=[I,contains(eventstemp,HeaderTextList{ii})];
        end
        I=sum(I,2)/length(HeaderTextList);
        NPreHeaderlines=find(I==1)-1;
        if ~isempty(NPreHeaderlines)
            HeaderTextExact = eventstemp{I==1};
            VarNameYN=1;
        else
            HeaderTextExact = 1;
            NPreHeaderlines=0;
            VarNameYN=0;
        end
        
        % Import Table (now that we know where the data starts)
        EpochEventT = readtable([directory fname],'Delimiter','\t','ReadVariableNames',VarNameYN,'HeaderLines',NPreHeaderlines,'FileType','text');
        
        numCols=size(EpochEventT,2);
        
        if numCols>1 % proper nsrr harmonized file
            
            clear EventStartTime
            I=strcmp(EpochEventT.instance,'.');
            EpochEventT.instance(I)={NaN};
            EpochEventT.Event=strcat(EpochEventT.class,' ',EpochEventT.instance);
            
            for jj=1:size(EpochEventT.Event,1)
                try
                    starttimetemp=mod(datenum(char(EpochEventT.start(jj)),'HH:MM:SS.FFF'),1)*86400;
                catch
                    starttimetemp=mod(datenum(char(EpochEventT.start(jj)),'HH:MM:SS'),1)*86400;
                end
                if starttimetemp<43200; starttimetemp=starttimetemp+86400; end
                try
                    endtimetemp=mod(datenum(char(EpochEventT.stop(jj)),'HH:MM:SS.FFF'),1)*86400;
                catch
                    endtimetemp=mod(datenum(char(EpochEventT.stop(jj)),'HH:MM:SS'),1)*86400;
                end
                if endtimetemp<43200; endtimetemp=endtimetemp+86400; end
                EpochEventT.StartTime(jj)=starttimetemp-StartTime;
                EpochEventT.EndTime(jj)=endtimetemp-StartTime;
                EpochEventT.Duration(jj)=EpochEventT.EndTime(jj)-EpochEventT.StartTime(jj);
                EpochEventT.EventStartTime(jj)=starttimetemp;
                EpochEventT.EventEndTime(jj)=endtimetemp;
            end
            
            %         EpochEventT.EventStartTime=EpochEventT.start+StartTime;
            %         EpochEventT.Duration=EpochEventT.stop-EpochEventT.start;
            
            StageLabelOptions = {'W','R','N1','N2','N3'}';
            StageCodeOptions = [4 3 2 1 -1] ; %can not be zeros here
            
            clear match1 match2
            EventName = EpochEventT.class;
            for m=1:length(StageLabelOptions)% each column is an event type
                match1(:,m)=strcmpi(string(EventName'),string(StageLabelOptions(m)));
                match2(:,m)=startsWith(string(EventName'),string(StageLabelOptions(m)));
            end
            match = 1*( match1>0 | match2>0);
            temp = match.*StageCodeOptions;
            temp(temp==0)=NaN;
            temp(temp==-1)=0;
            EpochEventT.EpochCodes=min(temp')'; %just in case two criteria are detected, will leave the lowest code value in list
            
            HypnogramT = EpochEventT(~isnan(EpochEventT.EpochCodes),:);
            HypnogramT.DurationDiff = [diff(HypnogramT.EventStartTime);NaN];
            if ~isnan(HypnogramT.Duration(end))
                HypnogramT.DurationDiff(end) = HypnogramT.Duration(end);
            else
                HypnogramT.DurationDiff(end) = 30;
            end
            HypnogramT.EventEndTime = HypnogramT.EventStartTime+HypnogramT.DurationDiff;
            
            Hypnogram_t = (HypnogramT.EventStartTime(1)-30:30:HypnogramT.EventEndTime(end))';
            tempH = [8;HypnogramT.EpochCodes;8];
            tempT = [HypnogramT.EventStartTime(1)-30;HypnogramT.EventStartTime;HypnogramT.EventEndTime(end)];
            Hypnogram = interp1(tempT,tempH,Hypnogram_t+15,'previous','extrap');
            Hypnogram(1)=[];
            Hypnogram_t(1)=[];
            
            starttimediff = Hypnogram_t(1)-StartTime;
            
            
            % respiratory events  arousals
            EventLabelOptions= {...
                'arousal',...
                'apneacentral',...
                'apneamixed',...
                'apneaobstructive',...
                'hypopnea',...
                };
            EventCodeOptions=[1 3 5 2 4];
            EventExact1orStartsWith2 = [2 1 2 2 2];
            EventList=EpochEventT.Event;
            Idx=zeros(length(EventList),1);
            for ii=1:length(EventLabelOptions)
                StageTF = startsWith(EventList,EventLabelOptions{ii},'IgnoreCase',true);
                Idx=Idx|StageTF;
            end
            EpochEventT.RespEvent=Idx; % just to double check
            EventName=EpochEventT.Event(EpochEventT.RespEvent==1);
            EventStart=EpochEventT.EventStartTime(EpochEventT.RespEvent==1);
            EventDuration=EpochEventT.Duration(EpochEventT.RespEvent==1);
            
            % end of nsrr events import
        else
            warning ("ASSUMING ANNOTATION HAS ONLY STAGING INFO. IF NOT PLEASE RECHECK ANNOTATION FILE")
            Nepochs = length(EpochEventT.Var1);
            Hypnogram=NaN*zeros(Nepochs,1);
            
            Hypnogram(strcmp(EpochEventT.Var1,'W')|strcmp(EpochEventT.Var1, 'wake'))=4;
            Hypnogram(strcmp(EpochEventT.Var1, 'NREM1') | strcmp(EpochEventT.Var1, 'N1'))=2;
            Hypnogram(strcmp(EpochEventT.Var1, 'NREM2')| strcmp(EpochEventT.Var1, 'N2'))=1;
            Hypnogram(strcmp(EpochEventT.Var1, 'NREM3')| strcmp(EpochEventT.Var1, 'N3'))=0;
            Hypnogram(strcmp(EpochEventT.Var1, 'NREM4')| strcmp(EpochEventT.Var1, 'N4'))=0;
            Hypnogram(strcmp(EpochEventT.Var1, 'REM')| strcmp(EpochEventT.Var1, 'R'))=3;
            Hypnogram(strcmp(EpochEventT.Var1, '?') | strcmp(EpochEventT.Var1, 'L'))=8;
            
        end
             
    case 'AnnotEannot'  % WSC
        
        directory = char(Filenames(5));
        fname = char(Filenames(3)); % .eannot for sleep stages
        SleepStagetemp = readtable(fname,'FileType','text','ReadVariableNames',false);
        Nepochs = length(SleepStagetemp.Var1);
        Hypnogram=NaN*zeros(Nepochs,1);
        
        Hypnogram(strcmp(SleepStagetemp.Var1,'W')|strcmp(SleepStagetemp.Var1, 'wake'))=4;
        Hypnogram(strcmp(SleepStagetemp.Var1, 'NREM1') | strcmp(SleepStagetemp.Var1, 'N1'))=2;
        Hypnogram(strcmp(SleepStagetemp.Var1, 'NREM2')| strcmp(SleepStagetemp.Var1, 'N2'))=1;
        Hypnogram(strcmp(SleepStagetemp.Var1, 'NREM3')| strcmp(SleepStagetemp.Var1, 'N3'))=0;
        Hypnogram(strcmp(SleepStagetemp.Var1, 'NREM4')| strcmp(SleepStagetemp.Var1, 'N4'))=0;
        Hypnogram(strcmp(SleepStagetemp.Var1, 'REM')| strcmp(SleepStagetemp.Var1, 'R'))=3;
        Hypnogram(strcmp(SleepStagetemp.Var1, '?') | strcmp(SleepStagetemp.Var1, 'L'))=8;
        
        fname = char(Filenames(2)); % .annot for resp. events/arousals
        EpochEventT = readtable(fname,'FileType','text','ReadVariableNames',false);
        EpochEventT.Properties.VariableNames={'Event','Misc1','Misc2','StartTime1','EndTime1','Comments'};
        for jj=1:size(EpochEventT.Event,1)
            starttimetemp=mod(datenum(char(EpochEventT.StartTime1(jj)),'HH:MM:SS'),1)*86400;
            if starttimetemp<43200; starttimetemp=starttimetemp+86400; end
            endtimetemp=mod(datenum(char(EpochEventT.EndTime1(jj)),'HH:MM:SS.FFF'),1)*86400;
            if endtimetemp<43200; endtimetemp=endtimetemp+86400; end
            EpochEventT.StartTime(jj)=starttimetemp-StartTime;
            EpochEventT.EndTime(jj)=endtimetemp-StartTime;
            EpochEventT.Duration(jj)=EpochEventT.EndTime(jj)-EpochEventT.StartTime(jj);
            EpochEventT.EventStartTime(jj)=starttimetemp;
            EpochEventT.EventEndTime(jj)=endtimetemp;
        end
        
        
        EventLabelOptions= {...
            'arousal_spontaneous',...
            'arousal_respiratory',...
            'arousal_lm',...
            'arousal',...
            'apnea_central',...
            'apnea_mixed',...
            'apnea_obstructive',...
            'hypopnea',...
            };
        EventCodeOptions=[1 1 1 1 3 5 2 4];
        EventExact1orStartsWith2 = [2 2 2 2 1 2 2 2];
        EventList=EpochEventT.Event;
        Idx=zeros(length(EventList),1);
        for ii=1:length(EventLabelOptions)
            StageTF = startsWith(EventList,EventLabelOptions{ii},'IgnoreCase',true);
            Idx=Idx|StageTF;
        end
        EpochEventT.RespEvent=Idx; % just to double check
        EventName=EpochEventT.Event(EpochEventT.RespEvent==1);
        EventStart=EpochEventT.EventStartTime(EpochEventT.RespEvent==1);
        EventDuration=EpochEventT.Duration(EpochEventT.RespEvent==1);
       
        
end % end of system specific events import

%% General Code: Make Events, Arousals, Epochs signals

%%
% Moved To Import Settings
% if isfield(settings, 'IncludeRERAs') && settings.IncludeRERAs==1
%     disp('Warning: IncludeRERAs is now obselete code, replace with ''settings.ReallocateEventCodes4 = 12'' ');
%     %Re-allocate event codes (12-RERA) to 4-obstructive hypopneas
%     if isfield(settings, 'ReallocateEventCodes4') && ~ismember(settings.ReallocateEventCodes4, 12)
%         settings.ReallocateEventCodes4 = [settings.ReallocateEventCodes4, 12];
%     end
% end

%Re-allocate chosen events as obstructive hypopneas
if isfield(settings,'ReallocateEventCodes4')
    for i=1:length(settings.ReallocateEventCodes4)
        temp = find(EventCodeOptions==settings.ReallocateEventCodes4(i));
        try disp(['Scored ' EventLabelOptions{temp} ' taken to be obstructive hypopneas [4]']); end
        EventCodeOptions(temp)=4;
    end
end


if strcmp(system, 'SpikeDise')
    EventStart=EventStart(:);
    EventDuration=EventDuration(:);
    if ~isfield(Evts,'source') % quick fix for data with no events data
        EventCodes = nan(0,0);
    else
        EventCodes = Evts.source;
    end
    
    clear Evts
    Evts.EventName = EventName;
    Evts.EventLabelOptions = EventLabelOptions;
    Evts.Table1 = table(EventName,EventStart,EventDuration);
    Evts.Table1.EventCodes = EventCodes;
    
    % might rename these depending on how cumbersome it will be later on
    Evts.RespT = Evts.Table1(Evts.Table1.EventCodes==1,:); % RespT = VOTE classifications
    Evts.SideSleep = Evts.Table1(Evts.Table1.EventCodes==3,:); % ArT = JawThrust
    Evts.MouthOpenClosed = Evts.Table1(Evts.Table1.EventCodes==2,:); % Site out of view marker
    
    Evts.ArT = Evts.Table1(Evts.Table1.EventCodes==2,:); % DUMMY VARIABLE
    Evts.Hypnogram = NaN;
    Evts.Hypnogram_t = NaN;
else
    EventStart=EventStart(:);
    EventDuration=EventDuration(:);
    EventName=EventName(:);
    
    clear Evts
    Evts.Table1 = table(EventName,EventStart,EventDuration);
    Evts.Table1 = sortrows(Evts.Table1,'EventStart');
    EventName = Evts.Table1.EventName;
    EventStart = Evts.Table1.EventStart;
    EventDuration = Evts.Table1.EventDuration;
    
    Evts.EventName = EventName;
    Evts.EventLabelOptions = EventLabelOptions;
    
    if ~isempty(EventName)
        clear match
        for m=1:length(EventLabelOptions)% each column is an event type
            if EventExact1orStartsWith2(m)==1 %exact
                match(:,m)=string(EventLabelOptions(m))==string(EventName');
            elseif EventExact1orStartsWith2(m)==2 %startswith
                match(:,m)=startsWith(string(EventName'),string(EventLabelOptions(m)));
            end
        end
        temp = match.*EventCodeOptions;
        temp(temp==0)=NaN;
        EventCodes=min(temp'); %just in case two criteria are detected, will leave the lowest code value in list
        
    else
        EventCodes=[];
        disp(['Events Empty']);
    end
    Evts.Table1.EventCodes=EventCodes(:);
    
    %check for events that should be in the list
    if ~isempty(EventName(isnan(EventCodes)))
        temp = string(EventName(isnan(EventCodes))); %converted to string as cell array was having irregular trouble.
        disp(['Events found but unused: ' char(strjoin(unique(temp),', '))]);
    else
        disp(['Events OK']);
    end
    
    Evts.RespT = Evts.Table1(Evts.Table1.EventCodes>1,:);
    Evts.ArT = Evts.Table1(Evts.Table1.EventCodes==1,:);
    
    %     I=isnan(Evts.Table1.EventCodes);
    %     Evts.RespT(I,:)=[];
    %     EventCodes(I)=[];
    %     EventStart(I)=[];
    %     EventDuration(I)=[];
    %
    %remove very short resp events
    I=isnan(Evts.RespT.EventDuration<2 & Evts.RespT.EventCodes(:)>1);
    Evts.RespT(I,:)=[];
    %
    %     I=EventDuration(:)<2&EventCodes(:)>1;
    %     EventCodes(I)=[];
    %     EventStart(I)=[];
    %     EventDuration(I)=[];
    
    if sum(Evts.ArT.EventCodes==1)==0
        displaytext=['Warning: No scored arousals'];
        disp(['Warning: No scored arousals']);
        if ~isempty(handletext), disp(displaytext); set(handletext,'String',displaytext); drawnow; end
        
    end
    
    
    if std(Evts.ArT.EventDuration)<1
        displaytext=['Warning: Arousals are all similar durations'];
        disp(['Warning: Arousals are all similar durations: ' num2str(mode(Evts.ArT.EventDuration)) ' s']);
        if ~isempty(handletext), disp(displaytext); set(handletext,'String',displaytext); drawnow; end
    end
    
    % if "EdfUnscored" is set, then processing gets to here then dies 
    % (i.e. drops out to the catch in the calling function)
    % it dies because there is no variable called 'Hypnogram'.
    % question is, do we need to exit cleanly with empty Evts, 
    % or is drop out an ok option? when it drops out, it checks for Evts, 
    % and if not exist, it makes empty evts as follows; Evts=struct();
    if 1 %set NaN, others, to 8
        Hypnogram(isnan(Hypnogram))=8; %no stage code
        Hypnogram(Hypnogram<0)=8; %no stage code
    end
    if ~exist('Hypnogram_t')
        Hypnogram_t = StartTime+starttimediff+(0:30:(30*(length(Hypnogram)-1)))';
    end
    
    Evts.Hypnogram = Hypnogram;
    Evts.Hypnogram_t = Hypnogram_t;
    
    if sum(Evts.Hypnogram==4)==length(Evts.Hypnogram)
        displaytext=['Warning: Entire hypnogram is wake'];
        if ~isempty(handletext), disp(displaytext); set(handletext,'String',displaytext); drawnow; end
    end
    
    if exist('PositionT')
        Evts.PositionT = PositionT; %Used in NKxdf, but can be written for other systems with pos lists.
    end
end