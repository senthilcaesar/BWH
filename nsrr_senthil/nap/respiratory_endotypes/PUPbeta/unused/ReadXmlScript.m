%readxml

filenamexml='F:\Read Xml Matlab\patient A'

S=xml2struct(filenamexml)

Nepochs = length(S.CMPStudyConfig.SleepStages.SleepStage);
clear Epochs
Epochs=NaN*zeros(Nepochs,1);
for i=1:Nepochs
Epochs(i) = str2num(S.CMPStudyConfig.SleepStages.SleepStage{1,i}.Text);
end %W=0,N1=1,N2=2,N3=3,R=5
figure(1); plot(Epochs);

Nevents = length(S.CMPStudyConfig.ScoredEvents.ScoredEvent);
clear Events
for i=1:Nevents
Events{i}.name = S.CMPStudyConfig.ScoredEvents.ScoredEvent{i}.Name.Text;
Events{i}.start = str2num(S.CMPStudyConfig.ScoredEvents.ScoredEvent{i}.Start.Text);
Events{i}.duration = str2num(S.CMPStudyConfig.ScoredEvents.ScoredEvent{i}.Duration.Text);
end %W=0,N1=1,N2=2,N3=3,R=5


Eventnames = {'Arousal (ARO RES)','Obstructive Apnea','Mixed Apnea','Hypopnea','Central Apnea'}

%find unique names
if 0
UniqueNames=[];
for i=1:Nevents
    unique=1;
    for n=1:length(UniqueNames)
        if strcmp(Events{i}.name.Text,UniqueNames{n})
            unique=0;
            break
        end
    end
    if unique
        UniqueNames{length(UniqueNames)+1}=Events{i}.name.Text
    end
end
end


Eventnames = {'Arousal (ARO RES)','Obstructive Apnea','Mixed Apnea','Hypopnea','Central Apnea'}


