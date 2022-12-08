function [DataEventHypnog_Mat,ChannelsList,ChannelsFs] = addWSIntensity(DataEventHypnog_Mat,ChannelsList,ChannelsFs,noscoredarinwake,mdlA,mdlNoise,RefTable)


logitinverse = @(p) 1./(1+exp(-p));
dt = DataEventHypnog_Mat(2,1) - DataEventHypnog_Mat(1,1);
Fs = 1./dt;
%%
dT = dt;

%%
dsfactor = dT/dt;
DataEventHypnog_Mat_ds = DataEventHypnog_Mat(1:round(dsfactor):end,:);
DataEventHypnog_Mat_ds = array2table(DataEventHypnog_Mat_ds);
DataEventHypnog_Mat_ds.Properties.VariableNames = ChannelsList;

dT = DataEventHypnog_Mat_ds.Time(2)-DataEventHypnog_Mat_ds.Time(1);

DataEventHypnog_Mat_ds.Epochs(DataEventHypnog_Mat_ds.Epochs<0|DataEventHypnog_Mat_ds.Epochs>4)=NaN;

if ~ismember(['Pbetalogfilt1'],DataEventHypnog_Mat_ds.Properties.VariableNames)
    disp('missing Pbetalogfilt1');
end

try
    SpO2off = 1*(DataEventHypnog_Mat_ds.SpO2<50);
    SpO2off = RemoveShortSegments(SpO2off,60,dT);
    
    SpO2starti = find(DataEventHypnog_Mat_ds.SpO2>50,1,'first');
    SpO2endi = find(DataEventHypnog_Mat_ds.SpO2>50,1,'last');
    SpO2ever = 1*SpO2off; SpO2ever(SpO2starti:SpO2endi)=0; SpO2ever=logical(SpO2ever);
catch me
    disp('missing SpO2 signal, skipping');
end

if sum(SpO2ever==1)>0.90*length(SpO2ever)
    disp('missing SpO2 signal, skipping');
end


%%
%DataEventHypnog_Mat_ds

Exclude = SpO2off==1 | isnan(DataEventHypnog_Mat_ds.Epochs) | DataEventHypnog_Mat_ds.Epochs<0 | DataEventHypnog_Mat_ds.Epochs>4;

%for judging methods
ExcludeREM = SpO2off | DataEventHypnog_Mat_ds.EventsAr==1 | DataEventHypnog_Mat_ds.Epochs<0 | DataEventHypnog_Mat_ds.Epochs>3 | isnan(DataEventHypnog_Mat_ds.Epochs);

REMvsNREM = (DataEventHypnog_Mat_ds.Epochs==3)*1;
REMvsNREM(ExcludeREM)=NaN;

%%

if ~noscoredarinwake
    ExcludeAR = Exclude | DataEventHypnog_Mat_ds.Epochs==4 & DataEventHypnog_Mat_ds.EventsAr==0 | DataEventHypnog_Mat_ds.Epochs~=4 & DataEventHypnog_Mat_ds.EventsAr==1;
else
    ExcludeAR = Exclude | DataEventHypnog_Mat_ds.Epochs~=4 & DataEventHypnog_Mat_ds.EventsAr==1;
end

dsfactor = round(3/dT);
I = (1:dsfactor:length(Exclude))';
clear WSpredlogit acc
for j=1:8
    try
        [~,acc(j)] = WSfromoneEEG(Exclude(I),mdlA,mdlNoise,RefTable,DataEventHypnog_Mat_ds(I,:),j,3,ExcludeAR(I));
    catch me
        
    end
end
[~,accmaxi] = max(acc);

[WSintensity,Acc,Exclude_,NoiseBinary] = WSfromoneEEG(Exclude,mdlA,mdlNoise,RefTable,DataEventHypnog_Mat_ds,accmaxi,dT,ExcludeAR);
%Exclude is updated here to include noise, caused trouble e.g.
%for n=53 RICCADSA

PrWake = logitinverse(WSintensity);
% 
% if sum(strcmp(ChannelsList,'WSintensity'))==0
%     DataEventHypnog_Mat = [DataEventHypnog_Mat WSintensity];
%     ChannelsList = [ChannelsList {'WSintensity'}];
%     ChannelsFs = [ChannelsFs;ChannelsFs(1)];
% else
%     idx = find(strcmp(ChannelsList,'WSintensity')==1);
%     DataEventHypnog_Mat(:,idx) = WSintensity;
% end

newsig = {'PrWake','WSintensity'};
for i=1:length(newsig)
if sum(strcmp(ChannelsList,newsig{i}))==0
    DataEventHypnog_Mat = [DataEventHypnog_Mat eval(newsig{i})];
    ChannelsList = [ChannelsList {newsig{i}}];
    ChannelsFs = [ChannelsFs;ChannelsFs(1)];
else
    idx = find(strcmp(ChannelsList,newsig{i})==1);
    DataEventHypnog_Mat(:,idx) = eval(newsig{i});
end
end


