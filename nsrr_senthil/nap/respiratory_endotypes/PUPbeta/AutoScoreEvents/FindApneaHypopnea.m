function [ApneaHypopneamat2]=FindApneaHypopnea(VE_peupnea,SigT,ChannelsFs,ChannelsList,minEventDuration)
global settings

    if ~exist('minEventDuration')==1 || isempty('minEventDuration')==1
        minEventDuration=10;
    end

    disp(['minEventDuration=' num2str(minEventDuration)]);
%% Find events
ApneaHypopneatemp = (VE_peupnea<70)*1;

%% merge closely-spaced dropouts, <2s
removeshorterthani = 2/(1/ChannelsFs);
ApneaHypopnea_=ApneaHypopneatemp;

I = diff([NaN;1-ApneaHypopnea_]);
I1 = find(I==1);
I2 = find(I==-1);
[I1N,I2N]=TidyStartEndEventList(I1,I2,length(I));
lengthsN = I2N-I1N;
removeX = lengthsN<removeshorterthani;
for j=1:length(lengthsN)
    if ~removeX(j)
        continue
    end
    ApneaHypopnea_(I1N(j):I2N(j))=1;
end

%% remove short events, <10s
removeshorterthani = minEventDuration/(1/ChannelsFs);
clear I I1 I2 I1N I2N
I = diff([NaN;ApneaHypopnea_]);
I1 = find(I==1);
I2 = find(I==-1);
[I1N,I2N]=TidyStartEndEventList(I1,I2,length(I));
lengthsN = I2N-I1N;
removeX = lengthsN<removeshorterthani;
for j=1:length(lengthsN)
    if ~removeX(j)
        continue
    end
    ApneaHypopnea_(I1N(j):I2N(j))=0;
end

%Repeat for apnea
Apnea = (VE_peupnea<10)*1;

removeshorterthani = 2/(1/ChannelsFs);
Apnea_=Apnea;

clear I I1 I2 I1N I2N
I = diff([NaN;1-Apnea_]);
I1 = find(I==1);
I2 = find(I==-1);
[I1N,I2N]=TidyStartEndEventList(I1,I2,length(I));
lengthsN = I2N-I1N;
removeX = lengthsN<removeshorterthani;
for j=1:length(lengthsN)
    if ~removeX(j)
        continue
    end
    Apnea_(I1N(j):I2N(j))=1;
end

removeshorterthani = minEventDuration/(1/ChannelsFs);
clear I I1 I2 I1N I2N
I = diff([NaN;Apnea_]);
I1 = find(I==1);
I2 = find(I==-1);
[I1N,I2N]=TidyStartEndEventList(I1,I2,length(I));
lengthsN = I2N-I1N;
removeX = lengthsN<removeshorterthani;
for j=1:length(lengthsN)
    if ~removeX(j)
        continue
    end
    Apnea_(I1N(j):I2N(j))=0;
end

ApneaHypopnea_(Apnea_==1)=2; % apnea=2; hypopnea=1

%% change hypopneas ending in apneas as apnea
clear I I1 I2 I1N I2N
I = diff([NaN;1*(ApneaHypopnea_>0)]);
I1 = find(I==1);
I2 = find(I==-1);
[I1N,I2N]=TidyStartEndEventList(I1,I2,length(I));
evntdur=(I2N-I1N)*(1/ChannelsFs);
ApneaHypopneamat2=ApneaHypopnea_;

for i=1:length(I1N)
    Ir = I1N(i):(I2N(i)-1);
    temp = ApneaHypopnea_(Ir);
    ApneaHypopneamat2(Ir)=max(temp);
    EventCodes(i)=max(temp);
end
figure (22);clf (figure (22)); plot(SigT.Time,ApneaHypopnea_);
hold on
plot(SigT.Time, 1*(ApneaHypopneamat2));

%% events array
EventStart=SigT{I1N,1};
EventDuration=(I2N-I1N)*(1/ChannelsFs);
EventCodes=EventCodes(:); % apnea=2; hypop=1

%% remove apneas in wake
Stages_=SigT.Epochs;
DistFromSleep=howfarawayfromsleepInf(1*(Stages_>3),ones(size(Stages_)))*(1/ChannelsFs(1));

figure(33); clf(figure(33));
ax(1)=subplot(3,1,1);
plot(SigT.Time,Stages_);
hold on;
plot(SigT.Time,ApneaHypopneamat2,'k'); % orginal apnea-hypopnea events
hold on;

temp=[]; SleepEvts=ones(length(I1N),1);
for ii=1:length(I1N)
    Ir = I1N(ii):(I2N(ii)-1);
    temp(ii) = max(DistFromSleep(Ir)); % changed to max; was min earlier
    if temp(ii)>60 % changed to >30; earlier it was <120.
    ApneaHypopneamat2(Ir)=0;
    SleepEvts(ii)=0;
    end
end
%% update events array
EventStart(SleepEvts==0)=[];
EventDuration(SleepEvts==0)=[];
EventCodes(SleepEvts==0)=[]; % apnea=2; hypop=1
if isfield(settings,'AutoScoredEventsAllCentral') && settings.AutoScoredEventsAllCentral==1
    ApneaHypopneamat2(ApneaHypopneamat2==1)=6; % hypopnea
    ApneaHypopneamat2(ApneaHypopneamat2==2)=3; % apnea
    EventCodes(EventCodes==1)=6;
    EventCodes(EventCodes==2)=3;
else
    ApneaHypopneamat2(ApneaHypopneamat2==1)=4; % hypopnea
    ApneaHypopneamat2(ApneaHypopneamat2==2)=2; % apnea
    EventCodes(EventCodes==1)=4;
    EventCodes(EventCodes==2)=2;
end

plot(SigT.Time,ApneaHypopneamat2,'Color',[1 0.5 0.6]); % events in sleep only
hold on;
legend('stage','all events','events in sleep');
ax(2)=subplot(3,1,2);
plot(SigT.Time,DistFromSleep);
legend('distance from sleep');

ax(3)=subplot(3,1,3);
yyaxis right
plot(SigT.Time,ApneaHypopneamat2,'Color',[1 0.5 0.6]); % events in sleep only
yyaxis left
plot(SigT.Time,SigT.Flow);
legend('flow','events in sleep');
xlabel('time(s)');
linkaxes(ax,'x');



