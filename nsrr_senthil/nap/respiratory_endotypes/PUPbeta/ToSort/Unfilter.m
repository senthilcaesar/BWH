%FilterUnfilter

Flow1 = Flow.values+0.005391;

Flow1(Flow1>0)=1.5*Flow1(Flow1>0);
dt = Time(2)-Time(1);

filter_HFcutoff_butter0 = 1/10; %4Hz=250ms
filter_order0 = 1;
[B_butter0,A_butter0] = butter(filter_order0,filter_HFcutoff_butter0/(1/dt/2),'high');
Flowf = filter(B_butter0,A_butter0,Flow1);


figure(1); clf(1)
plot(Time,Flow1);
hold('on');
plot(Time,Flowf);
noise = 0.01*randn(length(Time),1);
Flowuf = filter(A_butter0,B_butter0,Flowf+noise);
plot(Time,Flowuf);

%% Try with real data 
%fname = 'C:\Users\szs88\Dropbox (Partners HealthCare)\WtcOSATraits\061016dc-d.edf'
fname = 'C:\Users\szs88\Dropbox (Partners HealthCare)\WtcOSATraits\080220cd-d.edf'
fid = fopen(fname,'r');
fseek(fid,252,-1);
M  = sscanf(char(fread(fid,4,'char')'),'%d');
for m=1:M+12
    fseek(fid,(m-1)*16+256,-1);
    Label{m} = char(fread(fid,16,'char')'); %%%%%%%%%
end
fclose(fid); % Close file

ch=14;
eval(['[Signal,Fs,~,~,~,~,~,~,~] = readedfrev3(fname,ch-1,0,Inf);']); %the '-1' is because EDF channel numbers start from 0.

dt=1/Fs;
N=length(Signal);
Time = (0:dt:dt*(N-1))';

%%
figure(1); clf(1);
plot(Time,Signal);
hold('on');
plot(Time,0*Signal);

filter_HFcutoff_butter0 = 1/30; %4Hz=250ms
filter_order0 = 1;
[B_butter0,A_butter0] = butter(filter_order0,filter_HFcutoff_butter0/(1/dt/2),'high');
Flowuf = filter(A_butter0,B_butter0,Signal);
h=plot(Time,detrend(Flowuf));



