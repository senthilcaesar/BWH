clear all
close all

%[y,fs]=audioread('G:\Partners Healthcare Dropbox\SATP Group\20211202-1538_Recording_1.wav');

%[y,fs]=audioread('G:\Partners Healthcare Dropbox\SATP Group\20211202_163042.wav'); %no cal sig
%[y,fs]=audioread('G:\Partners Healthcare Dropbox\SATP Group\20211202_165337.wav'); biggerbydB = 94-94; %(clipping at 105.6) -dig gain = 1?
%[y,fs]=audioread('G:\Partners Healthcare Dropbox\SATP Group\20211202_165442.wav'); biggerbydB = 94-87.55; %(clipping at 111.1) -dig gain = 0.5 + foam?
%[y,fs]=audioread('G:\Partners Healthcare Dropbox\SATP Group\20211202_172127.wav'); biggerbydB = 94-100.4; %(clipping at 105.6) -dig gain = 0.5?
%[y,fs]=audioread('G:\Partners Healthcare Dropbox\SATP Group\20211202_173054.wav'); biggerbydB = 94-67.74; %(clipping at 123.2) -dig gain = 0.5 + foam?
%[y,fs]=audioread('G:\Partners Healthcare Dropbox\SATP Group\20211203_161653.wav'); biggerbydB = 94-71.66; %(clipping at 128) -dig gain = 1 + optimal foam?
%
%[y,fs]=audioreadforce('G:\Partners Healthcare Dropbox\SATP Group\20211202_173054.wav');
%
%[y,fs]=audioread('G:\Partners Healthcare Dropbox\SATP Group\20211203_161653.wav');

[Snore,fs]=audioread('G:\Partners Healthcare Dropbox\SATP Group\SnoreAtHome\1\before cal.wav'); biggerbydB = 114-92.02; %(clipping at 105.6) -dig gain = 1?
[Snore,fs]=audioread('G:\Partners Healthcare Dropbox\SATP Group\SnoreAtHome\1\after cal.wav'); biggerbydB = 114-95.56; %(clipping at 105.6) -dig gain = 1?
[Snore,fs]=audioread('G:\Partners Healthcare Dropbox\SATP Group\SnoreAtHome\1\dpw2192t1.wav'); biggerbydB = 114-95.56; %(clipping at 105.6) -dig gain = 1?


%%
% load 'G:\Partners Healthcare Dropbox\SATP Group\test.mat'
% y = Snore.values; %0.23921177
% fs=10000;
Time = [0:(1/fs):(length(Snore)-1)*(1/fs)]';
if 1
 Fcut = 1/(2*pi()*0.008);
 filter_order0 = 1;
 [B_butter0,A_butter0] = butter(filter_order0,[Fcut]/(fs/2),'high');
    Snore = filter(B_butter0,A_butter0,Snore);
else %bypass initial filtering

end

temp = 1; %will get overwritten 

if ~exist('biggerbydB')
    biggerbydB=94-94; %91.65 is ideal arbitrary size
end
temp = 10.^(biggerbydB/20);

Snore=Snore*10.055*temp; %calibrated Snore;

SnrPower = Snore.^2; %calibrated mean squared i.e. RMS^2 ~ sound power

Fcut = 1/(2*pi()*0.05);
 filter_order0 = 1;
 [B_butter0,A_butter0] = butter(filter_order0,[Fcut]/(fs/2),'low');
SnrPower = filter(B_butter0,A_butter0,SnrPower);
if 1
    clear yPower
end
SnrPower = downsample(SnrPower,25); %11025->441 Hz
fs2 = fs/25;
TimedB = [0:(1/fs2):(length(SnrPower)-1)*(1/fs2)]';
SnoreDB = 10*(log(SnrPower/(0.00002.^2))/log(10));

figure(1);
ax(1)=subplot(3,1,1); plot(Time,Snore);
ax(2)=subplot(3,1,2); plot(TimedB,SnrPower);
ax(3)=subplot(3,1,3); plot(TimedB,SnoreDB);
linkaxes(ax,'x')
%%
%clear Time TimedB
%save dpw2192t1 '-v7.3'
%save dpw2192t1_SnoreDB SnoreDB '-v7.3'