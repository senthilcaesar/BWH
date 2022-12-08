function[daqalltemp]=LabviewReadAllChannels(settings,n,Filenames,SynchSM)

LabViewFilesDir='LabViewFiles\';
Subjname=[Filenames{n}(1:end-4)];
Labviewfilenamesdir = [settings.workdir LabViewFilesDir Subjname];
lvmpattern=fullfile(Labviewfilenamesdir, '*.lvm');
dinfo = dir(lvmpattern);

LabViewAllChannels=[];

FsResampLabView=128; % 128 hz resampling

clear lvmstarttime lvmfilestarttime
for ii=1:size(dinfo,1)
    labfile1=[dinfo(ii).folder '\' dinfo(ii).name];
    fid2=fopen(labfile1);
%     daq=textscan(fid2,'%f %f %*f %*f %*f %*f %*f %*f %*f','HeaderLines',11);
%     linenum = 11;
    lvmstarttimetemp= textscan(fid2,'%4s%12q',1,'delimiter','','headerlines',10);
    lvmstarttime(ii,1)=mod(datenum(lvmstarttimetemp{1, 2},'HH:MM:SS.FFF'),1)*86400;
    lvmfilestarttimetemp=textscan(fid2,'%4s%12q',1,'delimiter',' ','headerlines',6);
    lvmfilestarttime(ii,1)=mod(datenum(lvmfilestarttimetemp{1, 2},'HH:MM:SS.FFF'),1)*86400;
%     synchdaq=cell2mat(daq);
    fclose(fid2);
end

TimeLabviewFileStart=find(lvmstarttime==mode(lvmstarttime),1,'first');

daqtemp=[];
for ii=TimeLabviewFileStart:size(dinfo,1)
    labfile1=[dinfo(ii).folder '\' dinfo(ii).name];
    fid2=fopen(labfile1);
    daq=textscan(fid2,'%f %f %f %f %f %f %f %f %f','HeaderLines',23);
    daqtemp=[daqtemp;cell2mat(daq)];
    fclose(fid2);
    clear daq
end

daqall=resample(daqtemp,FsResampLabView,1000);
Synch=SynchSM;
synchsandsamp=resample(Synch,FsResampLabView,256);
clear daqtemp;
daqsynch=daqall(:,2).*100;
synchsandsamp=synchsandsamp;  % extracting only the synch signal from sandman and dividing it by the Sandman amplification factor

% figure; plot(daqsynch); hold on; plot(synchsandsamp,'r');
Synchdelay=finddelay(synchsandsamp,daqsynch);

% keeping sandman in original length and adjusting labview for synchinh..
% if -ive sandman is leading, cut sandman
daqalltemp=daqall;
if settings.fnirs
    [daqfnirs]=fNirsRead(settings,n,Filenames,daqalltemp,FsResampLabView);
    daqfnirs=real(daqfnirs);
    daqall=[];
    daqall=daqfnirs;
    daqalltemp=daqall;
    clear daqfnirs
end
if Synchdelay >=0
    daqalltemp(1:Synchdelay,:)=[];
else
    daqalltemp=[NaN(abs(Synchdelay),(size(daqall,2)));daqall];
end

if size(daqalltemp,1)<length(synchsandsamp) % sandman is more than labview at the end
    daqalltemp(end:length(synchsandsamp),:)=NaN; % padding with nans
else  % Sandman stopped with before labview
    daqalltemp(length(synchsandsamp)+1:end,:)=[]; % discarding labview extras
   
end
figure(111); plot(daqalltemp(:,2)*100); hold on; plot(synchsandsamp,'r');

if settings.fnirs
    clf(figure(111));
    figure(111);
    plot(daqalltemp(:,6)); hold on; plot(daqalltemp(:,11),'r');
end
end





