function [daqfnirs]=fNirsRead(settings,n,Filenames,daqalltemp,FsResampLabView)
fnirsFilesDir='fNIRSFinalFiles\';
Subjname=[Filenames{n}(1:end-4)];
fnirsfilenamesdir = [settings.workdir fnirsFilesDir Subjname];
fnirspattern=fullfile(fnirsfilenamesdir, '*.mat');
dinfo1 = dir(fnirspattern); % fnirs files

% initialize labview colums for incoming fnirs data
daqalltemp(:,10:30)=NaN;
fnirstemp=[];
fnirsmarker=daqalltemp(:,6);   % extract the marker signal for fnirs on-off
a=round(fnirsmarker);  % round off the digits to get only 0s and 5s
ixT=find(diff(fnirsmarker)>3);
ixTD=find(diff(fnirsmarker)<-3);
% figure; plot(diff(fnirsmarker));
daqallbackup=daqalltemp;

try  
    
for jj=1:size(dinfo1,1)
   
    if jj==1
        clear fnirstemp temp
        labfile=[dinfo1(jj).folder '\' dinfo1(jj).name];
        temp = load(labfile);
        fnirstemp=[temp.fNIRS]; % 22 channels here
          
%         ix=find(diff(fnirsmarker)>3,1,'first'); % marker going up for 1st time
        
        ix=ixT(1); % marker going up for 1st time
%         if ix<360*FsResampLabView
%             ixt=find(ixT>(360*FsResampLabView),1,'first'); % giving a buffer time of 6min in the
%             % beginning to remove any noise or spikes
%             ix=ixT(ixt);
%         end
        
        daqallbackup(ix:ix+size(fnirstemp,1)-1,10:30)=fnirstemp(:,2:end);
        TotLen=ix+size(fnirstemp,1); 
        PrevEd=find(ixT>=TotLen,1,'first');
        PrevEdInd=ixT(PrevEd);
        NxtSt=find(ixTD>=PrevEdInd,1,'first');
        NxtStInd=ixTD(NxtSt);
        figure(333); clf(figure(333));
       plot(daqallbackup(:,1),daqallbackup(:,6)); 
        hold on; plot(daqallbackup(:,1),daqallbackup(:,11));
        
    else
        
        clear fnirstemp temp
        labfile=[dinfo1(jj).folder '\' dinfo1(jj).name];
        temp = load(labfile);
        fnirstemp=[temp.fNIRS]; % 22 channels here
        daqallbackup(NxtStInd:NxtStInd+size(fnirstemp,1)-1,10:30)=fnirstemp(:,2:end);
        TotLen=NxtStInd+size(fnirstemp,1); 
        PrevEd=find(ixT>=TotLen,1,'first');
        PrevEdInd=ixT(PrevEd);
        NxtSt=find(ixTD>=TotLen,1,'first');
        NxtStInd=ixTD(NxtSt);
       
    end
    
end

figure(333); clf(figure(333));
plot(daqallbackup(:,1),daqallbackup(:,6));
hold on; plot(daqallbackup(:,1),daqallbackup(:,11),'r');

daqfnirs=daqallbackup;

catch
    disp('fnirs synchronization failed');
    daqfnirs=daqalltemp;
end
