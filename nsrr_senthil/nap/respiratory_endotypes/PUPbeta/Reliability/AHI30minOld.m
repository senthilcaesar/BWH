clear all;
close all;
clc;
addpath(genpath('C:\Users\szs88\Partners Healthcare Dropbox\SATP Group\PUPbeta_git\PUPbeta\'));
g = groot;
global settings ChannelsList 
settings.AMasterSpreadsheet = 'D:\MESA\PUPStart\AMasterSpreadsheet';

% read spreadsheet (files worksheet)
[~,patients,~] = xlsread(settings.AMasterSpreadsheet,1,'AD4:AH10003');
% read spreadsheet (options worksheet)
[~,~,raw] = xlsread(settings.AMasterSpreadsheet,2,'C20:C57');
settings.Fs=raw{13};
settings.ignoreCPAPdata=logical(raw{14});
% Read spreadsheet (position codes worksheet)
[~,~,settings.poscodesdatabase] = xlsread(settings.AMasterSpreadsheet,3,'B2:J55');
settings.positioncodesout = PositionCodeLookup(settings.poscodesdatabase,'Output');
[~,~,settings.protocol] = xlsread(settings.AMasterSpreadsheet,1,'AF4:AF10003');
settings.positioncodes = PositionCodeLookup(settings.poscodesdatabase,settings.protocol{1});
settings.supinepositioncode = settings.positioncodes(1);

%
seg=30; % 30min windows

for n=1:size(patients,1)
    n
    settings.filename=[char(patients(n,1))];
    disp(['Subject: ', settings.filename])
    directoryn=char(patients(n,2));
    if directoryn(end)~=filesep
        directoryn=[directoryn filesep];
    end
    MATfilename=[directoryn char(patients(n,1))];
    load(MATfilename);
   
    
    % create odd and even windows
    width=seg*60; % window length in seconds
    t=(0:(1/settings.Fs):(size(DataEventHypnog_Mat,1)-1)*(1/settings.Fs))';
    D1=[0:2*width:size(DataEventHypnog_Mat,1)]; % duty cycle of pulse;off time and on time are same.
    D1=D1+width/2; % to account for duty cycle starting default at 1/2 of D1 from 0.
    winn_o = pulstran(t,D1,@rectpuls,width); % odd window
    winn_e=1-winn_o; % even window
    winn_all=ones(size(DataEventHypnog_Mat,1),1); % whole night
    %     figure; plot(t,winn_o); hold on; plot(t,winn_e,'r');
    clear t;
    
    
    % AHI for odd and even windows
    try
        if 0
        [AHIData_odd{n},Evts_odd{n}] = getAHI(DataEventHypnog_Mat,winn_o);
        [AHIData_even{n},Evts_even{n}] = getAHI(DataEventHypnog_Mat,winn_e);
        [AHIData_all{n},Evts_all{n}] = getAHI(DataEventHypnog_Mat,winn_all);
        else
           
                EpochsXHz=DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'Epochs')==1));
                Time=DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'Time')==1));
                [Position,PositionRaw] = getPos(DataEventHypnog_Mat,ChannelsList,settings);
                [AHI.Total,Evts,rowNames]=getAHIEvtSubset(Evts,EpochsXHz,Time,Position,winn_all);
                AHITable1=array2table(AHI.Total,'VariableNames',rowNames,'RowNames',{'AHITotal'});
                
                %Event subset, e.g. 4pc AHI
                desatlist=[3 4]; %note 4percent will not include arousal
                           
                [AHI.ThreePA,Evts,~] = getDesatArSubset(Evts,EpochsXHz,Time,Position,CPAPoff,desatlist(1));
                [AHI.FourP,Evts,rowNames] = getDesatArSubset(Evts,EpochsXHz,Time,Position,CPAPoff,desatlist(2));
                AHITable2=array2table(AHI.ThreePA,'VariableNames',rowNames,'RowNames',{'AHI3PA'});
                AHITable3=array2table(AHI.FourP,'VariableNames',rowNames,'RowNames',{'AHI4P'});
                AHIdata2{1}=[AHITable1;AHITable2;AHITable3];
        end
        
        [AllPositionsEventsTable,AllPositionsEventsTableN,SupineEventsTable,EvtDurTable] = AHIdata2Tbls(AHIData_all(n));
        
        AHIodd_total{n}=AHIData_odd{n}(58);
%         TSTodd_total(n)=AHIData_odd{n}(59)
        AHIeven_total{n}=AHIData_even{n}(58);
%          TSTeven_total(n)=AHIData_even{n}(59)
        AHIall_total{n}=AHIData_all{n}(58);
%          TSTall_total(n)=AHIData_all{n}(59)
    catch
        Nantmp1={NaN(1,size(AHIData_odd{1},2)) NaN(1,size(AHIData_odd{1},2))};
        [AHIData_odd{n},Evts_odd{n}] =  Nantmp1{:};
        [AHIData_even{n},Evts_even{n}] = Nantmp1{:};
        [AHIData_all{n},Evts_all{n}] = Nantmp1{:};
        
        AHIodd_total{n}=NaN;
        AHIeven_total{n}=NaN;
        AHIall_total{n}=NaN;
%         TSTodd_total(n)=NaN;
%         TSTeven_total(n)=NaN;
%         TSTall_total(n)=NaN;
         
    end
    
    
    %ODI
    try
        SaO2=DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'SpO2')==1));
        Include_all=1-(SaO2==0); % include only when sao2 is non-zero; use this for whole night
        SaO2(SaO2==0)=NaN;
        Include_o=zeros(size(DataEventHypnog_Mat,1),1);
        Include_e=zeros(size(DataEventHypnog_Mat,1),1);
        Include_o(Include_all==1&winn_o==1)=1; % odd window
        Include_e(Include_all==1&winn_e==1)=1; % even window
        dt = 1/settings.Fs;
        
        [ODI3_all{n},ODI4_all{n},~,~,~,~,~,~,~,~]=CalcODI(SaO2,dt,Include_all);
        [ODI3_odd{n},ODI4_odd{n},~,~,~,~,~,~,~,~]=CalcODI(SaO2,dt,Include_o);
        [ODI3_even{n},ODI4_even{n},~,~,~,~,~,~,~,~]=CalcODI(SaO2,dt,Include_e);
    catch
        Nantmp2={NaN NaN NaN NaN NaN NaN};
        [ODI3_all{n},ODI4_all{n},ODI3_odd{n},ODI4_odd{n},ODI3_even{n},ODI4_even{n}]=Nantmp2{:};
        
    end
end