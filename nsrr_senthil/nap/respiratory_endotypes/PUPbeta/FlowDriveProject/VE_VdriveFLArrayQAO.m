function [BB_time, EdiDrive, PesDrive, VE, Ar, win, notAr, Veup, hypnog, pcw, BB_Tot, BFLD, EType, ApneaB] =...
    VE_VdriveFLArrayQAO(pt,criteria,BreathFLData, BreathData)
%% build array of ventilation and ventilatory drive data specific to PUP format

NExcludepostAR = 2;
% BreathData    DataOut (complete listing at end of fn)
tcol=   2;      % 1  
t2col=  4;      % 3
Ecol=   13;     % 6
ARcol=  9;      % 10
hypcol= 12;     % 14
VEcol=  23;     % 15
Veupcol=24;     % 11
pescol= 14;     % 16
pmuscol=15;     % 17
edicol= 16;     % 19

EdiDrive=[]; %ventilatory drive, normalized locally, EDI
PesDrive= []; %ventilatory drive, normalized locally, Pes
VE=[]; %ventilation, normalized locally
BB_time=[];
BB_Tot=[];
Ar=[];
notAr=[];
win=[];
Veup=[];
hypnog=[];
pcw=[];
EType=[]; % event type, 5=Central 6=Obstructive 7=OHypopnea 15=Mixed
ApneaB=[]; % 1 if VEfromFlow scored this breath as Apnea

BFLD=[]; % breath FL data to be returned (as a table)

for w=1:length(BreathData{pt})
    if (size(BreathData{pt}{w},1)==1&&isnan(BreathData{pt}{w}))||...
            isempty(BreathData{pt}{w})||...
            isempty(BreathData{pt}{w})||...
            isempty(BreathFLData{pt}{w})||...
            criteria(w)==0
        continue
    end
    
    % convert BreathData table to mat
    BreathDataMat = table2array(BreathData{pt}{w});
    if size(BreathDataMat,2)>24
        ApneaB_Temp=BreathDataMat(:,25);
    else 
        ApneaB_Temp = NaN;
    end
    EType_Temp=BreathDataMat(:,Ecol);
    y1=BreathDataMat(:,VEcol);
    x1=BreathDataMat(:,edicol);
    x2=BreathDataMat(:,pescol);
    pcw1=BreathDataMat(:,pmuscol) - BreathDataMat(:,pescol);
    t1=BreathDataMat(:,tcol);
    t21=BreathDataMat(:,t2col);
    a1=ceil(BreathDataMat(:,ARcol)); %arousal
    hyp1=BreathDataMat(:,hypcol);
    % nota is ~a, but 1 to 0 change in a is delayed (by two) in nota
    nota1=1-a1;
    nota1(1:NExcludepostAR)=0;
    aoffset = -[NaN;diff(a1)];
    I=find(aoffset==1);
    if ~isempty(I)
        for i=1:length(I)
            li=I(i);
            ri=I(i)+NExcludepostAR-1;
            if ri>length(x1), ri=length(x1); end
            nota1(li:ri)=0;
        end
    end
    
    win1 = w+0*t1;
    veup1 = BreathDataMat(:,Veupcol);
    
    % BreathFLData
    BFLD_Temp=BreathFLData{pt}{w}(:,:);
    
    EdiDrive = [EdiDrive;x1];
    PesDrive = [PesDrive;x2];
    VE = [VE;y1];
    BB_time = [BB_time;t1];
    BB_Tot = [BB_Tot;t21];
    Ar = [Ar;a1];
    win = [win;win1];
    notAr=[notAr;nota1];
    Veup=[Veup;veup1];
    hypnog = [hypnog;hyp1];
    pcw = [pcw;pcw1];
    BFLD = [BFLD;BFLD_Temp];
    EType = [EType;EType_Temp];
    ApneaB = [ApneaB;ApneaB_Temp];
    
    if 0
        figure(100)
        stairs(t1,[x1 y1 2+0.3*a1]);
        hold('on')
    end
end

BB_Tot = BB_Tot-BB_time;

%% BreathData contains the following data:
%
% 1 'Time0'
% 2 'Time_start'
% 3 'Time_mid'
% 4 'Time_end'
% 5 'BB_i_start'
% 6 'BB_i_mid'
% 7 'BB_i_end'
% 8 'VI' (normalised)
% 9 'AR'
% 10 'spo2'
% 11 'pos_B'
% 12 'hypnog_B'
% 13 'Etype'
% 14 'DeltaPes'
% 15 'DeltaPmus'
% 16 'DeltaEdi'
% 17 'VIpes'
% 18 'VIedi'
% 19 'GGpeak'
% 20 'GGtonic'
% 21 'FlowPes_VI'
% 22 'FlowEdi_VI'
% 23 VE
% 24 Veup
% 25 ApneaB

%% DataOut contains the following:
%
% 1 Time(BB_i_start)
% 2 Time(BB_i_mid)
% 3 Time(BB_i_end)
% 4 VI'
% 5 Vdr_est
% 6 E1'
% 7 E_recover'
% 8 E_Terminate'
% 9 Error
% 10 AR
% 11 meanVIbeforenormalizing+0*AR
% 12 VAr_est
% 13 pos_B
% 14 hypnog_B
% 15 meanVIbeforenormalizing*VI'
% 16 DeltaPes
% 17 DeltaPmus
% 18 VIpes
% 19 DeltaEdi
% 20 VIedi
% 21 GGpeak
% 22 GGtonic
% 23 FlowPes_VI
% 24 FlowEdi_VI

