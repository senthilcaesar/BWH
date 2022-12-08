function [BB_time, EdiDrive, PesDrive, VE, Ar, win, notAr, Veup, hypnog, pcw, BB_Tot, BFLD, EType, ApneaB] =...
    VE_VdriveFLArray(pt,criteria,BreathFLData, BreathData, FeatureSize)
%% extract ventilation and ventilatory drive data from the BreathData tables
% this version works entirely within the tables (no conversion to mat)

%% settings
NExcludepostAR = 2;

%% variables
BB_time=[];
EdiDrive=[];    %ventilatory drive, normalized locally, EDI
PesDrive= [];   %ventilatory drive, normalized locally, Pes
VE=[];          %ventilation, normalized locally
Ar=[];
win=[];
notAr=[];
Veup=[];
hypnog=[];
pcw=[];
BB_Tot=[];
BFLD=[];        % breath FL data to be returned (as a table)
EType=[];       % event type, 5=Central 6=Obstructive 7=OHypopnea 15=Mixed
ApneaB=[];      % 1 if VEfromFlow scored this breath as Apnea

%% Anon function to keep code tidy
getcol = @(x,y) find(strcmp(x,y)==1);

%% process
for w=1:length(BreathData{pt})
    if (size(BreathData{pt}{w},1)==1&&isnan(BreathData{pt}{w}))||...
            isempty(BreathData{pt}{w})||...
            isempty(BreathData{pt}{w})||...
            isempty(BreathFLData{pt}{w})||...
            criteria(w)==0
        continue
    end
    
    % get the table variable names
    VarNames = BreathData{pt}{w}.Properties.VariableNames;
    
    if size(BreathData{pt}{w},2)<24
        ApneaB_Temp = NaN;
        disp('--- WARNING --- BreathData is small. Probable old version.');
    else 
        ApneaB_Temp=BreathData{pt}{w}{:,getcol(VarNames, 'ApneaB')};
    end
    EType_Temp=BreathData{pt}{w}{:,getcol(VarNames, 'Etype')};
    VE_temp=BreathData{pt}{w}{:,getcol(VarNames, 'VE')};
    x1=BreathData{pt}{w}{:,getcol(VarNames, 'DeltaEdi')};
    x2=BreathData{pt}{w}{:,getcol(VarNames, 'DeltaPes')};
    pcw1=BreathData{pt}{w}{:,getcol(VarNames, 'DeltaPmus')} - BreathData{pt}{w}{:,getcol(VarNames, 'DeltaPes')};
    t1=BreathData{pt}{w}{:,getcol(VarNames, 'Time_start')};
    t21=BreathData{pt}{w}{:,getcol(VarNames, 'Time_end')};
    a1=ceil(BreathData{pt}{w}{:,getcol(VarNames, 'ARei')}); %arousal
    hyp_Temp=BreathData{pt}{w}{:,getcol(VarNames, 'hypnog_B')};
    Veup_Temp = BreathData{pt}{w}{:,getcol(VarNames, 'Veup')};
    
    
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
    
    win_Temp = w+0*t1;
    
    % BreathFLData
    BFLD_Temp=BreathFLData{pt}{w}(:,:);
    if size(BFLD_Temp,2) ~= FeatureSize
        disp('!');
    end
    
    %% add to growing data
    EdiDrive = [EdiDrive;x1];
    PesDrive = [PesDrive;x2];
    VE = [VE;VE_temp];
    BB_time = [BB_time;t1];
    BB_Tot = [BB_Tot;t21];
    Ar = [Ar;a1];
    win = [win;win_Temp];
    notAr=[notAr;nota1];
    Veup=[Veup;Veup_Temp];
    hypnog = [hypnog;hyp_Temp];
    pcw = [pcw;pcw1];
    BFLD = [BFLD;BFLD_Temp];
    EType = [EType;EType_Temp];
    ApneaB = [ApneaB;ApneaB_Temp];
    
    if 0
        figure(100); clf(figure(100));
        stairs(t1,[x1/100 VE_temp 0.3*a1]);
        hold('on')
    end
end

BB_Tot = BB_Tot-BB_time;

%% BreathData contains 
% 'Time0','Time_start','Time_mid','Time_end','BB_i_start','BB_i_mid','BB_i_end',...
% 'VI','ARei','AReiF','ARie','ARieF','spo2','pos_B','hypnog_B','Etype','DeltaPes','DeltaPmus','DeltaEdi','VIpes','VIedi',...
% 'GGpeak','GGtonic','FlowPes_VI','FlowEdi_VI','VE','Veup','Pdelta','Ptheta','Palpha','Psigma','Pbeta','WakeSleep','ApneaB'};

