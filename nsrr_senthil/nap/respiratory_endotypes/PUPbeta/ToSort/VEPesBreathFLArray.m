function [t EdiDrive PesDrive VE a win nota veup hyp pcw t2 e BFLD] = VE_VdriveFLArray(DataOut,n,criteria,BreathFLData)
%%
%build array of ventilation and ventilatory drive data specific to PUP format

NExcludepostAR = 2;
VEcol = 15;
Veupcol = 11;
Errorcol = 9;
ARcol = 10;
tcol = 1;
t2col = 3;
hypcol = 14; 
edicol = 19; %16 [16:Pes, 17:Pmus, 18:VEpes, 19:Edi, 20:VEedi, 21:GGpeak]
pescol = 16;
pmuscol = 17; %17
Ecol = 6;

EdiDrive=[]; %ventilatory drive, normalized locally, EDI
PesDrive= []; %ventilatory drive, normalized locally, Pes
VE=[]; %ventilation, normalized locally
t=[];
t2=[];
a=[];
nota=[]; 
win=[];
veup=[];
hyp=[];
pcw=[];
e=[];
%subject = n;
BFLD=[]; % breath FL data to be returned 

for w=1:length(DataOut{n})   
    if (size(DataOut{n}{w},1)==1&&isnan(DataOut{n}{w}))||...
            isempty(DataOut{n}{w})||...
            isempty(BreathFLData{n}{w})||...
            criteria(w)==0
        continue
    end
    BFLD_Temp=BreathFLData{n}{w}(:,:);
    y1=DataOut{n}{w}(:,VEcol);
    e1=DataOut{n}{w}(:,Ecol);
    x1=DataOut{n}{w}(:,edicol);%-DataOut{n}{w}(:,Errorcol);
    x2=DataOut{n}{w}(:,pescol);
    pcw1=DataOut{n}{w}(:,pmuscol) - DataOut{n}{w}(:,pescol);%-DataOut{n}{w}(:,Errorcol);
    t1=DataOut{n}{w}(:,tcol);
    t21=DataOut{n}{w}(:,t2col);
    a1=ceil(DataOut{n}{w}(:,ARcol)); %arousal
    hyp1=DataOut{n}{w}(:,hypcol);
    % nota is ~a, but 1 to 0 change in a is delayed by two in nota
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
    
    win1 = w+0*x1;
    veup1 = DataOut{n}{w}(:,Veupcol);
    
    EdiDrive = [EdiDrive;x1];
    PesDrive = [PesDrive;x2];
    e = [e;e1];
    VE = [VE;y1];
    t = [t;t1];
    t2 = [t2;t21];
    a = [a;a1];
    win = [win;win1];
    nota=[nota;nota1];
    veup=[veup;veup1];
    hyp = [hyp;hyp1];
    pcw = [pcw;pcw1];
    BFLD = [BFLD;BFLD_Temp];
if 0
figure(100)
stairs(t1,[x1 y1 2+0.3*a1]);
hold('on')
end
end
t2 = t2-t;

%% DataOut =[
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
% 24 FlowEdi_VI];


