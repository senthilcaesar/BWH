function [bT,bSigT] = VEfromFlowT(time,Vflow,settings1,ploton)
global settings n winNum %updated substantially 2016-10-05 

% fakeIEratio=1.5
% if fakeIEratio~=1
%         Vflow(Vflow>0)=Vflow(Vflow>0)*(fakeIEratio^0.5);
%         Vflow(Vflow<0)=Vflow(Vflow<0)/(fakeIEratio^0.5);  
% end

if exist('settings1')
    settings=settings1;
end

%start attempting to handle some NaN data
nantracker = isnan(Vflow);
Vflow(nantracker)=0;

leaksignalfilt = 0*time;
if ~exist('Vflow_backup')
    Vflow_backup = Vflow;
end
if 0
    Vflow = Vflow_backup(:)';
end

minimum_figs=1-settings.plotbreathdetectionfigures;

Vflow=Vflow(:)';
time=time(:)';

N = length(Vflow);
dt=(time(end)-time(1))/(length(time)-1);

%% Offset such that the mean nasal pressure signal=0
detrend_flow=0;
if detrend_flow
    leak1=mean(Vflow);
    Vflow=Vflow-mean(Vflow);
else
    leak1=0;
end

%% low pass before sqrt transform
if 1
    filter_HFcutoff_butter0 = 6;
    filter_order0 = 2;
    [B_butter0,A_butter0] = butter(filter_order0,[filter_HFcutoff_butter0]/(1/dt/2),'low');
    Vflow = filtfilt(B_butter0,A_butter0,Vflow); %filtfilt, otherwise flow signal is right-shifted
end

%% Square root scaling of the pressure signal
% to make it a better approximate of a true flow signal:
if settings.sqrt_scaling==1
        [Vflow,leak,IEratio,IEratioEstimated]=sqrtscaling(time,Vflow,settings.scalingexponent,settings.plotfiguresqrtscaling);
        leak=leak+leak1;
        Vflow_backup2=Vflow;
elseif settings.sqrt_scaling==2 %already did drift removal previously
    
    
        if 0
        [~,~,IEratio_originalmethod,~]=sqrtscaling(time,Vflow,settings.scalingexponent,settings.plotfiguresqrtscaling);
        IEratio_originalmethod
        end
        
        
        Vflow(Vflow>0)=Vflow(Vflow>0).^settings.scalingexponent;
        Vflow(Vflow<0)=-((-Vflow(Vflow<0)).^settings.scalingexponent);
        
        %find IEratio
        medianlength=18;
        medianstepsize=1;
        Noverlaps=round(medianlength*settings.Fs);
        Nstepsize=round(medianstepsize*settings.Fs);
        
        FlowBuffer = buffer(Vflow,Noverlaps,Noverlaps-Nstepsize,'nodelay');
        I = FlowBuffer(:,end)==0; FlowBuffer(I,end)=NaN; %removes zeros in buffer on last window
        FlowBufferI=FlowBuffer;
        FlowBufferI(FlowBufferI<0)=NaN;
        FlowBufferE=FlowBuffer;
        FlowBufferE(FlowBufferE>0)=NaN;
        A = table();
        A.Ai = nanmean(FlowBufferI)';
        A.Ae = nanmean(-FlowBufferE)';
        A.Atot = A.Ai + A.Ae;
        I = A.Atot<nanmean(0.2*A.Atot);
        A{I==1,:}=NaN;
        
        if ploton
            figure(99)
            TimeDrift = time(1) + medianlength/2 + [0:medianstepsize:medianstepsize*(size(FlowBuffer,2)-1)]';
            subplot(4,1,1); plot(time,Vflow_backup,'k',time,0*Vflow,'r');
            subplot(4,1,2); plot(time,Vflow,'k',time,0*Vflow,'r');
            subplot(4,1,3); plot(TimeDrift,[A.Ai A.Ae]);
        end
        
        %nanmean(tempI>tempE)
        IEratioEstimated = nanmedian(A.Ai./A.Ae);
        
        %ConservativeIE 
        tempcentiles=[10:5:90];
        centilesx=prctile(A.Ai./A.Ae,tempcentiles);
        
        ConservativeIE=1;
        if centilesx(end)<1
            ConservativeIE=centilesx(end);
        elseif centilesx(1)>1
            ConservativeIE=centilesx(1);
        end
        IEratioEstimated=ConservativeIE;
        
        if ploton
        figure(56); plot(tempcentiles,centilesx)
        end
        
        expinspcorrectionlimit=2; 
        IEratio=IEratioEstimated;
        if IEratio>expinspcorrectionlimit
            IEratio=expinspcorrectionlimit;
        elseif IEratio<1/expinspcorrectionlimit
            IEratio=1/expinspcorrectionlimit;
        end
        
        Vflow(Vflow>0)=Vflow(Vflow>0)/(IEratio^0.5);
        Vflow(Vflow<0)=Vflow(Vflow<0)*(IEratio^0.5);        

        leak=0;
        Vflow_backup2=Vflow;
        if ploton
           figure(99);
           subplot(4,1,4); plot(time,Vflow,'k',time,0*Vflow,'r');
        end
else
    leak=leak1; IEratio=1;
    Vflow_backup2=Vflow;
    IEratioEstimated=NaN;
end

%% Gentle filter for analysis
if 1 %already filtered
    filter_HFcutoff_butter0 = 4;
    filter_order0 = 4;
    [B_butter0,A_butter0] = butter(filter_order0,[filter_HFcutoff_butter0]/(1/dt/2),'low');
    Vflow_backup2 = filtfilt(B_butter0,A_butter0,Vflow_backup2);
    %VFlow_filtered1 = filtfilt(B_butter0,A_butter0,Vflow);
end

if length(Vflow)<1000 && ~(isfield(settings, 'AllowLessThan10Secs') && settings.AllowLessThan10Secs==1)
    time=NaN;
    Vflow_out=NaN;
    BB_i_start=NaN;
    BB_i_mid=NaN;
    BB_i_end=NaN;
    BB_t=NaN;
    VI=NaN;
    VE=NaN;
    Ttot=NaN;
    leak=NaN;
    IEratio=NaN;
    VT=NaN;
    Vpeak=NaN;
    Vpeakmean=NaN;
    Apnea_B=NaN;
    leak_B=NaN;
    
    bSigT = table();
    bSigT.time=time(:);
    bSigT.Vflow=Vflow(:);
    bSigT.Vflow_out=Vflow_out(:);
    bT = table();
    bT.BB_i_start=BB_i_start(:);
    bT.BB_i_mid=BB_i_mid(:);
    bT.BB_i_end=BB_i_end(:);
    bT.BB_t=BB_t(:);
    bT.VI=VI(:);
    bT.VE=VE(:);
    bT.Ttot=Ttot(:);
    bT.leak=leak(:) + 0*Ttot(:);
    bT.IEratio=IEratio(:) + 0*Ttot(:);
    bT.VT=VT(:);
    bT.Vpeak=Vpeak(:);
    bT.Vpeakmean=Vpeakmean(:);
    bT.Apnea_B=Apnea_B(:);
    bT.leak_B=leak_B(:);
    
    return;
end

%Filter gently for use in analysis
% if 0
% filter_LFcutoff_butter1 = 1/10; %original 1/10
% filter_HFcutoff_butter1 = 3; %original 3
% filter_order1 = 4;
% [B_butter1,A_butter1] = butter(filter_order1,[filter_LFcutoff_butter1 filter_HFcutoff_butter1]/(1/dt/2));
% else
% filter_HFcutoff_butter1 = 3; %original 3
% filter_order1 = 4;
% [B_butter1,A_butter1] = butter(filter_order1,[filter_HFcutoff_butter1]/(1/dt/2),'low');
% end
% VFlow_filtered1 = filter(B_butter1,A_butter1,Vflow);

%% Aggressive low pass for apnea detection
% The Butterworth filter has a good frequency-response. Frequency response
% in the passband is as flat as possible. i.e. it has a sharp shoulder, and
% a steep roll-off. However, it has relatively poor time-delay performance,
% and will incur a time shift.
% In contrast, an equivalent order Bessel filter will have a more curved
% shoulder and roll-off in frequency-response, but its phase shift varies
% linearly with frequency; this is equivalent to a constant time delay.
% See Horowitz & Hill, The Art of Electronics, 2nd Ed, Ch 5, pp 269-272.
% So when we apply this (as we do in the processing steps) in the forward
% and reverse directions, it produces a better in-phase filtered signal.

if 1 %low pass filter signal
    FilterType = 'Butter'; %'Bessel';  %
    filter_HFcutoff_butter0 = 1.2; %0.75
    switch FilterType
        case 'Butter'
            filter_order0 = 2;
            [num,den] = butter(filter_order0,[filter_HFcutoff_butter0]/(1/dt/2),'low');
        case 'Bessel'
            [b,a] = besself(2,4*pi*filter_HFcutoff_butter0);
            [num,den] = bilinear(b,a,1/dt);
    end
    Vflow2 = filtfilt(num,den,Vflow_backup2(:)');
end

if ~minimum_figs
    figure(33); clf(figure(33));
    set(gcf,'color',[1 1 1]);
    ax1(1)=subplot(2,1,1); plot(time,cumsum(Vflow_backup2(:)'-leaksignalfilt(:)')*dt); box('off');
    ax1(2)=subplot(2,1,2); plot(time,Vflow,'r-',time,Vflow2,'b-'); box('off'); hold on;
    %     plot(time(up_indx),0, 'r^');
    %     plot(time(down_indx),0, 'rv');
    %     refline(0,0);
    linkaxes(ax1,'x');
end

vol=cumsum(Vflow_backup2(:)'-leaksignalfilt(:)')*dt;
Vflow = gradient(vol)/dt;
%Vflow_out = Vflow;

%% Loop
for loopnumber=1:2 %second loop subtracts extreme leaks
    if ~minimum_figs
        figure(33);
    end
    vol=cumsum(Vflow(:)'-leaksignalfilt(:)')*dt;
    Vflow = gradient(vol)/dt;
    %% detecting peaks (0.063440 seconds)
    
    % accompanying filtered flow trace
    %Vflow = [diff(vol)/dt]; Vflow=[Vflow(1) Vflow];
    
    [min_list,max_list] = peakdet(-vol,0.001*std(vol)); %in this configuration, max_list(i,1) > minlist(i,1).
    
    BB_i_start=min_list(:,1);
    BB_i_mid=max_list(:,1);
    if BB_i_start(1)==1
        BB_i_start(1)=[];
    end
    if BB_i_start(end)==length(vol)
        BB_i_start(end)=[];
    end
    if BB_i_start(1)>BB_i_mid(1)
        BB_i_mid(1)=[];
    end
    while BB_i_mid(end)>BB_i_start(end)
        BB_i_mid(end)=[];
    end
    
    if ~minimum_figs
        ax1(1)=subplot(1+1,1,1); plot(time,vol,time(BB_i_start),vol(BB_i_start),'r.-',time(BB_i_mid),vol(BB_i_mid),'k.-'); box('off');
        ax1(1+1)=subplot(1+1,1,1+1); plot(time,Vflow,time(BB_i_start),Vflow(BB_i_start),'r.',time(BB_i_mid),Vflow(BB_i_mid),'k.'); box('off');
    end
    
    %estimate normal tidal volume from the distance between upper and lower
    %vols.
    
    if ~minimum_figs
        ax1(1)=subplot(1+1,1,1); plot(time,vol,time(BB_i_start),vol(BB_i_start),'r.',time(BB_i_mid),vol(BB_i_mid),'k.'); box('off');
        ax1(1+1)=subplot(1+1,1,1+1); plot(time,Vflow,time,Vflow,time(BB_i_start),Vflow(BB_i_start),'r.',time(BB_i_mid),Vflow(BB_i_mid),'k.'); box('off');
    end
    
    BB_i_end = BB_i_start(2:end);
    BB_i_start = BB_i_start(1:end-1);
    
    VTi = vol(BB_i_mid) - vol(BB_i_start);
    VTe = vol(BB_i_mid) - vol(BB_i_end);
    Ti = (BB_i_mid-BB_i_start)'*dt;
    Te = (BB_i_end-BB_i_mid)'*dt;
    Ttot = Ti+Te;
    VT = (VTi.*Te+VTe.*Ti)./(Ttot);
    VTmean = sum(VT.*(Ttot))/sum(Ttot);
    criteriafornormal = VT>(0.2*VTmean)&Ttot<10;
    normalTtot=median(Ttot(criteriafornormal));
    %inverted data
    Te_previous = [NaN Te(1:end-1)];
    VTe_previous = [NaN VTe(1:end-1)];
    VTinverted = (VTi.*Te_previous+VTe_previous.*Ti)./(Ti+Te_previous);
    
    %% Remove extreme mini-breaths (0.000116 seconds)
    %heirarchical removal of mini breath patterns I-E or E-I (inverted VT), slow but very reliable, speed depends on noise in signal
    VT_thres = VTmean*0.05;
    while 1
        [minVTinspexp,minVTinspexp_i] = min(VT);
        [minVTinverted,minVTinverted_i] = min(VTinverted);
        [minVT,VTpattern] = min([minVTinspexp,minVTinverted]);
        
        if minVT>VT_thres
            break
        end
        
        if VTpattern==2
            i=minVTinverted_i;
            %time(BB_i_start(i))
            [~,maxi]=max([vol(BB_i_mid(i-1)),vol(BB_i_mid(i))]);
            BB_i_start(i)=[];
            BB_i_end(i-1)=[];
            BB_i_mid(i+1-maxi)=[]; %removes lower end-insp
        elseif VTpattern==1
            i=minVTinspexp_i;
            %time(BB_i_start(i))
            [~,mini]=min([vol(BB_i_start(i)),vol(BB_i_end(i))]);
            if (i+2-mini)>length(BB_i_start)||i==1
                BB_i_start(i)=[]; %min==1, keep lhs, remove rhs: index = i+1; min==2, keep rhs, remove lhs: index = i;
                BB_i_end(i)=[];
                BB_i_mid(i)=[]; %removes higher end-exp vol
            else
                BB_i_start(i+2-mini)=[]; %min==1, keep lhs, remove rhs: index = i+1; min==2, keep rhs, remove lhs: index = i;
                BB_i_end(i+2-mini-1)=[];
                BB_i_mid(i)=[]; %removes higher end-exp vol
            end
        end
        if ~minimum_figs&&1
            ax1(1)=subplot(1+1,1,1); plot(time,vol,time(BB_i_start),vol(BB_i_start),'r.',time(BB_i_mid),vol(BB_i_mid),'k.'); box('off');
            ax1(1+1)=subplot(1+1,1,1+1); plot(time,Vflow,time,Vflow,time(BB_i_start),Vflow(BB_i_start),'r.',time(BB_i_mid),Vflow(BB_i_mid),'k.'); box('off');
        end
        %recalculate
        VTi = vol(BB_i_mid) - vol(BB_i_start);
        VTe = vol(BB_i_mid) - vol(BB_i_end);
        Ti = (BB_i_mid-BB_i_start)'*dt;
        Te = (BB_i_end-BB_i_mid)'*dt;
        Ttot = Ti+Te;
        VT = (VTi.*Te+VTe.*Ti)./(Ti+Te);
        VTmean = sum(VT.*(Ttot))/sum(Ttot);
        VT_thres = VTmean*0.04;
        %inverted data
        Te_previous = [NaN Te(1:end-1)];
        VTe_previous = [NaN VTe(1:end-1)];
        VTinverted = (VTi.*Te_previous+VTe_previous.*Ti)./(Ti+Te_previous);
        %pause
    end
    criteriafornormal = VT>(0.2*VTmean)&Ttot<10;
    normalTtot=median(Ttot(criteriafornormal));
    normalTi=median(Ti(criteriafornormal));
    
    if ~minimum_figs
        ax1(1)=subplot(1+1,1,1); plot(time,vol,time(BB_i_start),vol(BB_i_start),'r.',time(BB_i_mid),vol(BB_i_mid),'k.'); box('off');
        ax1(1+1)=subplot(1+1,1,1+1); plot(time,Vflow,time,Vflow,time(BB_i_start),Vflow(BB_i_start),'r.',time(BB_i_mid),Vflow(BB_i_mid),'k.'); box('off');
    end
    
    %% If we find a short or small inspiration, presume it is part of the previous expiration (0.000073 seconds)
    VTi_thres = VTmean*0.05;
    Ti_thres = normalTtot/10;
    
    criteria=find((Ti<Ti_thres)|(VTi<VTi_thres*0.67)|((Ti<Ti_thres*2)&(VTi<VTi_thres*2)));
    %firstbreathisbad=sum(criteria==1);
    if ~isempty(criteria)
        firstbreathisbad=criteria(1)==1;
        while firstbreathisbad
            criteria(1)=[];
            BB_i_start(1)=[];
            BB_i_mid(1)=[];
            BB_i_end(1)=[];
            criteria=criteria-1;
            if isempty(criteria)
                break
            end
            firstbreathisbad=criteria(1)==1;
        end
        criteria2=criteria-1;
        BB_i_start(criteria)=[];
        BB_i_mid(criteria)=[];
        BB_i_end(criteria2)=[];
    end
    
    % recalculate breath info (0.001917 seconds)
    
    for i=1:length(BB_i_start)
        [~,tempi]=max(vol(BB_i_start(i):BB_i_end(i)));
        BB_i_mid(i,1)=BB_i_start(i)+tempi-1;
    end
    
    if ~minimum_figs
        ax1(1)=subplot(1+1,1,1); plot(time,vol,time(BB_i_start),vol(BB_i_start),'r.',time(BB_i_mid),vol(BB_i_mid),'k.'); box('off');
        ax1(1+1)=subplot(1+1,1,1+1); plot(time,Vflow,time(BB_i_start),Vflow(BB_i_start),'r.',time(BB_i_mid),Vflow(BB_i_mid),'k.'); box('off');
        linkaxes(ax1,'x');
    end
    VTi = vol(BB_i_mid) - vol(BB_i_start);
    VTe = vol(BB_i_mid) - vol(BB_i_end);
    Ti = (BB_i_mid-BB_i_start)'*dt;
    Te = (BB_i_end-BB_i_mid)'*dt;
    Ttot = Ti+Te;
    VT = (VTi.*Te+VTe.*Ti)./(Ttot);
    VTmean = sum(VT.*(Ttot))/sum(Ttot);
    criteriafornormal = VT>(0.2*VTmean)&Ttot<10;
    normalTtot=median(Ttot(criteriafornormal));
    normalTe=median(Te(criteriafornormal));
    
    
    %% If we find a small(er) expiration above FRC presume the breath is part of one continuous inspiration:
    
    %VolStartInsp = (VTi-VTe)/VTmean;
    FVTe = VTe./VTi;
    criteria=find(((FVTe<0.2)&Te<normalTe/2)); %VolStartInsp>0.33 %was leaving double peaks on smaller breaths with incl exhalation
    
    if ~isempty(criteria)
        NN=length(BB_i_start);
        lastbreathisbad=criteria(end)==NN;
        while lastbreathisbad
            criteria(end)=[];
            BB_i_start(end)=[];
            BB_i_mid(end)=[];
            BB_i_end(end)=[];
            NN=NN-1;
            if isempty(criteria)
                break
            end
            lastbreathisbad=criteria(end)==NN;
        end
        criteria2=criteria+1;
        BB_i_start(criteria2)=[];
        BB_i_mid(criteria)=[];
        BB_i_end(criteria)=[];
    end
    
    % recalculate breath info:
    
    for i=1:length(BB_i_start)
        [~,tempi]=max(vol(BB_i_start(i):BB_i_end(i)));
        BB_i_mid(i,1)=BB_i_start(i)+tempi-1;
    end
    
    if ~minimum_figs
        
        ax1(1)=subplot(1+1,1,1); plot(time,vol,time(BB_i_start),vol(BB_i_start),'r.',time(BB_i_mid),vol(BB_i_mid),'k.'); box('off');
        ax1(1+1)=subplot(1+1,1,1+1); plot(time,Vflow,time(BB_i_start),Vflow(BB_i_start),'r.',time(BB_i_mid),Vflow(BB_i_mid),'k.'); box('off');
        linkaxes(ax1,'x');
    end
    VTi = vol(BB_i_mid) - vol(BB_i_start);
    VTe = vol(BB_i_mid) - vol(BB_i_end);
    Ti = (BB_i_mid-BB_i_start)'*dt;
    Te = (BB_i_end-BB_i_mid)'*dt;
    Ttot = Ti+Te;
    VT = (VTi.*Te+VTe.*Ti)./(Ttot);
    VTmean = sum(VT.*(Ttot))/sum(Ttot);
    criteriafornormal = VT>(0.2*VTmean)&Ttot<10;
    normalTtot=median(Ttot(criteriafornormal));
    if sum(criteriafornormal)<0.25*length(Ttot);
        normalTtot=median(Ttot);
    end
    
    
    %% Improve detection of start of inspiration based on the time to 20% of peak Vflow (0.006290 seconds)
    if 0
        Fpeakflow = 0.2;
        BB_Ti_Xi = NaN*BB_i_start;
        
        %find time to 50 percent of peak flow for all breaths to provide a normal value
        for i=1:length(BB_i_start)
            try
                tempflow = (vol(BB_i_start(i):BB_i_mid(i))-vol(BB_i_start(i)));
                tempflow = tempflow/tempflow(end);
                BB_Ti_Xi(i) = find(tempflow>Fpeakflow,1);
            catch me
            end
        end
        
        BB_i_start_new = round(BB_i_mid-(BB_i_mid-BB_i_start-BB_Ti_Xi-1)/(1-Fpeakflow-0.10));
        
        if ~minimum_figs
            ax1(1)=subplot(1+1,1,1); plot(time,vol,time(BB_i_start),vol(BB_i_start),'r.',time(BB_i_mid),vol(BB_i_mid),'k.',time(BB_i_start+BB_Ti_Xi-1),vol(BB_i_start+BB_Ti_Xi-1),'g.',time(BB_i_start_new),vol(BB_i_start_new),'g*'); box('off');
            ax1(1+1)=subplot(1+1,1,1+1); plot(time,Vflow,time,Vflow_apnea0,time,Vflow_apnea,time(BB_i_start),Vflow(BB_i_start),'r.',time(BB_i_mid),Vflow(BB_i_mid),'k.',time(BB_i_start+BB_Ti_Xi-1),Vflow(BB_i_start+BB_Ti_Xi-1),'g.',time(BB_i_start_new),Vflow(BB_i_start_new),'g*'); box('off');
            linkaxes(ax1,'x');
        end
        
        
        Ti_X_upper = prctile(BB_Ti_X,90);
        
        Fearlyinspisapnea=0*BB_i_start;
        for i=1:length(BB_i_start)
            if isnan(BB_Ti_X(i))
                continue
            end
            Fearlyinspisapnea(i)=mean(Vflow_apnea(BB_i_start(i):(BB_i_start(i)+round(BB_Ti_X(i)/dt))));
            if BB_Ti_X(i)>Ti_X_upper&&Fearlyinspisapnea(i)>0.5 %at least half of early insp is apnea
                BB_i_start(i)=BB_i_start(i)+round((BB_Ti_X(i)-Ti_X_upper)/dt);
                if i>1
                    BB_i_end(i-1)=BB_i_start(i);
                end
            end
        end
        
        % recalculate breath info:
        
        if ~minimum_figs
            ax1(1)=subplot(1+1,1,1); plot(time,vol,time(BB_i_start),vol(BB_i_start),'r.',time(BB_i_mid),vol(BB_i_mid),'k.'); box('off');
            ax1(1+1)=subplot(1+1,1,1+1); plot(time,Vflow_apnea0,time,Vflow_apnea,time,Vflow,time,Vflow,time(BB_i_start),Vflow(BB_i_start),'r.',time(BB_i_mid),Vflow(BB_i_mid),'k.'); box('off');
            linkaxes(ax1,'x');
        end
        VTi = vol(BB_i_mid) - vol(BB_i_start);
        VTe = vol(BB_i_mid) - vol(BB_i_end);
        Ti = (BB_i_mid-BB_i_start)'*dt;
        Te = (BB_i_end-BB_i_mid)'*dt;
        Ttot = Ti+Te;
        VT = (VTi.*Te+VTe.*Ti)./(Ttot);
        VTmean = sum(VT.*(Ttot))/sum(Ttot);
        criteriafornormal = VT>(0.2*VTmean)&Ttot<10;
        normalTtot=median(Ttot(criteriafornormal));
        Vdot = VT./Ttot;
        
    end
    %% Find probable apneic periods using "Vflow2" (0.085613 seconds)
    
    if 1&&length(BB_i_start)>5
        if 0
            meanVflow = mean(Vflow(BB_i_start(1):BB_i_end(end)));
        else
            meanVflow = 0;
        end
        slidew = normalTtot*1.5;
        if isnan(normalTtot)
            normalTtot=4;
        end
        slidewi = round(slidew/dt);
        dTi=round(slidewi/10); %was 25, 50 halves the time...
        Nwindows = ceil((length(Vflow)-slidewi)/dTi)+1;
        normalInspFlow = VTmean/median(Ti(criteriafornormal));
        % Scotty changed to 0.05, was 0.08
        thres1 = normalInspFlow*0.067; %%%%%%%%%%%%%%%%%%%%%%%% low amplitude % Mesa version had this changed to 0.2
        thres2 = normalInspFlow*0.5; %%%%%%%%%%%%%%%%%%%%%%%% near flow = 0 (mean)
        Vflow_delta = zeros(Nwindows,1);
        Vflow_median = zeros(Nwindows,1);
        Vflow_apnea = 0*Vflow;
        temptime = zeros(Nwindows,1);
        
        for i=1:Nwindows %this loop is slow
            li = 1+(i-1)*dTi;
            ri = li+slidewi;
            if ri>length(Vflow),ri=length(Vflow); end
            Vflow_delta(i)=std(detrend(Vflow2(li:ri)));%
            Vflow_median(i)=median(Vflow2(li:ri));
            temptime(i) = median(time(li:ri));
            if Vflow_delta(i)<thres1&&(Vflow_median(i)-meanVflow)<thres2 &&(Vflow_median(i)-meanVflow)>-thres2
                Vflow_apnea((li+round(0.10*slidewi)-1):(li+round(0.90*slidewi))) = 1;
            end
        end
        
        % Some settings for DISE data DV 7/30/21
        if (isfield(settings,'DISEdata')&&settings.DISEdata)&&(isfield(settings,'modBBdetect')&&settings.modBBdetect)
            thres1 = 0.025;
            for i=1:Nwindows %this loop is slow
                li = 1+(i-1)*dTi;
                ri = li+slidewi;
                if ri>length(Vflow),ri=length(Vflow); end
                Vflow_delta(i)=std(detrend(Vflow2(li:ri)));%
                Vflow_median(i)=median(Vflow2(li:ri));
                temptime(i) = median(time(li:ri));
                if Vflow_delta(i)<thres1&&(Vflow_median(i)-meanVflow)<thres2 %&&(Vflow_median(i)-meanVflow)>-thres2
                    Vflow_apnea((li+round(0.10*slidewi)-1):(li+round(0.90*slidewi))) = 1;
                end
            end
        end
        
        Vflow_apnea = Vflow_apnea(1:length(Vflow));
        Vflow_apnea0 = Vflow_apnea;
        
        if ~minimum_figs
            ax1(1+1)=subplot(1+1,1,1+1); plot(time,Vflow_apnea0,time,Vflow_apnea,time,Vflow,time,Vflow2,time(BB_i_start),Vflow(BB_i_start),'r.',time(BB_i_mid),Vflow(BB_i_mid),'k.'); box('off');
            linkaxes(ax1,'x');
            hold on
            %temptime(i) = median(time(li:ri));
            plot(temptime,Vflow_delta,'g',temptime,Vflow_median,'b',temptime,thres1+0*temptime,'g:',temptime,thres2+0*temptime,'b:',temptime,-thres2+0*temptime,'b:')
            hold off
        end
    else
        Vflow_apnea = time*0;
        Vflow_apnea0 = time*0;
    end %if detect apneas
    
    
    %% In breaths within apneas - move the end-inspiration to the earliest maxima (0.003859 seconds)
    
    Fapnea = NaN*BB_i_start;
    peak_i_current = BB_i_mid-BB_i_start+1; %location of peak vol for each breath
    peak_i_earliest = peak_i_current;
    Earlierpeak = 0*BB_i_start;
    Fapneathres = 0.67;
    CurrentPeakinApnea = Vflow_apnea0(BB_i_mid)';
    
    for i=1:length(BB_i_start)
        Fapnea(i)=mean(Vflow_apnea0(BB_i_start(i):BB_i_end(i)));
        if Fapnea(i)<=Fapneathres
            continue
        end
        tempvol = vol(BB_i_start(i):BB_i_end(i));
        [max_list,~] = peakdet((tempvol),0.2*std(tempvol));
        if isempty(max_list)
            peak_i_earliest(i)=peak_i_current(i);
            Earlierpeak(i) = 0;
            continue
        end
        peak_i_earliest(i)=max_list(1,1);
        
        Earlierpeak(i) = peak_i_earliest(i)<peak_i_current(i);
        
        %figure(100); plot(tempvol); hold('on');
    end
    
    criteria = Fapnea>Fapneathres&Earlierpeak&Ttot(:)>normalTtot&CurrentPeakinApnea(:)==1;
    
    BB_i_mid(criteria)=BB_i_start(criteria)+peak_i_earliest(criteria)-1;
    
    
    %% Detect start of inspiration at the end of apnea based on the time to 50% of peak Vflow (0.006290 seconds)
    
    Fpeakflow = 0.5;
    BB_Ti_X = NaN*BB_i_start;
    
    %find time to 50 percent of peak flow for all breaths to provide a normal value
    for i=1:length(BB_i_start)
        try
            tempflow = Vflow(BB_i_start(i):BB_i_mid(i));
            temppeakflow = max(tempflow);
            BB_Ti_X(i) = dt*(find(tempflow>temppeakflow*Fpeakflow,1)-1);
        catch me
        end
    end
    
    Ti_X_upper = prctile(BB_Ti_X,90);
    % Why not 90th percentile (conservative)? Or 50th?
    Fearlyinspisapnea=0*BB_i_start;
    for i=1:length(BB_i_start)
        if isnan(BB_Ti_X(i))
            continue
        end
        Fearlyinspisapnea(i)=mean(Vflow_apnea(BB_i_start(i):(BB_i_start(i)+round(BB_Ti_X(i)/dt))));
        if BB_Ti_X(i)>Ti_X_upper&&Fearlyinspisapnea(i)>0.5 %at least half of early insp is apnea
            BB_i_start(i)=BB_i_start(i)+round((BB_Ti_X(i)-Ti_X_upper)/dt);
            if i>1
                BB_i_end(i-1)=BB_i_start(i);
            end
        end
    end
    
    % recalculate breath info:
    
    if ~minimum_figs
        ax1(1)=subplot(1+1,1,1); plot(time,vol,time(BB_i_start),vol(BB_i_start),'r.',time(BB_i_mid),vol(BB_i_mid),'k.'); box('off');
        ax1(1+1)=subplot(1+1,1,1+1); plot(time,Vflow_apnea0,time,Vflow_apnea,time,Vflow,time,Vflow,time(BB_i_start),Vflow(BB_i_start),'r.',time(BB_i_mid),Vflow(BB_i_mid),'k.'); box('off');
        linkaxes(ax1,'x');
    end
    VTi = vol(BB_i_mid) - vol(BB_i_start);
    VTe = vol(BB_i_mid) - vol(BB_i_end);
    Ti = (BB_i_mid-BB_i_start)'*dt;
    Te = (BB_i_end-BB_i_mid)'*dt;
    Ttot = Ti+Te;
    VT = (VTi.*Te+VTe.*Ti)./(Ttot);
    VTmean = sum(VT.*(Ttot))/sum(Ttot);
    criteriafornormal = VT>(0.2*VTmean)&Ttot<10;
    normalTtot=median(Ttot(criteriafornormal));
    Vdot = VT./Ttot;
    normalTi=median(Ti);
    
    %% Try again using time to peak vol
    
    BB_Ti_X = NaN*BB_i_start;
    
    %find time to 50 percent of peak flow for all breaths to provide a normal value
    
    Finspapnea=0*BB_i_start;
    for i=1:length(BB_i_start)
        
        Finspapnea(i)=mean(Vflow_apnea(BB_i_start(i):(BB_i_mid(i))));
        
        %abort if start insp is not within apnea i.e. mid insp pause -- 2018/09/25
        if ~Vflow_apnea(BB_i_start(i))
            continue
        end
        
        %move Ti to
        if Ti(i)>normalTtot&&Finspapnea(i)>0.1
            BB_i_start_new_temp=BB_i_mid(i)-round(1.1*normalTi/dt);
            if i>1
                if BB_i_start_new_temp<BB_i_mid(i-1)
                    BB_i_start_new_temp=BB_i_mid(i-1)+1; %can't go back to before last mid!
                end
                BB_i_start(i)=BB_i_start_new_temp;
                BB_i_end(i-1)=BB_i_start(i);
            end
        end
    end
    
    % recalculate breath info:
    
    if ~minimum_figs
        ax1(1)=subplot(1+1,1,1); plot(time,vol,time(BB_i_start),vol(BB_i_start),'r.',time(BB_i_mid),vol(BB_i_mid),'k.'); box('off');
        ax1(1+1)=subplot(1+1,1,1+1); plot(time,Vflow_apnea0,time,Vflow_apnea,time,Vflow,time,Vflow,time(BB_i_start),Vflow(BB_i_start),'r.',time(BB_i_mid),Vflow(BB_i_mid),'k.'); box('off');
        linkaxes(ax1,'x');
    end
    VTi = vol(BB_i_mid) - vol(BB_i_start);
    VTe = vol(BB_i_mid) - vol(BB_i_end);
    Ti = (BB_i_mid-BB_i_start)'*dt;
    Te = (BB_i_end-BB_i_mid)'*dt;
    Ttot = Ti+Te;
    VT = (VTi.*Te+VTe.*Ti)./(Ttot);
    VTmean = sum(VT.*(Ttot))/sum(Ttot);
    criteriafornormal = VT>(0.2*VTmean)&Ttot<10;
    normalTtot=median(Ttot(criteriafornormal));
    Vdot = VT./Ttot;
    
    %% Handle for apnea detected in long mid inspiratory pause (caused trouble, 2018/11, fixed with thresformovingmidleft=normalTi*1.5)
    handlemidinsppauseapnea=1;
    if handlemidinsppauseapnea
        
        %start insp before apnea, mid inside apnea:
        I = find(Vflow_apnea0(BB_i_start)==0 & Vflow_apnea0(BB_i_mid)==1);
        %assumes single apnea onset inside insp
        thresformovingmidleft=normalTi*1.5;
        for i=1:length(I)
            range = BB_i_start(I(i)):BB_i_mid(I(i));
            I2 = find(diff(Vflow_apnea0(range))==1,1); %indx of apnea onset in inspiration
            if length(I2)==1 %only one new apnea here in insp
                I2_ = I2 + 1 + range(1) - 1; %double correction due to use of diff
                differencetempi = BB_i_mid(I(i))-I2_;
                if differencetempi>(thresformovingmidleft/dt)
                    BB_i_mid(I(i)) = BB_i_mid(I(i)) - differencetempi;
                end
            end
        end
    end
    
    % recalculate breath info:
    
    if ~minimum_figs
        ax1(1)=subplot(1+1,1,1); plot(time,vol,time(BB_i_start),vol(BB_i_start),'r.',time(BB_i_mid),vol(BB_i_mid),'k.'); box('off');
        ax1(1+1)=subplot(1+1,1,1+1); plot(time,Vflow_apnea0,time,Vflow_apnea,time,Vflow,time,Vflow,time(BB_i_start),Vflow(BB_i_start),'r.',time(BB_i_mid),Vflow(BB_i_mid),'k.'); box('off');
        linkaxes(ax1,'x');
    end
    VTi = vol(BB_i_mid) - vol(BB_i_start);
    VTe = vol(BB_i_mid) - vol(BB_i_end);
    Ti = (BB_i_mid-BB_i_start)'*dt;
    Te = (BB_i_end-BB_i_mid)'*dt;
    Ttot = Ti+Te;
    VT = (VTi.*Te+VTe.*Ti)./(Ttot);
    VTmean = sum(VT.*(Ttot))/sum(Ttot);
    criteriafornormal = VT>(0.2*VTmean)&Ttot<10;
    normalTtot=median(Ttot(criteriafornormal));
    normalTi=median(Ti(criteriafornormal));
    Vdot = VT./Ttot;
    
    %% Breaths with low mean insp flow, occurring predominantly in detected apneas => apnea -- delete / merge with prior expiration
    %0.003667 seconds
    
    Apnea_B1 = 0*BB_i_start;
    for i=1:length(BB_i_start)
        if mean(Vflow_apnea0(BB_i_start(i):BB_i_mid(i)))>0.67&&((VTi(i)/Ti(i))/(VTmean/normalTi))<0.2
            Apnea_B1(i)=1;
        end
    end
    
    criteria = find(Apnea_B1==1);
    if ~isempty(criteria)
        firstbreathisbad=criteria(1)==1;
        while firstbreathisbad
            criteria(1)=[];
            BB_i_start(1)=[];
            BB_i_mid(1)=[];
            BB_i_end(1)=[];
            criteria=criteria-1;
            if isempty(criteria)
                break
            end
            firstbreathisbad=criteria(1)==1;
        end
        criteria2=criteria-1;
        BB_i_start(criteria)=[];
        BB_i_mid(criteria)=[];
        BB_i_end(criteria2)=[];
    end
    
    % recalculate breath info:
    
    if ~minimum_figs
        ax1(1)=subplot(1+1,1,1); plot(time,vol,time(BB_i_start),vol(BB_i_start),'r.',time(BB_i_mid),vol(BB_i_mid),'k.'); box('off');
        ax1(1+1)=subplot(1+1,1,1+1); plot(time,Vflow_apnea0,time,Vflow_apnea,time,Vflow,time,Vflow,time(BB_i_start),Vflow(BB_i_start),'r.',time(BB_i_mid),Vflow(BB_i_mid),'k.'); box('off');
        linkaxes(ax1,'x');
    end
    VTi = vol(BB_i_mid) - vol(BB_i_start);
    VTe = vol(BB_i_mid) - vol(BB_i_end);
    Ti = (BB_i_mid-BB_i_start)'*dt;
    Te = (BB_i_end-BB_i_mid)'*dt;
    Ttot = Ti+Te;
    VT = (VTi.*Te+VTe.*Ti)./(Ttot);
    VTmean = sum(VT.*(Ttot))/sum(Ttot);
    criteriafornormal = VT>(0.2*VTmean)&Ttot<10;
    Te = (BB_i_end-BB_i_mid)'*dt;
    normalTi=median(Ti(criteriafornormal));
    normalTe=median(Te(criteriafornormal));
    normalTtot=median(Ttot(criteriafornormal));
    
    %% Probable end-expiratory pause detector (also finds apneas fyi)
    
    if settings.modBB_i_start&&length(BB_i_start)>5

        meanVflow = 0;

        if isnan(normalTtot)
            normalTtot=4;
        end
        slidew = normalTtot*0.05;
        slidewi = round(slidew/dt);
        dTi=round(slidewi/10); %was 25, 50 halves the time...
        if dTi==0
            dTi=1;
        end
        Nwindows = ceil((length(Vflow)-slidewi)/dTi)+1;
        normalInspFlow = VTmean/median(Ti(criteriafornormal));
        % Scotty changed to 0.05, was 0.08
        thres1 = normalInspFlow*0.04; %%%%%%%%%%%%%%%%%%%%%%%% low amplitude % Mesa version had this changed to 0.2
        thres1B = normalInspFlow*0.9; %%%%%%%%%%%%%%%%%%%%%%%% low amplitude % Mesa version had this changed to 0.2
        thres2 = normalInspFlow*0.25; %%%%%%%%%%%%%%%%%%%%%%%% near flow = 0 (mean)
        
            Vflow2buffer = buffer(Vflow2,slidewi+1,slidewi+1-dTi,'nodelay');
                %handle zero in buffer in last col e.g. find(Vflow2buffer(:,end)==0)
            lastelementdelta = (size(Vflow2buffer,2)-1)*dTi + slidewi+1 - length(Vflow2);
            Vflow2buffer((end-lastelementdelta+1):end,end)=NaN;
            
            Nwindows = size(Vflow2buffer,2);
            Vflow_deltaB = zeros(Nwindows,1);
            Vflow_apnea = 0*Vflow;
            tic
            Vflow_delta = nanstd(Vflow2buffer);
            Vflow_median = nanmedian(Vflow2buffer);
            if ~minimum_figs
                Timebuffer = buffer(time,slidewi+1,slidewi+1-dTi,'nodelay');
                temptime = Timebuffer(1,:) + (slidewi+1)*dt/2;
            end
            
            X = [ones(slidewi+1,1) (1:slidewi+1)'];
            for i=1:Nwindows
                y = [Vflow2buffer(:,i)];
                if 0
                    %X=X(1:4:end,:);
                    y=y(1:4:end,:);
                end
                I = ~isnan(y);
                Vflow_deltaB(i)=std(y(I)-((X(I)\y(I))' * X(I)')');
                li = 1+(i-1)*dTi;
                if Vflow_deltaB(i)<thres1B&&Vflow_delta(i)<thres1&&(Vflow_median(i)-meanVflow)<thres2&&(Vflow_median(i)-meanVflow)>-thres2
                    Vflow_apnea((li+round(0*slidewi)):(li+round(1*slidewi))) = 1;
                end
            end
        
        Vflow_apnea = Vflow_apnea(1:length(Vflow));
        Vflow_apnea2 = Vflow_apnea;
        
        if ~minimum_figs
            ax1(1+1)=subplot(1+1,1,1+1); plot(time,Vflow_apnea2,time,Vflow_apnea0,time,Vflow_apnea,time,Vflow,time,Vflow2,time(BB_i_start),Vflow(BB_i_start),'r.',time(BB_i_mid),Vflow(BB_i_mid),'k.'); box('off');
            linkaxes(ax1,'x');
            hold on
            %temptime(i) = median(time(li:ri));
            plot(temptime,Vflow_delta,'g',temptime,Vflow_median,'b',temptime,thres1+0*temptime,'g:',temptime,thres2+0*temptime,'b:',temptime,-thres2+0*temptime,'b:')
            hold off
        end
        
    else
        Vflow_apnea = time*0;
        Vflow_apnea2 = time*0;
    end %if detect apneas
    
    %% Adjust Inspiratory start (BB_i_start) after end-expiratory pause
    if settings.modBB_i_start&&length(BB_i_start)>5
        Fpeakflow = 0.1;
        Fpeakvol = 0.2; %0.04;
        BB_i_start_new = NaN*BB_i_start;
        %find time to 50 percent of peak flow for all breaths to provide a normal value
        for i=1:length(BB_i_start)
            if Vflow_apnea2(BB_i_start(i))==1
                try
                    tempvol = (vol(BB_i_start(i):BB_i_mid(i))-vol(BB_i_start(i)));
                    tempvol = tempvol/tempvol(end);
                    new_start_vol_temp_delta = find(tempvol>Fpeakvol,1); %find where volume > threshold (Fpeakflow)
                    
                    
                    tempflow = Vflow(BB_i_start(i):BB_i_mid(i));
                    tempflow = tempflow/max(tempflow);
                    new_start_flow_temp_delta = find(tempflow>Fpeakflow,1); %find where volume > threshold (Fpeakflow)
                    
                    new_start_flowvol_temp_delta = min([new_start_vol_temp_delta new_start_flow_temp_delta]);
                    
                    %Only keep new start index if its within apnea2 range
                    new_start_silence_temp_delta = find(Vflow_apnea2(BB_i_start(i):BB_i_mid(i)) == 0, 1, 'first'); %is index of end silence
                    if isempty(new_start_silence_temp_delta)
                        continue
                    end
                    %                     if BB_i_start(i) + new_start_silence_temp_delta > BB_i_mid(i) %abort if algorithm classified entire breath as quiet and therefore fails
                    %                         continue
                    %                     end
                    %tempNewStartIdx = round(BB_i_mid(i)-(BB_i_mid(i)-BB_i_start(i)-new_start_silence_temp-1)/(1-0.25)); % had to do this to fix bug
                    %deltatemp = tempNewStartIdx - BB_i_start(i);
                    
                    if new_start_flowvol_temp_delta <= new_start_silence_temp_delta  %hit vol limit && tempNewStartIdx > 0
                        BB_i_start_new(i) = BB_i_start(i) + new_start_flowvol_temp_delta - 1;
                    else %end silence
                        BB_i_start_new(i) = BB_i_start(i) + new_start_silence_temp_delta - 1;
                        %                         BB_firstZero = new_start_silence_temp_delta - 1;
                        %                         BB_i_start_new(i) = round(BB_i_mid(i)-(BB_i_mid(i)-BB_i_start(i)-BB_firstZero-1));
                    end
                    
                catch me
                end
                
            end
        end
        
        nanIdx = ~isnan(BB_i_start_new);
        
        if ~minimum_figs
            ax1(1)=subplot(1+1,1,1); plot(time,vol,time(BB_i_start),vol(BB_i_start),'r.',time(BB_i_mid),vol(BB_i_mid),'k.'); box('off'); hold on
            plot(time(BB_i_start_new(nanIdx)),vol(BB_i_start_new(nanIdx)),'g*'); hold off
            ax1(1+1)=subplot(1+1,1,1+1); plot(temptime,Vflow_delta,'m',temptime,Vflow_median,'b',temptime,thres1+0*temptime,'m:',temptime,thres2+0*temptime,'b:',temptime,-thres2+0*temptime,'b:',temptime,3*Vflow_deltaB,'g',temptime,3*thres1B+0*temptime,'g:'); hold on
            plot(time,Vflow,time,Vflow2,time,Vflow_apnea2,time,Vflow_apnea0,time(BB_i_start),Vflow(BB_i_start),'rx',time(BB_i_mid),Vflow(BB_i_mid),'k.'); box('off');
            plot(time(BB_i_start_new(nanIdx)),Vflow(BB_i_start_new(nanIdx)),'g*'); hold off
            linkaxes(ax1,'x');
        end
        
        
        %         Ti_X_upper = prctile(BB_Ti_X,90);
        %
        %         Fearlyinspisapnea=0*BB_i_start;
        %         for i=1:length(BB_i_start)
        %             if isnan(BB_Ti_X(i))
        %                 continue
        %             end
        %             Fearlyinspisapnea(i)=mean(Vflow_apnea(BB_i_start(i):(BB_i_start(i)+round(BB_Ti_X(i)/dt))));
        %             if BB_Ti_X(i)>Ti_X_upper&&Fearlyinspisapnea(i)>0.5 %at least half of early insp is apnea
        %                 BB_i_start(i)=BB_i_start(i)+round((BB_Ti_X(i)-Ti_X_upper)/dt);
        %                 if i>1
        %                     BB_i_end(i-1)=BB_i_start(i);
        %                 end
        %             end
        %         end
        %
        %         % recalculate breath info:
        %
        %         if ~minimum_figs
        %             ax1(1)=subplot(1+1,1,1); plot(time,vol,time(BB_i_start),vol(BB_i_start),'r.',time(BB_i_mid),vol(BB_i_mid),'k.'); box('off');
        %             ax1(1+1)=subplot(1+1,1,1+1); plot(time,Vflow_apnea0,time,Vflow_apnea,time,Vflow,time,Vflow,time(BB_i_start),Vflow(BB_i_start),'r.',time(BB_i_mid),Vflow(BB_i_mid),'k.'); box('off');
        %             linkaxes(ax1,'x');
        %         end
        %% Update
        BB_i_start(nanIdx) = BB_i_start_new(nanIdx);
        BB_i_end(1:end-1) = BB_i_start(2:end);
        %         for i=1:length(BB_i_start)
        %             [~,tempi]=max(vol(BB_i_start(i):BB_i_end(i)-1));
        %             BB_i_mid(i,1)=BB_i_start(i)+tempi-1;
        %         end
        %
        if ~minimum_figs
            ax1(1)=subplot(1+1,1,1); plot(time,vol,time(BB_i_start),vol(BB_i_start),'r.',time(BB_i_mid),vol(BB_i_mid),'k.'); box('off');
            ax1(1+1)=subplot(1+1,1,1+1); plot(time,Vflow_apnea0,time,Vflow_apnea,time,Vflow,time,Vflow,time(BB_i_start),Vflow(BB_i_start),'r.',time(BB_i_mid),Vflow(BB_i_mid),'k.'); box('off');
            linkaxes(ax1,'x');
        end
        VTi = vol(BB_i_mid) - vol(BB_i_start);
        VTe = vol(BB_i_mid) - vol(BB_i_end);
        Ti = (BB_i_mid-BB_i_start)'*dt;
        Te = (BB_i_end-BB_i_mid)'*dt;
        Ttot = Ti+Te;
        VT = (VTi.*Te+VTe.*Ti)./(Ttot);
        VTmean = sum(VT.*(Ttot))/sum(Ttot);
        criteriafornormal = VT>(0.2*VTmean)&Ttot<10;
        normalTtot=median(Ttot(criteriafornormal));
        Vdot = VT./Ttot;
        
    end
    %%
    
    
    %% [OFF] Mid expiratory pause: If we find a small inspiration (with larger expiration) above FRC presume the breath is part of one continuous inspiration:
    %need to move to after apneas -- at least until after start insp is moved to end apnea.
    %VolStartInsp = (VTi-VTe)/VTmean;
    if 0 %causing trouble
        FVTe = VTe./VTi;
        
        VTe_prev = [NaN VTe(1:end-1)];
        VTi_prev = [NaN VTi(1:end-1)];
        
        FVT_prev = VTe_prev./VTi_prev;
        
        FVT_merged = (VTe+VTe_prev)./(VTi+VTi_prev);
        
        criteria=find(((FVT_merged<1.5)&(FVT_merged>0.67)&(FVTe>2)&(FVT_prev<0.67)&Te<(normalTe*0.75))); %VolStartInsp>0.33 %was leaving double peaks on smaller breaths with incl exhalation
        
        
        if ~isempty(criteria)
            firstbreathisbad=criteria(1)==1;
            while firstbreathisbad
                criteria(1)=[];
                BB_i_start(1)=[];
                BB_i_mid(1)=[];
                BB_i_end(1)=[];
                criteria=criteria-1;
                if isempty(criteria)
                    break
                end
                firstbreathisbad=criteria(1)==1;
            end
            criteria2=criteria-1;
            BB_i_start(criteria)=[];
            BB_i_mid(criteria)=[];
            BB_i_end(criteria2)=[];
        end
        
        % recalculate breath info:
        
        for i=1:length(BB_i_start)
            [~,tempi]=max(vol(BB_i_start(i):BB_i_end(i)));
            BB_i_mid(i,1)=BB_i_start(i)+tempi-1;
        end
        
        if ~minimum_figs
            
            ax1(1)=subplot(1+1,1,1); plot(time,vol,time(BB_i_start),vol(BB_i_start),'r.',time(BB_i_mid),vol(BB_i_mid),'k.'); box('off');
            ax1(1+1)=subplot(1+1,1,1+1); plot(time,Vflow,time(BB_i_start),Vflow(BB_i_start),'r.',time(BB_i_mid),Vflow(BB_i_mid),'k.'); box('off');
            linkaxes(ax1,'x');
        end
        VTi = vol(BB_i_mid) - vol(BB_i_start);
        VTe = vol(BB_i_mid) - vol(BB_i_end);
        Ti = (BB_i_mid-BB_i_start)'*dt;
        Te = (BB_i_end-BB_i_mid)'*dt;
        Ttot = Ti+Te;
        VT = (VTi.*Te+VTe.*Ti)./(Ttot);
        VTmean = sum(VT.*(Ttot))/sum(Ttot);
        criteriafornormal = VT>(0.2*VTmean)&Ttot<10;
        normalTtot=median(Ttot(criteriafornormal));
        if sum(criteriafornormal)<0.25*length(Ttot);
            normalTtot=median(Ttot);
        end
    end
    
    %% Long expirations are apneas (later broken up into smaller pieces) (0.000067 seconds unused)
    %But they should also contain detected apneas (Vflow_apnea0), otherwise slow large expirations get caught up
    
    FtotIsAnApnea = 1.0; % original setting was 1.0
    FteIsAnApnea = 0.5; % original setting was 0.5;
    TeThresholdForApnea = (FteIsAnApnea*normalTe+FtotIsAnApnea*normalTtot);
    Apnea_B = Te'>TeThresholdForApnea;
    Iapnea = find(Apnea_B);
    for i=length(Iapnea):-1:1
        Te_total = dt*(BB_i_end(Iapnea(i))-BB_i_mid(Iapnea(i)));
        Fapnea = mean(Vflow_apnea0(BB_i_mid(Iapnea(i)):BB_i_end(Iapnea(i))));
        if Fapnea==0 %no apnea detected
            Apnea_B(Iapnea(i))=0;
            continue
        end
        Te_est = normalTe;
        if (Te_total-Te_est)>FtotIsAnApnea*normalTtot     %if the new period is now too short to be called apnea, do not proceed
            BB_i_start_newapneabreath = BB_i_mid(Iapnea(i))+round(Te_est/dt);
            di=BB_i_end(Iapnea(i))-BB_i_start_newapneabreath;
            BB_i_mid_newapneabreath = BB_i_start_newapneabreath+round(normalTi/normalTtot*di);
            
            BB_i_start = [BB_i_start(1:Iapnea(i));BB_i_start_newapneabreath;BB_i_start(Iapnea(i)+1:end)];
            BB_i_end = [BB_i_end(1:Iapnea(i)-1);BB_i_start_newapneabreath;BB_i_end(Iapnea(i):end)];
            BB_i_mid = [BB_i_mid(1:Iapnea(i));BB_i_mid_newapneabreath;BB_i_mid(Iapnea(i)+1:end)];
            Apnea_B = [Apnea_B(1:Iapnea(i)-1);0;1;Apnea_B(Iapnea(i)+1:end)];
            
            if ~minimum_figs
                ax1(1)=subplot(2,1,1); plot(time,vol,...
                    time(BB_i_start),vol(BB_i_start),'r.',...
                    time(BB_i_mid),vol(BB_i_mid),'k.',...
                    time(BB_i_start_newapneabreath),vol(BB_i_start_newapneabreath),'ro',...
                    time(BB_i_mid_newapneabreath),vol(BB_i_mid_newapneabreath),'ko');
                box('off');
                ax1(1+1)=subplot(2,1,2);
                plot(time,Vflow_apnea0,...
                    time,Vflow_apnea,...
                    time,Vflow,...
                    time(BB_i_start),Vflow(BB_i_start),'r.',...
                    time(BB_i_mid),Vflow(BB_i_mid),'k.');
                box('off');
                linkaxes(ax1,'x');
            end
            
        else
            Apnea_B(Iapnea(i))=0;
        end
    end
    
    % recalculate breath info:
    
    % if ~minimum_figs
    %     ax1(1)=subplot(1+1,1,1); plot(time,vol,time(BB_i_start),vol(BB_i_start),'r.',time(BB_i_mid),vol(BB_i_mid),'k.');
    %     ax1(1+1)=subplot(1+1,1,1+1); plot(time,Vflow,time(BB_i_start),Vflow(BB_i_start),'r.',time(BB_i_mid),Vflow(BB_i_mid),'k.');
    %     linkaxes(ax1,'x');
    % end
    
    Ti = (BB_i_mid-BB_i_start)'*dt;
    Te = (BB_i_end-BB_i_mid)'*dt;
    Ttot = Ti+Te;
    
    
    %% Divide up apneas into Ttot-sized parcels (0.000086 seconds unused)
    Iapnea = find(Apnea_B);
    Nbreathsapnea=floor(Ttot(Iapnea)/normalTtot);
    for i=length(Iapnea):-1:1
        if Nbreathsapnea(i)>0
            di=round((BB_i_end(Iapnea(i))-BB_i_start(Iapnea(i)))/Nbreathsapnea(i));
            %     time(BB_i_start(Iapnea(i)))
            temp = (BB_i_start(Iapnea(i)):di:((BB_i_end(Iapnea(i)))+di/4))'; %di/5 just handles rounding issues
            BB_i_start = [BB_i_start(1:Iapnea(i));temp(2:(end-1));BB_i_start((Iapnea(i)+1):end)];
            BB_i_end = [BB_i_end(1:Iapnea(i)-1);temp(2:(end-1));BB_i_end(Iapnea(i):end)];
            BB_i_mid = [BB_i_mid(1:Iapnea(i)-1);temp(1:(end-1))+round((normalTi/normalTtot*di));BB_i_mid((Iapnea(i)+1):end)];
            Apnea_B = [Apnea_B(1:Iapnea(i));1+0*temp(2:(end-1));Apnea_B((Iapnea(i)+1):end)];
        end
    end
    
    % recalculate breath info:
    
    if ~minimum_figs
        ax1(1)=subplot(2,1,1);
        plot(time,vol,time(BB_i_start),vol(BB_i_start),'r.',...
            time(BB_i_mid),vol(BB_i_mid),'k.');
        box('off');
        ax1(1+1)=subplot(2,1,2);
        plot(time,Vflow,...
            time,Vflow_apnea0,...
            time,Vflow_apnea,...
            time(BB_i_start),Vflow(BB_i_start),'r.',...
            time(BB_i_mid),Vflow(BB_i_mid),'k.');
        hold('on'); stairs(time(BB_i_start),Apnea_B*max(Vflow)/2,'b');  hold('off');
        box('off');
        linkaxes(ax1,'x');
    end
    VTi = vol(BB_i_mid) - vol(BB_i_start);
    VTe = vol(BB_i_mid) - vol(BB_i_end);
    Ti = (BB_i_mid-BB_i_start)'*dt;
    Te = (BB_i_end-BB_i_mid)'*dt;
    Ttot = Ti+Te;
    VT = (VTi.*Te+VTe.*Ti)./(Ttot);
    VTmean = sum(VT.*(Ttot))/sum(Ttot);
    criteriafornormal = VT>(0.2*VTmean)&Ttot<10;
    normalTtot=median(Ttot(criteriafornormal));
    
    %% Load previously computed breath start and stop times
    % Dan Vena added this to avoid having to re-select manually edited
    % breath start and stop times
    if loopnumber == 2 && (isfield(settings, 'LoadPreviousBreathStartStop') && settings.LoadPreviousBreathStartStop)
        filename = [settings.AnalyzedDirectory,settings.savename,'_',num2str(n)];
        try
            load(filename,'BreathDataTable')

            if iscell(BreathDataTable)
                BrTable = BreathDataTable{winNum};
            else
                BrTable = BreathDataTable;
            end

            BB_i_start = BrTable.BB_i_start;
            BB_i_mid = BrTable.BB_i_mid;
            BB_i_end = BrTable.BB_i_end;
            Apnea_B = BrTable.ApneaB;

             % Recalculate breath info:
            VTi = vol(BB_i_mid) - vol(BB_i_start);
            VTe = vol(BB_i_mid) - vol(BB_i_end);
            Ti = (BB_i_mid-BB_i_start)'*dt;
            Te = (BB_i_end-BB_i_mid)'*dt;
            Ttot = Ti+Te;
            VT = (VTi.*Te+VTe.*Ti)./(Ttot);
            VTmean = sum(VT.*(Ttot))/sum(Ttot);
            criteriafornormal = VT>(0.2*VTmean)&Ttot<10;
            normalTtot=median(Ttot(criteriafornormal));
        catch me
            disp('Could not load previous breath start/stop times')
        end
    end
    
    %% Manually check and correct breath start and stop times 
    % Dan Vena added this for the DISE data which has few, but critical
    % breaths
    if loopnumber == 2 && (isfield(settings, 'ManuallyEditBreathTimes') && settings.ManuallyEditBreathTimes)
        Snore = evalin('caller','Snore');
        SnoreDB = evalin('caller','SnoreDB');
        BB_i_start_old = BB_i_start;

        GUIvars = struct('time',time','Vflow',Vflow','vol',vol',...
            'BB_i_start',BB_i_start,'BB_i_mid',BB_i_mid,'BB_i_end',BB_i_end,...
            'Snore',Snore,'SnoreDB',SnoreDB);
        hGui = BBStartStopAgain(GUIvars);
        % get the handle to the GUI and wait for it to be closed
        waitfor(hGui);
        disp('closed')
        % once GUI is closed, get updated version of GUIvars
        GUIVarsUpdate = get(0,'userdata');
        % continue with script
        
        BB_i_start = GUIVarsUpdate.BB_i_start;
        BB_i_mid = GUIVarsUpdate.BB_i_mid;
        BB_i_end = GUIVarsUpdate.BB_i_end;
        
        % fix Apnea_B
        Apnea_B_old = Apnea_B;
        Apnea_B = nan(size(BB_i_start));
        
        if length(BB_i_start) >= length(BB_i_start_old)
            same = ismember(BB_i_start,BB_i_start_old);
            Apnea_B(~same) = 0;
            sameold = ismember(BB_i_start_old,BB_i_start(same));
            Apnea_B(same) = Apnea_B_old(sameold);
        else
            same = ismember(BB_i_start_old,BB_i_start);
            samenew = ismember(BB_i_start,BB_i_start_old(same));
            Apnea_B(samenew) = Apnea_B_old(same);
            Apnea_B(~samenew) = 0;
        end
        
        if ~minimum_figs
            ax1(1)=subplot(2,1,1);
            plot(time,vol,time(BB_i_start),vol(BB_i_start),'r.',...
                time(BB_i_mid),vol(BB_i_mid),'k.');
            box('off');
            ax1(1+1)=subplot(2,1,2);
            plot(time,Vflow,...
                time,Vflow_apnea0,...
                time,Vflow_apnea,...
                time(BB_i_start),Vflow(BB_i_start),'r.',...
                time(BB_i_mid),Vflow(BB_i_mid),'k.');

            box('off');
            linkaxes(ax1,'x');
        end

        % Recalculate breath info:
        VTi = vol(BB_i_mid) - vol(BB_i_start);
        VTe = vol(BB_i_mid) - vol(BB_i_end);
        Ti = (BB_i_mid-BB_i_start)'*dt;
        Te = (BB_i_end-BB_i_mid)'*dt;
        Ttot = Ti+Te;
        VT = (VTi.*Te+VTe.*Ti)./(Ttot);
        VTmean = sum(VT.*(Ttot))/sum(Ttot);
        criteriafornormal = VT>(0.2*VTmean)&Ttot<10;
        normalTtot=median(Ttot(criteriafornormal));
    end
    
    %% Breaths with detected apneas have zero ventilation
    VT(Apnea_B==1)=0;
    if 0
        VTi(Apnea_B==1)=0;
        VTe(Apnea_B==1)=0;
    end
    
    %% Estimate leak (0.016668 seconds)
    %estimated residual leak
    estimateresidualleak=1;
    %estimateresidualleakfigs=1-minimum_figs;
    Nbreaths = length(VT); minbreaths = 5;
    if estimateresidualleak
        prclow = 10;
        if Nbreaths*prclow/100<minbreaths
            prclow = minbreaths/Nbreaths;
        end
        VTXtile = prctile(VT,prclow);
        prclow2 = 67;
        VTXtile2L = prctile(VT,prclow2);
        if VTXtile==0
            VTXtile=min(VT(VT>0));
        end
        if VTXtile2L==0
            VTXtile2L=min(VT(VT>0));
        end
        VTXtile2 = max(VT); %75
        VTbelowVTXtile=VT<=VTXtile;
        VTwithinVTXtile2=VT>=VTXtile2L&VT<=VTXtile2;
        
        VTlowestXtile = median(VT(VTbelowVTXtile));
        VTsecondlowestXtile = median(VT(VTwithinVTXtile2));
        
        iqrflow = prctile(Vflow,75)-prctile(Vflow,25);
        
        leak_res=iqrflow*[-3:0.05:3];
        FVTilowestXtile=0*leak_res;
        FVTisecondlowestXtile=FVTilowestXtile;
        projectedFVTizeroVT=FVTilowestXtile;
        
        for i=1:length(leak_res)
            %add leak back to tidal volumes
            VTi_lc = VTi + (leak+leak_res(i)).*Ti;
            VTe_lc = VTe - (leak+leak_res(i)).*Te;
            FVTi = (VTi_lc-VTe_lc)./Ttot;
            FVTilowestXtile(i) = median(FVTi(VTbelowVTXtile));
            FVTisecondlowestXtile(i) = median(FVTi(VTwithinVTXtile2));
            %yintercept
            projectedFVTizeroVT(i) = FVTilowestXtile(i)-VTXtile*(FVTisecondlowestXtile(i)-FVTilowestXtile(i))/(VTXtile2-VTXtile);
        end
        error = projectedFVTizeroVT;
        
        if 0
            figure(100)
            plot(leak_res,error,'.--');
        end
        
        for temp=1:1
            
            xintercept = interp1(error,leak_res,0,'linear','extrap');
            i = find(error>0,1);
            leak_res = xintercept;
            
            %from above for loop
            VTi_lc = VTi + (leak+leak_res).*Ti;
            VTe_lc = VTe - (leak+leak_res).*Te;
            FVTi = (VTi_lc-VTe_lc)./Ttot;
            
            FVTilowestXtile = median(FVTi(VTbelowVTXtile));
            FVTisecondlowestXtile = median(FVTi(VTwithinVTXtile2));
            projectedFVTizeroVT = FVTilowestXtile-VTXtile*(FVTisecondlowestXtile-FVTilowestXtile)/(VTXtile2-VTXtile); % output these?
            
            %error = abs(FVTilowestXtile)+abs(FVTisecondlowestXtile);
            error = projectedFVTizeroVT; %this should work for calculating leak, but not for minimising IE error effects
            
            if 0
                figure(100)
                plot(leak_res,error,'.--');
            end
        end
        
        if (~minimum_figs)&&0
            figure(12); set(gcf,'color',[1 1 1]);
            subplot(1,1,1);
            plot(VT,FVTi,'k.','markersize',15); hold('on'); box('off');
            xlabel('VT');ylabel('FVTi');
            plot([VTlowestXtile VTsecondlowestXtile],[FVTilowestXtile FVTisecondlowestXtile],'ko-'); box('off');
            hold('off');
        end
        
        leaktotal = -leak_res;
        leak_residual = leaktotal-leak;
        
        if 0
            figure(101)
            plot(time,Vflow+leak,time,0*Vflow+leaktotal,'--')
        end
        
    else
        leaktotal=NaN;
    end
    
    %% Additional parameters
    %time,VFlow_filtered1,BB_i_start,BB_i_mid,BB_i_end,BB_t,VI,VE,Ttot,leak,IEratio
    BB_t = time(BB_i_start);
    % Force V positive, just in case.
    if 1
        VT(VT<0)=0;
    end
    VE = VT./Ttot; %force these to be the same for now until LGfromFlow is updated
    VI = VT./Ttot;
    
    Vpeak = NaN*BB_i_start;
    Vpeakmean = NaN*BB_i_start;
    for i=1:length(BB_i_start)
        Vpeak(i) = max(Vflow(BB_i_start(i):BB_i_mid(i)));
        Vpeakmean(i) = mean(Vflow(BB_i_start(i):BB_i_mid(i)));
    end
    
    
    %% Find local leak using N=3 breaths
    if loopnumber==1
        leakBB=NaN*BB_i_start;
        maxleak=VTmean/normalTtot*0.1; % Scotty changed to 0.1, was 0.25
        leaksignal=0*time;
        for i=2:length(BB_i_start)-1
            leakBB(i) = (vol(BB_i_start(i+1))-vol(BB_i_start(i-1)))/(time(BB_i_start(i+1))-time(BB_i_start(i-1)));
            if leakBB(i)>maxleak
                rangei=BB_i_start(i):BB_i_end(i);
                leaksignal(rangei)=leakBB(i)-maxleak;
            elseif leakBB(i)<-maxleak
                rangei=BB_i_start(i):BB_i_end(i);
                leaksignal(rangei)=leakBB(i)+maxleak;
            end
        end
        leakBB(1)=leakBB(2);
        leakBB(length(BB_i_start))=leakBB(length(BB_i_start)-1);
        
        filter_HFcutoff_butter0 = 1;
        filter_order0 = 1;
        [B_butter0,A_butter0] = butter(filter_order0,[filter_HFcutoff_butter0]/(1/dt/2),'low');
        leaksignalfilt = filtfilt(B_butter0,A_butter0,leaksignal);
    end
    
end
%% rebuild leak-removed (transformed) flow from scratch without smoothing
% the scaling exponent here is set in Analysis.m depending on flow signal
Vflow_out = Vflow_backup-leak;
exponent = settings.scalingexponent;
Vflow_out(Vflow_out>0)=(Vflow_out(Vflow_out>0).^(exponent))/(IEratio^0.5);
Vflow_out(Vflow_out<0)=(-(-Vflow_out(Vflow_out<0)).^(exponent))*(IEratio^0.5);
Vflow_out = Vflow_out-leaksignalfilt(:);
if ~minimum_figs
    if 0
        figure(54); clf(figure(54));
        plot(time,Vflow_backup,'g', time,Vflow,'r',time,Vflow_out,'k');
    end
end
%Vflow_out
Vflow=Vflow(:);

%%
leak_B = leaksignalfilt(BB_i_mid);

%% NaN tracker
    Vflow_out(nantracker)=NaN;
    Vflow(nantracker)=NaN;
    nantracker_B = 0*BB_i_start;
    for i=1:length(BB_i_start)
        irange = BB_i_start(i):(BB_i_end(i)-1);
        nantracker_B(i) = max(nantracker(irange));
    end
    BB_i_start(nantracker_B==1)=[];
    BB_i_mid(nantracker_B==1)=[];
    BB_i_end(nantracker_B==1)=[];
    BB_t(nantracker_B==1)=[];
    VI(nantracker_B==1)=[];
    VE(nantracker_B==1)=[];
    Ttot(nantracker_B==1)=[];
    VT(nantracker_B==1)=[];
    Vpeak(nantracker_B==1)=[];
    Vpeakmean(nantracker_B==1)=[];
    Apnea_B(nantracker_B==1)=[];
    leak_B(nantracker_B==1)=[];

%%
%[VFlow_crossing_i; VFlow_crossingE_i]
% If less than 4 breaths detected in the period, set as artefact period.

% DLM added a special option here to allow windows with < 4 breaths.
% This is to suit some MRI data that only has very short respiratory traces
% This option, if not manually invoked as a special setting should have no 
% effect on usual processing.
if length(BB_i_start)<4 && ~(isfield(settings, 'AllowLessThan4Breaths') && settings.AllowLessThan4Breaths==1)
    time=NaN;
    Vflow_out=NaN;
    BB_i_start=NaN;
    BB_i_mid=NaN;
    BB_i_end=NaN;
    BB_t=NaN;
    VI=NaN;
    VE=NaN;
    Ttot=NaN;
    leak=NaN;
    IEratio=NaN;
    VT=NaN;
    Vpeak=NaN;
    Vpeakmean=NaN;
    Apnea_B=NaN;
    leak_B=NaN;
    
    bSigT = table();
    bSigT.time=time(:);
    bSigT.Vflow=Vflow(:);
    bSigT.Vflow_out=Vflow_out(:);
    bT = table();
    bT.BB_i_start=BB_i_start(:);
    bT.BB_i_mid=BB_i_mid(:);
    bT.BB_i_end=BB_i_end(:);
    bT.BB_t=BB_t(:);
    bT.VI=VI(:);
    bT.VE=VE(:);
    bT.Ttot=Ttot(:);
    bT.leak=leak(:) + 0*Ttot(:);
    bT.IEratio=IEratio(:) + 0*Ttot(:);
    bT.VT=VT(:);
    bT.Vpeak=Vpeak(:);
    bT.Vpeakmean=Vpeakmean(:);
    bT.Apnea_B=Apnea_B(:);
    bT.leak_B=leak_B(:);

    return;   
end

%%
% time,Vflow,BB_i_start,BB_i_mid,BB_i_end,BB_t,VI,VE,Ttot,leak,...
%    IEratio,VT,Vpeak,Vpeakmean,Apnea_B,Vflow_out,VTi,VTe,Ti,Te,leak_B,IEratioEstimated

bSigT = table();
bSigT.time=time(:);
bSigT.Vflow=Vflow(:);
bSigT.Vflow_out=Vflow_out(:);
bT = table();
bT.BB_i_start=BB_i_start(:);
bT.BB_i_mid=BB_i_mid(:);
bT.BB_i_end=BB_i_end(:);
bT.BB_t=BB_t(:);
bT.VI=VI(:);
bT.VE=VE(:);
bT.Ttot=Ttot(:);
bT.leak=leak(:) + 0*Ttot(:);
bT.IEratio=IEratio(:) + 0*Ttot(:);
bT.VT=VT(:);
bT.Vpeak=Vpeak(:);
bT.Vpeakmean=Vpeakmean(:);
bT.Apnea_B=Apnea_B(:);
bT.leak_B=leak_B(:);















