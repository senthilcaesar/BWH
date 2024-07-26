function [I,Vdot,VT,Ti,Te,VTi,VTe,VflowC_]=Vflowanalysis2CMM(Vflow,time,dt,minimum_figs,useactualzero,Nloops)
warning('off');
%% Respiratory Trace Analyse Using Volume as the primary method of breath ID
%[Pxy,F] = CPSD(X,Y,WINDOW,NOVERLAP,NFFT,Fs)
% Respiratory Trace Analyse
nicefig = @() set(gca,'box','off','tickdir','out','xcolor',[1 1 1],'xticklabel',[],'fontsize',8,'fontname','arial narrow');
nicefiglast = @() set(gca,'box','off','tickdir','out','fontsize',8,'fontname','arial narrow');

if useactualzero
    Nloops=1;
end

%transpose if needed:
if size(Vflow,2)<size(Vflow,1)
    Vflow=Vflow';
end
%transpose if needed:
if size(time,2)<size(time,1)
    time=time';
end

vol1=cumsum(Vflow)*dt;
leak = mean(Vflow);
if ~useactualzero
    voldetrend1 = leak*(time-time(1));
else
    leak = mean(Vflow);
    voldetrend1 = 0*(time-time(1));
end
vol2 = vol1-voldetrend1;


if ~minimum_figs
    figure(5-Nloops); set(gcf,'color',[1 1 1]); title('detect start and end of breaths');
    ax1(1)=subplot(Nloops+1,1,1); plot(time,vol1-voldetrend1); ylabel('vol');
    ax1(2)=subplot(Nloops+1,1,2); plot(time,vol2); ylabel('vol');
    linkaxes(ax1,'x');
end

%Filter volume gently for use in analysis
if 0
    filter_HFcutoff_butter1 = 10;
    filter_order1 = 2;
    [B_butterHcut,A_butterHcut] = butter(filter_order1,[filter_HFcutoff_butter1]/(1/dt/2),'low');
else
    if ~useactualzero
        filter_HFcutoff_butter1 = 10;
        filter_LFcutoff_butter1 = 1/15;
        filter_order1 = 2;
        [B_butterHcut,A_butterHcut] = butter(filter_order1,[filter_LFcutoff_butter1 filter_HFcutoff_butter1]/(1/dt/2));
        vol2_filtered1 = filtfilt(B_butterHcut,A_butterHcut,vol2);
        voldetrend2 = vol2_filtered1 - vol2;
    else
        filter_HFcutoff_butter1 = 10;
        filter_order1 = 2;
        [B_butterHcut,A_butterHcut] = butter(filter_order1,[filter_HFcutoff_butter1]/(1/dt/2),'low');
        vol2_filtered1 = filtfilt(B_butterHcut,A_butterHcut,vol2);
        voldetrend2 = vol2_filtered1 - vol2;
    end
end

%Vflow_filtered1 = filtfilt(B_butterHcut,A_butterHcut,Vflow);
Vflow_filtered1 = [diff(vol2_filtered1)/dt]; Vflow_filtered1=[Vflow_filtered1(1) Vflow_filtered1];
if ~minimum_figs
    ax1(2)=subplot(Nloops+1,1,2); plot(time,vol2_filtered1);
    linkaxes(ax1,'x');
end

%Aggressive low pass for apnea detection
if 1 %low pass filter signal
    filter_HFcutoff_butter0 = 1;
    filter_order0 = 3;
    [B_butter0,A_butter0] = butter(filter_order0,[filter_HFcutoff_butter0]/(1/dt/2),'low');
    Vflow2 = filtfilt(B_butter0,A_butter0,Vflow);
    if 1
        filter_LFcutoff_butter0 = 1/10;
        filter_order0 = 3;
        [B_butter0,A_butter0] = butter(filter_order0,[filter_LFcutoff_butter0]/(1/dt/2),'high');
        Vflow2 = filter(B_butter0,A_butter0,Vflow2);
    end
end


%% loop -- estimate start insp values, then detrend vol, then re-estimate
%start insp times...

for XX=1:Nloops
    
    %% detecting peaks (0.063440 seconds)
    
    % accompanying filtered flow trace
    %Vflow = [diff(vol)/dt]; Vflow=[Vflow(1) Vflow];
    vol=vol2_filtered1; %nomenclature fix needed
    
    [min_list,max_list] = peakdet(-vol,0.001*std(vol)); %in this configuration, max_list(i,1) > minlist(i,1).
    
    BB_i_start=min_list(:,1);
    BB_i_mid=max_list(:,1);
    if BB_i_start(1)==1;
        BB_i_start(1)=[];
    end
    if BB_i_start(end)==length(vol);
        BB_i_start(end)=[];
    end
    if BB_i_start(1)>BB_i_mid(1)
        BB_i_mid(1)=[];
    end
    while BB_i_mid(end)>BB_i_start(end)
        BB_i_mid(end)=[];
    end
    
    if ~minimum_figs
        
        ax1(1)=subplot(1+1,1,1); plot(time,vol,time(BB_i_start),vol(BB_i_start),'r.-',time(BB_i_mid),vol(BB_i_mid),'k.-');
        ax1(1+1)=subplot(1+1,1,1+1); plot(time,Vflow,time(BB_i_start),Vflow(BB_i_start),'r.',time(BB_i_mid),Vflow(BB_i_mid),'k.');
    end
    
    %estimate normal tidal volume from the distance between upper and lower
    %vols.
    
    if ~minimum_figs
        
        ax1(1)=subplot(1+1,1,1); plot(time,vol,time(BB_i_start),vol(BB_i_start),'r.',time(BB_i_mid),vol(BB_i_mid),'k.');
        ax1(1+1)=subplot(1+1,1,1+1); plot(time,Vflow,time,Vflow,time(BB_i_start),Vflow(BB_i_start),'r.',time(BB_i_mid),Vflow(BB_i_mid),'k.');
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
    VT_thres = VTmean*0.08;
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
            ax1(1)=subplot(1+1,1,1); plot(time,vol,time(BB_i_start),vol(BB_i_start),'r.',time(BB_i_mid),vol(BB_i_mid),'k.');
            ax1(1+1)=subplot(1+1,1,1+1); plot(time,Vflow,time,Vflow,time(BB_i_start),Vflow(BB_i_start),'r.',time(BB_i_mid),Vflow(BB_i_mid),'k.');
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
        ax1(1)=subplot(1+1,1,1); plot(time,vol,time(BB_i_start),vol(BB_i_start),'r.',time(BB_i_mid),vol(BB_i_mid),'k.');
        ax1(1+1)=subplot(1+1,1,1+1); plot(time,Vflow,time,Vflow,time(BB_i_start),Vflow(BB_i_start),'r.',time(BB_i_mid),Vflow(BB_i_mid),'k.');
    end
    
    %% If we find a  short or small inspiration, presume it is part of the previous expiration (0.000073 seconds)
    VTi_thres = VTmean*0.04;
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
        
        ax1(1)=subplot(1+1,1,1); plot(time,vol,time(BB_i_start),vol(BB_i_start),'r.',time(BB_i_mid),vol(BB_i_mid),'k.');
        ax1(1+1)=subplot(1+1,1,1+1); plot(time,Vflow,time(BB_i_start),Vflow(BB_i_start),'r.',time(BB_i_mid),Vflow(BB_i_mid),'k.');
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
        
        ax1(1)=subplot(1+1,1,1); plot(time,vol,time(BB_i_start),vol(BB_i_start),'r.',time(BB_i_mid),vol(BB_i_mid),'k.');
        ax1(1+1)=subplot(1+1,1,1+1); plot(time,Vflow,time(BB_i_start),Vflow(BB_i_start),'r.',time(BB_i_mid),Vflow(BB_i_mid),'k.');
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
        dTi=50; %was 25, 50 halves the time...
        slidewi = round(slidew/dt);
        Nwindows = ceil((length(Vflow)-slidewi)/dTi)+1;
        normalInspFlow = VTi/Ti;
        thres1 = normalInspFlow/0.75; %%%%%%%%%%%%%%%%%%%%%%%% low amplitude
        thres2 = normalInspFlow*0.5; %%%%%%%%%%%%%%%%%%%%%%%% near flow = 0 (mean)
        Vflow_delta = zeros(Nwindows,1);
        Vflow_median = zeros(Nwindows,1);
        Vflow_apnea = 0*Vflow;
        
        for i=1:Nwindows %this loop is slow
            li = 1+(i-1)*dTi;
            ri = li+slidewi;
            if ri>length(Vflow),ri=length(Vflow); end
            Vflow_delta(i)=max(Vflow2(li:ri))-min(Vflow2(li:ri));
            Vflow_median(i)=median(Vflow2(li:ri));
            if Vflow_delta(i)<thres1&&(Vflow_median(i)-meanVflow)<thres2&&(Vflow_median(i)-meanVflow)>-thres2
                Vflow_apnea((li+round(0.10*slidewi)-1):(li+round(0.90*slidewi))) = 1;
            end
        end
        
        Vflow_apnea = Vflow_apnea(1:length(Vflow));
        Vflow_apnea0 = Vflow_apnea;
        
        if ~minimum_figs
            ax1(1+1)=subplot(1+1,1,1+1); plot(time,Vflow_apnea0,time,Vflow_apnea,time,Vflow,time,Vflow2,time(BB_i_start),Vflow(BB_i_start),'r.',time(BB_i_mid),Vflow(BB_i_mid),'k.');
            linkaxes(ax1,'x');
        end
    else
        Vflow_apnea = time*0;
        Vflow_apnea0 = time*0;
    end %if detect apneas
    
    %% In breaths with apneas - move the end-inspiration to the earliest maxima (0.003859 seconds)
    
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
    
    criteria = Fapnea>Fapneathres&Earlierpeak&Ttot'>normalTtot&CurrentPeakinApnea==1;
    
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
        ax1(1)=subplot(1+1,1,1); plot(time,vol,time(BB_i_start),vol(BB_i_start),'r.',time(BB_i_mid),vol(BB_i_mid),'k.');
        ax1(1+1)=subplot(1+1,1,1+1); plot(time,Vflow_apnea0,time,Vflow_apnea,time,Vflow,time,Vflow,time(BB_i_start),Vflow(BB_i_start),'r.',time(BB_i_mid),Vflow(BB_i_mid),'k.');
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
    
    
    %% Shallow breaths / inspirations occurring predominantly in detected apneas with low ventilation is an apnea -- delete / merge with prior expiration
    %0.003667 seconds
    
    Apnea_B1 = 0*BB_i_start;
    for i=1:length(BB_i_start)
        if mean(Vflow_apnea0(BB_i_start(i):BB_i_end(i)))>0.67&&VT(i)<0.2*VTmean
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
        ax1(1)=subplot(1+1,1,1); plot(time,vol,time(BB_i_start),vol(BB_i_start),'r.',time(BB_i_mid),vol(BB_i_mid),'k.');
        ax1(1+1)=subplot(1+1,1,1+1); plot(time,Vflow_apnea0,time,Vflow_apnea,time,Vflow,time,Vflow,time(BB_i_start),Vflow(BB_i_start),'r.',time(BB_i_mid),Vflow(BB_i_mid),'k.');
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
    
    %% Long expirations are apneas (later broken up into smaller pieces) (0.000067 seconds unused)
    %But they should also contain detected apneas (Vflow_apnea0), otherwise slow large expirations get caught up
    
    FtotIsAnApnea = 1.0;
    TeThresholdForApnea = (0.5*normalTe+FtotIsAnApnea*normalTtot);
    Apnea_B = Te'>TeThresholdForApnea;
    Iapnea = find(Apnea_B);
    for i=length(Iapnea):-1:1
        Te_total = dt*(BB_i_end(Iapnea(i))-BB_i_mid(Iapnea(i)));
        clear temp
        temp = find(Vflow_apnea0(BB_i_mid(Iapnea(i)):BB_i_end(Iapnea(i))),1);
        Te_est1 = dt*(temp-1); %potentially 1. placing end-exp at the start of apnea...
        if isempty(temp)
            Apnea_B(Iapnea(i))=0;
            continue
        end
        Te_est2 = normalTe; %potentially 2. placing end-exp at time of normal expiration
        Te_est = max([Te_est1 Te_est2]); %use max of 1 and 2 to avoid the possiblity of apnea starting at end insp.
        
        if (Te_total-Te_est)>FtotIsAnApnea*normalTtot     %if the new period is now too short to be called apnea, do not proceed
            BB_i_start_newapneabreath = BB_i_mid(Iapnea(i))+round(Te_est/dt);
            di=BB_i_end(Iapnea(i))-BB_i_start_newapneabreath;
            BB_i_mid_newapneabreath = BB_i_start_newapneabreath+round(normalTi/normalTtot*di);
            
            BB_i_start = [BB_i_start(1:Iapnea(i));BB_i_start_newapneabreath;BB_i_start(Iapnea(i)+1:end)];
            BB_i_end = [BB_i_end(1:Iapnea(i)-1);BB_i_start_newapneabreath;BB_i_end(Iapnea(i):end)];
            BB_i_mid = [BB_i_mid(1:Iapnea(i));BB_i_mid_newapneabreath;BB_i_mid(Iapnea(i)+1:end)];
            Apnea_B = [Apnea_B(1:Iapnea(i)-1);0;1;Apnea_B(Iapnea(i)+1:end)];
            
            if ~minimum_figs
                ax1(1)=subplot(1+1,1,1); plot(time,vol,time(BB_i_start),vol(BB_i_start),'r.',time(BB_i_mid),vol(BB_i_mid),'k.',time(BB_i_start_newapneabreath),vol(BB_i_start_newapneabreath),'ro',time(BB_i_mid_newapneabreath),vol(BB_i_mid_newapneabreath),'ko');
                ax1(1+1)=subplot(1+1,1,1+1); plot(time,Vflow_apnea0,time,Vflow_apnea,time,Vflow,time,Vflow,time(BB_i_start),Vflow(BB_i_start),'r.',time(BB_i_mid),Vflow(BB_i_mid),'k.');
                linkaxes(ax1,'x');
            end
        else
            Apnea_B(Iapnea(i))=0;
        end
    end
    
    
    %% recalculate breath info:
    
    % if ~minimum_figs
    %
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
    
    %% recalculate breath info:
    
    if ~minimum_figs
        ax1(1)=subplot(1+1,1,1); plot(time,vol,time(BB_i_start),vol(BB_i_start),'r.',time(BB_i_mid),vol(BB_i_mid),'k.');
        ax1(1+1)=subplot(1+1,1,1+1); plot(time,Vflow,time,Vflow,time,Vflow_apnea0,time,Vflow_apnea,time(BB_i_start),Vflow(BB_i_start),'r.',time(BB_i_mid),Vflow(BB_i_mid),'k.');
        hold('on'); stairs(time(BB_i_start),Apnea_B*max(Vflow)/2,'b');  hold('off');
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
    
    
    %% Additional information
    % baselinedriftstd = std(vollower+voldetrend1);
    % Endexpvol = vol2_filtered1(BB_i_start);
    % Endexpvol_t = time(BB_i_start);
    vollower = interp1(time(BB_i_start),vol2_filtered1(BB_i_start),time,'linear');
    vollower(1:BB_i_start(1)-1)=vol2_filtered1(BB_i_start(1)); %extrapolation using nearest
    vollower(BB_i_start(end):end)=vol2_filtered1(BB_i_start(end)); %extrapolation using nearest
    temp = polyfit(time,vollower,1);
    
    leak = temp(1) + leak;
    % if ~minimum_figs
    % figure(99);
    % plot(time,voldetrend1+vol2_filtered1,time(BB_i_start),voldetrend1(BB_i_start)+vol2_filtered1(BB_i_start),'r.',time(BB_i_mid),voldetrend1(BB_i_mid)+vol2_filtered1(BB_i_mid),'k.',time,voldetrend1+volupper,'k:',time,voldetrend1+vollower,'r:');
    % title(['leak=' num2str(leak) ', baselinedriftstd=' num2str(baselinedriftstd)]);
    % set(gcf,'position',[50 100 500 500]);
    % end
    
    % vollower_record = - voldetrend2;
    %baselinedriftstd = std(vollower_record);
    
    % figure(900);
    % plot(time,vol2_filtered1+vollower_record,time(BB_i_start),vol2_filtered1(BB_i_start),'r.',time(BB_i_mid),vol2_filtered1(BB_i_mid),'k.',time,vollower_record,'k:');
    % set(gcf,'position',[50 100 500 500]);
    VflowC_ = Vflow_filtered1';
    %% detrend and re-run
    if XX<Nloops
        %vol2_filtered1_backup = vol2_filtered1;
        %vol2_filtered1
        if 1
            volarray = vol2_filtered1([BB_i_start;BB_i_end(end)]);
            %volarray = [volarray(1) volarray volarray(end)];
            %medfilt1
            %volarray2 = medfilt1(volarray,1)
            %volarray2 = volarray2(2:end-1);
            voltrend = interp1(time([BB_i_start;BB_i_end(end)]),volarray,time,'linear');
            voltrend(1:BB_i_start(1)-1)=vol2_filtered1(BB_i_start(1));
            voltrend(BB_i_end(end):end)=vol2_filtered1(BB_i_end(end));
            %figure(102);plot(time,[vol2_filtered1;voltrend]);
            vol2_filtered1=vol2_filtered1-voltrend;
            %     vollower_record = vollower_record + voltrend;
        else
            leak2 = mean(Vflow(BB_i_start(1):BB_i_end(end)));
            voltrend = leak2*(time-time(1));
            vol2_filtered1 = filtfilt(B_butterHcut,A_butterHcut,vol2);
            vol2_filtered1=vol2_filtered1-voltrend;
        end
    end
    
end
%% Force V positive, just in case.

VI = VTi./Ttot;
VE = VTe./Ttot;
Vdot = VT./Ttot;

VT(VT<0)=0;
Vdot(Vdot<0)=0;
VTi(VTi<0)=0;
VTe(VTe<0)=0;
VE(VE<0)=0;
VI(VI<0)=0;

%% Start/End breath indices in original indices
clear I;
I.starti = BB_i_start;
I.midi = BB_i_mid;
I.endi = BB_i_end;


