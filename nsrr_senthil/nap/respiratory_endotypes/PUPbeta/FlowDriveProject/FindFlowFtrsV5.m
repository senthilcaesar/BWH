% FlowFtrsALi = [];
% AAFs=round(1/mean(diff(time)));             % Fs, calc original Fs.
% for ii=1:length(BB_i_start)
%     AAFlow=Vflow2(BB_i_start(ii):BB_i_end(ii))-Vflow2(BB_i_start(ii)); % get flow for this breath; offset flow to zero, set as first sample
%     AAftrs=FindFlowFtrsV5(AAFlow,AAFs,BB_i_mid(ii)-BB_i_start(ii)+1);
%     AAvariableNames = fieldnames(AAftrs);
%     for jj=1:length(AAvariableNames) % add these results to the growing struct
%         FlowFtrsALi=setfield(FlowFtrsALi,{1},AAvariableNames{jj},{ii,':'},getfield(AAftrs, {1}, AAvariableNames{jj}));
%     end
% end
% % convert the struct into a matrix
% FlowFeaturesAli=[]; %NaN(size(FlowFtrsALi.D1,1),length(fieldnames(FlowFtrsALi)));
% AAvariableNames = fieldnames(FlowFtrsALi);
% for jj=1:length(AAvariableNames)
%     FlowFeaturesAli=[FlowFeaturesAli  getfield(FlowFtrsALi, {1}, AAvariableNames{jj})];
% end
% 

function [ftrs]=FindFlowFtrsV5(Flow,Fs_in,BrMid)

% Ali's processing uses the entire breath flow
% offset flow to zero, first sample
% BrMid is just the BB_i_mid for this breath



ValidFlag=1;
ftrs=[];
Fs=125;
Flow=resample(Flow,Fs,Fs_in); % DLM - why resampling to 250 ??
% Br.BrMid=round(Fs*(Br.BrMid/Fs_in));
% Br.BrEnd=length(Flow);

FlowUnfilterd=Flow;

% BrMid=Br.BrMid;
I_01=round(0.1*BrMid);
I_09=round(0.9*BrMid);

Time=(0:length(Flow)-1)/Fs;

Flutter=FindWaveLetCoefs(Flow,'db4',Fs,'Flow Noise');
Flow1=Flow-Flutter;

PeakFindThr=0.008*max(Flow1(1:BrMid));
if PeakFindThr<0
    ValidFlag=0;
end
Flow=Flow1;

%% find Ti, Te, Ttot, VTi, VTe, VT
vol=cumtrapz(Time,Flow); %figure();plot(vol);
VTi = vol(BrMid) - vol(1);
VTe = vol(BrMid) - vol(end);
Ti = (BrMid-1)/Fs;
Te = (length(vol)-BrMid)/Fs;
VT = (VTi*Te+VTe*Ti)/(Ti + Te);
ftrs.Ti_Te = Ti/Te;
ftrs.Ti_Ttot = Ti/(Ti+Te);
ftrs.VTi_VTe = VTi/VTe;
ftrs.VTi_VT = VTi/VT;

%% find peak flow
[pks,locs,w,p]=findpeaks(Flow(1:BrMid),'MinPeakProminence',0.03,'Annotate','extents');
%figure();plot(Flow(1:BrMid)); hold on; plot(locs,pks,'ro');
if ~isempty(locs)
%     I1_chnge=locs(1);
    if locs(1)<0.5*BrMid
        I1_chnge=locs(1);
    else
        I1_chnge=round(0.5*BrMid);
    end
    if locs(end)>0.8*BrMid
        I2_chnge=locs(end);
    else
        I2_chnge=round(0.8*BrMid);
    end
else
    I1_chnge=I_01;
    I2_chnge=I_09; 
end
if ~isempty(p)
    WidthMostProminentPeak=w(p==max(p))/Fs/Ti;
    PromMostProminentPeak=p(p==max(p));
else
    WidthMostProminentPeak=NaN;
    PromMostProminentPeak=NaN;
end
if length(locs)>=3
    fp1_idx=locs(1)/Fs;
    fp2_idx=locs(2)/Fs;
    fp3_idx=locs(3)/Fs;
    fp1=pks(1);
    fp2=pks(2);
    fp3=pks(3);
    
    pkProm1=p(1)/PromMostProminentPeak;
    pkProm2=p(2)/PromMostProminentPeak;
    pkProm3=p(3)/PromMostProminentPeak;
    
    
elseif length(locs)==2
    fp1_idx=locs(1)/Fs;
    fp2_idx=locs(2)/Fs;
    fp3_idx=NaN;
    fp1=pks(1);
    fp2=pks(2);
    fp3=NaN;
    
    pkProm1=p(1)/PromMostProminentPeak;
    pkProm2=p(2)/PromMostProminentPeak;
    pkProm3=NaN;
    
elseif length(locs)==1
    fp1_idx=locs(1)/Fs;
    fp2_idx=NaN;
    fp3_idx=NaN;
    fp1=pks(1);
    fp2=NaN;
    fp3=NaN;
    
    pkProm1=p(1)/PromMostProminentPeak;
    pkProm2=NaN;
    pkProm3=NaN;
    
else
    fp1_idx=NaN;
    fp2_idx=NaN;
    fp3_idx=NaN;
    fp1=NaN;
    fp2=NaN;
    fp3=NaN;
    
    pkProm1=NaN;
    pkProm2=NaN;
    pkProm3=NaN;
    
end

Vpeak=max(pks);
peakFlow_idx=locs(pks==max(pks))/Fs;
if isempty(Vpeak)
    [Vpeak,peakFlow_idx]=max(Flow);
end
    
Vmean=mean(Flow(1:BrMid));
VmeanE=-mean(Flow(BrMid+1:end));
[ExpFlowPeak,ExpFlowPeakIdx]=min(Flow1); % original 
%[ExpFlowPeak,ExpFlowPeakIdx]=min(Flow);   % without flutter removal
PeakExpFlowTimeFromBrMid=(ExpFlowPeakIdx-BrMid)/Fs;
VpeakE=-ExpFlowPeak;
TpeakE=PeakExpFlowTimeFromBrMid;

ftrs.Vpeak1_Vpeak=fp1/Vpeak;
ftrs.Vpeak2_Vpeak=fp2/Vpeak;
ftrs.Vpeak3_Vpeak=fp3/Vpeak;
ftrs.Vpeak_VpeakE=Vpeak/VpeakE;
ftrs.Vmean_VmeanE=Vmean/VmeanE;
ftrs.Vmean_Vpeak=Vmean/Vpeak;
ftrs.Vmean_VT=Vmean/VT;
ftrs.Vpeak_VT=Vpeak/VT;
ftrs.VpeakE_VT=VpeakE/VT;
ftrs.Tpeak_Ti=peakFlow_idx/Ti;
ftrs.Tpeak1_Ti=fp1_idx/Ti;
ftrs.Tpeak2_Ti=fp2_idx/Ti;
ftrs.Tpeak3_Ti=fp3_idx/Ti;
ftrs.TpeakE_Te=TpeakE/Te;
ftrs.Tpeak_TpeakE=peakFlow_idx/TpeakE;
ftrs.MostProminentPeakWidth_Ti=WidthMostProminentPeak;
ftrs.pkProm1_n=pkProm1;
ftrs.pkProm2_n=pkProm2;
ftrs.pkProm3_n=pkProm3;

%% find NED and ratio of flow at different times
NED=NaN;
NEDSpan=0.25*(BrMid-1);
if ~isnan(pks) & ValidFlag==1
    NEDMinValIdx=round(peakFlow_idx*Fs):round(BrMid/2+NEDSpan);
    if ~isempty(NEDMinValIdx)
        NED=100*(Flow(round(peakFlow_idx*Fs))-min(Flow(NEDMinValIdx)))/Vpeak;
    else
        NED=NaN;
    end
end
if isempty(NED)
    ftrs.NED=NaN;
else
    ftrs.NED=NED;
end

%% find flow derivative and greatest change in the flow between 0.25%-0.75% from start point
if ValidFlag==1
    % Calculating smooth derivative
    FlowI1_I2=Flow(I1_chnge:I2_chnge);
    TimeI1_I2=Time(I1_chnge:I2_chnge)';
    tt=-(I1_chnge):length(Flow)-(I1_chnge)-1;    
    ipt=findchangepts(FlowI1_I2,'MinThreshold' ,0.005,'Statistic','linear');   
    ipt2=[1;ipt;length(FlowI1_I2)];
    for ii=2:length(ipt2)
        clear X b
        X = [ones(size(TimeI1_I2(ipt2(ii-1):ipt2(ii)-1))) TimeI1_I2(ipt2(ii-1):ipt2(ii)-1)];
        b = regress(FlowI1_I2(ipt2(ii-1):ipt2(ii)-1),X);    % Removes NaN data
        slp(ii-1)=b(2);
        intercept=b(1);
        stFlow(ii-1)=intercept+TimeI1_I2(ipt2(ii-1))*slp(ii-1);
        endFlow(ii-1)=intercept+TimeI1_I2(ipt2(ii)-1)*slp(ii-1);   
    end
    
    rise=endFlow-stFlow;
    Disc=[0 stFlow(2:end)-endFlow(1:end-1)];
    riseSlp=abs(rise.*slp);
    rise=abs(rise);
    Disc=abs(Disc);   
    FlowChangeSrtd=NaN(5,1);
    riseSrtd=NaN(5,1);
    DiscSrtd=NaN(5,1);
    RiseSlope=NaN(5,1);
    if ~isempty(slp)
        [FlowChangeSrtd_temp,srtd_idx]=sort(abs(slp),'descend');
        rise_tmp=rise(srtd_idx);
        Disc_tmp=Disc(srtd_idx);
        riseSlp_tmp=sort(riseSlp,'descend');
        if length(FlowChangeSrtd_temp)<5
            FlowChangeSrtd_temp(length(FlowChangeSrtd_temp):5)=min(FlowChangeSrtd_temp);
            rise_tmp(length(rise_tmp):5)=rise_tmp(end);
            Disc_tmp(length(Disc_tmp):5)=Disc_tmp(end);
            riseSlp_tmp(length(riseSlp_tmp):5)=riseSlp_tmp(end);
        end
        for ii=1:length(FlowChangeSrtd_temp)
            if ii<=5
                FlowChangeSrtd(ii)=FlowChangeSrtd_temp(ii);
                riseSrtd(ii)=rise_tmp(ii);
                DiscSrtd(ii)=Disc_tmp(ii);
                RiseSlope(ii)=riseSlp_tmp(ii);
            else
                break;
            end 
        end
    end    
    ftrs.D1=FlowChangeSrtd(1)/VT; 
    ftrs.D2=riseSrtd(1)/VT;   
    ftrs.D3=DiscSrtd(1)/VT;  
    ftrs.D4=RiseSlope(1)/VT^2;      
else
    ftrs.D1=NaN;      
    ftrs.D2=NaN;   
    ftrs.D3=NaN;
    ftrs.D4=NaN;
end

%% Find inspiratory and expiratory fluttering power
% obviously, this must use unfiltered flow (ie without flutter removed)
if ValidFlag==1 
    FlowInsp=FlowUnfilterd(1:BrMid); 
    FlowExp=FlowUnfilterd(BrMid:end);
    
%     N=length(FlowInsp);
%     w=hanning(N)/rms(hanning(N));
%     X = fft(w.*FlowInsp);
%     Pxx = (X.*X')/N.^2;
    
    [Pxx1,f1] = periodogram(FlowInsp,[],[],Fs);
    [Pxx2,f2] = periodogram(FlowExp,[],[],Fs);
%     [Pxx1,f1]=pwelch(FlowInsp,[],[],[],Fs);
%     [Pxx2,f2]=pwelch(FlowExp,[],[],[],Fs);
    Pow1=bandpower(Pxx1,f1,[5,round(Fs/2)-1],'psd');
    Pow2=bandpower(Pxx2,f2,[5,round(Fs/2)-1],'psd');
    ftrs.InspFlutPow_ExpFlutPow=Pow1/Pow2;
    
    ftrs.InspFlutPow_Vpeak2=Pow1/Vpeak^2;
    ftrs.ExpFlutPow_VpeakE2=Pow2/VpeakE^2;
    ftrs.InspFlutPow_VT2=Pow1/VT^2;
    ftrs.ExpFlutPow_VT2=Pow2/VT^2;
else
    ftrs.InspFlutPow_ExpFlutPow=NaN;
    ftrs.InspFlutPow_Vpeak2=NaN;
    ftrs.ExpFlutPow_VpeakE2=NaN;
    ftrs.InspFlutPow_VT2=NaN;
    ftrs.ExpFlutPow_VT2=NaN;
end


%% Find Expiratory Flow Limitation Index
%ExpFlow sampled at 250Hz 
EFLI=CalcAsymmetry(Flow,BrMid,Fs);
ftrs.EFLI=EFLI;


%% Find Insp. flattness index
if ValidFlag==1
    if 0
        FlowData=Flow(1:BrMid);
        start=[max(FlowData)*2]; % A reasonable estimate of a starting physiological value for each parameter
        lower=[0]; % Minimum conveivable physiological value for each parameter
        upper=[max(FlowData)*20]; % Maximum conveivable physiological value for each parameter

        OPTIONS = optimset('Display','off','TolX',1e-4,'TolFun',1e-4,'MaxIter',4000,'MaxFunEvals',150,'Algorithm','interior-point');
        Parameters=start;
        [c,Error,~,~] = fmincon(@(Parameters) TheInvertedParabola(Parameters,FlowData),Parameters,[],[],[],[],lower,upper,[],OPTIONS); %@(Parameters) model_nonlinear_constraints() is now [],[]

        [~,y_insp,Rsquared]=TheInvertedParabola(c,FlowData);
    %     figure(1);subplot(211);plot(FlowData);hold on;plot(y,'r');hold off
        MeanIsnpFlow=mean(FlowData);

        InspFlatnessIndex = trapz((0:length(FlowData)-1)/Fs,abs(FlowData-y_insp'));

        ftrs.InspFlatnessIndex=InspFlatnessIndex/(MeanIsnpFlow*BrMid/Fs);
    else
        % at the top, the following were defined:
        %I_01=round(0.1*BrMid);
        %I_09=round(0.9*BrMid);
        MeanIsnpFlow=mean(Flow(1:BrMid));
        InspFlowTemp=Flow(I_01:I_09);
        InspFlatnessIndex = trapz((0:length(InspFlowTemp)-1)/Fs,abs(InspFlowTemp-MeanIsnpFlow)); 
        
        ftrs.InspFlatnessIndex=InspFlatnessIndex/(MeanIsnpFlow*BrMid/Fs);
    end
else
    ftrs.InspFlatnessIndex=NaN;
end


%% Find Exp. flattness index
if ValidFlag==1
    if 0
        clear FlowData OPTIONS Parameters c y

        FlowData=Flow(BrMid+1:end);

        start=[0.2 0.2]; % A reasonable estimate of a starting physiological value for each parameter
        lower=[0 0]; % Minimum conveivable physiological value for each parameter
        upper=[1 1]; % Maximum conveivable physiological value for each parameter

        OPTIONS = optimset('Display','off','TolX',1e-4,'TolFun',1e-4,'MaxIter',4000,'MaxFunEvals',150,'Algorithm','interior-point');
        Parameters=start;
        [c,Error,~,~] = fmincon(@(Parameters) SingleBreathModelExp(Parameters, FlowData ,Fs),Parameters,[],[],[],[],lower,upper,[],OPTIONS); %@(Parameters) model_nonlinear_constraints() is now [],[]

        [~,y,Rsquared]=SingleBreathModelExp(c,FlowData,Fs);

        MeanExpFlow=mean(FlowData);

        ExpFlatnessIndex = trapz((0:length(FlowData)-1)/Fs,abs(FlowData-y'));

        ftrs.ExpFlatnessIndex=ExpFlatnessIndex/(-MeanExpFlow*length(FlowData)/Fs);

        y=[y_insp y];

        time=(0:length(Flow)-1)/Fs;
        figure(1);plot(time,Flow);hold on;plot(time, y,'r');hold off
        text(time(round(BrMid/2)),0.75,num2str(ftrs.InspFlatnessIndex));
        text(time(round((BrMid+length(Flow))/2)),0.75,num2str(ftrs.ExpFlatnessIndex));
        ylim([min(y)-0.1 max(0.8,[max(y)+0.1])])
    else    
        MeanExpFlow=mean(Flow(BrMid:end));
        ExpFlowDur=length(Flow)-BrMid;
        ExpFlowTemp=Flow(BrMid+round(0.2*ExpFlowDur):BrMid+round(0.8*ExpFlowDur));
        ExpFlatnessIndex = trapz((0:length(ExpFlowTemp)-1)/Fs,abs(ExpFlowTemp-MeanExpFlow));

        ftrs.ExpFlatnessIndex=ExpFlatnessIndex/(-MeanExpFlow*ExpFlowDur/Fs);
    end
else
    ftrs.ExpFlatnessIndex=NaN;
end
ftrs.InspFlatnessIndex_ExpFlatnessIndex=ftrs.InspFlatnessIndex/ftrs.ExpFlatnessIndex;

end