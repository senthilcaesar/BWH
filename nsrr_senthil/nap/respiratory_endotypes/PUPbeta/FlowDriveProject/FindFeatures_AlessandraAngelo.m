% %% AlessandraAngelo flow features
% % Flow is the current Flow variable, but could do leak adjust as per Ali's
% Fs = round(1/mean(diff(time)));
% TemplateVol = NaN; % Need to know how template volume is made
% FlowFeaturesAA=NaN(length(BB_i_start),26); % 26 cols of features, n rows of breaths
% for bb = 1:length(BB_i_start)
%     % setting these (part) blindly...
%     Br.Start=BB_i_start(bb);
%     Br.Mid=BB_i_mid(bb);%-BB_i_start(ii)+1;
%     Br.End=BB_i_end(bb);%-BB_i_start(ii)+1;
%     Flow2 = Vflow2(Br.Start:Br.End);
%     Br.VT = VT(bb);
%     Br.VT95 = prctile(Br.VT,95); % confirm this is correct (see Sara's code)
%     [FlowFeaturesAA(bb,:),ftrsNames]=FindFeatures_AlessandraAngelo(Flow2,Fs,Br,TemplateVol);
% end

function [ftrs,ftrsNames]=FindFeatures_AlessandraAngelo(Flow,Fs,Br,TemplateVol)
showfigures=0;

ftrs=zeros(1,26);
ftrsNames=cell(1,26);

ftrsNames{1}='StartPoint';
ftrsNames{2}='Num Of Peaks';
ftrsNames{3}='Pf_1_3';
ftrsNames{4}='Pf_2_3';
ftrsNames{5}='Pf_3_3';
ftrsNames{6}='IsTerminalPeak';
ftrsNames{7}='ExpFlowPeak';
ftrsNames{8}='NED';
ftrsNames{9}='X';
ftrsNames{10}='X';
ftrsNames{11}='X';
ftrsNames{12}='InspFlutPow';
ftrsNames{13}='ExpFlutPow';
ftrsNames{14}='ExpFlowLimitationIndex';
ftrsNames{15}='InspFlatnessIndex';
ftrsNames{16}='ExpFlatnessIndex';
ftrsNames{17}='InspSkewness';
ftrsNames{18}='ExpSkewness';
ftrsNames{19}='InspKurt';
ftrsNames{20}='ExpKurt';
ftrsNames{21}='CorrIndexVol';
ftrsNames{22}='InspPeaksRatio';
ftrsNames{23}='X';
ftrsNames{24}='X';
ftrsNames{25}='X';
ftrsNames{26}='X';

ftrs(1,1)=Br.Start;

%figure(3); clf(figure(3)); plot(Flow(Br.Start:Br.End));hold on; % original

% Offset breath times to actual breath
BrMid=Br.Mid-Br.Start;
BrEnd=Br.End-Br.Start;

%plot(BrMid,Flow(Br.Mid),'rx');

ValidFlag=1;
I1=round(0.05*BrMid);
if I1==0
    I1=1;
end
I2=round(0.95*BrMid);

if Br.VT<0.02 || Br.VT>Br.VT95
    ValidFlag=0;
end

% Time=(0:length(Flow)-1)/Fs;

Flutter=FindWaveLetCoefs(Flow,'db4',Fs,'Flow Noise');
Flow1=Flow-Flutter;

PeakFindThr=0.01*max(Flow1);

if PeakFindThr<0
    ValidFlag=0;
end

%% ALI'S FEATURES OPTIMIZED

%% find peak flow
if ValidFlag==1

    [Flow_peaks,Flow_troughs] = peakdetOriginal(Flow1(1:BrMid),PeakFindThr); %%max e min of flow inspiration [index, value]

    Flow_troughs=[1 Flow1(1);Flow_troughs]; %%add start and end of insp as minima
    Flow_troughs=[Flow_troughs; BrMid Flow1(BrMid)];

    if ~isempty(Flow_peaks)
        while Flow_peaks(end,1)>=Flow_troughs(end,1)%%last peak must be before a min
            Flow_peaks(end,:)=[];
            if isempty(Flow_peaks)
                break
            end
        end
    end

    if ~isempty(Flow_peaks)
        while Flow_peaks(1,1)<=Flow_troughs(1,1) %%first peak must be after a min
            Flow_peaks(1,:)=[];
            if isempty(Flow_peaks)
                break
            end
        end

        jj_tr=[];
        for ii=2:size(Flow_troughs,1) %%find each peak before and after the first min (excluded the starting point of the breath)
            fltemp=find(Flow_peaks(:,1)>Flow_troughs(ii-1,1) & Flow_peaks(:,1)<Flow_troughs(ii,1),1);
            if isempty(fltemp)
                jj_tr=[jj_tr;ii-1];
            end
        end
        Flow_troughs(jj_tr,:)=[];
        pkTro1=Flow_peaks(:,2)-Flow_troughs(1:end-1,2);
        pkTro2=Flow_peaks(:,2)-Flow_troughs(2:end,2);

        pkTros=sort([pkTro1;pkTro2],'descend');

        if showfigures
            figure(3),plot(Flow1(1:BrMid));hold on;plot(Flow_peaks(:,1),Flow_peaks(:,2),'*g');plot(Flow_troughs(:,1),Flow_troughs(:,2),'*r');
            hold off;
        end

        if length(pkTros)>=8
            PeakFindThr = max(pkTros(8),0.03);
        else
            PeakFindThr = 0.03;
        end
        itr=0;
        while 1
            [minpkTro1,minpkTro1_i] = min(pkTro1);
            [minpkTro2,minpkTro2_i] = min(pkTro2);
            maxpkTro1 = max(pkTro1);
            maxpkTro2 = max(pkTro2);
            maxpkTro = max([maxpkTro1,maxpkTro2]);

            [minpkTro,Pkpattern] = min([minpkTro1,minpkTro2]);

            itr=itr+1;
            if isempty(minpkTro) || itr>100
                break
            else
                if minpkTro>PeakFindThr
                    break
                end
            end
            if maxpkTro<PeakFindThr
                ValidFlag=0;
                break
            end

            if Pkpattern==2
                Flow_peaks(minpkTro2_i,:)=[];
                Flow_troughs(minpkTro2_i+1,:)=[];
            elseif Pkpattern==1
                Flow_peaks(minpkTro1_i,:)=[];
                Flow_troughs(minpkTro1_i,:)=[];
            end
            %recalculate
            pkTro1=Flow_peaks(:,2)-Flow_troughs(1:end-1,2);
            pkTro2=Flow_peaks(:,2)-Flow_troughs(2:end,2);
        end
    end
end
fp1=NaN;fp1_idx=NaN;
fp2=NaN;fp2_idx=NaN;
fp3=NaN;fp3_idx=NaN;
peakFlow=NaN;peakFlow_idx=NaN;
NumberofPeaks=0;
IsTerminalPeak=0;

if ValidFlag==1
    if ~isempty(Flow_peaks) 
       i_temp=find(Flow_peaks(:,1)<=BrMid/3);
       if ~isempty(i_temp)
           NumberofPeaks=NumberofPeaks+1;
           [fp1,idx_temp]=max(Flow_peaks(i_temp,2));
           fp1_idx=Flow_peaks(i_temp(idx_temp),1)/Fs;
       end
       i_temp=find(Flow_peaks(:,1)>BrMid/3 & Flow_peaks(:,1)<=2*BrMid/3);
       if ~isempty(i_temp)
           NumberofPeaks=NumberofPeaks+1;
           [fp2,idx_temp]=max(Flow_peaks(i_temp,2));
           fp2_idx=Flow_peaks(i_temp(idx_temp),1)/Fs;
       end
       i_temp=find(Flow_peaks(:,1)>2*BrMid/3 & Flow_peaks(:,1)<3*BrMid/3);
       if ~isempty(i_temp)
           NumberofPeaks=NumberofPeaks+1;
           [fp3,idx_temp]=max(Flow_peaks(i_temp,2));
           fp3_idx=Flow_peaks(i_temp(idx_temp),1)/Fs;
       end
       [peakFlow,idx_temp]=max(Flow_peaks(:,2));
       peakFlow_idx=Flow_peaks(idx_temp,1)/Fs;
       if ~isnan(fp3) && fp3==peakFlow
           IsTerminalPeak=1;
       end    
    end
else
    NumberofPeaks=NaN;
    IsTerminalPeak=NaN;
end

if showfigures
    figure(3),plot(Flow1(1:BrMid));hold on;plot(Flow_peaks(:,1),Flow_peaks(:,2),'*g');plot(Flow_troughs(:,1),Flow_troughs(:,2),'*r');
    hold off;
end

ftrs(1,2)=NumberofPeaks;
ftrs(1,3)=fp1;
ftrs(1,4)=fp2;
ftrs(1,5)=fp3;
ftrs(1,6)=IsTerminalPeak;

%% Exp flow peak

if ValidFlag==1
    [ExpFlowPeak,ExpFlowPeakIdx]=min(Flow1);
    PeakExpFlowTimeFromBrMid=(ExpFlowPeakIdx-BrMid)/Fs;
    ftrs(1,7)=ExpFlowPeak;


else
    ftrs(1,7)=NaN;

end

%% find NED and ratio of flow at different times
%a=round([0.1:0.1:0.9]'*BrMid/2);
%b=BrMid-round([0.1:0.1:0.9]'*BrMid/2);
 
  NED=NaN;
    if (~isnan(fp1) && ValidFlag==1)
        
        NED=100*(fp1-Flow1(round(BrMid/2)))/fp1;
        
        if NED<0
            NED=0;
        end
    elseif (~isnan(fp2) && ValidFlag==1)
        
        NED=100*(fp2-Flow1(round(BrMid/2)))/fp2;
        
        if NED<0
            NED=0;
        end
    elseif (~isnan(fp3) && ValidFlag==1)
        
        NED=100*(fp3-Flow1(round(BrMid/2)))/fp3;
    
        if NED<0
            NED=0;
        end
        
    end
    
    if NED>100
       NED=100;
    end
    
    ftrs(1,8)=NED;

%% find flow derivative and greatest change in the flow between 0.05%-0.95% from start point

if ValidFlag==1 

    ftrs(1,9)=0; %% they were not necessary so setted to 0
    ftrs(1,10)=0;
    ftrs(1,11)=0;

else

    ftrs(1,9)=NaN;
    ftrs(1,10)=NaN;
    ftrs(1,11)=NaN;
  
end

%% Uncomment here otherwise to calculate them:

% if ValidFlag==1
%     VboxParams.L=5;
%     VboxParams.tetha=0.8;
%     VboxParams.H=0.1;
%     
%     % Calculating smooth rate of change
%     [~,~,~,FlowDiff]=Vbox(Flow,VboxParams,Fs);
%     flowDiffSmooth=FlowDiff(I1:I2);
%     
%     MinFlowChange=nanmin(flowDiffSmooth);
%     MedianFlowChange=nanmedian(flowDiffSmooth);
%     STDFlowChange=nanstd(flowDiffSmooth);
%     
%     ftrs(1,9)=MinFlowChange;
%     ftrs(1,10)=MedianFlowChange;
%     ftrs(1,11)=STDFlowChange;
% else
%     
%     ftrs(1,9)=NaN;
%     ftrs(1,10)=NaN;
%     ftrs(1,11)=NaN;
% end

%% Find inspiratory and expiratory fluttering power

% if ValidFlag==1
%      Flutter1=Flutter(1:BrMid); 
%      Flutter2=Flutter(BrMid:BrEnd);
%     
%     [Pxx1,f1]=pwelch(Flutter1,3,[],256,Fs);
%     [Pxx2,f2]=pwelch(Flutter2,3,[],[],Fs);
%  
%     Pow1=10*log10(bandpower(Pxx1,f1,[0 f1(end)],'psd')); % approx integral of PSD(Flutter)
%     Pow2=10*log10(bandpower(Pxx2,f2,[0 f2(end)],'psd'));
%     
%     ftrs(1,12)=Pow1;
%     ftrs(1,13)=Pow2;
% else
%     ftrs(1,12)=NaN;
%     ftrs(1,13)=NaN;
% end

if ValidFlag==1
    
    ftrs(1,12)=0;
    ftrs(1,13)=0;
else
    ftrs(1,12)=NaN;
    ftrs(1,13)=NaN;
end


%% Find Expiratory Flow Limitation Index
if ValidFlag==1
    ftrs(1,14)=IsExpFlowLimited(Flow1,1,BrEnd,Fs);
else
    ftrs(1,14)=NaN;
end

%% ALESSANDRA AND ANGELO'S FEATURES
% shift the Flow1 signal to start at the start of the breath
%Flow1 = Flow1(Br.Start:Br.End); % line added by DLM

% shift the struct times to start at zero
Br.End=Br.End-Br.Start;
Br.Mid=Br.Mid-Br.Start;
Br.Start=0;

middleInsp=round(Br.Mid/2);
middleExp=round((Br.End-Br.Mid)/2);
dt=1/Fs;

FlowInsp=Flow1(1:Br.Mid);
FlowExp=Flow1(Br.Mid:end);
Vol=cumsum(Flow1)*dt;
maxVol=max(abs(Vol));

FlowInsp(FlowInsp<0)=0; % force postive, or zero
FlowInsp=FlowInsp/max(FlowInsp);

FlowExp(FlowExp>0)=0; % force negative, or zero
FlowExp=FlowExp/min(FlowExp);
flat=0;

%% Inspiratory and Expiratory Flatness index
if ValidFlag
    flat(1)=sum(FlowInsp>0.9)/sum(FlowInsp>0.2);
    flat(2)=sum(FlowExp>0.85)/sum(FlowExp>0.4);
else
    flat = [NaN,NaN];
end
ftrs(1,15)=flat(1);
ftrs(1,16)=flat(2);

%% Inspiratory and Expiratory Skewness
if ValidFlag==1   
    if~isnan(FlowInsp) & ~isnan(middleInsp)
        AreaInspFlow=trapz(FlowInsp)*dt;
        areaL=trapz(FlowInsp(1:middleInsp))*dt;
        areaR=trapz(FlowInsp(middleInsp:end))*dt;
        ftrs(1,17)=(areaL/AreaInspFlow)-(areaR/AreaInspFlow);
    end
    if~isnan(FlowExp) & ~isnan(middleExp)
        AreaExpFlow=trapz(FlowExp)*dt;
        areaL=trapz(FlowExp(1:middleExp))*dt;
        areaR=trapz(FlowExp(middleExp:end))*dt;
        ftrs(1,18)=(areaL/AreaExpFlow)-(areaR/AreaExpFlow);
    end
else
    ftrs(1,17)=NaN;
    ftrs(1,18)=NaN;
end

%% Inspiratory and Expiratory Kurtosis
if ValidFlag==1
    if~isnan(FlowInsp)
        ftrs(1,19)=kurtosis(FlowInsp,0);
    end 
else  
    ftrs(1,19)=NaN;  
end

if ValidFlag==1   
    if~isnan(FlowExp)
        ftrs(1,20)=kurtosis(FlowExp,0);
    end    
else
    ftrs(1,20)=NaN;
end

%% Volume waveform correlation with template
if 0
    
    if~isnan(TemplateVol)
        
        if~isnan(Vol) && ~isnan(maxVol)
            Vol=Vol/norm(Vol,Inf);
            VolShape=interp1(linspace(1,1000, length(Vol)),(Vol),1:1000,'spline');
            figure();plot(VolShape);
            ftrs(1,21)=corr(VolShape',TemplateVol');
        end
        
    end
    
else
    
    ftrs(1,21)=NaN;
    
end

%% Inspiratory peaks ratio
if ValidFlag==1  
        P1P=fp1;
        P2P=fp2;
        P3P=fp3;      
        P=[P1P P2P P3P];
        P(isnan(P))=0; % replace NaN's with zeros
        %PmaxP=max(P);
        %P=P./PmaxP;            
        %ftrs(1,22)=sum(P)/3;
        ftrs(1,22)=sum(P./max(P))/3; % DLM addedd this line to replace above
else   
        ftrs(1,22)=NaN;
end

%% Ali features about discontinuity points
if ValidFlag==1  
    ftrs(1,23)=0; %% they were not necessary so setted to 0
    ftrs(1,24)=0;
    ftrs(1,25)=0;
    ftrs(1,26)=0;
else
    ftrs(1,23)=NaN;
    ftrs(1,24)=NaN;
    ftrs(1,25)=NaN;
    ftrs(1,26)=NaN;  
end

end