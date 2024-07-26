function [PrUpright,Fnoise2All,Info]=FlowInvertedDetectTool(FlowSig,TimeSig,Info)



global settings
Time = [TimeSig(1):(1./settings.Fs):TimeSig(end)]';
Flow = interp1(TimeSig,FlowSig,Time);
clear TimeSig FlowSig

% clear all;
% close all;
% clc;

% addpath(genpath('C:\Users\rma56\Dropbox (Partners HealthCare)\PUPbeta_git\PUPbeta\'))
% load('D:\MrOS\Visit1\Converted\mros-visit1-aa0003_XHz.mat')

% load model
load('FlowInvertedDetector.mat')

% settings required
%need to have these defined before here, let's not overwrite the global variable in here
if 0
%     settings.plotbreathdetectionfigures=0;
%     settings.sqrt_scaling=1;
%     settings.plotfiguresqrtscaling=0;
%     settings.modBB_i_start = logical(1); % modifybreathstartusingflowsilence
%     settings.scalingexponent=0.67;
%     settings.plotfigure=0;
end

PrctileList = [12.5:12.5:87.5];

% FlowSig = -DataEventHypnog_Mat(:,2);
% TimeSig = DataEventHypnog_Mat(:,1);

%overall noise
noisewavAll = FlowSignalToNoise(Time,Flow,settings.plotfigure); %replace with Flow and Time
Fnoise2All=sum(noisewavAll>=2)/length(noisewavAll);

figure(898);clf(898);
set(gcf,'color',[1 1 1]);

clear Tone
Tone.diffcurr = nan(length(PrctileList),1);
Tone.diffprev = nan(length(PrctileList),1);
Tone.COVVTi = nan(length(PrctileList),1);

for i=1:length(PrctileList)
  figure(898);  
    try
        lt = prctile(Time,PrctileList(i));
        I = find(Time>lt & Time<lt+420);
        
        
        %plot(TimeSig(I),FlowSig(I));
        tempx = Time(I);
        temp = Flow(I);
        temp = (temp-nanmean(temp))/nanstd(temp);
        subplot(length(PrctileList),1,i);
        plot(tempx,temp);
        set(gca,'box','off','ylim',[-4 4],'xlim',[min(tempx) max(tempx)],'xtick',[],'xcolor',[1 1 1],'ytick',[],'ycolor',[1 1 1]);
        if i==1
            h=ylabel('Example signals','fontweight','normal','color',[0 0 0]);
            if isfield(settings,'fname')
                title(settings.fname)
            end
        end
        
        %% adding code for SNR breath and window
        
        if 1
            TimeWin=Time(I);
            FlowWin=Flow(I);
            [~,~,BB_i_start,BB_i_mid,BB_i_end,~,~,~,~,leak,...
                IEratio,VT,~,~,Apnea_B,~,VTi,VTe,~,~,leak_B,IEratioEstimated] =...
                VEfromFlow(TimeWin,FlowWin);
            
            Time0 = Time(1)+0*BB_i_start;
            IEratio = IEratio+0*BB_i_start;
            IEratioEstimated = IEratioEstimated+0*BB_i_start;
            leak_A = leak+0*BB_i_start;
            
            BreathDataList = {'Time0','Time_start','Time_mid','Time_end','BB_i_start','BB_i_mid','BB_i_end',...
                'IEratio','IEratioEstimated','leak','leak2'};
            
            BreathDataTable=table(Time0, TimeWin(BB_i_start), TimeWin(BB_i_mid), TimeWin(BB_i_end), BB_i_start, BB_i_mid, BB_i_end, ...
                IEratio(:),IEratioEstimated(:),leak_A(:),leak_B(:),'VariableNames', BreathDataList);
            
%             starttime =BreathDataTable.Time0(1);
%             endtime=starttime+(settings.windowlength*60);
%             starti= find(Time==starttime); % using Time here and not TimeFlow--error in SignalToNoiseWin.m line 123
%             endi= find(Time>=endtime,1,'first')-1; % some time issues in windows

            BreathDataTable = SignalToNoiseWin(BreathDataTable,FlowWin,TimeWin);
            
            
            IEratioW(i,1)=nanmean(IEratio);
            IEratioW2(i,1)=nanmean(10.^abs(log10(IEratio)));
            IEratioEstimatedW(i,1)=nanmean(IEratioEstimated);
            IEratioEstimatedW2(i,1)=nanmean(10.^abs(log10(IEratioEstimated)));
%             SNRBreath(i,1)=nanmean(BreathDataTable.SNRbreath);
            SNRwindow(i,1)=nanmedian(BreathDataTable.SNRwindow);
%             leakW(i,1)=nanmean(BreathDataTable.leak);
%             leakBW(i,1)=nanmean(BreathDataTable.leak2);
            
        end
        
        %add code to detct if flow is there (e.g. Noise code, uses fft)
        noisewav = FlowSignalToNoise(tempx,temp,settings.plotfigure);
        
        Fnoise2(i,1)=sum(noisewav>=2)/length(noisewav);
        if Fnoise2(i,1)>=0.1
            %display(['Skipping ' num2str(i) ':' ', noisy or absent flow signal']);
            Tone.diffcurr(i,1)=NaN;
            Tone.diffprev(i,1)=NaN;
            Tone.COVVTi(i,1)=NaN;
        else
            %             [~,~,~,~,~,~,~,~,~,~,...
            %                 ~,VT,~,~,Apnea_B,~,VTi,VTe,Ti,Te] =...
            %                 VEfromFlow(TimeSig(I),FlowSig(I));
            
            [Tone.diffcurr(i,1),Tone.diffprev(i,1),Tone.COVVTi(i,1)] = FlowInvertedParameters(VTi,VTe,VT,Apnea_B);
        end
    catch me
    end
end


Tone = struct2table(Tone);
Tone.diffcurroverprev = Tone.diffcurr./Tone.diffprev;

Tone.PappearsInverted = 1-predict(FlowInvertedDetector.mdlUpright,Tone);
Tone.Perror = predict(FlowInvertedDetector.mdlError,Tone);

Ttwo=table([]);
MeanPw = nanmean(Tone.PappearsInverted.*(1-Tone.Perror))/nanmean(1-Tone.Perror);
Ttwo=table(MeanPw);
PrUpright = predict(FlowInvertedDetector.mdlSubjUpright,Ttwo);

%% Flow noise metrics--one value for subject
Info.FlowQ.IEratio = nanmedian(IEratioW);
Info.FlowQ.IEratio2 = nanmedian(IEratioW2);
Info.FlowQ.IEratioEstimated = nanmedian(IEratioEstimatedW);
Info.FlowQ.IEratioEstimated2 = nanmedian(IEratioEstimatedW2);
Info.FlowQ.SNRwindow = nanmedian(SNRwindow);
Info.FlowQ.SNRwindowNoiseMean = 10*log10(1./nanmean(1./(10.^([SNRwindow]/10))));
Info.FlowQ.FSNRover15 = nansum(SNRwindow>15)/nansum(SNRwindow>-Inf); %fixed 8/4/2021
Info.FlowQ.FSNRover20 = nansum(SNRwindow>20)/nansum(SNRwindow>-Inf);
Info.FlowQ.FSNRover25 = nansum(SNRwindow>25)/nansum(SNRwindow>-Inf);
Info.FlowQ.FSNRover30 = nansum(SNRwindow>30)/nansum(SNRwindow>-Inf);


%%
if isfield(settings,'savefigure')&& settings.savefigure==1
    savefigdir=[settings.workdir 'FlowQCPlots\'];
    if ~(exist(savefigdir, 'dir') == 7)
        mkdir(savefigdir);
    end
    savefigname=[savefigdir extractBefore(settings.fname,'.edf') '.png'];
    saveas(figure(898),savefigname);
    
    savefigdir2=[settings.workdir 'FlowQCPlots\MatFigs\'];
    if ~(exist(savefigdir2, 'dir') == 7)
        mkdir(savefigdir2);
    end
    savefigname2=[savefigdir2 extractBefore(settings.fname,'.edf') '.fig'];
    
    saveas(figure(898),savefigname2);
end






