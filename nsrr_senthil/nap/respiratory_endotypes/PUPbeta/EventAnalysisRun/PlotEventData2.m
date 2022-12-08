function EAinfo = PlotEventData2(Boxes,Ensembles,EAinfo,criteria)
%%
global settings

ArSignal='EventsAr';
if isfield(settings,'useWSanalysisToReplaceAr') && settings.useWSanalysisToReplaceAr == 1 %0=use original scoring, 1=use best EEG, 2=use "predicted best" EEG.
    ArSignal='EventsArWS';
elseif isfield(settings,'useWSanalysisToReplaceAr') && settings.useWSanalysisToReplaceAr == 2
    ArSignal='EventsArWSB';
end

if ~exist('criteria')
    criteria=true(size(Boxes.Time,1),1)
end

%PtName=EAinfo.FileName;
EnTime=Ensembles.Time(:)';
Xlims = [-45 30];
EvMatrix= Boxes.EventsResp(criteria,:)';

if isfield(settings,'PlotEventData2Option')
    switch settings.PlotEventData2Option
        case 1
            plotchanneltext = { ...
                {'VI','VdriveEdiNorm'}, ...
                {'GGpeak','GGtonic'}, ...
                {ArSignal} ...
                }; %,'WakeSleepVe' ,{'DeltaPes'}     {'Dsat'} ...
            C = {  {[0.2 0.2 1],[0.2 0.2 0.3],[0.95 0.5 0.1],[0.4 0.95 0.3]},  {[0.95 0.5 0.1],[0.95 0.5 0.1]},  ...
                {[0.4 0.95 0.3],[0.95 0.5 0.1],[0.2 0.2 1]} ,  {[0.4 0.95 0.3],[0.95 0.4 0.3],[0.3 0.4 0.95]} ,  {[0.4 0.95 0.3],[0.95 0.4 0.3],[0.3 0.4 0.95]}  };
            
        case 2
            plotchanneltext = { ...
                {'VI'}, ...
                {'FlowDrive'}, ...
                {'SpO2'}, ...
                {'HR'} ...
                {ArSignal} ...
                };
            C = {  {[0.2 0.2 1],[0.2 0.2 0.3],[0.95 0.5 0.1],[0.4 0.95 0.3]},  {[0.95 0.5 0.1],[0.95 0.5 0.1]},  ...
                {[0.4 0.95 0.3],[0.95 0.5 0.1],[0.2 0.2 1]} ,  {[0.4 0.95 0.3],[0.95 0.4 0.3],[0.3 0.4 0.95]} ,  {[0.4 0.95 0.3],[0.95 0.4 0.3],[0.3 0.4 0.95]}  };
            
        case 3 %DEFAULT
            plotchanneltext = { ...
                {'VI','VdriveEdiNorm'}, ...
                {'SpO2'}, ...
                {'HRfilt','Pulse'}, ...
                {'FlowDrive'}, ...
                {ArSignal,'ArPr'},...
                {'WPrLogit','ArPrLogit'} ...
                {'SpO2'} ...
                };
            C = {  {[0.2 0.2 1],[0.2 0.2 0.3],[0.95 0.5 0.1],[0.4 0.95 0.3]},  {[0.95 0.5 0.1],[0.95 0.5 0.1]},  ...
                {[0.4 0.95 0.3],[0.95 0.5 0.1],[0.2 0.2 1]} ,  {[0.4 0.95 0.3],[0.95 0.4 0.3],[0.3 0.4 0.95]} ,  {[0.4 0.95 0.3],[0.95 0.4 0.3],[0.3 0.4 0.95]}  };
            
        case 4
            plotchanneltext = { ...
                {'VI','VdriveEdiNorm'}, ... %,'VdriveEdiNorm'
                {'SpO2'}, ...
                {'ArPr'},...
                {'GGpeak','GGtonic'} %'WPrLogit','ArPrLogit'} ...
                }; %{ArSignal,'WPr','ArPr'} ...,'WakeSleepVe' ,    {'Dsat'} ...
            C = {  {[0.2 0.2 1],[0.2 0.2 0.3],[0.95 0.5 0.1],[0.4 0.95 0.3]},  {[0.95 0.5 0.1],[0.95 0.5 0.1]},  ...
                {[0.4 0.95 0.3],[0.95 0.5 0.1],[0.2 0.2 1]} ,  {[0.4 0.95 0.3],[0.95 0.4 0.3],[0.3 0.4 0.95]} ,  {[0.4 0.95 0.3],[0.95 0.4 0.3],[0.3 0.4 0.95]}  };
            
        case 5
            plotchanneltext = { ...
                {'VI','VdriveEdiNorm','VdrivePesNorm'} ,... %,'VdriveEdiNorm'
                {'EventsAr'} ,...
                {'GGpeakCalibrated','GGtonicCalibrated'} ... %'WPrLogit','ArPrLogit'} ...
                }; %,'WakeSleepVe' ,{'DeltaPes'}     {'Dsat'} ...
            C = {  {[0.2 0.2 1],[0.2 0.2 0.3],[0.95 0.5 0.1],[0.4 0.95 0.3]}, {[0.95 0.5 0.1],[0.95 0.5 0.1]},  ...
                {[0.4 0.95 0.3],[0.95 0.5 0.1],[0.2 0.2 1]} ,  {[0.95 0.4 0.3],[0.95 0.4 0.3],[0.3 0.4 0.95]} ,  {[0.4 0.95 0.3],[0.95 0.4 0.3],[0.3 0.4 0.95]}  };
        case 6
            plotchanneltext = { ...
                {'VI','VdriveEdiNorm'}, ... %,'VdriveEdiNorm'
                {'GGpeak','GGtonic'} %'WPrLogit','ArPrLogit'} ...
                }; %,'WakeSleepVe' ,{'DeltaPes'}     {'Dsat'} ...   {'ArPr'} ...
            C = {  {[0.2 0.2 1],[0.2 0.2 0.3],[0.95 0.5 0.1],[0.4 0.95 0.3]},   ...
                {[0.95 0.4 0.3],[0.95 0.4 0.3],[0.3 0.4 0.95]} , {[0.4 0.95 0.3],[0.95 0.5 0.1],[0.2 0.2 1]}, {[0.4 0.95 0.3],[0.95 0.4 0.3],[0.3 0.4 0.95]}  };
            
        case 7
            plotchanneltext = { ...
                {'VI','VdrivePesNorm','VdrivePmusNorm'}, ... %,'VdriveEdiNorm'
                {'DeltaPes','DeltaPmus'}, ...
                {'FlowDrive'}, ...
                {'SpO2'}, ...
                {'HRfilt','Pulse'}, ...
                {ArSignal,'WPr','ArPr'} ...
                }; %,'WakeSleepVe' ,{'DeltaPes'}     {'Dsat'} ...
            C = { {[0.2 0.2 1],[0.2 0.2 0.3],[0.95 0.5 0.1]},  {[0.2 0.2 0.3],[0.95 0.5 0.1]}, {[0.4 0.95 0.3]}, ...
                {[0.95 0.5 0.1]} ,  {[0.4 0.95 0.3]} ,  {[0.4 0.95 0.3],[0.95 0.4 0.3],[0.3 0.4 0.95]} ,  {[0.4 0.95 0.3]}  };
            
            
        case 8
            plotchanneltext = { ...
                {'DeltaThorax','DeltaAb'},...
                {'VdrivePesNorm'}, ...
                {'VI','VdriveEdiNorm'}, ...
               % {'VdrivePepiNorm','VdrivePesNorm'}...
             %   {'VdriveEdiNorm'}...
                %{'HR'} ...
               % {ArSignal} ...
                };  
                C ={ ...
                %{[0.2 0.2 1]},{[0.4 0.95 0.3],[0.2 0.2 1],[0.95 0.5 0.1],[0.2 0.2 0.3]},...
                 {[0.4 0.95 0.3]},... 
                {[0.2 0.7 0.7],[0.95 0.5 0.1],[0.95 0.5 0.1],[0.4 0.95 0.3]},...
                {[0.2 0.2 0.3],[0.95 0.5 0.1],[0.4 0.95 0.3]},...         
                {[0.2 0.7 0.7],[0.95 0.5 0.1],[0.2 0.2 0.3]},...
           {[0.2 0.2 0.3],[0.95 0.5 0.1],[0.4 0.95 0.3]}};
            %    
        case 9
            plotchanneltext = { ...
                {'VI'},... %,'PdriveK','VdrivePdriveKrawNorm'}, ... %,'VdriveEdiNorm'
                {'GGpeak'},...  %  %''GGpeak'',''GGpeak''} ...
                {'VdriveEdiNorm','VdrivePesNorm'}...
                }; %{ArSignal,'WPr','ArPr'} ...,'WakeSleepVe' ,    {'Dsat'} ...
            C = { ...
                {[0.2 0.2 1],[0.4 0.95 0.3],[0.2 0.2 1]},...
                {[0.4 0.95 0.3],[0.2 0.2 1],[0.95 0.5 0.1],[0.2 0.2 0.3]},...
                {[0.2 0.2 0.3],[0.95 0.5 0.1],[0.4 0.95 0.3]},...
                {[0.2 0.2 0.3],[0.95 0.5 0.1],[0.4 0.95 0.3]}};
            %                 C = {  {[0.2 0.2 1],[0.2 0.2 0.3],[0.95 0.5 0.1],[0.4 0.5 0.3]},  {[0.95 0.5 0.1],[0.95 0.5 0.1]},  ...
            %     {[0.4 0.95 0.3],[0.95 0.5 0.1],[0.2 0.2 1]} ,  {[0.4 0.95 0.3],[0.95 0.4 0.3],[0.3 0.4 0.95]} ,  {[0.4 0.95 0.3],[0.95 0.4 0.3],[0.3 0.4 0.95]}  };
            
        case 11
            plotchanneltext = { ...
                {'VI'} ,... %,'PdriveK','VdrivePdriveKrawNorm'}, ... %,'VdriveEdiNorm'
                {'ArPr',ArSignal},...
                {'FlowDrive'},...  %  %''GGpeak'',''GGpeak''} ...
                {'SpO2'},...
                {'HRfilt','Pulse'},...
                {'AA_ExpFlat90Te_O'},...
                {'Tpeak1_Ti_T'},...
                {'VI'}... 
                }; %{ArSignal,'WPr','ArPr'} ...,'WakeSleepVe' ,    {'Dsat'} ...
            C = { ...
                 {[0.4 0.95 0.3],[0.2 0.2 1],[0.95 0.5 0.1],[0.2 0.2 0.3]}...
                ,{[0.4 0.95 0.3],[0.2 0.2 1],[0.95 0.5 0.1],[0.2 0.2 0.3]}...
                ,{[0.2 0.2 0.3],[0.95 0.5 0.1],[0.4 0.95 0.3]}...
                ,{[0.4 0.95 0.3],[0.2 0.2 1],[0.95 0.5 0.1],[0.2 0.2 0.3]}...
                ,{[0.4 0.95 0.3],[0.2 0.2 1],[0.95 0.5 0.1],[0.2 0.2 0.3]}...
                ,{[0.4 0.95 0.3],[0.2 0.2 1],[0.95 0.5 0.1],[0.2 0.2 0.3]}...
                ,{[0.4 0.95 0.3],[0.2 0.2 1],[0.95 0.5 0.1],[0.2 0.2 0.3]}...
                ,{[0.4 0.95 0.3],[0.2 0.2 1],[0.95 0.5 0.1],[0.2 0.2 0.3]}...
            };
            %                 C = {  {[0.2 0.2 1],[0.2 0.2 0.3],[0.95 0.5 0.1],[0.4 0.5 0.3]},  {[0.95 0.5 0.1],[0.95 0.5 0.1]},  ...
            %     {[0.4 0.95 0.3],[0.95 0.5 0.1],[0.2 0.2 1]} ,  {[0.4 0.95 0.3],[0.95 0.4 0.3],[0.3 0.4 0.95]} ,  {[0.4 0.95 0.3],[0.95 0.4 0.3],[0.3 0.4 0.95]}  };
            
        case 12
            plotchanneltext = { ...
                {'VI'} ,... %,'PdriveK','VdrivePdriveKrawNorm'}, ... %,'VdriveEdiNorm'
                {'ArPr',ArSignal},...
                {'FlowDrive'},...  %  %''GGpeak'',''GGpeak''} ...
                {'SpO2'},...
                {'HRfilt','Pulse'},...
                {'DeltaPepi'},...
                {'YFOTmeanI','YFOTmeanE'},...
                {'VI'}... 
                }; %{ArSignal,'WPr','ArPr'} ...,'WakeSleepVe' ,    {'Dsat'} ...
                
            C = {
                {[0.4 0.95 0.3],[0.2 0.2 1],[0.95 0.5 0.1],[0.2 0.2 0.3]}...
                ,{[0.4 0.95 0.3],[0.2 0.2 1],[0.95 0.5 0.1],[0.2 0.2 0.3]}...
                ,{[0.2 0.2 0.3],[0.95 0.5 0.1],[0.4 0.95 0.3]}...
                ,{[0.4 0.95 0.3],[0.2 0.2 1],[0.95 0.5 0.1],[0.2 0.2 0.3]}...
                ,{[0.4 0.95 0.3],[0.2 0.2 1],[0.95 0.5 0.1],[0.2 0.2 0.3]}...
                ,{[0.4 0.95 0.3],[0.2 0.2 1],[0.95 0.5 0.1],[0.2 0.2 0.3]}...
                ,{[0.4 0.95 0.3],[0.2 0.2 1],[0.95 0.5 0.1],[0.2 0.2 0.3]}...
                ,{[0.4 0.95 0.3],[0.2 0.2 1],[0.95 0.5 0.1],[0.2 0.2 0.3]}... %G,B,R,black
                ,{[0.4 0.95 0.3],[0.2 0.2 1],[0.95 0.5 0.1],[0.2 0.2 0.3]}...
            };
            
            %                 C = {  {[0.2 0.2 1],[0.2 0.2 0.3],[0.95 0.5 0.1],[0.4 0.5 0.3]},  {[0.95 0.5 0.1],[0.95 0.5 0.1]},  ...
            %     {[0.4 0.95 0.3],[0.95 0.5 0.1],[0.2 0.2 1]} ,  {[0.4 0.95 0.3],[0.95 0.4 0.3],[0.3 0.4 0.95]} ,  {[0.4 0.95 0.3],[0.95 0.4 0.3],[0.3 0.4 0.95]}  };
                
        otherwise
            plotchanneltext = { ...
                {'VI'}, ...
                {'SpO2'}, ...
                {'HR'} ...
                {ArSignal} ...
                };
            
            C = {  {[0.2 0.2 1],[0.2 0.2 0.3],[0.95 0.5 0.1],[0.4 0.95 0.3]},  {[0.95 0.5 0.1],[0.95 0.5 0.1]},  ...
                {[0.4 0.95 0.3],[0.95 0.5 0.1],[0.2 0.2 1]} ,  {[0.4 0.95 0.3],[0.95 0.4 0.3],[0.3 0.4 0.95]} ,  {[0.4 0.95 0.3],[0.95 0.4 0.3],[0.3 0.4 0.95]}  };
            
    end
end

%define more precise event times using ensemble VE
VEmean = nanmean(Boxes.VI(criteria,:));

thres=1;


temp = VEmean>thres;
I=find(diff(temp)==1);
[~,i]=min(abs(EnTime(I))); %find the surge in VI to >thres closest to time=0
I2=I(i);

deltaT = interp1(VEmean([I2 I2+1]),EnTime([I2 I2+1]),thres);


%find start event time
I3=find(diff(temp)==-1); % added 1 to get the next second, within the event
I3(I3>I2)=[]; %remove post event "start events"
if ~isempty(I3)
    I3=I3(end); %get last index
    deltaT2 = interp1(VEmean([I3 I3+1]),EnTime([I3 I3+1]),thres);
else
    deltaT2=-100;
    disp('warning no 100 detected');
end

if 1
%find start event time
I3=find(diff(VEmean>0.9)==-1); % added 1 to get the next second, within the event
I3(I3>I2)=[]; %remove post event "start events"
if ~isempty(I3)
    I3=I3(end); %get last index
    deltaT90 = interp1(VEmean([I3 I3+1]),EnTime([I3 I3+1]),0.9);
    deltaT2estT90 = deltaT90*1.1;
    if deltaT2<deltaT2estT90*1.1 
        deltaT2=deltaT2estT90*1.1;
        disp('using 90');
    end
else
    deltaT2=-100;
    disp('warning no 90 detected');
end

I3=find(diff(VEmean>0.8)==-1); % added 1 to get the next second, within the event
I3(I3>I2)=[]; %remove post event "start events"
if ~isempty(I3)
    I3=I3(end); %get last index
    deltaT80 = interp1(VEmean([I3 I3+1]),EnTime([I3 I3+1]),0.8);
    deltaT2estT80 = deltaT80*1.22;
    if deltaT2<deltaT2estT80*1.1 
        deltaT2=deltaT2estT80*1.1;
        disp('using 80');
    end
else
    deltaT2=-20;
    disp('warning no 80 detected');
end




%%
%deltaT2=-EAinfo.meanDuration
if 1
    EnTime = EnTime - deltaT;
    deltaT2=deltaT2-deltaT;
    deltaT=0;
    EAinfo.EnTime = EnTime;
    EAinfo.deltaT2 = deltaT2;
end


%% Figure
figure(15); clf(15);
set(gcf,'Position',[95 95 375 860])
Nsubplots = length(plotchanneltext) + 1;
set(gcf,'color',[1 1 1]);
for i=1:length(plotchanneltext)
% currently causing problems:
%     if i==1 & settings.PlotEventData2Option~=[8 9 2]
%         ax(i)=subplot(Nsubplots,1,i:(i+1)); %takes up two rows at this level, makes signal larger
%     elseif i==4 & settings.PlotEventData2Option==[8]
%         
%         ax(i)=subplot(Nsubplots,1,i:(i+1));
%         
%     else
%         ax(i)=subplot(Nsubplots,1,i);
%     end
    ax(i)=subplot(Nsubplots,1,i);
    box('off');
    J=length(plotchanneltext{i});
    
    for j=1:J
        try
            ydatastr = plotchanneltext{i}{j};
            ydata=getfield(Boxes,ydatastr);
            
            ydata = ydata(criteria,:);
            ndata = sum(~isnan(ydata));
            % Special edits
            if settings.PlotEventData2Option==1
                switch i
                    case 1
                        ydata=ydata*100;
                    case 4
                        ydata=ydata*100;
                    case 5
                        ydata=ydata*100;
                end
            end   
            if settings.PlotEventData2Option==9
                switch i
                    
                    case 1
                        ydata=ydata*100;
                    case 3
                        ydata=ydata*100;
                end
            end
            
            ydatamean = nanmean(ydata);
            
            %Default color
            c = [0.5 0.5 0.5];
            try
                c = C{i}{j};
            end
            plot(EnTime,ydatamean,'-','color',c.^10,'linewidth',1.5);
            hold('on');
            
            ydataupper95CI = ydatamean + 1.96*nanstd(ydata)./sqrt(ndata);
            ydatalower95CI = ydatamean - 1.96*nanstd(ydata)./sqrt(ndata);
            h=fill([EnTime fliplr(EnTime)],[ydataupper95CI fliplr(ydatalower95CI)],c);
            set(h,'FaceAlpha',0.5,'EdgeAlpha',0);
            set(gca,'tickdir','out');
            
        end
    end
    
    
    box('off');
     set(gca,'fontsize',10)
    % hold('off');
    if i<length(plotchanneltext)
        set(gca,'xtick',[],'xcolor',[1 1 1])
    end
    
    
    
    
    
    %     posA = get(gca,'position');
    %     posA2 = posA;
    %     posA2(1)=posA(1) - 0.025;
    %     posA2(2)=posA(2) - 0.025;
    %     posA2(3)=1-(1-posA(3)-posA(1))*0.10-posA(1);
    %     posA2(4)=posA(4) + 0.035;
    %     set(gca,'position',posA2);
    
    % Special edits
    switch settings.PlotEventData2Option
        case {1,5}
            switch i
                case 1
                    temp = get(gca,'ylim');
                    temp(1)=0;
                    set(gca,'ylim',temp);
                    
                    ylim([0 400]); 
                case 3
                    %  ylim([0 1]);
                    
                    %                 case 2
                    %                     ylim([92 110]);
            end
        case 6
            switch i
                case 2
                    ylim([0 350]);
            end
      
    end
    
    ylimtemp=get(gca,'ylim');
    plot(0*[1 1],ylimtemp,'-','color',[1 0.5 0.5]);
    hold on
    plot(deltaT2*[1 1],ylimtemp,'--','color',[1 0.5 0.5]);
    chH = get(gca,'Children');
    set(gca,'Children',flipud(chH));
    
    box('off');
    
    xlim(Xlims);
    linkaxes(ax,'x');
    
    
    
end

end
