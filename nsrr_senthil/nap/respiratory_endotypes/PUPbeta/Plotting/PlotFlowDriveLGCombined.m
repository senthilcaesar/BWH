% 2021-11-26 17:21:22 Friday EAS added this have flow:drive and LG on same
% plot
function [] = PlotFlowDriveLGCombined(SigT,BreathDataTable,Vflow_out, VI, Veupnea, E1, AR3, TimeB, Ecentralapnea, Ttot_ser, Ttot_ser_prev, AR1, E, BB_i_start, BB_i_mid, BB_i_end, BB_t, VraOn, Ab, Th, Time, hypnog, spo2, Arousal, position_mode, Percent_position)

global n winNum settings ChannelsList

%% MAIN PLOT
while 1
    Ventilation_EventScoring_data=[VI(:)-Veupnea E1(:) AR3(1:length(E1)) TimeB(:) Ecentralapnea(:) [Ttot_ser(:); 0] Ttot_ser_prev(:)];
    
    %find arousal onset
    ARonset=AR1;
    ARonset(1:end)=0;
    for i=2:1:length(AR1)
        if AR1(i)==1&&AR1(i-1)==0
            ARonset(i)=1;
        end
    end
    
    %counting events, note this only counts obstructive events.
    N_events=0;
    for k=2:length(E)
        if E(k)==0&&E(k-1)==1
            N_events=N_events+1;
        end
    end
    
    
    modelpolyfitorder=3;
    temp_i_i=10;
    if sum(AR3((temp_i_i+1):end)==0)==0
        VraOn=0;
    end
    [Error,Vdr_est,VAr_est,LoopGainParameters, Optimal_vent_control_parameters,Optimal_MSE,FitQual,i_1,i_end,...
        lowerSEM,upperSEM,CI_parameters] = ...
        FindBestModelParametersCI(Ventilation_EventScoring_data,VraOn,Veupnea,modelpolyfitorder);
    %[Error,Vdr_est,VAr_est,LoopGainParameters,Optimal_vent_control_parameters,Optimal_MSE,FitQual,i_1,i_end] = ...
    %   FindBestModelParameters(Ventilation_EventScoring_data,VraOn,Veupnea,modelpolyfitorder);
    
    if VraOn==0
        LoopGainParameters(8)=NaN; %VRA
    end
    
    %Vdr_est=Vdr_est_X
    
    if isnan(LoopGainParameters(1))==1
        WinInfo(1:22)=NaN;
        disp('exiting 6');
        return
    end
    
    Arthressample=Vdr_est+Veupnea;
    Arthressample([1:i_1 i_end:end])=[];
    ARonset([1:i_1 i_end:end])=[];
    Arthressample(ARonset==0)=[];
    Arthres=mean(Arthressample);
    N_arousals_used=sum(ARonset);
    
    % calculate this outside of plot if statement
    meanSEM = nanmean(upperSEM(i_1:end)-lowerSEM(i_1:end));
    
    %% 
    RIPinfo='';
    
    %% Plot
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    %If plotfigure=1, create plot
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    if settings.plotfigure==1
        LGfig = 333;
        figure(LGfig); clf(LGfig); LGfigID = gcf;
        LGfigID.Color = [1 1 1];
        LGfigID.Units = 'Inches';
        %LGfigID.Position = [8 3 11 8]; %removed, irritating :)
        
        if 1
            plotflow=Vflow_out; plotflow_t=Time;
        else
            plotflow=Vflow; plotflow_t=time;
        end
        if 1
            RIP_Thorax_plot=Th;
            RIP_Abdo_plot=Ab;
            plotfilteredRIP=0;
            if plotfilteredRIP
                RIP_Thorax_plot = filtfilt(B_butter1,A_butter1,RIP_Thorax_plot);
                RIP_Abdo_plot = filtfilt(B_butter1,A_butter1,RIP_Abdo_plot);
            end
        end



        set(gcf,'Color',[1 1 1]);
        ax1(1)=subplot('Position',[0.131 0.782 0.775 0.071]); % hypnogram plot
        
        plot(Time,hypnog,'k');
        
        %plot(BB_t(AR1>0),AR1(AR1>0),'r.',BB_t(E1<1),E1(E1<1),'k.',BB_t(E<1),E(E<1)+1,'c.',BB_t(Ecentralhypop>0),Ecentralhypop(Ecentralhypop>0)-2,'b.',...
        %    BB_t(Ecentralapnea>0),Ecentralapnea(Ecentralapnea>0)-3,'g.');
        set(gca,'XColor',[1 1 1],'YColor',[0 0 0],'ylim',[-0.1 5.1],...
            'FontSize',7,...
            'YTickLabel',{'3','2','1','R','W'},...
            'YTick',[0:4],...
            'TickDir','out',...
            'TickLength',[0 0],...
            'box','off',...
            'FontName','Arial');
        %line([Time(1) Time(1)],[-0.1 5.1]);
        xlim([min(Time) max(Time)]);
        title(['Patient ' num2str(n) ', Window ' num2str(winNum) ' ' RIPinfo],'fontname','Arial Narrow','fontsize',12);

        %% signals plot
        ax1(2)=subplot('Position',[0.131 0.455 0.775 0.321]); 
        
        box('off');
        xlim([min(Time) max(Time)]);
        hold('on')
        
        alpha_=98.5; % Mesa version set this to 99.5
        nflow=1.5*(prctile(plotflow,alpha_)-prctile(plotflow,100-alpha_));    mflow=mean(plotflow);
        nTH=2*(prctile(RIP_Thorax_plot,alpha_)-prctile(RIP_Thorax_plot,100-alpha_));    mTH=mean(RIP_Thorax_plot);
        nAB=2*(prctile(RIP_Abdo_plot,alpha_)-prctile(RIP_Abdo_plot,100-alpha_));     mAB=mean(RIP_Abdo_plot);
        minspo2=min(spo2);
        if minspo2>90
            minspo2=90;
        end
        nspo2=1.5*(100-minspo2);
        mspo2=mean([100,minspo2]);
        
        %plot EEG and arousals scored over the top in green
        if max(Arousal)>0
            tempds=4; tempsize=0.2; tempoffset=1.2-0.0;
            tempx=downsample(Time,tempds);
            tempy=tempsize*downsample(Arousal,tempds)+tempoffset;
            tempxx=[tempx;flipud(tempx)];
            tempyy=[tempy;tempoffset+zeros(length(tempy),1)];
            fill(tempxx,tempyy,[0 1 0],'EdgeColor','none','FaceAlpha',0.75);
            %fill(tempxx,tempyy,[0 1 0],'EdgeColor','none');
        end
        if exist('ArousalOrig') && max(ArousalOrig)>0 && 1 %compare scoring, plot original scoring higherup
            tempds=4; tempsize=0.2; tempoffset=1.2 - 0.0 + 0.1;
            tempx=downsample(Time,tempds);
            tempy=tempsize*downsample(ArousalOrig,tempds)+tempoffset;
            tempxx=[tempx;flipud(tempx)];
            tempyy=[tempy;tempoffset+zeros(length(tempy),1)];
            fill(tempxx,tempyy,[0.05 0.2 0.9],'EdgeColor','none','FaceAlpha',0.75);
            %fill(tempxx,tempyy,[0 1 0],'EdgeColor','none');
        end
        if sum(strcmp(ChannelsList,'EEG')==1)&&1
            EEGplot = SigT.EEG;
            plot(Time,1.1+0.1*EEGplot/(prctile(EEGplot,99.5)-prctile(EEGplot,0.05)),'color',[0.1 0.5 0.3],'linewidth',0.25);
        end
        
        if sum(strcmp(ChannelsList,'WakeSleep')==1)
            WakeSleep = SigT.WakeSleep;
            plot(Time,1.2+0.2*WakeSleep,'k');
        end
        if sum(strcmp(ChannelsList,'WPr')==1)
            WPr = SigT.WPr;
            ArPr = SigT.ArPr;
            WAPr = max([WPr ArPr]')';
            plot(Time,1.2+0.2*WAPr,'r');
        end
        if sum(strcmp(ChannelsList,'WPrB')==1)
            WPr = SigT.WPrB;
            ArPr = SigT.ArPrB;
            WAPr = max([WPr ArPr]')';
            plot(Time,1.2+0.2*WAPr,'r');
        end
        
        plotbreathgrid=1;
        flowplotpos=0.5;
        if plotbreathgrid
            gridheight=0.5;
            xgridlist = Time(BB_i_start)';
            plot([xgridlist;xgridlist],flowplotpos+gridheight*[-1 +1],'color',1-(1-0.9*[1 0.5 0.5]).^2);
            xgridlist = Time(BB_i_mid)';
            plot([xgridlist;xgridlist],flowplotpos+gridheight*[-1 +1],'color',1-(1-0.9*[1 0.5 0.5]).^4);
            %zero flow baseline:
            plot(plotflow_t,flowplotpos+0*plotflow_t,'color',0.8*[1 0.5 0.5]);
        end
        
        %back to plotting flow, RIP, SpO2 (over EEG)
        plot(plotflow_t,0.5+(plotflow-mflow)/nflow,'k',...
            Time,-0.25+(RIP_Thorax_plot-mTH)/nTH,'k',...)
            Time,-1+(RIP_Abdo_plot-mAB)/nAB,'k',...
            Time,-1.75+(spo2-mspo2)/nspo2,'r');
        plotpos = -1.75;
        
        set(gca,'YTick',[-1.75 -1 -0.25 0.5 1.2],...
            'XColor',[1 1 1],...
            'TickDir','out',...
            'TickLength',[0 0],...
            'FontName','Arial','FontSize',7);
        box('off');
        ylim([-2.35 1.6]);
        set(gca,'YDir','normal','FontSize',7);
        set(gca,'YTickLabel',{'SpO2', 'Abdomen', 'Thorax','Flow','EEG'});
        
        optionalchannelslist = {'Edi','Pes','GGpmax'};
        if 1 %clip Pes upper
            try
                Pes(Pes>Pesbaselineest+20)=NaN; %overwrite just for plot, careful not to analyze more below
            catch me
            end
        end
        for aa=1:length(optionalchannelslist)
            try
            channeltoplot = optionalchannelslist{aa};
            if sum(strcmp(ChannelsList,channeltoplot)==1)
                hold('on')
                tempupper=prctile(eval(channeltoplot),alpha_);
                templower=prctile(eval(channeltoplot),100-alpha_);
                tempdelta=tempupper-templower;
                tempupper=tempupper-tempdelta*0.1;
                templower=templower+tempdelta*0.1;
                tempdelta=tempupper-templower;
                tempmid = templower+tempdelta*0.5;
                %nX=1.5*(prctile(eval(channeltoplot),alpha_)-prctile(eval(channeltoplot),100-alpha_)); mX=nanmean(eval(channeltoplot));
                plotpos=plotpos-0.75;
                plot(plotflow_t,plotpos+(eval(channeltoplot)-tempmid)/tempdelta/2.2,'k')
                hold('on')
                temp = get(gca,'YTick');
                set(gca,'YTick',[temp(1)-0.75 temp]);
                temp = get(gca,'Ylim');
                set(gca,'Ylim',[temp(1)-0.75 temp(2)]);
                temp = get(gca,'YTickLabel');
                set(gca,'YTickLabel',[channeltoplot;temp]);
            end
            catch
            end
        end
        
        if 1
            %plot final event breaths over the top of flow in blue
            hold('on')
            dt_new=0.25;
            startT=ceil(BB_t(1));
            endT=ceil(BB_t(length(BB_t)));
            time_dt=startT:dt_new:endT;
            E1_rs=0*time_dt;
            for i=1:length(time_dt)
                E1_rs(i) = 1-E1(find(Time(BB_i_start)<=time_dt(i),1,'last'));
            end
            C1_rs=0*time_dt;
            for i=1:length(time_dt)
                C1_rs(i) = Ecentralapnea(find(Time(BB_i_start)<=time_dt(i),1,'last'));
            end
            E1_rs=E1_rs&(~C1_rs);
            if sum(E1_rs)>0
                tempsize=0.25; tempoffset=0.375;
                tempy=tempsize*(E1_rs)+tempoffset;
                tempxx=[time_dt fliplr(time_dt)];
                tempyy=[tempy tempoffset+zeros(1,length(tempy))];
                fill(tempxx,tempyy,[0.1 0.1 1],'EdgeColor','none','FaceAlpha',0.25);
            end
            hold('on')
            if sum(C1_rs)>0
                tempsize=0.25; tempoffset=0.375;
                tempy=tempsize*(C1_rs)+tempoffset;
                tempxx=[time_dt fliplr(time_dt)];
                tempyy=[tempy tempoffset+zeros(1,length(tempy))];
                fill(tempxx,tempyy,[1 0.1 0.1],'EdgeColor','none','FaceAlpha',0.25);
            end
            %figure(980); plot(time_dt,E1_rs,'r',BB_t,E1,'r.')
        end

        %% flow:drive plot
        ax1(3)=subplot('Position',[0.131 0.325 0.775 0.125]); hold on;
        FlowDrive = BreathDataTable{:, 'FlowDrive'};
        Apnea_B = BreathDataTable{:, 'ApneaB'};
        Apnea_B = logical(Apnea_B);
        FlowDriveOutPlot = FlowDrive*100;
        FlowDriveOutPlot(FlowDriveOutPlot>105)=105;
        FlowDriveOutPlot(FlowDriveOutPlot<0)=0;
        FlowDriveOutPlot2 = FlowDriveOutPlot;
        FlowDriveOutPlot2(Apnea_B)=0;
        FlowDriveOutPlot3 = FlowDriveOutPlot;
        FlowDriveOutPlot3(isnan(FlowDriveOutPlot3))=0;        
        plot(Time,100+0*Time,'k--');
        plot(Time,50+0*Time,'k--');
        plot(Time,0+0*Time,'k--');
        stairs(BB_t(1:i_end),FlowDriveOutPlot3, 'Color', [1 0 1]); % plot all NaN flow:drive in pink
        stairs(BB_t(1:i_end),FlowDriveOutPlot2, 'Color', [1 0 0]); % plot all NaN flow:drive apnea breaths in red
        stairs(BB_t(1:i_end),FlowDriveOutPlot, 'Color', [0 0 1]); % plot all flow:drive in blue solid
        xlim([min(Time) max(Time)]);
        ylim([-5 106]);
        set(gca,'box','off','tickdir','out', 'fontname','Arial','fontsize',7);
        set(gca,'XColor',[1 1 1]); % make x axis white (hide it)
        ylabel('Flow:Drive (%)');

        %% ventilation plot
        ax1(4)=subplot('Position',[0.131 0.100 0.775 0.21]); 
        hold('on');
        %plot event breaths over the top of flow in blue
        tempsize=0.5; tempoffset=0.75;
        tempy=tempsize*(E1_rs)+tempoffset;
        tempxx=[time_dt fliplr(time_dt)];
        tempyy=[tempy tempoffset+zeros(1,length(tempy))];
        fill(tempxx,tempyy,[0.1 0.1 1],'EdgeColor','none','FaceAlpha',0.25);
        tempsize=0.5; tempoffset=0.75; %updated 6/2014
        tempy=tempsize*(C1_rs)+tempoffset;
        tempxx=[time_dt fliplr(time_dt)];
        tempyy=[tempy tempoffset+zeros(1,length(tempy))];
        fill(tempxx,tempyy,[1 0.1 0.1],'EdgeColor','none','FaceAlpha',0.25);
        
        if 0
            try
                realVdrive = FlowEdi_VI/(meanVIbeforenormalizing);
                stairs(BB_t(1:i_end),realVdrive(1:i_end),'color',[1 0.6 0.6]);
                realVdrive2 = FlowPes_VI/(meanVIbeforenormalizing);
                stairs(BB_t(1:i_end),realVdrive2(1:i_end),'color',[0.6 1 1]);
            catch me
            end
        end
        
        stairs(BB_t(1:i_end),Ventilation_EventScoring_data(1:i_end,1)+Veupnea,'color',[0.4 0.4 0.7]);
        plot(BB_t(1:i_end),[zeros(1,length(BB_t(1:i_end))); ones(1,length(BB_t(1:i_end)))],'k:');
        stairs(BB_t(i_1:i_end)+0.3,Vdr_est(i_1:end)'+Veupnea+VAr_est(i_1:i_end)','color',[0.1 0.8 0.1]);
        plot(BB_t(i_1:i_end),Vdr_est(i_1:end)'+Veupnea,'color',[0 0 0]);
        
        % meanSEM = nanmean(upperSEM(i_1:end)-lowerSEM(i_1:end));
        
        %         if exist('lowerSEM')&&meanSEM<1&&0
        %         plot(BB_t(i_1:i_end),lowerSEM(i_1:end)'+Veupnea,'color',[0.6 0.6 0.6]);
        %         plot(BB_t(i_1:i_end),upperSEM(i_1:end)'+Veupnea,'color',[0.6 0.6 0.6]);
        %         end
        %plot(BB_t(i_1:i_end),W(i_1:end)','r');
        %plot(BB_t(i_1:i_end),VAr_est(i_1:i_end),'b');
        %xlabel(num2str(LGplusinfo(4),2))
        %LGplusinfo=[Gtot_est tau1_est tau2_est LG180_est T180_est LG60_est LG30_est delay_est VRA_est];
        if ~isnan(Arthres)
            arthrestext = [' ArThr=' num2str(Arthres,2)];
        else
            arthrestext=[];
        end
        %xlabel([ 'LG1=' num2str(LoopGainParameters(6),2) ' Tn=' num2str(LoopGainParameters(5),2) ' LGn=' num2str(LoopGainParameters(4),2) ' delay=' num2str(LoopGainParameters(3),2) arthrestext ' RsqTotal=' num2str(FitQual(2),2)])
        xlabel([ 'Pos=', num2str(position_mode,1), '(', num2str(Percent_position,0), ')',...
            ' meanSEM=', num2str(meanSEM,2),...
            ' LG1=', num2str(LoopGainParameters(6),2),'(', num2str(CI_parameters(1),2), ')',...
            ' Tn=', num2str(LoopGainParameters(5),2), ' LGn=', num2str(LoopGainParameters(4),2), ...
            ' delay=', num2str(LoopGainParameters(3),2), '(', num2str(CI_parameters(5),2), ')',...
            arthrestext,' RsqTotal=', num2str(FitQual(2),2)]);
        
        hold('off');
        set(gca,'FontSize',7,...
            'TickLength',[0.0 0.0],...
            'TickDir','out',...
            'xtick',[Time(1) max(Time)],'xticklabel',{datestr(Time(1)/86400,'HH:MM:SS'),datestr(max(Time)/86400,'HH:MM:SS')},...
            'FontName','Arial');
        box('off');
        
        xlim([min(Time) max(Time)]);
        
        linkaxes(ax1,'x');

        %% manual scoring touchups
        if settings.manualscoringtouchups
            disp('click left of data if no touch ups needed, or left then right to add events in range, or right then left to remove events in range');
            [t(1),~,~]=ginput(1);
            if t(1)<min(Time)
                manualscoringtouchupsincomplete=0;
            else
                [t(2),~,button]=ginput(1);
                manualscoringtouchupsincomplete=1;
                if button==1
                    if t(1)<t(2)
                        E1(BB_t>t(1)&BB_t<t(2))=0; %add events...
                    elseif t(1)>t(2)
                        E1(BB_t>t(2)&BB_t<t(1))=1; %remove events...
                    end
                elseif button>1 %right click
                    if t(1)<t(2)
                        Ecentralapnea(BB_t>t(1)&BB_t<t(2))=1; %add events...
                    elseif t(1)>t(2)
                        Ecentralapnea(BB_t>t(2)&BB_t<t(1))=0; %remove events...
                    end
                end
            end
        else
            manualscoringtouchupsincomplete=0;
        end
    else
        manualscoringtouchupsincomplete=0;
    end
    
    %%
    if manualscoringtouchupsincomplete==0
        break
    end
    
end %end while manualscoringtouchups loop

%% save loop gain plot
if settings.saveplots
    % old
    if 0
        savefig(gcf,[settings.OutputDataDirectory '\' settings.savename '_n=' num2str(n) '_w=' num2str(winNum)],'compact');
    else 
        % save in folder for each patient
        saveloc=[settings.OutputDataDirectory settings.savename filesep settings.filename];
        if ~(exist(saveloc, 'dir') == 7)
            mkdir(saveloc);
        end
        % save LGfigID in the file format specified.
        switch settings.savefigas
            case 'saveasTIFF'
                print(LGfigID, [saveloc, '\LG_window=',num2str(winNum)], '-dtiff', '-r300');               
            case 'saveasPNG'
                saveas(LGfigID, [saveloc, '\LG_window=',num2str(winNum)], 'png');                
            case 'saveasFIG' 
                savefig(LGfigID, [saveloc, '\LG_window=',num2str(winNum)],'compact');                
            % 2021-11-02 12:15:33 Tuesday EAS added this option to save both .fig and .png plots
            case 'saveasFIGandPNG'            
                LGfigID.WindowState = 'maximized';
                saveloc1 = [saveloc '\lgfigs'];
                saveloc2 = [saveloc '\lgpngs'];
                if ~(exist(saveloc1, 'dir') == 7)
                    mkdir([saveloc1]);
                end                    
                if ~(exist(saveloc2, 'dir') == 7)
                    mkdir([saveloc2]);
                end      
                savefig(LGfigID, [saveloc1, '\LG_window=',num2str(winNum)],'compact');
                % saveas(LGfigID, [saveloc, '\LG_window=',num2str(winNum)], 'png');
                % exportgraphics is better than saveas because there's no
                % blank space around the figure
                exportgraphics(LGfigID, [saveloc2, '\LG_window=',num2str(winNum),'.png']);
        end
    end
end

end
