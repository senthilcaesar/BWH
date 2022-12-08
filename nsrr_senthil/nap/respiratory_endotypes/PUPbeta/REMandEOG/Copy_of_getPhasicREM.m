function [DataEventHypnog_Mat,ChannelsList,ChannelsFs]=getPhasicREM(DataEventHypnog_Mat,ChannelsList,ChannelsFs)
% tic
%             Options = [];
%             %ArSignal = DataEventHypnog_MatT.EventsAr;
%             %ArSignal = logitinverse(WSpredlogit)>0.25;
%             ArSignal = 1*(DataEventHypnog_MatT.EventsAr==1 | logitinverse(WSpredlogit)>0.75);
%             
%             widthrms=120;
%             BufferI = buffer(1:length(WSpredlogit),round(widthrms*1/dT),round((widthrms-0.5)*1/dT),'nodelay');
%                     BufferI(:,sum(BufferI==-8)>0)=[];       
%             t = DataEventHypnog_MatT.Time(BufferI);% buffer(DataEventHypnog_MatT.Time,round(widthrms*1/dT),round((widthrms-0.5)*1/dT),'nodelay');
%             t_ = t(1,:)' + widthrms/2;
%             WSslow = nanmedian(WSpredlogit(BufferI))';
%             WSslow = interp1(t_,WSslow,DataEventHypnog_MatT.Time,'linear');
%             
%             ArSignal = 1*(DataEventHypnog_MatT.EventsAr==1 | logitinverse(WSpredlogit)>0.75) | WSslow>0.75;
%             
            %%
            addpath(genpath('G:\Dropbox (Personal)\PUPbeta_git\PUPbeta20190629'))
            
            if exist('DataEventHypnog_Mat')
                DataEventHypnog_MatT = array2table(DataEventHypnog_Mat);
            end
            
            DataEventHypnog_MatT.Properties.VariableNames = ChannelsList;
            %clear DataEventHypnog_Mat;
            dT = DataEventHypnog_MatT.Time(2)-DataEventHypnog_MatT.Time(1);
            
            ArSignal = 1*(DataEventHypnog_MatT.EventsAr==1);
            
            Exclude = DataEventHypnog_MatT.SpO2<50 | DataEventHypnog_MatT.Epochs>4 | isnan(DataEventHypnog_MatT.Epochs);   
            
            Options = [];
            
            load mdlREM
            [PrREM,Tbl,EOGrmsConj_p50,ymtaAr180b,t_,rmssumbyrbyar,rmssumbyr,rmssum,Thres,LOCfilt,ROCfilt] = REMfromEOG(DataEventHypnog_MatT.LOC,DataEventHypnog_MatT.ROC,DataEventHypnog_MatT.Time,ArSignal,Exclude,Options,mdlREM);
            
            %%
            n=1; 
            Exp=1;
            plotfig1=1;
            
            if 0
                PrREM(NoiseBinary==1)=NaN;
            end
            
            if 1
                ExcludeREM = DataEventHypnog_MatT.SpO2<50 | DataEventHypnog_MatT.EventsAr==1 | DataEventHypnog_MatT.Epochs<0 | DataEventHypnog_MatT.Epochs>3 | isnan(DataEventHypnog_MatT.Epochs);
                REMvsNREM = (DataEventHypnog_MatT.Epochs==3)*1;
                REMvsNREM(ExcludeREM)=NaN;
                balance = nansum(REMvsNREM)/sum(~isnan(REMvsNREM));
            end
                        
            % Plot
            if 1
                Tbl.Time = DataEventHypnog_MatT.Time;
                Tbl.Subj = ones(length(REMvsNREM),1)*n+Exp*10000;
                Tbl.REMvsNREM = REMvsNREM;
                Tbl.ExcludeREM = ExcludeREM;
                %Tbl.WSslowPr = logitinverse(WSslow);
                %Tbl.WSslowLogOdds = WSslow;
            end
            %%
            %ThresBackup=Thres;
            %Thres=Thres+2;
            %%
            ThresF = 100; %1
            FthresForNTonic=25; %25
            %TonicREM_ = rmssumbyr<ThresF*(10.^Thres) & rmssum<(ThresF*50*(10.^Thres));
            PhasicREM_ = DataEventHypnog_MatT.Epochs==3 & rmssumbyr.*(1-ArSignal)>ThresF*(10.^Thres);
            TonicREM_ = DataEventHypnog_MatT.Epochs==3 & ~PhasicREM_ & ~ArSignal;
            
            TotalPowerSurgeInTonicREM = TonicREM_ & rmssum>ThresF*(10.^Thres)*FthresForNTonic;
            PercentREMinartifact = sum(TotalPowerSurgeInTonicREM)/sum(DataEventHypnog_MatT.Epochs==3)*100
            
            TonicREM_ = DataEventHypnog_MatT.Epochs==3 & ~PhasicREM_ & ~TotalPowerSurgeInTonicREM;
            
            %TonicREM = DataEventHypnog_MatT.Epochs==3 & ~PhasicREM;
            %PhasicREM = DataEventHypnog_MatT.Epochs==3 & (1-RemoveShortSegments(1-PhasicREM_,2,dT,0));
            if 0
            TonicREM = DataEventHypnog_MatT.Epochs==3 & (RemoveShortSegments(TonicREM_,3,dT,0));
            
            PhasicREM = DataEventHypnog_MatT.Epochs==3 & ~TonicREM;
            PhasicREM = (RemoveShortSegments(PhasicREM,2,dT,0));
            TonicREM = TonicREM & ~ArSignal;
            TonicREM = (RemoveShortSegments(TonicREM,6,dT,0));
            end
            
            temp = DataEventHypnog_MatT.Epochs==3 & (RemoveShortSegments(1-PhasicREM_,3,dT,0));
            %First we remove gaps between phasicREM* of <3 s (joins adjacent phasic REM)
            PhasicREM = (1-temp);
            PhasicREM = (RemoveShortSegments(PhasicREM,2,dT,0));
            %Then remove phasic REM <2 s.
            %Remove arousals again if needed.
            PhasicREM = PhasicREM & ~ArSignal & DataEventHypnog_MatT.Epochs==3; %just in case
            
            %Tonic REM is defined by absence of phasic REM, absence of arousal, and must be longer than 3 s. 
            TonicREM = TonicREM_ & ~PhasicREM & ~ArSignal;
            TonicREM = RemoveShortSegments(TonicREM,3,dT,0);
            
            %Find phasic REM, remove gaps between phasic of <3 sec.
            %Remove phasic REM <2 secs
            
            %Find tonic REM, remove short tonic <3 sec 
            %Remove tonic scored in arousal
            %remove tonic <6 seconds
            %
            if plotfig1
                figure(1); clf(1); set(gcf,'color',[1 1 1]);
                ax(1)=subplot(4,1,1);
                plot(DataEventHypnog_MatT.Time,[DataEventHypnog_MatT.LOC LOCfilt]);
                set(gca,'xtick',[],'box','off','xcolor',[1 1 1]);
                pos = get(gca,'position'); pos2 = [pos(1)*0.5 pos(2) pos(3)*1.15 pos(4)*1.2]; set(gca,'position',pos2);
                ax(2)=subplot(4,1,2);
                plot(DataEventHypnog_MatT.Time,[DataEventHypnog_MatT.ROC ROCfilt]);
                set(gca,'xtick',[],'box','off','xcolor',[1 1 1]);
                pos = get(gca,'position'); pos2 = [pos(1)*0.5 pos(2) pos(3)*1.15 pos(4)*1.2]; set(gca,'position',pos2);
                ax(3)=subplot(4,1,3);
                plot(DataEventHypnog_MatT.Time,DataEventHypnog_MatT.Epochs);
                hold on
                plot(DataEventHypnog_MatT.Time,[DataEventHypnog_MatT.EventsAr] + 4.5); 
                %plot(DataEventHypnog_MatT.Time,logitinverse(WSpredlogit) + 4.5);
                %plot(DataEventHypnog_MatT.Time,logitinverse(WSslow) + 4.5,'color','r');
                plot(DataEventHypnog_MatT.Time,[TonicREM_+5.9 TonicREM+6 PhasicREM_+7.4 PhasicREM+7.5]);
                %plot(DataEventHypnog_MatT.Time,[ArSignal*1+8.75]);
                plot(Tbl.Time,PrREM,'linewidth',1.5);
                set(gca,'xtick',[],'box','off','xcolor',[1 1 1]);
                pos = get(gca,'position'); pos2 = [pos(1)*0.5 pos(2) pos(3)*1.15 pos(4)*1.2]; set(gca,'position',pos2);
                ax(4)=subplot(4,1,4);
                %all signals on sqrt scaling
                plot(DataEventHypnog_MatT.Time,[rmssum.^0.5 -rmssum.^0.5],'color',[0.6 0.6 0.6]);
                hold on
                plot(DataEventHypnog_MatT.Time,[-rmssumbyr.^0.5],'color',[0.2 0.2 0.2]);
                plot(DataEventHypnog_MatT.Time,[rmssumbyrbyar.^0.5],'color',[0.9 0.1 0.2]);
                plot(DataEventHypnog_MatT.Time,[rmssumbyr.^0.5].*(1-ArSignal),'g:');
                plot(DataEventHypnog_MatT.Time,[(ThresF*(10.^Thres)) (10.^ymtaAr180b) ThresF*(10.^Thres)*25].^0.5,'color',[0.4 0.3 0.2]);
                pos = get(gca,'position'); pos2 = [pos(1)*0.5 pos(2) pos(3)*1.15 pos(4)*1.2]; set(gca,'position',pos2);
                set(gca,'box','off','tickdir','out');
                %         figure(8)
                %         plot(-1:0.01:1,logitinverse(((-1:0.01:1)-0.5)*10))
                linkaxes(ax,'x');
            end
            
            %%
            
            
            DataEventHypnog_MatT.TonicREM = TonicREM;
            DataEventHypnog_MatT.PhasicREM = PhasicREM;            
            DataEventHypnog_MatT.EOGconjugate = rmssumbyr;
            DataEventHypnog_MatT.EOGconjugateFThres = log10(rmssumbyr ./ (ThresF*(10.^Thres)) );
            
            ChannelsFs = [ChannelsFs;repmat(ChannelsFs(1),4,1)];
            ChannelsList = DataEventHypnog_MatT.Properties.VariableNames;
            DataEventHypnog_Mat = DataEventHypnog_MatT{:,:};
            
            