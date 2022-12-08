function [BreathSigT,BreathBox] = EventFlowShape(BreathDataTable2,Evts,SigT,settings)
% Generates ensemble average flow time series for a patient's events

% [BreathDataTable2,~,BreathFLDataTable2,~,BreathDataTableFulls]=GetNonOvlappedVE(BreathDataTable,BreathFLDataTable);
            
            %%
            dt = 1/settings.Fs;
            ploton=0;
            EvtEnd = Evts.RespT.EventEnd(Evts.RespT.EventCodes==4);
            
            EvtStart = Evts.RespT.EventStart(Evts.RespT.EventCodes==4);
            EvStIdx = nan(length(EvtStart),1);
            EvEndIdx = nan(length(EvtStart),1);
            Window = nan(length(EvtStart),1);
            
            for ii=1:length(EvtStart)
                if EvtStart(ii) < BreathDataTable2.Time_start(1) ||... %if outside whole table then skip
                        EvtStart(ii) > BreathDataTable2.Time_start(end)
                    excluderow = ii;
                    continue
                end
                
                % find event start index and time
                TimeDiff = BreathDataTable2.Time_start - EvtStart(ii);
                TimeDiff(TimeDiff<0 | BreathDataTable2.FDuplicated2~=0) = nan;
                [~,EvStIdx(ii,1)] = min(TimeDiff);
                
                Windowtemp = BreathDataTable2.Time0(EvStIdx(ii,1));
                Window(ii)=Windowtemp;
                
                % find event end index and time - FIX THIS NEEDS TO BE IN SAME
                % WINDOW AS ABOVE
                TimeDiff = BreathDataTable2.Time_end - EvtEnd(ii);
                TimeDiff(TimeDiff>0 | BreathDataTable2.Time0 ~= Windowtemp) = nan;
                [~,EvEndIdx(ii,1)] = max(TimeDiff);
                if isnan(max(TimeDiff))
                    EvEndIdx(ii,1)=NaN;
                end
                
            end
            
            EvEndIdx(isnan(EvEndIdx))=[];
            
            %% Actual Plot
            Nbreaths=15;
            %to do, handle errors if there is not enough data to make a
            %breath
            for M=1:(Nbreaths*2+1)
                ilocal = EvEndIdx-Nbreaths+M-1; %breath (row) selection 
                ilocal(ilocal>height(BreathDataTable2))=[]; %remove if it doesn't exist
                ilocal(ilocal<1)=[]; %remove if it doesn't exist
                tempT = BreathDataTable2(ilocal,:);
                
                %to do, add check that the window is the same as that for
                %the end-event window
                
                Tbreath = nanmedian(tempT.Ttot);
                Tbreathi = round(Tbreath/dt);
                TbreathInsp(M) = nanmedian(tempT.Ti);
                 BreathSigT.TbreathInspi(M) = round(TbreathInsp(M)/dt);
                
                TbreathExp(M) = nanmedian(tempT.Te);
                 BreathSigT.TbreathExpi(M) = round(TbreathExp(M)/dt);
                
                invertflow=0; %fix **************************************
                
                SigIn = ((-1)^settings.Pnasaldownisinsp)*((-1)^invertflow)*SigT.Flow; %%% added invert flow
                %SigIn = DataEventHypnog_Mat(:,find(strcmp(ChannelsList,'SpO2')==1));
                BreathBox{M}=nan(height(tempT), BreathSigT.TbreathInspi(M)+ BreathSigT.TbreathExpi(M));
                
                if ploton
                    figure(98); clf(98)
                end
                
                for i=1:height(tempT)
                    
                    Time_starti=round((tempT.Time_start(i)-SigT.Time(1))/dt)+1;
                    I = Time_starti:(Time_starti+round(tempT.Ttot(i)/dt)-1);
                    %I = Time_starti:(Time_starti+round(tempT.Ti(i)/dt)-1);
                    
                    if 1
                        exponent = settings.scalingexponent;
                        leak = tempT.leak(i);
                        leak2 = tempT.leak2(i);
                        IEratio = tempT.IEratio(i);
                        Vflow_out2 = SigIn(I)-leak;
                        Vflow_out2(Vflow_out2>0)=(Vflow_out2(Vflow_out2>0).^(exponent))/(IEratio^0.5);
                        Vflow_out2(Vflow_out2<0)=(-(-Vflow_out2(Vflow_out2<0)).^(exponent))*(IEratio^0.5);
                        Vflow_out2 = Vflow_out2-leak2;
                        Vflow_out2 = Vflow_out2./tempT.Veup(i); %normalize the whole flow signal by local eupnea
                        sigout = Vflow_out2;
                    else
                        sigout=SigIn(I); %single breath linearized flow
                    end
                    
                    
                    Xin2 = 0:1:([tempT.Ti(i)/dt]-1);
                    Xin = Xin2 / (([tempT.Ti(i)/dt-1])/( BreathSigT.TbreathInspi(M)-1));
                    Xin3 = 0:1:(BreathSigT.TbreathInspi(M)-1);
                    
                    sigout2 = interp1(Xin,sigout(1:length(Xin)),Xin3); %insp
                    
                    Xin2 = 0:1:(tempT.Te(i)/dt-1);
                    Xin = Xin2 / (([tempT.Te(i)/dt-1])/( BreathSigT.TbreathExpi(M)-1));
                    Xin3 = 0:1:(BreathSigT.TbreathExpi(M)-1);
                    
                    sigout3 = interp1(Xin,sigout(end-length(Xin)+1:end),Xin3); %exp
                    
                    sigoutS = [sigout2 sigout3]; %joined
                    if ploton
                        figure(98)
                        subplot(3,1,1);
                        plot([sigout2 sigout3])
                        hold on
                    end
                    BreathBox{M}(i,1:length(sigoutS))=(sigoutS(:))'; %keep it
                end
                
                %removes outliers (e.g. 3 times); edit to make sure there
                %are at least X breaths (e.g. X=10)
                for jj=1:3
                    if ploton
                        figure(98)
                        subplot(3,1,2);
                        plot(nanmean(BreathBox{M}))
                        hold on
                    end
                    temp = BreathBox{M};
                    temp2 = temp - nanmean(temp);
                    %temp2 = temp - prctile(abs(temp),25);
                    %temp3 = nanmedian(abs(temp2),2);
                    temp3 = nanmean(abs(temp2),2);
                    temp3 = nanmean(temp2.^2,2);
                    temp3 = max(abs(temp2)')';
                    [~,sort_] = sort(temp3,'descend');
                    percentile = min([50]); %use this or use some other value that would leave you with at least X breaths.
                    BreathBox{M}(temp3>prctile(temp3,percentile),:)=NaN; %select on 50%
                    %best breaths
                    %BreathBox{M}(sort_<31,:)=NaN; %find 30 best events
                    if ploton
                        figure(98)
                        % plot(nanmean(BreathBox{M}));
                        subplot(3,1,3);
                        plot(BreathBox{M}');
                        if 0
                            pause(0.5)
                        end
                    end
                end
            end
            
            sig=[];
            lengthBreathBox = zeros(Nbreaths*2+1,1);
            for M=1:(Nbreaths*2+1)
                
                %     sig_breath = nanmean(BreathBox{M}); % temp hack, can be removed
                %     sig_breath(isnan(sig_breath))=[];
                sig = [sig nanmean(BreathBox{M})];
                lengthBreathBox(M) = length(BreathBox{M});
            end
            ilastevent = sum(lengthBreathBox(1:(Nbreaths+1))); %add 1 or 2?
            
            BreathSigT.Time=[0:dt:(length(sig)-1)*dt]'-ilastevent*dt;
            BreathSigT.Flow=sig(:);
            if 1 %fillnan
                I = ~isnan(BreathSigT.Flow);
                BreathSigT.Flow=interp1(BreathSigT.Time(I),BreathSigT.Flow(I),BreathSigT.Time,'pchip');
            end
            
            figure(15);
            hold off;
            try 
                xlims = get(gca,'xlim');
            end
            plot(BreathSigT.Time,BreathSigT.Flow,'k');
            set(gca,'box','off');
            try 
                set(gca,'xlim',xlims);
            end

            