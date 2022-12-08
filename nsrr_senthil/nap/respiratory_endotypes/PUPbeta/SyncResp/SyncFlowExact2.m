function Lagaverage2=SyncFlowExact(FlowSpikeSync,FlowNoxSync,TimeNoxRIP1,Rthres)

%Later use Inductance_Thorax = interp1(TimeNoxRIP1+Lagaverage2,Inductance_Thorax,TimeNoxRIP1,'linear','extrap');
dtNox = TimeNoxRIP1(2)-TimeNoxRIP1(1);

ploton=1;
ploton2=1;
if ploton
  figure(21);
        ax21(1)=subplot(2,1,1); plot(TimeNoxRIP1,FlowSpikeSync); set(gca,'xtick',[]);
        ax21(2)=subplot(2,1,2); plot(TimeNoxRIP1,FlowNoxSync); 
        linkaxes(ax21,'x');
end

%dtNox = TimeNoxRIP1(2)-TimeNoxRIP1(1);
N2 = length(TimeNoxRIP1);


        %Rthres=0.25;

%% Xcorr1
        maxlags=round(1/dtNox)*60; %60
        Twindow=300; %180
       
        
        %[c,lags]=xcorr(VolNoxSync,RIPvol,maxlags);
        if ploton2
        [c,lags]=crosscorr(FlowNoxSync,FlowSpikeSync,maxlags); %positive "lag" indicates Spike signal is to the right.
        [maxcT,inxT]=max(c);
        figure(2); plot(lags*dtNox,c); hold('on');
        leadi=inxT-(maxlags+1);
        end
        leadi=0;
        %% max range for lags
        
        clear centertime maxc inx
        %bestlag_ = (inxT+lags(Lilag)-1)*dtNox;
        
        %maxlags = abs(lagmaxT);
        
        %M=round(N2/25/300); %600 is 10 min blocks
        %segmentlength=round(N2/M);
        
        Foverlap=0.75;
        
        Nlength = round(Twindow/dtNox);
        Nstep = round(Nlength*(1-Foverlap));
        M = floor((N2-Nlength)/Nstep);        
        
        maxlags2=max([abs(maxlags+leadi) abs(maxlags-leadi)]);
        includelagslist=(-maxlags:maxlags)+leadi;
        includelagslisti=includelagslist+maxlags2+1;
        %actuallagslist=-maxlags2:maxlags2;
        
        inx=NaN*(zeros(1,M));
        maxc=NaN*(zeros(1,M));
        centertime=NaN*(zeros(1,M));
        figure(2); clf(2);
        for m=2:M-1
            Li = 1+(m-1)*Nstep;
            Ri = Li+Nlength-1;
            
            if Ri>N2, Ri=N2; end
            xdata = FlowSpikeSync(Li:Ri);
            ydata = FlowNoxSync(Li:Ri);
            centertime(m) = mean(TimeNoxRIP1(Li:Ri));
            [c,lags]=crosscorr(ydata,xdata,maxlags2);
            temp=NaN*c; temp(includelagslisti)=c(includelagslisti); c=temp;
                [maxc(m),inx(m)]=max(c);
            if inx(m)>1&&inx(m)<length(c)
                [inx(m),maxc(m)] = PeakFitQuadratic(((inx(m)-1):(inx(m)+1)),c((inx(m)-1):(inx(m)+1)));
            end
            if ploton
            figure(2); 
            subplot(3,1,1); 
            try
                plot(lags*dtNox,c,'color',[0.5*m/M 0.5*m/M 0.5*m/M]); %hold('on');
                ylim([-0.5 1]);
                axx(1)=subplot(3,1,2);
                ydata_=(inx+lags(1)-1)*dtNox;
                plot(centertime,ydata_,'.'); %hold('on');
                ylim([prctile(ydata_,1)-0.01 prctile(ydata_,99)+0.01]);
                axx(2)=subplot(3,1,3);
                plot(centertime,maxc,'.'); %hold('on');
            catch me
            end
            %pause(0.00001)
            end
        end
        I=isnan(maxc);
        maxc(I==1)=[];
        inx(I==1)=[];
        centertime(I==1)=[];
        
        %figure(3); plot(centertime,maxc)
        %m=1; P1 = -2.1017e-05; -2.1200e-05; -1.9884e-05; -2.0469e-05;
        %m=3; P1 = -2.4848e-05;
        bestlag = (inx+lags(1)-1)*dtNox;    
        maxc2 = maxc;
        
        bestlag2 = bestlag;
        centertime2 = centertime;
        bestlag2(maxc<Rthres)=[];
        centertime2(maxc<Rthres)=[];
        maxc2(maxc<Rthres)=[];
        
           
%% fit linear time loss
        
        
        subplot(3,1,2); 
        plot(centertime2,bestlag2,'.');
        
        [P,~] = polyfit(centertime2,bestlag2,1);
        bestlagmodel = polyval(P,centertime2);
        
        subplot(3,1,2); plot(centertime2,bestlag2,'.',centertime2,bestlagmodel,'r:');
        ylabel('Delay, Nox versus Spike, (s)')
        
        subplot(3,1,3); plot(centertime2,maxc2,'.');
        ylabel('Correlation')
        
        %% fit moving-time-average time loss
        maxt = 300;
            Lagaverage=NaN*centertime2;
            for j=1:length(centertime2)
                temp2 = abs(centertime2(j)-centertime2); %list of times
                weights = (maxt-temp2).*(maxc2); weights(weights<0)=0; 
                weights = weights/sum(weights);
                Lagaverage(j)=nansum(bestlag2(weights>0).*weights(weights>0));
            end
            Lagaverage1=0*TimeNoxRIP1;
            if length(centertime2)>1
                Lagaverage1 = interp1(centertime2,Lagaverage,TimeNoxRIP1,'linear');
            end
                Lagaverage1(TimeNoxRIP1>=centertime2(end))=Lagaverage(end);
                Lagaverage1(TimeNoxRIP1<=centertime2(1))=Lagaverage(1);
                
        
        subplot(3,1,2); 
        hold('on');
        plot(TimeNoxRIP1,Lagaverage1,'g--');
      
        TimeNoxRIP4=TimeNoxRIP1+Lagaverage1;
        %confirminorder = sum(sort(TimeNoxRIP4)-TimeNoxRIP4)==0; %must be zero
        
        
        %% BSR1
        
        if length(centertime2)>1
            Nbp = [-1:20];
            clear Fsse SSE
            Fsse = NaN*Nbp;
            SSE = NaN*Nbp;
            for i=1:length(Nbp)
                try
                if Nbp(i)>0
                    ab = BrokenStickRegression(centertime2,bestlag2,Nbp(i));
                    if sum(isnan(sum(ab,2)))>1
                        predy = bestlag2*NaN;
                    else
                        %ab(isnan(sum(ab,2)),:)=[]; %catch for when BrokenStickRegression yields NaN for yvalues of a breakpoint
                        predy = interp1(ab(:, 1),ab(:, 2),centertime2,'linear','extrap');
                    end
                elseif Nbp(i)==0        
                    [P,S] = polyfit(centertime2,bestlag2,1);
                    predy = polyval(P,centertime2);
                elseif Nbp(i)==-1
                    predy = nanmean(bestlag2)+0*centertime2; 
                end
               
                SSE(i) = sum((predy-bestlag2).^2);
                CumminSSE(i) = min(SSE);
                if i>1
                    Fsse(i)=SSE(i)/CumminSSE(i-1);
                end
                catch me
                end 
            end
            
            i=find(Fsse<0.8&~isnan(Fsse),1,'last');
            Nbp=Nbp(i);
            if isempty(Nbp), Nbp=-1; end
            
                if Nbp>0
                    ab = BrokenStickRegression(centertime2,bestlag2,Nbp);
                        ab(isnan(sum(ab,2)),:)=[]; %catch for when BrokenStickRegression yields NaN for yvalues of a breakpoint
                    predy = interp1(ab(:, 1),ab(:, 2),centertime2,'linear','extrap');
                elseif Nbp==0        
                    [P,~] = polyfit(centertime2,bestlag2,1);
                    predy = polyval(P,centertime2);
                elseif Nbp==-1
                    predy = nanmean(bestlag2)+0*centertime2;
                end
            Lagaverage2 = interp1(centertime2,predy,TimeNoxRIP1,'linear','extrap');
        else
            Lagaverage2 = Lagaverage1;
        end

        
        subplot(3,1,2); 
        hold('on');
        plot(TimeNoxRIP1,Lagaverage2,'r-');
        
        