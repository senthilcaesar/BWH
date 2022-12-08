%SpO2DriveOne

plotfig=1;
%%




% if plotfig
%     
%     ax11(1)=subplot(3,1,[1 2]);
%     temprand=0.4*(rand(sum(criteria_),1)-0.5);
%     scatter(T1.spo2(criteria_)+temprand,yval,8,[0.8 0.3 0.1],'filled','markerfacealpha',0.4);
% end
%%

clear bSlope Rsq
kval=-3:14;
figure(896); clf(896);
criteria_ = T1.E1==1&T1.AR3==1; %not in events; and in arousal
temprand=0.4*(rand(sum(criteria_),1)-0.5);

for kk=1:length(kval)
    k=kval(kk);
    
    T1.SpO2shift = circshift(T1.spo2,-k);
    T1.Time0shift = circshift(T1.Time0,-k); 
    Iexclude = T1.Time0shift~=T1.Time0;
    T1.SpO2shift(Iexclude)=NaN;
    
    xval = 100-T1.SpO2shift(criteria_);
    yval = 100*T1.VI(criteria_);
    Amatrix = [xval];

    subplot(4,5,kk);
    J=2;
    for jj=1:J
        [b,~,stats]=glmfit(Amatrix,yval);
        [b stats.p stats.se];
        bSlope(kk)=b(2);
        yvalpred = Amatrix*b(2:end)+b(1);
        Error = abs(yvalpred-yval);
        I=Error>prctile(Error,100*(0.975^(1/J))); %will leave best 90% of data in four steps: 100*(0.9^(1/4))
        yval(I)=NaN;
    end
    
    
    if plotfig
        
        hold('on');
        scatter(100-xval+temprand,yval,8,[0.0 0.3 0.8],'filled','markerfacealpha',0.2);
    end
    [b,~,stats]=glmfit(Amatrix,yval);
    [b stats.p stats.se];
    bSlope(kk)=b(2);
    
    yvalpred = Amatrix*b(2:end)+b(1);
    if plotfig
        hold('on');
        plot(100-xval,yvalpred,'-','color',[0.1 0.3 0.8]);
    end
    
    SSres = nansum((yval - yvalpred).^2);
    SStot = nansum((yval - nanmean(yval)).^2);
    Rsq(kk) = 1-SSres/SStot;
    
    
    
    
    
    %SpO2arthres=nanmedian(100-xval);
    %SpO2arthres_q1=prctile(100-xval,25);
    %SpO2arthres_q3=prctile(100-xval,75);
    
%     if plotfig
%         %plot(100-xbins,VEmedianbin);
%         xvals2=min(100-xval):1:max(100-xval);
%         ax11(2)=subplot(3,1,3);
%         h=hist(100-xval,xvals2);
%         h=h/nansum(h);
%         bar(xvals2,h,'edgecolor','none','barwidth',1,'facealpha',0.5); box('off');
%         hold('on');
%         set(gca,'tickdir','out');
%         plot(SpO2arthres*[1 1],[0 h(xvals2==SpO2arthres)],'k-');
%         plot(SpO2arthres_q1*[1 1],[0 h(xvals2==SpO2arthres_q1)],'k--');
%         plot(SpO2arthres_q3*[1 1],[0 h(xvals2==SpO2arthres_q3)],'k--');
%         hold('off');
%     end
    
end

if plotfig
    figure(13)
    subplot(2,1,1); plot(kval,Rsq,'.-');
    subplot(2,1,2); plot(kval,bSlope,'.-'); 
end

[VESpO2slope,kkmax]=max(bSlope);

%% recalculate with best delay
 T1.SpO2shift = circshift(T1.spo2,-kval(kkmax));
    T1.Time0shift = circshift(T1.Time0,-kval(kkmax)); 
    Iexclude = T1.Time0shift~=T1.Time0;
    T1.SpO2shift(Iexclude)=NaN;
    
    xval = 100-T1.SpO2shift(criteria_);
    yval = 100*T1.VI(criteria_);
    Amatrix = [xval];

    J=2;
    for jj=1:J
        [b,~,stats]=glmfit(Amatrix,yval);
        [b stats.p stats.se];
        %bSlope(kk)=b(2);
        yvalpred = Amatrix*b(2:end)+b(1);
        Error = abs(yvalpred-yval);
        I=Error>prctile(Error,100*(0.975^(1/J))); %will leave best 90% of data in four steps: 100*(0.9^(1/4))
        yval(I)=NaN;
    end    
    T1.SpO2Drive = b(2)*(100-T1.SpO2shift) + b(1);
    
    %% Once Best delay is known, do fun things
    
    Iwake = criteria_; %criteria_ = T1.E1==1&T1.AR3==1; %not in events and in arousal
    figure(452); clf(452);
        hold('on');
        subplot(1,2,1);
        scatter(T1.SpO2Drive(Iwake),100*T1.VI(Iwake),5,[0.0 0.3 0.8],'filled','markerfacealpha',0.1);
        ylim([0 300]);
        subplot(1,2,2);
        scatter(T1.SpO2shift(Iwake),100*T1.VI(Iwake),5,[0.0 0.3 0.8],'filled','markerfacealpha',0.1);
        ylim([0 300]);

    Isleep = T1.hypnog_B<4 & T1.AR3==0; %in sleep and not in arousal
    figure(453); clf(453);
        hold('on');
        subplot(1,2,1);
        scatter(T1.SpO2Drive(Isleep),100*T1.VI(Isleep),5,[0.0 0.3 0.8],'filled','markerfacealpha',0.1);  hold('on');
        ylim([0 300]);
        subplot(1,2,2);
        scatter(T1.SpO2shift(Isleep),100*T1.VI(Isleep),5,[0.0 0.3 0.8],'filled','markerfacealpha',0.1);  hold('on');
        ylim([0 300]);

    %% bindata
%     xbins = prctile(T1.SpO2shift(Isleep),[0:10:100])
%     clear VEBin SpO2ShiftBin SpO2DriveBin NBin
%     for i=1:length(xbins)-1
%         Ii = T1.SpO2shift>=xbins(i) & T1.SpO2shift<xbins(i+1) & Isleep==1;
%         VEBin(i)=nanmedian(yval(Ii));
%         SpO2ShiftBin(i)=nanmedian(T1.SpO2shift(Ii));
%         SpO2DriveBin(i)=nanmedian(T1.SpO2shift(Ii));
%         NBin(i)=sum(T1.SpO2shift(Ii));
%     end
%       
%     to do, plot bin data over the scatter
%         subplot(1,2,1);
%         
%         scatter(T1.SpO2Drive(Isleep),100*T1.VI(Isleep),5,[0.0 0.3 0.8],'filled','markerfacealpha',0.1);
%         ylim([0 300]);
%         subplot(1,2,2);
%         scatter(T1.SpO2shift(Isleep),100*T1.VI(Isleep),5,[0.0 0.3 0.8],'filled','markerfacealpha',0.1);
%         ylim([0 300]);
%     if plotfig
%         hold('on');
%         plot(100-SpO2medianbin,VEmedianbin,'-','color',[0 0 0]);
%     end










