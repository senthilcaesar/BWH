
criteria = Poscriteria==1 & nremXcriteria==1 & FREM<=maxFREM & FREM>=minFrem & TimeOfNightCriteria==1 & N_events>=minNevents & minwake<=maxwakethres;

sumnremXcriteria=sum(nremXcriteria==1 & FREM<=maxFREM & FREM>=minFrem);

SpO2VESyncData=[];
clear spo2_w VI_w E_w T_w Ar_w Er1_w hyp_w Er2_w
for w=1:W
    try
        if criteria(w)==0
            continue
        end
        TableTemp=BreathDataTable{1}{w};
        VI_w = TableTemp.VI;
        spo2_w = TableTemp.spo2;
        hyp_w = TableTemp.hypnog_B;
        if k>0
            spo2_w = [spo2_w(k+1:end); ones(k,1)*NaN ];
        elseif k<0
            spo2_w = [ones(-k,1)*NaN ; spo2_w(1:end+k)];
        end
        E_w = DataOut{1}{w}(:,6);
        Ar_w = DataOut{1}{w}(:,10);
        T_w = DataOut{1}{w}(:,1);
        Er1_w = DataOut{1}{w}(:,7);
        Er2_w = [NaN ; Er1_w(1:end-1)];
        Alpha=FVeupnea(spo2_w,spo2wakemean10to90);
        if 0
            VI_w=VI_w*Alpha;
        end
        SpO2VESyncData = [SpO2VESyncData;[spo2_w VI_w E_w T_w Ar_w Er1_w hyp_w Er2_w];[NaN NaN NaN NaN NaN NaN NaN NaN]];
    catch me
    end
end
spo2_w=SpO2VESyncData(:,1);
VI_w=SpO2VESyncData(:,2);
E_w=SpO2VESyncData(:,3);
T_w=SpO2VESyncData(:,4);
Ar_w=SpO2VESyncData(:,5);
Er1_w=SpO2VESyncData(:,6);
hyp_w=SpO2VESyncData(:,7);
Er2_w=SpO2VESyncData(:,8);
%%
if plotfig
    figure(12); clf(12);
    ax(1)=subplot(2,1,1);stairs(T_w,VI_w);
    hold('on');stairs(T_w,Ar_w);
    ax(2)=subplot(2,1,2);stairs(T_w,spo2_w);
    linkaxes(ax,'x');
end

criteria_ = spo2_w<prctile(spo2_w,50)&E_w==1;
criteria_ = Ar_w==1;
%criteria_ = Er1_w==1&Ar_w==1;
criteria_ = Er1_w==1&Ar_w==1;

%criteria = (0*spo2_w+1)>0; %great for anatomy?

%criteria = E_w==1; %not in events

xval = 100-spo2_w(criteria_);
yval = VI_w(criteria_);
Amatrix = [xval ];

if plotfig
    figure(11); clf(11);
    ax11(1)=subplot(3,1,[1 2]);
    temprand=0.4*(rand(sum(criteria_),1)-0.5);
    scatter(spo2_w(criteria_)+temprand,yval,8,[0.8 0.3 0.1],'filled','markerfacealpha',0.4);
end

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
    scatter(spo2_w(criteria_)+temprand,yval,8,[0.0 0.3 0.8],'filled','markerfacealpha',0.2);
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

%bindata
xbins = min(xval):1:max(xval);
clear VEmedianbin VEmedianbin_N
for i=1:length(xbins)
    VEmedianbin(i)=nanmedian(yval(xval==xbins(i)));
    VEmedianbin_N(i)=sum(xval==xbins(i));
end
VEmedianbin(VEmedianbin_N<10)=NaN;

SpO2arthres=nanmedian(100-xval);
    SpO2arthres_q1=prctile(100-xval,25);
    SpO2arthres_q3=prctile(100-xval,75);
    
if plotfig
    %plot(100-xbins,VEmedianbin);
    xvals2=min(100-xval):1:max(100-xval);
    ax11(2)=subplot(3,1,3);
    h=hist(100-xval,xvals2);
    h=h/nansum(h);
    bar(xvals2,h,'edgecolor','none','barwidth',1,'facealpha',0.5); box('off');
    hold('on');
    set(gca,'tickdir','out');
    plot(SpO2arthres*[1 1],[0 h(xvals2==SpO2arthres)],'k-');
    plot(SpO2arthres_q1*[1 1],[0 h(xvals2==SpO2arthres_q1)],'k--');
    plot(SpO2arthres_q3*[1 1],[0 h(xvals2==SpO2arthres_q3)],'k--');
    hold('off');
end