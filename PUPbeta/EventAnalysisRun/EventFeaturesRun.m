function [Ftrs]=EventFeaturesRun(Boxes,Ensembles,RespT,TotalAHI,Models)

%% Quantify different patterns of sleep apnea (Sinusoidal, Exponential, Square)
ploton=0;
[~,~,~,~,T1] = WaveformBreakdown(Ensembles.VI,Ensembles.Time,Models,ploton);
I = find(Ensembles.Time>-1.67*T1&Ensembles.Time<1.67*T1);
[Ftrs.VE_FsinEst,Ftrs.VE_FexpEst,Ftrs.VE_FsquEst,~,~] = WaveformBreakdown(Ensembles.VI(I),Ensembles.Time(I),Models,ploton);


%% Quantify inter-interval variability between events
PostEv=Boxes.EventsResp(:,Ensembles.Time>=1 & Ensembles.Time<=1.75*T1);
PostEv=PostEv';
TnextEv=nan(size(Boxes.EventsResp,1),1);
for ii=1:size(PostEv,2)
    IdxUp=find([diff(PostEv(:,ii))<0; 0]);
    if ~isempty(IdxUp)
        TnextEv(ii)=IdxUp(1);
    end
end

Ftrs.Ev_TnextEv=nanmean(TnextEv);
Ftrs.Ev_NnextEv=sum(~isnan(TnextEv));
Ftrs.Ev_PnextEv=sum(~isnan(TnextEv))/length(TnextEv);
Ftrs.Ev_CV=100*nanstd(TnextEv)/nanmean(TnextEv);

%% Quantify sleep apnea burdens
% plotfig=0;
% %[VeFtrs.OSBurden,VeFtrs.MinMax]=FindBurden(VE1,VE_mean,TimeEnsemble,TotalAHI,plotfig);
% [VeFtrs.Burden]=FindBurden(repmat(max(VE1),size(VE1,1),1)-VE1,max(VE_mean)-VE_mean,TimeEnsemble,TotalAHI,plotfig);
% [SatFtrs.Burden,SatFtrs.MinMax]=FindBurden(100-Sat1,100-Sat_mean,TimeEnsemble,TotalAHI,plotfig);
% [HRFtrs.Burden,HRFtrs.MinMax]=FindBurden(HR1,HR_mean,TimeEnsemble,TotalAHI,plotfig);
% [ArFtrs.Burden,ArFtrs.MinMax]=FindBurden(AR1,AR_mean,TimeEnsemble,TotalAHI,plotfig);
% [ArScFtrs.Burden,ArScFtrs.MinMax]=FindBurden(ArScore1,ArScore_mean,TimeEnsemble,TotalAHI,plotfig);
% [ArIntFtrs.Burden,ArIntFtrs.MinMax]=FindBurden(ArInt1,ArInt_mean,TimeEnsemble,TotalAHI,plotfig);

%% Quantify periodicity of VE, Sat, and ArInt
plotfig=0;
[VE_Periodicity]=QuantifyPeriodicity(Ensembles.VI,Ensembles.Time,T1,plotfig);
mergestructs = @(x,y) cell2struct([struct2cell(x);struct2cell(y)],[fieldnames(x);fieldnames(y)]);
Ftrs = mergestructs(Ftrs,VE_Periodicity);
end


function [Burden,MinMax]=FindBurden(EventMatrix,EventMean,TimeEnsemble,TotalAHI,plotfig)

    Burden.PreEvMinBsln=NaN;
    Burden.PostEvMinBsln=NaN;
    Burden.PreEvAvgBsln=NaN;
    Burden.NoBsln=NaN;
    Burden.AllEvsMinBsln=NaN; %The legendary hypoxic burden
    MinMax.maxAtZero=NaN;
    MinMax.maxAtZeroT=NaN;
    MinMax.minPreEvt=NaN;
    MinMax.minPreEvtT=NaN;
    MinMax.minPostEvt=NaN;
    MinMax.minPostEvtT=NaN;
    MinMax.maxPreEvt=NaN;
    MinMax.maxPreEvtT=NaN;
    MinMax.maxPostEvt=NaN;
    MinMax.maxPostEvtT=NaN;

    if sum(isnan(EventMean))<10
        [maxS,minS] = peakdetOriginal(EventMean,0.1,TimeEnsemble);

        %%%Max response occurs somewhere between event end -20 and event end
        %%%+20 seconds
        MaxWin=20;
        if ~isempty(maxS) & ~isempty(minS)
            max_resp=maxS(maxS(:,1)>=-MaxWin & maxS(:,1)<=MaxWin,:);
            max_resp=max_resp(max_resp(:,2)==max(max_resp(:,2)),:);

            if ~isempty(max_resp)
                MinMax.maxAtZero=max_resp(1,2);
                MinMax.maxAtZeroT=max_resp(1,1);

                min_pre=minS(minS(:,1)<max_resp(1,1),:);
                if ~isempty(min_pre)
                    min_pre=min_pre(end,:);
                    MinMax.minPreEvt=min_pre(1,2);
                    MinMax.minPreEvtT=min_pre(1,1);
                end

                min_post=minS(minS(:,1)>max_resp(1,1),:);
                if ~isempty(min_post)
                    min_post=min_post(1,:);
                    MinMax.minPostEvt=min_post(1,2);
                    MinMax.minPostEvtT=min_post(1,1);
                end

                max_pre=maxS(maxS(:,1)<max_resp(1,1),:);
                if ~isempty(max_pre)
                    max_pre=max_pre(end,:);
                    MinMax.maxPreEvt=max_pre(1,2);
                    MinMax.maxPreEvtT=max_pre(1,1);
                end

                max_post=maxS(maxS(:,1)>max_resp(1,1),:);
                if ~isempty(max_post)
                    max_post=max_post(1,:);
                    MinMax.maxPostEvt=max_post(1,2);
                    MinMax.maxPostEvtT=max_post(1,1);
                end

                if plotfig
                    figure(14);plot(TimeEnsemble,EventMean);hold on;
                    plot(min_pre(:,1),min_pre(:,2),'*',min_post(:,1),min_post(:,2),'*');
                    plot(max_resp(:,1),max_resp(:,2),'^',max_pre(:,1),max_pre(:,2),'^',max_post(:,1),max_post(:,2),'^');
                    hold off;
                end

                if ~isempty(min_pre) && ~isempty(min_post)
                    MeanEvResp=EventMean(TimeEnsemble>=min_pre(1,1) & TimeEnsemble<=min_post(1,1));

                    %%%Referencing with pre-event minimum
                    MeanEvResp1=MeanEvResp-min_pre(1,2);
                    MeanEvResp1(MeanEvResp1<0)=0;

                    %%%Referencing with post-event minimum
                    MeanEvResp2=MeanEvResp-min_post(1,2);
                    MeanEvResp2(MeanEvResp2<0)=0;

                    %%%Referencing with pre-event average
                    MeanEvResp3=MeanEvResp-mean(EventMean(TimeEnsemble<min_pre(1,1)));
                    MeanEvResp3(MeanEvResp3<0)=0;

                    Burden.PreEvMinBsln=nansum(MeanEvResp1)*TotalAHI/60;
                    Burden.PostEvMinBsln=nansum(MeanEvResp2)*TotalAHI/60;
                    Burden.PreEvAvgBsln=nansum(MeanEvResp3)*TotalAHI/60;
                    Burden.NoBsln=nansum(MeanEvResp)*TotalAHI/60;

                    AllEvResp=EventMatrix(TimeEnsemble>=min_pre(1,1) & TimeEnsemble <= min_post(1,1),:);
                    AllBsLns=repmat(min(AllEvResp),size(AllEvResp,1),1); %baseline is minimum in search window for all individual events, thus result always positive
                    AllEvMinBsLns=AllEvResp-AllBsLns;

                    Burden.AllEvsMinBsln=nanmean(nansum(AllEvMinBsLns))*TotalAHI/60; %warning nansum assumes dT = 1 (multiply by dT if dT~=1)
                end
            end
        end
    end
    
end


function [features]=QuantifyPeriodicity(EventMean,TimeEnsemble,T0,plotfig)

Fs=1/mean(diff(TimeEnsemble));
EventMean=EventMean-mean(EventMean);
[autocor,lags] = xcorr(EventMean,round((TimeEnsemble(end)-TimeEnsemble(1))*Fs),'coeff');

[pks,lcs] = findpeaks(autocor,lags, ...
    'MinPeakDistance',0.6*T0,'MinPeakheight',0.1);
[trs,trlcs] = findpeaks(-autocor,lags, ...
    'MinPeakDistance',0.6*T0,'MinPeakheight',0.1);
T1 = mean(diff(lcs))/Fs;
T2 = mean(diff(trlcs))/Fs;
T=mean([T1 T2]);
pks(lcs<=0)=[];
lcs(lcs<=0)=[];

trs=-trs;
trs(trlcs<=0)=[];
trlcs(trlcs<=0)=[];

if plotfig
    figure(21);clf(21);
    plot(lags/Fs,autocor)
    xlabel('Lag (secs)')
    ylabel('Autocorrelation')
    hold on
    plpk=plot(lcs/Fs,pks,'or');
    plot(trlcs/Fs,trs,'og');
    hold off
    legend(plpk,[repmat('Period: ',[1 1]) num2str(T,3)])
    axis([0 2*TimeEnsemble(end) -1.1 1.1])
end

if ~isempty(T)
    features.PeriodACorr=T;
else
    features.PeriodACorr=NaN;
end

if length(pks)>=2
    features.Pk1ToPk0=pks(1);
    features.Pk2ToPk0=pks(2);
elseif length(pks)==1
    features.Pk1ToPk0=pks(1);
    features.Pk2ToPk0=0;
else
    features.Pk1ToPk0=0;
    features.Pk2ToPk0=0;
end

if length(trs)>=2
    features.Tr1ToPk0=trs(1);
    features.Tr2ToPk0=trs(2);
elseif length(trs)==1
    features.Tr1ToPk0=trs(1);
    features.Tr2ToPk0=0;
else
    features.Tr1ToPk0=0;
    features.Tr2ToPk0=0;
end

end
   
    
    
    




