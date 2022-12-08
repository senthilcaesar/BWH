function [WSpredlogit,acc,ExcludeEEG,NoiseBinary,EEGRef,Tselect] = WSfromoneEEG(Exclude,mdlA,mdlNoise,RefTable,DataEventHypnog_Mat_ds,j,dT,ExcludeAR)

    NoisePredThres=0.5;
    MinNoiseTime=60;

    Tselect=table();
    Tselect.Pbeta = eval(['DataEventHypnog_Mat_ds.Pbetalogfilt' num2str(j)]);
    Tselect.Palpha = eval(['DataEventHypnog_Mat_ds.Palphalogfilt' num2str(j)]);
    Tselect.Ptheta = eval(['DataEventHypnog_Mat_ds.Pthetalogfilt' num2str(j)]);
    Tselect.Pdelta = eval(['DataEventHypnog_Mat_ds.Pdeltalogfilt' num2str(j)]);

    ExcludeEEG = Exclude | Tselect.Pbeta==-Inf;

    [WSpredlogit,Tselect,EEGRef] = PredWakeSleep(Tselect,mdlA,RefTable,ExcludeEEG);

    Tselect.Total = Tselect.Pbeta_ref + Tselect.Palpha_ref + Tselect.Ptheta_ref + Tselect.Pdelta_ref;

    NoiseBinary = RemoveShortSegments(1*(predict(mdlNoise,Tselect)>NoisePredThres),MinNoiseTime,dT);

    %noise based on j=last EEG signal only
    Tselect.Total = Tselect.Pbeta_ref + Tselect.Palpha_ref + Tselect.Ptheta_ref + Tselect.Pdelta_ref;
    Tselect.WSpredlogit = WSpredlogit;

    NoiseBinary = RemoveShortSegments(1*(predict(mdlNoise,Tselect)>NoisePredThres),MinNoiseTime,dT);

    ExcludeEEG = NoiseBinary==1 | Exclude==1 | Tselect.Pbeta==-Inf;

    [WSpredlogit,Tselect,EEGRef] = PredWakeSleep(Tselect,mdlA,RefTable,ExcludeEEG);
    
    WakeNoAR = 1*(DataEventHypnog_Mat_ds.Epochs==4);
    
    try
        performance = PredictiveValue(1*(WakeNoAR(~ExcludeAR)>0.5),1*(WSpredlogit(~ExcludeAR)>0),WakeNoAR(~ExcludeAR));
        %     acc=performance.Acc_sem_chance_p(1);
        acc=performance.Acc(1);
    catch me
        acc=NaN;
    end
    
    if 0 %debug
    figure(2)
    logitinverse = @(p) 1./(1+exp(-p));
    plot(DataEventHypnog_Mat_ds.Time,[DataEventHypnog_Mat_ds.Epochs DataEventHypnog_Mat_ds.EventsAr logitinverse(WSpredlogit) ExcludeAR-0.1 WakeNoAR-0.2]);
    end



