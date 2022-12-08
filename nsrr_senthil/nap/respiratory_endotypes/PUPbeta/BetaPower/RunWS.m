function [WPr,WArPr,WPrB,WArPrB,WSinfo]=RunWS(SigT,ChannelsList,noscoredarinwake,mdlA,mdlNoise,RefTable,mdlAcc,mdlAR)
%Arousal probability is now directly taken from WArPr and WArPrB 
%(later called ArPr and ArPrB)

global settings
    %% Arousal Scoring

logit = @(p) log(p./(1-p));
logitinverse = @(p) 1./(1+exp(-p));

    dt = SigT.Time(2)-SigT.Time(1);
    Fs = 1./dt;
    dT = dt;

    dsfactor = dT/dt;
    DataEventHypnog_Mat_ds = SigT;

    dT = DataEventHypnog_Mat_ds.Time(2)-DataEventHypnog_Mat_ds.Time(1);

    DataEventHypnog_Mat_ds.Epochs(DataEventHypnog_Mat_ds.Epochs<0|DataEventHypnog_Mat_ds.Epochs>4)=NaN;

    if ~ismember(['Pbetalogfilt1'],DataEventHypnog_Mat_ds.Properties.VariableNames)
        disp('missing Pbetalogfilt1');
    end

    SpO2off = 1*(DataEventHypnog_Mat_ds.SpO2<50);
    
    if isfield(settings,'ChangeSpO2offinRunWS') 
        % no setting or 0 is no change (default), 
        % 1 is ForceSpO2 to all bad, 
        % 2 is ForceSpO2 to all good,
        % 3 is check and decide based on badness frequency
        switch settings.ChangeSpO2offinRunWS
            case 0
                 % do nothing, change nothing, process as per usual
            case 1
                % 1 is ForceSpO2 to all bad,
                disp('Warning: Forcing SpO2 signal for WakeSleep analysis to ALL BAD');
                SpO2off = ones(length(SpO2off),1);
            case 2
                % 2 is ForceSpO2 to all good,
                disp('Warning: Forcing SpO2 signal for WakeSleep analysis to ALL GOOD');
                SpO2off = zeros(length(SpO2off),1);
            case 3
                % 3 is check and decide based on signal 'badness'.  
                % this is unfinished, untested, and potentially dangerous.                
                % some data has VERY frequent SpO2 dropouts, i.e. on average, 
                % every eight seconds, it drops out for one sample or second. 
                % This was super annoying for the old RemoveShortSegments code...
                % so here we are testing to see if we can do this on the fly
                % and force the SpO2off signal to Good or Bad accordingly. 
                I = diff([NaN;SpO2off]); I1 = find(I==1); I2 = find(I==-1);
                [I1Ones,I2Ones]=TidyStartEndEventList(I1,I2,length(I));
                lengthsOnes = I2Ones-I1Ones; %length of segments of ones (i.e. bad data)
                %median(lengthsOnes) mean(lengthsOnes)
                if (nnz(lengthsOnes<Fs)/nnz(lengthsOnes<Inf)) > 0.95
                    % if the vast majority (i.e. >95%) is very short bad segments
                    % just force these really short periods to Good, and
                    % then process as per usual.
                    disp('Warning: Forcing lots of v.short BAD segments in signal for WakeSleep analysis to GOOD');
                    BadDataToBeExcluded = find(lengthsOnes<Fs);
                    for i = 1:length(BadDataToBeExcluded) 
                        SpO2off(I1Ones(BadDataToBeExcluded(i)):I2Ones(BadDataToBeExcluded(i)))=0; 
                    end
                else
                    % given the conditions upon which we get here, it may be safe to do nothing, and process normally     
                end                                
                %[I1Zeros,I2Zeros]=TidyStartEndEventList(I2,I1,length(I)); %Zeros are zeros
                %lengthsZeros = I2Zeros-I1Zeros; % length of segments of zeros (i.e. good data)
                %median(lengthsZeros) mean(lengthsZeros)
        end                
    end

    SpO2off = RemoveShortSegments(SpO2off,60,dT);
    SpO2starti = find(DataEventHypnog_Mat_ds.SpO2>50,1,'first');
    SpO2endi = find(DataEventHypnog_Mat_ds.SpO2>50,1,'last');
    SpO2ever = 1*SpO2off; SpO2ever(SpO2starti:SpO2endi)=0; SpO2ever=logical(SpO2ever);
    if sum(SpO2ever==1)>0.90*length(SpO2ever)
        disp('missing SpO2 signal, skipping'); % DLM asks: skipping what exactly?
    end
 
    
    %% Start
    UnknownStage = isnan(DataEventHypnog_Mat_ds.Epochs) | DataEventHypnog_Mat_ds.Epochs<0 | DataEventHypnog_Mat_ds.Epochs>4; 
    if sum(UnknownStage==1)<0.9*length(UnknownStage)
        Exclude = SpO2off==1 | UnknownStage;
        %for judging methods
        ExcludeREM = SpO2off | DataEventHypnog_Mat_ds.EventsAr==1 | DataEventHypnog_Mat_ds.Epochs<0 | DataEventHypnog_Mat_ds.Epochs>3 | isnan(DataEventHypnog_Mat_ds.Epochs);
        REMvsNREM = (DataEventHypnog_Mat_ds.Epochs==3)*1;
        REMvsNREM(ExcludeREM)=NaN;
        
        if ~noscoredarinwake
            %scored arousals in wake: excluding arousals in sleep and excluding sleep onset in wake to judge WS accuracy etc
            ExcludeAR = Exclude | DataEventHypnog_Mat_ds.Epochs==4 & DataEventHypnog_Mat_ds.EventsAr==0 | DataEventHypnog_Mat_ds.Epochs~=4 & DataEventHypnog_Mat_ds.EventsAr==1;
        else
            %no scored arousals in wake: excluding arousals in sleep to judge WS accuracy etc
            ExcludeAR = Exclude | DataEventHypnog_Mat_ds.Epochs~=4 & DataEventHypnog_Mat_ds.EventsAr==1;
        end
    
    else
        Exclude = SpO2off==1;
        ExcludeAR = Exclude;
        disp('Warning, predominantly unscored; ignoring staging entirely');
    end
   

    dsfactor = round(3/dT);
    I = (1:dsfactor:length(Exclude))';
    clear WSpredlogit Acc Excludei
    for j=1:12 %up to 8 EEGs - may fail with one EEG?
        try
            [WSpredlogiti(:,j),Acc(1,j),Excludei(:,j)] = WSfromoneEEG(Exclude(I),mdlA,mdlNoise,RefTable,DataEventHypnog_Mat_ds(I,:),j,3,ExcludeAR(I));
        catch me

        end
    end
    TAccFeatures = AccFeatures(WSpredlogiti,Excludei); %- may fail with one EEG?
    AccPred = predict(mdlAcc,TAccFeatures);
    AccPred(AccPred<0)=0; AccPred(AccPred>1)=1;

    [~,accmaxi] = max(Acc);
    [~,accmaxiPred] = max(AccPred);

    %Slow:
    [WSpredlogit,acc,~,~,EEGRef,Tselect] = WSfromoneEEG(Exclude,mdlA,mdlNoise,RefTable,DataEventHypnog_Mat_ds,accmaxi,dT,ExcludeAR);
    WPr = logitinverse(WSpredlogit);
    
    if accmaxi==accmaxiPred
        WSpredlogitB=WSpredlogit;
        accB = acc;
        EEGRefB=EEGRef;
        TselectB=Tselect;
    else
        [WSpredlogitB,accB,~,~,EEGRefB,TselectB] = WSfromoneEEG(Exclude,mdlA,mdlNoise,RefTable,DataEventHypnog_Mat_ds,accmaxiPred,dT,ExcludeAR);
    end
    WPrB = logitinverse(WSpredlogitB);

    %% Arousal Scoring
 
    if settings.WSArVersion==1
        
    [ArPr,WArPr]=RunWASar(WSpredlogit,mdlAR,dT);
    
    [ArPrB,WArPrB]=RunWASar(WSpredlogitB,mdlAR,dT);
    
    elseif settings.WSArVersion==2
       
        T = Tselect(:,{'Pbeta_ref','Palpha_ref','Ptheta_ref','Pdelta_ref'});
        
        T{:,{'Ptheta_ref'}} = T{:,{'Ptheta_ref'}}*0 + 1.2263;
        T{:,{'Pdelta_ref'}} = T{:,{'Pdelta_ref'}}*0 + 1.4154;
        
        DataEventHypnog_Mat_ds.WSBalance = Tselect.WSpredlogit; %logit(predict(mdlA,DataEventHypnog_Mat_ds));
        DataEventHypnog_Mat_ds.WakeIntensity = logit(predict(mdlA,T));
        DataEventHypnog_Mat_ds.SleepIntensity =  DataEventHypnog_Mat_ds.WakeIntensity - DataEventHypnog_Mat_ds.WSBalance;

    [ArPr,WArPr]=RunWASar2(DataEventHypnog_Mat_ds(:,{'WakeIntensity','SleepIntensity','WSBalance'}),mdlAR,dT);
    %ArSig = max([WPr ArPr]')' > 0.5;
    
     T = TselectB(:,{'Pbeta_ref','Palpha_ref','Ptheta_ref','Pdelta_ref'});
        
        T{:,{'Ptheta_ref'}} = T{:,{'Ptheta_ref'}}*0 + 1.2263;
        T{:,{'Pdelta_ref'}} = T{:,{'Pdelta_ref'}}*0 + 1.4154;
        
        %note overwrite here:
        DataEventHypnog_Mat_ds.WSBalance = Tselect.WSpredlogit; %logit(predict(mdlA,DataEventHypnog_Mat_ds));
        DataEventHypnog_Mat_ds.WakeIntensity = logit(predict(mdlA,T));
        DataEventHypnog_Mat_ds.SleepIntensity =  DataEventHypnog_Mat_ds.WakeIntensity - DataEventHypnog_Mat_ds.WSBalance;
        
        [ArPrB,WArPrB]=RunWASar2(DataEventHypnog_Mat_ds(:,{'WakeIntensity','SleepIntensity','WSBalance'}),mdlAR,dT);
    end
    
   WSinfo.Acc = Acc(:)';
   WSinfo.AccPred = AccPred(:)';
   WSInfo.EEGRef = EEGRef;
   WSInfo.EEGRefB = EEGRefB;
   WSInfo.EEGChannel = accmaxi;
   WSInfo.EEGChannelB = accmaxiPred;
   
