disp('Starting ArInt Adjustment');
                
                thresh=2; % threshold of 2 or 3
                
                % 1: get arousal from ArIntCon
                clear I
                I=Sigs.SleepStage<=0 | Sigs.SleepStage>5; % remove arousals in wake
                Sigs.ArIntCont(I)=0;
                
                clear I I1 I2 I1N I2N ArousalSig_
                ArousalSig_=Sigs.ArIntCont>=thresh;
                
                % merge closely-spaced dropouts, <2.9s %
                % 2.9s
                removeshorterthani = 2.9/(1/Fs);
                I = diff([NaN;1-ArousalSig_]);
                I1 = find(I==1);
                I2 = find(I==-1);
                [I1N,I2N]=TidyStartEndEventList(I1,I2,length(I));
                lengthsN = I2N-I1N;
                removeX = lengthsN<removeshorterthani;
                for j=1:length(lengthsN)
                    if ~removeX(j)
                        continue
                    end
                    ArousalSig_(I1N(j):I2N(j))=1;
                end
                
                % remove arousals <3s
                removeshorterthani = 3/(1/Fs);
                removelongerthani = 15/(1/Fs);
                clear I I1 I2 I1N I2N
                I = diff([NaN;ArousalSig_]);
                I1 = find(I==1);
                I2 = find(I==-1);
                [I1N,I2N]=TidyStartEndEventList(I1,I2,length(I));
                lengthsN = I2N-I1N;
                removeX = lengthsN<removeshorterthani;
                for j=1:length(lengthsN)
                    if ~removeX(j)
                        continue
                    end
                    ArousalSig_(I1N(j):I2N(j))=0;
                end
                
                clear I I1 I2 I1N I2N
                I = diff([NaN;ArousalSig_]);
                I1 = find(I==1);
                I2 = find(I==-1);
                [I1N,I2N]=TidyStartEndEventList(I1,I2,length(I));
                lengthsN = I2N-I1N;
                N_ArCon = length(I1); % final number of arousals after all criteria (wake/<3s)
                
                if ~isfield(Sigs,'ArAdj')
                    
                    %2. adjust arousal duration/create new arousal signal
                    Arousal_Adjusted=zeros(length(Sigs.Arousal),1);
                    clear I I1
                    I1 = diff([NaN;Sigs.Arousal]);
                    I = find(I1==1);
                    N_ArManual=length(I);
                    
                    ArAdjInd=[];
                    for j=1:length(I)
                        endi=find(I2N>I(j),1,'first');
                        I1temp= I1N(endi);
                        I2temp=I2N(endi);
                        rangei=I1temp-5*Fs; % checking for manual Ar marking within 5s of Arcont
                        if I(j)>=rangei
                            if 0 % using manual scoring as starting time
                                ArDur_Adjusted=I2temp-I(j);
                                if  ArDur_Adjusted>=3*Fs && ArDur_Adjusted <=15*Fs
                                    Arousal_Adjusted(I(j):I2temp)=1;
                                elseif ArDur_Adjusted >15*Fs
                                    Arousal_Adjusted(I(j):I(j)+15*Fs)=1;
                                elseif ArDur_Adjusted <=3*Fs
                                    Arousal_Adjusted(I1temp:I2temp)=1; % assigning based on cont ar if dur<3s
                                    
                                end
                                ArAdjInd=[ArAdjInd;j];
                                
                            else % SS suggested using autoscoring as the starting time
                                ArDur_Adjusted=I2temp-I1temp;
                                if  ArDur_Adjusted <=15*Fs
                                    Arousal_Adjusted(I1temp:I2temp)=1;
                                elseif ArDur_Adjusted >15*Fs
                                    Arousal_Adjusted(I1temp:I1temp+15*Fs)=1;
                                end
                                ArAdjInd=[ArAdjInd;j];
                            end
                            
                        end
                    end
                    N_ArAdj=length(ArAdjInd);
                    
                    figure(111); clf(111);plot(Sigs.Time,Sigs.ArIntCont); hold on;
                    plot(Sigs.Time,ArousalSig_,'r'); hold on;
                    plot(Sigs.Time,Sigs.Arousal,'k');
                    plot(Sigs.Time,Arousal_Adjusted,'g --');
                    
                    % plot
                    figure(333); clf(333);
                    ax1(1)=subplot(3,1,1);
                    plot(Sigs.Time,Sigs.ArIntCont,'Color',[0.3 0.74 0.93]);hold on;plot(Sigs.Time,ArousalSig_,'Color',[0.85 0.33 0.09]);
                    plot(Sigs.Time,Sigs.Arousal,'k','LineWidth',2);hold on;plot(Sigs.Time,Arousal_Adjusted,'g');
                    legend('ArInContinuous','ArSig-EEG','Manual Scoring','Arousals-Adjusted')
                    ax1(2)=subplot(3,1,2);
                    plot(Sigs.Time,EEG(:,1));hold on;plot(Sigs.Time,EEG(:,2)+200);
                    legend('EEG1','EEG2');
                    ax1(3)=subplot(3,1,3);
                    plot(Sigs.Time,Sigs.SleepStage);legend('SleepStage')
                    xlabel('Time(s)')
                    linkaxes(ax1,'x')
                    
                    
                    SensitivityT=[SensitivityT; table(SubId,thresh,N_ArCon,N_ArManual,N_ArAdj)];
                    
                    % 3. re-assign Sigs.Arousal
                    clear  Sigs.Arousal
                    Sigs.Arousal=Arousal_Adjusted;
                    
                    % 4. re-calculate ArInt
                    if size(EEG,2)>=1
                        ArousalIntensity=NaN(length(Sigs.Time),size(EEG,2));
                        for jj=1:size(EEG,2)
                            try
                                [ArousalIntensitytmp,~]=ArousalIntensityRun(EEG(:,jj),Sigs.Arousal,Sigs.SleepStage,Sigs.Time,TimeEEG);
                                ArousalIntensity(:,jj)=ArousalIntensitytmp;
                            end
                            
                        end
                    else
                        ArousalIntensity=NaN(size(Sig.Time));
                    end
                    Sigs.ArInt=max(ArousalIntensity,[],2);
                    Sigs.ArAdj=1;
                    
                    save([resDir2 filename '-Sigs.mat'],'Sigs','Fs')
                else
                    disp('Already had ArInt Adjustment--using existing adjusted signal');
                    %                     ArousalSig_=Sigs.ArIntCont;
                end
       
        