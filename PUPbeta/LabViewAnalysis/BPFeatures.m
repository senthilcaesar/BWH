function [Systolic_int,Diastolic_int,Map_int]=BPFeatures(BPWave,Time,ChannelsList,Channels_Fs)

I=find(strcmp(ChannelsList,'BPWave'));
if ~isempty(I)
    BPBackup=BPWave; % for debugging only
    %     BPWave=BPBackup;
    %     BPWave=DataEventHypnog_Mat(:,I);
    
    %% artifact rejection
    BPWave(BPWave<=0.45)=NaN; % deleting anything less than 40mmhg
    BPWave(BPWave>1.5*prctile(BPWave,95))=NaN;
    
    %% create time based on sampling freq of Blood pressure
    %     TimeBP=DataEventHypnog_Mat(:,1);
    Fs_BP=Channels_Fs(find(strcmp(ChannelsList,'BPWave')));
    N_BP=length(BPWave);
    N_timeBP = round((N_BP/Fs_BP*Fs_BP));
    TimeBP=(Time(1):(1/Fs_BP):Time(1)+(N_timeBP-1)*(1/Fs_BP))';
    
    %% finding peaks and valleys
    [pks_BP,BP_peak_locs]=findpeaks(BPWave,'MinPeakProminence',0.3,'MinPeakDistance',45,'MinPeakHeight',0.45);
    figure; plot(TimeBP,BPWave); hold on
    plot(TimeBP(BP_peak_locs),BPWave(BP_peak_locs),'*r');
    
    systolic_locs=BP_peak_locs';
    systolic_val=pks_BP';
    
    for i=1:length(BP_peak_locs)-1
        [diastolic_val(i),I]= min(BPWave((BP_peak_locs(i)+1):BP_peak_locs(i+1)));
        diastolic_locs(i)=I+BP_peak_locs(i);
    end
    %% Interpolation of peaks and valleys
    %     Systolic_int=interp1(systolic_locs,systolic_val,[1:length(BPWave)],'linear');
    Systolic_int=interp1(systolic_locs,systolic_val,[1:length(BPWave)],'pchip',nanmean(BPWave));
    hold on;
    plot(TimeBP,Systolic_int,'g');
    Diastolic_int=interp1(diastolic_locs,diastolic_val,[1:length(BPWave)],'pchip',nanmean(BPWave));
    
    hold on;
    plot(TimeBP,Diastolic_int,'r');
    %% removing calibration and gaps with no data
    for i11  = 1:length(systolic_locs)
        time_val(1,i11) = TimeBP(systolic_locs(1,i11));
        time_locs(1,i11) = find(time_val(1,i11)==TimeBP);
    end
    t_1 = zeros(1,length(systolic_locs)); % first location of systolic peak
    t_2 = zeros(1,length(systolic_locs)); % second location of systolic peak
    new_width = zeros(1,length(TimeBP));
    combined_clipped_BP_full = BPWave;
    
    for j1 = 1:length(time_val)-1
        
        dist = time_val(1,j1+1) - time_val(1,j1);
        
        if abs(dist)>5
            t_1(1,j1) = systolic_locs(1,j1);
            t_2(1,j1) = systolic_locs(1,j1+1);
            combined_clipped_BP_full(t_1(1,j1)+25:t_2(1,j1)-25) = NaN; % to preserve the peaks added 25 points
            %         t1_loc = time_locs(1,j1);
            %         t2_loc = time_locs(1,j1+1);
            %         a = length(t1_loc+25:t2_loc-25);
            %         new_width(1,(t1_loc+25:t2_loc-25)) = 3*ones(1,a);
        end
    end
    
    
%     figure;
%     plot(TimeBP,BPWave),hold on,plot(TimeBP(systolic_locs),systolic_val,'*r')
%     plot(TimeBP,combined_clipped_BP_full,'g')
    
    IndRem=[nonzeros(t_1(:)), nonzeros(t_2(:))];
    
    % case 1: several consecutive gaps represent overall noisy
    % section-->delete the entire section
    IndRem(:,3)=[1;(IndRem(2:end,1)-IndRem(1:end-1,2))];
    IndRem(:,4)=(IndRem(:,3)==0);
    I = diff([IndRem(:,4);NaN]);
    I1 = find(I==1);
    I2 = find(I==-1);
    [I1N,I2N]=TidyStartEndEventList(I1,I2,length(I));
    I1I2N=[I1N,I2N];
    IndRemF=IndRem(:,1:2);
    for kk=1:size(I1I2N,1)
        tempind1=I1I2N(kk,1);
        tempind2=I1I2N(kk,2);
        IndRemF(tempind1,2)=IndRem(tempind2,2);
        IndRemF(tempind1+1:tempind2,:)=[NaN];
    end
    
    IndRemF(isnan(IndRemF(:,1)),:)=[];
    
    % case 2: nans in the beginning
    if systolic_locs(1)>3*Fs_BP % presence of nans
        IndRemF=[[1,systolic_locs(1)-1];IndRemF];
    end
    
    % case 3: nans in the end
    if (TimeBP(end)-TimeBP(systolic_locs(end)))>=3
        IndRemF=[IndRemF;[systolic_locs(end)+1,N_timeBP]];
    end
    
    %% remove the gaps from systolic and diastolic trace
    for kk=1:size(IndRemF,1)
        tempind1=IndRemF(kk,1);
        tempind2=IndRemF(kk,2);
        Systolic_int(tempind1:tempind2)=NaN;
        Diastolic_int(tempind1:tempind2)=NaN;
    end
    
    %% Calculate MAP
    Map_int=[((1/3)*Systolic_int)+((2/3)*Diastolic_int)];
  
    
    %% convert from voltage to mmHg
    if nanmean(Map_int)<2 % means range of signal is in voltage 0-2 etc
        Map_int=Map_int.*100;
        Systolic_int=Systolic_int*100;
        Diastolic_int=Diastolic_int*100;
    end
    figure;
    plot(TimeBP,BPWave*100); hold on
    plot(TimeBP,Systolic_int,'r');
    plot(TimeBP,Diastolic_int,'g');
    plot(TimeBP,Map_int,'k');
end