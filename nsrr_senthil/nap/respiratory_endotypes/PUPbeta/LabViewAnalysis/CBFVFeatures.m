function [Systolic_int,Diastolic_int,Map_int]=CBFVFeatures(CBFVWave,Time,ChannelsList,Channels_Fs)

I=find(strcmp(ChannelsList,'CBFVWave'));
if ~isempty(I)
    CBFVBackup=CBFVWave; % for debugging only
    %     BPWave=BPBackup;
    %     BPWave=DataEventHypnog_Mat(:,I);
    filt_CBFV=conv(CBFVWave, ones(7,1)/7, 'same');
%     filt_CBFV = highpass(filt_CBFV_tmp,75,128);
%     figure;plot(filt_CBFV)

 
    %% create time based on sampling freq of Blood pressure
    %     TimeBP=DataEventHypnog_Mat(:,1);
    Fs_CBFV=Channels_Fs(find(strcmp(ChannelsList,'CBFVWave')));
    N_CBFV=length(CBFVWave);
    N_timeCBFV = round((N_CBFV/Fs_CBFV*Fs_CBFV));
    TimeCBFV=(Time(1):(1/Fs_CBFV):Time(1)+(N_timeCBFV-1)*(1/Fs_CBFV))';
   
%      [CBFVwavtmp,lCBFVwavtmp,I]=Spikedetect(TimeCBFV,filt_CBFV,Fs_CBFV);

    %% finding peaks and valleys
    [pks_CBFV,CBFV_peak_locs]=findpeaks(CBFVWave,'MinPeakProminence',0.10,'MinPeakDistance',20,'MinPeakHeight',0.15,'MinPeakWidth',10);
    figure; plot(TimeCBFV,CBFVWave); hold on
    plot(TimeCBFV(CBFV_peak_locs),CBFVWave(CBFV_peak_locs),'*r');
    
    systolic_locs_CBFV=CBFV_peak_locs';
    systolic_val_CBFV=pks_CBFV';
    
    for i=1:length(CBFV_peak_locs)-1
        [diastolic_val(i),I]= min(filt_CBFV((CBFV_peak_locs(i)+1):CBFV_peak_locs(i+1)));
        diastolic_locs(i)=I+CBFV_peak_locs(i);
    end
    %% Interpolation of peaks and valleys
    %     Systolic_int=interp1(systolic_locs,systolic_val,[1:length(BPWave)],'linear');
    Systolic_int=interp1(systolic_locs_CBFV,systolic_val_CBFV,[1:length(CBFVWave)],'pchip',nanmean(CBFVWave));
    hold on;
    plot(TimeCBFV,Systolic_int,'g');
    Diastolic_int=interp1(diastolic_locs,diastolic_val,[1:length(CBFVWave)],'pchip',nanmean(CBFVWave));
    
    hold on;
    plot(TimeCBFV,Diastolic_int,'r');
    %%
    %%
%test1
% check if the values are not good
tmp_low_location  = find(filt_CBFV<10==1);
tmp1 = 1;
tmp_start = [];
tmp_end = [];
t2_tmp = 1;

for t2 = 1:length(tmp_low_location)-1
    
    tmp_diff = tmp_low_location(t2+1)-tmp_low_location(t2);
    
    if tmp_diff > 2
        
    tmp_start(1,t2_tmp) = tmp_low_location(t2+1);
    tmp_end(1,t2_tmp) = tmp_low_location(t2+1)-1;
    t2_tmp = t2_tmp+1;
    
    end
    
end


test2_tmp_filt_CBFV = filt_CBFV;
for t1 = 1:length(tmp_start)
    
    if tmp_start(t1)>100 && (tmp_start(t1)+100)<length(test2_tmp_filt_CBFV)
    
    tmp = [];
    tmp = find((test2_tmp_filt_CBFV(tmp_start(t1)-99:tmp_start(t1)+100)<0.10)==1);
    
%     if length(tmp)>100 && ~isempty(nonzeros(tmp))
    if ~isempty(nonzeros(tmp))
%     tmp_ip(tmp1,1:length(tmp)) = tmp;
%     tmp1 = tmp1+1;
    
    test2_tmp_filt_CBFV((tmp_start(t1)-29:tmp_start(t1)+29)) = NaN;
    
    end
    
    end
    
end

tmp_nan = find(isnan(test2_tmp_filt_CBFV)==1);
tmp_start_nan = [];
tmp_end_nan = [];
nan_tmp = 1;

for t3 = 1:length(tmp_nan)-1
    
    tmp_start =  tmp_nan(t3+1) - tmp_nan(t3);
    
    if tmp_start >1
    tmp_start_nan(1,nan_tmp) = tmp_nan(t3+1);
    tmp_end_nan(1,nan_tmp) = tmp_nan(t3);
%     if tmp_start_nan(1,nan_tmp) - tmp_end_nan(1,nan_tmp)<3200
%     tmp_filt_CBFV(tmp_end_nan(1,nan_tmp):tmp_start_nan(1,nan_tmp)) = NaN;
%     end
    nan_tmp = nan_tmp+1;
    end    
    
end

tmp_start_nan =  [tmp_nan(1) tmp_start_nan];


for o2 = 1:length(tmp_start_nan)-1
    
%     if (tmp_end_nan(o2) - tmp_start_nan(o2))<3200
    test2_tmp_filt_CBFV(tmp_start_nan(o2):tmp_end_nan(o2)) = NaN;
%     end
    
end

[test1_pks_CBFV,test1_CBFV_peak_locs]=findpeaks(test2_tmp_filt_CBFV,'MinPeakProminence',0.10,'MinPeakDistance',20,'MinPeakHeight',0.15,'MinPeakWidth',10);
test1_systolic_locs_CBFV=test1_CBFV_peak_locs';
test1_systolic_val_CBFV=test1_pks_CBFV';
%--------------------Calculating CBF troughs(filtered)--------------------
test1_new_Systolic_int_CBFV=interp1(test1_systolic_locs_CBFV,test1_systolic_val_CBFV,[1:length(test2_tmp_filt_CBFV)],'pchip',nanmean(test2_tmp_filt_CBFV));

%--------------------Calculating CBF troughs(filtered)--------------------
for i=1:length(test1_systolic_locs_CBFV)-1
    [test1_filt_CBF_trough(i),I]=min(test2_tmp_filt_CBFV((test1_systolic_locs_CBFV(i)+1):test1_systolic_locs_CBFV(i+1)));
    test1_filt_CBF_trough_locs(i)=I+test1_systolic_locs_CBFV(i);
end 

test1_filt_CBF_trough_int=interp1(test1_filt_CBF_trough_locs,test1_filt_CBF_trough,[1:length(test2_tmp_filt_CBFV)],'pchip',nanmean(test2_tmp_filt_CBFV));


%%
% test2
% test with the distance
tmp_systolic_locs_CBFV = test1_systolic_locs_CBFV;
tmp_systolic_val_CBFV = test1_systolic_val_CBFV;

for i =  1:length(test1_systolic_locs_CBFV)-1
    
    tmp = []; 
    tmp = test1_systolic_locs_CBFV(i+1) - test1_systolic_locs_CBFV(i);
    
    if tmp> 100
        
    tmp_systolic_locs_CBFV(1,i) = NaN;
    tmp_systolic_val_CBFV(1,i) = NaN;
    
    end
end

tmp_nan_systolic = find(isnan(tmp_systolic_locs_CBFV)==1);
tmp_filt_CBFV =test2_tmp_filt_CBFV ;

tmp_start = [];
tmp_end =[];

for i1 = 2:length(tmp_nan_systolic)
    
    tmp_start = tmp_nan_systolic(i1)-1;
    tmp_end = tmp_nan_systolic(i1);
    tmp_filt_CBFV(tmp_start:tmp_end) = NaN;
    
    
end

[test1_pks_CBFV,test1_CBFV_peak_locs]=findpeaks(tmp_filt_CBFV,'MinPeakProminence',0.10,'MinPeakDistance',20,'MinPeakHeight',0.15,'MinPeakWidth',10);
test2_systolic_locs_CBFV=test1_CBFV_peak_locs';
test2_systolic_val_CBFV=test1_pks_CBFV';
%--------------------Calculating CBF peaks(filtered)--------------------
test2_new_Systolic_int_CBFV=interp1(test2_systolic_locs_CBFV,test2_systolic_val_CBFV,[1:length(tmp_filt_CBFV)],'pchip',nanmean(tmp_filt_CBFV));

%--------------------Calculating CBF troughs(filtered)--------------------
for i=1:length(test2_systolic_locs_CBFV)-1
    [test2_filt_CBF_trough(i),I]=min(tmp_filt_CBFV((test2_systolic_locs_CBFV(i)+1):test2_systolic_locs_CBFV(i+1)));
    test2_filt_CBF_trough_locs(i)=I+test2_systolic_locs_CBFV(i);
end 

test2_filt_CBF_trough_int=interp1(test2_filt_CBF_trough_locs,test2_filt_CBF_trough,[1:length(tmp_filt_CBFV)],'pchip',nanmean(tmp_filt_CBFV));
 
%%
% rise in amplitude
test3_systolic_val_CBFV = test2_systolic_val_CBFV;
test3_systolic_locs_CBFV = test2_systolic_locs_CBFV;

for o = 1:length(test3_systolic_val_CBFV)-1
       
    tmp = abs((test3_systolic_val_CBFV(o+1) - test3_systolic_val_CBFV(o))/test3_systolic_val_CBFV(o));
    
    if tmp > 0.20
        test3_systolic_val_CBFV(o+1) = NaN;
        test3_systolic_locs_CBFV(o+1) = NaN;
    end
        
end

tmp_nan = [];
tmp_nan = find(isnan(test3_systolic_locs_CBFV)==1);
tmp_start_nan = [];
tmp_end_nan = [];
nan_tmp = 1;

tmp_start = [];

for t3 = 1:length(tmp_nan)-1
    
    tmp_start =  tmp_nan(t3+1) - tmp_nan(t3);
    
    if tmp_start >1
    tmp_start_nan(1,nan_tmp) = tmp_nan(t3+1);
    tmp_end_nan(1,nan_tmp) = tmp_nan(t3);
%     if tmp_start_nan(1,nan_tmp) - tmp_end_nan(1,nan_tmp)<3200
%     tmp_filt_CBFV(tmp_end_nan(1,nan_tmp):tmp_start_nan(1,nan_tmp)) = NaN;
%     end
    nan_tmp = nan_tmp+1;
    end    
    
end

tmp_start_nan =  [tmp_nan(1) tmp_start_nan];

test3_systolic_locs_CBFV = test3_systolic_locs_CBFV(~isnan(test3_systolic_locs_CBFV));
test3_systolic_val_CBFV = test3_systolic_val_CBFV(~isnan(test3_systolic_val_CBFV));

new_Systolic_int_CBFV=interp1(test3_systolic_locs_CBFV,test3_systolic_val_CBFV(~isnan(test3_systolic_val_CBFV)),[1:length(test2_tmp_filt_CBFV)],'pchip',nanmean(tmp_filt_CBFV));
tmp_new_Systolic_int_CBFV = new_Systolic_int_CBFV;

% [filt_CBF_peak,filt_CBF_peak_locs]=findpeaks(test2_tmp_filt_CBFV,'MinPeakProminence',0.1,'MinPeakDistance',50,'MinPeakHeight',0.15,'MinPeakWidth',10);
% filt_CBF_peak_int= interp1(filt_CBF_peak_locs,filt_CBF_peak,[1:length(test2_tmp_filt_CBFV)],'spline');

for i=1:length(test3_systolic_locs_CBFV)-1
    [test3_filt_CBF_trough(i),I]=min(tmp_filt_CBFV((test3_systolic_locs_CBFV(i)+1):test3_systolic_locs_CBFV(i+1)));
    test3_filt_CBF_trough_locs(i)=I+test3_systolic_locs_CBFV(i);
end 

filt_CBF_trough_int=interp1(test3_filt_CBF_trough_locs,test3_filt_CBF_trough,[1:length(test2_tmp_filt_CBFV)],'pchip',nanmean(tmp_filt_CBFV));

for o1= 1:length(tmp_start_nan)-1
    
    tmp_new_Systolic_int_CBFV(tmp_start_nan(o1):tmp_end_nan(o1)) = NaN;
    filt_CBF_trough_int(tmp_start_nan(o1):tmp_end_nan(o1)) = NaN;
    
end

    %% removing calibration and gaps with no data
    for i11  = 1:length(test3_systolic_locs_CBFV)
        time_val(1,i11) = TimeCBFV(test3_systolic_locs_CBFV(1,i11));
        time_locs(1,i11) = find(time_val(1,i11)==TimeCBFV);
    end
    t_1 = zeros(1,length(test3_systolic_locs_CBFV)); % first location of systolic peak
    t_2 = zeros(1,length(test3_systolic_locs_CBFV)); % second location of systolic peak
    new_width = zeros(1,length(TimeCBFV));
    combined_clipped_CBFV_full = tmp_filt_CBFV;
    
    for j1 = 1:length(time_val)-1
        
        dist = time_val(1,j1+1) - time_val(1,j1);
        
        if abs(dist)>5
            t_1(1,j1) = test3_systolic_locs_CBFV(1,j1);
            t_2(1,j1) = test3_systolic_locs_CBFV(1,j1+1);
            combined_clipped_CBFV_full(t_1(1,j1)+25:t_2(1,j1)-25) = NaN; % to preserve the peaks added 25 points
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
    if test3_systolic_locs_CBFV(1)>3*Fs_CBFV % presence of nans
        IndRemF=[[1,test3_systolic_locs_CBFV(1)-1];IndRemF];
    end
    
    % case 3: nans in the end
    if (TimeCBFV(end)-TimeCBFV(test3_systolic_locs_CBFV(end)))>=3
        IndRemF=[IndRemF;[test3_systolic_locs_CBFV(end)+1,N_timeCBFV]];
    end
    
    %% remove the gaps from systolic and diastolic trace
    for kk=1:size(IndRemF,1)
        tempind1=IndRemF(kk,1);
        tempind2=IndRemF(kk,2);
        tmp_new_Systolic_int_CBFV(tempind1:tempind2)=NaN;
        filt_CBF_trough_int(tempind1:tempind2)=NaN;
    end
    
    %% Calculate MAP
    Map_int=[((1/3)*tmp_new_Systolic_int_CBFV)+((2/3)*filt_CBF_trough_int)];
  
    
    %% convert from voltage to mmHg
    if nanmean(Map_int)<2 % means range of signal is in voltage 0-2 etc
        Map_int=Map_int.*100;
        tmp_new_Systolic_int_CBFV=tmp_new_Systolic_int_CBFV*100;
        filt_CBF_trough_int=filt_CBF_trough_int*100;
    end
    
%     combined_clipped_CBFV_full_temp1 = [];
%     combined_clipped_CBFV_full_temp2 = [];
%     combined_clipped_CBFV_full_temp1 = combined_clipped_CBFV_full;
%     combined_clipped_CBFV_full_temp2(isnan(combined_clipped_CBFV_full_temp1)) = 0;
%     combined_clipped_CBFV_full_temp = lowpass(combined_clipped_CBFV_full_temp2,10,128);
%     figure;
%     plot(combined_clipped_CBFV_full_temp)
    
    figure;
    plot(TimeCBFV,combined_clipped_CBFV_full); hold on
    plot(TimeCBFV,tmp_new_Systolic_int_CBFV,'r');
    plot(TimeCBFV,filt_CBF_trough_int,'g');
    plot(TimeCBFV,Map_int,'k');
end