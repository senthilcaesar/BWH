function [Systolic_int,Diastolic_int,Map_int] = new_CBFV_feature(CBFVWave,Time,ChannelsList,ChannelsFs)
CBFVBackup=CBFVWave; % for debugging only
filt_CBFV=conv(CBFVWave, ones(7,1)/7, 'same');

 %% create time based on sampling freq of Blood pressure
    %     TimeBP=DataEventHypnog_Mat(:,1);
    Fs_CBFV=128;
    N_CBFV=length(filt_CBFV);
    N_timeCBFV = round((N_CBFV/Fs_CBFV*Fs_CBFV));
    TimeCBFV=(Time(1):(1/Fs_CBFV):Time(1)+(N_timeCBFV-1)*(1/Fs_CBFV))';
    
%% finding peaks and valleys
    [pks_CBFV,CBFV_peak_locs]=findpeaks(filt_CBFV,'MinPeakProminence',0.10,'MinPeakDistance',20,'MinPeakHeight',0.15,'MinPeakWidth',10);
%     figure; plot(TimeCBFV,filt_CBFV); hold on
%     plot(TimeCBFV(CBFV_peak_locs),filt_CBFV(CBFV_peak_locs),'*r');
    
    figure; plot(filt_CBFV); hold on
    plot(CBFV_peak_locs,pks_CBFV,'*r');

    systolic_locs_CBFV=CBFV_peak_locs';
    systolic_val_CBFV=pks_CBFV';
    
    for i=1:length(CBFV_peak_locs)-1
        [diastolic_val(i),I]= min(filt_CBFV((CBFV_peak_locs(i)+1):CBFV_peak_locs(i+1)));
        diastolic_locs(i)=I+CBFV_peak_locs(i);
    end
    
%% Interpolation of peaks and valleys
    %     Systolic_int=interp1(systolic_locs,systolic_val,[1:length(BPWave)],'linear');
    Systolic_int=interp1(systolic_locs_CBFV,systolic_val_CBFV,[1:length(filt_CBFV)],'pchip',nanmean(filt_CBFV));
    hold on;
%     plot(TimeCBFV,Systolic_int,'g');
plot(Systolic_int,'g');
    Diastolic_int=interp1(diastolic_locs,diastolic_val,[1:length(filt_CBFV)],'pchip',nanmean(filt_CBFV));
    
    hold on;
%     plot(TimeCBFV,Diastolic_int,'r');
     plot(Diastolic_int,'r');
%%
window = 1;
filt_set2_CBFV = filt_CBFV;

for i2 = 1:window*128:length(filt_set2_CBFV)-window*128
    
    tmp_signal = [];
    tmp_signal = filt_set2_CBFV(i2:i2+(window*128)-1);
    tmp_rise = [];
    
    if (~isempty(tmp_signal(~isnan(tmp_signal))))
        
    max_val = abs(max(tmp_signal));
    min_val = abs(min(tmp_signal));
    tmp_rise = (max_val - min_val)/((max_val + min_val)/2);
    
    if tmp_rise >1.0
        
        filt_set2_CBFV(i2:i2+(window*128)-1) = NaN;
        
    end
    
    end
end 
% filt_set2_CBFV = [filt_CBFV(isnan(filt_CBFV));filt_set2_CBFV];
figure;plot(filt_CBFV)
hold on;
plot(filt_set2_CBFV)

filt_set3_CBFV = filt_set2_CBFV;
filt_set3_CBFV(filt_set3_CBFV<10) = NaN;
filt_set3_CBFV(filt_set3_CBFV>110) = NaN;
tmp_nan_loc = find(isnan(filt_set3_CBFV));

for i3 = 1:2:length(tmp_nan_loc)-1
    
    tmp_diff = tmp_nan_loc(i3+1)-tmp_nan_loc(i3);
    
    if tmp_diff < 128*5
        
        filt_set3_CBFV(tmp_nan_loc(i3):tmp_nan_loc(i3+1)) = NaN;
        
    end
end


 %% finding peaks and valleys
    [new_pks_CBFV,new_CBFV_peak_locs]=findpeaks(filt_set3_CBFV,'MinPeakProminence',0.70,'MinPeakDistance',30,'MinPeakHeight',27,'MinPeakWidth',30);
    figure;plot(TimeCBFV,filt_CBFV),hold on,plot(TimeCBFV,filt_set3_CBFV);
    plot(TimeCBFV(new_CBFV_peak_locs),filt_set2_CBFV(new_CBFV_peak_locs),'*r');
    
    new_systolic_locs_CBFV=new_CBFV_peak_locs';
    new_systolic_val_CBFV=new_pks_CBFV';
    
    for i=1:length(new_systolic_locs_CBFV)-1
        [new_diastolic_val(i),I]= min(filt_set2_CBFV((new_systolic_locs_CBFV(i)+1):new_systolic_locs_CBFV(i+1)));
        new_diastolic_locs(i)=I+new_systolic_locs_CBFV(i);
    end

    nan_loc =find(isnan(filt_set3_CBFV));
    %% Interpolation of peaks and valleys
    %     Systolic_int=interp1(systolic_locs,systolic_val,[1:length(BPWave)],'linear');
    new_Systolic_int=interp1(new_systolic_locs_CBFV,new_systolic_val_CBFV,[1:length(filt_set3_CBFV)],'pchip',nanmean(filt_set3_CBFV));
    new_Systolic_int(nan_loc(1:end)) = NaN;
    
    hold on;
    plot(TimeCBFV,new_Systolic_int,'g');
    new_Diastolic_int=interp1(new_diastolic_locs,new_diastolic_val,[1:length(filt_set3_CBFV)],'pchip',nanmean(filt_set3_CBFV));
    
    hold on;
    new_Diastolic_int(nan_loc(1:end)) = NaN;
    plot(TimeCBFV,new_Diastolic_int,'r');
    
%     filt_set3_CBFV = filt_set2_CBFV;

    Map_int=[((1/3)*new_Systolic_int)+((2/3)*new_Diastolic_int)];








