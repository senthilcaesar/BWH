function Score=ScoreArousal(ArStruct,fs)
%[Score]=ScoreArousal(ArStruct,fs)
%ArStruct: Arousal struct with the following fields 
%     ArStruct.sig: EEG signal during an arousal. 
%         It has to be prefiltered (lowpass filter with cut-off frequency of 50Hz, notch filter with cut-off frequency of 60 Hz)
%     ArStruct.presig: EEG signal before the arousal. The length of this
%         signal has to be the same as that of ArStruct.sig
%     ArStruct.training_set: Training data set with 542 rows (observations)
%         and 16 columns. First 15 columns are the wavelet features and las column (16) contains the visual scoring 
%fs: EEG sampling frequency
%Score: Arousal score, a number between 0 and 9

ar_sig=ArStruct.sig;
ar_pre_sig=ArStruct.presig;
TrainMtx=ArStruct.training_set;

decomLevel=5; %%% Wavelet decomposition level
wavefun='db4'; %%% Wavelet function

%%%This function calculates the 15 features for the input signal.
%%%"baseline_ftrs" and "ar_ftrs" are two 1x15 vectors  
[baseline_ftrs]=FindWaveFtrs(ar_pre_sig,decomLevel,wavefun,fs);
[ar_ftrs]=FindWaveFtrs(ar_sig,decomLevel,wavefun,fs);

%%%This function classifies the input arousal segment based on its feature
%%%vector both during and before arousal and generates a score between 0
%%%and 9 depending on the intensity of the arousal
Score=ClassifyArousal(ar_ftrs,baseline_ftrs,TrainMtx);

end
   

