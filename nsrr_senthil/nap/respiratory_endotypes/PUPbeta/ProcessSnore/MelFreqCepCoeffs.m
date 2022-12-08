function MelFreqCepCoeffs(signal,Fs)
bar(1:5,coeffs(1:5,sampStart1),'r')

[cepstra,aspectrum,pspectrum] = melfcc(d, sr, 'wintime',0.2,'hoptime',0.05);
  % Append deltas and double-deltas onto the cepstral vectors
 del = deltas(cep2);
 % Double deltas are deltas applied twice with a shorter window
 ddel = deltas(deltas(cep2,5),5);
 % Composite, 39-element feature vector, just like we use for speech recognition
 cepDpDD = [cep2;del;ddel];
 
 %% Recreate mfcc output from matlab - works but absolute values are off
 [audioIn,fs] = audioread("Rainbow-16-8-mono-114secs.wav");
 d=audioIn;
 sr=fs;
 [cepstra,aspectrum,pspectrum] = melfcc(d, sr, 'wintime',0.1280,'hoptime',0.0640);
 
nbins = 60;
for i = 1:size(cepstra,2)
    figure
    histogram(cepstra(i,:),nbins,"Normalization","pdf")
    title(sprintf("Coefficient %d",i-1))
end
 
%% 
d=FilteredSnd;
sr=Fs;
[d,sr] = audioread('sm1_cln.wav');
 % Look at its regular spectrogram
 ax1=subplot(511);
 specgram(d, 1025, sr);
 
 % Calculate basic RASTA-PLP cepstra and spectra
 [cep1, spec1] = rastaplp(d, sr);
 % .. and plot them
 ax2=subplot(512);
 imagesc(10*log10(spec1)); % Power spectrum, so dB is 10log10
 axis xy
 ax3=subplot(513);
 imagesc(cep1)
 axis xy
 % Notice the auditory warping of the frequency axis to give more 
 % space to low frequencies and the way that RASTA filtering 
 % emphasizes the onsets of static sounds like vowels

 % Calculate 12th order PLP features without RASTA
 [cep2, spec2] = rastaplp(d, sr, 0, 12);
 % .. and plot them
 ax4=subplot(514);
 imagesc(10*log10(spec2));
 axis xy
 
 ax5=subplot(515);
 imagesc(cep2);
 axis xy
 
 linkaxes([ax1,ax2,ax3,ax4,ax5],'x');
 
 % Notice the greater level of temporal detail compared to the 
 % RASTA-filtered version.  There is also greater spectral detail 
 % because our PLP model order is larger than the default of 8

  % Append deltas and double-deltas onto the cepstral vectors
 del = deltas(cep2);
 % Double deltas are deltas applied twice with a shorter window
 ddel = deltas(deltas(cep2,5),5);
 % Composite, 39-element feature vector, just like we use for speech recognition
 cepDpDD = [cep2;del;ddel];

