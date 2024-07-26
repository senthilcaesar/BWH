function [FtrMtx]=FindWaveFtrs2(sig,WaveLetLevel,WaveFun,fs)
%[FtrMtx]=FindWaveLetFtrs(sig,WaveLetLevel,WaveFun,TimeRes)
%sig: input signal
%WaveLetLevel: wavelet decomposition level
%WaveFun: Wavelet basis function
%FtrMtx : Extracted Features for every timeres second

[C,L] = wavedec(sig,WaveLetLevel,WaveFun);
D = detcoef(C,L,1:5);
D{6} = appcoef(C,L,WaveFun,5);%%Approximate Coeff

D2=nan(length(D{1}),length(D));
D2(:,1)=D{1};
D2(1:length(D{2}),2)=D{2};
D2(1:length(D{3}),3)=D{3};
D2(1:length(D{4}),4)=D{4};
D2(1:length(D{5}),5)=D{5};
D2(1:length(D{6}),6)=D{6};

FtrMtx(:,1:2:12)=nanmean(abs(D2));
FtrMtx(:,2:2:12)=nanvar(D2);

FtrMtx(:,13)=FtrMtx(:,1)./FtrMtx(:,3);
FtrMtx(:,14)=FtrMtx(:,3)./FtrMtx(:,5);
FtrMtx(:,15)=FtrMtx(:,5)./FtrMtx(:,7);
FtrMtx(:,16)=FtrMtx(:,7)./FtrMtx(:,9);
FtrMtx(:,17)=FtrMtx(:,9)./FtrMtx(:,11);

FtrMtx(:,18)=FtrMtx(:,1)./FtrMtx(:,5);
FtrMtx(:,19)=FtrMtx(:,3)./FtrMtx(:,7);
FtrMtx(:,20)=FtrMtx(:,5)./FtrMtx(:,9);
FtrMtx(:,21)=FtrMtx(:,7)./FtrMtx(:,11);

FtrMtx(:,22)=FtrMtx(:,1)./FtrMtx(:,7);
FtrMtx(:,23)=FtrMtx(:,3)./FtrMtx(:,9);
FtrMtx(:,24)=FtrMtx(:,5)./FtrMtx(:,11);

FtrMtx(:,25)=FtrMtx(:,1)./FtrMtx(:,9);
FtrMtx(:,26)=FtrMtx(:,3)./FtrMtx(:,11);

FtrMtx(:,27)=FtrMtx(:,1)./FtrMtx(:,11);
FtrMtx(:,28:28+length(D)-1)=nansum(abs(diff(D2)));

FtrMtx(:,34)=length(sig)/fs;
FtrMtx(:,[7 8 13 15 16 17 20 24 26 27])=[];


