function [FtrMtx]=WaveletFtrsRun(D2)
%[FtrMtx]=FindWaveLetFtrs(D2)
%FtrMtx : Extracted Features for every timeres second

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


