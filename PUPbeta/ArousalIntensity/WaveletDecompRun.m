function [D2]=WaveletDecompRun(sig,decomLevel,wavefun)

[C,L] = wavedec(sig,decomLevel,wavefun);
D = detcoef(C,L,1:decomLevel);
D{6} = appcoef(C,L,wavefun,decomLevel);%%Approximate Coeff

D2=nan(length(D{1}),length(D));
D2(:,1)=D{1};
D2(1:length(D{2}),2)=D{2};
D2(1:length(D{3}),3)=D{3};
D2(1:length(D{4}),4)=D{4};
D2(1:length(D{5}),5)=D{5};
D2(1:length(D{6}),6)=D{6};







