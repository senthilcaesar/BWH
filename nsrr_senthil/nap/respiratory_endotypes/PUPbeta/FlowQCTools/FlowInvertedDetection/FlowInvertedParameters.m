function [diffcurr,diffprev,COVVTi] = FlowInvertedParameters(VTi,VTe,VT,Apnea_B)

VTeprev = VTe(1:end-1);
VTicurr = VTi(2:end);
VTecurr = VTe(2:end);
temp = Apnea_B(2:end);
VTeprev(temp==1)=[];
VTicurr(temp==1)=[];
VTecurr(temp==1)=[];
diffcurr = nanmedian(abs((VTicurr-VTecurr)))/nanmean(VT);
diffprev = nanmedian(abs((VTicurr-VTeprev)))/nanmean(VT);
COVVTi = nanstd(VTicurr(2:end)-VTicurr(1:end-1))/nanmean(VTicurr);
