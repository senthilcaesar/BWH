function [VT,VTi,VTe,VI,VE] = VEfromFlowWithoutBreathDetection(Flow, BB_i_start,BB_i_mid,BB_i_end,Ti,Te,Apnea_B,dt)

I = [BB_i_start; BB_i_mid; BB_i_end(end)]; %zero flow indices
Flow2 = Flow - nanmedian(Flow(I)); % remove leak
vol=cumsum(Flow2)*dt;
%figure(89); plot(SigT.Time,SigT.FlowOral);
VTi = vol(BB_i_mid(:)) - vol(BB_i_start(:)); 
VTe = vol(BB_i_mid(:)) - vol(BB_i_end(:));
%Ti = (BB_i_mid-BB_i_start)'*dt;
%Te = (BB_i_end-BB_i_mid)'*dt;
VT = (VTi.*Te(:)+VTe.*Ti(:))./(Ti(:)+Te(:));
VTe(VTe<0)=0;
VTi(VTi<0)=0;
VT(VT<0)=0;
VI = VTi./(Ti(:)+Te(:));
VE = VTe./(Ti(:)+Te(:));
VTe(Apnea_B(:)==1)=0;
VTi(Apnea_B(:)==1)=0;
VI(Apnea_B(:)==1)=0;
VE(Apnea_B(:)==1)=0;
VT(Apnea_B(:)==1)=0;