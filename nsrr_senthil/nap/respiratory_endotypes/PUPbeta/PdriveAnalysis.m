function     [BreathDataTable2,alpha,scalealpha,scalealpha30,scalealpha10]= PdriveAnalysis(BreathDataTable2,RespT,settings)
        %[coeff,~,~]=xlsread('E:\Dropbox (Partners HealthCare)\PhenotypeDrive2018\PUPstart\AMasterSpreadsheet','BO4:BO63');
        alpha2=29;
        %find coeff
       % BreathDataTable2.DeltaEdi2=abs(BreathDataTable2.PeakEdi)
       % [BreathDataTable2,DriveExponent] = VdriveFromEdi(BreathDataTable2,'VE','DeltaEdi2',1);
        [~,EnsemblesX]=EventAnalysisRun([],BreathDataTable2,RespT,settings.BreathEnsembleMethod,{'DeltaPes','VT','DeltaEdi','VdriveEdiNorm','PeakEdi','VdriveDeltaEdi2Norm'})
        
        Pes=EnsemblesX.DeltaPes;
        Edi=EnsemblesX.VdriveEdiNorm;
        VT=EnsemblesX.VT;
        Time=EnsemblesX.Time;
        win1=1*(Time>-15&Time<0); win1(win1==0)=NaN;
        [~,a]=min(VT.*win1);
        win2=1*(Time>0&Time<15); win2(win2==0)=NaN;
        [~,b]=max(VT.*win2);
        Time([a b])
        alpha=(Pes(a)*Edi(b)-Pes(b)*Edi(a))./(Edi(a)*VT(b)-VT(a)*Edi(b));
        alpha(alpha<0.33*alpha2)=0.33*alpha2;
        alpha(alpha>3*alpha2)=3*alpha2;
        scalealpha = Edi(a)./(Pes(a)+alpha*VT(a));
        scalealpha30 = Edi(a)./(Pes(a)+30*VT(a));
        scalealpha10 = Edi(a)./(Pes(a)+10*VT(a));
        
        
        
         BreathDataTable2.PdriveK=(BreathDataTable2.DeltaPes+alpha.*BreathDataTable2.VT)*scalealpha;  
        BreathDataTable2.Pdrive30=(BreathDataTable2.DeltaPes+30.*BreathDataTable2.VT)*scalealpha30;
        BreathDataTable2.Pdrive10=(BreathDataTable2.DeltaPes+10.*BreathDataTable2.VT)*scalealpha10;
  
        BreathDataTable2.PdriveKraw=(BreathDataTable2.DeltaPes+alpha.*BreathDataTable2.VT);
       
        % [BreathDataTable2,DriveExponent] = VdriveFromEdi(BreathDataTable2,'VE','PeakEdi',1);
       % Edi=abs(EnsemblesX.PeakEdi);
        %alpha=(Pes(a)*Edi(b)-Pes(b)*Edi(a))./(Edi(a)*VT(b)-VT(a)*Edi(b));
%         alpha(alpha<0.33*alpha2)=0.33*alpha2;
%         alpha(alpha>3*alpha2)=3*alpha2;
        
        
        [BreathDataTable2] = VdriveFromEdi(BreathDataTable2,'VE','DeltaPes',0);
        [BreathDataTable2] = VdriveFromEdi(BreathDataTable2,'VE','Pdrive30',0);
        [BreathDataTable2] = VdriveFromEdi(BreathDataTable2,'VE','PdriveKraw',0);
   
        [~,EnsemblesX]=EventAnalysisRun([],BreathDataTable2,RespT,settings.BreathEnsembleMethod);
   
       
    figure(102); clf(102); set(gcf,'color',[1 1 1]);
    Nsubs=5;
   ax102(1)=subplot(Nsubs,1,1); 
     plot(EnsemblesX.Time,EnsemblesX.DeltaPes);
     hold on
     plot(EnsemblesX.Time([a b]),EnsemblesX.DeltaPes([a b]),'.','markersize',20);
     plot(EnsemblesX.Time,EnsemblesX.PdriveKraw);
     plot(EnsemblesX.Time([a b]),EnsemblesX.PdriveKraw([a b]),'.','markersize',20);
     box off
   ax102(2)=subplot(Nsubs,1,2); 
     plot(EnsemblesX.Time,EnsemblesX.VT);
     hold on
     plot(EnsemblesX.Time([a b]),EnsemblesX.VT([a b]),'.','markersize',20);
     box off
   ax102(3)=subplot(Nsubs,1,3); 
     plot(EnsemblesX.Time,EnsemblesX.VdriveEdiNorm);
     hold on
     plot(EnsemblesX.Time([a b]),EnsemblesX.VdriveEdiNorm([a b]),'.','markersize',20);
     box off  
     plot(EnsemblesX.Time,EnsemblesX.PdriveK);
     
   ax102(4)=subplot(Nsubs,1,4); 
     plot(EnsemblesX.Time,EnsemblesX.PdriveK);
     hold on
     plot(EnsemblesX.Time([a b]),EnsemblesX.PdriveK([a b]),'.','markersize',20);
     box off  
   ax102(5)=subplot(Nsubs,1,5); 
     plot(EnsemblesX.Time,EnsemblesX.VdriveEdiNorm);
     hold on
     plot(EnsemblesX.Time([a b]),EnsemblesX.VdriveEdiNorm([a b]),'.','markersize',20);
     plot(EnsemblesX.Time,EnsemblesX.PdriveK);
     plot(EnsemblesX.Time([a b]),EnsemblesX.PdriveK([a b]),'.','markersize',20);
     plot(EnsemblesX.Time,EnsemblesX.VdrivePdriveKrawNorm);
     plot(EnsemblesX.Time([a b]),EnsemblesX.VdrivePdriveKrawNorm([a b]),'.','markersize',20);
     
%      plot(EnsemblesX.Time,EnsemblesX.VdrivePeakEdiNorm);
%      plot(EnsemblesX.Time([a b]),EnsemblesX.VdrivePeakEdiNorm([a b]),'.','markersize',20);
%      plot(EnsemblesX.Time,EnsemblesX.VdrivePesNorm);
%      plot(EnsemblesX.Time([a b]),EnsemblesX.VdrivePesNorm([a b]),'.','markersize',20);
     box off  
   
   end