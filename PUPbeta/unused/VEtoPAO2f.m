function [F,PAO2]=VEtoPAO2f(Parameters,VI_rs,PO2_,dt,iindex,FVeupnea,Peupnea)

%dt=0.008;

%FVeupnea=Parameters(4);%-0.05

Vdeadspace = Parameters(4)*mean(VI_rs);

VI_rs=VI_rs-Vdeadspace;
VI_rs(VI_rs<0)=0;

error1=Parameters(3);%0;
k2=Parameters(2);%40;
k1=Parameters(1);%15;

Veupnea=mean(VI_rs)*(1-FVeupnea);
%Peupnea=Parameters(5);
Pexpected=150-(150-Peupnea)*(1-FVeupnea);

PAO2=0*VI_rs;
PAO2(1)=Pexpected+error1;

ODC=1;
if ODC==1
% p50 = 27;
% p90 = 60; %60 / 2.2350
% beta = (9.*(p50.^3)-p90.^3)./(p90-(9.*p50)); 
% alpha = p50.^3+beta.*p50;
alpha=23400;
beta=150;
elseif ODC==3 %Hill
    alpha = 50;
    beta=5;
end

AaDO2=0;

if 0 
figure(12); 
Ptemp=50:140;
Stemp=P_O2toS_O2(ODC,Ptemp-AaDO2,alpha,beta);
plot(Ptemp,Stemp)
end

Seupnea=P_O2toS_O2(ODC,Peupnea-AaDO2,alpha,beta);

%QbonQdot=60; Fa=0.1; SaSvdelay = Fa*QbonQdot; SaSvdelayi = round(SaSvdelay/dt); tauSv = (1-Fa)*QbonQdot;
%DSv_O2(1)=P_O2toS_O2(ODC,Pexpected-AaDO2,alpha,beta)-Seupnea;

for i=1:length(VI_rs)
    if 1
        if i<length(VI_rs)
            PAO2(i+1)=PAO2(i)+dt*(k1*(VI_rs(i)-Veupnea)-k2*(PAO2(i)-Peupnea));
        end
    elseif 0
        %Sa_O2(i)=P_O2toS_O2(1,PAO2(i),23400,150);
        if i<length(VI_rs)
            PAO2(i+1)=PAO2(i)+dt*((k1/50)*(VI_rs(i)-Veupnea)*(150-PAO2(i)) - k2*(PAO2(i)-Peupnea));  %k (or k=50) is VL*P0; tau = VL/(Qdot*beta);
        end
    elseif 0
        Sa_O2(i)=P_O2toS_O2(ODC,PAO2(i)-AaDO2,alpha,beta);
        if i<length(VI_rs)
            if i>SaSvdelayi
                DSv_O2(i+1)=DSv_O2(i)+dt/tauSv*(Sa_O2(i-SaSvdelayi)-Seupnea-DSv_O2(i));
                %Sa_O2_mixed(i+1)=Sa_O2_mixed(i)+dt/tauSv*(Sa_O2(i-SaSvdelayi));
            else
                DSv_O2(i+1)=DSv_O2(1);
            end
            PAO2(i+1)=PAO2(i)+dt*((k1/50)*(VI_rs(i)-Veupnea)*(150-PAO2(i)) - k2*(Sa_O2(i)-Seupnea+DSv_O2(i)));  
        end 
    else
        Sa_O2(i)=P_O2toS_O2(ODC,PAO2(i)-AaDO2,alpha,beta);
        if i<length(VI_rs)
            PAO2(i+1)=PAO2(i)+dt*((k1)*(VI_rs(i)-Veupnea)*(150-PAO2(i)) - k2*(Sa_O2(i)-Seupnea));  
        end
    end
end
%Sa_O2=P_O2toS_O2(ODC,PAO2-AaDO2,alpha,beta);
%Sa_O2_delayed=([NaN*zeros(delay_i,1);Sa_O2(1:(end-delay_i))]);
%Sa_O2_delayed = [Sa_O2(1)*ones(delay_i,1);Sa_O2(1:(end-delay_i))];

%delay_maxi = round(30/dt);
%F=sum((Sa_O2_delayed(delay_maxi:end)-SaO2(delay_maxi:end)).^2)/length((Sa_O2_delayed(delay_maxi:end)));

Err = PO2_(iindex)-PAO2(iindex);
F=sum((Err(1:end)).^2);




