function P_O2 = S_O2toP_O2(ODC,S_O2,alpha,beta)

if ODC==1
%Saturation to Partial Pressure
%sat_alpha=23400;
%sat_beta=150;
%beta=150;
%alpha=23400;
temp=((-108.*S_O2.*(alpha)+12*(3^0.5).*(4*(beta)^3.*S_O2.^2-800*(beta)^3.*S_O2+40000*(beta)^3+27.*S_O2.^2.*(alpha)^2).^0.5).*(S_O2-100).^2).^(1/3);
P_O2=1./(6.*(S_O2-100)).*temp+(-2*(beta).*(S_O2-100))./temp;
elseif ODC==2
        P_O2=10.^(alpha + beta*log10(S_O2./(100-S_O2)));  
end