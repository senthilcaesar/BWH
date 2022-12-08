function S_O2 = P_O2toS_O2(ODC,P_O2,alpha,beta)
%%%%%% Oxygen Saturation %%%%%%
if ODC==1
%sat_alpha=23400;
%sat_beta=150;
S_O2 = 100*(alpha.*(P_O2.^3+beta.*P_O2).^-1+1).^-1;
%Severinghaus, J. W. (1979). "Simple, accurate equations for human blood O2 dissociation computations." J Appl Physiol 46(3): 599-602. 
%http://www.health.adelaide.edu.au/paed-anaes/javaman/Respiratory/oxygen/O2Equations.html
elseif ODC==2
%k1=alpha, k3=beta
temp = ((P_O2/(10.^alpha)).^(1/beta));
S_O2 = 100*(temp./(temp+1));
clear temp;     
%Battaglia
elseif ODC==3
    S_O2 = 100*((P_O2./alpha).^beta)./(1+(P_O2./alpha).^beta); 
end