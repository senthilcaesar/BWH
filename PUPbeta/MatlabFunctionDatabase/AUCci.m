function [ci,se,p] = AUCci(AUC,N1,N2)
% Written by S Sands 2016
% See Hanley & McNeil, 1982; Cortex & Mohri, 2004
% Confidence intervals
invertAUC = AUC<0.5;
AUC(invertAUC)=1-AUC(invertAUC);
alpha = 0.05;
z = norminv(1-alpha/2);
Q1 = AUC./(2-AUC);
Q2 = (2*AUC.^2)./(1+AUC);
se = sqrt((AUC.*(1-AUC) + (N1-1).*(Q1-AUC.^2) + (N2-1).*(Q2-AUC.^2))./(N1.*N2));
%AUC(invertAUC)=1-AUC(invertAUC); 
ci = [AUC(:)-z*se(:) AUC(:)+z*se(:)];
p = [1-normcdf(abs(0.5-AUC)./se,0,1)];