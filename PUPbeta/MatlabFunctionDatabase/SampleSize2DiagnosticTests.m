function SampleSize2DiagnosticTests(varargin)

%written by S.Sands 2017-02-17 based on
%http://www.sciencedirect.com/science/article/pii/S1532046414000501
%assumes normal binomial distribution

P1=0.6; %sensitivity or specificity of current diagnostic test.
deltaP=0.2; %degree of improvement that is expected of new diagnostic test
N=40;

alpha=0.05;
beta=0.8;
zbeta=norminv(beta); %corresponds to 80% power
zalpha=norminv(1-alpha/2); %1.96 corresponds to 95% confidence

%calculate 
P2=P1+deltaP; %sensitivity or specificity of new diagnostic test
Pave = (P1+P2)/2; %average
n = (zalpha*(2*Pave*(1-Pave))^0.5 + zbeta*(P1*(1-P1)+P2*(1-P2))^0.5)/(deltaP^2);

zBeta = (N*(deltaP^2) - zalpha*(2*Pave*(1-Pave))^0.5)/((P1*(1-P1)+P2*(1-P2))^0.5);
Beta=normcdf(zBeta);