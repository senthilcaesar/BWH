function [y]= fVpassiveTmesa(x)

%setup
tol=0.5;
kVpassive=0.33;
maxV = (1-(1-((100-tol)/100)).^kVpassive) / ((100-tol)/100);
F1.Link = @(x) 100*(1-(1-x/100).^(kVpassive))/maxV;
%x is Vpassive (0-100)
x(x>=100-tol) = 100-tol;
x(x<=tol) = tol;
x = x/100; %x is clipped Vpassive 0.005-0.995
%y = F1.Link(x*100)/100;
y = F1.Link(x*100);