function [Xlimited,I] = limitoutliers(X,coef1,coef2,ploton)
% 
% X = Ts.ArTH_bmi_BL;
% coef1=3;
% coef2=84.13; %1 SD from mean equivalent
% ploton = 1;

%y = cdf(makedist('Normal','mu',0,'sigma',1),[-2,-1,0,1,2])

A = prctile(X,[(100-coef2) 50 coef2]);

upper = (A(3)-A(2))*coef1 + A(2);
lower =  A(2) - (A(2)-A(1))*coef1;

I = find(X<lower | X>upper);
%length(I)


if ploton>0
figure(ploton);
histogram(X,30);
end

X(X<lower) = lower;
X(X>upper) = upper;
Xlimited = X;

if ploton>0
hold('on');
histogram(Xlimited,30);
end

% figure(1)
% histogram(X,30)





