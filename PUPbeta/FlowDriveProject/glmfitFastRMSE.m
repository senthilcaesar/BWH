function [RMSE]=glmfitFastRMSE(X,y,w)

%written by SS
%due credit to https://www.mathworks.com/matlabcentral/fileexchange/34394-fast---detailed-multivariate-ols-regression
%...and author of Matlab's glmfit
if 0
nanx = sum(isnan(X),2)>0;
nany = isnan(y);
X(nanx|nany,:)=[];
y(nanx|nany)=[];
w(nanx|nany)=[];
w=w/mean(w);
end

[n,k] = size(X);
X = [ones(n,1),X]; % Add column of ones 
dfe = n - k - 1; % 

beta = ((X' * (w.*X)) \ X') * (w.*y); %checked
sw = w.^0.5;

Estimates = (beta' * X')';
SSE = sum((sw.*(Estimates - mean(y))) .^ 2);
SST = sum((sw.*(y - mean(y))).^2);
SSR = SST - SSE; %"dev" in glmfit
Rsq = SSE / SST;
MSE = SSR ./ dfe;
RMSE = MSE .^ 0.5; %"sfit" in glmfit
% F = (SSE ./ k) ./ (SSR ./ dfe);
% Fp = 1 - fcdf(F,dfe,1);
% covb = MSE * (inv((sw.*X)' * (sw.*X)));
% se = sqrt(diag(covb));
% t = beta./se;
% p = 2*tcdf(-abs(t),dfe);
% 
% if 0 %compare with Matlab function
% [b,dev,stats]=glmfit(X(:,2:end),y,'normal','weights',w,'estdisp','off');
% end
%





