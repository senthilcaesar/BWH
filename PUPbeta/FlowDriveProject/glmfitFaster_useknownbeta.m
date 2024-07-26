function [Rsq,p,MSE,ypred]=glmfitFaster_useknownbeta(X,y,w,remnans, beta)

%written by SS
%due credit to https://www.mathworks.com/matlabcentral/fileexchange/34394-fast---detailed-multivariate-ols-regression
%...and author of Matlab's glmfit

% this version does not determine beta, instead using the beta passed in.

verbose = 0;

if remnans
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

%beta = ((X'*(w.*X))\X')*(w.*y); %checked 
sw = w.^0.5;

ypred = (beta' * X')';
SSE = sum((sw.*(ypred - mean(y))) .^ 2);
SST = sum((sw.*(y - mean(y))).^2);
SSR = SST - SSE; %"dev" in glmfit
SSR_2 = sum((sw.*(y - ypred)).^2); % alternate method of calculating this value, then check equivalence
if ~isequal(round(SSR,3), round(SSR_2,3))&&verbose; disp('Unequal SSresiduals in fastfit'); end
Rsq = SSE / SST;
Rsq_2 = 1-(SSR / SST); % alternate method of calculating this value, then check equivalence
if ~isequal(round(Rsq,3), round(Rsq_2,3))&&verbose; disp('Unequal Rsquared in fastfit'); end
MSE = SSR ./ dfe;
%RMSE = MSE .^ 0.5; %"sfit" in glmfit
F = (SSE ./ k) ./ (SSR ./ dfe);
Fp = 1 - fcdf(F,dfe,1);
covb = MSE * (inv((sw.*X)' * (sw.*X)));
se = sqrt(diag(covb));
t = beta./se;
p = 2*tcdf(-abs(t),dfe);

if 0 %compare with Matlab function
[b,dev,stats]=glmfit(X(:,2:end),y,'normal','weights',w,'estdisp','off');
end

%
% DLM
% weighted coefficient of determination
%     SSres_w = sum(weights.*(Gtest_All - predyL1O_array(:,NumFtrs)).^2);         % 
%     SStot_w = sum(weights.*(Gtest_All - mean(Gtest_All)).^2);                   % SST
%     SSreg_w = sum(weights.*(predyL1O_array(:,NumFtrs)-mean(Gtest_All)).^2);     % SSE
%     Rsquared_1_w = 1 - (SSres_w/SStot_w)
%     Rsquared_2_w = (SSreg_w/SStot_w)
%     SStot2_w = SSres_w+SSreg_w;
% 




