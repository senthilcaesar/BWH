function [AR,beta,constant,variance,residuals] = arx0(OBJ, YData, XData)
%ARX0 Initial coefficient estimates of univariate ARX models
%
% Syntax:
%
%   [AR, beta, constant, variance, residuals] = arx0(OBJ, Y, X)
%
% Description:
%
%   Estimate the coefficients of the stable seasonal and non-seasonal
%   autoregressive polynomials associated with univariate ARX models via 
%   ordinary least squares (OLS). The models may include both seasonal and 
%   non-seasonal integration effects as well as exogenous regression 
%   components; moving average components, if any, are ignored. If necessary, 
%   an additive constant is also estimated.
%
% Input Arguments:
% 
%   OBJ - An ARIMA model whose seasonal and non-seasonal AR coefficients 
%     are estimated. An additive constant and any coefficients associated 
%     with an exogenous regression component are also estimated.
%
%   Y - Time series vector of interest, assumed to follow an ARX process 
%     specified in OBJ. The last element contains the most recent 
%     observation. 
%
%   X - Matrix of exogenous data used to include a regression component. Each
%     column of X is a separate time series, and the last row of each contains 
%     the most recent observation of each series. If missing, the conditional 
%     mean will have no regression component regardless of the presence of 
%     any regression coefficients found in the model.
%
% Output Arguments:
%
%   AR - Vector of coefficients associated with the stable autoregressive 
%     polynomials known to the optimizer. These coefficients include all 
%     non-zero autoregressive coefficients included in the model at lags 
%     greater than zero (i.e., coefficients at positive lags of the 'AR' 
%     and 'SAR' polynomials, with the 'AR' coefficients followed by the 'SAR' 
%     coefficients). Both seasonal and non-seasonal integration effects are 
%     accounted for by data adjustment, and coefficients associated with 
%     integration are NOT included in AR because the optimizer knows 
%     nothing about them.
%
%   beta - Vector of coefficients associated with the columns of the input
%     exogenous data matrix X (see above).
%
%   constant - Model constant. 
%
%   variance - Variance of OLS residuals.
%
%   residuals - Residuals from the OLS fit. These are returned to support
%     subsequent analysis, such as estimating a moving average component.
%
% Notes: 
%
%  o  Regression components included in the model are based on the presence 
%     of the exogenous data matrix X, and the number of elements in the 
%     output regression coefficient vector (see beta above) is equal to the
%     number of columns in the data matrix X.
%
%     Provided the input data matrix X is not empty, any non-NaN coefficients 
%     found in OBJ.Beta indicate existing initial regression coefficient
%     estimates. Such coefficients are not estimated, but simply copied 
%     directly into the output beta vector. NaN-valued coefficients, however, 
%     are estimated in the presence of any initial estimates found in the 
%     model OBJ.
%
%   o Any non-zero and non-NaN coefficients at positive lags specified in 
%     OBJ.AR and OBJ.SAR indicate existing initial estimates. Such 
%     coefficients are not estimated, but simply copied directly into the 
%     output AR vector. NaN-valued coefficients, however, are estimated in 
%     the presence of any initial estimates found in the model OBJ.
%
%   o For autoregressive components, if an OLS coefficient cannot be uniquely 
%     associated with a single coefficient, then the corresponding element 
%     of the AR output vector is set to zero (see detailed comments below).
%
%   o If the autoregressive estimates indicate one or more unstable 
%     polynomials, then all estimated coefficients of the output vectors AR 
%     and beta are set to zero; user-specified values, however, are honored.
%
%   o Although user-specified initial estimates are honored for an additive 
%     constant, for all AR and SAR coefficients, and for all coefficients 
%     associated with a regression component, a user-specified initial 
%     estimate of the model variance is not; in all cases the output variance 
%     is the sample variance of the OLS residuals.

% Copyright 2013 The MathWorks, Inc.

%
% Filter Y(t) to remove the effects of integration.
%

I     = getLagOp(OBJ, 'Integrated Non-Seasonal') * getLagOp(OBJ, 'Integrated Seasonal');
YData = I(YData);

%
% Since the exogenous data array takes precedence over any regression
% coefficients, the OBJ.Beta property might be empty. The following ensures 
% that the regression coefficient vector is of correct length.
%

isRegressionIncluded = ~isempty(XData);

if isRegressionIncluded
   if isempty(OBJ.Beta)
      Beta = nan(1,size(XData,2));   % When empty, estimate all of them
   else
      Beta = OBJ.Beta;
      Beta = Beta(:)';               % Guarantee a row vector
   end
%
%  Adjust the number of observations of X(t) to match those of the filtered Y(t).
%
   XData = XData((end - size(YData,1) + 1):end,:);
else
   Beta = [];
end

%
% Access the non-seasonal/seasonal AR polynomials, and determine the
% non-zero lags.
% 

LagOpAR        = getLagOp(OBJ, 'AR'); 
LagOpSAR       = getLagOp(OBJ, 'SAR');

LagsNonZeroAR  = LagOpAR.Lags;                       % Non-zero lags of AR
LagsNonZeroSAR = LagOpSAR.Lags;                      % Non-zero lags of SAR
LagsNonZeroAR  = LagsNonZeroAR(LagsNonZeroAR > 0);   % Exclude lag 0
LagsNonZeroSAR = LagsNonZeroSAR(LagsNonZeroSAR > 0); % Exclude lag 0

AR = LagOpAR * LagOpSAR;          % Product polynomial
L  = AR.Lags;                     % Retain all non-zero lags of the product
AR = cell2mat(toCellArray(AR))';  % Convert to a vector
AR = AR(L + 1);                   % Remove intermediate zeros, and avoid index 0

%
% Estimate any unknown AR and regression coefficients by OLS.
%
% Form the design matrix X(t) with the lags of Y(t) associated with unknown 
% AR coefficients and unknown regression coefficients. Then form a new time 
% series as the linear combination of the input Y(t) and any lags of Y(t) 
% associated with known AR coefficients and known regression coefficients.
%

if isnan(OBJ.Constant)
%
%  The constant is estimated, so allow for it in the design matrix X(t).
%
   Y = lagmatrix(YData, L(~isnan(AR))) * AR(~isnan(AR));

   if isRegressionIncluded
      X = [ones(numel(YData),1)  lagmatrix(YData, L(isnan(AR)))  XData(:,isnan(Beta))];
      Y = Y - XData(:,~isnan(Beta)) * Beta(1,~isnan(Beta))';
   else
      X = [ones(numel(YData),1)  lagmatrix(YData, L(isnan(AR)))];
   end

else
%
%  The constant is known, so allow for it in the adjusted time series matrix.
%
   if isRegressionIncluded
      X = [lagmatrix(YData, L(isnan(AR)))   XData(:,isnan(Beta))];
      Y = lagmatrix(YData, L(~isnan(AR))) * AR(~isnan(AR)) - OBJ.Constant - XData(:,~isnan(Beta)) * Beta(1,~isnan(Beta))';
   else
      X = lagmatrix(YData, L(isnan(AR)));
      Y = lagmatrix(YData, L(~isnan(AR))) * AR(~isnan(AR)) - OBJ.Constant;
   end
end

Last  = L(end);                    % Row of the last NaN found in X(t) or Y(t)
[Q,R] = qr(X((Last + 1):end,:),0);
b     = R \ (Q' * Y((Last + 1):end));

residuals         = Y - X * b;
residuals(1:Last) = 0;

%
% Pack the regression coefficients and additive constant into the appropriate outputs.
%

beta              = Beta(:);
beta(isnan(Beta)) = b((end - sum(isnan(Beta)) + 1):end);

if isnan(OBJ.Constant)
   constant = b(1);
   b        = b(2:end);            % Exclude the constant from the coefficient vector
else
   constant = OBJ.Constant;
end

%
% Assign the coefficients from the estimated regression vector to the 
% appropriate lags of the AR and SAR polynomials.
%

L          =  L(isnan(AR));               % Retain only lags associated with estimated coefficients
AR         =  LagOpAR.Coefficients;
AR         = -[AR{1:LagOpAR.Degree}];     % Returns positive lags 1, 2, 3, ...
SAR        =  LagOpSAR.Coefficients;
SAR        = -[SAR{1:LagOpSAR.Degree}];   % Returns positive lags 1, 2, 3, ...
LagsNanAR  =  find(isnan(AR));            % AR  lags assigned to 
LagsNanSAR =  find(isnan(SAR));           % SAR lags assigned to 
[~,  iAR]  =  intersect(L, LagsNanAR);    % AR  indices assigned from the regression vector
[~, iSAR]  =  intersect(L, LagsNanSAR);   % SAR indices assigned from the regression vector

AR (LagsNanAR)  = b(iAR);
SAR(LagsNanSAR) = b(iSAR);

%
% Determine which lags of the regression vector are uniquely associated with 
% a single coefficient:
%
%   o If the non-seasonal and seasonal autoregressive polynomials have any 
%     non-zero coefficients associated with positive lags in common (i.e., 
%     the intersection of non-zero coefficient lags of LagOpAR and LagOpSAR 
%     is a non-empty set), then these common lags will have an estimated 
%     regression coefficient which the sum the polynomial coefficients of 
%     interest; furthermore, the coefficients at these common lags will also 
%     appear as a multiplicative regression coefficient at twice the common 
%     lags.
%
%   o If the combined non-zero lags of the non-seasonal and seasonal 
%     polynomials (i.e., the union of non-zero coefficient lags of LagOpAR 
%     and LagOpSAR) have any lags in common with the higher-order lags 
%     produced as a result of their multiplication (i.e., the intersection 
%     of the union of the lags associated non-zero coefficients of LagOpAR 
%     and LagOpSAR and the higher-order lags associated with the polynomial 
%     multiplication is a non-empty set), then these common lags will have 
%     estimated regression coefficients which are the sum of various 
%     polynomial coefficients of interest.
%
% If either of the above intersections is a non-empty set, then the common 
% lags of each intersection are associated with polynomial coefficients not 
% uniquely associated with a single coefficient and the corresponding output
% coefficients are set to zero.
%

LagsUnion   = union(LagsNonZeroAR, LagsNonZeroSAR);
LagsProduct = zeros(numel(LagsNonZeroSAR), numel(LagsNonZeroAR));

if ~isempty(LagsProduct)
   for i = 1:numel(LagsNonZeroSAR)
       LagsProduct(i,:) = LagsNonZeroSAR(i) + LagsNonZeroAR;  % Lags due to multiplication
   end
end

LagsIntersect1  = intersect(LagsNonZeroAR, LagsNonZeroSAR);
LagsIntersect2  = intersect(LagsUnion, LagsProduct(:)');
LagsUnion       = union(LagsIntersect1, LagsIntersect2);
LagsToZero      = intersect(LagsUnion, LagsNanAR);          % Honor user-specified values
AR (LagsToZero) = 0;
LagsToZero      = intersect(LagsUnion, LagsNanSAR);         % Honor user-specified values
SAR(LagsToZero) = 0;

%
% Ensure stationarity of autoregressive processes. 
%
% If the autoregressive estimates indicate one or more unstable processes, 
% then simply assume all estimated coefficients are zero while honoring 
% user-specified values; otherwise, pack the coefficients into the appropriate 
% element of the output vector.
%
% Notice that the stationarity test is offset from one by a small tolerance, 
% consistent with the stability tolerance of the LagOp/isStable method.
%

if any(abs([roots([1 -AR]) ; roots([1 -SAR])]) >= (1 - 10*eps))
   AR0               =  LagOpAR.Coefficients;
   AR0               = -[AR0{1:LagOpAR.Degree}];
   SAR0              =  LagOpSAR.Coefficients; 
   SAR0              = -[SAR0{1:LagOpSAR.Degree}];
   AR                =  [AR0(LagsNonZeroAR) SAR0(LagsNonZeroSAR)]';
   AR(isnan(AR))     =  0;
   beta(isnan(Beta)) =  0;  
   variance          =  var(Y(~isnan(Y)));
else
   AR       = [AR(LagsNonZeroAR) SAR(LagsNonZeroSAR)]';
   variance = var(residuals((Last + 1):end));
end

end