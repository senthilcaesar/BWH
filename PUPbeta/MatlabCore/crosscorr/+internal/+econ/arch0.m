function [K, ARCH] = arch0(e, Q)
%ARCH0 Initial parameter estimates of univariate ARCH(Q) processes
%
% Syntax:
%
%   [Constant, ARCH] = arch0(e, Q)
%
% Description:
%
%   Compute initial parameter estimates of a univariate ARCH(Q) model via 
%   ordinary least squares,
%
%      s^2(t) = c + A1*e^2(t-1) + ... + Aq*e^2(t-Q).
%
% Input Arguments:
%
%   e - Time series vector of interest. e is a stochastic process assumed
%     to follow an ARCH(Q) form. The last element contains the most recent
%     observation. 
%
%   Q - Nonnegative integer, the degree of the ARCH(Q) model.
%
% Output Arguments:
%
%   Constant - Estimate of the constant of an ARCH(Q) model. 
%
%   ARCH - Q-by-1 vector of estimates of the Q coefficients of an ARCH(Q) 
%     model. The first element of ARCH is an estimate of the coefficient of 
%     e(t-1)^2, the second is an estimate of the coefficient of e(t-2)^2, 
%     and so forth.
%
% References:
%
%   [1] Box, G. E. P., G. M. Jenkins, and G. C. Reinsel. Time Series
%       Analysis: Forecasting and Control. 3rd edition. Upper Saddle River,
%       NJ: Prentice-Hall, 1994.
%
%   [2] Greene, W. H. Econometric Analysis, Prentice Hall, 3rd Edition, 1997
%
%   [3] Hamilton, J. D. Time Series Analysis. Princeton, NJ: Princeton
%       University Press, 1994.

% Copyright 1999-2011 The MathWorks, Inc.

%
% Square the input residuals and compute lagged values.
%

e2    = e.^2;
e2Lag = lagmatrix(e2, 1:Q);

%
% Compute initial ARCH(Q) model coefficients via OLS.
%

X = [ones(size(e2,1) - Q,1) e2Lag((Q + 1):end,:)];
b = X \ e2((Q + 1):end);

%
% Now update via FGLS (see Greene, page 571-572).
%

f = X * b;
X = [1./f   e2Lag((Q + 1):end,:) ./ f(:,ones(1,Q))];
y = (e2((Q + 1):end) ./ f) - 1;
b = b + (X \ y);

%
% Ensure positivity of all coefficients.
%

b(b < 0) = 0;

% 
% Ensure positivity of the constant K and stationarity of ARCH coefficients. 
%
% If the sum of the ARCH coefficients is not less than 1, assume the sum 
% A1 + A2 + ... + AQ = 0.95 and normalize the coefficients.
%

K    = max(b(1), 1e-6);                 % Ensure the constant is > 0.
ARCH = b(2:(Q + 1));

if sum(ARCH) >= 1
   ARCH =  0.95 * (ARCH / sum(ARCH));   % This maintains proportionality.
end

end