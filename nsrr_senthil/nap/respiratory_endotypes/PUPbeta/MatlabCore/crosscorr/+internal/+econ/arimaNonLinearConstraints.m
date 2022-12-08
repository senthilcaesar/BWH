function [c,ceq] = arimaNonLinearConstraints(X,LagsAR,LagsSAR,LagsMA,LagsSMA,tolerance)
%ARIMANONLINEARCONSTRAINTS Nonlinear parameter constraints for ARIMA models
%
% Syntax:
%
%   [C,CEQ] = arimaNonLinearConstraints(X,LagsAR,LagsSAR,LagsMA,LagsSMA,tolerance)
%
% Description:
%
%   Enforce any nonlinear constraints needed for maximum likelihood
%   parameter estimation of conditional mean models with ARIMA components. 
%   This function serves as the nonlinear constraint function required by
%   the Optimization Toolbox function FMINCON. The nonlinear constraints
%   enforced are the stationarity and invertibility requirement of the
%   autoregressive and moving average polynomials, respectively. 
%
% Input Arguments:
%
%   X - Column vector of process parameters associated with fitting 
%     conditional mean and variance specifications to an observed series y.
%
%   LagsAR - Vector of lags associated with the non-zero coefficients of the
%     non-seasonal AR polynomial.
%
%   LagsSAR - Vector of lags associated with the non-zero coefficients of the
%     seasonal AR polynomial.
%
%   LagsMA - Vector of lags associated with the non-zero coefficients of the
%     non-seasonal MA polynomial.
%
%   LagsSMA - Vector of lags associated with the non-zero coefficients of the
%     seasonal MA polynomial.
%
%   tolerance - Positive numerical tolerance, or offset, applied to all 
%     constraints. Tolerance parameter specifies how close estimated 
%     coefficients may get to a corresponding constraint, and is related to 
%     the maximum constraint violation tolerance parameter (TolCon) of the 
%     Optimization Toolbox.
%
% Output Arguments:
%
%   C - 4-element column vector of nonlinear inequality constraints 
%     associated with the roots of each of the autoregressive and moving 
%     average polynomials. These constraints ensure that the magnitudes of 
%     the largest eigenvalues of each polynomial lie inside the unit circle.
%
%   CEQ - Empty matrix placeholder for future compatibility.
%
% Note:
%
%   This function is needed only for optimization, and is not meant to be 
%   called directly. No error checking is performed.

% Copyright 2013 The MathWorks, Inc.

%
% Determine the indices of the input vector X associated with each of the
% four polynomials.
%

nCoefficients = numel(LagsAR) + numel(LagsSAR) + numel(LagsMA) + numel(LagsSMA);

iAR  = 2:(2 + numel(LagsAR) - 1);
iSAR = (numel(iAR) + 2):(numel(iAR) + numel(LagsSAR) + 1);
iMA  = (numel(iAR) + numel(iSAR) + 2):(numel(iAR) + numel(iSAR) + numel(LagsMA) + 1);
iSMA = (numel(iAR) + numel(iSAR) + numel(iMA) + 2):(nCoefficients + 1);

%
% Now extract the relevant parameters from the input vector X and place
% them into the appropriate element of the polynomial vectors.
%

AR           = zeros(1,max(LagsAR)); 
AR(LagsAR)   = X(iAR);
SAR          = zeros(1,max(LagsSAR)); 
SAR(LagsSAR) = X(iSAR);
MA           = zeros(1,max(LagsMA)); 
MA(LagsMA)   = X(iMA);
SMA          = zeros(1,max(LagsSMA)); 
SMA(LagsSMA) = X(iSMA);

%
% Nonlinear inequality constraints ensure that the magnitudes of the 
% largest eigenvalues of the autoregressive and moving average polynomials 
% lie inside the unit circle.
%

AR   = roots([1  -AR]);
SAR  = roots([1 -SAR]);
MA   = roots([1   MA]);
SMA  = roots([1  SMA]);

c    = [max(abs(AR)).^2 ; max(abs(SAR)).^2 ; max(abs(MA)).^2 ; max(abs(SMA)).^2] - (1 - tolerance);
ceq  = [];

end