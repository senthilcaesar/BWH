function [AR,MA,constant,variance] = arma0(y, LagOpAR, LagOpMA)
%ARMA0 Initial parameter estimates of univariate ARMA(P,Q) processes
%
% Syntax:
%
%   [AR,MA,constant,variance] = arma0(y, LagOpAR, LagOpMA)
%
% Description:
%
%   Compute initial estimates of the autoregressive (AR) and moving
%   average (MA) coefficients of a stationary/invertible univariate
%   ARMA(P,Q) model of y(t),
%
%      y(t) = c + A1*y(t-1) + ... + Ap*y(t-P) + e(t) 
%               + M1*e(t-1) + ... + Mq*e(t-Q).
%
%   Estimates of the model constant c and the variance of the noise 
%   process e(t) are also provided.
%
% Input Arguments:
%
%   y - Time series vector of interest. y is a stochastic process assumed
%     to follow a stationary/invertible ARMA(P,Q) model (see notes). The 
%     last element contains the most recent observation. 
%
%   LagOpAR - The autoregressive (AR) lag operator polynomial of degree P 
%     associated with a stationary/invertible ARMA(P,Q) model (see notes). 
%
%   LagOpMA - The moving average (MA) lag operator polynomial of degree Q 
%     associated with a stationary/invertible ARMA(P,Q) model (see notes). 
%
% Output Arguments:
%
%   AR - P-by-1 vector of initial estimates of the autoregressive parameters
%     of an ARMA(P,Q) model. The first element of AR is an estimate of the 
%     coefficient of the first lag of y, the second is an estimate of the 
%     coefficient of the second lag of y, and so forth. Only elements of AR 
%     at lags associated with NaN coefficients found in LagOpAR are estimated; 
%     all others are held fixed.
%      
%   MA - Q-by-1 vector of initial estimates of the moving-average parameters 
%     of an ARMA(P,Q) model. The first element of MA is an estimate of the 
%     coefficient of the first lag of the innovations process, the second 
%     is an estimate of the coefficient of the second lag of the innovations 
%     process, and so forth. Only elements of AR at lags associated with NaN 
%     coefficients found in LagOpMA are estimated; all others are held fixed.
%
%   constant - Estimate of the constant included in a general ARMA model. 
%
%   variance - Estimate of the unconditional variance of the noise process.
%
% Notes:
%
%   o The algorithm, outlined in Appendix A6.2 of [1], is a 2-step process 
%     in which the AR coefficients are first found by solving the modified 
%     Yule-Walker equations; the MA coefficients are then found by an 
%     iterative technique applied to the filtered time series derived from
%     the AR coefficients estimated in step one.
%
%   o If the input autoregressive lag operator polynomial, LagOpAR, is
%     associated with an ARIMA(P,D,Q) model, then all non-seasonal and 
%     seasonal integration must be excluded from the input polynomial
%     LagOpAR as well as the input series. In this case, the input time 
%     series must be filtered to remove the effects of integration, and
%     LagOpAR represents the product of stable non-seasonal and seasonal 
%     autoregressive polynomials. Similarly, LagOpMA represents the product 
%     of invertible non-seasonal and seasonal moving average polynomials.
%
% References:
%
%   [1] Box, G. E. P., G. M. Jenkins, and G. C. Reinsel. Time Series
%       Analysis: Forecasting and Control. 3rd edition. Upper Saddle River,
%       NJ: Prentice-Hall, 1994.
%
%   [2] Hamilton, J. D. Time Series Analysis. Princeton, NJ: Princeton
%       University Press, 1994.

% Copyright 1999-2011 The MathWorks, Inc.

%
% Get the degree of each polynomial.
%

P = LagOpAR.Degree;
Q = LagOpMA.Degree;

%
% Check for the case of no ARMA model. In this case, compute the sample
% mean and variance and exit.
%

if (P + Q) == 0
   AR       = [];
   MA       = [];
   constant = mean(y);
   variance = var(y,1);
   return
end

%
% Set the polynomial stability tolerance to be consistent with the stability 
% tolerance of the LagOp/isStable method.
%

tolerance = 10 * eps;

%
% Determine if any AR or MA coefficients are held fixed as user-specified 
% initial guesses or exclusion constraints.
%

AR0      = toCellArray(reflect(LagOpAR));
AR0      = cell2mat(AR0(2:end))';          % Exclude coefficient at lag zero
iFixedAR = ~isnan(AR0);                    % Identify lags of coefficients held fixed

MA0      = toCellArray(LagOpMA);
MA0      = cell2mat(MA0(2:end))';          % Exclude coefficient at lag zero
iFixedMA = ~isnan(MA0);                    % Identify lags of coefficients held fixed

%
% Estimate the AR coefficients of a general ARMA(P,Q) model:
%

if P > 0

%
% Compute the autocovariance sequence of the y(t) process. The variance
% (autocovariance at lag 0) is found in the first element.
%

   correlation = autocorr(y,P + Q);      % autocorrelation sequence
   covariance  = correlation * var(y,1); % autocovariance sequence

   if any(~isfinite(covariance))
%
%     The presence of INF or NaN is found, so just exit gracefully.
%
      AR           = zeros(P,1);
      MA           = zeros(Q,1);
      AR(iFixedAR) = AR0(iFixedAR);    % Restore coefficients held fixed
      MA(iFixedMA) = MA0(iFixedMA);    % Restore coefficients held fixed
      constant     = mean(y);
      variance     = var(y,1);
      return
   end

%
%  In each case below, the matrix C of covariances is given in equation 
%  A6.2.1 (page 220) of [1].
%

   if Q > 0

%
%     For ARMA processes, the matrix C of covariances derived from the 
%     estimated autocovariance sequence is Toeplitz, but nonsymmetric. The
%     AR coefficients are found by solving the modified Yule-Walker equations.
%
      i = abs(Q:-1:Q-P+1) + 1;        % covariance(k) = covariance(-k)
      C = toeplitz(covariance(Q+1:Q+P), covariance(i));
      d = covariance(Q+2:Q+P+1);

   else

%
%     For AR processes, the matrix C of covariances derived from the estimated 
%     auto-covariance sequence is Toeplitz and symmetric. The AR coefficients 
%     are found by solving the Yule-Walker equations.
%
      C = toeplitz(covariance(1:P));
      d = covariance(2:P+1);

   end

%
%  Estimate AR coefficients by solving a least-squares problem. 
%
   
   if sum(iFixedAR) == 0
%
%     No coefficients are held fixed, so solve the least-squares problem with 
%     a simple backslash operation.
%
      AR = C \ d;
      
   else
%
%     Some coefficients are held fixed, so impose equality constraints by 
%     solving a constrained linear least-squares problem with LSQLIN.
%
      Aeq               = zeros(sum(iFixedAR), P);
      Aeq(:,iFixedAR)   = eye(sum(iFixedAR));
      beq               = AR0(iFixedAR);
      options           = optimoptions('lsqlin', 'Algorithm', 'interior-point', 'Display', 'off'); 
      [AR,~,~,exitFlag] = lsqlin(C, d, [], [], Aeq, beq, [], [], [], options);
%
%     As a safeguard, revert to the traditional backslash solution if necessary.
%
      if exitFlag <= 0
         AR           = C \ d;
         AR(iFixedAR) = AR0(iFixedAR);  % Restore coefficients held fixed
      end
      
   end
%
%  Test the AR process for stationarity. 
%
%  In the event the AR polynomial is unstable, set all estimated coefficients 
%  to 0 while retaining those held fixed.
%
   eigenValues = roots([1 ; -AR]);

   if any(abs(eigenValues) >= (1 - tolerance))
      AR            = zeros(P,1);
      AR(iFixedAR)  = AR0(iFixedAR);  % Restore coefficients held fixed
   end

else

   AR = [];

end

%
% Filter the ARMA(P,Q) input series y(t) with the estimated AR coefficients 
% to obtain a pure MA process. If the input moving-average model order Q is
% zero, then the filtered output is just a pure innovations process (i.e., 
% an MA(0) process); in this case the innovations variance estimate is just 
% the sample variance of the filtered output. If Q > 0, then compute the 
% autocovariance sequence of the MA process and continue.
%

x        = filter([1 -AR'], 1, y);
constant = mean(x);

if Q == 0
   variance = var(x,1);
   MA       = [];
   return
end

c = autocorr(x,Q) * var(x,1); % Covariance of an MA(Q) process

%
% Estimate the variance of the white noise innovations process e(t) and the 
% MA coefficients of a general ARMA(P,Q) model. The method of computation 
% is outlined in equation A6.2.4 (page 221) of [1].
%

MA           = zeros(Q,1);             % Initialize MA coefficients
MA(iFixedMA) = MA0(iFixedMA);          % Copy coefficients held fixed
MA1          = ones(Q,1);              % Saved MA coefficients from previous iteration
counter      = 1;                      % Iteration counter
tol          = 0.01;                   % Convergence tolerance

while ((norm(MA - MA1) > tol) && (counter < 100))

    MA1 = MA;
%
%   Estimate the variance of the innovations process e(t):
%
    variance = c(1)/([1 ; MA]'*[1 ; MA]);
%
%   Ensure invertibility of the MA polynomial. If necessary, set all 
%   estimated coefficients to 0 while retaining those held fixed.
%
    eigenValues = roots([1 ; MA]);
    
    if (variance <= 0) || any(abs(eigenValues) >= (1 - tolerance))
       MA           = zeros(Q,1);
       MA(iFixedMA) = MA0(iFixedMA);
       variance     = var(x,1);
       return
    end

    for j = Q:-1:1
%
%       Estimate the moving-average coefficients. The MA coefficients are
%       the negative of those appearing in equation A6.2.4 (page 221) of [1]. 
%
%       The following test enforces user-specified initial guesses and exclusion 
%       constraints associated with zero-valued coefficients (i.e., coefficients 
%       associated with lags missing from the MA polynomial).
%
        if isnan(MA0(j))    
           MA(j) = [c(j+1) ; -MA(1:Q-j)]' * [1/variance ; MA(j+1:Q)];
        end

    end
    
    counter = counter + 1;    % Update the counter

end

%
% Once again, ensure invertibility of the MA polynomial in the event the
% above WHILE loop terminated due to excessive iterations rather than
% convergence.
%

eigenValues = roots([1 ; MA]);
    
if any(abs(eigenValues) >= (1 - tolerance))
   MA           = zeros(Q,1);         % Set estimated coefficients to zero
   MA(iFixedMA) = MA0(iFixedMA);      % Restore coefficients held fixed
   variance     = var(x,1);
end

end