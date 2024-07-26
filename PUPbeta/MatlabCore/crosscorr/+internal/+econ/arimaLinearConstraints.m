function constraints = arimaLinearConstraints(LagsAR,LagsSAR,LagsMA,LagsSMA,beta,...
                       isDistributionT,isVarianceConstant,solve,x0,tolerance)
%ARIMALINEARCONSTRAINTS Linear parameter constraints for ARIMA models
%
% Syntax:
%
%   constraints = arimaLinearConstraints(LagsAR,LagsSAR,LagsMA,LagsSMA,...
%              beta,isDistributionT,isVarianceConstant,solve,x0,tolerance)
%
% Description:
%
%   Given information about an ARIMA model, create a data structure with 
%   lower and upper bounds and linear equality and inequality constraints. 
%
% Input Arguments:
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
%   beta - Vector of regression coefficients.
%
%   isDistributionT - Boolean flag indicating whether or not the residual 
%     process is fit to a standardized t distribution. TRUE indicates the 
%     residuals are t-distributed; FALSE indicates that the residuals are 
%     not t-distributed.
%
%   isVarianceConstant - Boolean flag indicating whether or not the variance
%     is constant. TRUE indicates the variance is constant; FALSE indicates
%     the variance ARIMA model contains a conditional variance model.
%
%   solve - Logical column vector of flags that indicate which parameters 
%     are estimated and which are held fixed throughout the optimization 
%     process. TRUE in the i-th element indicates that the i-th parameter 
%     is estimated; FALSE indicates that the corresponding parameter is 
%     held fixed at some user-specified initial value. 
%
%   x0 - Column vector of initial parameter values included in the 
%     optimization. This vector is used to construct the output equality 
%     constraints Aeq and beq. Elements of x0 correspond to model 
%     coefficients ordered as follows:
%
%      o Constant
%      o non-zero  AR coefficients
%      o non-zero SAR coefficients 
%      o non-zero  MA coefficients 
%      o non-zero SMA coefficients
%      o variance (constant-variance models only)
%      o DoF (t distribution only)
%
%   tolerance - Positive numerical tolerance, or offset, applied to all 
%     constraints. Tolerance parameter specifies how close estimated 
%     coefficients may get to a corresponding constraint, and is related to 
%     the maximum constraint violation tolerance parameter (TolCon) of the 
%     Optimization Toolbox.
%
% Output Arguments:
%
%   constraints - A data structure of constraint information with the
%     following fields:
%
%     lb - Vector of coefficient lower bounds such that lb <= x.
%
%     ub - Vector of coefficient upper bounds such that x <= ub.
%
%     A - Linear inequality constraint matrix of the form A*x <= b. Each 
%       row corresponds to a linear inequality constraint, and each column
%       corresponds to an estimated coefficient.
%
%     b - Linear inequality constraint vector of the form A*x <= b. Each
%       element corresponds to a linear inequality constraint.
%
%     Aeq - Linear equality constraint matrix of the form Aeq*x = beq. Each
%       row corresponds to a linear equality constraint, and each column
%       corresponds to an estimated coefficient.
%
%     beq - Linear equality constraint vector of the form Aeq*x = beq.
%       Each element corresponds to a linear inequality constraint.
%
% Note:
%
%   o All elements of constraint vectors and columns of constraint matrices
%     correspond to the order of the model coefficients found in x0.
%

% Copyright 2013 The MathWorks, Inc.

%
% Preallocate the lower & upper bound constraint vectors.
%

nAR   = numel(LagsAR) + numel(LagsSAR);
nMA   = numel(LagsMA) + numel(LagsSMA);
nBeta = numel(beta);

nParameters = 1 + nAR + nMA + nBeta + isVarianceConstant + isDistributionT;
LB          = zeros(nParameters,1);
UB          = zeros(nParameters,1);
hundred     = 100;

isRegressionIncluded = (nBeta > 0);

%
% Set the following lower and upper bounds for ARIMA models.
%
%   Coefficient  Lower Bound  Upper Bound
%   -----------  -----------  -----------
%    Constant        -10          10
%     AR/SAR         see notes below
%     MA/SMA         see notes below
%      Beta         -100         100
%    Variance         0          100
%

i1 = 1;                                            % First index of ARIMA parameters (Constant/Intercept)
i2 = i1 + nAR + nMA + nBeta + isVarianceConstant;  % Last  index of ARIMA parameters (Variance for constant-variance models only!)

LB(i1:i2) = [-10 ; repmat(-1 + tolerance, nAR + nMA, 1) ; repmat(-hundred, nBeta, 1) ; tolerance(isVarianceConstant)];
UB(i1:i2) = [ 10 ; repmat( 1 - tolerance, nAR + nMA, 1) ; repmat( hundred, nBeta, 1) ; hundred(isVarianceConstant)];

%
% Adjust the lower & upper bounds for AR, SAR, MA, and SMA coefficients
% using an ad-hoc approach:
%
%   o If the polynomial degree is greater than 1 AND the first lag is included 
%     in the estimation, then set the bounds of the first lag to +/-2.
%
%   o If the polynomial degree is greater than 2 AND the second lag is included 
%     in the estimation, then set the bounds of the second lag to +/-1.5.
% 
%   o The coefficient bounds at all other lags are +/-1, which are already
%     set above.
%

iCount =  2;

limit  =  2.0;
lower  = -1.5;
upper  =  1.5;

if numel(LagsAR) > 0
   degree  = LagsAR(end);
   for L = LagsAR
       if (L == 1) && (degree > 1)
          LB(iCount) = -limit + tolerance;
          UB(iCount) =  limit - tolerance;
       elseif (L == 2) && (degree > 2)
          LB(iCount) = lower;
          UB(iCount) = upper;
       end
       iCount = iCount + 1;
   end
end

if numel(LagsSAR) > 0
   degree = LagsSAR(end);
   for L = LagsSAR
       if (L == 1) && (degree > 1)
          LB(iCount) = -limit + tolerance;
          UB(iCount) =  limit - tolerance;
       elseif (L == 2) && (degree > 2)
          LB(iCount) = lower;
          UB(iCount) = upper;
       end
       iCount = iCount + 1;
   end
end

if numel(LagsMA) > 0
   degree  = LagsMA(end);
   for L = LagsMA
       if (L == 1) && (degree > 1)
          LB(iCount) = -limit + tolerance;
          UB(iCount) =  limit - tolerance;
       elseif (L == 2) && (degree > 2)
          LB(iCount) = lower;
          UB(iCount) = upper;
       end
       iCount = iCount + 1;
   end
end

if numel(LagsSMA) > 0
   degree = LagsSMA(end);
   for L = LagsSMA
       if (L == 1) && (degree > 1)
          LB(iCount) = -limit + tolerance;
          UB(iCount) =  limit - tolerance;
       elseif (L == 2) && (degree > 2)
          LB(iCount) = lower;
          UB(iCount) = upper;
       end
       iCount = iCount + 1;
   end
end

%
% If a t distribution is selected, then there is the additional constraint
%
%   DoF > 2
%
% so bound 2 < DoF <= 200.
%

if isDistributionT
   LB(i2 + 1) = 2 + tolerance;   % Strictly greater than 2
   UB(i2 + 1) = 200;             % Limit beyond which a t distribution is Gaussian
end

%
% In the event that an element of the initial estimate vector x0 falls 
% outside the corresponding lower/upper bound, then reset the bound to the 
% initial estimate. 
%

LB = min(LB,x0);
UB = max(UB,x0);

%
% Further widen the lower and upper bounds on the constant/intercept and the
% regression coefficients, and the upper bound of the variance (constant-variance 
% models only) by an arbitrary 20 percent.
%

LB(1) = min(LB(1), 1.2 * x0(1));
UB(1) = max(UB(1), 1.2 * x0(1));

if isRegressionIncluded
   iBeta     = (1 + nAR + nMA + 1):(1 + nAR + nMA + nBeta);
   LB(iBeta) = min(LB(iBeta), 1.2 * x0(iBeta));
   UB(iBeta) = max(UB(iBeta), 1.2 * x0(iBeta));
end

if isVarianceConstant
   UB(i2) = max(UB(i2), 1.2 * x0(i2));
end

%
% Set user-specified linear equality constraints of the form Aeq*x = beq.
%

Fix = ~solve;

if any(Fix)

   i   = find(Fix); 
   Aeq = zeros(length(i),nParameters);

   for j = 1:length(i)
       Aeq(j,i(j)) = 1;
   end

   beq = x0(Fix);

else

   Aeq = [];
   beq = [];

end

%
% Pack all constraints into a data structure.
%

constraints.lb  =  LB;
constraints.ub  =  UB;
constraints.A   =  [];
constraints.b   =  [];
constraints.Aeq = Aeq;
constraints.beq = beq;

end