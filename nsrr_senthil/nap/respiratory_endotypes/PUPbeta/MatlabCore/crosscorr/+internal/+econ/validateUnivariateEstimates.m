function value = validateUnivariateEstimates(name, value, numElements, lagsToKeep)
%VALIDATEUNIVARIATEESTIMATES Validate initial estimates for univariate models
%
% Syntax:
%
%   value = validateUnivariateEstimates(name,value,numElements,lagsToKeep)
%
% Description:
%
%   Ensure that a scalar coefficient or vector of coefficients is of correct 
%   length. Additionally, once the initial validation is complete, any
%   vector-valued input is re-sized as necessary to enforce any zero-valued
%   exclusion constraints. If validation fails, an error is thrown.
%
% Input Arguments:
%
%   name - Character string indicating the name of the input coefficient or
%     parameter to validate.
%
%   value - Scalar or vector coefficient or parameter to validate.
%
%   numElements - A positive integer indicating the number of elements
%     (i.e., length) that value is expected to contain.
%
%   lagsToKeep - A vector of positive integers indicating the lags of value
%     retained in the output. In the event that value is a vector of 
%     coefficients, the validated output value will contain the lags in 
%     lagsToKeep of the input value.
%
% Output Arguments:
%
%   value - Scalar or vector of validated coefficients or parameters. If no
%     error occurs, the output value will always have the same number of 
%     elements as the number of lags found in lagsToKeep.

% Copyright 1999-2011 The MathWorks, Inc.

%
% Ensure the input coefficient is a vector of correct length.
%

if (numel(value) ~= numElements) || (~isvector(value) && ~isempty(value))
   error(message('econ:internal:econ:validateUnivariateEstimates:IncorrectLength', name, numElements))
end

%
% Now ensure only the relevant lags are retained. The following check
% associates a NaN with the input lags as "not applicable", which means  
% the coefficient is a constant having no lags associated with it.
%

if ~isnan(lagsToKeep)
   value = value(lagsToKeep);
end

end
