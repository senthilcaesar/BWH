function Z = simulateStandardizedVariates(OBJ, numRows, numColumns)
%SIMULATESTANDARDIZEDVARIATES Generate i.i.d. standardized residuals.
%
% Syntax:
%
%   Z = simulateStandardizedVariates(OBJ,numRows,numColumns)
%
% Description:
%
%   Given the required number of rows (observations) and columns (paths),
%   generate an array of mean-zero, unit-variance, i.i.d. random variates 
%   from the probability distribution found in the input model OBJ. 
%
% Input Arguments:
%
%   OBJ - Model specification object, encapsulating the probability 
%     distribution and any necessary parameters required to simulate the
%     i.i.d. observations.
%
%   numRows - The required number of observations (rows).
%
%   numColumns - The required number of paths (columns).
%
% Output Arguments:
%
%   Z - A numRows-by-numColumns array of i.i.d. random variates consistent 
%     with the probability distribution found in the input specification 
%     OBJ. 

% Copyright 2018 The MathWorks, Inc.

%
% Determine the distribution and generate standardized noise.
%

switch upper(OBJ.Distribution.Name)
    
   case 'GAUSSIAN'
       
      Z = randn(numRows,numColumns); % i.i.d. N(0,1) variates

   case 'T'
       
      DoF = OBJ.Distribution.DoF;
      if isnan(DoF)
         error(message('econ:LagIndexableTimeSeries:simulateStandardizedVariates:UnspecifiedDoF'))
      end

      Z = randn(numRows,numColumns);

      if DoF <= 200   % i.i.d. standardized t(0,1) variates
         Z = Z .* sqrt(DoF./(randg(DoF/2,numRows,numColumns).*2)) / sqrt(DoF/(DoF-2));
      end

   otherwise
      error(message('econ:LagIndexableTimeSeries:simulateStandardizedVariates:InvalidDistribution'))
      
end

end % Simulate standardized residuals.