classdef LagIndexableTimeSeries
%LAGINDEXABLETIMESERIES Create a lag-indexable time series object
%
% Description:
%
%   LagIndexableTimeSeries is the abstract superclass for all univariate 
%   time series classes whose coefficients associated with lag operator 
%   polynomials are indexable by lag. The class provides subscripted 
%   referencing and assignment methods, publishes method and property 
%   templates, and serves as a repository for some common static methods.
%

% Copyright 1999-2011 The MathWorks, Inc.

properties (Access = public)
  Distribution = struct('Name', 'Gaussian')
end

properties (GetAccess = public, SetAccess = private, Abstract)
  P
  Q
end

properties (Abstract)
  Constant
end

properties (Access = protected, Abstract)
  LHS       % Cell array of properties of subclasses whose polynomials are not reflected
  RHS       % Cell array of properties of subclasses whose polynomials are reflected
end

methods (Access = public, Abstract)
  varargout  = simulate(OBJ, numObs, varargin)
  varargout  = estimate(OBJ, Y, varargin)
  varargout  = forecast(OBJ,Y0, varargin)
  varargout  = infer   (OBJ, Y, varargin)
  polynomial = getLagOp(OBJ, name)
  OBJ        = setLagOp(OBJ, name, polynomial)
  OBJ        = validateModel(OBJ, varargin)
end

methods

%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

  function OBJ = set.Distribution(OBJ, Distribution)

%
%  Impose consistency on the probability distribution.
%

   if ischar(Distribution) || ( isstring(Distribution) && isscalar(Distribution) )

      if strcmpi(Distribution, 'Gaussian')
         Distribution = struct('Name', 'Gaussian');
      elseif strcmpi(Distribution, 't')
         Distribution = struct('Name', 't', 'DoF', NaN);
      else
         error(message('econ:LagIndexableTimeSeries:Distribution:InvalidDistributionString'))
      end

   elseif isstruct(Distribution)

      fields = fieldnames(Distribution);

      if ~any(strcmp('Name', fields))
         error(message('econ:LagIndexableTimeSeries:Distribution:MissingName'))
      end

      if strcmpi(Distribution.Name, 'Gaussian')

         if numel(fields) > 1
            error(message('econ:LagIndexableTimeSeries:Distribution:InvalidGaussian'))
         end
         Distribution.Name = 'Gaussian';

      elseif strcmpi(Distribution.Name, 't')

         if ~any(strcmp('DoF', fields)) || (numel(fields) > 2)
            error(message('econ:LagIndexableTimeSeries:Distribution:InvalidT'))
         end

         if ~isscalar(Distribution.DoF) || ~strcmp(class(Distribution.DoF), 'double') ...
                                        || (Distribution.DoF <= 2)
            error(message('econ:LagIndexableTimeSeries:Distribution:InvalidDoF'))
         end

         Distribution.Name = 't';

      else
         error(message('econ:LagIndexableTimeSeries:Distribution:InvalidDistributionName'))
      end

   end

   OBJ.Distribution = Distribution;

  end

end

methods (Access = public)

%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

function varargout = subsref(OBJ, args)
%SUBSREF Subscripted referencing of LagIndexableTimeSeries objects

switch args(1).subs

  case [OBJ.LHS  OBJ.RHS]   % Is the reference to a coefficient indexed by lag?

    if numel(args) > 1
%
%      When indexing into a property associated with a Lag Operator
%      Polynomial, the indices must be either a double vector of
%      non-negative integer lags or a colon (':').
%
       if ~(isempty(args(2).subs) || isa(args(2).subs{:}, 'double') || strcmp(args(2).subs,':'))
          error(message('econ:LagIndexableTimeSeries:subsref:InvalidReference'))
       end
%
%      Get the underlying polynomial.
%
       polynomial = getLagOp(OBJ, args(1).subs);

       if any(strcmp(args(1).subs, OBJ.LHS))  % Is it a LHS polynomial?
          isLagOpLHS = true;
       else
          isLagOpLHS = false;
       end

       if strcmpi(args(2).type, '{}')

          if strcmp(args(2).subs,':')
             args(2).subs{1} = (1:polynomial.Degree);  % Assume lags 1, 2, ...
          end
%
%         Let the contained polynomial handle the subscripted referencing 
%         to ensure consistency.
%
          args(1).subs = 'Coefficients'; 
          [varargout{1:nargout}] = subsref(polynomial, args);

          if ~isempty(varargout)
             if isLagOpLHS
%
%               Since {} cell indexing has been requested, negate all coefficients
%               associated with non-zero lags to allow for the difference equation
%               sign convention.
%
                for i = 1:max(nargout,1)
                    if args(2).subs{:}(i) ~= 0
                       varargout{i} = -varargout{i};
                    end
                end
             end
          end

       elseif strcmpi(args(2).type, '()')

          if nargout > 1
             error(message('econ:LagIndexableTimeSeries:subsref:InvalidParenthesisReference'))
          end

          if isempty(args(2).subs)           % Handle "()" references.
             error(message('econ:LagIndexableTimeSeries:subsref:NoLagIndices'))
          end

          C = LagOp(subsref(polynomial.Coefficients, args(2:end)));
          
          if ischar(args(2).subs{:})         % Handle "(:)" references.
             args(2).subs{1} = 1:C.Degree;   % Get all lags >= 1
             N     = numel(args(2).subs{:});
             array = cell(N, N > 0);         % Ensure "(:)" returns a column array.
          else
             N     = numel(args(2).subs{:});
             array = cell(N > 0 ,N);         % Ensure a row array.
          end
%
%         Negate coefficients associated with non-zero lags.
%
          if isLagOpLHS
             C = reflect(C);
          end
%
%         Extract the actual coefficients.
%
          C = C.Coefficients;

          for i = 1:N
              array{i} = C{args(2).subs{1}(i)};
          end
          varargout = {array};
       end

    else

        varargout = {OBJ.(args.subs)};
    end

  otherwise  

    if ismethod(OBJ, args(1).subs)                      % Is it a method?
%
%      The reference is to a method, so just call the built-in.
%
       [varargout{1:nargout}] = builtin('subsref', OBJ, args);

    elseif any(strcmp(args(1).subs, properties(OBJ)))   % Is it a property?
%
%      The reference is to a property.
%
       if isa(OBJ.(args(1).subs), 'internal.econ.LagIndexableTimeSeries')
%
%         The syntax references a property which is itself an LagIndexableTimeSeries object. 
%
          if numel(args) == 1
%
%            The contained object is referenced: OBJ.Variance
%
             varargout = {builtin('subsref', OBJ, args)};
          else
             if ischar(args(2).subs)
%
%               The reference is to a method or property of contained 
%               LagIndexableTimeSeries object, with no lags specified.
%
                if ismethod(OBJ.(args(1).subs), args(2).subs)
%
%                  A method of a contained LagIndexableTimeSeries object is 
%                  referenced, such as:
%
%                  [...] = OBJ.Variance.simulate(...)
%
                   [varargout{1:nargout}] = subsref(OBJ.(args(1).subs), args(2:end));
                else
                   if numel(args) == 2
%
%                     The syntax references a property of a contained LagIndexableTimeSeries 
%                     object, with no lags specified. For example, consider a 
%                     GARCH model stored in the Variance property of an ARIMA model. 
%
%                     The following handles the syntaxes:
%
%                     OBJ.Variance.Garch
%                     OBJ.Variance.Arch
%
                      varargout = {builtin('subsref', OBJ.(args(1).subs), args(2:end))};
                   else
%
%                     The syntax references a property of a contained LagIndexableTimeSeries 
%                     object, with lags specified. For example, consider a GARCH 
%                     model stored in the Variance property of an ARIMA model. 
%
%                     The following handles the syntaxes:
%
%                     OBJ.Variance.Garch([1 2])
%                     OBJ.Variance.Arch{2}
%
                      [varargout{1:nargout}] = subsref(OBJ.(args(1).subs), args(2:end));
                   end
               end
            else
%
%              The reference is to a property inaccessible via lag-based indexing.
%              However, lags were specified and so the command makes no sense. 
%
               error(message('econ:LagIndexableTimeSeries:subsref:InvalidLagReference'))
            end
         end
       else
%
%        The reference is to a property that is not an LagIndexableTimeSeries object.
%
         if ischar(args(end).subs)
%
%           The last subscript is a character string, and so might be a
%           valid reference such as:
%
%           OBJ.Distribution
%           OBJ.Distribution.DoF
%           OBJ.Variance.Distribution
%           OBJ.Variance.Distribution.DoF
%
            varargout = {builtin('subsref', OBJ, args)};
         else
%
%           The last subscript is not a string, and so is of the form:
%
%           OBJ.Distribution(...)
%           OBJ.Variance.Distribution.DoF{...}
%
%           and references a property inaccessible via lag-based indexing. 
%           However, lags were specified and so the command makes no sense. 
%
            error(message('econ:LagIndexableTimeSeries:subsref:InvalidLagReference'))
         end
       end

     else

       error(message('econ:LagIndexableTimeSeries:subsref:InvalidReference'))

     end

end

end  % SUBSREF

%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

function OBJ = subsasgn(OBJ,args,value)
%SUBSASGN Subscripted assignment of LagIndexableTimeSeries objects

try

   switch args(1).subs

      case [OBJ.LHS OBJ.RHS]   % Lag Operator Polynomials

%
%        When assigning to properties associated with Lag Operator Polynomials,
%        the right-hand side value must be either a cell vector or a double
%        matrix (direct assignment of LagOp objects is not allowed, which is
%        consistent with the corresponding constructors).
%
         isCell = iscell(value);

         if ~(isa(value,'double') || isCell)
            error(message('econ:LagIndexableTimeSeries:subsasgn:InvalidAssignmentType'))
         end
%
%        Get the underlying polynomial to update.
%
         polynomial = getLagOp(OBJ, args(1).subs);
%
%        Check that the assigned coefficients (i.e., value) are specified as 
%        a 1-by-N or n-by-1 cell vector.
%
         if isCell
%
%           The following check excludes assignments such as:
%
%           OBJ.LHS/RHS = {0.2 -0.4 ; 0.1 0.07}
%
%           which seems somewhat ambiguous for polynomials of any dimensionality.
%
            if ~isempty(value) && (size(value,1) ~= 1) && (size(value,2) ~= 1)
                error(message('econ:LagIndexableTimeSeries:subsasgn:InvalidCoefficientCellFormat'))
            end
         end
% 
%        Reflect coefficients associated with LHS polynomials.
%
         if any(strcmp(args(1).subs, OBJ.LHS))     % Is it a LHS polynomial?
            for i = 1:numel(value)
                if (numel(args) == 1) || ((i <= numel(args(2).subs{1})) && (args(2).subs{1}(i) ~= 0))
                   if isCell
                      value{i} = -value{i};
                   else
                      value(i) = -value(i);
                   end
                end
            end
         end
%
%        Since model coefficients are re-directed to the underlying OBJ.LHS/RHS 
%        lag operator polynomial, save the name of the original coefficient 
%        and overwrite the first subscript reference to access the appropriate 
%        coefficients of the polynomial.
%
         coefficient  = args(1).subs;   
         args(1).subs = 'Coefficients';
%
%        Handle the following syntaxes:
%
%        OBJ.LHS = {G1, G2, ... GP} or OBJ.RHS  = {A1, A2, ... AQ}
%
%        in which no lag references appear on the left-hand side of the 
%        assignment and cell vectors appear on the right-hand side. The 
%        coefficients are assumed to correspond to lags 1, 2, ... to the 
%        degree of the polynomial.
%
%        This syntax is consistent with the default behavior of the constructor
%        and supports the common "reduced-form" representation in which the 
%        polynomial coefficient at lag 0 is assumed known (usually 1 or 0), and
%        is saved.
%
         if (numel(args) == 1)

            polynomial.Coefficients = polynomial.Coefficients(0);
            args(2).type            = '()'; 

            if isempty(value)
%
%              The right-hand side cell array is empty, so overwrite it with 
%              a cell array whose only element is a matrix of zeros. This 
%              allows the subsequent conversion to a Lag Operator Polynomial
%              to determine the dimensionality of the underlying polynomial.
%
%              The following line of code converts the assignment from
%              OBJ.LHS/RHS = {} to OBJ.LHS/RHS = {0}.
%
               value = {zeros(polynomial.Dimension)};
            end

            args(2).subs = {1:numel(value)}; % Assume lags 1, 2, ...

         end
%
%        Update the underlying polynomial and then update the object.
%
         polynomial = subsasgn(polynomial, args, value);
         OBJ        = setLagOp(OBJ, coefficient, polynomial);

      otherwise

         if ismethod(OBJ, args(1).subs)
            error(message('econ:LagIndexableTimeSeries:subsasgn:InvalidMethodAssignment'))
         end

         if isa(OBJ.(args(1).subs), 'internal.econ.LagIndexableTimeSeries') && (numel(args) > 1)
%
%           The value is assigned to a property which is itself a 
%           LagIndexableTimeSeries object. 
%
%           For example, consider the situation in which a GARCH model is 
%           stored in the Variance property of an ARIMA model. 
%
%           The following would handle assignments like the following:
%
%           OBJ.Variance.Garch([1 2]) = {0.4 0.2}
%           OBJ.Variance.Arch{2}      = 0.03
%
%           OBJ.Variance.Distribution = 'Gaussian'
%           OBJ.Variance.Distribution = struct('Name', 't', 'DoF', pi)
%
            OBJ.(args(1).subs) = subsasgn(OBJ.(args(1).subs), args(2:end), value);
%
%           The 'Distribution' property is handled as a special case to
%           enforce consistency between a contained LagIndexableTimeSeries 
%           object and its container, also an LagIndexableTimeSeries object.
%
            if strcmpi('Distribution', args(2).subs)
               OBJ = builtin('subsasgn', OBJ, args(2:end), value);
            end

         else

            if numel(args) > 1
               if strcmpi('Distribution', args(1).subs) && ischar(args(end).subs)
%
%                 For clarity, the 'Distribution' property is handled as a 
%                 special case. The following check disallows direct assignments
%                 to specific fields of the underlying distribution
%                 property, such as 
%
%                 OBJ.Distribution.Name = 't'
%                 OBJ.Distribution.DoF  =  10

                  error(message('econ:LagIndexableTimeSeries:subsasgn:InvalidDistributionAssignment'))
               else
                  error(message('econ:LagIndexableTimeSeries:subsasgn:InvalidLagAssignment'))
               end
            else
               OBJ = builtin('subsasgn', OBJ, args, value);

               if strcmpi('Distribution', args(1).subs)
%
%                 The 'Distribution' property is handled as a special case
%                 to enforce consistency between a contained LagIndexableTimeSeries object 
%                 and its container, also an LagIndexableTimeSeries object.
%
                  OBJ = validateModel(OBJ, 'Distribution',  OBJ.Distribution);
               else
                  OBJ = validateModel(OBJ);
               end
            end

         end
   end

catch exception

   exception.throwAsCaller();

end

end  % SUBSASGN

end  % METHODS Block


methods (Hidden)
 
%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

function varargout = cat(OBJ, varargin)
  error(message('econ:LagIndexableTimeSeries:cat:CannotConcatenate'))
end

%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

function varargout = horzcat(OBJ, varargin)
  error(message('econ:LagIndexableTimeSeries:horzcat:CannotConcatenate'))
end

%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

function varargout = vertcat(OBJ, varargin)
  error(message('econ:LagIndexableTimeSeries:vertcat:CannotConcatenate'))
end

%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

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


end % Methods (Hidden)


methods (Static, Hidden) % Error Checking Utilities

%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

function empty
  error(message('econ:LagIndexableTimeSeries:empty:InvalidMethod'))
end

%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

function outputArray = checkPresampleData(outputArray, inputName, inputArray, nRows)
%CHECKPRESAMPLEDATA Check and allocate user-specified pre-sample data
%
% Syntax:
%
%	Y = OBJ.checkPresampleData(Y, name, Y0, numRows)
%
% Description:
%
%   Given the required number of rows (observations) and columns (paths),
%   check user-specified pre-sample data arrays. If no errors are found, 
%   then update the output pre-sample data array with the relevant 
%   pre-sample data.
%
% Inputs:
%   Y - A matrix of presample data pre-allocated to the correct size. This 
%     matrix will have a sufficient number of rows to store all required 
%     presample observations and a number of columns equal to the number of 
%     sample paths.
%
%   name - A character string containing the name of the presample
%     data array being checked for errors. This string is placed into the
%     output error message for additional information if necessary.
%
%   Y0 - The user-specified presample data matrix to check. This must be a
%     column vector or a matrix with at least as many columns as Y. In 
%     either case, it must have at least numRows rows.
%
%   numRows - The required number of presample observations (i.e., rows).
%
% Outputs:
%   Y - If no errors occur, this will be an updated version of the input 
%     matrix Y, in which case the last numRows rows will contain the
%     required pre-sample data.
%

%
% Test for errors in the pre-sample data.
%

if size(inputArray,1) < nRows

   error(message('econ:LagIndexableTimeSeries:checkPresampleData:InvalidPreSampleRows', inputName, nRows))

end

nColumns  =  size(outputArray,2);

if nRows > 0     % Allow an empty matrix

   if (size(inputArray,2) ~= 1) && (size(inputArray,2) < nColumns)
      if nColumns == 1
         error(message('econ:LagIndexableTimeSeries:checkPresampleData:InvalidPreSampleColumnVector', inputName))
      else
         error(message('econ:LagIndexableTimeSeries:checkPresampleData:InvalidPreSampleMatrix', inputName, nColumns))
      end
   end

%
%  Re-format the error-free pre-sample data.
%

   inputArray  =  inputArray((end - nRows + 1):end , :);  % Retain no more rows than necessary.

%
%  If the user-specified input presample array is a single 
%  column vector, replicate it across all columns (i.e., paths).
%

   if size(inputArray,2) == 1
      outputArray((end - nRows + 1):end , :)  =  repmat(inputArray , 1 , nColumns);
   else
      outputArray((end - nRows + 1):end , :)  =  inputArray(: , 1:nColumns);
   end        

end

end

%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

function varargout = listwiseDelete(varargin)
%LISTWISEDELETE Remove missing observations by listwise deletion.
%
% Syntax:
%
%   [series1,series2,...] = listwiseDelete(OBJ,series1,series2,...)
%
% Description:
%
%   This static method synchronizes and removes missing observations from 
%   any number of time series. Missing observations are indicated by NaN's. 
%   When removing missing data, listwise deletion is applied across all 
%   time series as if the inputs are, taken collectively, a single 
%   composite series.
%
% Input Arguments:
%
%   series1, series2, ... - A variable-length list of time series arrays
%     in which rows correspond to time-tagged observations and columns to 
%     paths or variables. The last row of each series contains the most 
%     recent observation. Missing observations are indicated by NaN's. The 
%     inputs do not need to have the same number of rows nor the same 
%     number of columns. However, when removing missing data, all inputs 
%     are merged into a single, composite series, and listwise deletion is 
%     applied across each row of the composite. The data is synchronized 
%     such that the last (most recent) observation of each series is 
%     assumed to occur at the same time.
%
% Output Arguments:
%
%   series1, series2, ... - A list of time series arrays with missing
%     observations removed by listwise deletion. The number of outputs is
%     indicated by the user, but the outputs are always ordered such that
%     the first output series corresponds to the first input series, the
%     second output series corresponds to the second input series, and so
%     forth. Each output series is synchronized such that the last (most 
%     recent) observation is placed in the last row upon return.
%

%
% Pre-allocate arrays and determine sizes of input time series.
%

nSeries = numel(varargin);
nRows   = zeros(nSeries, 1);
nCols   = zeros(nSeries, 1);

for i = 1:nSeries
    [nRows(i), nCols(i)] = size(varargin{i});
end

X = zeros(max(nRows), sum(nCols));     % The composite array

%
% Pack each input time series into the composite data array (X), and  
% synchronize the data such that the last (most recent) observation of each
% series is placed in the last row of the composite array.
%

nCols  = cumsum([1 ; nCols]);
nRowsX = size(X,1);

for i = 1:nSeries
    rows             = (nRowsX - nRows(i) + 1):nRowsX;
    columns          =  nCols(i):(nCols(i + 1) - 1);
    X(rows, columns) =  varargin{i};
end

%
% Identify only those rows of the composite data array with no missing 
% observations (i.e., rows with no NaN's).
%

iX = all(~isnan(X),2);

%
% Now extract the remaining observations of each time series. 
%
% Each series in the output list is a cleansed version of the corresponding 
% series in the input list with all NaN's removed. Although the number of 
% columns of each series is unchanged, the output will likely have fewer 
% rows than its corresponding input. 
%
% The actual data found in a given output series will contain only those
% observations of its input for which no missing observations were found in
% ANY input series.
%

varargout =  cell(1,nargout);  % Only return what they ask for ...

for i = 1:nargout
    iSeries      = [false(nRowsX - nRows(i), 1) ; all(~isnan(varargin{i}),2)];
    rows         = iX & iSeries;                   % logicals
    columns      = nCols(i):(nCols(i + 1) - 1);    % doubles
    varargout{i} = X(rows, columns);
end

end % End of Listwise Deletion Method

%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

function [covariance, scores, numEvals] = errorCovarianceOpg(nLogL, x0, solve, lb, ub)
%ERRORCOVARIANCEOPG Error covariance matrix of MLE estimates by OPG method
%
% Syntax:
%
%   covariance = errorCovarianceOpg(nLogL,x0,solve,lb,ub)
%
% Description:
%
%   Estimate the error covariance matrix of model parameters using the outer 
%   product of gradients (OPG) method. The gradient of the loglikelihood 
%   function is approximated by numerically calculating the scores using 
%   finite central differences. Valid only for parameters estimated by 
%   maximum likelihood, the method assumes that the peak of the loglikelihood 
%   objective function has been found within the interior of the allowable 
%   parameter space, such that no boundary constraints are actively enforced.
%
% Input Arguments:
%
%   nLogL - Negative loglikelihood objective to evaluate (function handle).
%
%   x0 - Column vector of optimized maximum likelihood parameter estimates.
%
%   solve - Boolean vector the same length as x0. True elements indicate 
%     that the coefficient in the corresponding element of x0 has been 
%     estimated; false elements indicate that the coefficient was held 
%     fixed at the value in the corresponding element of x0.
%
%   lb - Vector of coefficient lower bounds the same length as x0.
%
%   ub - Vector of coefficient upper bounds the same length as x0.
%
% Output Arguments:
%
%   Cov - Covariance matrix of the estimation errors associated with
%     parameter estimates obtained by maximum likelihood. The standard 
%     errors of the individual parameter estimates are the square root of 
%     the diagonal elements. 
%
% Note:
%
%   o This is an hidden, undocumented method. No error checks are
%     performed, and users should not invoke it directly.
%
% Reference:
%
%   [1] Hamilton, J. D. Time Series Analysis. Princeton, NJ: Princeton
%       University Press, 1994, pp. 142-144, 427, 660-661.
%
%   [2] Numerical Recipes in C. The Art of Scientific Computing.

x0                  = x0(:);       % Guarantee a column vector
numberOfVariables   = numel(x0);
[~,nLogLikelihoods] = nLogL(x0);   % Called to determine number of observations

%
% Create structure of flags for finite differences.
%

OPTIONS                     = [];
defaultopt                  = optimset('fmincon');
finDiffOpts.DiffMinChange   = optimget(OPTIONS,'DiffMinChange',defaultopt,'fast');
finDiffOpts.DiffMaxChange   = optimget(OPTIONS,'DiffMaxChange',defaultopt,'fast');
finDiffOpts.TypicalX        = ones(size(x0));

%
% Determine whether we use central or forward finite differences.
%
% Note: 
%
%   Central finite differences: finDiffOpts.FinDiffRelStep = eps^(1/3)
%   Forward finite differences: finDiffOpts.FinDiffRelStep = sqrt(eps)
%

finDiffOpts.FinDiffType     = 'central';
finDiffOpts.FinDiffRelStep  = repmat(eps^(1/3), numberOfVariables, 1);    % Assumes CENTRAL differences
finDiffFlags.fwdFinDiff     = strcmpi(finDiffOpts.FinDiffType,'forward'); % Set the Boolean flag

finDiffFlags.scaleObjConstr = false; % No scaling in this algorithm
finDiffFlags.chkFunEval     = false; % Don't validate function values
finDiffFlags.chkComplexObj  = true;  % Check whether objective function values are complex if chkFunEval is true
finDiffFlags.isGrad         = false; % Compute objective gradient, not Jacobian of a system
finDiffFlags.hasLBs         = false(numberOfVariables,1);
finDiffFlags.hasUBs         = false(numberOfVariables,1);

if ~isempty(lb)
    finDiffFlags.hasLBs = isfinite(lb); % Check for lower bounds
end
if ~isempty(ub)
    finDiffFlags.hasUBs = isfinite(ub); % Check for upper bounds
end

sizes.xShape = size(x0);                % input to finitedifferences()

%
% Call finite differences to calculate the scores.
%

N = 1:numberOfVariables;

[~, scores, ~, numEvals] = finitedifferences(x0, [], @dummyConstFunc, lb, ub,...
    [], nLogLikelihoods', zeros(0,1),N(solve), finDiffOpts, sizes, [], ...
    repmat(nLogLikelihoods,numel(x0),1), [], finDiffFlags, []);

covariance      =  zeros(length(x0));
j               =  solve;
covariance(j,j) =  pinv(scores(j,:)*scores(j,:)');
covariance      = (covariance + covariance')/2;

%
% Nested function.
%

  function [c,ceq] = dummyConstFunc(x)
    [~,c] = nLogL(x);
    ceq   = [];
  end

end % Error Covariance Method


end % Methods (Static, Hidden)

end % Class definition