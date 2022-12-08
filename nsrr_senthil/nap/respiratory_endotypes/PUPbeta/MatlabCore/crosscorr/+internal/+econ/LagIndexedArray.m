classdef (Sealed) LagIndexedArray
%LAGINDEXEDARRAY Create a LagIndexedArray utility object.
%
% Syntax:
%	OBJ = LagIndexedArray(coefficients, lags)
%
% Description:
%   Create a lag indexed array object by specifying the various coefficients and
%   the corresponding lags with which each coefficient is associated. Note that 
%   LagIndexedArray objects are not meant to be instantiated directly, but
%   rather are designed to support the creation and manipulation of lag operator 
%   polynomial (LagOp) objects.
%
% Input Arguments:
%   coefficients - The coefficients of a lag operator polynomial specified as a 
%     cell array of square matrices.
%
%   lags - Vector of integer lags associated with the polynomial coefficients. 
%     The number of elements of lags must be the same as the number of 
%     coefficients.
%
% Output Arguments:
%   OBJ - Lag indexed array (LagIndexedArray) object.
%
% See also LagOp.

% Copyright 2009 The MathWorks, Inc.

properties (Access = private)
  Coefficients
  Lags
  Dimension
end

methods (Access = public, Hidden)

function OBJ = LagIndexedArray(Coefficients, Lags, Dimension)

%
% Perform some limited error checking.
%

if ~(isempty(Coefficients) && isempty(Lags))

   LagOp.checkCoefficients(Coefficients);
   LagOp.checkLags(Lags);

   if numel(Coefficients) ~= numel(Lags)
      error(message('econ:LagIndexedArray:LagIndexedArray:InputSizeMismatch'))
   end

end

%
% Note that the coefficients are NOT sorted by lag index, which is needed to
% support assignments of, for example, the following:
%
%     A.Coefficients([4 0 1]) = A.Coefficients([1 3 0]);
%

OBJ.Coefficients = Coefficients;
OBJ.Lags         = Lags;

if nargin > 2
   if ~isempty(Coefficients) && (size(Coefficients{1}, 1) ~= Dimension)
      error(message('econ:LagIndexedArray:LagIndexedArray:InputDimensionMismatch'))
   end
   OBJ.Dimension = Dimension;
else
   OBJ.Dimension = size(Coefficients{1}, 1);
end

end

%
% * * * * * * * * * LagIndexedArray Utility Methods * * * * * * * * *
%

function N = numel(OBJ, varargin) 
  if isempty(varargin)
     N = 0;
  else
     if strcmp(':', varargin)
        N = numel(getLags(OBJ));
     else
        N = numel(varargin{1});
     end
  end
end

function lastLag = end(OBJ, ~, ~)
  lastLag = max(OBJ.Lags);
end

function lags = getLags(OBJ)
  lags = OBJ.Lags;
end

function dimension = getDimension(OBJ)
  dimension = OBJ.Dimension;
end

%
% * * * * * * * * * LagIndexedArray-to-LagOp Converter * * * * * * * * * *
%

function OBJ = LagOp(OBJ, varargin)

  if isempty(varargin(1:2:end))
%
%    No optional inputs, so just pass it off to the LagOp constructor.
%
     OBJ = LagOp(OBJ.Coefficients, 'Lags', OBJ.Lags);

  else  % Optional inputs are specified.

     iLags = find(strcmpi('Lags', varargin(1:2:end)));

     if any(iLags)
%
%       The user specified specific lags to extract from the original object, so
%       get the coefficients at these specific lags ONLY. 
%
%       Please note that, aside from some minor near-zero-coefficient filtering 
%       via the optional "Tolerance" parameter, many of the same "extraction 
%       style" LagIndexedArray-to-LagOp conversions may also be performed via
%       subscripted "()" assignments. 
%
%       For example, given some lag operator polynomial A(L), the following two 
%       statements produce the same results:
%
%       B = LagOp(A.Coefficients, 'Lags', [0 4]);  % Returns a LagOp object.
%       C = A.Coefficients([0 4]);                 % Returns a LagIndexedArray object.
%
        coefficients = getCoefficients(OBJ, varargin{2 * iLags});
        iTolerance   = find(strcmpi('Tolerance', varargin(1:2:end)));

        if any(iTolerance)   % Did they also specify a near-zero tolerance?
           OBJ = LagOp(coefficients, 'Lags'     , varargin{2 * iLags}, ...
                                     'Tolerance', varargin{2 * iTolerance});
        else
           OBJ = LagOp(coefficients, 'Lags', varargin{2 * iLags});
        end

     else
        OBJ = LagOp(OBJ.Coefficients, 'Lags', OBJ.Lags, varargin{:});
     end
  end

end

%
% * * * * * * Subscripted Referencing of LagIndexedArray Objects * * * * * * * 
%

function varargout = subsref(OBJ, args)
%SUBSREF Subscripted referencing of LagIndexedArray objects.
%
% For a LagIndexedArray object B and lag L, valid subscripted references are:
% 
%  B(...)    - Get a subset of the LagIndexedArray object.
%  B{...}    - Get a subset of n-by-n coefficients matrices.
%  B{L}(i,j) - Get elements (i,j) of a coefficient matrix at lag L.
%

nSubscripts = internal.econ.LagIndexedArray.checkSubscripts(args);

try
%
% When indexing into a LagIndexedArray, determine if the user requested all 
% coefficients via the shorthand ":" syntax so explicit error checking can be 
% avoided.
%
  if isnumeric(args(1).subs{:})
     lags = args(1).subs{:};      % Allow for explicit lag indexing.
     LagOp.checkLags(lags);       % Ensure lags are non-negative integers.
  else
     lags = getLags(OBJ);         % Allow for ":" syntax (i.e., get all user-specified lags).
  end

  C = getCoefficients(OBJ, lags); % Get the coefficients at specified lags.

  switch args(1).type

    case '{}'   % "{}" indexing returns a coefficient matrix (or part of one).
       if nSubscripts > 1
          [varargout{1:nargout}] = builtin('subsref', C{:}, args(end));
       else
          if nargout > numel(C)
             error(message('econ:LagIndexedArray:subsref:InsufficientLHSArgumentList'));
          else
             varargout = C;
          end
       end

    case '()'   % "()" indexing returns a LagIndexedArray object (preserving type).
       varargout = {internal.econ.LagIndexedArray(C, lags)};

    otherwise
       error(message('econ:LagIndexedArray:subsref:InvalidAccessType'));
  end

catch exception

  exception.throwAsCaller();

end

end    % SUBSREF

%
% * * * * * * Subscripted Assignment of LagIndexedArray Objects * * * * * * 
%

function OBJ = subsasgn(OBJ, args, value)
%SUBSASGN Subscripted assignment of LagIndexedArray objects.
%
% For a LagIndexedArray object B and lag L, valid assignments are:
%
%  B         =  {...} - Assign all coefficients of a LagIndexedArray.
%  B(...)    =  {...} - Assign subset of coefficients of a LagIndexedArray.
%  B         = B(...) - Assign all coefficients of a LagIndexedArray from 
%                       another LagIndexedArray.
%  B(...)    = B(...) - Assign subset of coefficients of a LagIndexedArray from 
%                       another LagIndexedArray.
%  B{L}      =  [...] - Assign subset of n-by-n coefficient matrices.
%  B{L}(i,j) =  [...] - Assign elements of a coefficient matrix.
%

nSubscripts = internal.econ.LagIndexedArray.checkSubscripts(args);

try

   switch args(1).type

     case '{}'

        if ~isnumeric(value)    % B{L} = {...} assignments are errors.
           error(message('econ:LagIndexedArray:subsasgn:InvalidAccessType'));
        end

        if nSubscripts > 1
%
%          Convert a partial matrix into a full/square coefficient matrix
%          appropriate for assignment of the form A{L}(i,j) = [ ].
%
           C     = getCoefficients(OBJ, args(1).subs{:});
           value = builtin('subsasgn', C{:}, args(end), value);

        end

    case '()'
%
%       LagIndexedArray object "()" assignment requires the RHS to be either a 
%       cell array of coefficients or another LagIndexedArray object. So, for a 
%       LagIndexedArray object B(L), the following assignments are acceptable:
%
%            (1) B(...) =  {...}
%            (2) B(...) = B(...)
%
        if ~(isa(value, 'internal.econ.LagIndexedArray') || iscell(value))
           error(message('econ:LagIndexedArray:subsasgn:InvalidAccessType'));
        end

        if isa(value, 'internal.econ.LagIndexedArray')
           value = value.Coefficients;
        end

     otherwise
        error(message('econ:LagIndexedArray:subsasgn:InvalidAccessType'));

   end

   LHS_Lags  = args(1).subs{:}; 
   [C, lags] = setCoefficients(OBJ, value, LHS_Lags);
   OBJ       = internal.econ.LagIndexedArray(C, lags);

catch exception

   exception.throwAsCaller();

end

end    % SUBSASGN

%
% * * * * * * * * * * Display of LagIndexedArray objects * * * * * * * * *
%

function disp(OBJ)

%
% Get information about the specific lags with which the object was created as
% well as the order in which the lags were specified.
%

D      = OBJ.Dimension;  % Get the Dimension.
lags   = OBJ.Lags;       % Explicit lags (including the order!) used to create the object.
nTotal = numel(lags);    % Total # of lags used to create the object.

if nTotal == 0
   S1 = sprintf('%d-D Lag-Indexed Cell Array with', D);   % OBJ created w/out any non-zero coefficients
elseif nTotal <= 12
   s  = repmat(' %d', 1, nTotal);
   S1 = sprintf(['%d-D Lag-Indexed Cell Array Created at Lags [' s(2:end) '] with'], D, lags);
elseif sum(diff(lags)) == (nTotal - 1)
   S1 = sprintf('%d-D Lag-Indexed Cell Array Created at Lags [%d ... %d] with', D, lags(1), lags(end));
else
   S1 = sprintf('%d-D Lag-Indexed Cell Array Created at %d Lags with', D, nTotal);
end

%
% Get information about the lags with non-zero coefficient matrices.
%

indices = true(numel(lags),1);

for L = 1:numel(indices)
     indices(L) = any( (abs(OBJ.Coefficients{L}(:)) > 0) | (isnan(OBJ.Coefficients{L}(:))) );  % Assume zero-tolerance here.
end

nonZeroLags = unique(lags(indices),'legacy');                    % This sorts the unique lags too.
nNonZero    = numel(nonZeroLags);

if nNonZero == 0
   S2 = ' Zero Coefficients at All Lags.\n';
elseif nNonZero <= 12
   s  = repmat(' %d', 1, nNonZero);
   S2 = sprintf([' Non-Zero Coefficients at Lags [' s(2:end) '].\n'], nonZeroLags);
elseif sum(diff(nonZeroLags)) == (nNonZero - 1)
   S2 = sprintf(' Non-Zero Coefficients at Lags [%d ... %d].\n', nonZeroLags(1), nonZeroLags(end));
else
   S2 = sprintf(' Non-Zero Coefficients at %d Lags.\n', nNonZero);
end

%
% Now summarize the information.
%

fprintf(['    '  S1  '\n'  '   '  S2]);

end

end    % Methods (Access = public, Hidden)


methods (Access = private)

%
% * * * * * * Get Coefficients of LagIndexedArray Objects * * * * * * 
%

function coefficients = getCoefficients(OBJ, RHS_Lags)

%
% Handles the syntax C = A.Coefficients([...]).
%

currentLags  = OBJ.Lags;
coefficients = cell(1, numel(RHS_Lags));
zeroMatrix   = zeros(OBJ.Dimension);

%
% At this level, traditional MATLAB logical indexing is used (array elements are 
% sequentially accessed by the corresponding 1-based element indices 1,2,...).
%

for L = 1:numel(RHS_Lags)
    iMatch = find(currentLags == RHS_Lags(L));
    if any(iMatch)
       coefficients{L} = OBJ.Coefficients{iMatch(1)};
    else
       coefficients{L} = zeroMatrix;    % Assign the implicit 0 matrix.
    end
end

end    % getCoefficients

%
% * * * * * * Set Coefficients of LagIndexedArray Objects * * * * * * 
%

function [coefficients, newLags] = setCoefficients(OBJ, Coefficients, LHS_Lags)

currentLags  = OBJ.Lags;                     % Lags ALREADY specified in OBJ.
newLags      = union(currentLags, LHS_Lags,'legacy'); % Lags that WILL BE specified in OBJ after exit.
coefficients = cell(1, numel(newLags));

if isa(Coefficients, 'double')
   Coefficients = {Coefficients};
end

if numel(Coefficients) == 1                  % Perform scalar expansion if necessary.
   Coefficients = Coefficients(ones(1,numel(LHS_Lags)));
end

if numel(Coefficients) ~= numel(LHS_Lags)
   error(message('econ:LagIndexedArray:setCoefficients:InvalidAssignment'));
end

%
% At this level, traditional MATLAB logical indexing is used (array elements are 
% sequentially accessed by the corresponding 1-based element indices 1,2,...).
%

for L = 1:numel(newLags)

    iMatch = find(LHS_Lags == newLags(L));

    if any(iMatch)
%
%      In the event of multiple matches, notice that the LAST matching lag is 
%      used in the assignment; this mimics the behavior of core MATLAB cell
%      arrays.
%
       coefficients{L} = Coefficients{iMatch(end)};
    else
       coefficients{L} = OBJ.Coefficients{currentLags == newLags(L)};
    end

end

end    % setCoefficients

end    % Methods (Access = private)

%
% * * * * * * * * * * * Error Checking Utilities * * * * * * * * * * *
% 

methods (Static, Hidden)

function [nSubscripts, errorStruct] = checkSubscripts(args) 

errorStruct.message    = 'Invalid subscript or access syntax.';
errorStruct.identifier = 'econ:LagIndexedArray:checkSubscripts:InvalidAccessType';

%
% Ensure multiple references are of acceptable forms.
%

nSubscripts = length(args);

if nSubscripts > 3
   error(message('econ:LagIndexedArray:checkSubscripts:InvalidAccessType'));
end

if ~any(strcmp(args(1).type, {'()' '{}'}))
   error(message('econ:LagIndexedArray:checkSubscripts:InvalidAccessType'));
end

end    % Error Checking.

end    % Methods (Static, Hidden)

end    % Class Definition
