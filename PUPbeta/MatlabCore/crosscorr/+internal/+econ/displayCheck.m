function OK = displayCheck(display)
%DISPLAYCHECK Check value of 'display' parameter for univariate estimation.
%
% Syntax:
%
%   OK = displayCheck(display)
%
% Description:
%
%   Check that the value of 'display' parameter (used to display model 
%   estimation information) is valid.
%
% Input Argument:
%
%   display - String or cell vector of strings indicating what information
%             to display in the command window. 
%
% Output Argument:
%
%   OK - Logical flag. If TRUE, then the input display is valid. If 
%        validation fails, an error is thrown.
%

if ~isvector(display)

   error(message('econ:internal:econ:displayCheck:DisplayFlagNonVector'))

elseif isnumeric(display) || (iscell(display) && any(cellfun(@isnumeric,display)))

   error(message('econ:internal:econ:displayCheck:DisplayFlagNumeric'))

elseif ~all(ismember(lower(display),{'off','params','iter','diagnostics','full'}))

   error(message('econ:internal:econ:displayCheck:DisplayFlagInvalid'))

else

   OK = true;

end

end