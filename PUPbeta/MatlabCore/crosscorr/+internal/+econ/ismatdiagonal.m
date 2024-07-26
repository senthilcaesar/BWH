function ok = ismatdiagonal(A)
%ISMATDIAGONAL - Check if a matrix is a square diagonal matrix.
%   This undocumented function may be modified or removed in a future release.
%
%	ok = internal.econ.ismatdiagonal(A);
%
% Inputs:
%	A - A matrix.
%
% Outputs:
%	ok - 1 if matrix is square and exactly diagonal, 0 otherwise.
%
% Comments:
%	1) This function is strict about the input matrix being a diagonal matrix.
%   2) If the input matrix is not square, the returned value is 0. If the input matrix is not a
%      valid matrix, the returned value is also 0.

% Copyright 2008 The MathWorks, Inc.  

if nargin < 1 || isempty(A) || ~isnumeric(A) || ndims(A) > 2
	ok = false;
elseif size(A,1) ~= size(A,2)
    ok = false;
elseif size(A,1) == 1
    ok = true;
else
	B = A - diag(diag(A));
	ok = ~any(B(:));
end
