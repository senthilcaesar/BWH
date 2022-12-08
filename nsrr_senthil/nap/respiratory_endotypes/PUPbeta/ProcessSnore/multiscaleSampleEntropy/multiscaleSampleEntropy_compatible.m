function [ e, A, B ] = multiscaleSampleEntropy_compatible( x, m, r, tau )
%MULTISCALESAMPLEENTROPY
%
% Based on "Multiscale entropy analysis of biological signals"
% By Madalena Costa, Ary L. Goldberger, and C.-K. Peng
% Published on 18 February 2005 in Phys. Rev. E 71, 021906.
%
% This code was implemented by John Malik on 26 April 2017.
% Contact: john.malik@duke.edu

switch nargin
    case 1
        m = 2;
        r = 0.15;
        tau = 1;
    case 2
        r = 0.15;
        tau = 1;
    case 3
        tau = 1;
end

% coarse signal
x = x(:);
u = ceil(length(x)/tau);
y = zeros(tau, u);
y(1:length(x)) = x;
y = mean(y, 1);

% (m+1)-element sequences
N = length(y);
X = zeros(N - m, m + 1);
for i = 1:m + 1
    X(:, i) = y((i - 1) + 1:N - m + (i - 1));
end

% matching (m+1)-element sequences
A = sum(pdist(X, 'chebychev') < r * nanstd(x, 1));

% matching m-element sequences
X = X(:, 1:m);
B = sum(pdist(X, 'chebychev') < r * nanstd(x, 1));

% take log
if A == 0 || B == 0
    e = NaN;
    return
end
e = log(B / A);

end

