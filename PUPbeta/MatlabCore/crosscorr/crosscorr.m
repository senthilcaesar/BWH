function varargout = crosscorr(varargin)
%CROSSCORR Sample cross-correlation
%
% Syntax:
%
%   [xcf,lags,bounds] = crosscorr(y1,y2)
%   [xcf,lags,bounds] = crosscorr(y1,y2,param,val,...)
%   [xcf,lags,bounds,h] = crosscorr(...)
%   [xcf,lags,bounds,h] = crosscorr(ax,...)
%   crosscorr(...)
%
% Description:
%
%   Compute the sample cross-correlation function (XCF) between univariate,
%   stochastic time series y1 and y2. CROSSCORR optionally plots the XCF
%   sequence with confidence bounds.
%
% Input Arguments:
%
%   y1 - First vector of observations of a univariate time series. The last
%        element of y1 contains the most recent observation.
%
%   y2 - Second vector of observations of a univariate time series. The
%        last element of y1 contains the most recent observation.
%
%   ax - Axes object in which to plot. If unspecified, CROSSCORR plots to
%        the current axes (gca).
%
% Optional Input Parameter Name/Value Pairs:
%
%  'NumLags' Positive integer that determines the number of lags at which the 
%            XCF is computed. The lags used to compute the XCF are 0, +/-1, 
%            +/-2, ..., +/-NumLags. The default is min[20,min(N1,N2)-1], 
%            where N1 and N2 are the effective sample sizes of y1 and y2.
%
%  'NumSTD'  For computing confidence bounds, a nonnegative scalar multiple
%            specifying an interval of +/-(NumSTD) times the standard error
%            computed under the assumption that y1 and y2 are uncorrelated.
%            The default is 2 (approximate 95% confidence).
%
% Output Arguments:
%
%   xcf - Sample XCF. Vector of length 2*NumLags+1 of values computed at
%       lags 0, +/-1, +/-2, ..., +/-NumLags. The center element of xcf
%       contains the 0-lag cross-correlation.
%
%   lags - Vector of lag numbers of length 2*NumLags+1 used to compute xcf.
%
%   bounds - Two-element vector of approximate upper and lower confidence
%       bounds, assuming that the input series are uncorrelated.
%
%   h - Vector of handles to plotted graphics objects. CROSSCORR plots the
%       XCF when the number of output arguments is 0 or 4.
%
% Example:
%
%   % Create a random sequence of 100 Gaussian deviates and a delayed
%   % version lagged by 4 periods. Observe the XCF peak at the 4th lag:
%
%   x = randn(100,1);    % 100 Gaussian deviates ~ N(0,1)
%   y = lagmatrix(x,4);  % Delay x by 4 periods
%   y(isnan(y)) = 0;     % Replace NaNs with zeros
%   crosscorr(x,y)       % Observe peak at the 4th lag
%
% Reference:
%
%   [1] Box, G. E. P., G. M. Jenkins, and G. C. Reinsel. Time Series
%       Analysis: Forecasting and Control. 3rd edition. Upper Saddle River,
%       NJ: Prentice-Hall, 1994.
%
% See also AUTOCORR, PARCORR, FILTER.

% Copyright 2018 The MathWorks, Inc.   

% Preprocess varargin for target axes:

try
    
    [ax,args] = internal.econ.axesparser(varargin{:});
    
catch ME
    
    throw(ME)
    
end

% This function produces a single plot:

if ~isempty(ax) && ~isscalar(ax)
    
    error(message('econ:internal:econ:axesparser:InvalidParent'));
    
end

% Parse inputs and set defaults:

if numel(args) < 2
    
   error(message('econ:crosscorr:UnspecifiedInput'))
   
end

if (numel(args) > 2) && isnumeric(args{3}) % Deprecated positional inputs
    
    % Positional input syntax:
    % [...] = crosscorr(y1,y2,numLags,numSTD)

    y1 = args{1};
    y2 = args{2};
    args = args(3:end);
   
    try
       
        parser = checkPositionalInputs(y1,y2,args{:});
     
    catch ME
        
        throwAsCaller(ME)
     
    end
   
    numLags = parser.Results.NumLags;
   
else % Name-value pair inputs
    
    parser = inputParser;
    parser.addRequired ('y1'     ,     @(x) validateattributes(x, {'double'}, {'vector'}, '', 'first input series'));
    parser.addRequired ('y2'     ,     @(x) validateattributes(x, {'double'}, {'vector'}, '', 'second input series'));
    parser.addParameter('NumLags', 20, @(x) validateattributes(x, {'double'}, {'scalar' 'integer' '>'  0}, '', 'number of lags'));
    parser.addParameter('NumSTD' ,  2, @(x) validateattributes(x, {'double'}, {'scalar'           '>=' 0}, '', 'number of standard deviations'));
   
    y1 = args{1};
    y2 = args{2};
    args = args(3:end);
    
    try
        
        parser.parse(y1,y2,args{:});
     
    catch ME
        
        throwAsCaller(ME)
     
    end
    
    % Set lags:

    numLags = parser.Results.NumLags;

    if any(strcmpi('NumLags',parser.UsingDefaults)) % Default lags
       
        numLags = min(numLags,min(length(y1),length(y2))-1);
      
    else % User-specified lags
       
        if numLags > (min(length(y1),length(y2))-1)
          
            error(message('econ:crosscorr:InputTooLarge'))
         
        end
      
    end
   
end % Input parse

% Preprocess validated inputs:

y1        = parser.Results.y1;
y2        = parser.Results.y2;
rowSeries = (size(y1,1) == 1);
y1 = y1(:);
y2 = y2(:);
N = min(length(y1),length(y2)); % Sample size

numSTD    = parser.Results.NumSTD;

% Compute XCF:

y1 = y1-mean(y1);
y2 = y2-mean(y2);
L1 = length(y1);
L2 = length(y2);

if L1 > L2
    y2(L1) = 0;
elseif L1 < L2
    y1(L2) = 0;
end

% FFT is significantly faster than filtering for large data sets:

nFFT = 2^(nextpow2(max([L1 L2]))+1);
F = fft([y1(:) y2(:)],nFFT); 

ACF1 = ifft(F(:,1).*conj(F(:,1)));
ACF2 = ifft(F(:,2).*conj(F(:,2)));

xcf = ifft(F(:,1).*conj(F(:,2)));
xcf = xcf([(numLags+1:-1:1) (nFFT:-1:(nFFT-numLags+1))]);
xcf = real(xcf)/(sqrt(ACF1(1))*sqrt(ACF2(1)));

lags = (-numLags:numLags)';
bounds = [numSTD;-numSTD]/sqrt(N);

% Perform nargout-dependent operations:

nargoutchk(0,4)

if (nargout == 0) || (nargout == 4) % Create plot
    
    % Plot to gca if no parent axes is specified:

    if isempty(ax)
    
        ax = gca;
    
    end

    % Store NextPlot flag (and restore on cleanup):

    next = get(ax,'NextPlot');
    cleanupObj = onCleanup(@()set(ax,'NextPlot',next));

    % Plot the sample XCF:

    hPlot = stem(ax,lags,xcf,'filled','r-o','MarkerSize',4,'Tag','XCF');
   
    % Plot confidence bounds under the hypothesis that y is uncorrelated:

    x1 = ax.XLim(1);
    x2 = ax.XLim(2);
	set(ax,'NextPlot','add')
    hBounds = plot(ax,[x1 x1; x2 x2],[bounds([1 1]) bounds([2 2])],'-b',...
                   'Tag','Confidence Bound');
	hXAxis = plot(ax,[x1 x2],[0 0],'-k','Tag','X Axis');
    
    % Return "plot object":

    h = [hPlot;hBounds;hXAxis];
    
    % Modify axes properties conditional on NextPlot flag:
    
    ax.Tag = 'XCFPlot';

    switch next
    
        case {'replace','replaceall'}
    
            grid(ax,'on')
            xlabel(ax,'Lag')
            ylabel(ax,'Sample Cross Correlation')
            title(ax,'Sample Cross Correlation Function')
            
        case {'replacechildren','add'}
        
            % Do not modify axes properties
    
    end

end

if nargout > 0

    % Re-format outputs to conform to row/column orientation of y:

    if rowSeries

      xcf = xcf';
      lags = lags';
      bounds = bounds';

    end

end

% Suppress assignment to ans:

if (nargout > 0) && (nargout < 4)

    varargout = {xcf,lags,bounds};
    
elseif nargout == 4
    
    varargout = {xcf,lags,bounds,h};
    
end

% -------------------------------------------------------------------------
function parser = checkPositionalInputs(y1,y2,numLags,numSTD)

% Ensure the sample data are vectors:

[rows,columns] = size(y1);

if (rows ~= 1) && (columns ~= 1)
    
   error(message('econ:crosscorr:NonVectorSeries1'))
   
end

[rows,columns] = size(y2);

if (rows ~= 1) && (columns ~= 1)
    
   error(message('econ:crosscorr:NonVectorSeries2'))
   
end

N = min(length(y1),length(y2)); % Sample size

% Ensure numLags is a positive integer or set default:

if (nargin >= 3) && ~isempty(numLags)
    
   if numel(numLags) > 1
       
      error(message('econ:crosscorr:NonScalarLags'))
      
   end
   
   if (round(numLags) ~= numLags) || (numLags <= 0)
       
      error(message('econ:crosscorr:NonPositiveInteger'))
      
   end
   
   if numLags > (N-1)
       
      error(message('econ:crosscorr:InputTooLarge'))
      
   end
   
else
    
   numLags = min(20,N-1); % Default
   
end

% Ensure numSTD is a positive scalar or set default:

if (nargin >= 4) && ~isempty(numSTD)
    
   if numel(numSTD) > 1
       
      error(message('econ:crosscorr:NonScalarSTDs'))
      
   end
   
   if numSTD < 0
       
      error(message('econ:crosscorr:NegativeSTDs'))
      
   end
   
else
    
   numSTD = 2; % Default
   
end

parser.Results.y1      = y1;
parser.Results.y2      = y2;
parser.Results.NumLags = numLags;
parser.Results.NumSTD  = numSTD;