function [Rsquared, R] = UnivariateStats(X, y, w)
%% This function simply calculates R and R-squared,
% with, or without weighting


%% R, Correlation Coefficient
if nargin==2 % unweighted
    % covariance method
    r_ = cov(X,y) ./ (std(X).*std(y));
    
    % alternate, long-hand method
    %Xterm = (X-mean(X))./std(X); 
    %yterm = (y-mean(y))./std(y);
    %r_longhand = (1./(length(X)-1)).*sum((Xterm).*(yterm));
     
    % alternate, built-in method
    %r_corr = corr(X,y); % built-in method
    
    % alternate, built-in but Spearman
    % r_spearman = corr(X,y,'type','Spearman'); % rank correlation
    
else % use weights
    r_ = weightedcorrs([X,y], w); % use external fn
end
R = (r_(1,2));

%% Rsquared, Coefficient of Determination
if nargin==2
    SSR = nansum(((X-nanmean(y))).^2);  
    SSE = nansum(((y - X)).^2);
    SST = nansum(((y - nanmean(y))).^2);
else
    % checked
    SSR =    nansum( w.*(X-nanmean(y)).^2);  
    SSE =    nansum( w.*(y - X) .^2);
    SST =    nansum( w.*(y - nanmean(y)) .^2); 
end

%alternate methods for calculating Rsq
%Rsq1 = SSR / SST; 
%Rsq2 = 1-(SSE / SST); % 
%Rsq3 = 1-(SSR / SST); % this is equivalent to calc in glmfastFit ?
%[Rsq1, Rsq2, Rsq3]
Rsquared = 1-(SSE / SST);

if 0
    figure(); scatter(X,y);
    hold on; 
    lsline;
    x_lim = xlim();
    plot([x_lim(1), x_lim(2)], [mean(y), mean(y)],'k-', 'linewidth', 1.5);  
end