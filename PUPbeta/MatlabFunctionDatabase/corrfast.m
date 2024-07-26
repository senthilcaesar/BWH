function [r,p] = corrfast(X,Y,Fkeep)
    
if ~exist('Fkeep')
    Fkeep=0.5;
end

    I = sum(isnan(X)|isnan(Y))>(Fkeep*size(X,1)); %remove if half of data are NaN
    
    % Mean subtract cols
    if 0
    X = bsxfun(@minus,X,nansum(X,1)./size(X,1));
    Y = bsxfun(@minus,Y,nansum(Y,1)./size(Y,1));
    else
    X = bsxfun(@minus,X,nansum(X,1)./sum(~isnan(X)));
    Y = bsxfun(@minus,Y,nansum(Y,1)./sum(~isnan(X)));    
    end
    
    % Normalize by the L2-norm (Euclidean) of Rows:
    X = X.*repmat(sqrt(1./max(eps,nansum(abs(X).^2,1))),[size(X,1),1]); 
    Y = Y.*repmat(sqrt(1./max(eps,nansum(abs(Y).^2,1))),[size(Y,1),1]);
    
    % Calc R
    r = nansum(X.*Y);
    
    r(I)=NaN;
    
    % p-values if required
    if nargout==2
        t = (r.*(r1-2).^0.5)./(1-r.^2).^0.5;
        p = 2*tcdf(abs(t),(r1-2),'upper');
        p(I)=NaN;
    end
    
    %inspired by:
    %https://www.mathworks.com/matlabcentral/fileexchange/63082-fast_corr,
    %fast_corr, version 1.1.0.0 (2.09 KB) by Elliot Layden
    
end