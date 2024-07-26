%Q: does including uncorrelated covariates reduce power?

%generate correlated variables

M=500;

Rsq = 0.2;
varLG = (0.1)^2;
    LGSD = (varLG)^0.5;
beta = 20; %
varNoise = (varLG*beta^2)*(1/Rsq - 1);
    NoiseSD = (varNoise)^0.5;

N=36;
Out = NaN*zeros(M,1);
for m=1:M

    LG = 0.5 + LGSD*randn(N,1);

    Response = 0.5 + beta.*LG + NoiseSD*randn(N,1);

    UncorrelatedX = randn(N,1);

    %mdl = fitglm([LG],Response);
    mdl = fitglm([LG UncorrelatedX],Response);
    mdl.Rsquared.Ordinary;
    Out(m)=mdl.Coefficients.pValue(2);
end

Power = sum(Out<0.05)/length(Out)

%%
%Q: does including uncorrelated covariates reduce power?

%generate correlated variables
tic
M=500;

% Rsq = 0.2;
% varLG = (0.1)^2;
    LGSD = 0.12; %LGSD = 0.14
% beta = 20; %
% varNoise = (varLG*beta^2)*(1/Rsq - 1);
%     NoiseSD = (varNoise)^0.5;

 % LGresponders=0.5308, LGNresponders=0.6477
 % LGrespondersSD=0.1099, LGNrespondersSD=0.168
 
N=36;
Out = NaN*zeros(M,1);

for m=1:M
    
    N1 = round(N/2);
    N2 = round(N/2);
    if 0
    LG1 = 0.53 + 0.110*randn(N1,1);
    LG2 = 0.65 + 0.168*randn(N2,1);
    else
    LG1 = 0.53 + LGSD*randn(N1,1);
    LG2 = 0.65 + LGSD*randn(N2,1);
    end
    Response1 = ones(N1,1);%0.5 + beta.*LG + NoiseSD*randn(N,1);
    Response2 = zeros(N2,1);
    
    LG = [LG1;LG2];
    Response = [Response1;Response2];
    
    UncorrelatedX = randn(N,1);
    UncorrelatedY = randn(N,1);
    UncorrelatedZ = randn(N,1);

    mdl = fitglm([LG],Response,'Distribution','binomial');
    
    %mdl = fitglm([LG UncorrelatedX],Response,'Distribution','binomial');
    %mdl = fitglm([LG UncorrelatedX UncorrelatedY UncorrelatedZ],Response,'Distribution','binomial');
    mdl.Rsquared.Ordinary;
    Out(m)=mdl.Coefficients.pValue(2);
    
end
toc
Power = sum(Out<0.025)/length(Out)





