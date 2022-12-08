function [FlowLimCertOut,PrOut2] = computeFLcertainty(FinalModelTable,BreathFLDataTable)

logit = @(p) log(p./(1-p));
logitinverse = @(p) 1./(1+exp(-p)); % also called sigmoid

FlowCertIn = nan(size(BreathFLDataTable,1), size(FinalModelTable,1)-2);

for ii = 3:size(FinalModelTable,1)
    varname = FinalModelTable.ftrName{ii};
    trExp = str2num([varname(end-2),'.',varname(end)]); %transform exponent
    
    % limit feature term to within upper and lower limits
    if 1
        if 0 % this is how it was done for FlowDrive
            % we apply limits to untransformed data, then apply transform
            % I'm not 100% sure this is correct
            temp = BreathFLDataTable.(varname(1:end-4));
            upper = FinalModelTable.upper(ii);
            lower = FinalModelTable.lower(ii);
            temp(temp<lower)=lower;
            temp(temp>upper)=upper;
            FlowCertIn(:,ii-2) = sign(temp).*...
                abs(temp).^trExp.*...
                FinalModelTable.betas(ii);
        else % this is how i've done it for FLcert
            % do the transform, then apply limits, becuase the limits were
            % got from data that was already transformed.
            % also, in mnrval, we don't times by beta
            temp = sign(BreathFLDataTable.(varname(1:end-4))).*...
                abs(BreathFLDataTable.(varname(1:end-4))).^trExp;
            upper = FinalModelTable.upper(ii);
            lower = FinalModelTable.lower(ii);
            temp(temp<lower)=lower;
            temp(temp>upper)=upper;
            FlowCertIn(:,ii-2) = temp;
        end
    else % no limit applied to feture terms
        % also, do not multiply by beta for mnrval
        FlowCertIn(:,ii-2) = sign(BreathFLDataTable.(varname(1:end-4))).*...
            abs(BreathFLDataTable.(varname(1:end-4))).^trExp;
    end
end

%[betas]=mnrfit(Amatrix2_flow(:,topFtrsfromSFFS(:,1)),HS_cat, 'model', 'ordinal', 'link','logit');
[pihat] = mnrval(FinalModelTable.betas,FlowCertIn, 'model', 'ordinal', 'link','logit');
[~,cat_f] = max(pihat(:,:),[],2);
FlowLimCertOutA = categorical(cat_f);

%% newly added code to match paper output
sumval = sum(pihat(:,1:2),2); %Pr of FL (possible + certain) direct from model
out = logit(1-sumval); % 1-Pr

%[kappa, agree, T1, T2, AUC, accur, T3] = PerformanceOfFtrGivenLabel(HS, out)

%note these are in the log odds space predicting normal
T1 = -3.0650; %threshold separating Normal and possible from certain
T2 = 0.0483; %threshold separating Normal from possible and certain
T3 = -1.0677; %threshold separating Normal from certain (possible removed)

out2 = out - T2; % adjust at optimal two-class threshold, >0 indicates not normal.
out2(out2<-10)=-10;
out2 = -out2;
PrOut2 = logitinverse(out2); %>0.5 indicates not normal
PrcutoffCertain = logitinverse(T2-T1); %>0.957 (this number) indicates certainly not normal CFL
Prcutoff2states = logitinverse(T2-T3); %>0.753 gives a 2-state decision, N vs CFL
Thres = [0.5 PrcutoffCertain];
FlowLimCertOut = sum(PrOut2<Thres,2)+1;


end


















