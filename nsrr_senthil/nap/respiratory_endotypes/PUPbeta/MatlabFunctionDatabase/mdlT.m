function mdlCurrentT = mdlT(mdl2,remdownstreamterms) 
warning('off');

T = mdl2.Variables;
Yvar = mdl2.ResponseName;
varlistO = mdl2.Coefficients.Properties.RowNames(2:end)';
varlist=varlistO;

%swap 'X^2' for 'X:X'
temp = find(endsWith(varlist,'^2')==1);
for i=temp
    temp2 = varlist{i}(1:end-2)
    varlist{i} = [temp2 ':' temp2];
end

varlistCurrent = varlist;
LRp = nan(length(varlist),1);

if ~exist('remdownstreamterms')
    remdownstreamterms=1;
end

%swap 'X^2' for 'X:X'
temp = find(endsWith(varlistCurrent,'^2')==1);
for i=temp
    temp2 = varlistCurrent{i}(1:end-2)
    varlistCurrent{i} = [temp2 ':' temp2];
end

    clear mdlReduced
    for i=1:length(varlistCurrent)
        varlist2 = varlistCurrent;
        varCurrent = varlist2(i);
        idxrow = find(varlist==string(varCurrent));
        
        
            labelset = varlistCurrent;
            if remdownstreamterms
            DownstreamTerm = find(contains(labelset,labelset(idxrow)) & (contains(labelset,':') | contains(labelset,'^'))==1);
            remterms = unique([i DownstreamTerm]);
            else
            remterms =  i;
            end
            
        varlist2(remterms)=[];    
        mdlReduced{idxrow} = fitglm(T,[Yvar ' ~ ' listtoeqn(varlist2)],'Distribution',mdl2.Distribution.Name);
        
         
        
        try
            [~,LRp(idxrow,1)] = lratiotest(mdl2.LogLikelihood,mdlReduced{idxrow}.LogLikelihood,mdl2.NumEstimatedCoefficients-mdlReduced{idxrow}.NumEstimatedCoefficients);
        catch me
           LRp(idxrow,1) = NaN;
        end
    end

mdlCurrentT = mdl2.Coefficients;
Ilink = sum((mdl2.Coefficients.Properties.RowNames == string(varlistO)).*repmat(1:length(varlistO),height(mdlCurrentT),1),2);
mdlCurrentT.pValueLL = [NaN;LRp(Ilink(2:end))];
mdlCurrentT.Lower = mdlCurrentT.Estimate - 1.96*mdlCurrentT.SE;
mdlCurrentT.Upper = mdlCurrentT.Estimate + 1.96*mdlCurrentT.SE;

try
tempdata = T{:,varlist(Ilink(2:end))};
tempcatvars = nanmean(tempdata==1 | tempdata==0)==1;
TwoSD = 2*nanstd(tempdata)';
TwoSD(tempcatvars)=1;
mdlCurrentT.TwoSD = [NaN;TwoSD];
mdlCurrentT.Beta2SD = mdlCurrentT.Estimate.*mdlCurrentT.TwoSD;
mdlCurrentT.SE2SD = mdlCurrentT.SE.*mdlCurrentT.TwoSD;
mdlCurrentT.Lower2SD = mdlCurrentT.Beta2SD - 1.96*mdlCurrentT.SE2SD;
mdlCurrentT.Upper2SD = mdlCurrentT.Beta2SD + 1.96*mdlCurrentT.SE2SD;
if mdl2.Distribution.Name=="Binomial"
    mdlCurrentT.OR2SD = exp(abs(mdlCurrentT.Beta2SD)).*sign(mdlCurrentT.Beta2SD);
    mdlCurrentT.ORLower2SD = exp(abs(mdlCurrentT.Beta2SD)- 1.96*mdlCurrentT.SE2SD).*sign(mdlCurrentT.Beta2SD);
    mdlCurrentT.ORUpper2SD = exp(abs(mdlCurrentT.Beta2SD)+ 1.96*mdlCurrentT.SE2SD).*sign(mdlCurrentT.Beta2SD);
end
end