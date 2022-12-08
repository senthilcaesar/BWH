function y = logit(x,invflag)

if ~exist('invflag')
    invflag=0;
end

if invflag==1
y = 1./(1 + exp(-x));
else
y = log(x./(1-x));
end