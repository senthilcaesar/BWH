function yvalues_filtered=filter121(yvalues,fold)
yvalues_filtered=yvalues;
if ~(exist('fold')==1)
    fold=1;
end

for ii=1:fold
    for i=1:length(yvalues)
        range=[(i-1) i i (i+1)];
        if i==1
            range(1)=[];
        elseif i==length(yvalues)
            range(length(range))=[];
        end
        yvalues_filtered(i)=nanmean(yvalues(range));
    end
    yvalues=yvalues_filtered;
end

return