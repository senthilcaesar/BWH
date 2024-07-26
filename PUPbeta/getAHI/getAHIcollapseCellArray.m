%% For AHIdata in a single cell, e.g. when using AHIonly=1

%%
function Array = getAHIcollapseCellArray(AHIdata)
N=size(AHIdata,2);
Empty = zeros(N,1);
for i=1:N
    if isempty(AHIdata{1,i})
        Empty(i)=1;
    end
end
I=find(Empty==0);
Nrows = NaN*zeros(length(I),1);
for i=1:length(I)
    Nrows(i) = length(AHIdata{1,I(i)})
end
M = max(Nrows);
Array = NaN*zeros(N,M);
for i=1:length(I)
    Array(I(i),1:Nrows(i)) = AHIdata{1,I(i)};
end