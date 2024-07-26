function ExpArrayOut = normalizeExp(ExpArrayIn)

ExpArrayOut = nan(size(ExpArrayIn));
for ii = 1:size(ExpArrayIn,1)
    minExp = min(ExpArrayIn(ii,:));
    ExpArrayOut(ii,:) = ExpArrayIn(ii,:)./abs(minExp);
end