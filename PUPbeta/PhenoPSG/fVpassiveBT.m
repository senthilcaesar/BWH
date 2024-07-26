function Vpassive = fVpassiveBT(VpassiveT)
temp = 1-(VpassiveT/100);
Vpassive=100*(1-(temp.^2));

