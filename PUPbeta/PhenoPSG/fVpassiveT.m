function VpassiveT = fVpassiveT(Vpassive)
temp = 1-(Vpassive/100);
temp(temp<0)=0;
VpassiveT=100*(1-(temp.^0.5));

