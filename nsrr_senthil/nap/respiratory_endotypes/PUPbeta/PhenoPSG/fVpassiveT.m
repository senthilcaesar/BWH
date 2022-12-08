function VpassiveT = fVpassiveT(Vpassive,Inverse)
temp = 1-(Vpassive/100);
temp(temp<0)=0;
VpassiveT=100*(1-(temp.^0.5));

if exist('Inverse') && Inverse==1 %to check
    temp = 1-(Vpassive/100);
    VpassiveT=100*(1-(temp.^2));
end