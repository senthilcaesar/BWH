function ArThresT = fArThresT(ArThres,Inverse)
temp = (ArThres/100)-1;
temp(temp<0)=0;
ArThresT=100*(1+(temp.^0.5));

if exist('Inverse') && Inverse==1
    temp=100*(((ArThres/100-1).^2)+1);
    ArThresT=temp;  
end

