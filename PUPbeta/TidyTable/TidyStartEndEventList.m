function [I,I2] = TidyStartEndEventList(I,I2,N)

if isempty(I) && isempty(I2) %does not handle all night being an event or "good", not enough information provided, needs an exception    
else
    if isempty(I) 
        I=1;
    end
    if isempty(I2)
        I2=N;
    end
    if I(1)>I2(1)
        I=[1;I];
    end
    if I(end)>I2(end)
        I2=[I2;N];
    end
end
