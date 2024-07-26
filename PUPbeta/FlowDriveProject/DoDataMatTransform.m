function [Amatrix2, Labels] = DoDataMatTransform(Amatrix, FeatureNames, addextratransform)

I1=Amatrix>=0;
I2=Amatrix<0;

AmatrixSQRT=NaN*Amatrix;
AmatrixSQRT(I1) = Amatrix(I1).^0.5;
AmatrixSQRT(I2) = -(-Amatrix(I2)).^0.5;

AmatrixSQ=NaN*Amatrix;
AmatrixSQ(I1) = Amatrix(I1).^2;
AmatrixSQ(I2) = -(-Amatrix(I2)).^2;
if addextratransform
    AmatrixX=NaN*Amatrix;
    AmatrixX = exp(-Amatrix);
    AmatrixY=NaN*Amatrix;
    AmatrixY = exp(Amatrix);
end

% Apply appropriate names for ftrs with transforms
temp = FeatureNames.Name;
clear temp1 temp2 temp3 tempX tempY
for i=1:length(temp)
    temp1{i}=[temp{i} '_1p0'];
    temp2{i}=[temp{i} '_0p5'];
    temp3{i}=[temp{i} '_2p0'];
    if addextratransform
        tempX{i}=[temp{i} '_0pX'];
        tempY{i}=[temp{i} '_0pY'];
    end
end

% Make large matrix including square-root and square transformed data where possible
Amatrix2 = [Amatrix AmatrixSQRT AmatrixSQ];
Labels = [temp1 temp2 temp3]';

if addextratransform
    Amatrix2 = [Amatrix2 AmatrixX AmatrixY];
    Labels = [Labels;tempX';tempY'];
end

end