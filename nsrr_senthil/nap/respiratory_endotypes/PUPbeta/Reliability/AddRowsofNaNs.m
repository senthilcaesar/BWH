function Out=AddRowsofNaNs(In,finalNrows)

Out=In;
Nr = size(Out,1);

if istable(Out)
    Out{Nr+1:finalNrows,:}=NaN;
else
    Out(Nr+1:finalNrows,:)=NaN;
end
