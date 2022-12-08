function Out=AddRowsofZeros(In,finalNrows)

Out=In;
Nr = size(Out,1);

if istable(Out)
    Out{Nr+1:finalNrows,:}=0;
else
    Out(Nr+1:finalNrows,:)=0;
end
