function [BinaryVector]=FindVecs(Segs,len)

%%% Find the sample number of rising and falling edge of a binary vector
BinaryVector=zeros(len,1);
for ii=1:size(Segs,2)
    if Segs(1,ii)==0
        Segs(1,ii)=1;
    end
    if Segs(2,ii)==0
        Segs(2,ii)=1;
    end
    if Segs(1,ii)>len
        Segs(1,ii)=len;
    end
    if Segs(2,ii)>len
        Segs(2,ii)=len;
    end
    
    BinaryVector(Segs(1,ii):Segs(2,ii))=1;
end
