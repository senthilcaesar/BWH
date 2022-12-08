function [Segs]=FindSegs(BinaryVector)

%%% Find the sample number of rising and falling edge of a binary vector
if size(BinaryVector,2)==1
    BinaryVector=BinaryVector';
end
Sig=BinaryVector;
Sig=[0 Sig 0];
SigTmp=Sig+[0 Sig(1:end-1)];
IdxTemp=find(SigTmp==1);
Starts=IdxTemp(1:2:end);
Ends=IdxTemp(2:2:end);
Segs(1,:)=Starts-1;
Segs(2,:)=Ends-2;
