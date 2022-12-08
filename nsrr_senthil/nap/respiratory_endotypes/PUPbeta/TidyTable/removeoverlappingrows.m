function Data1 = removeoverlappingrows(Data,cols,addNaNgaps)
%cols = [1 5 8]

w = Data(:,cols(1));
tL = Data(:,cols(2));
tR = tL+Data(:,cols(3));

tol=0.1;
Alltimes = unique(sort([tL;tR]));
delta = [NaN;diff(Alltimes)];
Alltimes(delta<tol)=[];
freq = Alltimes*0;
w_ = cell(length(freq),1);
for i=1:length(tL)
    a=tL(i);
    b=tR(i);
    Ir1=find(Alltimes>=(a-tol),1,'first');
    Il2=find(Alltimes<=(b+tol),1,'last');
    irange = Ir1:(Il2-1);
    for j=1:length(irange)
        w_{irange(j),1}(length(w_{irange(j)})+1)=w(i);
    end
end

Wmid = 0*Alltimes;
for i=1:length(w_)
    Wmid(i) = ceil(median(w_{i}));
end
keep=0*tL;
for i=1:length(tL)
    a=tL(i);
    b=tR(i);    
    Ir1=find(Alltimes>=(a-tol),1,'first');
    Il2=find(Alltimes<=(b+tol),1,'last');
    irange = Ir1:(Il2-1);
    keep(i)=w(i)==mode(Wmid(irange));
end

Data1 = Data;
Data1(keep==0,:)=[];

if 1 %occasionally, some rows were out of order
Data1 = sortrows(Data1,cols(2));
end

if addNaNgaps
tol2=1;
i=1;
M=size(Data1,2);
while i<(size(Data1,1)-1)
    if (Data1(i,cols(2))+Data1(i,cols(3))+tol2)<Data1(i+1,cols(2))
        Data1 = [Data1(1:i,:); NaN*ones(1,M); Data1((i+1):size(Data1,1),:)];
        %keyboard
        i=i+1;
    end
    i=i+1;
end
end
% figure();
% stairs(Data1(:,5),Data1(:,6));
