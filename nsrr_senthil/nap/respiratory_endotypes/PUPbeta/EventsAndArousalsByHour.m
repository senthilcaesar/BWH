function T2 = EventsAndArousalsByHour(Evts)

%%
T = Evts.EpochT;
Isleep1 = find(T.Epochs<4,1);
T(1:Isleep1-1,:)=[];

clear T2
i=1
while 1
    li = 1 + (i-1)*120;
    ri1 = li + 120-1;
    ri=ri1;
    ri(ri>height(T))=height(T);
    Irange = li:ri;
    T2.Nepochs(i,1)=length(Irange);
    T2.TST(i,1) = sum(T.Epochs(Irange)<4)*0.5;
    T2.Narousals(i,1) = sum(T.NArousals(Irange));
    T2.Nevents(i,1) = sum(T.Nevents(Irange));
    if ri1>height(T)
        break
    end
    i=i+1;
end
T2 = struct2table(T2);
T2.AHI = T2.Nevents./T2.TST*60;
T2.ArI = T2.Narousals./T2.TST*60;

minsleep=10;
T2.AHI(T2.TST<minsleep)=NaN;
T2.ArI(T2.TST<minsleep)=NaN;
