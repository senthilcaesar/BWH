function [DeltaPes_]=AlignWithBreaths(DeltaPes,Time,Iflow,Ipes)

%% Make continuous data
% VE_rs = NaN*Time;
showfigures=0;
PesDelta_rs = NaN*Time;

for i=1:length(Time)
    %     try
    %     VE_rs(i) = VE(find(Time(Iflow.starti)<=Time(i),1,'last'));
    %     catch me
    %     end
    try
        PesDelta_rs(i) = DeltaPes(find(Time(Ipes.starti)<=Time(i),1,'last'));
    catch me
    end
end

if showfigures
    figure(100);
    %     ax100(1)=subplot(211);
    %     stairs(Time(Iflow.starti),VE); hold('on');
    %     plot(Time,VE_rs,'k'); hold('off');
    stairs(Time(Ipes.starti),DeltaPes); hold('on');
    plot(Time,PesDelta_rs,'k'); hold('off');
end

for i=1:length(Iflow.starti)
    DeltaPes_(i,1)=nanmedian(PesDelta_rs(Iflow.starti(i):Iflow.endi(i)));
    if mean(isnan(PesDelta_rs(Iflow.starti(i):Iflow.endi(i))))>0.5 %if half of breath is NaN then use NaN
        DeltaPes_(i,1)=NaN;
    end
end

if showfigures
    hold('on');
    stairs(Time(Iflow.starti),DeltaPes_,'r--');
    hold('off');
end



