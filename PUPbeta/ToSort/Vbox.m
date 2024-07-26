function [disLoc,disHeight,b_LH,Flchange]=Vbox(Flow)
disLoc=zeros(size(Flow));
L=5;
tetha=0.8;
H=0.1;

b_LH=zeros(size(Flow));
maxmin1=zeros(size(Flow));
maxmin2=zeros(size(Flow));
maxmin3=zeros(size(Flow));
for ii=2*L+1:length(Flow)-L
    clear BoxIdx1 BoxIdx2 BoxIdx3
    BoxIdx1=find(Flow>=Flow(ii-L)-H & Flow<=Flow(ii-L)+H);
    BoxIdx1=BoxIdx1(BoxIdx1>=ii-2*L & BoxIdx1<ii-L);
    maxmin1(ii)=median(Flow(BoxIdx1));%max(Flow(BoxIdx1))-min(Flow(BoxIdx1));
    if isnan(maxmin1(ii))
        if ~isempty(BoxIdx1)
            maxmin1(ii)=Flow(BoxIdx1);
        else
            maxmin1(ii)=mean(Flow);
        end
    end
    
    BoxIdx2=find(Flow>=Flow(ii)-H & Flow<=Flow(ii)+H);
    BoxIdx2=BoxIdx2(BoxIdx2>=ii-L & BoxIdx2<ii);
    maxmin2(ii)=median(Flow(BoxIdx2));%max(Flow(BoxIdx2))-min(Flow(BoxIdx2));
    if isnan(maxmin2(ii))
        if ~isempty(BoxIdx2)
            maxmin2(ii)=Flow(BoxIdx2);
        else
            maxmin2(ii)=mean(Flow);
        end
    end
    
    b_LH(ii)=length(BoxIdx2);
    BoxIdx3=find(Flow>=Flow(ii+L)-H & Flow<=Flow(ii+L)+H);
    BoxIdx3=BoxIdx3(BoxIdx3>=ii & BoxIdx3<ii+L);
    if isnan(maxmin3(ii))
        if ~isempty(BoxIdx3)
            maxmin3(ii)=Flow(BoxIdx3);
        else
            maxmin3(ii)=mean(Flow);
        end
    end
    
    maxmin3(ii)=median(Flow(BoxIdx3));%max(Flow(BoxIdx3))-min(Flow(BoxIdx3));
    
end
b_LH(1:L)=L;

Flchange=abs(maxmin3-maxmin2)+abs(maxmin2-maxmin1);
disLoc(b_LH>0 & b_LH<tetha*L)=1;

disHeight=0;

if 0
    Fs = 125;
    Time=(0:length(Flow)-1)/Fs;
    figure(2); clf(figure(2));
    ax(1)=subplot(211);
    plot(Time,Flow);ylabel('Flow(L/S)');
    hold on;
    stairs(Time,disLoc,'r');
    hold off;
    ylim([-2 2])
    
    ax(2)=subplot(212);
    plot(Time,b_LH);
    axis tight;
    ylim([-1 20])
    linkaxes(ax,'x')
end
end
