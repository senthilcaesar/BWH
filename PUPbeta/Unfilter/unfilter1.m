function Y2 = unfilter1(Flow,Fs,w1,order,ploton,newhighpass,newhighpassorder)

N=length(Flow);
Time = (0:1/Fs:(N-1)/Fs)';
% w1 = 0.1;
Wn = w1/(Fs/2);
[a,b] = butter(order,Wn,'high');
%[a,b] = besself(order,Wn,'high');
Y = filter(b,a,Flow);
Y2 = detrend(Y);

if ~isempty(newhighpass)
[b,a] = butter(newhighpassorder,newhighpass/(Fs/2),'high');
%[b,a] = besself(newhighpassorder,newhighpass/(Fs/2),'high');
%[b,a] = cheby2(newhighpassorder,newhighpass/(Fs/2),'high');
 %BUTTORD, BESSELF, CHEBY1, CHEBY2, ELLIP, FREQZ, FILTER.
Y2 = filter(b,a,Y2);
end

if ploton
ax(1) = subplot(2,1,1);
plot(Time,Flow)
ax(2) = subplot(2,1,2);
plot(Time,Y2)
linkaxes(ax,'x')
end


