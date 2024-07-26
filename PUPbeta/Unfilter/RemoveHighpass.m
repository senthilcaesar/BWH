Time = DataEventHypnog_Mat(:,1);
Fs = 1/(Time(2)-Time(1));
Flow = DataEventHypnog_Mat(:,2);
w1 = 0.1;
n = 1;

Wn = w1/(Fs/2);
[a,b] = butter(n,Wn,'high');

Y = filter(b,a,Flow);
Y2 = detrend(Y);

if 1            
ax(1) = subplot(2,1,1);
plot(Time,Y2)
ax(2) = subplot(2,1,2);
plot(Time,Flow)
linkaxes(ax,'x')
end
