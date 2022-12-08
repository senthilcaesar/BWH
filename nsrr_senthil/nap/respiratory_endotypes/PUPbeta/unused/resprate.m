function Ttotest = resprate(y_data,dt)

N = length(y_data);
T0 = N*dt;
df = 1/T0;

tic;
%Detect respiratory frequency using volume trace
%Highpass filter
LF = 1/15; %original 1/10
HF = 3; %original 1/10
PYYupperi=round(15/df);
filter_order1 = 2;
[B_butter1,A_butter1] = butter(filter_order1,[LF HF]/(1/dt/2));
y_data_filtered1 = filtfilt(B_butter1,A_butter1,detrend(y_data));

Yfft = fft(y_data_filtered1);
PYY = Yfft.*conj(Yfft);
PYY(PYYupperi+1:end)=[];
F=0:df:df*(length(PYY)-1);
%center frequency approach
temp = cumsum(PYY); 
Ttotest = 1/F(find((temp/temp(end))>0.5,1));
if 0
    figure(11); plot(0:dt:dt*(N-1),y_data_filtered1-mean(y_data_filtered1))
    figure(12);
    subplot(2,1,1); plot(F,PYY);
    xlim([1/15 1/1.5]);
    subplot(2,1,2); plot(F,temp/temp(end));
    xlim([1/15 1/1.5]);
end
toc;