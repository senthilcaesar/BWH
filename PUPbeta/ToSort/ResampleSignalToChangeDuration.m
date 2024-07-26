%% This is a test script for artificially changing the duration of breaths

% newFsrate > Fs will stretch the signal duration, while maintaining y values
% newFsrate < Fs will shrink the signal duration, while maintaining y values
newFsrate = 200; % change this value and re-run (F5), e.g. 50, 100, 125, 200, 250, 500

% set up original signal
Fs = 125;
dt = 1/Fs;
secs = 3;
Time = 0:dt:Fs*secs;
Flow = sin(Time);
oldlength=length(Time);
[~, firstpeak] = max(Flow(1:Fs*secs));

% plot original signal
figure(1); clf(figure(1));
plot(Time, Flow, 'b'); hold on;
plot(Time(firstpeak), Flow(firstpeak), 'b^');
ylim([-1.1 1.1]); xlim([0 secs*2.1]);
ylabel('Signal (e.g. Flow)'); xlabel('Time (Seconds)');

% this section does the work. it was copied directly from MIFL
if newFsrate~=0
    FsFactor = newFsrate / Fs;
    Flow=interp1(linspace(1,length(Flow)*FsFactor, length(Flow)),Flow,1:length(Flow)*FsFactor,'spline'); 
    Time = 0:dt:length(Flow)/Fs;                % make new "Time"
    firstpeak = round(firstpeak .* FsFactor);   % realign the breath timing  
end

% handle odd lengths for plot
if length(Time) < oldlength
    Time(end:oldlength) = NaN;
    Flow(end:oldlength) = NaN;
end

% plot resampled signal
plot(Time(1:oldlength), Flow(1:oldlength), 'r');
plot(Time(firstpeak), Flow(firstpeak), 'r^');
legend('original', 'orig peak', 'adjusted', 'adj peak'); 
refline(0,0);
