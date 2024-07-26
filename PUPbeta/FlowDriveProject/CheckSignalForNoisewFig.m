%% Check if respiratory signal is all noise
function [f_scale, y] = CheckSignalForNoisewFig(Time, Signal, Fs)
% Test the signal for the predominant frequency
% If outside of (very wide) respiratory rate limits, reject as noise
% current lower limit set as 15 sec / cycle (or 0.0667 Hz)
% current upper limit set as 0.05 sec / cycle (or 20 Hz)
% This is good at picking up long periods of 50 Hz noise

%
% ToDo: 
% 1. run this in windows across full study
%   - option one, call this Fn as is in prior step with window length data
%   - option two, do the windowing here...
% 2. change the upper and lower limits (this is what Scotty is using)
%   - lower set as 10 sec / cycle (0.1 Hz)
%   - upper set as 1 sec / cycle (1 Hz)
%

if 0
Time = timewav;
Signal = respwav;
Fs = 125;
end

verbose = 0;
showplot = 0;
TF = 1; % defaut set to 1 to continue operation, reset to 0 if fault
%Fs_ = 1/(Time(2)-Time(1))
x = Signal;
x = x - nanmean(x); 
x(isnan(x))=0; % set zero for NaNs
N=length(x);
nfft = 2^nextpow2(N); % next larger power of 2
y = fft(x,nfft); % Fast Fourier Transform
y = abs(y.^2); % raw power spectrum density
y = y(1:1+nfft/2); % half-spectrum
[~,k] = max(y); % find maximum
f_scale = (0:nfft/2)* Fs/nfft; % frequency scale
if f_scale(k)<0.0667 || f_scale(k)>20
%if f_scale(k)>20
    % if predominant frequency is less than 0.0667 Hz (i.e 15 sec cycle) or
    % greater then 20 Hz (i.e. 0.05 sec cycle), then it is noise
    % the range in AA flutpower is [5,round(Fs/2)-1] which is [5 62]
    if verbose
        disp('more noise than signal');
    end
    if showplot
        figure(1); clf(figure(1));
        ax1(1)=subplot(3,2,1:4);
        plot(f_scale, BasicFilt(y,9)); hold on; % was 100, then 3, now 9
        axis('tight');grid('on');
        title('Dominant Frequency in Signal');
        xlim([0 5]);
        xlabel('Frequency (Hz)'); ylabel('Amplitude');
        ax1(2)=subplot(3,2,5:6);
        plot(Time, Signal, 'b');
        xlabel('Time'); ylabel('Signal');
        hold off;
        
        if 0
            figure(2451); clf(figure(2451)); fig = gcf;
            fig.Color = [1 1 1]; fig.Units = 'inches';
            fig.Position = [2  2   4   5];
            plot(f_scale, BasicFilt(y,2)); hold on; % was 100, then 3, now 9
            set(gca,'box','off','tickdir','out', 'FontSize', 12, 'FontName','Arial Narrow');
            ylabel('Amplitude', 'FontSize', 18, 'FontName', 'Arial Narrow');
            xlabel('Frequency (Hz)', 'FontSize', 18, 'FontName', 'Arial Narrow');
            xlim([0 20]);
            str = ['..\Figures\FFT_DemoFig']; % 
            %print(fig, str, '-dtiff', '-r1000');
        end
    end
    
    %figure(2);clf(figure(2));
    %obw(x,Fs);
    
    TF = 0; % return indicating this in noisy
    
    % pause;
    % could prompt user what to do
    % then set TF flag accordingly
end
 if 0
            figure(2451); clf(figure(2451)); fig = gcf;
            fig.Color = [1 1 1]; fig.Units = 'inches';
            fig.Position = [12  4.5   4   5];
            plot(f_scale, BasicFilt(y,2)); hold on; % was 100, then 3, now 9
            set(gca,'box','off','tickdir','out', 'FontSize', 12, 'FontName','Arial Narrow');
            ylabel('Amplitude', 'FontSize', 18, 'FontName', 'Arial Narrow');
            xlabel('Frequency (Hz)', 'FontSize', 18, 'FontName', 'Arial Narrow');
            xlim([0 20]);
            str = ['..\Figures\FFT_DemoFig']; % 
            %print(fig, str, '-dtiff', '-r1000');
        end

end
