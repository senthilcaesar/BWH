function [TAA]=CalcTAA(DataIn)
%% This function calculates waveform correlation
% as the point by point sum of products of the two waveforms, of set window
% length, which produces a single value. Initially, the reference waveform
% shifted left by the defined offset. The point by point sum of products is
% repeated after sliding the reference waveform one data point to the right.
% This is repeated for each point in the slide distance, and the resulting 
% vector is normalised (-1 to +1). The value at mid-slide (which occurs
% when the two waveforms are naturally aligned) is taken as the value for 
% the current step. This is repeated for every step in the entire waveform.
% This mimics the Spike method of calculating waveform correlation.
% However, the smoothing/filtering used here is not performed in Spike.

%%
% Inputs
%  DataIn
%   (:,1) - Time vector
%   (:,2) - Thor data
%   (:,3) - Abdo data
%
% Outputs
%  TAA(:,1) - Thoracoabdominal asynchrony, of length(Time)
% 

%% Extra thoughts
% It would be nice to do a 'bake-off' between TAA calculated by DLM script
% within Spike, and the TAA calculated here.

%% pre-processing 
% if one or both RIP bands are missing, just return with NaN
% ToDo: expand this to include cases were one or both bands are just noise.
if all(DataIn(:,3)==-1) || all(DataIn(:,2)==-1)
    TAA = NaN(size(DataIn(:,3),1),1);
    %disp('RIP data missing, skipping TAA calculation');
    return
end

%% options
ShowPerBreathFigures = 0;   % this is only possible if not run within parfor
ShowSummaryFigures = 0;
FilterFactor = 0.08;        % number of points in MA filter is FilterFactor*fs
Verbose = 0;                % set to 1 for processing details, or zero for quiet operation
%% data
dt = DataIn(2,1)-DataIn(1,1);
fs = 1 / dt;
abdo_filt = normalise_minus1to1(BasicFilt(detrend(DataIn(:,2)),round(FilterFactor*fs)));
thor_filt = normalise_minus1to1(BasicFilt(detrend(DataIn(:,3)),round(FilterFactor*fs)));

%% check for noise in data
AbdoNoise = CheckSignalForNoise(DataIn(:,1), abdo_filt, fs); 
ThorNoise = CheckSignalForNoise(DataIn(:,1), thor_filt, fs);
if (AbdoNoise==0) || (ThorNoise==0)
    TAA = NaN(size(DataIn(:,3),1),1);
    if Verbose; disp('RIP data noisy, skipping TAA calculation'); end
    return
end


%% settings
% notes on output, for specified window, offset, slide and step settings
% W, O, S, S - comments
% 2, 1, 2, 1 - looks noisy, higher frequency changes
% 3, 1, 2, 1 - pretty good
% 4, 1, 2, 1 - probably the best
% 5, 1, 2, 1 - settings used in spike version
% 6, 2, 4, 1 - overly smoothed
W = 4;
O = 1;
S = 2;
St= 1;
window_length = single(W*fs);   % average Ttot
offset = single(O*fs);          % 1/4 average Ttot
slide = single(S*fs);           % 1/2 average Ttot
step = single(St*fs);           % Should not be affected by Ttot 

% ToDo: test if this improves with settings determined by pt's own Ttot.

%% processing
abdo = [abdo_filt;zeros(offset,1)]; % pad out with zeros, because of offest
thor = [zeros(offset,1);thor_filt;zeros(slide,1)]; % pad out with zeros
data_length=round(length(abdo_filt));
numsteps = round((data_length-window_length)/step);
corr_val = NaN(numsteps,1);
if Verbose; disp('Doing TAA calculation'); end
parfor n=1:1:numsteps
    win_st = n*step-step+1;
    win_nd = win_st+window_length;
    abdo_ = (abdo(win_st:win_nd));
    thor_ = (thor(win_st:win_nd+slide));
    sum_val = NaN(slide,1); 
    for m=1:1:slide
        thor__ = thor_(m:m+window_length);
        sum_val(m)= sum(abdo_.*thor__);
    end
    
    sum_val=normalise_minus1to1(sum_val);
    corr_val(n) = sum_val(ceil(length(sum_val)/2));
    
%     if ShowPerBreathFigures
%         figure(1); clf(figure(1));
%         ax(1)=subplot(2,1,1);
%         plot(sum_val); hold on;
%         
%         ax(2)=subplot(2,1,2);
%         plot(abdo_,'r-');hold on;
%         plot(thor_(offset:end-offset),'b-');
%         box off
%     end

end

%% normalise and pad the start and end with zeros
corr_val = [zeros(round(offset/fs),1);normalise_minus1to1(corr_val);zeros(round(offset/fs),1)];

%% interp to full length and smooth by fs
if any(isnan(corr_val))
    TAA = NaN(data_length,1);
else
    TAA=BasicFilt(interp1(linspace(1,data_length, length(corr_val)),(corr_val),1:data_length,'spline'),fs);
end

if ShowSummaryFigures
    figure(2); clf(figure(2));
    ax(1)=subplot(2,1,1);
    plot(TAA); hold on;
    refline(0,0);
    ax(2)=subplot(2,1,2);
    plot(abdo_filt,'r-'); hold on;
    plot(thor_filt,'b-');
    linkaxes(ax, 'x');
end

end