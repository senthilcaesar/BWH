%% TTransFromVolume
% Find Transition period using Volume trace
% This will find only one transition period:
%   (1) between exp(N) and insp(N+1).
%
% The basic process is as follows:
% Find BB timing in Volume 
%   use VentilationFromVolume function.
% For 1 to length(BB) 
%   find AUC between BB_start and BB_end
%   calculate required size as threshold*AUC
%   find best elapsed time to reach required size
%   note BB_thres_start and BB_thres_end
% 
% For 1 to length(BB)-1 % can't find a Ttran on the last breath, becuase it
% requires knowledge of the next breath (doesn't exist in this sequence).
%   Ttran(N) = BB_thres_start(N+1) - BB_thres_end(N)
%

function [ I_mod, TTrans] = TTransFromVol( time, I, Vol, range)
%%   BreathTimeAUC Determines breath time based on area under curve.
% Input:
% time - the time base of the Flow data
% I - original BB timing (I.i_start, I.i_end, I.e_start, I.e_end)
%   this method uses both start and end times for both insp and exp
% Vol - the volume waveform
% range - nominated range to find volume for

% Outputs:
% I_mod - index of start and end times that reflect the nominated range
% TTrans - the inter-breath transition period

verbose = 1; % switch to display error messages throughout code execution
ShowFigs = 0;

BB_i_start = I.starti;
BB_i_end = I.endi;
if ~isnan(BB_i_start)
    N_Breaths = length(BB_i_start);   
    % preallocated variables
    BBAUC_i_start = zeros(N_Breaths,1);
    BBAUC_i_end = zeros(N_Breaths,1);
    TotalAUC = zeros(N_Breaths,1);
    ThresAUC = zeros(N_Breaths,1);
    for i = 1:N_Breaths     
        ydata = Vol(BB_i_start(i):BB_i_end(i));
        xdata = time(BB_i_start(i):BB_i_end(i));
                
        if (~isempty(ydata) && ~isempty(xdata))
            % find interval of inspiratory signal
            dt_i=(xdata(end)-xdata(1))/(length(xdata));
                    
            % set up vector of inspiratory integral data (to determine best region)
            AUCdata = zeros(length(xdata),1);
            for j = 1:(length(AUCdata))
                AUCdata(j) = dt_i*ydata(j);
            end
            
            % find Total inspiratory AUC
            TotalAUC(i) = dt_i*sum(ydata);
                
            % find intended Threshold AUC, based on insp area
            ThresAUC(i) = TotalAUC(i) * range;         
            
            % find shortest time to reach this volume
            [time_part, Tstart, Tend] = FindTimeForVolume(AUCdata, ThresAUC(i), dt_i);
            
            % if no data back from FindTimeForVolume, return NaNs
            if (isempty(time_part) || isempty(Tstart) || isempty(Tend))
                BBAUC_i_start(i) = NaN;
                BBAUC_i_end(i) = NaN;
            else
                % reset timings
                % set new inspiratory start and end index
                BBAUC_i_start(i) = BB_i_start(i) + Tstart;
                BBAUC_i_end(i) = BB_i_start(i) + Tend;
 
            end
        else % if no data in, return NaNs
            if verbose
                sprintf(['No data at breath ' num2str(i) ]);
            end
            BBAUC_i_start(i) = NaN;
            BBAUC_i_end(i) = NaN;
        end
    end
else % if no data in, return NaNs
    if verbose
        sprintf(['No data.']);
    end
    BBAUC_i_start = NaN;
    BBAUC_i_end = NaN;
end

I_mod.starti = BBAUC_i_start;
I_mod.midi = I.midi; % this is unchanged from data passed into function
I_mod.endi = BBAUC_i_end;

% calculate TTrans
TTrans = diff([time(BBAUC_i_end(1:end-1)); time(BBAUC_i_start(2:end))],1,1);
TTrans = [TTrans mean(TTrans)];
% to do - a better way to handle this last value
% add mean value as last value, seems like a better option than just
% appending a zero or NaN

if ShowFigs
    figure();
    plot(time, Vol); hold on;
    plot(time(BBAUC_i_start), Vol(BBAUC_i_start), 'g^');
    plot(time(BBAUC_i_end), Vol(BBAUC_i_end), 'g^');
end

end

