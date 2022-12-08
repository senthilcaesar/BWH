%% TTransFromFlow
% Find Transition period using Flow trace
% This will find two transition periods (TiTrans and TeTrans:
%   (1) one between insp(N) and exp(N) (usually very short)
%   (2) one between exp(N) and insp(N+1) (longer, known as the main TTrans)
%
% The basic process is as follows:
% Previously, find BB timing in Flow using VentilationFromFlow function.
% For 1 to length(BB) 
%   find AUC between BB_start and BB_end
%   calculate required size as threshold*AUC
%   find best elapsed time to reach required size
%   note BB_thres_start and BB_thres_end
% 
% For 1 to length(BB)-1 % can't find a Ttran on the last breath, becuase it
% requires knowledge of the next breath (doesn't exist in this sequence).
%   Ttran(N) = BB_thres_start(N+1) - BB_thres_end(N)

function [BB_Ttrans,TiTrans,TeTrans] = TTransFromFlow( time, BBs, Flow, range, Apnea_B, dt)
%%   BreathTimeAUC Determines breath time based on area under curve.
% Input:
% time - the time base of the Flow data
% BBs - original breath timing 
%   BB_original = [BB_i_start BB_i_mid BB_i_end];
%   this method uses both start and end times for both insp and exp
% Flow - the airflow waveform
% range - nominated range to find volume for
% Apnea_B - breaths added by VEfromFlow during apnea, or no airflow 

% Outputs:
% BB_Trans - index of start and end modified breath times
%   BB_Ttrans = [ BBAUC_i_start BBAUC_i_end BBAUC_e_start BBAUC_e_end];
% TiTrans - the inspiratory transition period
% TeTrans - the expiratory transition period

verbose = 0; % switch to display messages throughout code execution
ShowFigures = 0; % summary figure(s)
ShowExtraFigures = 0; % per breath figure(s)

BB_i_start = BBs(:,1);
BB_i_end = BBs(:,2);
BB_e_start = BBs(:,2);
BB_e_end = BBs(:,3);

if ~isnan(BB_i_start)
    N_Breaths = length(BB_i_start);   
    % preallocated variables
    BBAUC_i_start = zeros(N_Breaths,1);
    BBAUC_i_end = zeros(N_Breaths,1);
    BBAUC_e_start = zeros(N_Breaths,1);
    BBAUC_e_end = zeros(N_Breaths,1);
    TotalAUC_I = zeros(N_Breaths,1);
    TotalAUC_E = zeros(N_Breaths,1);
    ThresAUC_I = zeros(N_Breaths,1);
    ThresAUC_E = zeros(N_Breaths,1);    
    for i = 1:N_Breaths
        if Apnea_B(i)
            % do not attempt to find values for these breaths
            % default back to original breath timing               
            BBAUC_i_start(i) = BBs(i,1);
            BBAUC_i_end(i) = BBs(i,2);
            BBAUC_e_start(i) = BBs(i,2);
            BBAUC_e_end(i) = BBs(i,3);
            continue
        end
        ydata_insp = Flow(BB_i_start(i):BB_i_end(i));
        ydata_exp = -Flow(BB_e_start(i):BB_e_end(i));
        xdata_insp = time(BB_i_start(i):BB_i_end(i));
        xdata_exp = time(BB_e_start(i):BB_e_end(i));
        if ShowExtraFigures
            figure(10); clf(figure(10));
            subplot(2,2,1);
            plot(xdata_insp, ydata_insp); ylabel('insp'); hold on;
            refline(0,0);
            subplot(2,2,2);
            plot(xdata_exp, -ydata_exp); ylabel('exp'); hold on;
            refline(0,0);
            str = ['Breath number: ' num2str(i)]; suptitle(str);
        end
        if 0 % debugging plot
            figure(100); clf(figure(100));
            plot(time(BB_i_start(i-2):BB_e_end(i+2)), Flow(BB_i_start(i-2):BB_e_end(i+2))); hold on;
            
            plot(time(BB_i_start(i-2)), Flow(BB_i_start(i-2)), 'r^');
            plot(time(BB_i_end(i-2)), Flow(BB_i_end(i-2)), 'rd');
            plot(time(BB_e_end(i-2)), Flow(BB_e_end(i-2)), 'k.');
            
            plot(time(BB_i_start(i-1)), Flow(BB_i_start(i-1)), 'r^');
            plot(time(BB_i_end(i-1)), Flow(BB_i_end(i-1)), 'rd');
            plot(time(BB_e_end(i-1)), Flow(BB_e_end(i-1)), 'k.');
            
            plot(time(BB_i_start(i)), Flow(BB_i_start(i)), 'r^');
            plot(time(BB_i_end(i)), Flow(BB_i_end(i)), 'rd');
            plot(time(BB_e_end(i)), Flow(BB_e_end(i)), 'k.');
            
            plot(time(BB_i_start(i+1)), Flow(BB_i_start(i+1)), 'r^');
            plot(time(BB_i_end(i+1)), Flow(BB_i_end(i+1)), 'rd');
            plot(time(BB_e_end(i+1)), Flow(BB_e_end(i+1)), 'k.');
            
            plot(time(BB_i_start(i+2)), Flow(BB_i_start(i+2)), 'r^');
            plot(time(BB_i_end(i+2)), Flow(BB_i_end(i+2)), 'rd');
            plot(time(BB_e_end(i+2)), Flow(BB_e_end(i+2)), 'k.');
            
            plot(time(BBAUC_i_start(i-2)), Flow(BBAUC_i_start(i-2)), 'g^');
            plot(time(BBAUC_i_end(i-2)), Flow(BBAUC_i_end(i-2)), 'g^');
            plot(time(BBAUC_e_start(i-2)), Flow(BBAUC_e_start(i-2)), 'gv');
            plot(time(BBAUC_e_end(i-2)), Flow(BBAUC_e_end(i-2)), 'gv');
            
            plot(time(BBAUC_i_start(i-1)), Flow(BBAUC_i_start(i-1)), 'g^');
            plot(time(BBAUC_i_end(i-1)), Flow(BBAUC_i_end(i-1)), 'g^');
            plot(time(BBAUC_e_start(i-1)), Flow(BBAUC_e_start(i-1)), 'gv');
            plot(time(BBAUC_e_end(i-1)), Flow(BBAUC_e_end(i-1)), 'gv');
           
            plot(time(BBAUC_i_start(i)), Flow(BBAUC_i_start(i)), 'g^');
            plot(time(BBAUC_i_end(i)), Flow(BBAUC_i_end(i)), 'g^');
        end
        
        if (~isempty(ydata_insp) && ~isempty(ydata_exp) && ~isempty(xdata_insp) && ~isempty(xdata_exp))
            % find interval of inspiratory signal
            %dt_i=(xdata_insp(end)-xdata_insp(1))/(length(xdata_insp)-1);
            dt_i = dt;
            % find interval of expiratory signal
            %dt_e=(xdata_exp(end)-xdata_exp(1))/(length(xdata_exp)-1);
            dt_e = dt;
            if 0 % old method
                % set up vector of inspiratory integral data (used to determine best region)
                AUCdata_insp = zeros(length(xdata_insp),1);
                for j = 1:(length(AUCdata_insp))
                    AUCdata_insp(j) = dt_i*ydata_insp(j);
                end
                % set up vector of expiratory integral data (used to determine best region)
                AUCdata_exp = zeros(length(xdata_exp),1);
                for j = 1:(length(AUCdata_exp))
                    AUCdata_exp(j) = dt_e*ydata_exp(j);
                end
                % find Total inspiratory AUC
                TotalAUC_I(i) = dt_i*sum(ydata_insp);
                % find Total expiratory AUC
                TotalAUC_E(i) = dt_e*sum(ydata_exp);
            else % new method
                if length(ydata_insp)==1
                    AUCdata_insp = 0;
                else
                    AUCdata_insp = cumtrapz(xdata_insp, ydata_insp);
                end  
                if length(ydata_exp)==1
                    AUCdata_exp = 0;
                else
                    AUCdata_exp = cumtrapz(xdata_exp, ydata_exp);
                end  
                AUCdata_exp = cumtrapz(xdata_exp, ydata_exp);
                TotalAUC_I(i) = max(AUCdata_insp);
                TotalAUC_E(i) = max(AUCdata_exp);
            end
            
            % find intended Threshold AUC, based on insp area
            ThresAUC_I(i) = TotalAUC_I(i) * range;         
          
            % find Threshold AUC, based on exp area
            ThresAUC_E(i) = TotalAUC_E(i) * range;

            % find shortest exp time to reach this volume
            [exp_time_part, ExpTstart, ExpTend] = FindTimeForVolume(AUCdata_exp, ydata_exp, ThresAUC_E(i), dt_e);
            
            % find shortest insp time to reach this volume
            [insp_time_part, InspTstart, InspTend] = FindTimeForVolume(AUCdata_insp, ydata_insp, ThresAUC_I(i), dt_i);
            
            % if no data back from FindTimeForVolume, return NaNs
            if (isempty(insp_time_part) || isempty(exp_time_part) || isempty(ExpTstart) || ...
                    isempty(ExpTend) || isempty(InspTstart) || isempty(InspTend))
                if verbose
                    str = ['No result returned for breath ', num2str(i) ];
                    disp(str);
                end
                % if no Ttran timing result for this breath, default 
                % back to original breath timing               
                BBAUC_i_start(i) = BBs(i,1);
                BBAUC_i_end(i) = BBs(i,2);
                BBAUC_e_start(i) = BBs(i,2);
                BBAUC_e_end(i) = BBs(i,3);
            else
                % reset timings
                % set new inspiratory start and end index
                BBAUC_i_start(i) = BB_i_start(i) + InspTstart;
                BBAUC_i_end(i) = BB_i_start(i) + InspTend;
                
                % set new expiratory start and end index
                BBAUC_e_start(i) = BB_e_start(i) + ExpTstart;
                BBAUC_e_end(i) = BB_e_start(i) + ExpTend; 
                
                if ShowExtraFigures
                    figure(10); % build upon existing figure
                    subplot(2,2,1); hold on;
                    plot(xdata_insp, ydata_insp); ylabel('insp'); hold on;
                    refline(0,0);
                    plot(time(BBAUC_i_start(i)),Flow(BBAUC_i_start(i)),'r^');
                    plot(time(BBAUC_i_end(i)),Flow(BBAUC_i_end(i)),'rv');
                    subplot(2,2,2); hold on;
                    plot(xdata_exp, -ydata_exp); ylabel('exp'); hold on;
                    refline(0,0);
                    plot(time(BBAUC_e_start(i)),Flow(BBAUC_e_start(i)),'r^');
                    plot(time(BBAUC_e_end(i)),Flow(BBAUC_e_end(i)),'rv');
                    subplot(2,2,3);
                    plot(AUCdata_insp);
                    subplot(2,2,4);
                    plot(AUCdata_exp);
                end
            end
        else % if no data in, return NaNs
            if verbose
                sprintf(['No signal data for breath ' num2str(i) ]);
            end
            % if no Ttran timing result for this breath, default
            % back to original breath timing
            BBAUC_i_start(i) = BBs(i,1);
            BBAUC_i_end(i) = BBs(i,2);
            BBAUC_e_start(i) = BBs(i,2);
            BBAUC_e_end(i) = BBs(i,3);
        end
    end
else % if no data in, return NaNs
    if verbose
        sprintf(['No data.']);
    end
    % no data in, default back to original breath timing
    BBAUC_i_start = BBs(:,1);
    BBAUC_i_end = BBs(:,2);
    BBAUC_e_start = BBs(:,2);
    BBAUC_e_end(i) = BBs(:,3);
end

%BB_Ttrans = [ BBAUC_i_start BBAUC_i_end BBAUC_e_start BBAUC_e_end];

% calculate TiTrans
if sum(isnan([(BBAUC_i_end); (BBAUC_e_start)])) == 0
    TiTrans = diff([time(BBAUC_i_end); time(BBAUC_e_start)], 1, 1);
else
    for i = 1:N_Breaths
        if sum(isnan([(BBAUC_i_end(i)); (BBAUC_e_start(i))])) == 0
            TiTrans(i) = diff([time(BBAUC_i_end(i)); time(BBAUC_e_start(i))], 1, 1);
        else
            TiTrans(i) = NaN;
        end
    end
end   
%TiTrans(TiTrans==0)=0.008; % set the minimum Tip as one sample

% calculate TeTrans
if sum(isnan([(BBAUC_e_end); (BBAUC_i_start)])) == 0
    TeTrans = diff([time(BBAUC_e_end(1:end-1)); time(BBAUC_i_start(2:end))],1,1);
else
    for i = 1:N_Breaths-1
        if sum(isnan([(BBAUC_i_end(i)); (BBAUC_e_start(i+1))])) == 0
            TeTrans(i) = diff([time(BBAUC_e_end(i)); time(BBAUC_i_start(i+1))], 1, 1);
        else
            TeTrans(i) = NaN;
        end
    end
end

% ToDo: a better way to handle this last value. Could add average value as 
% the last value, not sure that's a better option than just appending a NaN
%TeTrans = [TeTrans NaN];
if size(TeTrans,2)>size(TeTrans,1) % data is row, add at end
    TeTrans = [TeTrans nanmedian(TeTrans)];
else %size(TeTrans,1)>size(TeTrans,2) % data is vertical, add at bottom
    TeTrans = [TeTrans; nanmedian(TeTrans)];
end

if 0 % don't do this, keep as much data as possible
TiTrans(isnan(TeTrans))=NaN;
BBAUC_i_start(isnan(TeTrans))=NaN;
BBAUC_i_end(isnan(TeTrans))=NaN;
BBAUC_e_start(isnan(TeTrans))=NaN;
BBAUC_e_end(isnan(TeTrans))=NaN;
end

BB_Ttrans = [ BBAUC_i_start BBAUC_i_end BBAUC_e_start BBAUC_e_end];

if ShowFigures
    exclude=isnan(BBAUC_i_start);
    i_start = BBAUC_i_start; i_start(isnan(i_start))=[];
    i_end = BBAUC_i_end; i_end(isnan(i_end))=[];
    e_start = BBAUC_e_start; e_start(isnan(e_start))=[];
    e_end = BBAUC_e_end; e_end(isnan(e_end))=[];
    Tran = TeTrans; Tran(isnan(Tran))=-1;
    figure(10);clf(figure(10));
    plot(time, Flow); hold on;
    plot(time(i_start), Flow(i_start), 'g^');
    plot(time(i_end), Flow(i_end), 'gv');
    plot(time(e_start), Flow(e_start), 'rv');
    plot(time(e_end), Flow(e_end), 'r^');
    refline(0,0);
    for n = 1 : length(Tran)
        if exclude(n)
            continue
        end
        x = time(BBAUC_e_end(n));
        y = Flow(BBAUC_e_end(n));
        str = ([{num2str(n)},{num2str(TeTrans(n))}]);
        text(x,y,str);
    end
    
end

end
