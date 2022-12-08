function [vOut, ivOut] = BBSelectManual(Time,vol,Flow,Snore,SnoreDB,ivIn1,vIn1,ivIn2,fsFlow)
%%
% The code below allows you to interactively correct any misplaced
% troughs in the volume signal.  If a trough is mismarked, then left click
% just to the left of it and then right click so that the vertical line for
% the cross hairs is located where the trough should be.  See figure.
l = 4000;
n = 1;
x=[];
% xx = [];
% yy = [];
Button = [];
plotSnore=1;
while n < length(vol)
    while 1% plots until no new input      
        interactive_vmin3(Time,vol,Flow,Snore,SnoreDB,ivIn1,ivIn2,n,l,plotSnore);
        [xt,y,button] = ginput;
        x = xt*fsFlow;
        
        % make deletions
        xdel = x(button==1);
        for ii = 1:length(xdel)
            ind = find(ivIn1>round(xdel(ii)), 1, 'first'); %returns the larger index
            ivIn1(ind) = [];
            vIn1(ind) = [];
        end
        
        % make adds
        xadd = x(button==3);
        for ii = 1:length(xadd)
            if xadd(ii) < ivIn1(end)
                ind = find(ivIn1>round(xadd(ii)), 1, 'first'); %returns the larger index
                ivIn1 = [ivIn1(1:ind-1)' round(xadd(ii)) ivIn1(ind:end)']'; %ivmin must be a column vector
                vIn1 = [vIn1(1:ind-1)' vol(round(xadd(ii))) vIn1(ind:end)']';
            else
                ivIn1 = [ivIn1' round(xadd(ii))]';
                vIn1 = [vIn1' vol(round(xadd(ii)))]';
            end
        end
        close
        
        % break out of while loop once there is no input (i.e. x is empty)
        if isempty(x)
            break
        end
    end
%     
%     Button = [Button button'];
%     xx = [xx x'];
%     yy = [yy y'];
    n = n + l;
    close
end

% return these
vOut = vIn1;
ivOut = ivIn1;