% Call the function "get_vmin2" to find the troughs of the volume signal
ti = 2; te=2;
[vmin,ivmin] = get_vmin2(vol,ti,te,1.4,fsFlow);
tvmin = ivmin/fsFlow;

% Find the peaks of the volume signal
vmax = [];
ivmax = [];
for n = 2:length(ivmin)
    [maxvol,indmax] = max(vol(ivmin(n-1):ivmin(n)));
    indmax = (ivmin(n-1)+indmax)-1;
    vmax = [vmax maxvol];
    ivmax = [ivmax indmax];
end
vmax = vmax';
ivmax = ivmax';
tvmax = ivmax/fsFlow;
clear indmax maxvol n;

%%
% The code below allows you to interactively correct any misplaced
% troughs in the volume signal.  If a trough is mismarked, then left click
% just to the left of it and then right click so that the vertical line for
% the cross hairs is located where the trough should be.  See figure.
l = 9002;
n = 1;
x=[];
% xx = [];
% yy = [];
Button = [];
plotSnore=1;
while n < length(vol)
    while 1% plots until no new input      
        interactive_vmin3(Time,vol,Flow,Snore,SnoreDB,ivmin,ivmax,n,l,plotSnore);
        [xt,y,button] = ginput;
        x = xt*fsFlow;
        
        % make deletions
        xdel = x(button==1);
        for ii = 1:length(xdel)
            ind = find(ivmin>round(xdel(ii)), 1, 'first'); %returns the larger index
            ivmin(ind) = [];
            vmin(ind) = [];
        end
        
        % make adds
        xadd = x(button==3);
        for ii = 1:length(xadd)
            if xadd(ii) < ivmin(end)
                ind = find(ivmin>round(xadd(ii)), 1, 'first'); %returns the larger index
                ivmin = [ivmin(1:ind-1)' round(xadd(ii)) ivmin(ind:end)']'; %ivmin must be a column vector
                vmin = [vmin(1:ind-1)' vol(round(xadd(ii))) vmin(ind:end)']';
            else
                ivmin = [ivmin' round(xadd(ii))]';
                vmin = [vmin' vol(round(xadd(ii)))]';
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
clear n x y yy button;

%%
tvmin = ivmin/fsFlow;

clear ivmax tvmax vmax

% re-find the peaks of the volume signal after you have corrected the
% troughs
vmax = [];
ivmax = [];
for n = 2:length(ivmin)
    [maxvol,indmax] = max(vol(ivmin(n-1):ivmin(n)));
    indmax = (ivmin(n-1)+indmax)-1;
    vmax = [vmax maxvol];
    ivmax = [ivmax indmax];
end
tvmax = ivmax/fsFlow;
clear indmax maxvol n;



% Correct any tvmax's that equal tvmin's
for n = 1:length(ivmax)
    if ivmax(n) == ivmin(n)
        ivmax(n) = ivmin(n) + round((ivmin(n+1)-ivmin(n))/2); %put it halfway between the ivmins
    end
end

ivmax = ivmax';
vmax = vmax';
%% Now check and adjust vmax
% The code below allows you to interactively correct any misplaced
% troughs in the volume signal.  If a trough is mismarked, then left click
% just to the left of it and then right click so that the vertical line for
% the cross hairs is located where the trough should be.  See figure.
l = 9002;
n = 1;
x=[];
% xx = [];
% yy = [];
Button = [];
plotSnore=1;
while n < length(vol)
    while 1% plots until no new input      
        interactive_vmin3(Time,vol,Flow,Snore,SnoreDB,ivmax,ivmin,n,l,plotSnore);
        [xt,y,button] = ginput;
        x = xt*fsFlow;
        
        % make deletions
        xdel = x(button==1);
        for ii = 1:length(xdel)
            ind = find(ivmax>round(xdel(ii)), 1, 'first'); %returns the larger index
            ivmax(ind) = [];
            vmax(ind) = [];
        end
        
        % make adds
        xadd = x(button==3);
        for ii = 1:length(xadd)
            if xadd(ii) < ivmax(end)
                ind = find(ivmax>round(xadd(ii)), 1, 'first'); %returns the larger index
                ivmax = [ivmax(1:ind-1)' round(xadd(ii)) ivmax(ind:end)']'; %ivmin must be a column vector
                vmax = [vmin(1:ind-1)' vol(round(xadd(ii))) vmax(ind:end)']';
            else
                ivmax = [ivmax' round(xadd(ii))]';
                vmax = [vmax' vol(round(xadd(ii)))]';
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
clear n x y yy button;
