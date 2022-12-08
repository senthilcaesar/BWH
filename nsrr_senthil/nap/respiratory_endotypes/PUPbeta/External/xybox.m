function [xmin2,xmax2,ymin2,ymax2] = xybox(xmin,xmax,ymin,ymax,m,c)

temp = [xmin xmax].*m+c;
temp2 = ([ymin ymax]-c)/m;

if m<0
    if temp(2)<ymin
        xmax2 = (ymin-c)/m;
    else
        xmax2 = xmax;
    end
    if temp(1)>ymax
        xmin2 = (ymax-c)/m;
    else
        xmin2 = xmin;
    end
elseif m>0
    if temp(1)<ymin
        xmin2 = (ymin-c)/m;
    else
        xmin2 = xmin;
    end
    if temp(1)>ymax
        xmax2 = (ymax-c)/m;
    else
        xmax2 = xmax;
    end
end
ymin2 = min([xmin2 xmax2].*m+c);
ymax2 = max([xmin2 xmax2].*m+c);

% figure(1)
% plot([xmin,xmax],[ymin,ymax],':');
% hold('on')
% plot([xmin,xmax],[xmin,xmax]*m+c,'r');
% plot([xmin2,xmax2],[xmin2,xmax2]*m+c,'k');