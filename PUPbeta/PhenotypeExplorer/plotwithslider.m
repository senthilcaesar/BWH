function plotwithslider(timeminmax)
date2PUPtime = @(x) (datenum(x)-693961)*86400;
PUPtime2date = @(x) datetime(x/86400,'ConvertFrom','excel');

global tminmax P h1 valuelast xlimlast sliderstep xvalues yvalues range ax2 xvalues_rm yvalues_rm
sliderstep = [0.025 0.8]; % [0.25 1]
tminmax = timeminmax;
%range = 2; % visible x range
if 1
    P = [tminmax(1)+range/2 tminmax(2)-range/2]; % min and max mean(xlim Positions)
else
    P = [tminmax(1) tminmax(2)]; % min and max mean(xlim Positions)
end

   if isa(tminmax,'datetime')
       if diff(P)<0, P=PUPtime2date((date2PUPtime(tminmax(1))+date2PUPtime(tminmax(2)))/2*[1 1]+[-0.5 +0.5]); end
   else
       if diff(P)<0, P=(tminmax(1)+tminmax(2))/2*[1 1]+[-0.5 +0.5]; end
   end


axis tight; xlim([tminmax(1) tminmax(1)+range])
%plot(x,y);

set(gca,'ylimmode','auto');
set(gca,'ylimmode','manual');
h1=uicontrol('style','slider','units','normalized','position',[0 0 1 0.05],'KeyPressFcn',@key_pressed_fcn,...
    'callback',@slider1);
valuelast=get(h1,'Value');
xlimlast=get(gca,'xlim');
Fsstep = sliderstep*range/(P(2)-P(1));
Fsstep(Fsstep>1)=1;
set(h1,'sliderstep',Fsstep);
uicontrol(...
    'Units','normalized',...
    'Position',[0.01,0.99-0.025,0.15,0.025],...
    'String','Select',...
    'KeyPressFcn', @key_pressed_fcn,...
    'Callback',@selectdata);
uicontrol(...
    'Units','normalized',...
    'Position',[0.17,0.99-0.025,0.15,0.025],...
    'String','Remove',...
    'KeyPressFcn', @key_pressed_fcn,...
    'Callback',@removedata);
uicontrol(...
    'Units','normalized',...
    'Position',[0.17+0.16,0.99-0.025,0.15,0.025],...
    'String','ClrLastSelect',...
    'KeyPressFcn', @key_pressed_fcn,...
    'Callback',@ClrLastSelect);
uicontrol(...
    'Units','normalized',...
    'Position',[0.17+2*0.16,0.99-0.025,0.15,0.025],...
    'String','ClrLastRemove',...
    'KeyPressFcn', @key_pressed_fcn,...
    'Callback',@ClrLastRemove);
% uicontrol(...
%     'Units','normalized',...
%     'Position',[0.17+0.16,0.99-0.025,0.15,0.025],...
%     'String','ZoomInX',...
%     'KeyPressFcn', @key_pressed_fcn,...
%     'Callback',@zoominX);
% uicontrol(...
%     'Units','normalized',...
%     'Position',[0.17+2*0.16,0.99-0.025,0.15,0.025],...
%     'String','ZoomOutX',...
%     'KeyPressFcn', @key_pressed_fcn,...
%     'Callback',@zoomoutX);

set(gcf,'KeyPressFcn', @key_pressed_fcn);
end

function slider1(hObj,event)
date2PUPtime = @(x) (datenum(x)-693961)*86400;
PUPtime2date = @(x) datetime(x/86400,'ConvertFrom','excel');
global P valuelast xlimlast sliderstep range tminmax
% called when slider was used
global ax2
for i=1:length(ax2)
    set(ax2(i),'ylimmode','manual');
end
valtemp=get(hObj,'value');
attemptedwindowstepnorm=(valtemp-valuelast)*(P(2)-P(1))/range;

xlimnew=get(gca,'xlim');
range=xlimnew(2)-xlimnew(1);
P = [tminmax(1)+range/2 tminmax(2)-range/2]; % min and max mean(xlim Positions)
if isa(tminmax,'datetime')
    if diff(P)<0, P=PUPtime2date((date2PUPtime(tminmax(1))+date2PUPtime(tminmax(2)))/2*[1 1]+[-0.5 +0.5]); end
else
    if diff(P)<0, P=(tminmax(1)+tminmax(2))/2*[1 1]+[-0.5 +0.5]; end
end
Fsstep = sliderstep*range/(P(2)-P(1));
Fsstep(Fsstep>1)=1;
set(hObj,'sliderstep',Fsstep);

if sum(xlimnew==xlimlast)==0 %xrange moved not via slider
    mxl = xlimnew(1)+range/2;
    val = (mxl-P(1))/(P(2)-P(1)); %first move slider to appropriate location
    %set(hObj,'value',val); %seems needless
    %         switch round(1000*attemptedwindowstepnorm)/1000 %minor rounding issues
    %             case sliderstep(2) %large step forward
    %             val=val+range*sliderstep(2)/(P(2)-P(1));
    %             case sliderstep(1) %small step forward
    %             val=val+range*sliderstep(1)/(P(2)-P(1));
    %             case -sliderstep(2) %large step backward
    %             val=val-range*sliderstep(2)/(P(2)-P(1));
    %             case -sliderstep(1) %small step backward
    %             val=val-range*sliderstep(1)/(P(2)-P(1));
    %         end
    val=val+range*attemptedwindowstepnorm/(P(2)-P(1));
    set(hObj,'value',val);
else %move x limits simply according to slider
    val = get(hObj,'value');
end
mxl = val*(P(2)-P(1))+P(1); % mean xlim, denormalize
xlim([mxl-range/2 mxl+range/2]);
xlimlast=get(gca,'xlim');
valuelast=val;
end

function selectdata(varargin)
global xvalues yvalues
[tempx,tempy] = ginput(2);
xvalues = [xvalues;tempx'];
yvalues = [yvalues;tempy'];
outputlive = [xvalues yvalues]
hold(gca,'on');
plot([xvalues(end,1) xvalues(end,2) xvalues(end,2) xvalues(end,1) xvalues(end,1)],[yvalues(end,1) yvalues(end,1) yvalues(end,2) yvalues(end,2) yvalues(end,1)],'color',[0 1 0]);
hold(gca,'off');
end %inputCallback

function ClrLastSelect(varargin)
global xvalues yvalues
hold(gca,'on');
plot([xvalues(end,1) xvalues(end,2) xvalues(end,2) xvalues(end,1) xvalues(end,1)],[yvalues(end,1) yvalues(end,1) yvalues(end,2) yvalues(end,2) yvalues(end,1)],'color',[1 0.8 0.8]);
hold(gca,'off');
xvalues(end,:)=[];
yvalues(end,:)=[];
outputlive = [xvalues yvalues]
end %inputCallback


function removedata(varargin)
global xvalues_rm yvalues_rm
[tempx_rm,tempy_rm] = ginput(2);
xvalues_rm = [xvalues_rm;tempx_rm'];
yvalues_rm = [yvalues_rm;tempy_rm'];
% outputlive = [xvalues_rm yvalues_rm]
hold(gca,'on');
plot([xvalues_rm(end,1) xvalues_rm(end,2) xvalues_rm(end,2) xvalues_rm(end,1) xvalues_rm(end,1)],[yvalues_rm(end,1) yvalues_rm(end,1) yvalues_rm(end,2) yvalues_rm(end,2) yvalues_rm(end,1)],'color',[1 0 0]);
hold(gca,'off');
end %inputCallback

function ClrLastRemove(varargin)
global xvalues_rm yvalues_rm
hold(gca,'on');
plot([xvalues_rm(end,1) xvalues_rm(end,2) xvalues_rm(end,2) xvalues_rm(end,1) xvalues_rm(end,1)],[yvalues_rm(end,1) yvalues_rm(end,1) yvalues_rm(end,2) yvalues_rm(end,2) yvalues_rm(end,1)],'color',[1 0.8 0.8]);
hold(gca,'off');
xvalues_rm(end,:)=[];
yvalues_rm(end,:)=[];
% outputlive = [xvalues_rm yvalues_rm]
end %inputCallback

function zoominX(varargin)
date2PUPtime = @(x) (datenum(x)-693961)*86400;
PUPtime2date = @(x) datetime(x/86400,'ConvertFrom','excel');
global range h1 sliderstep P tminmax
xlimnew=get(gca,'xlim');
xlimdelta=xlimnew(2)-xlimnew(1);
xlimnew(2)=xlimnew(1)+xlimdelta/2;
set(gca,'xlim',xlimnew);

range=xlimnew(2)-xlimnew(1);
P = [tminmax(1)+range/2 tminmax(2)-range/2]; % min and max mean(xlim Positions)
if isa(tminmax,'datetime')
    if diff(P)<0, P=PUPtime2date((date2PUPtime(tminmax(1))+date2PUPtime(tminmax(2)))/2*[1 1]+[-0.5 +0.5]); end
else
    if diff(P)<0, P=(tminmax(1)+tminmax(2))/2*[1 1]+[-0.5 +0.5]; end
end
Fsstep = sliderstep*range/(P(2)-P(1));
Fsstep(Fsstep>1)=1;
set(h1,'sliderstep',Fsstep);

mxl = xlimnew(1)+range/2;
val = (mxl-P(1))/(P(2)-P(1)); %first move slider to appropriate location
val(val<0)=0; %was needed, bug fix
set(h1,'value',val);
valuelast=val; %unused
end %inputCallback

function zoomoutX(varargin)
date2PUPtime = @(x) (datenum(x)-693961)*86400;
PUPtime2date = @(x) datetime(x/86400,'ConvertFrom','excel');
global range h1 sliderstep P valuelast tminmax
xlimnew=get(gca,'xlim');
xlimdelta=xlimnew(2)-xlimnew(1);
xlimnew(2)=xlimnew(1)+xlimdelta*2;
set(gca,'xlim',xlimnew);

range=xlimnew(2)-xlimnew(1);
P = [tminmax(1)+range/2 tminmax(2)-range/2]; % min and max mean(xlim Positions)
if isa(tminmax,'datetime')
    if diff(P)<0, P=PUPtime2date((date2PUPtime(tminmax(1))+date2PUPtime(tminmax(2)))/2*[1 1]+[-0.5 +0.5]); end
else
    if diff(P)<0, P=(tminmax(1)+tminmax(2))/2*[1 1]+[-0.5 +0.5]; end
end
Fsstep = sliderstep*range/(P(2)-P(1));
Fsstep(Fsstep>1)=1;
set(h1,'sliderstep',Fsstep);

mxl = xlimnew(1)+range/2;
val = (mxl-P(1))/(P(2)-P(1)); %first move slider to appropriate location
val(val<0)=0; %if needed
set(h1,'value',val);
valuelast=val;
end %inputCallback


function optimize(varargin)
global ax2
for i=1:length(ax2)
    %set(ax2(i),'ylimmode','auto');
    try
        xlims = get(ax2(i),'xlim');
        h = get(ax2(i), 'Children');
        clear ylims
        for j=1:length(h)
            temp = (h(j).XData > xlims(1) & h(j).XData < xlims(2));
            ylims(j,1) = min(h(j).YData(temp));
            ylims(j,2) = max(h(j).YData(temp));
        end
        ylims = [min(ylims(:,1)) max(ylims(:,2))];
        set(ax2(i),'ylim',ylims);
    catch
    end
end
% figure(gcf);
% set(gca,'ylimmode','manual');
end %inputCallback


function key_pressed_fcn(varargin)
global h1 range sliderstep P valuelast tminmax
switch varargin{2}.Key
    case 'c'
        zoomoutX();
        slider1(h1);
    case 'x'
        zoominX();
        slider1(h1);
    case 's'
        selectdata();
    case 'z'
        optimize();
    case 'rightarrow'
        %'rightarrow'
        val=get(h1,'value');
        val=val+range*sliderstep(1)/(P(2)-P(1));
        val(val<0)=0; %if needed
        set(h1,'value',val);
        slider1(h1);
        valuelast=val;
    case 'leftarrow'
        %'rightarrow'
        val=get(h1,'value');
        val=val-range*sliderstep(1)/(P(2)-P(1));
        val(val<0)=0; %if needed
        set(h1,'value',val);
        slider1(h1);
        valuelast=val;
    case 'pageup'
        %'rightarrow'
        val=get(h1,'value');
        val=val-range*sliderstep(2)/(P(2)-P(1));
        val(val<0)=0; %if needed
        set(h1,'value',val);
        slider1(h1);
        valuelast=val;
    case 'pagedown'
        %'rightarrow'
        val=get(h1,'value');
        val=val+range*sliderstep(2)/(P(2)-P(1));
        val(val<0)=0; %if needed
        set(h1,'value',val);
        slider1(h1);
        valuelast=val;
    case 'q' %pageup
        %'rightarrow'
        val=get(h1,'value');
        val=val-range*sliderstep(2)/(P(2)-P(1));
        val(val<0)=0; %if needed
        set(h1,'value',val);
        slider1(h1);
        valuelast=val;
    case 'w' %pagedown
        %'rightarrow'
        val=get(h1,'value');
        val=val+range*sliderstep(2)/(P(2)-P(1));
        val(val<0)=0; %if needed
        set(h1,'value',val);
        slider1(h1);
        valuelast=val;
    case '3' %set 30s page width, anchored at left
        delta = 15;
        if isdatetime(tminmax)
            delta = seconds(15);
        end
        xlims=get(gca,'xlim');
        xlimnew(1) = mean(xlims)-delta;
        xlimnew(2) = mean(xlims)+delta;
        set(gca,'xlim',xlimnew);
        range=xlimnew(2)-xlimnew(1);
        P = [tminmax(1)+range/2 tminmax(2)-range/2]; % min and max mean(xlim Positions)
        if diff(P)<0, P=(tminmax(1)+tminmax(2))/2*[1 1]+[-0.5 +0.5]; end
        %set(h1,'sliderstep',sliderstep*range/(P(2)-P(1)));
        
        Fsstep = sliderstep*range/(P(2)-P(1));
        Fsstep(Fsstep>1)=1;
        set(h1,'sliderstep',Fsstep);
        
        mxl = xlimnew(1)+range/2;
        val = (mxl-P(1))/(P(2)-P(1)); %first move slider to appropriate location
        val(val<0)=0; %if needed
        set(h1,'value',val);
        valuelast=val;
        slider1(h1);
end
end