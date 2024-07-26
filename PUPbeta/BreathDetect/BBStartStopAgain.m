function varargout = BBStartStopAgain(varargin)
% BBStartStopAgain MATLAB code for BBStartStopAgain.fig
%      BBSTARTSTOPAGAIN, by itself, creates a new BBSTARTSTOPAGAIN or raises the existing
%      singleton*.
%
%      H = BBSTARTSTOPAGAIN returns the handle to a new BBSTARTSTOPAGAIN or the handle to
%      the existing singleton*.
%
%      BBSTARTSTOPAGAIN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BBSTARTSTOPAGAIN.M with the given input arguments.
%
%      BBSTARTSTOPAGAIN('Property','Value',...) creates a new BBSTARTSTOPAGAIN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before BBStartStopAgain_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to BBStartStopAgain_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help BBStartStopAgain

% Last Modified by GUIDE v2.5 08-May-2020 15:59:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @BBStartStopAgain_OpeningFcn, ...
                   'gui_OutputFcn',  @BBStartStopAgain_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before BBStartStopAgain is made visible.
function handles = BBStartStopAgain_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to BBStartStopAgain (see VARARGIN)

% Choose default command line output for BBStartStopAgain
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using BBStartStopAgain.
    if strcmp(get(hObject,'Visible'),'off')
        handles.GUIvars = varargin{1};
        guidata(hObject, handles); % updates handles structure

        PlotSubPlotFigure(hObject,handles);
        MakeInspStartStopTable(hObject,handles);
        
        % in case nothing gets changed, send out the same GUIvars to avoid
        % errors
        GUIvars = handles.GUIvars;
        set(0,'userdata',GUIvars);
    end

% UIWAIT makes BBStartStopAgain wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = BBStartStopAgain_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pushbutton1. ADJUST INSP START
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
GUIvars = handles.GUIvars;

time = GUIvars.time;
flow = GUIvars.Vflow;
vol = GUIvars.vol;
snore = GUIvars.Snore;
snoreDB = GUIvars.SnoreDB;
BB_i_start = GUIvars.BB_i_start;
BB_i_mid = GUIvars.BB_i_mid;
BB_i_end = GUIvars.BB_i_end;

Fs = 1/(time(2)-time(1));

% keeping this variable because the rest of the code uses it
ivmin = [BB_i_start;BB_i_end(end)];
% ivmax = BB_i_mid;

[xt,y,button] = ginput;
x = (xt-time(1))*Fs;

% make deletions
xdel = x(button==1);
for ii = 1:length(xdel)
    ind = find(ivmin>round(xdel(ii)), 1, 'first'); %returns the larger index
    ivmin(ind) = [];
%     vmin(ind) = [];
end

% make adds
xadd = x(button==3);
for ii = 1:length(xadd)
    if xadd(ii) < ivmin(end)
        ind = find(ivmin>round(xadd(ii)), 1, 'first'); %returns the larger index
        ivmin = [ivmin(1:ind-1)' round(xadd(ii)) ivmin(ind:end)']'; %ivmin must be a column vector
%         vmin = [vmin(1:ind-1)' vol(round(xadd(ii))) vmin(ind:end)']';
    else
        ivmin = [ivmin' round(xadd(ii))]';
%         vmin = [vmin' vol(round(xadd(ii)))]';
    end
end

% rename ivmin
BB_i_start = ivmin(1:end-1); % start inspiration
BB_i_end = ivmin(2:end);

% repack updated GUIvars
GUIvars = struct('time',time,'Vflow',flow,'vol',vol,...
            'BB_i_start',BB_i_start,'BB_i_mid',BB_i_mid,'BB_i_end',BB_i_end,...
            'Snore',snore,'SnoreDB',snoreDB);
        
handles.GUIvars = GUIvars; % smuggle GUIvars out in handles
guidata(hObject, handles); % update handles

%% plot again
% get axes data
axesdata = findall(gcf,'type','axes');
xlimits = get(axesdata(1),'XLim');
ylim1 = get(axesdata(1),'YLim');
ylim2 = get(axesdata(2),'YLim');
ylim3 = get(axesdata(3),'YLim');
ylim4 = get(axesdata(4),'YLim');

PlotSubPlotFigure(hObject,handles,xlimits,ylim1,ylim2,ylim3,ylim4)
MakeInspStartStopTable(hObject,handles)

% --- Executes on button press in pushbutton2. ADJUST INSP STOP
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
GUIvars = handles.GUIvars;

time = GUIvars.time;
flow = GUIvars.Vflow;
vol = GUIvars.vol;
snore = GUIvars.Snore;
snoreDB = GUIvars.SnoreDB;
BB_i_start = GUIvars.BB_i_start;
BB_i_mid = GUIvars.BB_i_mid;
BB_i_end = GUIvars.BB_i_end;

Fs = 1/(time(2)-time(1));

% keeping this variable because the rest of the code uses it
ivmax = BB_i_mid;

[xt,y,button] = ginput;
x = (xt-time(1))*Fs;

% make deletions
xdel = x(button==1);
for ii = 1:length(xdel)
    ind = find(ivmax>round(xdel(ii)), 1, 'first'); %returns the larger index
    ivmax(ind) = [];
%     vmax(ind) = [];
end

% make adds
xadd = x(button==3);
for ii = 1:length(xadd)
    if xadd(ii) < ivmax(end)
        ind = find(ivmax>round(xadd(ii)), 1, 'first'); %returns the larger index
        ivmax = [ivmax(1:ind-1)' round(xadd(ii)) ivmax(ind:end)']'; %ivmin must be a column vector
%         vmax = [vmax(1:ind-1)' vol(round(xadd(ii))) vmax(ind:end)']';
    else
        ivmax = [ivmax' round(xadd(ii))]';
%         vmax = [vmax' vol(round(xadd(ii)))]';
    end
end

% rename ivmax
BB_i_mid = ivmax;

% repack updated GUIvars
GUIvars = struct('time',time,'Vflow',flow,'vol',vol,...
            'BB_i_start',BB_i_start,'BB_i_mid',BB_i_mid,'BB_i_end',BB_i_end,...
            'Snore',snore,'SnoreDB',snoreDB);
        
handles.GUIvars = GUIvars; % smuggle GUIvars out in handles
guidata(hObject, handles); % update handles

%% plot again
% get axis limits
axesdata = findall(gcf,'type','axes');
xlimits = get(axesdata(1),'XLim');
ylim1 = get(axesdata(1),'YLim');
ylim2 = get(axesdata(2),'YLim');
ylim3 = get(axesdata(3),'YLim');
ylim4 = get(axesdata(4),'YLim');

PlotSubPlotFigure(hObject,handles,xlimits,ylim1,ylim2,ylim3,ylim4)
MakeInspStartStopTable(hObject,handles)

%[REMOVE] - only 1 minute windows
% --- Executes on button press in pushbutton4. ZOOM VIEW BUTTON
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
GUIvars = handles.GUIvars;

time = GUIvars.time;
flow = GUIvars.Vflow;
vol = GUIvars.vol;
snore = GUIvars.Snore;
snoreDB = GUIvars.SnoreDB;
BB_i_start = GUIvars.BB_i_start;
BB_i_mid = GUIvars.BB_i_mid;
BB_i_end = GUIvars.BB_i_end;

Fs = 1/(time(2)-time(1));
n=1;
l=length(flow)-1;

% get axes data % don't need this anymore because 1 minute at a time
axesdata = findall(gcf,'type','axes');
% 
% xlimits = get(axesdata(1),'XLim');
% xlim_new = [xlimits(1) xlimits(1)+7500/Fs];
% xlim_i = round(xlim_new.*Fs);
% if xlim_i(1) == 0
%     xlim_i(1) = 1;
% end

xlim_i = [1 length(time)];
xlim_new = [time(1) time(end)];
ylim1 = [min(vol(xlim_i(1):xlim_i(2)))-0.1 max(vol(xlim_i(1):xlim_i(2)))+.1];
ylim2 = [min(flow(xlim_i(1):xlim_i(2)))-0.1 max(flow(xlim_i(1):xlim_i(2)))+.1];
ylim3 = [min(snoreDB(xlim_i(1):xlim_i(2))) max(snoreDB(xlim_i(1):xlim_i(2)))+2];
ylim4 = [min(snore(xlim_i(1):xlim_i(2)))-0.1 max(snore(xlim_i(1):xlim_i(2)))+.1];

% reset axes
set(axesdata(1),'XLim',xlim_new);
set(axesdata(1),'YLim',ylim1);
set(axesdata(2),'YLim',ylim2);
set(axesdata(3),'YLim',ylim3);
set(axesdata(4),'YLim',ylim4);

% --- Executes on button press in pushbutton5. RESET VIEW BUTTON
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Load data
GUIvars = handles.GUIvars;

time = GUIvars.time;
flow = GUIvars.Vflow;
vol = GUIvars.vol;
snore = GUIvars.Snore;
snoreDB = GUIvars.SnoreDB;
BB_i_start = GUIvars.BB_i_start;
BB_i_mid = GUIvars.BB_i_mid;
BB_i_end = GUIvars.BB_i_end;

% indices
n=1;
l=length(flow)-1;

% reset axes values
xlimits = [time(n) time(n+l)];
ylim1 = [min(vol(n:n+l))-0.1 max(vol(n:n+l))+.1];
ylim2 = [min(flow(n:n+l))-0.1 max(flow(n:n+l))+.1];
ylim3 = [min(snoreDB(n:n+l)) max(snoreDB(n:n+l))+2];
ylim4 = [min(snore(n:n+l))-0.1 max(snore(n:n+l))+.1];

% get axes data
axesdata = findall(gcf,'type','axes');

% reset axes
set(axesdata(1),'XLim',xlimits);
set(axesdata(1),'YLim',ylim1);
set(axesdata(2),'YLim',ylim2);
set(axesdata(3),'YLim',ylim3);
set(axesdata(4),'YLim',ylim4);

%[REMOVED] - data is windowed one minute at a time
% --- Executes on button press in pushbutton7. BACK 1 min
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
GUIvars = handles.GUIvars;

time = GUIvars.time;
flow = GUIvars.Vflow;
vol = GUIvars.vol;
snore = GUIvars.Snore;
snoreDB = GUIvars.SnoreDB;
BB_i_start = GUIvars.BB_i_start;
BB_i_mid = GUIvars.BB_i_mid;
BB_i_end = GUIvars.BB_i_end;

Fs = 1/(time(2)-time(1));
n=1;
l=length(flow)-1;

% get axes data
axesdata = findall(gcf,'type','axes');

xlimits = get(axesdata(1),'XLim');
xlim_new = [xlimits(1)-(7500-1250)/Fs xlimits(1)+1250/Fs];
if xlim_new(1) < 0
    xlim_new(1)=0;
    xlim_new(2)=7500/Fs;
end
if xlim_new(2) > time(end)
    xlim_new(2) = time(end);
    xlim_new(1) = time(end)-7500/Fs;
end

xlim_i = round(xlim_new.*Fs);
if xlim_i(1)==0
    xlim_i(1)=1;
end

ylim1 = [min(vol(xlim_i(1):xlim_i(2)))-0.1 max(vol(xlim_i(1):xlim_i(2)))+.1];
ylim2 = [min(flow(xlim_i(1):xlim_i(2)))-0.1 max(flow(xlim_i(1):xlim_i(2)))+.1];
ylim3 = [min(snoreDB(xlim_i(1):xlim_i(2))) max(snoreDB(xlim_i(1):xlim_i(2)))+2];
ylim4 = [min(snore(xlim_i(1):xlim_i(2)))-0.1 max(snore(xlim_i(1):xlim_i(2)))+.1];

% reset axes
set(axesdata(1),'XLim',xlim_new);
set(axesdata(1),'YLim',ylim1);
set(axesdata(2),'YLim',ylim2);
set(axesdata(3),'YLim',ylim3);
set(axesdata(4),'YLim',ylim4);

%[REMOVED] - data is windowed one minute at a time
% --- Executes on button press in pushbutton8. FORWARD 1 min
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
GUIvars = handles.GUIvars;

time = GUIvars.time;
flow = GUIvars.Vflow;
vol = GUIvars.vol;
snore = GUIvars.Snore;
snoreDB = GUIvars.SnoreDB;
BB_i_start = GUIvars.BB_i_start;
BB_i_mid = GUIvars.BB_i_mid;
BB_i_end = GUIvars.BB_i_end;

Fs = 1/(time(2)-time(1));
n=1;
l=length(flow)-1;

% get axes data
axesdata = findall(gcf,'type','axes');

xlimits = get(axesdata(1),'XLim');
xlim_new = [xlimits(2)-(1250/Fs) xlimits(2)+(7500-1250)/Fs];
if xlim_new(1) < 0
    xlim_new(1)=0;
    xlim_new(2)=7500/Fs;
end
if xlim_new(2) > time(end)
    xlim_new(2) = time(end);
    xlim_new(1) = time(end)-7500/Fs;
end

xlim_i = round(xlim_new.*Fs);
if xlim_i(1)==0
    xlim_i(1)=1;
end

ylim1 = [min(vol(xlim_i(1):xlim_i(2)))-0.1 max(vol(xlim_i(1):xlim_i(2)))+.1];
ylim2 = [min(flow(xlim_i(1):xlim_i(2)))-0.1 max(flow(xlim_i(1):xlim_i(2)))+.1];
ylim3 = [min(snoreDB(xlim_i(1):xlim_i(2))) max(snoreDB(xlim_i(1):xlim_i(2)))+2];
ylim4 = [min(snore(xlim_i(1):xlim_i(2)))-0.1 max(snore(xlim_i(1):xlim_i(2)))+.1];

% reset axes
set(axesdata(1),'XLim',xlim_new);
set(axesdata(1),'YLim',ylim1);
set(axesdata(2),'YLim',ylim2);
set(axesdata(3),'YLim',ylim3);
set(axesdata(4),'YLim',ylim4);

% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file = uigetfile('*.fig');
if ~isequal(file, 0)
    open(file);
end

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
     set(hObject,'BackgroundColor','white');
end

set(hObject, 'String', {'plot(rand(5))', 'plot(sin(1:0.01:25))', 'bar(1:.5:10)', 'plot(membrane)', 'surf(peaks)'});

function PlotSubPlotFigure(hObject,handles,varargin)
GUIvars = handles.GUIvars;

time = GUIvars.time;
flow = GUIvars.Vflow;
vol = GUIvars.vol;
snore = GUIvars.Snore;
snoreDB = GUIvars.SnoreDB;
BB_i_start = GUIvars.BB_i_start;
BB_i_mid = GUIvars.BB_i_mid;
BB_i_end = GUIvars.BB_i_end;

ivmin = [BB_i_start;BB_i_end(end)];
ivmax = BB_i_mid;
n = 1;
l = length(flow)-1;

if isempty(varargin)
    xlimits = [time(n) time(n+l)];
    ylim1 = [min(vol(n:n+l))-0.1 max(vol(n:n+l))+.1];
    ylim2 = [min(flow(n:n+l))-0.1 max(flow(n:n+l))+.1];
    ylim3 = [min(snoreDB(n:n+l)) max(snoreDB(n:n+l))+2];
    ylim4 = [min(snore(n:n+l))-0.1 max(snore(n:n+l))+.1];
else
    xlimits = varargin{1};
    ylim1 = varargin{2};
    ylim2 = varargin{3}; 
    ylim3 = varargin{4}; 
    ylim4 = varargin{5}; 
end

ax1=subplot(4,1,1);
plot(time, snore)
if ~isempty(ivmin)
    hold on
    plot([time(ivmin) time(ivmin)],[min(snore(n:n+l))-0.1 max(snore(n:n+l))+.1],'r')
    plot([time(ivmax) time(ivmax)],[min(snore(n:n+l))-0.1 max(snore(n:n+l))+.1],'Color',[0.7 0.7 0.7]), hold off
end
ylabel('Snore')
xlim(xlimits)
ylim(ylim4)
ax1.Position = ax1.Position - [0.085 0 0.04 0];

ax2=subplot(4,1,2);
plot(time, snoreDB)
if ~isempty(ivmin)
    hold on
    plot([time(ivmin) time(ivmin)],[min(snoreDB(n:n+l))-0.1 max(snoreDB(n:n+l))+.1],'r')
    plot([time(ivmax) time(ivmax)],[min(snoreDB(n:n+l))-0.1 max(snoreDB(n:n+l))+.1],'Color',[0.7 0.7 0.7]), hold off
end
ylabel('SnoreDB')
xlim(xlimits)
ylim(ylim3)
ax2.Position = ax2.Position - [0.085 0 0.04 0];

nfig=2;
totalfig=4;

ax3=subplot(totalfig,1,nfig+1);
plot(time, flow)
if ~isempty(ivmin)
    hold on
    plot([time(ivmin) time(ivmin)],[min(flow(n:n+l))-0.1 max(flow(n:n+l))+.1],'r')
    plot([time(ivmax) time(ivmax)],[min(flow(n:n+l))-0.1 max(flow(n:n+l))+.1],'Color',[0.7 0.7 0.7]), hold off
end
ylabel('Flow')
xlim(xlimits)
ylim(ylim2)
ax3.Position = ax3.Position - [0.085 0 0.04 0];

ax4=subplot(totalfig,1,nfig+2);
plot(time, vol)
if ~isempty(ivmin)
    hold on
    plot([time(ivmin) time(ivmin)],[min(vol(n:n+l))-0.1 max(vol(n:n+l))+.1],'r')
    plot([time(ivmax) time(ivmax)],[min(vol(n:n+l))-0.1 max(vol(n:n+l))+.1],'Color',[0.7 0.7 0.7]), hold off
end
ylabel('Vol')
xlabel('Time')
xlim(xlimits)
ylim(ylim1)
ax4.Position = ax4.Position - [0.085 0 0.04 0];

linkaxes([ax1 ax2 ax3 ax4],'x');
set(gcf,'toolbar','figure');
set(gcf,'menubar','figure');

% repack GUIvars - I don't think this is necessary, but makes it flexibible
GUIvars = struct('time',time,'Vflow',flow,'vol',vol,...
            'BB_i_start',BB_i_start,'BB_i_mid',BB_i_mid,'BB_i_end',BB_i_end,...
            'Snore',snore,'SnoreDB',snoreDB);
        
handles.GUIvars = GUIvars; % smuggle GUIvars out in handles
guidata(hObject, handles); % update handles

% this is where GUIdata gets sent out to the main function
set(0,'userdata',GUIvars);

function MakeInspStartStopTable(hObject,handles)
% Load data into table
% unpack GUIvars
GUIvars = handles.GUIvars;
time = GUIvars.time;
flow = GUIvars.Vflow;
vol = GUIvars.vol;
snore = GUIvars.Snore;
snoreDB = GUIvars.SnoreDB;
BB_i_start = GUIvars.BB_i_start;
BB_i_mid = GUIvars.BB_i_mid;
BB_i_end = GUIvars.BB_i_end;

ivmin = [BB_i_start;BB_i_end(end)];
ivmax = BB_i_mid;
Fs = 1/(time(2)-time(1));

if length(ivmax) > length(ivmin)
    diffL = length(ivmax) - length(ivmin);
    ivmin(end+1:end+diffL) = nan;
elseif length(ivmax) < length(ivmin)
    diffL = length(ivmin) - length(ivmax);
    ivmax(end+1:end+diffL) = nan;
end

vmin_t = ivmin./Fs;
breathtimes = (ivmax-ivmin)./Fs;
dat = [round(vmin_t,1),round(ivmin,0),round(ivmax,0),round(breathtimes,1)];

set(handles.uitable1, 'Data',dat, 'ColumnFormat',{'numeric'});
set(handles.uitable1, 'ColumnName',{'Time','StartI','StopI','Ti'});
set(handles.uitable1, 'ColumnWidth',{38,40,40,48})

% repack GUIvars - I don't think this is necessary, but makes it flexibible
GUIvars = struct('time',time,'Vflow',flow,'vol',vol,...
            'BB_i_start',BB_i_start,'BB_i_mid',BB_i_mid,'BB_i_end',BB_i_end,...
            'Snore',snore,'SnoreDB',snoreDB);
        
handles.GUIvars = GUIvars; % smuggle GUIvars out in handles
guidata(hObject, handles); % update handles

% --- Executes on key press with focus on figure1 or any of its controls.
function figure1_WindowKeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

switch eventdata.Key
  case 'a'
    pushbutton7_Callback(hObject, eventdata, handles)
  case 'd'
    pushbutton8_Callback(hObject, eventdata, handles)
  ...
end

% [REMOVED]
% --- Executes on button press in pushbutton10. Calculate vmin (i.e.
% inspiratory start)
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

vol = evalin('base','vol');
ti = evalin('base','vol');
te = evalin('base','vol');
time = evalin('base','Time');
Fs = 1/(time(2)-time(1));

% Find initial peaks and troughs
% Call the function "get_vmin2" to find the troughs of the volume signal
ti = 2; te=2;
[vmin,ivmin] = get_vmin2(vol,ti,te,1.4,Fs);
% tvmin = ivmin/Fs;

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
% tvmax = ivmax/Fs;

assignin('base','vmax',vmax)
assignin('base','ivmax',ivmax)
assignin('base','vmin',vmin)
assignin('base','ivmin',ivmin)

% update plot
PlotSubPlotFigure()
MakeInspStartStopTable(handles)

%[REMOVED]
% --- Executes on button press in pushbutton11. % LOAD IVMIN AND IVMAX
function pushbutton11_Callback(hObject, eventdata, handles)
% % hObject    handle to pushbutton11 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% settings = evalin('base','settings');
% patients = evalin('base','patients');
% ptnum = evalin('base','ptnum');
% 
% % load data
% fpath = [settings.workdir,'Analyzed\'];
% loadname = [fpath, patients{ptnum,1}(1:end-4) '_A.mat'];
% load(loadname,'ivmin','ivmax')
% 
% assignin('base','ivmin',ivmin)
% assignin('base','ivmax',ivmax)
PlotSubPlotFigure(handles)
MakeInspStartStopTable(handles)
handles.GUIvars = GUIvars;

% [REMOVED]
% --- Executes on button press in pushbutton12. % Compute Breath Tables
function pushbutton12_Callback(hObject, eventdata, handles)
% close gui
closereq(); 

% % hObject    handle to pushbutton12 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% ivmin = evalin('base','ivmin');
% ivmax = evalin('base','ivmax');
% vol = evalin('base','vol');
% Time = evalin('base','Time');
% Flow = evalin('base','Flow');
% DataEventHypnog_Mat = evalin('base','DataEventHypnog_Mat');
% SnoreChannelsList = evalin('base','SnoreChannelsList');
% SnoreInterpStruct = evalin('base','SnoreInterpStruct');
% Fs = 1/(Time(2)-Time(1));
% 
% BB_i_start = ivmin(1:end-1); % start inspiration
% BB_i_mid = ivmax(1:end);% end inspiration/start expiration
% BB_i_end = ivmin(2:end);
% VTi = vol(BB_i_mid) - vol(BB_i_start);
% VTe = vol(BB_i_mid) - vol(BB_i_end);
% Ti = (BB_i_mid-BB_i_start)*(1/Fs);
% Te = (BB_i_end-BB_i_mid)*(1/Fs);
% Ttot = Ti+Te;
% VT = (VTi.*Te+VTe.*Ti)./(Ttot);
% % VTmean = sum(VT.*(Ttot))/sum(Ttot);
% % criteriafornormal = VT>(0.2*VTmean)&Ttot<10;
% % normalTtot=median(Ttot(criteriafornormal));
% VE = VT./Ttot; %force these to be the same for now until LGfromFlow is updated
% 
% %% Construct breath data table
% BB_Times(:,1) = BB_i_start; 
% BB_Times(:,2) = BB_i_mid; 
% BB_Times(:,3) = BB_i_end; 
% 
% Apnea_B = zeros(1,length(BB_Times)); % can also base this on ventilation
% BB_Ttrans = []; TiTrans = nan(length(BB_i_start),1); 
% TeTrans = nan(length(BB_i_start),1); TAA_ = [];
% Time0 = zeros(length(BB_i_start),1);
% 
% BreathDataTable=table(Time0, Time(BB_i_start), Time(BB_i_mid),...
%     Time(BB_i_end), BB_i_start, BB_i_mid, BB_i_end, VE(:));
% 
% [BreathFLDataTable] = ComputeBreathFeatures(Time, Flow, BB_Times, BB_Ttrans, TiTrans, TeTrans, TAA_, Apnea_B,[0 1 0]); %[0: downsampleHz; 1=original timing, no resample]
% 
% [BreathSnoreTable] = ComputeSnoreFeatures(Time,DataEventHypnog_Mat,SnoreInterpStruct, SnoreChannelsList,BB_Times);
% 
% % send variables to workspace
% assignin('base','BB_Times',BB_Times)
% assignin('base','BreathDataTable',BreathDataTable)
% assignin('base','BreathFLDataTable',BreathFLDataTable)
% assignin('base','BreathSnoreTable',BreathSnoreTable)

% --- Executes on button press in pushbutton13. DONE
function pushbutton13_Callback(hObject, eventdata, handles)

% close gui
closereq(); 

% % hObject    handle to pushbutton13 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% ivmin = evalin('base', 'ivmin');
% vmin = evalin('base', 'vmin');
% ivmax = evalin('base', 'ivmax');
% vmax = evalin('base', 'vmax');
% BreathDataTable = evalin('base', 'BreathDataTable');
% BreathFLDataTable = evalin('base', 'BreathFLDataTable');
% BreathSnoreTable = evalin('base', 'BreathSnoreTable');
% settings = evalin('base','settings');
% patients = evalin('base','patients');
% ptnum = evalin('base','ptnum');
% 
% % Save data
% datatosave = struct('ivmin',ivmin, 'vmin',vmin,'ivmax',ivmax,'vmax',vmax,...
%     'BreathDataTable',BreathDataTable,'BreathFLDataTable',BreathFLDataTable,...
%     'BreathSnoreTable',BreathSnoreTable);
% cd([settings.workdir,'Analyzed\']);
% savename = [patients{ptnum,1}(1:end-4) '_A.mat'];
% save(savename,'-struct','datatosave')
% disp(['Saved: ', patients{ptnum,1}(1:end-4)])
