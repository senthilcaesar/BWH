function fog3D(ah,Command,varargin)
% FOG3D adds 3D fog to the axes to emphasize depth of the 3D plot
% 
% FOG3D adds 3D fog to the axes simulated by series of semi-transparent
% planes oriented perpendicularly to the CAMERAPOSITION-CAMERATERGET line.
% The planes cover whole 3D scene and tract 3D plot rotation making an
% illusion of 3D depth fog.
% 
% Only small subset of the fog planes are visible during rotation of the 3D
% plot to make user interaction smooth and fast.
% As soon as the 3D rotation is finished (when user releases mouse button)
% all fog planes are shown.
% 
% Syntax
%       fog3D(ah,Command)
%       fog3D(ah,Command, Param1Name,Param1Value, ... )
% 
%       ah          handle of axes to add fog
% 
%       Command     string command with desired action:
%           'add'       to add fog
%                           does nothing if fog already added
%           'remove'    to remove fog completely
%                           does nothing if no fog was already added
%           'hide'      to hide fog
%                           works only if fog was already added
%           'show'      to show previously hidden fog
%                           works only if fog was already added
%           'toggle'    to toggle fog visibility
%                           works only if fog was already added
% 
%           If no COMMAND parameter is provided, 'add' is assumed.
% 
%       Available paameters are:
% 
%           Visible             logical     default: true
%                                   controls whether fog is initially visible
% 
%           FollowAxesRotation  logical     default: true
%                                   controls whether fog follows rotation of 3D plot
%                                   If FOLLOWAXESROTATION is false fog becomes
%                                   invisible when 3D plot is rotated so the
%                                   CAMERAPOSITION-CAMERATERGET line is
%                                   parallel to fog planes.
% 
%           SimpleWhenRotated   logical     default: true
%                                   controls whether fewer fog planes are
%                                   displayed when user rotates 3D plot
% 
%           NPlanes             number      default: 101
%                                   controls number of fog planes
%                                   More planes make fog visually better.
% 
%           NSimplePlanes       number      default: 11
%                                   controls how many of fog planes are visible during 3D plot rotation
% 
%           Color               1x3 vector  default: [ 1.0 1.0 1.0 ];
%                                   controls fog color
% 
%           Origin              1x3 vector  default: get(ah,'CameraTarget');
% 
% Examples
% 
%   logo;
%   fog3D(gca,'remove');
% 
%       In this case figure color is black while axes and fog color is white.
%       In order to create better image set figure color also to white.
% 
%   logo;
%   set(gcf,'Color','w')
%   fog3D(gca,'add');
% 
% 
%   Author: Andriy Nych (nych.andriy@gmail.com)
% Veersion: 1.0 (first public release)


hFig            = GetParentFigure(ah);

%% Utility functions
function SafelyDeleteAppData(h,name)
    if isappdata(h,name)
        rmappdata(h,name);
    end
end

function res = logical2OnOff(v)
    if v
        res = 'On';
    else
        res = 'Off';
    end
end

function res = fogAlreadyAdded(ah)
    res = isappdata(ah,'fog3d_props') | isappdata(ah,'hR3D');
end

function RotateFogToNewView(FogProps,newView)
    for k=1:FogProps.NPlanes
        set( FogProps.hFog(k), ...
            'XData',FogProps.XData0{k}, ...
            'YData',FogProps.YData0{k}, ...
            'ZData',FogProps.ZData0{k} );
        rotate( FogProps.hFog(k), [1 0 0],-newView(2)+00, FogProps.Origin );
        rotate( FogProps.hFog(k), [0 0 1],+newView(1)+00, FogProps.Origin );
    end
end

if ~isappdata(hFig,'fogAlreadyAdded')
    setappdata(hFig,'fogAlreadyAdded',@fogAlreadyAdded);
end
if ~isappdata(hFig,'RotateFogToNewView')
    setappdata(hFig,'RotateFogToNewView',@RotateFogToNewView);
end

%% Functions
function FogProps = CreateFogPlanes(ah,FogProps)
        XL_             = xlim(ah);
        YL_             = ylim(ah);
        ZL_             = zlim(ah);
        FogPlaneSize    = max( abs( [ XL_ YL_ ZL_] ) ) * 2; %4 % fog planes should be big enough cover whole scene
        axView          = get(ah,'View');
        FogProps.hFog   = zeros(FogProps.NPlanes,1);
        xx              = linspace( -FogProps.Thickness/2 , +FogProps.Thickness/2 , FogProps.NPlanes );
        for k=1:FogProps.NPlanes
            cc              = FogProps.Color(1,1:3);
            if size(FogProps.Color,1)==FogProps.NPlanes
            cc              = FogProps.Color(k,1:3);
            end
            FogProps.hFog(k) = patch( ...
                FogProps.Origin(1)+[+1 +1 +1 +1]*xx(k), ...
                FogProps.Origin(2)+[-1 -1 +1 +1]*FogPlaneSize/2, ...
                FogProps.Origin(3)+[-1 +1 +1 -1]*FogPlaneSize/2, ...
                cat(3, ...
                    ones(1,4)*1.00, ...
                    ones(1,4)*1.00, ...
                    ones(1,4)*1.00 ), ...
                'AmbientStrength', 1.0, ...
                'FaceLighting',     'none', ...
                'EdgeLighting',     'none', ...
                'FaceColor',        cc, ...
                'EdgeColor',        FogProps.EdgeColor, ...
                'FaceAlpha',        FogProps.PlaneAlpha, ...
                'EdgeAlpha',        1.0, ...
                'Parent',           ah ...
                );
            rotate( FogProps.hFog(k), [0 0 1],90);
        end
        FogProps.hFogSimple = FogProps.hFog( unique(round(linspace(1,FogProps.NPlanes,FogProps.NSimplePlanes))) );
        FogProps.XData0     = cell(FogProps.NPlanes,1);
        FogProps.YData0     = cell(FogProps.NPlanes,1);
        FogProps.ZData0     = cell(FogProps.NPlanes,1);
        for k=1:FogProps.NPlanes
        FogProps.XData0{k}  = get(FogProps.hFog(k),'XData');
        FogProps.YData0{k}  = get(FogProps.hFog(k),'YData');
        FogProps.ZData0{k}  = get(FogProps.hFog(k),'ZData');
        rotate( FogProps.hFog(k), [1 0 0],-axView(2)+00, FogProps.Origin );
        rotate( FogProps.hFog(k), [0 0 1],+axView(1)+00, FogProps.Origin );
        end
end

if nargin<2
    Command         = 'Add';
end

switch lower(Command)

    case 'toggle'
        if fogAlreadyAdded(ah)
            FogProps = getappdata(ah,'fog3d_props');
            FogProps.Visible = ~FogProps.Visible;
            set(FogProps.hFog, 'Visible',logical2OnOff(FogProps.Visible));
            setappdata(ah,'fog3d_props',FogProps);
        else
            return;
        end

    case 'hide'
        if fogAlreadyAdded(ah)
            FogProps = getappdata(ah,'fog3d_props');
            if FogProps.Visible
            FogProps.Visible = false;
            else
            end
            set(FogProps.hFog, 'Visible',logical2OnOff(FogProps.Visible));
            setappdata(ah,'fog3d_props',FogProps);
        else
            return;
        end

    case 'show'
        if fogAlreadyAdded(ah)
            FogProps = getappdata(ah,'fog3d_props');
            if FogProps.Visible
            else
            FogProps.Visible = true;
            end
            set(FogProps.hFog, 'Visible',logical2OnOff(FogProps.Visible));
            setappdata(ah,'fog3d_props',FogProps);
        else
            return;
        end

    case 'remove'
        if fogAlreadyAdded(ah)
            FogProps = getappdata(ah,'fog3d_props');
            delete(FogProps.hFog);
            SafelyDeleteAppData(ah,   'hR3D'           );
            SafelyDeleteAppData(ah,   'fog3d_props'    );
            SafelyDeleteAppData(ah,   'fog3d_rotating' );
            SafelyDeleteAppData(hFig, 'ActiveFogAxes'  );
        else
            return;
        end

    case 'add'
        if fogAlreadyAdded(ah)
            return;
        else
            XL                          = xlim(ah);
            YL                          = ylim(ah);
            ZL                          = zlim(ah);
            FogProps.Visible            = true;
            FogProps.FollowAxesRotation = false;
            FogProps.SimpleWhenRotated  = true;
            FogProps.NPlanes            = 101;
            FogProps.NSimplePlanes      = 11;
            FogProps.Thickness          = mean( [ diff(XL) diff(YL) diff(ZL) ] );
            FogProps.PlaneAlpha         = 1-nthroot(0.5,FogProps.NPlanes);       %%% fog strength
            FogProps.Color              = [ 1.0 1.0 1.0 ] * 1.0;
            FogProps.Origin             = get(ah,'CameraTarget');

            FogProps                    = parseargs_modified(FogProps,varargin{:});

            FogProps.EdgeColor          = [ 1.0 1.0 1.0 ] * 0.5;
            FogProps.EdgeColor          = 'none';
            FogProps                    = CreateFogPlanes(ah,FogProps);
            hR3D                        = rotate3d(ah);
            set(hR3D,'RotateStyle','orbit','Enable','on');
            setappdata(ah,   'hR3D',           hR3D);
            setappdata(ah,   'fog3d_props',    FogProps);
            setappdata(ah,   'fog3d_rotating', false );
            setappdata(hFig, 'ActiveFogAxes',  []    );
            set(hR3D,'ActionPreCallback',@myprecallback);
            set(hR3D,'ActionPostCallback',@mypostcallback);
            set(hFig, 'WindowButtonMotionFcn',@fig_WindowButtonMotionFcn_callback );
            set(hR3D,'Enable','off');
        end
end

end

function myprecallback(obj,evd) %#ok

if isappdata(gcbf,'fogAlreadyAdded')
    fogAlreadyAdded_fcn = getappdata(gcbf,'fogAlreadyAdded');
end

if fogAlreadyAdded_fcn(evd.Axes)
    setappdata(gcbf, 'ActiveFogAxes',evd.Axes );
    FogProps = getappdata(evd.Axes,'fog3d_props');
    setappdata(gca, 'fog3d_rotating',true );
    if FogProps.FollowAxesRotation
        if FogProps.SimpleWhenRotated
            set( FogProps.hFog,       'Visible','off');
            set( FogProps.hFogSimple, 'Visible','on');
            % we have to adjust transparency of the fog planes in
            % interactive mode according to number of planes shown
            set( FogProps.hFog,       'FaceAlpha',FogProps.PlaneAlpha*length(FogProps.hFog)/length(FogProps.hFogSimple) );
        end  % if FogProps.SimpleWhenRotated
    end  % if FogProps.FollowAxesRotation
end  % if fogAlreadyAdded

end

function mypostcallback(obj,evd) %#ok

if isappdata(gcbf,'fogAlreadyAdded')
    fogAlreadyAdded_fcn = getappdata(gcbf,'fogAlreadyAdded');
end
if isappdata(gcbf,'RotateFogToNewView')
    RotateFogToNewView_fcn = getappdata(gcbf,'RotateFogToNewView');
end

if fogAlreadyAdded_fcn(evd.Axes)
    setappdata(gcbf, 'ActiveFogAxes',[] );
    FogProps = getappdata(evd.Axes,'fog3d_props');
    setappdata( gca, 'fog3d_rotating',false );
    if FogProps.FollowAxesRotation
        newView = get(evd.Axes,'View');
        RotateFogToNewView_fcn(FogProps,newView);
    end  % if FogProps.FollowAxesRotation
    if FogProps.SimpleWhenRotated
        % now we have to restore adjusted transparency of the fog planes
        % interactive mode
        set( FogProps.hFog,       'Visible','on');
        set( FogProps.hFog,       'FaceAlpha',FogProps.PlaneAlpha);
    end  % if FogProps.SimpleWhenRotated
end  % if fogAlreadyAdded(evd.Axes)

end

function fig_WindowButtonMotionFcn_callback(obj,evd) %#ok

if isappdata(gcbf,'fogAlreadyAdded')
    fogAlreadyAdded_fcn = getappdata(gcbf,'fogAlreadyAdded');
end
if isappdata(gcbf,'RotateFogToNewView')
    RotateFogToNewView_fcn = getappdata(gcbf,'RotateFogToNewView');
end

ActiveFogAxes = getappdata(gcbf, 'ActiveFogAxes' );
if isempty(ActiveFogAxes)
    return;
end

if fogAlreadyAdded_fcn(ActiveFogAxes)
    if getappdata(gca,'fog3d_rotating')
        newView = round(get(ActiveFogAxes,'View'));
        FogProps = getappdata(gca,'fog3d_props');
        if FogProps.FollowAxesRotation
            RotateFogToNewView_fcn(FogProps,newView);
        end  % if FogProps.FollowAxesRotation
    end
end
end

function h = GetParentFigure(h)
while ~strcmp(get(h,'Type'),'figure')
    h = get(h,'Parent');
end
end

function X = parseargs_modified(X,varargin)
%PARSEARGS_MODIFIED - Parses name-value pairs
%
% PARSEARGS original:   Behaves like setfield, but accepts multiple name-value pairs and provides
% PARSEARGS original:   some additional features:
% PARSEARGS original:   1) If any field of X is an cell-array of strings, it can only be set to
% PARSEARGS original:      one of those strings.  If no value is specified for that field, the
% PARSEARGS original:      first string is selected.
% PARSEARGS original:   2) Where the field is not empty, its data type cannot be changed
% PARSEARGS original:   3) Where the field contains a scalar, its size cannot be changed.
%
% code was taken from PARSEARGS and stripped a little to remove feature #1

remaining = nargin-1; % number of arguments other than X
count = 1;
fields = fieldnames(X);
modified = zeros(size(fields));
% Take input arguments two at a time until we run out.
while remaining>=2
    fieldname = varargin{count};
    fieldind = find(strcmp(fieldname,fields));
    if ~isempty(fieldind)
        oldvalue = getfield(X,fieldname); %#ok
        newvalue = varargin{count+1};
        if ~isempty(oldvalue)
            % The caller isn't allowed to change the data type of a non-empty property,
            % and scalars must remain as scalars.
            if ~strcmp(class(oldvalue),class(newvalue))
                error(sprintf('Cannot change class of field "%s" from "%s" to "%s"',...
                    fieldname,class(oldvalue),class(newvalue))); %#ok
            elseif numel(oldvalue)==1 & numel(newvalue)~=1 %#ok
                error(sprintf('New value for "%s" must be a scalar',fieldname));  %#ok
            end
        end
        X = setfield(X,fieldname,newvalue); %#ok
        modified(fieldind) = 1;
    else
        error(['Not a valid field name: ' fieldname]);
    end
    remaining = remaining - 2;
    count = count + 2;
end
% Check that we had a value for every name.
if remaining~=0
    error('Odd number of arguments supplied.  Name-value pairs required');
end
end
