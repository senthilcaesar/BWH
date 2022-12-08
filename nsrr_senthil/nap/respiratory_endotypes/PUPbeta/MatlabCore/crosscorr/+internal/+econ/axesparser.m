function [ax,args] = axesparser(varargin)

% AXESPARSER extracts an (optional) axes handle ax from the initial
% position in an input argument list, and checks for validity. Use
% AXESPARSER as a preprocessor for plotting functions, before input
% parsing on args. If isempty(ax), plot to gca.

ax = gobjects(0);
args = varargin;

if ~isempty(args)
    
    arg1 = args{1};
    
    if isa(arg1,'matlab.graphics.Graphics')

        if isgraphics(arg1,'axes')

            args = args(2:end);
            ax = arg1;

        else

            error(message('econ:internal:econ:axesparser:InvalidParent'));

        end

    end

end