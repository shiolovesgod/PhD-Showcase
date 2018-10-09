classdef mobileLine < handle
%Last modified: 10/17/2013 to include a user defined range of mobility
    properties
        handle = [];
        hAx = [];
        hFig = [];
        thickness = [];
        callback = [];
    end
    
    properties(SetAccess = private)
        isHorizontal =  [];
        isEnabled = true;
        api = [];
    end
    
    
    methods
        
        function obj = mobileLine(varargin)
            
            [props, obj.api, obj.handle] = mobileLineAPI(varargin{:});
            
            obj.thickness = props.thickness;
            obj.hAx = props.hAx;
            obj.isHorizontal = props.isHorizontal;
            obj.hFig = props.hFig;
        end
        
        function updateVariable(obj, varName, value)
            obj.api.updateVariable(varName, value);
        end
        
        function setPosition(obj, varargin)
            
            if nargin > 1
                obj.api.setPosition(varargin{1});
            else
                obj.api.setPosition();
            end
            
        end
        
        function [startz varargout] = getPosition(obj)
            hLine = obj.handle;
            xdata = get(hLine, 'XData');
            ydata = get(hLine, 'YData'); 
            
            if obj.isHorizontal
                endpoints = [ydata(1), get(hLine,'BaseValue')];
                startz = min(endpoints);
                endz = max(endpoints);
            else
                startz = min(xdata);
                endz = max(xdata);
            end
            
            if nargout >1
                varargout{1} = endz;
            end
            
        end
        
        function set.thickness(obj,value)
            updateVariable(obj, 'thickness', value);
            obj.thickness = value;
            updateLine(obj);
            runCallback(obj);
        end
        
        
        
        function set.callback(obj,value)
            updateVariable(obj, 'callback', value)
            obj.callback = value;
        end
        
        function initializeAxis(obj)
            hAxes = obj.hAx;
            
            xlims = [0 1000]; ylims = [0 1024];
            
            set(hAxes, 'XLimMode','manual', 'YLimMode', 'manual', 'XLim', xlims, 'YLim', ylims);
            hold on
        end
        
        function updateLine(obj)
            obj.api.setPosition(); %update the line
        end
        
        function disable(obj)
            set(obj.handle,'HitTest','off');
            set(obj.handle, 'Visible','off')
            obj.isEnabled = false;
        end
        
        function enable(obj)
            set(obj.handle,'HitTest','on');
            set(obj.handle, 'Visible', 'on');
            obj.isEnabled = true;
            
            if obj.isHorizontal
%                 set(obj.hAx, 'CLim',[0 1]); %commented 10/26/2012 
            end
        end
        
        function runCallback(obj)
            if isa(obj.callback, 'function_handle')
                try
                    feval(obj.callback, obj.handle);
                catch
                    error('Error Running Callback')
                end
            end
        end
        
        function setEdgeColor(obj, color)
            if ishandle(obj.handle)
                try
                    set(obj.handle, 'EdgeColor', color)
                end
            end
        end
        
        function setFaceColor(obj, color)
            if ishandle(obj.handle)
                try
                    set(obj.handle, 'FaceColor', color)
                end
            end
        end
        
        function delete(obj)
            
            if ishandle(obj.handle)
                delete(obj.handle);
            end
            
            if isvalid(obj)
                delete(obj);
            else
                return
            end
            
            
        end
        
    end
    
    methods (Static)
        
        
    end
    
end

function [objProps api hLine] = mobileLineAPI(varargin)
%The user inputs should be: hAx, thickness, isHorizontal, limits of line

%These are the default values
hAx = gca;
thickness = 10; %z-stack
isHorizontal = true;
userLimits = [];

if nargin > 3
    hAx = varargin{1};
    thickness = varargin{2};
    isHorizontal = varargin{3};
    userLimits = varargin{4};
elseif nargin>2    
    hAx = varargin{1};
    thickness = varargin{2};
    isHorizontal = varargin{3};
elseif nargin>1
    hAx = varargin{1};
    thickness = varargin{2};
elseif nargin == 1
    hAx = varargin{1};
end

hFig = ancestor(hAx, 'figure');
xlims = get(hAx, 'Xlim'); ylims = get(hAx, 'Ylim'); clims = get(hAx, 'CLim');

%SET THE BOUDARIES AS LIMITS ON THE GRAPH

y_limits = ylims;
x_limits = xlims;

if isempty(userLimits)
    if isHorizontal
        userLimits = ylims;
    else
        userLimits = xlims;
    end
end

%DRAW THE LINE

%Set current axis first
axes(hAx);

hold on
if isHorizontal
    startPos = sum(userLimits)/2; %find the midpoint of the user defined area
    baseValue = startPos - thickness;
    hLine = area(xlims, [startPos startPos], baseValue);
    set(hAx, 'CLim',clims);
    y_limits = userLimits;
else
    startPos = sum(userLimits)/2; %find the midpoint of the user defined area
    baseValue = startPos - thickness;
    hLine = area([baseValue, baseValue, startPos, startPos], zeros(1,4)+ ylims(2), ylims(1));
    x_limits = userLimits;
end
hold off

set(hLine, 'FaceColor','none', 'EdgeColor', 'g','PickableParts','visible');

%Drawing the hLine makes the image darker, this is corrected by adj. clim

%Set Cursor action
iptPointerManager(hFig); 
iptSetPointerBehavior(hLine, @(f, cp) set(f, 'Pointer', 'fleur'));

%Set Button Down Functions
set(hLine, 'ButtonDownFcn', @startDrag);

dragFcn = [];

api.setPosition                 = @setPosition;
api.getPosition                 = @getPosition;
api.delete                      = @deleteLine;
api.setColor                    = @setColor;
api.updateVariable              = @updateVariable;
api.updateLine                  = @dragMotion;
api.isHorizontal                = isHorizontal;

setappdata(hLine, 'api', api);

%set output
objProps.thickness = thickness;
objProps.hAx = hAx;
objProps.isHorizontal= isHorizontal;
objProps.hFig = hFig;

% Initialize drag variables at function scope.
[callback ,start_x,start_y,...
    drag_motion_callback_id,drag_up_callback_id] = deal([]);



    function startDrag(varargin) %multiple inputs because it's a callback
        if strcmp(get(hFig, 'SelectionType'), 'normal')
            
            iptPointerManager(hFig, 'disable');
            
            drag_motion_callback_id = iptaddcallback(hFig, ...
                'WindowButtonMotionFcn', ...
                @dragMotion);
            
            drag_up_callback_id = iptaddcallback(hFig, ...
                'WindowButtonUpFcn', ...
                @stopDrag);
            
            dragFcn = @dragMotion;
            
        end %end if user single clicks
        
    end

    function dragMotion(varargin)
        
        if ~isempty(varargin)
            %get mouse location
            p = get(hAx, 'CurrentPoint');
        elseif isempty(start_x) || isempty(start_y)
            return
        else
            p = [start_x + 0.01, start_y + 0.01]; %offset so the function doesn't exit
        end
        
        if isHorizontal
            %original line postion
            start_y = get(hLine, 'YData'); start_y = start_y(1);
            
            new_y = p(1,2);
            
            delta_y = new_y - start_y;
            
            if delta_y == 0
                return
            end
            
            proposedVal = start_y + delta_y; %round(start_y + delta_y);
            
        else
            start_x = get(hLine, 'XData'); start_x = start_x(end);
            
            new_x = p(1,1);
            
            delta_x = new_x - start_x;
            
            if delta_x == 0
                return
            end
            
            proposedVal = start_x + delta_x;%round(start_x + delta_x);
        end
        
        setPosition(proposedVal);
        
    end

    function stopDrag(varargin)
        dragFcn();
        
        iptremovecallback(hFig, 'WindowButtonMotionFcn', ...
            drag_motion_callback_id);
        iptremovecallback(hFig, 'WindowButtonUpFcn', ...
            drag_up_callback_id);
        
        % Enable figure's pointer manager.
        iptPointerManager(hFig, 'enable');
        
        %Execute callback
        if ~isempty(callback)
            try
                feval(callback, hLine)
            end
        end
    end

    function [acceptedVal baseValue] = checkProposal(proposedVal)
        
        %         proposedVal = round(proposedVal); %make sure it's an integer
        
        if isHorizontal
            %Make sure the new y function is between the limits
            %minVal = floor(y_limits(1))+1; maxVal = floor(y_limits(2));
            minVal = y_limits(1); maxVal = y_limits(2);
        else
            %minVal = floor(x_limits(1))+1; maxVal = floor(x_limits(2));
            minVal = x_limits(1); maxVal = x_limits(2);
        end
        
        %Make sure the new y function is between the limits
        if proposedVal < minVal
            acceptedVal = minVal;
            baseValue= acceptedVal + thickness;
        elseif proposedVal + thickness > maxVal
            acceptedVal = maxVal - thickness;
            baseValue = maxVal;
        else
            acceptedVal = proposedVal;
            baseValue = acceptedVal + thickness;
        end
    end

    function setPosition(varargin)
        persistent pos bValue
        
        if nargin > 0
            proposedPosition = varargin{1};
        elseif isempty(pos) || isempty(bValue)
            return
        else
            proposedPosition = pos;
        end
        
        [pos bValue] = checkProposal(proposedPosition);
        
        
        %Make sure the right axis is the current one
        axes(hAx);
        
        if isHorizontal
            set(hLine, 'YData', [pos pos], 'BaseValue', bValue);
        else
            set(hLine, 'XData', [bValue, bValue, pos, pos])
        end
        axis tight;
        drawnow();
        
        
    end

    function updateVariable(varName, value)
        eval([varName, '= value;']);
    end


end