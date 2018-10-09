
function makeLineDraggable(hPlot, isVertical, lineData)


switch nargin 
    case 1
        isVertical = true;
        lineData = get(hPlot, 'YData');
    case 2
        if isVertical
            lineData = get(hPlot,'YData');
        else
            lineData = get(hPlot, 'XData');
        end
end
  
%Set the y data
setappdata(hPlot, 'lineData', lineData);
setappdata(hPlot, 'isVertical', isVertical)

%Set object pointer
iptPointerManager(ancestor(hPlot,'figure'));
iptSetPointerBehavior(hPlot, @(f,cp) set(f, 'Pointer', 'hand'));

%set button down fcn
set(hPlot, 'ButtonDownFcn', @startDrag);

%--------------------------------------------------------------------------
%OFFSET SUB-FUNCTIONS
%--------------------------------------------------------------------------
function startDrag(hObject, ~) %multiple inputs because it's a callback

hFig = ancestor(hObject, 'figure');
if strcmp(get(hFig, 'SelectionType'), 'normal')
    
    
    iptPointerManager(hFig, 'disable');
    set(hFig, 'Pointer','fleur');
    
    drag_motion_callback_id = iptaddcallback(hFig, ...
        'WindowButtonMotionFcn', ...
        @(hObj, edata) dragMotion(hObject)); %the line
    
    drag_up_callback_id = iptaddcallback(hFig, ...
        'WindowButtonUpFcn', @stopDrag);
    
    setappdata(hFig, 'motionFcnID', drag_motion_callback_id);
    setappdata(hFig, 'buttonUpFcnID', drag_up_callback_id);
end %end if user single clicks

function dragMotion(hLine, ~)
persistent xcoord mouse_yLoc

if ishandle(hLine)
    %get mouse location
    hAx = ancestor(hLine, 'axes');
    p = get(hAx, 'CurrentPoint');
elseif isempty(xcoord) || isempty(mouse_yLoc)
    return
else
    %     p = [xcoord + 0.01, mouse_yLoc + 0.01]; %offset so the function doesn't exit
    return
end

isVertical = getappdata(hLine, 'isVertical');

if isVertical
    xcoord = p(1);
    xdata = get(hLine, 'XData'); ydata = getappdata(hLine, 'lineData');
    [~,minI] = min(abs(xdata-xcoord));
    
    graph_yLoc = ydata(minI);
    mouse_yLoc = p(1,2);
    
    delta_y = mouse_yLoc - graph_yLoc;
    
    if delta_y == 0
        return
    end
    
    newYData = ydata+delta_y;
    ylims = get(hAx, 'YLim');
    set(hLine, 'YData', newYData);
    
    if max(newYData) > ylims(2) || min(newYData)<ylims(1)
        set(hAx, 'YLim', ylims)
    end
    
    newLoc = newYData(1);
else
    
    xcoord = p(1);
    currentxdata = get(hLine, 'XData'); xdata = getappdata(hLine, 'lineData');
    [~,minI] = min(abs(currentxdata-xcoord));
    
    graph_xLoc = xdata(minI);
    mouse_xLoc = p(1,1);
    
    delta_x = mouse_xLoc - graph_xLoc;
    
    if delta_x == 0
        return
    end
    
    newXData = xdata+delta_x;
    xlims = get(hAx, 'XLim');
    set(hLine, 'XData', newXData);
    
    if max(newXData) > xlims(2) || min(newXData)<xlims(1)
        set(hAx, 'XLim', xlims)
    end
    
	newLoc = newXData(1);
end

%if there is any TEXT associated with the line, move that too
hg = ancestor(hLine,'hggroup'); 
hgText = findobj('Parent', hg, 'Type', 'text');
if ~isempty(hgText)
    moveText(hgText, newLoc, isVertical)
end

function stopDrag(hFig, ~)
dragMotion('go');

drag_motion_callback_id = getappdata(hFig, 'motionFcnID');
drag_up_callback_id = getappdata(hFig, 'buttonUpFcnID');

iptremovecallback(hFig, 'WindowButtonMotionFcn', ...
    drag_motion_callback_id);
iptremovecallback(hFig, 'WindowButtonUpFcn', ...
    drag_up_callback_id);

set(hFig, 'Pointer','arrow')
% Enable figure's pointer manager.
iptPointerManager(hFig, 'enable')


function moveText(hgText, newLoc, isVertical)

for iHandle = hgText
    
    xyz = get(iHandle, 'Position'); x = xyz(1); y = xyz(2);
    
    if isVertical
        set(iHandle, 'Position', [x, newLoc])
    else
        set(iHandle, 'Position', [newLoc, y])
    end
    
end
