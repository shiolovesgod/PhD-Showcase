function varargout = manualBiometry(varargin)
% MANUALBIOMETRY MATLAB code for manualBiometry.fig
%      MANUALBIOMETRY, by itself, creates a new MANUALBIOMETRY or raises the existing
%      singleton*.
%
%      H = MANUALBIOMETRY returns the handle to a new MANUALBIOMETRY or the handle to
%      the existing singleton*.
%
%      MANUALBIOMETRY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MANUALBIOMETRY.M with the given input arguments.
%
%      MANUALBIOMETRY('Property','Value',...) creates a new MANUALBIOMETRY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before manualBiometry_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to manualBiometry_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help manualBiometry

% Last Modified by GUIDE v2.5 07-Jun-2016 13:32:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @manualBiometry_OpeningFcn, ...
                   'gui_OutputFcn',  @manualBiometry_OutputFcn, ...
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


% --- Executes just before manualBiometry is made visible.
function manualBiometry_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to manualBiometry (see VARARGIN)

% Choose default command line output for manualBiometry
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


if isempty(varargin)
    imgInfo = loadImages;
else
    imgInfo = varargin{1};
end

if isempty(imgInfo)
    return
end

if ~isfield(imgInfo,'surfLoc') || isempty(imgInfo.surfLoc)
    imgInfo.surfLoc = estimateSurfaceLocation(imgInfo);
end

imgInfo.tol = 0; %set tolerance for center A-Line
hDummyFig = figure; colormap (hObject, gray); delete(hDummyFig);

%Initialize Plot
initializeFigurePlots(handles, imgInfo);
setappdata(hObject, 'imgInfo', imgInfo);


%Initialize plot tools
hToggles = findobj(handles.panel_toolbar, 'Style', 'togglebutton');
defaultImgs = get(hToggles,'UserData');

%convert images to true color
trueColorImgs = defaultImgs;
for i = 1:numel(hToggles)
    for j = 1:2 
        jImg = defaultImgs{i}{j};
        
        if size(jImg, 3) < 3
            jImg = repmat(jImg,[1,1,3]);
        end
        
        if max(jImg(:)) == 1
            jImg = uint8(jImg*255);
        elseif ~strcmpi(class(jImg),'uint8')
            jImg = uint8(jImg);
        end 
        trueColorImgs{i}{j} = jImg; 
    end
end

arrayfun(@(iToggle) iptaddcallback(iToggle,'Callback',...
    @(hObj, edata) setToggleImgs(hToggles, trueColorImgs)), hToggles);



% UIWAIT makes manualBiometry wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = manualBiometry_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function imgInfo = loadImages

%Initialize output
imgInfo = [];

%Refractive index
n_default = [1,...air
    1.387,...cornea
    1.342,...aqueous
    1.415,...lens
    1.341,...vitreous
    1.380]; %retina
n = n_default(2:end);

imgInfo.n = n;


%Let the user select three images (or two)
[fNames, nImgsPerScan] = uiloadimgdir;

if isempty(fNames)
    return
end

fNames = fNames(1:nImgsPerScan);


sigma = 6; kSize = sigma*2;
H = fspecial('gaussian',kSize, sigma); %img Filter
img1 = imread(fNames{1}); imgSize = size(img1);
imgsIn = zeros([size(img1),3]);

for i = 1:nImgsPerScan
    iRaw = imread(fNames{i});
    iFiltered = imfilter(iRaw,H);
    imgsIn(:,:,i) = iFiltered;
end

stitchOutput = stitchMarcoOCT(imgsIn);
fullImg = stitchOutput.stitchedImg;

%CORNEA AND LENS
img2End = stitchOutput.absoluteImgEnd(2);
img3Range = stitchOutput.absoluteStitchStart(3,:);
topImg = fullImg(1:img2End, :);
bottomImg = fullImg(img3Range(1):img3Range(2),:,1);

imgInfo.topImg = topImg;
imgInfo.bottomImg = bottomImg;
imgInfo.imgF = imgsIn;
imgInfo.nImgsPerScan = nImgsPerScan;

%Estimate center a line and coords
aLineSum = sum(imgsIn(:,:,1),1);
midX = round(imgSize(2)/2); halfTol = round(imgSize(2)/8);
midALine = aLineSum; midALine([1:midX-halfTol,halfTol+midX:end]) = 0;
[~,cntrLineIdx] = max(midALine);

yRes = imgSize(1)/10.43;

imgInfo.cntrIdx = cntrLineIdx;
imgInfo.stitchOutput = stitchOutput;
imgInfo.yRes = yRes;
[~, imgInfo.ID] = fileparts(fNames{1});

allDistances = estimateSurfaceLocation(imgInfo);

imgInfo.surfLoc = allDistances;

function allDistances = estimateSurfaceLocation(imgInfo)
%This function guesses initial values for each surface.
%INPUTS
    %imgsIn --> (m x n x nImgs) filtered images
    %stitchOutput --> structure containing: absoluteImgStart &  frameShift
        %fields
    %yRes --> yResolution of image

imgsIn = imgInfo.imgF;
cntrLineIdx = imgInfo.cntrIdx;
stitchOutput = imgInfo.stitchOutput;
yRes = imgInfo.yRes;
nImgsPerScan = imgInfo.nImgsPerScan;

%Guess-timate every other surface by average distances
%CCT ~= 0.8; ACD ~=3.824; LT ~=5.349; RT ~=0.1

cntrAScan1 = smooth(double(imgsIn(:,cntrLineIdx,1)),15);
cntrAScan1(1:100)= 0;

mpd1 = round(0.5*yRes);
corneaIdx = findALinePeaks(cntrAScan1, mpd1);
allDistances = round(force1D(corneaIdx, 2)); 


if nImgsPerScan > 1
    cntrAScan2 = smooth(double(imgsIn(:,cntrLineIdx,2)),15);
    cntrAScan2(1:50)= 0;
    
    mpd2 = round(4*yRes);
    lensIdx = round(force1D(findALinePeaks(cntrAScan2, mpd2),2));

    absImgStart = stitchOutput.absoluteImgStart;
    frameShift = absImgStart - 1;
    if stitchOutput.isFlipImg(2)
        lensIdx = flipud(lensIdx);
        shift2 =  -(size(imgsIn, 1) + frameShift(2));
    else
        shift2 = frameShift(2);
    end
    
    %Shift it in accordance with stitch

    lensIdx = abs(lensIdx + shift2);
    allDistances = [allDistances; lensIdx];
end


if nImgsPerScan > 2
    cntrAScan3 = smooth(double(imgsIn(:,cntrLineIdx,3)),15);
    cntrAScan5(1:50)= 0;
    
    mpd3 = round(0.01*yRes);
    retinaIdx = findALinePeaks(cntrAScan3, mpd3);
    
    
    if stitchOutput.isFlipImg(3)
        shift3 = -size(imgsIn,1);
        retinaIdx = flipud(abs(retinaIdx + shift3));
    end
    
     allDistances = [allDistances; retinaIdx];
end

function initializeFigurePlots(handles, imgInfo)
tol = imgInfo.tol;
cntrIdx = imgInfo.cntrIdx; surfLoc = imgInfo.surfLoc;
nImgs = imgInfo.nImgsPerScan;
idxRange = getIdxRange(imgInfo.topImg, cntrIdx, tol);


%..........................................................................
%Display the top and bottom images
%..........................................................................
topImgDisplay = mat2gray(imgInfo.topImg);
topImgDisplay = imadjust(topImgDisplay, [0.2, 1], [0, 1], 1);

bottomImgDisplay = mat2gray(imgInfo.bottomImg);
bottomImgDisplay = imadjust(bottomImgDisplay, [0.45, 1],[0, 1], 1);

hAxTop = handles.axes_top; hAxBottom = handles.axes_bottom;
hAxInt = [handles.axes_int, handles.axes_int2]; 
hAxDer = [handles.axes_der, handles.axes_der2];

hImgTop = imagesc(topImgDisplay, 'Parent', hAxTop);
hImgBottom = imagesc(bottomImgDisplay, 'Parent', hAxBottom);

hAllTopAxes = [hAxTop, hAxInt(1), hAxDer(1)];
hAllBottomAxes = [hAxBottom, hAxInt(2), hAxDer(2)];

linkaxes(hAllTopAxes, 'y'); linkaxes(hAllBottomAxes,'y');

setappdata(hAxTop, 'img', topImgDisplay); 
setappdata(hAxTop, 'hImg', hImgTop);
setappdata(hAxBottom, 'img', bottomImgDisplay);
setappdata(hAxBottom, 'hImg', hImgBottom);



set([hAxTop, hAxBottom],'XTick',cntrIdx, 'YTick',[]);


%..........................................................................
%Display lines for intensity and derivative
%..........................................................................
intTop = smooth(mean(imgInfo.topImg(:,idxRange),2),15); 
derTop = smooth(diff(intTop),15); derTop = [derTop(1);derTop];
intBottom = mean(imgInfo.bottomImg(:,idxRange),2);
derBottom = smooth(diff(intBottom),15); derBottom = [derBottom(1);derBottom];

yDataTop = (1:size(topImgDisplay,1))'; 
yDataBottom = (1:size(bottomImgDisplay,1))';
hP(1) = plot(hAxInt(1), intTop, yDataTop,'g');
hP(2) = plot(hAxInt(2), intBottom, yDataBottom,'g'); 
setappdata(hAxInt(1), 'hPlot', hP(1)); setappdata(hAxInt(2), 'hPlot', hP(2));

hP(3) = plot(hAxDer(1), derTop, yDataTop,'g'); 
hP(4) = plot(hAxDer(2), derBottom, yDataBottom,'g'); 
setappdata(hAxDer(1), 'hPlot', hP(3)); setappdata(hAxDer(2), 'hPlot', hP(4));

axis([hAxInt, hAxDer], 'tight'); axis([hAxInt, hAxDer], 'ij'); 
set([hAxInt, hAxDer], 'XTick',[], 'YTick',[],'Color','k'); %
set([hAxDer, hAxInt], 'XLimMode','auto');

%..........................................................................
%CREATE THE MOBILE LINE
%..........................................................................
%SIOBHAN MOD: 6/27/2016
thickness = tol; isHorizontal = false;
hLineObj(1) = mobileLine(hAxTop, thickness, isHorizontal);

if nImgs > 2
    hLineObj(2) = mobileLine(hAxBottom, thickness, isHorizontal);
end

hLine = [hLineObj.handle];

%Set the line object as appdata to the line itself
arrayfun(@(iLine) setappdata(iLine, 'lineObj', hLineObj), hLine);
arrayfun(@(iLine) setappdata(iLine, 'hLine', hLine), hLine);


%Whenever one line moves, update the other line
hLineLink = linkprop([hLineObj.handle],'XData');
hLineObj(1).callback = @(iLine) axis(getappdata(hLine), 'tight');

if nImgs > 2
    hLineObj(2).callback = @(iLine) axis(getappdata(hLine), 'tight');
end

setPosition(hLineObj(1), cntrIdx);

arrayfun(@(lineH) iptaddcallback(lineH, 'ButtonDownFcn',...
    @(iLine, edata) lineObjCallback(iLine, handles)), hLine);


%Disable interactivity of line
arrayfun(@(iObj) iObj.disable, hLineObj);
   
%Change the color and style of the line
lineColor = [1,0.5 0.5];
set(hLine, 'EdgeColor', lineColor,'FaceColor', lineColor, 'LineStyle','-.',...
    'Visible','on','HitTest','on');

%Change pointer
%Set Cursor action
iptPointerManager(ancestor(hAxTop, 'figure')); 
iptSetPointerBehavior(hLine, @(f, cp) set(f, 'Pointer', 'arrow'));

%..........................................................................
%CREATE SURFACE POINTS
%..........................................................................
pointLabels = {'AC','PC','AL','PL','AR', 'PR'};
rectConstraintFcn1 = makeConstrainToRectFcn('impoint',[cntrIdx, cntrIdx],...
    [1, size(topImgDisplay,1)]);
rectConstraintFcn2 = makeConstrainToRectFcn('impoint',[cntrIdx, cntrIdx],...
    [1, size(bottomImgDisplay,1)]);

intData_top = [intTop,yDataTop]; intData_bottom = [intBottom,yDataBottom];
derData_top = [derTop,yDataTop]; derData_bottom =  [derBottom, yDataBottom];

hText = zeros(nImgs*2,3); hPointCell = cell(nImgs*2,3);
hPointHandles = cell(nImgs*2,1);



for i = 1:nImgs*2

    if i < 5 
        iAxImg = hAxTop; iAxInt = hAxInt(1); iAxDer = hAxDer(1);
        iRectCons = rectConstraintFcn1;
        iConsFcn = @(pos) intData_top(knnsearch(intData_top, pos),:);
        dConsFcn = @(pos) derData_top(knnsearch(derData_top, pos),:);
        xInt = intData_top(cntrIdx,1); xDer = derData_top(cntrIdx,1);

    else
        iAxImg = hAxBottom; iAxInt = hAxInt(2); iAxDer = hAxDer(2);
        iRectCons = rectConstraintFcn2;
        iConsFcn = @(pos) intData_bottom(knnsearch(intData_bottom, pos),:);
        dConsFcn = @(pos) derData_bottom(knnsearch(derData_bottom, pos),:);
        xInt = intData_bottom(cntrIdx,1); xDer = derData_bottom(cntrIdx,1);
    end
    
    iAxes = [iAxImg, iAxInt, iAxDer];
    cFcn = {iRectCons, iConsFcn, dConsFcn};
    xPoints = [cntrIdx, xInt, xDer];
    offset = [5, 10, 3];
    for j = 1:3
        jAx = iAxes(j); hold(jAx, 'on');
        jPoint = impoint(jAx, [xPoints(j), surfLoc(i)]); 
        setPositionConstraintFcn(jPoint, cFcn{j});     
        setConstrainedPosition(jPoint, getROIPosition(jPoint));
        jPos = getROIPosition(jPoint);
        
        hText(i,j) = text(jPos(1)+offset(j), jPos(2), pointLabels{i},...
            'Parent', jAx,'Color','w');
        hold(jAx, 'off');
        addNewPositionCallback(jPoint,...
            @(pos) pointMoveCallback(jPoint, pos));
        
        hPointHandles{i} = [hPointHandles{i}; get(jPoint, 'Children')];
        
       hPointCell{i,j} = jPoint;
    end
    
end

hPoint = reshape([hPointCell{:}], [nImgs*2, 3]);



title(hAxTop, imgInfo.ID,'Interpreter','none'); 
title(hAxBottom, imgInfo.ID,'Interpreter','none'); 
title(hAxInt(1), 'Intensity');  title(hAxInt(2), 'Intensity'); 
title(hAxDer(1), 'Derivative'); title(hAxDer(2), 'Derivative');


%..........................................................................
%Delete bottom axes if necessary
%..........................................................................

if nImgs < 3
    arrayfun(@cla,hAllBottomAxes); set(hAllBottomAxes, 'XTick',[], 'YTick',[]);
    set(handles.tool_switchView, 'Enable','off','TooltipString','View Retina');
    set(handles.tool_ar, 'Enable','off');
end

setappdata(handles.tool_switchView, 'corneaFlag',true)
axis([handles.axes_top, handles.axes_bottom],'tight');

%..........................................................................
%Save important appdata
%..........................................................................
hFig = ancestor(hAxTop, 'figure');
setappdata(hFig, 'hPoints', hPoint);
setappdata(hFig, 'hPointHandles', hPointHandles);
setappdata(hFig, 'hText', hText);
setappdata(hFig, 'lineObj', hLineObj);
setappdata(hFig, 'cntrIdx', cntrIdx);

function peaksOut = findALinePeaks(aLine, mpd)

[val, idx] = findpeaks(aLine,'MinPeakDistance', mpd,...
    'MinPeakHeight', mean(aLine) + 1.5*std(aLine));

if numel(idx) < 1
    
    [~, peaksOut] = max(aLine); peaksOut = peaksOut-15;
    peaksOut(2) = peaksOut + mpd;

elseif numel(idx) < 2
    peaksOut = val;
    peaksOut(2) = peaksOut + mpd;
else
    %nothing too far apart
    
    sortedPks = sortrows(round([idx, val]),1); %sort by peak value
    %take first and last
    peaksOut = sortedPks([1, end],:);   
    peaksOut = peaksOut(:,1);
    
end

function lineObjCallback(iLine, handles)


hFig = ancestor(iLine, 'figure');

%Only if the user double clicks, continue
if ~strcmpi(get(hFig,'SelectionType'),'open')
    return
end

lineObj = getappdata(iLine, 'lineObj'); %both line objects
hLine = getappdata(iLine, 'hLine'); %both line handles

hPoints = getappdata(hFig, 'hPoints');
hText = getappdata(hFig, 'hText');
isEnabled = lineObj(1).isEnabled;
hPoints = getappdata(hFig, 'hPointHandles');

if iscell(hPoints{1})
    hPointHandles = cell2mat(hPoints{1});
else
    hPointHandles = hPoints{1};
end

lineIdx = ~getappdata(handles.tool_switchView, 'corneaFlag'); lineIdx = lineIdx+1;

if isEnabled %disable it
    %update center index
    imgInfo= getappdata(gcbf, 'imgInfo');
    tol = imgInfo.tol;
    
    linePos = round(lineObj(lineIdx).getPosition);
    cntrIdx = linePos + tol; 
    setPosition(lineObj(lineIdx), linePos);
    axis([handles.axes_top, handles.axes_bottom],'tight');
    setappdata(hFig, 'cntrIdx', cntrIdx);
    
    %update the X Position of all points;
    updatePlots(hFig, cntrIdx);

    %Disable interactivity of line 
    arrayfun(@(iObj) iObj.disable, lineObj);    
    
    %Change the color and style of the line
    lineColor = [1,0.5 0.5];
    set(hLine, 'EdgeColor', lineColor,'FaceColor', lineColor, 'LineStyle','-.',...
        'Visible','on','HitTest','on');
    
    set(hPointHandles, 'Visible','on');
    set(hText, 'Visible', 'on');

    %Set Cursor action
    iptPointerManager(hFig);
    iptSetPointerBehavior(hLine, @(f, cp) set(f, 'Pointer', 'arrow'));

else %enable it
    %hide all points and text
    set(hPointHandles, 'Visible','off');
    set(hText, 'Visible', 'off');
    
    
    %Change the color and style of the line
    lineColor = 'g';
    set(hLine, 'EdgeColor', lineColor,'FaceColor', lineColor, 'LineStyle','-');
    
    arrayfun(@(iObj) iObj.enable, lineObj);
      
    %Set Cursor action
    iptPointerManager(hFig);
    iptSetPointerBehavior(hLine, @(f, cp) set(f, 'Pointer', 'fleur'));

end

function pointMoveCallback(iPoint, pos)
%This function updates the y position of all text and points associated
%with a surface

hFig = ancestor(get(iPoint, 'Children'), 'figure'); hFig = hFig{1};
hPoints = getappdata(hFig, 'hPoints');
hText = getappdata(hFig, 'hText');
newY = round(pos(2));
handles = guidata(hFig);
[surfIdx,~] = find(hPoints == iPoint);
surfText = hText(surfIdx,:);
surfPoints = hPoints(surfIdx,:);


surfText(~ishandle(surfText)) = [];
surfPoints(~ishandle(surfText)) = [];

%..........................................................................
%SNAP TO VERTEX OPTION
%..........................................................................
if ~isempty(handles) && get(handles.tool_snap2vertex,'Value')
    [~, col] = find(hPoints == iPoint );
    
    if col == 3
        %snap to peak of derivative
        hBothAx = [handles.axes_der, handles.axes_der2];
        hPanel = [handles.panel_top, handles.panel_bottom];
        hActiveAx = hBothAx(strcmpi(get(hPanel, 'Visible'),'on'));
        
        hP = getappdata(hActiveAx, 'hPlot');
        x = get(hP, 'XData'); xMask = false(size(x));
        xThresh = mean(abs(x)) +std(abs(x));
        xMask(newY-30:newY+30) = true; xMask = xMask*max(x);
        x(~xMask) = 0;
        
        if any(abs(x)>xThresh)
            %snap to the vertex
            [~,newY] = max(abs(x));
        end      
        
    end
else isempty(handles)
    [~, col] = find(hPoints == iPoint );
    
    if col == 2
        %snap to peak of derivative
        hAxInt = getappdata(hFig, 'hAxInt');
        
        hP = getappdata(hAxInt, 'hPlot');
        x = get(hP, 'XData'); xMask = false(size(x));
        xThresh = mean(abs(x)) +std(abs(x));
        xMask(newY-30:newY+30) = true; xMask = xMask*max(x);
        x(~xMask) = 0;
        
        if any(abs(x)>xThresh)
            %snap to the vertex
            [~,newY] = max(abs(x));
        end      
        
    end
    
end
%..........................................................................


%Update all Points
textOffset = [5 10 1];
for j = 1:numel(surfPoints)
    jPoint = surfPoints(j);
    
    
    jText = surfText(j); jKids = get(jPoint, 'Children');
    jAx = ancestor(jKids(1), 'axes'); hPlot = getappdata(jAx, 'hPlot');
    oldPos = getROIPosition(jPoint); newPos = oldPos;
    posConFcn = getPositionConstraintFcn(jPoint);
    
    
    if ~isempty(hPlot) && ishandle(hPlot)
        x = get(hPlot,'XData');
        newPos = [x(newY), newY];
    else
        newPos(2) = newY; 
    end 
    
    %Set position for point
    newPos = posConFcn(newPos);
    set(jKids, 'XData', newPos(1), 'YData', newPos(2));
    jTextPos = newPos; jTextPos(1) = jTextPos(1) + textOffset(j); %pause(0.01);
    %Set the text
    set(jText,'Position', jTextPos);
    
end

%Shift only the second coordinate of all text

textPos = cell2mat(get(surfText,'Position'));
textPos(:,2) = pos(2);
set(surfText, {'Position'}, num2cell(textPos, 2));

%Make that point an active color
arrayfun(@(iPoint) setColor(iPoint, 'b'), hPoints);
arrayfun(@(iPoint) setColor(iPoint, 'm'), hPoints(surfIdx,:));



function pointMoveCallback_mod(iPoint, pos)
%This function updates the y position of all text and points associated
%with a surface

hFig = ancestor(get(iPoint, 'Children'), 'figure'); hFig = hFig{1};
hPoints = getappdata(hFig, 'hPoints');
hText = getappdata(hFig, 'hText');
newY = round(pos(2));
handles = guidata(hFig);
[surfIdx,~] = find(hPoints == iPoint);
surfText = hText(surfIdx,:);
surfPoints = hPoints(surfIdx,:);


surfText(~ishandle(surfText)) = [];
surfPoints(~ishandle(surfText)) = [];

%..........................................................................
%SNAP TO VERTEX OPTION
%..........................................................................
if false %get(handles.tool_snap2vertex,'Value')
    [~, col] = find(hPoints == iPoint );
    
    if col == 3
        %snap to peak of derivative
        hBothAx = [handles.axes_der, handles.axes_der2];
        hPanel = [handles.panel_top, handles.panel_bottom];
        hActiveAx = hBothAx(strcmpi(get(hPanel, 'Visible'),'on'));
        
        hP = getappdata(hActiveAx, 'hPlot');
        x = get(hP, 'XData'); xMask = false(size(x));
        xThresh = mean(abs(x)) +std(abs(x));
        xMask(newY-30:newY+30) = true; xMask = xMask*max(x);
        x(~xMask) = 0;
        
        if any(abs(x)>xThresh)
            %snap to the vertex
            [~,newY] = max(abs(x));
        end      
        
    end
end
%..........................................................................


%Update all Points
textOffset = [5 10 1];
for j = 1:numel(surfPoints)
    jPoint = surfPoints(j);
    
    
    jText = surfText(j); jKids = get(jPoint, 'Children');
    jAx = ancestor(jKids(1), 'axes'); hPlot = getappdata(jAx, 'hPlot');
    oldPos = getROIPosition(jPoint); newPos = oldPos;
    posConFcn = getPositionConstraintFcn(jPoint);
    
    
    if ~isempty(hPlot) && ishandle(hPlot)
        x = get(hPlot,'XData');
        newPos = [x(newY), newY];
    else
        newPos(2) = newY; 
    end 
    
    %Set position for point
    newPos = posConFcn(newPos);
    set(jKids, 'XData', newPos(1), 'YData', newPos(2));
    jTextPos = newPos; jTextPos(1) = jTextPos(1) + textOffset(j); %pause(0.01);
    %Set the text
    set(jText,'Position', jTextPos);
    
end

%Shift only the second coordinate of all text

textPos = cell2mat(get(surfText,'Position'));
textPos(:,2) = pos(2);
set(surfText, {'Position'}, num2cell(textPos, 2));

%Make that point an active color
arrayfun(@(iPoint) setColor(iPoint, 'b'), hPoints);
arrayfun(@(iPoint) setColor(iPoint, 'm'), hPoints(surfIdx,:));


function updatePlots(hFig, cntrIdx)

imgInfo = getappdata(hFig, 'imgInfo');
hPoints = getappdata(hFig, 'hPoints');
hText = getappdata(hFig, 'hText');
handles= guidata(hFig);

tol = imgInfo.tol;
idxRange = getIdxRange(imgInfo.topImg, cntrIdx, tol);


%..........................................................................
%Change the A-lines
%..........................................................................

intTop = smooth(mean(imgInfo.topImg(:,idxRange),2),15); 
derTop = smooth(diff(intTop),15); derTop = [derTop(1);derTop];

hP_int =getappdata(handles.axes_int, 'hPlot');
hP_der =getappdata(handles.axes_der, 'hPlot');

set(hP_int, 'XData', intTop); set(hP_der, 'XData', derTop);
axis([handles.axes_int, handles.axes_der], 'tight');
set([handles.axes_int, handles.axes_der], 'XLimMode','auto')
nImgs = imgInfo.nImgsPerScan;


if nImgs > 2 %Plot the bottom
    intBottom = mean(imgInfo.bottomImg(:,idxRange),2);
    derBottom = smooth(diff(intBottom),15); derBottom = [derBottom(1);derBottom];

    hP_int2 =getappdata(handles.axes_int2, 'hPlot');
    hP_der2 =getappdata(handles.axes_der2, 'hPlot');
    set(hP_int2, 'XData', intBottom); 
    set(hP_der2, 'XData', derBottom);
    axis([handles.axes_int2, handles.axes_der2], 'tight');
    set([handles.axes_int2, handles.axes_der2], 'XLimMode','auto')
end

set([handles.axes_top, handles.axes_bottom],'XTick', cntrIdx);


%..........................................................................
%Update the points
%..........................................................................

for j = 1:numel(hPoints)
    jPoint = hPoints(j); jText = hText(j); 
    jAx = ancestor(jText, 'axes'); hPlot = getappdata(jAx, 'hPlot');
    newPos = getROIPosition(jPoint);
    
    xlims = get(jAx, 'XLim'); offset = ceil(abs(diff(xlims))/80);
    if ~isempty(hPlot)
        x = get(hPlot,'XData'); y = get(hPlot, 'YData'); plotData = [x;y]';
        oldY = round(newPos(2));
        newPos = [x(oldY), y(oldY)]; 
        
        %update position constraint fcn
        cFcn = @(pos) plotData(knnsearch(plotData, pos),:);
        setPositionConstraintFcn(jPoint, cFcn);
        
        
    else
        jKids = get(jPoint, 'Children'); jParentAx = ancestor(jKids(1),'axes');
        yLims = get(jParentAx, 'YLim'); yLims = [ceil(yLims(1)), floor(yLims(2))];
        
        %update position constraint fcn
        cFcn = makeConstrainToRectFcn('impoint',[cntrIdx, cntrIdx], yLims);
        setPositionConstraintFcn(jPoint, cFcn);
        
        newPos(1) = cntrIdx; 
    end
    
    %Set position for point
    setConstrainedPosition(jPoint, newPos);
    jTextPos = newPos; jTextPos(1) = jTextPos(1) + offset;
    %Set the text
    set(jText,'Position', jTextPos);
end

function idxRange = getIdxRange(img, cntrIdx, tol)

idxRange = max([1, cntrIdx-tol]):min([size(img, 2), cntrIdx+tol]);
%==========================================================================
%Toolbar callbacks
%==========================================================================

function tool_zoomIn_Callback(hObject, eventdata, handles)
% hObject    handle to tool_zoomIn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

state = {'off','on'};

hFig = ancestor(hObject, 'figure'); hZoom = zoom(hFig);
set(hZoom, 'Direction', 'in', 'Enable', state{get(hObject, 'Value')+1});

set([handles.tool_zoomOut, handles.tool_datacursor, handles.tool_pan], 'Value', 0);

function tool_zoomOut_Callback(hObject, eventdata, handles)
% hObject    handle to tool_zoomOut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

state = {'off','on'};

% Hint: get(hObject,'Value') returns toggle state of tool_zoomOut
hFig = ancestor(hObject, 'figure'); hZoom = zoom(hFig); 
set(hZoom, 'Direction', 'out', 'Enable', state{get(hObject, 'Value')+1});

set([handles.tool_zoomIn, handles.tool_datacursor, handles.tool_pan], 'Value', 0);

function tool_pan_Callback(hObject, eventdata, handles)
% hObject    handle to tool_pan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

state = {'off','on'};

pan(state{get(hObject, 'Value')+1});
set([handles.tool_zoomOut, handles.tool_datacursor, handles.tool_zoomIn], 'Value', 0);

function tool_datacursor_Callback(hObject, eventdata, handles)
% hObject    handle to tool_datacursor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

state = {'off','on'};

datacursormode(state{get(hObject, 'Value')+1});

set([handles.tool_zoomIn, handles.tool_zoomOut, handles.tool_pan], 'Value', 0);

function tool_switchView_Callback(hObject, eventdata, handles)
% hObject    handle to tool_switchView (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

corneaFlag = getappdata(handles.tool_switchView, 'corneaFlag'); %is cornea visible
hTop = handles.panel_top;
hBottom = handles.panel_bottom;

if isempty(corneaFlag) || corneaFlag %switch to lense
    topMode = 'off'; bottomMode = 'on'; toolStr = 'View Cornea';
else  %switch to cornea
    topMode = 'on'; bottomMode = 'off'; toolStr = 'View Lens';
end

set(hTop, 'Visible', topMode); set(hBottom,'Visible', bottomMode);
setappdata(handles.tool_switchView, 'corneaFlag',~corneaFlag);
set(handles.tool_switchView, 'TooltipString',toolStr);

function tool_ar_Callback(hObject, eventdata, handles)
% hObject    handle to tool_ar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hFig = ancestor(hObject, 'figure');
hPointHandles = getappdata(hFig, 'hPointHandles');
hText = getappdata(hFig, 'hText');

hPoints_ar = hPointHandles{5}; hText_ar = hText(5,:);
state = {'off','on'};
set([hPoints_ar; hText_ar'], 'Visible',state{get(hObject, 'Value')+1});

function tool_export_Callback(hObject, eventdata, handles)
% hObject    handle to tool_export (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hFig = ancestor(hObject, 'figure');
hDummyFig = getappdata(hFig, 'hDummyFig');

biometryOut = calculateBiometry(handles);

if isempty(hDummyFig) || ~ishandle(hDummyFig)
    uInput = questdlg('Where would you like to export the data?', 'Export Data',...
       'Excel', 'Print'); %,Worksapce

uInput = 'excel'; %for now
    switch lower(uInput)

        case 'workspace'
            assignin('base','biometryOut', biometryOut);
        case 'excel'
            
            biometryXL =  num2cell([biometryOut.centerALine,...
                biometryOut.optDist', biometryOut.optEyeLength,...
                biometryOut.geoDist', biometryOut.geoEyeLength]);
            
            mat2xl(biometryXL);
        case 'print'
    end    
else
    
    
end

function tool_viewBiometry_Callback(hObject, eventdata, handles)
% hObject    handle to tool_viewBiometry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


hFig = ancestor(hObject, 'figure');
biometryOut = calculateBiometry(handles);
plotBiometryResults(getappdata(hFig, 'imgInfo'), biometryOut);

function setToggleImgs(hToggles, imgData)

for i = 1:numel(hToggles)
    iToggle = hToggles(i); iData = imgData{i};
    set(iToggle, 'CData', iData{get(iToggle, 'Value')+1});
end
%==========================================================================
%Calculate biometry
%==========================================================================

function biometryOut = calculateBiometry(handles)

hFig = ancestor(handles.axes_top, 'figure');

imgInfo = getappdata(hFig, 'imgInfo'); 
hPoints = getappdata(hFig, 'hPoints');
cntrIdx = getappdata(hFig, 'cntrIdx');

%Get all the points 
hMainPoints = hPoints(:,1);

pointLocationOrig = cell2mat(arrayfun(@(iPoint) getROIPosition(iPoint),...
    hMainPoints,'UniformOutput', false));
pointLocation = pointLocationOrig;


stitchOutput = imgInfo.stitchOutput;
nImgs = imgInfo.nImgsPerScan;
yRes =  imgInfo.yRes;
n = imgInfo.n;


%shift the retina point accordingly 
if nImgs > 2
    absImgStart = stitchOutput.absoluteImgStart;
    frameShift = absImgStart - 1;
    
    %don't worry about inversion, image is already right side up
    
    pointLocation(5:6,2) = pointLocation(5:6,2) + frameShift(3);
end

coord_mm = pointLocation./yRes;
yLoc = pointLocation(:,2); yLoc_mm = yLoc./yRes;

if nImgs > 2
    if get(handles.tool_ar, 'Value')
        lineIdx = [1 2 3 4 6 5];
        
    else
        lineIdx = [1 2 3 4 6];
       
    end
    yLoc = yLoc(lineIdx);
    yLoc_mm = yLoc_mm(lineIdx);
    pointLocation = pointLocation(lineIdx,:);
    coord_mm = coord_mm(lineIdx,:);
end

%calculate distances with and without anterior retina

opticalD_mm = abs(diff(yLoc_mm));
geoD_mm = opticalD_mm./n(1:numel(opticalD_mm))';

biometryOut.optDist = opticalD_mm;
biometryOut.geoDist = geoD_mm;
biometryOut.coord_pxl = pointLocation;
biometryOut.coord_mm = coord_mm;
biometryOut.geoEyeLength = sum(geoD_mm(1:3));
biometryOut.optEyeLength = sum(opticalD_mm(1:3));
biometryOut.pupilWidth = [];
biometryOut.centerALine = cntrIdx;
biometryOut.img = stitchOutput.stitchedImg;
biometryOut.rawCoord_pxl = pointLocationOrig;

function plotBiometryResults(imgInfo, bStruct)

optDistances_mm = bStruct.optDist;
geoDistances_mm = bStruct.geoDist;
retinaFlag = false;
handles = guidata(gcbf);

if imgInfo.nImgsPerScan > 2
    
    
    if get(handles.tool_ar,'Value')
    uInput = questdlg('What would you like to display?','Display Biometry',...
        'Cornea+Lens', 'Retina','All','Cornea+Lens');
    end
    
    if isempty(uInput)
        return
    end
    
    switch lower(uInput)
        case 'cornea+lens'
            optDistances_mm = optDistances_mm(1:3);
            geoDistances_mm = geoDistances_mm(1:3);
            img = imgInfo.topImg;
        case 'retina'
            optDistances_mm = optDistances_mm(end);
            geoDistances_mm = geoDistances_mm(end);
            img = imgInfo.bottomImg;
            retinaFlag = true;
        otherwise
            img = bStruct.img;
    end
else
    img = imgInfo.topImg;
end

figure; hAx = gca;
imagesc(img), colormap gray; set(gca,'XTick',[], 'YTick',[]);
hold on;

nDistances = numel(geoDistances_mm);

if retinaFlag
    coord_pxl = bStruct.rawCoord_pxl(5:end,:);
else
   coord_pxl = bStruct.coord_pxl;
end

xCenterLine = coord_pxl(:,1); xCenterLine = xCenterLine(1:nDistances+1);
yCenterLine = coord_pxl(:,2); yCenterLine = yCenterLine(1:nDistances+1);

allXCenterLines = [xCenterLine, circshift(xCenterLine, [-1,0])]';
allXCenterLines(:,end) = [];


if nDistances == 5
    %shift the intra-retinal distances for plot reasons
    allXCenterLines(:,end) = allXCenterLines(:,end) + 5;
end

allYCenterLines = [yCenterLine, circshift(yCenterLine, [-1,0])]';
allYCenterLines(:,end) = [];

optDistances2Plot = optDistances_mm(1:nDistances);
geoDistances2Plot = geoDistances_mm(1:nDistances);


arrowLocationX = max(allXCenterLines,[],1)+ 10;
arrowLocationY = round(mean(allYCenterLines,1));

hold(hAx, 'on'); hP = plot(hAx, allXCenterLines, allYCenterLines);
labelColor = 'ygrcm'; 
arrayfun(@(i) set(hP(i), 'Color', labelColor(i)), 1:nDistances);

arrayfun(@(i) text(arrowLocationX(i),arrowLocationY(i),...
    ['\leftarrow', sprintf(' d = %2.3fmm (%2.3fmm)',...
    optDistances2Plot(i), geoDistances2Plot(i))],...
    'HorizontalAlignment', 'Left',...
    'Color', labelColor(i), 'Parent',hAx), 1:nDistances);

title(hAx, imgInfo.ID);

function pos = getROIPosition(hPoint)

hKids = get(hPoint, 'Children'); 
pos = [get(hKids(1), 'XData'), get(hKids(1), 'YData')];


% --- Executes on button press in tool_previewFull.
function tool_previewFull_Callback(hObject, eventdata, handles)
% hObject    handle to tool_previewFull (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Initialize vars
cntrIdx = getappdata(handles.figure1, 'cntrIdx'); %210
imgInfo = getappdata(handles.figure1, 'imgInfo');

%%
%get full stitched image
stitchOutput= imgInfo.stitchOutput;
stitchImg = stitchOutput.stitchedImg;
topImg = imgInfo.topImg; bottomImg = imgInfo.bottomImg;


topImgDisplay = mat2gray(topImg);
topImgDisplay = imadjust(topImgDisplay, [0.2, 1], [0, 1], 1);
bottomImgDisplay = mat2gray(bottomImg);
bottomImgDisplay = imadjust(bottomImgDisplay, [0.45, 1],[0, 1], 1);

img2End = stitchOutput.absoluteImgEnd(2);
img3Range = stitchOutput.absoluteStitchStart(3,:);

stitchImg(1:img2End,:) = topImgDisplay; %topImg; 
stitchImg(img3Range(1):img3Range(2),:,1) = bottomImgDisplay; %bottomImg;

%get intensity profile
iProfile = stitchImg(:,cntrIdx);

%%
%Set figure aesthetics
hFig21 = figure(21); clf
set(hFig21, 'Color','w')
set(hFig21, 'Units', 'characters');
set(hFig21, 'Position',[103.8 2.846, 220.8, 59.923]);
hAxImg = subplot(1,2,1); hAxInt = subplot(1,2,2);
set(hAxImg, 'Units', 'characters');
set(hAxImg, 'Position', [9.2, 5.53, 136.2, 51]);



%%
%show full stitched image
imagesc(stitchImg,'Parent', hAxImg); colormap gray;

%plot central A-Line
axes(hAxImg); hold on; 
hLine = plot([cntrIdx, cntrIdx], get(hAxImg, 'YLim'),'--r','Parent', hAxImg);
hold off;

%show the intensity Profile
hPI = plot(iProfile,1:numel(iProfile),...
    'Color', 'g','Parent',hAxInt); axis tight
set(hAxInt, 'Units', 'characters');
set(hAxInt, 'Position', [159, 5.53, 36.8, 51]);
set(hAxInt, 'YDir', 'reverse')
set(hAxInt, 'Color','k')


setappdata(hAxInt, 'hPlot', hPI);
set([hAxImg, hAxInt],'XTick',[],'YTick',[]);
%%

%..........................................................................
%CREATE SURFACE POINTS
%..........................................................................
pointLabels = {'AC','PC','AL','PL','AR', 'PR'};
rectConstraintFcn1 = makeConstrainToRectFcn('impoint',[cntrIdx, cntrIdx],...
    [1, size(stitchImg,1)]);
rectConstraintFcn2 = makeConstrainToRectFcn('impoint',[cntrIdx, cntrIdx],...
    [1, size(stitchImg,1)]);



intTop = smooth(imgInfo.topImg(:,cntrIdx),15); 
derTop = smooth(diff(intTop),15); derTop = [derTop(1);derTop];
intBottom = imgInfo.bottomImg(:,cntrIdx);
derBottom = smooth(diff(intBottom),15); derBottom = [derBottom(1);derBottom];

yDataTop = (1:size(stitchImg,1))'; 
yDataBottom = (1:size(stitchImg,1))';
intDataFull = [iProfile,yDataTop]; 
nAxes = 2; %iamge and int

nImgs = imgInfo.nImgsPerScan;
hText = zeros(nImgs*2,nAxes); hPointCell = cell(nImgs*2,nAxes);
hPointHandles = cell(nImgs*2,1);



hPoints_old = getappdata(handles.figure1, 'hPoints');
if isempty(hPoints_old)
    surfLoc = imgInfo.surfLoc; surfLoc(end-1) = 7000; surfLoc(end) = 7200;
else
    surfLocCell = arrayfun(@getROIPosition, hPoints_old(:,1),...
        'UniformOutput',false);
    surfLocMat = cat(1, surfLocCell{:});
    surfLoc = surfLocMat(:,2);
    surfLoc(end) = surfLoc(end)+stitchOutput.absoluteImgStart(end);
end


%%
%Plot first set of markers
for i = 1:nImgs*2
    iAxImg = hAxImg; iAxInt = hAxInt(1);
    iRectCons = rectConstraintFcn1;
    iConsFcn = @(pos) intDataFull(knnsearch(intDataFull, pos),:);
    xInt = intDataFull(cntrIdx,1);

    
    iAxes = [iAxImg, iAxInt];
    cFcn = {iRectCons, iConsFcn};
    xPoints = [cntrIdx, xInt];
    offset = [5, mean(iProfile)*1, 3];
    
    for j = 1:numel(iAxes)
        jAx = iAxes(j); hold(jAx, 'on');
        jPoint = impoint(jAx, [xPoints(j), surfLoc(i)]); 
        setPositionConstraintFcn(jPoint, cFcn{j});     
        setConstrainedPosition(jPoint, getROIPosition(jPoint));
        jPos = getROIPosition(jPoint);
        
        hText(i,j) = text(jPos(1)+offset(j), jPos(2), pointLabels{i},...nII
            'Parent', jAx,'Color','w');
        hold(jAx, 'off');
        addNewPositionCallback(jPoint,...
            @(pos) pointMoveCallback(jPoint, pos));
        
        hPointHandles{i} = [hPointHandles{i}; get(jPoint, 'Children')];
        
       hPointCell{i,j} = jPoint;
    end
    
end

%%
hPoint = reshape([hPointCell{:}], [nImgs*2, nAxes]);

hFig21 = ancestor(hAxImg, 'figure');
setappdata(hFig21, 'hPoints', hPoint);
setappdata(hFig21, 'hPointHandles', hPointHandles);
setappdata(hFig21, 'hText', hText);
% setappdata(hFig, 'lineObj', hLineObj);
setappdata(hFig21, 'cntrIdx', cntrIdx);
setappdata(hFig21, 'handles', handles);
setappdata(hFig21, 'hAxInt', hAxInt);

%%
%Hide Ant Retina
if nImgs > 2
    set(hPoint(5,1),'Visible','off');
    set(hPoint(5,2),'Visible','off');
    
    set(hText(5,:),'Visible','off');
end