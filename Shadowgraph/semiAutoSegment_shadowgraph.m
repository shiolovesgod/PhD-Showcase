function varargout = semiAutoSegment_shadowgraph(varargin)
% SEMIAUTOSEGMENT MATLAB code for semiAutoSegment.fig
%      SEMIAUTOSEGMENT, by itself, creates a new SEMIAUTOSEGMENT or raises the existing
%      singleton*.
%
%      H = SEMIAUTOSEGMENT returns the handle to a new SEMIAUTOSEGMENT or the handle to
%      the existing singleton*.
%
%      SEMIAUTOSEGMENT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SEMIAUTOSEGMENT.M with the given input arguments.
%
%      SEMIAUTOSEGMENT('Property','Value',...) creates a new SEMIAUTOSEGMENT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before semiAutoSegment_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to semiAutoSegment_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help semiAutoSegment

% Last Modified by GUIDE v2.5 03-Feb-2017 11:33:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @semiAutoSegment_shadowgraph_OpeningFcn, ...
    'gui_OutputFcn',  @semiAutoSegment_shadowgraph_OutputFcn, ...
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


% --- Executes just before semiAutoSegment is made visible.
function semiAutoSegment_shadowgraph_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to semiAutoSegment (see VARARGIN)

% Choose default command line output for semiAutoSegment
handles.output = hObject;
img_splash = get(hObject, 'UserData');
isUseSplash = false;

if ~isempty(img_splash) && isUseSplash && numel(varargin) > 1
    hSplash = splash(img_splash);
else
    hSplash = -1;
end

setappdata(hObject, 'hSplash',hSplash); %save the handle somewhere

%parse inputs
if isempty(varargin)
    userData = loadImg(); %working on the retina right now
    
    if isempty(userData)
        delete(hObject);
        return
    end
    
    
else %THIS NEEDS TO BE MODIFIED ONCE I CALL THIS PROGRAM FROM AN M-FILE
    %fields of userData needed: fileLocation, xyRes, scanWidth, imgSize,
    %imgIn {raw OCT image}
    userData = varargin{1}; %raw img
    
    if userData.scanWidth > 10 && userData.imgType < 5 %not shadowgraph imgs
        
        [imgOut, ~] = cropImg(userData.imgIn, userData.xyRes, 8.5);
        userData.imgIn = imgOut;
    end
end

%if it is an exvivo lens, initialize rotations and variants
if any(userData.imgType == [4, 5]) %OCT or shadowgraph lens ex-vivo
    userData = initialize_exVivoLens(userData);
    
    if numel(varargin) > 0
        default.fitWindow = userData.params.fitWindow;
    end    
end

if numel(varargin) > 1
    userData.runSilent = varargin{2};
else
    userData.runSilent = false;
end


%..........................................................................
%Standard Default Settings
%..........................................................................

userData.variantIdx = 3; %Green chanel  (1-RGB, 2R,3G,4B,5R-G,6R-B)

if exist('default','var') && isfield(default, 'fitWindow')
    %do nothing
else
    default.fitWindow = 6; %6mm
end

default.nRecursions = 5;
default.fitTol = 0.5; %start fit tolerance width in mm
default.fitDivFactor = 1.5;

%set default colors
highlight1 = [0.078, 0.169, 0.549]; %dark blue
highlight2 = [0.957, 0.455, 0]; %yellow/gold
highlight3 = [0.275, 0.753, 0.933]; %light blue
inactive1= [0.314, 0.314, 0.314]; %dark grey
    inactive2 = [0.8, 0.8, 0.8]; %light grey

colors.inactivePanel = inactive1;
colors.inactiveText = inactive2;
colors.inactiveTitle = highlight2;
colors.activePanel = highlight1;
colors.activeText = [1 1 1]; %white
colors.activeTitle = highlight3;

default.colors = colors;

%..........................................................................

%POSSIBLE Previous Settings
[userData, hasSettings]= loadParameters(userData);
userData.presetSettings = hasSettings;

userData.default = default;

if ~hasSettings
    set(hObject, 'Name', ancestorDir(fileparts(userData.fileLocation),0))
elseif isfield(userData, 'ID')
    set(hObject, 'Name', userData.ID);
end


switch userData.imgType
    case {1, 2}
        %         hasSettings = false;
        try
            [userData, ~, ~] = autoSegment(userData.imgIn, userData, hasSettings);
        catch
            [userData, ~, ~] = autoSegment(userData.imgIn, userData, false);
        end
    
    case 3
        hSplash.dispose;
        delete(hObject);
        errordlg('Retina Segmentation GUI is still in development.');
        return
%         [userData, ~, ~] = autoSegmentRetina(userData.imgIn, userData, hasSettings);
    case 4
        userData = autoSegment(userData.imgIn, userData, hasSettings);
    case 5
        
        userData = autoSegment(userData.imgIn, userData, hasSettings);
    case 6
        %we're not here yet....coronal

end

if ~isfield(userData,'zoomRegion')
    userData.zoomRegion = []; %where to focus the zoom
end

setappdata(hObject, 'userData', userData);


%No need to prep figure if running silent
if userData.runSilent
    % Update handles structure
    guidata(hObject, handles);
    return
    
end

initalizeSettingValues(handles, userData);


%This makes the code visible
updateThumbnail(handles, userData, {'1','2','3','4','5','6','7','8','9','10'});
%Update the plots
initializeFigure(handles, userData);

%Create a timer object and set it into figure appdata
hoverTime = 1; %seconds
tObj = timer('StartDelay', hoverTime,...
    'TimerFcn',@(x,y) hoverCallback(x,y,handles), 'ExecutionMode','singleShot');
setappdata(hObject, 'timerObj',tObj);
set(hObject, 'WindowButtonMotionFcn',...
    @(hObj, edata) mouseMoveCallback(hObj, edata, handles))


%delete splash screen
if ishandle(hSplash)
    hSplash.dispose;
end

set(hObject, 'CloseRequestFcn', @(hObj, edata) figure1_CloseRequestFcn(hObj, edata, guidata(hObj)));

% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = semiAutoSegment_shadowgraph_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
if isstruct(handles)
    varargout{1} = handles.output;
end

if ~ishandle(hObject)
    return
end

%Kil the splash screen if necessary
hSplash = getappdata(hObject, 'hSplash');
if ~isempty(hSplash) && ishandle(hSplash)
    hSplash.dispose;
end

hDummyFig = figure('Visible','off','IntegerHandle','off','HandleVisibility','off');
setappdata(hObject, 'hInvisibleFig', hDummyFig);
setappdata(hDummyFig, 'exportFlag', false);

if nargout > 0
    userData = getappdata(hObject,'userData');
    
    if isstruct(userData) && isfield(userData, 'runSilent') && userData.runSilent
        %just output the auto segment
        setappdata(hDummyFig, 'exportFlag', true);

    else
        set(hObject,'Visible','on');
        uiwait_mod(hDummyFig);
    end
    
    try
        isApply2All = get(handles.check_apply2all, 'Value');
        isOverwrite = get(handles.check_overwrite, 'Value') && isApply2All;
        
    catch
        isApply2All = false;
        isOverwrite = false;
    end
    
    if ishandle(hDummyFig)
        exportFlag  = getappdata(hDummyFig, 'exportFlag');
        wasHandle = true;
        delete(hDummyFig);        
    else
        wasHandle = false;
        exportFlag = false;
    end
    
    if exportFlag
        userData = getappdata(hObject,'userData');
        if userData.imgType > 4
            %Shadowgraph Output
            %Generate Output as a table
            try
                varargout{1} = generateTableOutput(userData);
                varargout{2} = [isApply2All, isOverwrite];
            catch
                varargout{1} = []; %there was a problem
                varargout{2} = [];
            end
        else
            %OCT Software Output
            varargout{1} = userData;
            varargout{2} = [isApply2All, isOverwrite];
        end
        delete(hObject);
    else
        varargout{1} = [];
        varargout{2} = [];
        
        if wasHandle
            delete(hObject);
        end
    end
else
    set(hObject,'Visible','on');
end


function [userData, hasSettings] = loadParameters(userData)
% userData.presetSettings = false;

%Initialize Output
hasSettings = false;

if ~isfield(userData, 'params')
    return
end


switch userData.imgType
    case {4,5} %Loaded from shadowgraph software
        hasSettings = userData.presetSettings;
        
        %Find the right variant Image
        if isfield(userData, 'variantFcn') && hasSettings %has preset settings
            %figure out what the variant Idx should be
            ImgVariants = userData.ImgVariants;
            fcnHandles = {ImgVariants.fcn_cdata};
            fcnStr = cellfun(@func2str, fcnHandles, 'UniformOutput', false);
            
            
            variantIdx = find(strcmpi(fcnStr, userData.variantFcn),1);
            if ~isempty(variantIdx)
                userData.variantIdx = variantIdx;
                userData = updateImgVariant(userData, variantIdx);
            end
            
        end
        
    otherwise
        %Loaded from OCT software
        hasSettings = true;
        
        userData = rmfield(userData, {'img', 'allArtifactIdx', 'intMask',...
            'artifactMask', 'segPoints', 'allGroups', 'validPixelMask',...
            'allPoints', 'finalFit'});
        
        
        %variable parameters
        imgF =  filterImage(userData.imgIn, userData.params.sigma);
        userData.filteredI = imgF;
        userData.imgF = imgF;
        userData.allArtifactIdx = findArtifactIdx(imgF, userData);
        [userData.img, userData.artifactMask] = removeArtifact(imgF,...
            userData.artifactIdx, userData.params.artifactTol);
        
end

function userData = loadImg(imgType)
%load an image
%imgType: 1=cornea, 2=lens, 3=retina, 4=ex-vivo OCT lens, 5=shadowgraph
    %lens


if nargin < 1
    imgType = []; %default this program to shadowgraph images
end

if isempty(imgType)
    [imgType, isOk] = listdlg('PromptString', 'Please select the image type:',...
        'ListString',{'Cornea','Lens', 'Retina','Ex-Vivo Lens', 'Shadowgraph Lens'},...
        'SelectionMode','single','Name', 'Image Type', 'InitialValue',1,...
        'ListSize',[160 75]);

    if ~isOk
        userData = [];
        return
    end
end


[img, fileLocation] = uigetimg('','off');

if isempty(img)
    userData = [];
    return
end

%..........................................................................
%Get scan width
%..........................................................................
[ySize, xSize] = size(img); info = imfinfo(fileLocation);

if isfield(info, 'ImageDescription')
    description = info.ImageDescription;
    resStr = regexp(description, '\d+\.\d+', 'match');
    res =str2double(resStr);
else
    res = repmat(1/0.0064102564,[1,2]); %pixels/mm
end

scanWidth = xSize/res(1); %[width, height]
%..........................................................................

if imgType == 5 %if it is NOT OCT
    xyRes = res;
else
    scanDepth = 10.43; 
    xyRes = [xSize/scanWidth, ySize/scanDepth];
end


userData.fileLocation = fileLocation;
userData.scanWidth = scanWidth;
userData.imgSize = size(img);
userData.xyRes = xyRes;
userData.imgIn = img;

userData.imgType = imgType; %1-cornea, 2-lens, 3-retina, 4, ex-vivo lens, 5-shadowgraph lens

function initializeFigure(handles, userData)


%Sliders
maxTol = round(size(userData.img, 1)/40)*10;

%tag, min, max, stepsize (integer)
sliderInfo = {'fitWindow', [1 8];...
    'nRecursion', [1 15]};

% set(handles.slider_fitWindow, 'SliderStep', repmat(1/diff([1 8])*4,[1 2]));


%Set max, min and range on slider functions
arrayfun(@(i) set(handles.(['slider_',sliderInfo{i,1}]),...
    'Min', sliderInfo{i,2}(1), 'Max', sliderInfo{i,2}(2),...
    'SliderStep', repmat(1/abs(diff(sliderInfo{i,2})),[1,2])),...
    1:size(sliderInfo,1));

function updateThumbnail(handles, userData, thumbnailNo, hAx)
%thumbNailNo = char

if nargin < 3
    %update all
    thumbnailNo = {'1','2','3','4','5','6','7','8','9','10'};
end

iThumbnail = thumbnailNo;

if ~iscell(iThumbnail)
    iThumbnail = cellstr(iThumbnail);
end

isShowMasks = get(handles.radio_showMasks,'Value');

while ~isempty(iThumbnail)
    
    if ~isempty(strfind(iThumbnail{1},'b'))
        iSurfNo = 2;
    else
        iSurfNo = 1;
    end
    
    if nargin > 3
        iAx = hAx; cla(iAx);
        try rmappdata(iAx, 'hImg'); end
        iCallback = '';
    else
        iAx = handles.(['axes', iThumbnail{1}]); cla(iAx);
        iCallback = @(varargin) thumbnailCallback(handles, iAx, iSurfNo);
    end
    
    
    
    switch iThumbnail{1}
        case '1' %rotated image
            iImg = userData.ImgRotated(1).cdata;
            updateImgDisplay(iAx, iImg, iCallback, userData);
            
        case '2' %image variants
            iImg = userData.selectedVariant{1};
            updateImgDisplay(iAx, iImg, iCallback, userData);
            
        case '3' %intensity threshold?
            %nothing for right now
            set(iAx, 'Color','k');
            
        case '4' %canny filter
            iImg = userData.imgEdge;
            
            iCallback = ''; %no callback for now
            updateImgDisplay(iAx, iImg, iCallback, userData);
            
        case '5' %image ROI
            iMask = userData.imgRoi;
            iImg = userData.selectedVariant{1};
            updateImgDisplay(iAx, iImg, iCallback, userData);
            hImg = getappdata(iAx, 'hImg');
            maxAlpha = 255; set(iAx, 'ALim',[0 maxAlpha]);
            transparency = maxAlpha*0.2;
            adata = ones(size(iImg))*transparency; adata(iMask) = maxAlpha;
            set(hImg, 'AlphaData', adata, 'AlphaDataMapping', 'direct')
            
            
            plotRoiRegion(iAx, iMask);
        
            
        case '6' %show edge points
            
            iImg = userData.selectedVariant{1};
            updateImgDisplay(iAx, iImg, '', userData); %iCallback
            edgePoints = userData.edgePoints;
            %show edge points
            
                
            hold(iAx, 'on');
            if ~isempty(edgePoints) 
                plot(iAx, edgePoints(:,1), edgePoints(:,2), '.m');
            end
           
            hold(iAx, 'off');
            
%             updateImgDisplay(iAx, userData.derMask{iSurfNo}, iCallback);
            
        case '7' %cosine fit
            
            
            if isShowMasks
                iImg = userData.imgEdge;
            else
                iImg = userData.selectedVariant{1};
            end
            
            
            if ~isfield(userData, 'cosFit_orig')
                iThumbnail(1) = [];
                continue
            end
            
            %show the image in mm
            updateImgDisplay_mm(iAx, iImg, '', userData);  
            
            
            %unshifted points
            cosFit_orig = userData.cosFit_orig; 
            segPoints = cosFit_orig.segPoints;
            xyFit = cosFit_orig.xyFit;     
          
            hold(iAx, 'on')
            %show the points            
            if ~isempty(segPoints) 
                plot(iAx, segPoints(:,1), segPoints(:,2), '.m','MarkerSize',0.5);
            end
            
            %show the fit
            if ~isempty(xyFit) 
                plot(iAx, xyFit(:,1), xyFit(:,2),'LineWidth',2);
            end
            
            hold(iAx, 'off');
            
        case '8' %window select
            
            if ~isfield(userData, 'cosFit_orig')
                iThumbnail(1) = [];
                continue
            end
            
            %show the image in mm
            if isShowMasks
                iImg = userData.imgEdge;
            else
                iImg = userData.selectedVariant{1};
            end
            
            %this function clears the axis first
            updateImgDisplay_mm(iAx, iImg, iCallback, userData);  
            
            %display all the points within the window
            %show: all points, inlierPts
            segPts_all = userData.cosFit_orig.segPoints;
            isInWindow = userData.cosFit_corr.isInWindow{1}; %not a mistake, use corrected
            segPts = segPts_all(isInWindow, :);
            
            
            hold(iAx, 'on');
            hPoints(1) = plot(iAx, segPts_all(:,1), segPts_all(:,2),'.r');
            hPoints(2) = plot(iAx, segPts(:,1), segPts(:,2),'.b');
            
            setappdata(iAx, 'hPoints', hPoints); %just encase you need to delete
            
            hold(iAx, 'off');
            %--------------------------------------------------------------
            %OBSOLETE CODE
            %--------------------------------------------------------------
            
            iThumbnail(1) = [];
            continue %for now
            
            %fit with tolerance range
            updateImgDisplay(iAx, userData.img, iCallback);
            
            allPoints = userData.allPoints{iSurfNo};
            inlierPoints = userData.segPoints{iSurfNo};
            xyFit = userData.finalFit{iSurfNo};
            if any([isempty(allPoints), isempty(inlierPoints), isempty(xyFit)])
                iThumbnail(1) = [];
                continue
            end
            
            params = userData.params;
           
            
            plotFitTolerance(iAx, allPoints, inlierPoints, xyFit, tolRange_pxls);
        case '9' %recursive fit
            
            if ~isfield(userData, 'quadFit')
                iThumbnail(1) = [];
                continue
            end
            
            %show the image in mm
            if isShowMasks
                iImg = userData.imgEdge;
            else
                iImg = userData.selectedVariant{1};
            end
            
            %this function clears the axis first
            updateImgDisplay_mm(iAx, iImg, iCallback, userData);  
            
            %show the original unshifted segmentation points
            fitNo = 1;
            
            %original points
            allPts_ant = userData.quadFit(fitNo).inputXY{1};
            allPts_pos = userData.quadFit(fitNo).inputXY{2};
            
            %inliers
            inlierPts_ant = userData.quadFit(fitNo).inlierXY{1};
            inlierPts_pos = userData.quadFit(fitNo).inlierXY{2};
            
            %fit
            inlierFit_ant = userData.quadFit(fitNo).inlierFit{1};
            inlierFit_pos = userData.quadFit(fitNo).inlierFit{2};
            
            params = userData.params;
            tolRange_mm = (params.startFitTol)...
                ./(params.fitDivFactor*params.nRecursions);
            
            %anterior
            plotFitTolerance(iAx, allPts_ant, inlierPts_ant,...
                inlierFit_ant, tolRange_mm);
            
            %posterior
            plotFitTolerance(iAx, allPts_pos, inlierPts_pos,...
                inlierFit_pos, tolRange_mm);
            
        case '10' %conic Fit
            
            if ~isfield(userData, 'conicFit')
                iThumbnail(1) = [];
                continue
            end
            
            %show the image in mm
            if isShowMasks
                iImg = userData.imgEdge;
            else
                iImg = userData.selectedVariant{1};
            end
            
            %this function clears the axis first
            updateImgDisplay_mm(iAx, iImg, '', userData);
            
            %show the points and fit
            fitNo = 1;
            
            %points
            ant_xyPts = userData.conicFit(fitNo).inputXY{1};
            pos_xyPts = userData.conicFit(fitNo).inputXY{2};
            
            %fit
            ant_xyFit = userData.conicFit(fitNo).inlierFit{1};
            pos_xyFit = userData.conicFit(fitNo).inlierFit{2};
            
            
            hPoints = [];
            hold(iAx, 'on'); 
            
            %plot points
            hPoints(1) = plot(iAx, ant_xyPts(:,1), ant_xyPts(:,2),'.r');
            hPoints(2) = plot(iAx, pos_xyPts(:,1), pos_xyPts(:,2),'.r');
            
            %plot fit
            hPoints(3) = plot(iAx, ant_xyFit(:,1), ant_xyFit(:,2),'b',...
                'LineWidth', 2);
            hPoints(4) = plot(iAx, pos_xyFit(:,1), pos_xyFit(:,2),'b',...
                'LineWidth',2);
            
            hold(iAx, 'off');
            
            setappdata(iAx, 'hPoints', hPoints);
            
            allR = cat(2,userData.conicFit.R)'; %{ant1, pos1;ant2,pos2}
            allP = cat(2,userData.conicFit.p)';
            tableData = [allR(:,1), allP(:,1), allR(:,2), allP(:,2)];
            set(handles.table_fitVals,'Data', tableData);
            
    end
    iThumbnail(1) = [];
end

resetAxesTitles(handles);

function resetAxesTitles(handles)
%Set the title on every axis
titleInfo = {'axes1', 'Rotate Image';...
    'axes2', 'Select Image';...
    'axes3', 'Intensity Threshold';...
    'axes4', 'Edge Filter';...
    'axes5', 'Image ROI';...
    'axes6', 'Select Edges';...
    'axes7', 'Straighten Lens';...
    'axes8', 'Select Window';...
    'axes9', 'Remove Outliers';...
    'axes10', 'Conic Fit'};

arrayfun(@(i) title(handles.(titleInfo{i,1}), titleInfo{i,2},...
    'FontWeight', 'bold'), 1:size(titleInfo,1));

function userData = initialize_exVivoLens(userData)

imgIn = userData.imgIn;

%Rotate the Images
ImgRotated(1).label = '0 Degrees';
ImgRotated(1).cdata = imgIn;
ImgRotated(1).nRotations = 0;

ImgRotated(2).label = '90 Degrees';
ImgRotated(2).cdata = imrotate(imgIn, 90);
ImgRotated(2).nRotations = 1;

ImgRotated(3).label = '180 Degrees';
ImgRotated(3).cdata = imrotate(imgIn, 180);
ImgRotated(3).nRotations = 2;

ImgRotated(4).label = '270 Degrees';
ImgRotated(4).cdata = imrotate(imgIn, 270);
ImgRotated(4).nRotations = 3;

%Image Select Options

if size(imgIn, 3) > 1
    imgVariantSettings = {'Grayscale Image', @(I) rgb2gray(I);...
        'Red Channel (RC)', @(I) I(:,:,1);...
        'Green Channel (GC)', @(I) I(:,:,2);...
        'Blue Channel (BC)', @(I) I(:,:,3);...
        'RC - GC', @(I) I(:,:,1) - I(:,:,2);...
        'RC - BC', @(I) I(:,:,1) - I(:,:,3);...
        'GC - BC', @(I) I(:,:,2) - I(:,:,3)};
    
    for i = 1:size(imgVariantSettings,1)
        iSettings = imgVariantSettings(i,:);
        ImgVariants(i).label = iSettings{1};
        ImgVariants(i).fcn_cdata = iSettings{2};
    end
    
    userData.hasVariants = true;
    userData.ImgVariants = ImgVariants;
else
    userData.hasVariants = false;
end

userData.ImgRotated = ImgRotated;



%==========================================================================
%MOUSE HOVER
%==========================================================================

function mouseMoveCallback(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Let me know how long you've been hovering

tObj = getappdata(hObject, 'timerObj');

isThumbPanelVisible = strcmpi(get(handles.panel_thumbnails, 'Visible'),'on');

if isempty(tObj)
    return
end

%Restart Timer
if ~isThumbPanelVisible
    stop(tObj); set(handles.panel_hoverAxes, 'Visible','off');
elseif ~strcmpi(tObj.Running,'on')
    stop(tObj); start(tObj);
    set(handles.panel_hoverAxes, 'Visible','off'); %was: handles.axes_hover
end
%..........................................................................

function hoverCallback(tObj,edata, handles)
%This is called for a timer function; edata -->edata.Type, edata.Data.time

% disp('Timer Event');

%Find out if you're hovering an axis, and if yes, display the popup axis

%Find out if I'm hovering anything
hHoverAx = handles.axes_hover; hFig = ancestor(hHoverAx, 'figure');
userData = getappdata(hFig, 'userData');
hThumbnailAx = getIntersectingThumbnail(handles);

if ~isempty(hThumbnailAx)
    %Move the hovering axis
    repositionThumbnail(handles, hThumbnailAx)
    
    %clear axes and copy objects
    thumbTag = get(hThumbnailAx,'Tag');
    
    if isempty(thumbTag)
        thumbTag = getappdata(hThumbnailAx, 'lastValidTag');
    end
    
    updateThumbnail(handles, userData, {thumbTag(5:end)}, hHoverAx);
    set(hHoverAx, 'XTick',[], 'YTick',[],'Visible','on');
    
    thumbTitleStr = get(get(hThumbnailAx, 'Title'),'String');
    if isempty(thumbTitleStr) && hThumbnailAx ~= handles.axes10
        thumbTitleStr = get(get(handles.([thumbTag(1:5),'a']), 'Title'),'String');
    end
    
    title(hHoverAx, thumbTitleStr, 'Color','w', 'FontWeight', 'bold');
    
end

function repositionThumbnail(handles, hThumb)

hHoverAx = handles.panel_hoverAxes; %, handles.axes_hover];
hFig = ancestor(hHoverAx, 'figure');
hPanel= handles.panel_thumbnails;

allObjects = [hFig, hPanel, hHoverAx, hThumb];
originalUnits = get(allObjects, 'Units');

set(allObjects, 'Units','pixels');
panelPos = get(hPanel, 'Position');
thumbPos = get(hThumb, 'Position');

if strcmpi(getappdata(hThumb,'lastValidTag'),'axes10')
    thumbPos_absolute = thumbPos;
else
    thumbPos_absolute = [thumbPos(1:2) + panelPos(1:2), thumbPos(3:4)];
end

%left, bottom, right, top
thumbBounds = [thumbPos_absolute(1:2),thumbPos_absolute(1:2) + thumbPos_absolute(3:4)];
offsetAmount = 20; %80pixles = 10 characters

hoverPos = get(hHoverAx, 'Position');

%Determine where to put the hovering axis
%pg 185, OBC notebook 1
thumbX = thumbBounds(1); thumbY = thumbBounds(2);
leftOffset = thumbX - (hoverPos(3) + offsetAmount);
rightOffset = thumbX + (thumbPos_absolute(3) + offsetAmount);
topAlign = thumbY + (thumbPos_absolute(4) - hoverPos(4));
middleAlign = thumbY + (0.5*(thumbPos_absolute(4) - hoverPos(4)));
bottomAlign = thumbY;

newPos = hoverPos;
switch getappdata(hThumb,'lastValidTag')
    case {'axes1', 'axes2'} %align with top of axis, right offset
        newPos(1:2) = [rightOffset, topAlign];
    case {'axes3'} %align with top of axis, left offset
        newPos(1:2) = [leftOffset, topAlign];
    case {'axes4', 'axes5'} %align with middle of axis, right offset
        newPos(1:2) = [rightOffset, middleAlign];
    case {'axes6','axes10'} %align with middle of axis, left offset
        newPos(1:2) = [leftOffset, middleAlign];
    case {'axes7','axes8'} %align with bottom of axis, right offset
        newPos(1:2) = [rightOffset, bottomAlign];
    case {'OBSOLETE-bottom left'} %align with bottom of axis, left offset
        newPos(1:2) = [leftOffset, bottomAlign];
    case {'axes9'}
        newPos(1:2) = [leftOffset, bottomAlign];
end
set(hHoverAx, 'Position', newPos, 'Visible','on');
arrayfun(@(i) set(allObjects(i), 'Units', originalUnits{i}),1:numel(allObjects));

function axOut = getIntersectingThumbnail(handles)
hFig = ancestor(handles.axes_hover,'figure');
axesIDs = {'1','2','3','4','5','6','7','8','9','10'};
axesTags = cellfun(@(iID) ['axes', iID], axesIDs, 'UniformOutput', false);

axisHandles = force1D(cellfun(@(iTag) handles.(iTag), axesTags,...
    'UniformOutput', false));
axisHandles = cat(1, axisHandles{:});
axisHandles(~ishandle(axisHandles)) = [];

oldAxUnits = get(axisHandles, 'Units');
oldFigUnits = get(hFig, 'Units');
set([hFig;axisHandles],'Units', 'pixels');

%NOTE: THIS DOESN'T TAKE INTO ACCOUNT PANEL OFFSET

%[left, bottom, width height]
axisPositions = get(axisHandles, 'Position');
axisPositions = cat(1, axisPositions{:});

%[left, bottom, right, top]
axisBounds = axisPositions;
axisBounds(:,3:4) = axisPositions(:, 1:2) + axisPositions(:,3:4);

cursorPos = get(hFig, 'CurrentPoint');

isHovering = sum([cursorPos(1)>axisBounds(:,1),... %leftBound vs x
    cursorPos(2)>axisBounds(:,2),... %bottomBound vs y
    cursorPos(1)<axisBounds(:,3),... %rightBound vs x
    cursorPos(2)<axisBounds(:,4)], 2)==4;... %rightBound vs x
    
axOut = axisHandles(isHovering);
%RESET UNITS
set(axisHandles,'Units',oldAxUnits{1});
set(hFig, 'Units', oldFigUnits)


%==========================================================================
%THUMBNAIL CALLBACKS
%==========================================================================

function thumbnailCallback(handles, hThumbAx, imgNo)
%hAx --> thumbnail axes corresponding to the callback
%imgNo --> 1 = anterior; 2 = posterior

hFig = ancestor(hThumbAx,'figure');
userData = getappdata(hFig, 'userData');
hRadioMasks = handles.radio_showMasks;
isShowMasks = get(hRadioMasks , 'Value');
set(handles.radio_showImgs, 'Value',1);

hAxTag = getappdata(hThumbAx, 'lastValidTag');
panelTag = ['panel_', hAxTag(1:5)];

if isfield(handles, panelTag)
    hActivePanel = handles.(panelTag);
else
    hActivePanel = handles.panel_axesX;
end

configurePanel(handles, hActivePanel, @(varargin) updatePanelDisplay(handles, hActivePanel));

hMainPanel  = ancestor(hActivePanel, 'uipanel','toplevel');
%Make only the main panel visible
set([handles.panel_thumbnails,...
    handles.panel_singleAxis, handles.panel_doubleAxis, handles.panel_hoverAxes],'Visible','off')
set([hMainPanel, handles.toolbar_plot], 'Visible','on');

hDoneButton = findobj(hMainPanel,'Style','pushbutton','String','Done');
setappdata(hDoneButton,'userData',userData);
setappdata(hDoneButton,'hThumbnail',hThumbAx);
setappdata(hDoneButton, 'panelFlag', false);
setappdata(hMainPanel, 'imgNo', imgNo);
setappdata(hMainPanel, 'hDone', hDoneButton);
setappdata(hActivePanel, 'hThumbnail', hThumbAx);

iAx = [handles.axes_settings2a, handles.axes_settings2b];


switch hAxTag
    case 'bZoom'
        iAx = handles.axes_settings1;
        updateThumbnail(handles, userData, '2', iAx); %show variantImg
        setappdata(hDoneButton,'panelFlag', true); %always update, even if no change
        
    case {'axes1', 'axes2', 'axes4'} %single
        iAx = handles.axes_settings1;
        updateThumbnail(handles, userData, hAxTag(5:end), iAx);
    case 'axes3' %artifact removal
        %show the filtered image with the artifacts
        artifactTol = userData.params.artifactTol;
        updateImgDisplay(iAx(1), userData.imgF,'');
        usedArtifacts = plotArtifacts(handles.listbox_artifactIdx, iAx(1));
        
        %show the modified image
        [imgOut, artifactMask] = removeArtifact(userData.imgF, usedArtifacts, artifactTol);
        updateImgDisplay(iAx(2), imgOut,'');
        
        %Outputs
        setappdata(hActivePanel, 'artifactMask', artifactMask);
        setappdata(hActivePanel, 'img', imgOut);
        
    case 'axes5' %single axis
        iAx = handles.axes_settings1;
        setappdata(hActivePanel, 'hAx', iAx);
        updatePanelDisplay(handles, hActivePanel);
        
        
    case 'axesNONE' %int mask --> double axis: panelX
        %show the original
        imgF = userData.imgF; thresh = userData.intThresh;
        filterSize = userData.params.intFilterSize;
        setappdata(hActivePanel, 'img', imgF);
        
        %show the mask
        updateImgDisplay(iAx(1), imgF,'');
        initializePanelX(handles, thresh, filterSize, imgF);
        updatePanelDisplay(handles, hActivePanel);
        
        set(handles.title_panelX, 'String', 'Int. Threshold Settings');
    case 'axes6' %axes INACTIVE% %double: X + a/b
        %show the original
        img = userData.imgDer{imgNo}; thresh = userData.derThresh(imgNo);
        filterSize = userData.params.derFilterSize(imgNo,:);
        setappdata(hActivePanel, 'img', img);
        
        %show the mask
        updateImgDisplay(iAx(1), img,'');
        initializePanelX(handles, thresh, filterSize, img);
        updatePanelDisplay(handles, hActivePanel);
        
        switch imgNo
            case 1
                set(handles.title_panelX, 'String', 'Ant. Derivative Settings');
            case 2
                set(handles.title_panelX, 'String', 'Pos. Derivative Settings');
        end
        
    case 'axes7' %axes INACTIVE% %single: a/b
        iAx = handles.axes_settings1;
        updateImgDisplay(iAx, userData.img,'');
        
        allGroups = userData.allGroups{imgNo};
        
        gIdx = 1:numel(allGroups); gNo = gIdx([allGroups.isBest]);
        [graphLabels, hPlots] = plotGroups(allGroups, iAx, gNo);
        hListbox = handles.listbox_selectBestGroup;
        
        %populate select group listbox
        set(hListbox,'String', graphLabels,'Value', gNo,'Min',1, 'Max',1);
        setappdata(hListbox, 'hGroupPlots', hPlots);
        setappdata(hListbox, 'gNo', gNo);
        
    case {'axes8', 'axes9'} %single: a/b
        updateThumbnail(handles, userData, hAxTag(5:end), iAx(1));
end

setappdata(hActivePanel, 'hAx', iAx);
set(hRadioMasks , 'Value', isShowMasks);
%--------------------------------------------------------------------------
function initalizeSettingValues(handles, userData)
%This function sets values for all sliders, edit boxes etc

params = userData.params;

%PANEL1: Rotate Image
%nothing 

%PANEL2: Select Image Variant
if userData.hasVariants
    variantLabels = {userData.ImgVariants.label}; 
    variantIdx = userData.variantIdx;
    set(handles.popup_variantSelect, 'String', variantLabels, 'Value',...
        variantIdx);
else
    %do nothing
end


% %PANEL8: Recursive Settings
%nRecursions
set([handles.slider_nRecursion, handles.edit_nRecursion],'Value', params.nRecursions,...
    'String',params.nRecursions);


%Fit Tolerance
set(handles.edit_startWidth, 'String', params.startFitTol, 'Value', params.startFitTol);
set(handles.edit_divFactor, 'String', params.fitDivFactor, 'Value', params.fitDivFactor);


%PANEL9: Fit Window
set([handles.edit_fitWindow, handles.slider_fitWindow], 'String', params.fitWindow,...
    'Value', params.fitWindow);

%PANEL3: Artifact Idx
%roi Tol
% set([handles.slider_artifactTol, handles.edit_artifactTol],...
%     'Value', params.artifactTol, 'String',params.artifactTol);
% 
% %artifact listbox
% hListbox = handles.listbox_artifactIdx; allArtifacts = userData.allArtifactIdx;
% listboxStr = arrayfun(@num2str, allArtifacts,'UniformOutput',false);
% nArtifacts = numel(listboxStr); listboxIdx = 1:nArtifacts;
% isSelected = arrayfun(@(iIdx) any(iIdx == allArtifacts), userData.artifactIdx);
% set(hListbox, 'String',listboxStr, 'Min',0, 'Max',nArtifacts, 'Value',...
%     listboxIdx(isSelected));
% setappdata(hListbox, 'allArtifactIdx', allArtifacts);
% 
% %PANEL7: Select Best Group
% %do nothing
% 
% 


params = userData.params; %this is just to contain the code

function initializePanelX(handles, thresh, filterSize, img)

set([handles.edit_thresholdX, handles.slider_thresholdX],'Value', thresh,...
    'String',thresh);

set(handles.edit_medianFilterX,'Value', filterSize(2),'String',filterSize(2));
set(handles.edit_medianFilterY,'Value', filterSize(1), 'String',filterSize(1));

minThresh = min(img(:));
[counts, grayVal]= imhist(uint8(img));
grayIdx = find(counts > 100,1,'last');
maxThresh = grayVal(grayIdx);

stepSize = 1/double(maxThresh-minThresh);
set(handles.slider_thresholdX, 'Min', minThresh, 'Max', maxThresh,...
    'SliderStep',[stepSize, stepSize]);

function configurePanel(handles, hActiveSubPanel, hPanelCallback)

%Set all controls to inactive state
hParentPanel = get(hActiveSubPanel, 'Parent');

hPanels = findobj(hParentPanel, 'Type','uipanel');
hTitles = findobj(hParentPanel, 'Type','uicontrol','Style','text',...
    '-regexp', 'Tag','title*');
hControls = findobj(hParentPanel, 'Type','uicontrol','-not','Style','text');
hLabels = findobj(hParentPanel, 'Type','uicontrol','-regexp','Tag','label*');
hSliderLabels = findobj(hParentPanel, 'Style','text','-regexp','Tag','edit_*');

uData= getappdata(ancestor(hParentPanel,'figure'),'userData');
default = uData.default; col = default.colors;


set([hPanels;hLabels; hTitles], 'BackgroundColor', col.inactivePanel,'ForegroundColor',...
    col.inactiveText);
set(hTitles, 'ForegroundColor', col.inactiveTitle);
set(hControls, 'Enable', 'off');
set(hSliderLabels,'BackgroundColor', col.inactiveText);

%Set active panel
hActiveTitle = findobj(hActiveSubPanel, 'Type','uicontrol','Style','text',...
    '-regexp', 'Tag','title*');
hActiveControls = findobj(hActiveSubPanel, 'Type','uicontrol','-not','Style','text');
hActiveLabels = findobj(hActiveSubPanel, 'Type','uicontrol','-regexp','Tag','label*');
hActiveSliderLabels = findobj(hActiveSubPanel, 'Style','text','-and','-regexp','Tag','edit_*');
hActiveSliderAndEdits = findobj(hActiveSubPanel, 'Style','edit','-or','Style','Slider');

set([hActiveSubPanel;hActiveLabels; hActiveTitle], 'BackgroundColor', col.activePanel,'ForegroundColor',...
    col.activeText);
set(hActiveTitle, 'ForegroundColor', col.activeTitle);
set(hActiveControls, 'Enable', 'on');
set(hActiveSliderLabels,'BackgroundColor', col.activeText);

if nargin > 2
    %remove these when the user clicks ok
    setappdata(hParentPanel,'hEditControls', hActiveControls);
    setappdata(hParentPanel,'hPanelCallback', hPanelCallback);
    
    if ~isempty(hActiveControls)
        
        hasCallback= arrayfun(@(iHandle) ~isempty(getappdata(iHandle, 'callbackID')),...
            hActiveControls);
        
        %MODIFYING hActiveControls, DO NOTTTT use after the if statement
        hActiveControls = hActiveControls(~hasCallback);
        
        %add callbacks if necessary
        callBackIds = arrayfun(@(iHandle) iptaddcallback(iHandle,'Callback',hPanelCallback),...
            hActiveControls);
        arrayfun(@(i) setappdata(hActiveControls(i),'callbackID', callBackIds(i)),...
            1:numel(callBackIds));
    end
end
set(hParentPanel, 'BackgroundColor', get(get(hParentPanel,'Parent'), 'BackgroundColor'))

function updatePanelDisplay(handles, hPanel)

panelName = get(hPanel, 'Tag');
hMainPanel = ancestor(hPanel,'uipanel','toplevel'); %panel_doubleAxis
imgNo = getappdata(hMainPanel, 'imgNo');

if isempty(panelName)
    panelName = getappdata(hPanel, 'lastValidTag');
end

hAx = getappdata(hPanel, 'hAx');
hMainPanel = ancestor(hPanel,'uipanel','toplevel');
hDoneButton = findobj(hMainPanel,'Style','pushbutton','String','Done');
setappdata(hDoneButton, 'panelFlag', true);
userData = getappdata(hDoneButton,'userData');


switch panelName
    
    case {'panel_axes1', 'panel_axes2'} %rotate image & selected image variant
        
        cla(hAx); %make sure you clear it before rotating
        updateThumbnail(handles, userData, panelName(end), hAx(1));
        
    case 'OBSOLETE - panel_axes3' %intensity thresh? code UNTOUCHED for now
        
        %show the filtered image with the artifacts
        artifactTol = get(handles.slider_artifactTol, 'Value');
        usedArtifacts = plotArtifacts(handles.listbox_artifactIdx, hAx(1));
        
        %show the modified image
        [imgOut, artifactMask] = removeArtifact(userData.imgF, usedArtifacts, artifactTol);
        updateImgDisplay(hAx(2), imgOut,'');
        
        %Recalculate image derivative
        imgDer = getDerivative(double(imgOut) - mean(imgOut(:)), 3);
        [posDer, negDer] = deal(imgDer);
        negDer(negDer>0) = 0; posDer(posDer<0) = 0;
        
        
        userData.imgDer{1} = abs(posDer);
        userData.imgDer{2} = abs(negDer);
        
        userData.artifactTol = artifactTol;
        userData.artifactMask = artifactMask;
        userData.artifactIdx = usedArtifacts;
        userData.img = imgOut;
        
    case 'panel_axes5' %ROI
        %if user chooses to display both
        hAx = handles.axes_settings1;
        iImg = userData.selectedVariant{1};
        iRoi = userData.imgRoi;
        iEdge = userData.imgEdge;
        
        isShowEdge = false;
        switch get(handles.popup_RoiDisplay,'Value')
            case 1 %Image Only
                plotImg = iImg;
            case 2 %Edge Only
                plotImg = iEdge;
            case 3 %Image+Edge
                plotImg = iImg;
                isShowEdge = true;

        end
        
        updateImgDisplay(hAx, plotImg);
        hImg = getappdata(hAx, 'hImg');
        maxAlpha = 255; set(hAx, 'ALim',[0 maxAlpha]);
        transparency = maxAlpha*0.1;
        adata = ones(size(iImg))*transparency; adata(iRoi) = maxAlpha;
        set(hImg, 'AlphaData', adata, 'AlphaDataMapping', 'direct')
        
        plotRoiRegion(hAx, iRoi);
        
        if isShowEdge
            edgeMask = iEdge; edgeMask(~iRoi) = false;
            plotEdgePoints(hAx, edgeMask);
        else
            hEdges = getappdata(hAx, 'hEdges');
            if ~isempty(hEdges) && ishandle(hEdges)
                delete(hEdges);
            end
        end
        
    case 'panel_axes8'
        userData.params.fitWindow = get(handles.edit_fitWindow,'Value');
        userData = updateWindow(userData);
        updateThumbnail(handles, userData, '8', hAx(1));  
        
    case 'panel_axes9' %recursion
        %get the recursion plots
        %replot accepted and discarded
        
        userData.params.startFitTol = get(handles.edit_startWidth,'Value');
        userData.params.fitDivFactor = get(handles.edit_divFactor,'Value');
        userData.params.nRecursions = get(handles.slider_nRecursion,'Value');
        
        %show both axes
        %         dbs= dbstack;
        %         disp({dbs.name});
        %         disp('................................');
        
        %re-run quadfit        
        quadFit = runStateNine(userData);
        userData.quadFit = quadFit;
        updateThumbnail(handles, userData, '9', hAx(1));        
        
    case 'OBSOLETE - panel_axesX'
        img = getappdata(hPanel, 'img');
        thresh = get(handles.slider_thresholdX,'Value');
        threshMin = get(handles.slider_thresholdX, 'Min'); 
        threshMax = get(handles.slider_thresholdX, 'Max');
        if thresh < threshMin
            set(handles.slider_thresholdX, 'Value', threshMin);
        elseif thresh > threshMax
            set(handles.slider_thresholdX, 'Value', threshMax);
        end
        filterSize = [get(handles.edit_medianFilterY,'Value'),...
            get(handles.edit_medianFilterX,'Value')];
        hThumbAx = getappdata(hPanel, 'hThumbnail');
        imgMask = medfilt2(img > thresh, filterSize);
        imagesc(imgMask, 'Parent',handles.axes_settings2b); 
        grayMap(handles.axes_settings2b);
        set(handles.axes_settings2b, 'XTick',[],'YTick',[]);
        
        %update userData structure
        if strcmpi(getappdata(hThumbAx,'lastValidTag'), 'axes5')
            userData.intMask = imgMask;
            userData.intThresh = thresh;
            userData.params.intFilterSize = filterSize;
        else
            userData.derMask{imgNo} = imgMask;
            userData.derThresh(imgNo) = thresh;
            userData.params.derFilterSize(imgNo,:) = filterSize;
        end
end

setappdata(hDoneButton,'userData', userData);


%==========================================================================
%PLOT/IMG DISPLAYS
%==========================================================================

function updateImgDisplay(hAx, img, callback, userData)

hImg = getappdata(hAx, 'hImg');

if isempty(hImg) || ~ishandle(hImg)
    hImg = imagesc(img,'Parent', hAx); axes(hAx),  grayMap(hAx);
    set(hAx, 'XTick',[], 'YTick',[],'Visible','on');
else
    set(hImg, 'CData', img);
end

if nargin > 2
    
    if ~isempty(callback)
        %Reset the callback
        set(hImg, 'ButtonDownFcn', callback);
        
        %when user hovers make it a hand
        iptPointerManager(ancestor(hAx,'figure'));
        iptSetPointerBehavior(hImg, @(f, cp) set(f, 'Pointer', 'hand'));
    end
    
end

%zoom as necessary

hFig = ancestor(hAx, 'figure');
if nargin < 4
    userData = getappdata(hFig,'userData');
end

if ~isempty(userData) && ~isempty(userData.zoomRegion)
    zoom(hFig, 'reset');
    set(hAx, {'XLim', 'YLim'}, userData.zoomRegion);
    drawnow;
end

setappdata(hAx,'hImg', hImg);
setappdata(hAx, 'img', img);

function updateImgDisplay_mm(hAx, img, callback, userData)  

hFig = ancestor(hAx, 'figure');

if nargin < 4
    userData = getappdata(hFig,'userData');    
end

xScale = userData.cosFit_orig.xScale;
yScale = userData.cosFit_orig.yScale;


%show the image
cla(hAx);
hImg = imagesc(xScale, yScale, img,'Parent', hAx);
grayMap(hAx);
set(hAx, 'XTick',[], 'YTick',[],'Visible','on');

%add callback
if nargin > 2
    
    if ~isempty(callback)
        %Reset the callback
        set(hImg, 'ButtonDownFcn', callback);
        
        %when user hovers make it a hand
        iptPointerManager(ancestor(hAx,'figure'));
        iptSetPointerBehavior(hImg, @(f, cp) set(f, 'Pointer', 'hand'));
    end
    
end

%zoom as necessary
if ~isempty(userData) && ~isempty(userData.zoomRegion)
    imgRes = userData.xyRes(1); shiftAmount = userData.cosFit_orig.shiftAmount;
    zoomRegion = userData.zoomRegion; zoomRegion_mm = zoomRegion; 
    zoomRegion_mm{1} = (zoomRegion{1}./imgRes) - shiftAmount(1);
    zoomRegion_mm{2} = (zoomRegion{2}./imgRes) - shiftAmount(2);
    
    zoom(hFig, 'reset');
    set(hAx, {'XLim', 'YLim'}, zoomRegion_mm);
    axis(hAx,'manual'); %to make sure it stays zoomed
    drawnow;
end


%tilt axis
tilt_deg = userData.cosFit_corr.tilt * (180/pi);
[az, el] = view(hAx);
view(hAx, az+tilt_deg,el);

setappdata(hAx,'hImg', hImg);
setappdata(hAx, 'img', img);



function hPlot = plotRoiRegion(hAx, roiMask)


if ~any(roiMask(:))
    return
end
hold(hAx, 'on');

roiBound = bwmorph(roiMask, 'remove');
[yVals, xVals] = find(roiBound);

hPlot = plot(hAx, xVals, yVals, '.m','MarkerSize',1);

hold(hAx, 'off');

setappdata(hAx, 'hROI', hPlot);

function hPlot = plotEdgePoints(hAx, edgeMask)
%This function plots all of the 1's in a logical mask as a point

hold(hAx, 'on');

[yVals, xVals] = find(edgeMask);
hPlot = plot(hAx, xVals, yVals, '.r', 'MarkerSize',2);
hold(hAx, 'off');

setappdata(hAx, 'hEdges', hPlot);

function plotFitTolerance(hAx, allPoints, inlierPoints, xyFit, tolRange_pxls)
%tolRange_pxls --> final tolerance range after all recursions

hPlots = [];

halfTol = tolRange_pxls/2;

hold(hAx, 'on');

xFit = xyFit(:,1); yFit = xyFit(:,2);
hPlots(1) = plot(hAx, allPoints(:,1), allPoints(:,2),'xr'); %all Points
hPlots(2) = plot(hAx, inlierPoints(:,1), inlierPoints(:,2), 'xg'); %inliers
hPlots(3) = plot(hAx, xFit, yFit,'y'); %fit

%Plot tolerance ranges
[upperBound, lowerBound] = deal(yFit);
upperBound = upperBound + halfTol;
lowerBound = lowerBound - halfTol;

hPlots(4) = plot(hAx, xFit, upperBound, '--m'); plot(hAx, xFit, lowerBound, '--m')

hold(hAx, 'off');


hPoints = getappdata(hAx, 'hPoints');

if isempty(hPoints)
    hPoints = hPlots;
else
    hPoints(~ishandle(hPoints)) = [];
    hPoints = [hPoints, hPlots];
end

setappdata(hAx, 'hPoints', hPoints);

function plotLensCoords(hAx, lensCoords)

hOld = getappdata(hAx, 'hLensCooords');

if ~isempty(hOld) && ishandle(hOld)
    delete(hOld);
end

hold(hAx, 'on'); hP = plot(hAx, lensCoords(:,1),lensCoords(:,2),...
    'MarkerSize',6, 'MarkerFaceColor','m','Marker','o',...
    'Color','m','LineStyle','none');
setappdata(hAx, 'hLensCoords',hP);
hold(hAx,'off');

function [graphLabels, hPlot] = plotGroups(groups, hAx, gNo)

if nargin < 2
    figure; hAx = gca;
    gNo = [];
elseif nargin < 3
    gNo = [];
end


nGroups = numel(groups);

if ~isempty(gNo)
    iColor = repmat([1 1 0],[nGroups,1]);
    iColor(gNo,:) = [0 1 0];
else
    iColor = rand(nGroups,3);
end

hold(hAx, 'on');
hPlot = zeros(1,nGroups);
for i = 1:nGroups
    
    if isfield(groups, 'global')
        iData = groups(i).global.inliers; %data;
    elseif isfield(groups,'data')
        iData = groups(i).data;
    else
        continue
    end
    
    if isempty(iData)
        continue
    end
    %     xMin = min(iData(:,1)); xMax= max(iData(:,1));
    %     iAvg = groups(i).avg(2);
    
    hPlot(i) = plot(hAx, iData(:,1), iData(:,2), 'x','Color', iColor(i,:));
end

%
% graphLabels = cellfun(@(n) {sprintf('Group %d',n); sprintf('Group %d',n)}, num2cell(1:numel(groups)),...
%     'UniformOutput',false);

graphLabels = arrayfun(@(n) sprintf('Group %d',n), 1:numel(groups),...
    'UniformOutput',false);

if nargin < 2
    legend(cat(1,graphLabels{:}), 'Location', 'EastOutside');
    axis ij;
end

hold(hAx, 'off')

function usedArtifacts = plotArtifacts(hListbox, hAx)

hOldPlot = getappdata(hAx, 'hArtifact');
if ~isempty(hOldPlot)
    hOldPlot = hOldPlot(ishandle(hOldPlot));
    delete(hOldPlot);
end

allArtifacts = getappdata(hListbox, 'allArtifactIdx');
usedArtifacts = allArtifacts(get(hListbox,'Value'));

%use full artifact range
artifactTol = get(findobj('Tag','slider_artifactTol'), 'Value');
img = getappdata(hAx, 'img'); ySize = size(img,1);
artifactRange = unique(cell2mat(arrayfun(@(x) (x-artifactTol:x+artifactTol)',...
    usedArtifacts,'UniformOutput',false)));

artifactRange(artifactRange<1) = [];
artifactRange(artifactRange>ySize) = [];

xLim = get(hAx, 'XLim')';
artifactPlot = (repmat(artifactRange, [1,2]))';
[xPlot, yPlot] = meshgrid(xLim,artifactPlot);
hold(hAx, 'on'); hPlot = plot(hAx, xPlot', yPlot', 'r'); hold(hAx, 'off');

setappdata(hAx, 'hArtifact', hPlot);



%--------------------------------------------------------------------------
%BUTTON CALLBACKS
%--------------------------------------------------------------------------

% --- Executes on button press in bZoom.
function bZoom_Callback(hObject, eventdata, handles)
% hObject    handle to bZoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%This has nothing to do with this function, I just needed to see it
% userData = getappdata(ancestor(hObject,'figure'),'userData');
% generateTableOutput(userData);

thumbnailCallback(handles, hObject, 1);


function b_rotateImg_Callback(hObject, eventdata, handles)
% hObject    handle to b_rotateLeft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hPanel = ancestor(hObject,'uipanel');

hMainPanel = ancestor(hPanel,'uipanel','toplevel');
hDoneButton = getappdata(hMainPanel, 'hDone');
userData = getappdata(hDoneButton, 'userData');

ImgRotated = userData.ImgRotated;

%shift left
imgRoi = userData.imgRoi; %If you rotate the image, rotate the imgRoi
switch hObject
    case handles.b_rotateLeft %ccw
        shift = -1; 
        userData.imgRoi = imrotate(imgRoi,90);  
    case handles.b_rotateRight %cw
        shift = 1;
        userData.imgRoi = imrotate(imgRoi,-90);  
end

ImgRotated = circshift(ImgRotated, [0,shift]);
    
%DISPLAY: updated in the updatePanelDisplay Function

userData.ImgRotated = ImgRotated;
userData.zoomRegion = []; %reset the zoom

setappdata(hDoneButton, 'userData', userData);
setappdata(hDoneButton, 'panelFlag', true);

function popup_variantSelect_Callback(hObject, eventdata, handles)
% hObject    handle to popup_variantSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hPanel = ancestor(hObject,'uipanel');

hMainPanel = ancestor(hPanel,'uipanel','toplevel');
hDoneButton = getappdata(hMainPanel, 'hDone');
userData = getappdata(hDoneButton, 'userData');

%Update variant image in userData
variantIdx = get(hObject, 'Value');
userData = updateImgVariant(userData, variantIdx);

setappdata(hDoneButton, 'userData', userData);
setappdata(hDoneButton, 'panelFlag', true);

function b_selectIrisROI_Callback(hObject, eventdata, handles)
% hObject    handle to b_selectIrisROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hPanel = ancestor(hObject,'uipanel');
hAx = getappdata(hPanel, 'hAx'); img = getappdata(hAx, 'img');
hLensCoords = getappdata(hAx, 'hLensCoords');

hMainPanel = ancestor(hPanel,'uipanel','toplevel');
userData = getappdata(getappdata(hMainPanel, 'hDone'), 'userData');


if ~isempty(hLensCoords) && ishandle(hLensCoords)
    lensCoords = [get(hLensCoords, 'XData'); get(hLensCoords,'YData')]';
    delete(hLensCoords);
else
    lensCoords = userData.lensCoords;
end


xlims = [1, size(img,2)]; ylims = [1, size(img,1)];
fcn = makeConstrainToRectFcn('impoly', xlims, ylims);
hLine = imline(hAx,lensCoords, 'PositionConstraintFcn',fcn);
wait(hLine);

if hLine.isvalid
    
    lensCoords = round(hLine.getPosition);
    plotLensCoords(hAx,lensCoords);
    
    userData.lensCoords = lensCoords;
    if userData.imgType == 2 || userData.imgType == 4
        userData.fitLimits = lensCoords(:,1);
    end
    setappdata(getappdata(hMainPanel, 'hDone'), 'userData', userData);
    delete(hLine);
end

function b_clearIrisROI_Callback(hObject, eventdata, handles)
% hObject    handle to b_clearIrisROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hPanel = ancestor(hObject,'uipanel');
hAx = getappdata(hPanel, 'hAx');
hLensCoords = getappdata(hAx, 'hLensCoords');

delete(hLensCoords);

function b_clearImgROI_Callback(hObject, eventdata, handles)
% hObject    handle to b_clearImgROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


hPanel = ancestor(hObject,'uipanel');
hAx = getappdata(hPanel, 'hAx'); img = getappdata(hAx, 'img');
hOld = getappdata(hAx, 'hROI');

%just let the user draw a new one
if ~isempty(hOld) && ishandle(hOld)
    delete(hOld);
end

roiMask = true(size(img));
plotRoiRegion(hAx, roiMask);

%Save changes to user data structure
hMainPanel = ancestor(hPanel,'uipanel','toplevel');
hDoneButton = getappdata(hMainPanel, 'hDone');
userData = getappdata(hDoneButton, 'userData');

userData.imgRoi = roiMask;

setappdata(hDoneButton, 'userData', userData);
setappdata(hDoneButton, 'panelFlag', true);

function modImgROI_Callback(hObject, eventdata, handles)
% hObject    handle to b_addImgROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%This function modifies the image ROI

hPanel = ancestor(hObject,'uipanel');
hAx = getappdata(hPanel, 'hAx'); img = getappdata(hAx, 'img');

hMainPanel = ancestor(hPanel,'uipanel','toplevel');
hDoneButton = getappdata(hMainPanel, 'hDone');
userData = getappdata(hDoneButton, 'userData');

oldMask = userData.imgRoi;


%Let the user define the new region
axes(hAx)
currentMask = roipoly();


switch hObject
    case handles.b_addImgROI
        
        if ~all(oldMask(:))
            newMask = oldMask | currentMask;
        else
            newMask = currentMask;
        end
    case handles.b_removeImgROI
        newMask = oldMask;
        newMask(currentMask) = false;
    case handles.b_clearImgROI
        newMask = true(size(img));
end


hROI = getappdata(hAx, 'hROI');
%remove old plot
if ~isempty(hROI) && ishandle(hROI)
    delete(hROI);
end

hEdges = getappdata(hAx, 'hEdges');

if ~isempty(hEdges) && ishandle(hEdges)
    delete(hEdges);
end

% 
% updateImgDisplay(hAx, userData.img); hImg = getappdata(hAx, 'hImg');
% maxAlpha = 255; set(hAx, 'ALim',[0 maxAlpha]); transparency = maxAlpha*0.2;
% adata = ones(size(img))*transparency; adata(newMask) = maxAlpha;
% set(hImg, 'AlphaData', adata, 'AlphaDataMapping', 'direct')
% 
% plotRoiRegion(hAx, newMask);

userData.imgRoi = newMask;

setappdata(hDoneButton, 'userData', userData);
setappdata(hDoneButton, 'panelFlag', true);


function popup_RoiDisplay_Callback(hObject, eventdata, handles)
% hObject    handle to popup_RoiDisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_RoiDisplay contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_RoiDisplay




function b_DoneCallback(hObject, eventdata, handles)
% hObject    handle to b_doneSettings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Get the user data structure, send it back
hFig = ancestor(hObject, 'figure');
userData = getappdata(hObject, 'userData');

set(hFig, 'Pointer', 'watch');

hThumbAx = getappdata(hObject,'hThumbnail');%Update the thumbnail
hAxTag = getappdata(hThumbAx, 'lastValidTag');
startState = str2double(cell2mat(regexp(hAxTag, '\d*','match')));

if isnan(startState) && hThumbAx == handles.bZoom    
    zoomRegion = get(handles.axes_settings1, {'XLim', 'YLim'});
    userData.zoomRegion = zoomRegion;
    startState = 10;
    updateThumbnail(handles, userData,...
        {'1','2','3','4','5','6','7','8','9','10'});
end

%re-run everything

%Determine if any changes were made
if getappdata(hObject, 'panelFlag');
    userData = reprocessSegmentation(handles, userData, startState);
    setappdata(hFig,'userData', userData);
else 
    updateThumbnail(handles, userData, hAxTag(5:end));
end

%Make the thumbnails panel visible
set([handles.panel_singleAxis, handles.panel_doubleAxis], 'Visible','off');
set(handles.panel_thumbnails, 'Visible','on');


%Remove added callbacks
hParentPanel = ancestor(hObject, 'uipanel');
hActiveSliderAndEdits = getappdata(hParentPanel,'hEditControls');
hPanelCallback =  getappdata(hParentPanel,'hPanelCallback');

if ~isempty(hActiveSliderAndEdits)
    iptremovecallback(hActiveSliderAndEdits,'Callback',hPanelCallback);
end

set(hFig, 'Pointer', 'arrow');
turnOffPlotTools(handles);


function b_CancelCallback(hObject, eventdata, handles)
% hObject    handle to b_cancelSettings1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Make the thumbnails panel visible
set([handles.panel_singleAxis, handles.panel_doubleAxis], 'Visible','off');
set(handles.panel_thumbnails, 'Visible','on');
turnOffPlotTools(handles);

% --- Executes on button press in b_exportData.
function b_exportData_Callback(hObject, eventdata, handles)
% hObject    handle to b_exportData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uInput = questdlg('Are you satisfied with the segmentation?',....
    'Segmentation Complete?','Yes','No','Yes');

if ~strcmpi(uInput, 'yes')
    return
end

%Send the layer points for segmentation
%if it's a lens, send the lens coords
hFig = ancestor(hObject, 'figure');
hInvisibleFig = getappdata(hFig, 'hInvisibleFig');
setappdata(hInvisibleFig, 'exportFlag', true);

uiresume(hInvisibleFig);

% Hint: delete(hObject) closes the figure
delete(hObject);






% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hFig = hObject;
hInvisibleFig = getappdata(hFig, 'hInvisibleFig');

if ~isempty(hInvisibleFig) && ishandle(hInvisibleFig) && ~strcmpi(get(hInvisibleFig, 'waitstatus'),'inactive')
    setappdata(hInvisibleFig, 'exportFlag', false);
    uiresume(hInvisibleFig);
else
    if ishandle(hInvisibleFig)
        delete(hInvisibleFig);
    end
    delete(hObject);
end

function turnOffPlotTools(handles)
hFig = ancestor(handles.axes1,'figure');
zoom(hFig, 'off'); pan(hFig, 'off');
dcm_obj = datacursormode(hFig); set(dcm_obj, 'Enable', 'off');

set(handles.toolbar_plot, 'Visible','off');

%--------------------------------------------------------------------------
%EDIT CALLBACKS
%--------------------------------------------------------------------------

function defaultSliderEditCallback(hObject, edata, handles)

objStyle = get(hObject, 'Style');

switch objStyle
    case 'slider'
        objTag = get(hObject, 'Tag');
        
        if isempty(objTag)
            objTag = getappdata(hObject,'lastValidTag');
        end
        objLabel = ['edit',objTag(numel(objStyle)+1:end)];
        
        value = round(get(hObject, 'Value'));
        hAllObj = [hObject, handles.(objLabel)];
    case 'edit'
        
        if (strcmpi(get(hObject, 'Tag'),'edit_startWidth') ||...
                strcmpi(get(hObject, 'Tag'),'edit_divFactor') ||...
                strcmpi(get(hObject,'Tag'),'edit_window'));
            value = str2double(get(hObject, 'String'));
        else
            value = round(str2double(get(hObject, 'String')));
        end
        
        hAllObj = hObject;
end

set(hAllObj, 'String', value, 'Value', value);



%--------------------------------------------------------------------------
%LISTBOX CALLBACKS
%--------------------------------------------------------------------------
function listbox_artifactIdx_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_artifactIdx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function listbox_selectBestGroup_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_selectBestGroup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function panel_showMasks_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in panel_showMasks
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

userData = getappdata(ancestor(hObject, 'figure'), 'userData');
updateThumbnail(handles, userData, {'2','4','5','6','7','8','9','10'});


%--------------------------------------------------------------------------
%CREATE FUNCTIONS
%--------------------------------------------------------------------------

function listbox_selectBestGroup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_selectBestGroup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function slider_DefaultCreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_artifactTol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function listbox_artifactIdx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_artifactIdx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_DefaultCreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_medianFilterY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%--------------------------------------------------------------------------

%==========================================================================
%SEGMENTS - split the lens into segments
%==========================================================================

function allSegmentationPts = extractEdgePoints(edgeMask)
%this function uses the ROI mask to divide the lens into 
occurence = 'first';
img = edgeMask;

[yInd, xInd] = find(edgeMask);
minY = min(yInd); maxY = max(yInd); midY = round(minY+maxY/2);
minX = min(xInd); maxX = max(xInd); midX = round(minX+maxX/2);
roiBounds = [minX, midY;...
    maxX, midY;...
    midX, minY;...
    midX, maxY];


%We will use the x-axis only to find the center coordinate and to divide
%top and bottom halves
xAxis = roiBounds(1:2,:); yAxis = roiBounds(3:4,:);

%Get the horizontal and vertical limits
leftBoundX = xAxis(1,1); rightBoundX = xAxis(2,1);
topBoundY = yAxis(1,2); bottomBoundY = yAxis(2,2);


lensCntrX = round(sum(roiBounds(1:2,1))/2); lensCntrY = round(sum(roiBounds(3:4,2))/2);
imgSize = size(img);
nCols = size(img,2);

%separate into thirds
nSegments = 7;
dx = round( (rightBoundX - leftBoundX) / nSegments );
leftSegmentX = leftBoundX + dx;
rightSegmentX= leftBoundX + dx*(nSegments-1);

%..........................................................................
%DIVIDE THE IMAGE INTO L,R, TOP & BOTTOM HALVES ACCORDING TO THE LENS CNTR
%..........................................................................
idxMat = reshape(1:numel(img), size(img));



%Left Half
params(1).indx = sub2ind(imgSize, idxMat(topBoundY:bottomBoundY,1:leftSegmentX));
params(1).t = @(A) flipud(A'); %ccw 90deg
params(1).tinv = @(A) flipud(A)'; %cw 90deg
params(1).label = 'left';

%Right Half
params(2).indx = sub2ind(imgSize, idxMat(topBoundY:bottomBoundY,rightSegmentX+1:nCols));
params(2).t = @(A) fliplr(A');
params(2).tinv = @(A) fliplr(A)';
params(2).label = 'right';

%Top Half
params(3).indx = sub2ind(imgSize, idxMat(topBoundY:lensCntrY, leftSegmentX:rightSegmentX)); 
params(3).t = @(A) flipud(A);
params(3).tinv = @(A) flipud(A);
params(3).label = 'top';

%Bottom Half
params(4).indx = sub2ind(imgSize, idxMat(lensCntrY+1:bottomBoundY, leftSegmentX:rightSegmentX)); 
params(4).t = @(A) A; % do nothing
params(4).tinv = @(A) A;
params(4).label = 'bottom';



[segmentation, areaImg] = analyzeQuadrant(img, params, occurence);
% 
% figure; imagesc(img); colormap(gray)
% hold on;
% cellfun(@(x) plot(x(:,1), x(:,2),'x', 'Color', rand(1,3)),...
%     segmentation ,'UniformOutput', false);
% legend({'Left', 'Right', 'Top', 'Bottom'});

allSegmentationPts = cat(1,segmentation{:});

function [output, areaImg] = analyzeQuadrant(img, params, occurence)
%This sub-function finds the edges in a given quadrant by transposing it so
%that the surface in the center is facing the top of the image.  It outputs
%the segmentation points in the original coordinate system
%..........................................................................
%INPUTS:
%..........................................................................
%img: a 2D matrix containing the image to be segmented
%
%params: a structure with the following fields:
%   indx - indices of sub image in the original image
%   t - anonymous function that rotates the image
%   tinv - anonymous function that brings image back to original
%   coordinate system
%
%derType: Scalar representing ype of derivative:
%    0 = positive  (first top edge)
%    1 = negative (first bottom edge)
%    2 = both (first top/bottom edge)
%
%..........................................................................
%OUTPUTS:
%..........................................................................
%
% output = cell array same size as input with x and y points of
% segmentation

boolImg = false(size(img));
[nRows, nCols] = size(img);
imgXInd = repmat(1:nCols,[nRows,1]);
imgYInd = repmat((1:nRows)',[1,nCols]);
output = cell(1,numel(params));
areaImg = boolImg;

for i = 1:numel(params)
    iIdx = params(i).indx; iT = params(i).t; iTInv = params(i).tinv;
    iImg = img(iIdx); 
    
    %Rotate the image
    iImg_rot = iT(iImg); [ySize, xSize] = size(iImg_rot);
    
    %Get only unique values
    [iYValue,iXValue] = find (iImg_rot); 
    iXYValues = [iXValue, iYValue];
    iXYValues_sorted = sortrows(iXYValues, [1 2]);
    iXSorted = iXYValues_sorted(:,1);
    
    [~,uniqueIdx,~] = unique(iXSorted,'stable');
    iXYValues_unique = iXYValues_sorted(uniqueIdx,:);
    iXUnique = iXYValues_unique(:,1); iYUnique = iXYValues_unique(:,2);
   
    
    iXInterp = min(iXValue):max(iXValue);
    if numel(iXInterp) < numel(uniqueIdx)
        
        isMissing = arrayfun(@(x) iXInterp == x, iXUnique);
        missingIdx = iXInterp(isMissing);
        notMissingIdx = iXInterp(~isMissing);
        
        extras = numel(missingIdx);
        for j = 1:numel(missingIdx)
            jX= missingIdx(j);
            
            if missingIdx < numel(missingIdx)
                jY = notMissingIdx(missingIdx(j+1));
            else
                jY = notMissingIdx(missingIdx(j-1));
            end
            extras(j,:) = [jX, jY];

        end
        
        iXYValues_unique = sortrows([iXYValues_unique;extras],1);
        iXUnique = iXYValues_unique(:,1); iYUnique = iXYValues_unique(:,2);

%         warning('There are %d values unaccounted for in perimeter calculation',...
%             numel(uniqueIdx) - numel(iXInterp));
    end
    
    iBounds_rot = iXYValues_unique;
    
    %......................................................................
    %AREA CALCULATIONS
    %......................................................................
    %Find points within the lens
    iYIdx = repmat((1:ySize)',[1,xSize]);
    thresholdArray = zeros(1,xSize);
    thresholdArray(iXUnique) = iYUnique;
    thresholdMatrix = repmat(thresholdArray, [ySize,1]);
    isWithinBoundary = iYIdx < thresholdMatrix;
    
    iIsArea_rot = isWithinBoundary;
    iIsArea =  iTInv(iIsArea_rot); %reverse rotation
    areaImg(iIdx) = iIsArea;
    %......................................................................
    
    %Return to the original coordinate system
    iIsImg_rot = false(size(iImg_rot)); %Boolean representative of rotated image
    iIsImg_rot(sub2ind(size(iImg_rot),iBounds_rot(:,2),iBounds_rot(:,1))) = true;

    iIsImg = iTInv(iIsImg_rot); %upright subimg in boolean form
    iFullImg = boolImg; iFullImg(iIdx) = iIsImg; %upright full img in boolean
    
    iXOut = imgXInd(iFullImg); iYOut = imgYInd(iFullImg);
    output{i} = [force1D(iXOut),force1D(iYOut)];
    
end


%--------------------------------------------------------------------------

%==========================================================================
%CORE FUNCTIONS - auto segmentation
%==========================================================================

%MAIN FUNCTIONS
function userData = autoSegment(imgIn, userData, hasSettings)

default = userData.default; 

if hasSettings
    params = userData.params; 
    params.xyRes = userData.xyRes;
else
    params.xyRes = userData.xyRes;
    params.nRecursions = default.nRecursions;
    params.fitType = @quadfit;
    params.fitDivFactor = default.fitDivFactor;
    params.startFitTol = default.fitTol;
    params.fitWindow = default.fitWindow;
end

userData.params = params;
imgType = userData.imgType;
xyRes = userData.xyRes;


%1: Rotate The Image
%pre-process the image (find the ROIs of lens and orientation)
img_raw = userData.ImgRotated(1).cdata;
[lensMask, ~] = estimateShadowgraphROIs(img_raw, false); %isPlot = false


if ~isempty(lensMask) 
    nRotations = getShadowgraphLensOrientation(lensMask);
    
else
    lensMask = true(size(img_raw(:,:,1)));
    nRotations = 0; %don't rotate, you don't know how
    userData.zoomRegion = [];
end


if hasSettings %overwrite rotations if you have it
    nRotations = params.nRotations;
end        

%Rotate the image
if nRotations > 0
    %rotate the image using circshift
    ImgRotated = userData.ImgRotated;
    userData.ImgRotated = circshift(ImgRotated, [0,-nRotations]);
    lensMask = imrotate(lensMask, nRotations*90);
end


%Apply the zoom region
if ~all(lensMask(:)) && (~hasSettings || isempty(userData.zoomRegion))
    %guess zoom region
    props = regionprops(lensMask,'BoundingBox'); tol = 100; %pxls
    bBox = props.BoundingBox;
    xlim = [bBox(1)-tol, bBox(1)+bBox(3)+tol];
    ylim = [bBox(2)-tol, bBox(2)+bBox(4)+tol];
    
    userData.zoomRegion = {xlim, ylim};
end

imgR = userData.ImgRotated(1).cdata; %the first is always the current

%2: Select Image
if userData.hasVariants
    variantIdx = userData.variantIdx;
    fcn_variant = userData.ImgVariants(variantIdx).fcn_cdata;
    variantImg = fcn_variant(imgR);
    variantLabel = userData.ImgVariants(variantIdx).label;
else
    variantImg = imgR;
    variantLabel = 'Original Image';
end

%3: Select Intensity Threshold
% [intMask, ~, intThresh] = extractEdgeMask(variantImg, imgType); %NOT SURE WHY I NEEED THIS....

%4: Canny Filter
[~, thresh]= edge(variantImg, 'canny');
imgEdge = edge(variantImg, 'canny', [3*thresh(1), 3*thresh(2)], 5);


%5: ROIs
%eliminate within the lens as well
if hasSettings && isfield(userData, 'imgRoi') && ~all(userData.imgRoi(:))
    imgRoi = userData.imgRoi; %useful in reloading previous segmentation    
elseif ~all(lensMask(:))
    imgRoi = lensMask & ~imerode(lensMask,strel('disk',75));
else 
    imgRoi = lensMask;
end

%if no Roi, return
if all(imgRoi(:))
    userData.img = variantImg; %almost equivalent to imgIn
    userData.selectedVariant = {variantImg, variantLabel};
    userData.imgEdge = imgEdge;
    userData.imgRoi = imgRoi;
    userData.edgePoints = [];
    userData.params.cannyThresh = thresh;
    userData.params.nRotations = nRotations;
    return
end

%6: Select Edges
edgeMask = imgEdge; 
edgeMask(~imgRoi) = false;
edgePoints = extractEdgePoints(edgeMask);

%7: Cosine Fit
userData.selectedVariant = {variantImg, variantLabel}; %need it to straighten lens
userData.imgRoi = imgRoi; %need it to straighten lens
userData = straightenLens(userData, edgePoints);

%8: Window Select
userData = updateWindow(userData); 

%9: Fit Top and Bottom
%Split the points into anterior, posterior and run a recursive quad fit
quadFit_out = runStateNine(userData);
userData.quadFit = quadFit_out;

%10: Conic Fit
userData = runStateTen(userData);

%update various user data fields
userData.img = variantImg; %almost equivalent to imgIn
userData.selectedVariant = {variantImg, variantLabel};
userData.imgEdge = imgEdge;
userData.imgRoi = imgRoi;
userData.edgePoints = edgePoints;
userData.params.cannyThresh = thresh;
userData.params.nRotations = nRotations;

function userData = reprocessSegmentation(handles, userData, startState)

if startState > 1
    state = getNextState(startState);
else 
    state = 1;
end

while ~isempty(state)
    thumbnailID = num2str(state);
    
    switch state
        
        case {1,3} %filtered image
            %do nothing but update thumbnail
            
        case 2 %select image variant
            %update the variant (i.e. based on rotation)
            userData = updateImgVariant(userData);
            
        
        case 4 %reprocess edge with filter parameters
            thresh = userData.params.cannyThresh;
            variantImg = userData.selectedVariant{1};
            imgEdge = edge(variantImg, 'canny', [3*thresh(1), 3*thresh(2)], 5);
            
            userData.imgEdge = imgEdge;
        
        case 5 %image ROI
            %if you rotated the image, reset the ROI
            currentRoi = userData.imgRoi;
            currentImg = userData.selectedVariant{1};
            
            if ~all(size(currentRoi) == size(currentImg))
               userData.imgRoi = true(size(currentImg)); %reset ROI? or rotate ROI?
            end
            
        case 6 % select edge points
            
            iEdge = userData.imgEdge;
            imgRoi = userData.imgRoi; 
            
            edgeMask = iEdge; edgeMask(~imgRoi) = false;
            userData.edgePoints = extractEdgePoints(edgeMask);
            
        case 7 %cosine fit
            
            userData = straightenLens(userData);
            
        case 8 %update window
            userData = updateWindow(userData);
            
        case 9 %recursion for outliers
            quadFit = runStateNine(userData);
            userData.quadFit = quadFit;
        case 10 %conic Fit
            userData = runStateTen(userData);
        otherwise
            %do nothing, the thumbnail will be updated
    end
    
    state = getNextState(state);    
    updateThumbnail(handles, userData, thumbnailID);
end

function nextState = getNextState(state)
%This currently looks like a redundant diagram but a state diagram will
%help remove unecessary re-processing steps

switch state
    case num2cell(1:9)
        nextState = state+1;        
    otherwise
        nextState = [];
end

function userData = OBSOLETE_runStateSeven(userData, iSurfNo)
%This function reselcts valid pixel points based on image points and
%selects the best group by setting the "isBest" field in allGroups to true
%INPUTS
%imgNo = 1 --> anterior, imgNo = 2--> posterior


stateID = str2double(sprintf('%d',[userData.imgType, iSurfNo]));

artifactMask = userData.artifactMask;
derMask = userData.derMask{iSurfNo}; intMask = userData.intMask;
imgROI = userData.imgRoi;

if ~any(imgROI(:))
    imgROI = true(size(imgROI));
end

basicMask = ~artifactMask & imgROI & derMask & intMask;

switch stateID
    case {11,21,41} %ant cornea, ant lens, ex-vivo ant lens
        validPixelMask = basicMask;
    case {12, 42} %pos cornea, ex-vivo pos lens
        %remove the part close to the previous surface
        xyRes = userData.params.xyRes;
        allGroups = userData.allGroups{1};
        bestGroup = allGroups([allGroups.isBest]);
        
        if ~isempty(bestGroup)
            img = userData.img;
            xyFit = bestGroup.global.fit;
            xFit = xyFit(:,1); yFit = xyFit(:,2);
            
            %Discard everything above anterior surface
            if stateID == 12 %remove .4mm
                dy = 0.4; %mm
            else %%remove 1.2mm
                dy = 0.4; %mm
            end
            
            
            idx2Exclude = segmentMatrix([force1D(xFit), force1D(round(yFit+dy*xyRes(2)))], img);
            posteriorROI = ones(size(img)); posteriorROI(idx2Exclude) = 0;
            
            %discard everything a below maximum distance away
            if stateID == 12
                bottomIdx2Exclude = segmentMatrix([force1D(xFit), force1D(round(yFit+dy*xyRes(2)*6))], img);
                posteriorROI(~bottomIdx2Exclude) = 0;
            end
            validPixelMask = posteriorROI & basicMask;
            
            if ~any(validPixelMask(:))
                validPixelMask = basicMask;
            end
        else
            validPixelMask = basicMask;
        end
        
    case 22 %pos lens (pos when inverted...it's technically anterior lens)
        lensROI = userData.lensRoi;
        validPixelMask = basicMask & lensROI;
    
end


aScanDMasked = userData.imgDer{iSurfNo}; aScanDMasked(~validPixelMask) = 0;

validPoints = findValidPoints(aScanDMasked, false, userData);
if isempty(validPoints)
    return
end

rawAntX = validPoints(:,1); rawAntY = validPoints(:,2);

%BestGroup
isConvex = (iSurfNo == 1) || userData.imgType == 1;
[~, allGroups] = selectBestGroup(rawAntX, rawAntY, @quadfit,...
    isConvex, userData);

userData.allGroups{iSurfNo} = allGroups;
userData.validPixelMask{iSurfNo} = validPixelMask;
userData.allPoints{iSurfNo} = validPoints;

function userData = OBSOLETE_runStateEight(userData, imgNo)

fitLimits = userData.fitLimits;

allGroups = userData.allGroups{imgNo};

if isfield(allGroups, 'isBest')
    bestGroup = allGroups([allGroups.isBest]);
else
    bestGroup = allGroups(1);
end

if isempty(bestGroup)
    xyData = userData.allPoints{imgNo};
else
    xyData = bestGroup.global.inliers;
end

[~, xyFit, xyPoints, ~,~] = runRecursion_mm(userData.img,userData.allPoints{imgNo},...
    xyData, fitLimits, userData.params);
userData.segPoints{imgNo} = xyPoints;
userData.finalFit{imgNo} = xyFit;

function quadFit_out = runStateNine(userData)
%This function splits the points into two curves, and recursively selects
%the inliers.  It mus be ran after the cosine fit & the window is selected
%INPUT: cosFit structure with isInWindow field
%OUTPUT: structure with 4 elements and fields:
	%inputXY - input into the quad fit recursion funcion
    %inlierXY - inler outputs of the recursion function (input to the conicFit)
    %inlierFit - inlier fit [xFit, yFit]
    %p_mm - coefficients of the quadratic fit
    %inputDescripton
        %the 4 elements are: cosFit_ptsOrig, cosFit_fitOrig, cosFit_ptsCor,
            %cosFit_fitCorr

            
set(gcbf, 'Pointer', 'watch')
drawnow;

%..........................................................................
%1: Split the points into two surfaces
%..........................................................................
%I'm still undecided on which points to use, the fit or original points

%corrected points
cosFit_corr = userData.cosFit_corr;
output_corr = splitPts(cosFit_corr); %split the points and return only those in the window
    % DESCRIPTION of output_corr
    % antPts_corr = output_corr{1,1};
    % posPts_corr = output_corr{2,1};
    % antFit_corr = output_corr{1,2};
    % posFit_corr = output_corr{2,2};

%original points
cosFit_orig = userData.cosFit_orig;
output_orig = splitPts(cosFit_orig);
    % DESCRIPTION OF output_orig
    % antPts_orig = output_orig{1,1};
    % posPts_orig = output_orig{2,1};
    % antFit_orig = output_orig{1,2};
    % posFit_orig = output_orig{2,2};

%..........................................................................
%2: Run recursive fit
%..........................................................................

fitWindow = userData.params.fitWindow; dw = fitWindow/2;
xFit = (-dw:0.05:dw)'; %should be the same for all surfaces
params = userData.params;

xyIn = [output_orig, output_corr];
    % xyIn = {antPts_orig, antFit_orig, antPts_corr, antFit_corr;...
    %         posPts_orig, posFit_orig, posPts_corr, posFit_corr};
    
nPairs = size(xyIn,2);
inputLabels = {'cosFit_ptsOrig', 'cosFit_fitOrig',...
    'cosFit_ptsCorr', 'cosFit_fitCorr'};

S = struct();
for iCol = 1:nPairs
    
    iX_ant = xyIn{1,iCol}(:,1); iY_ant = xyIn{1,iCol}(:,2); 
    iX_pos = xyIn{2,iCol}(:,1); iY_pos = xyIn{2,iCol}(:,2); 
    
    [iPts_ant, iFitY_ant, iPmm_ant] = runRecursion_mm(iX_ant, iY_ant,...
        iX_ant, iY_ant, xFit, params);
    
    [iPts_pos, iFitY_pos, iPmm_pos] = runRecursion_mm(iX_pos, iY_pos,...
        iX_pos, iY_pos, xFit, params);
    
    S(iCol).inputXY = {[iX_ant, iY_ant]; [iX_pos, iY_pos]};
    S(iCol).inlierXY = {iPts_ant; iPts_pos};
    S(iCol).inlierFit = {[xFit, iFitY_ant]; [xFit, iFitY_pos]};
    S(iCol).p_mm = {iPmm_ant; iPmm_pos};
    S(iCol).inputDescription = inputLabels{iCol}; 
end

quadFit_out = S;

% your computation
set(gcbf, 'Pointer', 'arrow')
drawnow;

function userData = runStateTen(userData)
%This function runs the conic fit on all members of the quadFit structure
%It is dependent on the function "runStateNine" which does a quadFit
%It looks for userData.quadFit(i).inlierXY



quadFit = userData.quadFit;
S = struct([]);
for i = 1:numel(quadFit)
    iXY = quadFit(i).inlierXY;
    iXFit = quadFit(i).inlierFit{1}(:,1);
    antXY = iXY{1}; posXY = iXY{2};
    antX = antXY(:,1); antY = antXY(:,2);
    posX = posXY(:,1); posY = posXY(:,2);
    
    try
        [iCurve_ant, iGOF_ant] = newConicFit(antX, antY, 3, true);
        [iCurve_pos, iGOF_pos] = newConicFit(posX, posY, 4, true);
    catch
        try
            [iCurve_ant, iGOF_ant] = newConicFit(antX, antY, 5, true);
            [iCurve_pos, iGOF_pos] = newConicFit(posX, posY, 6, true);
        catch
            return
        end
    end
        %3=ant lens; 4=pos lens true - isCorrected
    
	iCoeffVals_ant = coeffvalues(iCurve_ant);
    iCoeffVals_pos = coeffvalues(iCurve_pos);
    
    iYfit_ant = feval(iCurve_ant, iXFit);
    iYfit_pos = feval(iCurve_pos, iXFit);
    
    S(i).inputXY = iXY;
    S(i).inlierFit = {[iXFit, iYfit_ant]; [iXFit, iYfit_pos]};
    S(i).fitObj = {iCurve_ant, iCurve_pos};
    S(i).R = {1/iCoeffVals_ant(1); 1/iCoeffVals_pos(1)};
    S(i).p = {iCoeffVals_ant(2); 1/iCoeffVals_ant(2)};
    S(i).gof = {iGOF_ant; iGOF_pos};
end

userData.conicFit = S;


%--------------------------------------------------------------------------

function validPoints = findValidPoints(aScanDMasked, isFirstPeak, userData)
%if the separation is too small intra layers will be segmented together
switch userData.imgType %modified 3/6/14
    case 1 %invivo cornea
        separationInMM = 2;
    case 2 %invivo lens
        separationInMM = 3;
    otherwise
        separationInMM = 3;
end

u_yRes = userData.xyRes(2);
nSurfaces = 4;
imgCols = size(aScanDMasked,2);

imageGradientCell = num2cell(aScanDMasked, 1);
boundLocation = cellfun(@(iCol) (findMaximums(iCol,nSurfaces,u_yRes,...
    isFirstPeak, separationInMM)),...
    imageGradientCell,'UniformOutput', false);
boundLocation = [boundLocation{:}]';


hasDataArray = ~cellfun(@isempty, boundLocation);
xVals = (1:imgCols)';  isBothZero = sum(hasDataArray,2) == 0;
xVals(isBothZero) = []; boundLocation(isBothZero,:) = [];

isSurfaceZero = sum(hasDataArray,1) == 0;
boundLocation(:,isSurfaceZero) = [];

nSurfacesDetected = size(boundLocation,2);

allData = cell(1,nSurfacesDetected);

for i = 1:nSurfacesDetected
    iY = boundLocation(:,i);
    iX = xVals;
    isValidPartOfImage = ~cellfun(@isempty,iY);
    iX(~isValidPartOfImage) = [];
    allData(i) = {[iX, ([iY{:}]')]};
end

if isempty(allData)
    warndlg('Enlarge your ROI');
    validPoints =[];
    return
end

validPoints = double(sortrows(cat(1,allData{:}),1));

%State 1: FILTER THE IMAGE
function imgOut = filterImage(imgIn, sigma)

filterType = 'gaussian';
kSize = sigma*2;
H = fspecial(filterType, kSize, sigma);
imgOut = imfilter(imgIn, H, 'symmetric');

%State 2: SEGMENT THE IRIS
function lensRoiMask = getLensRoi(img, lensCoords, roiTol)

mc = polyfit(lensCoords(:,1), lensCoords(:,2),1);
xFit = 1:size(img,2);%min(lensX)-5:max(lensX)+5;
yFit = polyval(mc, xFit);
yFit_top = yFit - roiTol; yFit_bottom = yFit + roiTol;

xLeft = max([1, lensCoords(1,1)-round(roiTol/2)]);
xRight = min([lensCoords(2,1)+round(roiTol/2), size(img,2)]);
lensRoiMask = poly2mask([xLeft, xRight, xRight, xLeft, xLeft],...
    [yFit_top(1), yFit_top(end), yFit_bottom(end), yFit_bottom(1), yFit_top(1)],...
    size(img,1), size(img,2));


%State 3: ARTIFACT REMOVAL
function artifactIdxOut = findArtifactIdx(imgIn, userData)

%userData needs: xyRes,
xyRes  = userData.xyRes; u_yRes = xyRes(2);
maxDy_mm = 0.15; %mm
maxDy = maxDy_mm*u_yRes;
imgType = userData.imgType;


%FIND EDGES
img = imadjust(mat2gray(imgIn),[0.16, 1],[0,1], 0.5);
[imgRows, imgCols] = size(img);
yVals = repmat((1:imgRows)', [1,imgCols]);
xVals = repmat((1:imgCols), [imgRows,1]);

[intMask, gMask] = extractEdgeMask(imgIn, imgType);
validPixelMask = bitand(intMask, gMask);

xData = xVals(validPixelMask); yData = yVals(validPixelMask);
allGroups = groupByHistogram(xData, yData, maxDy, false, userData);

if isempty(allGroups)
    artifactIdxOut = [];
    return
end

gNo = 1:numel(allGroups);
gCount = [allGroups.n];

switch imgType
    case 2
        artifactThresh = 15;
        dataThresh = 25;
    otherwise
        artifactThresh = 25;
        dataThresh = 30;
        
end

%make standard deviation and number of members criteria for artifacts
isFirstOrLastGroup = false(size(allGroups)); isFirstOrLastGroup([1 end]) = true;
isArtifact = bitand([allGroups.stdev] < artifactThresh,...
    ([allGroups.hasArtifact] | isFirstOrLastGroup)); %%bitand:
artifactIdx = gNo(isArtifact);
isData = bitand([allGroups.stdev] > dataThresh, [allGroups.n] > mean(gCount)/2);
dataIdx = gNo(isData);


preDataArtifactIdx = max(artifactIdx(artifactIdx<dataIdx(1)));

if ~isempty(preDataArtifactIdx)
    preDataArtifactRange = (1:max(allGroups(preDataArtifactIdx).data(:,2)))';
else
    preDataArtifactRange = [];
end

artifactRange = arrayfun(@(iIdx) (min(allGroups(iIdx).data(:,2)):max(allGroups(iIdx).data(:,2)))',...
    artifactIdx,'UniformOutput', false);

artifactIdxOut  = unique(cat(1, artifactRange{:}, preDataArtifactRange, cat(1,allGroups.artifactIdx)));

function [imgOut, artifactMask] = removeArtifact(imgIn, artifactIdx, artifactTol)

img = imadjust(mat2gray(imgIn),[0.16, 1],[0,1], 0.5); %img = imgIn;

[ySize, ~] = size(imgIn);
artifactMask = false(size(img));

% if isempty(artifactIdx)
%     imgOut = img;
%     return
% end


%artifact + or - 2 pixels changed to +/-6 pixels modified: 11/4/2013
artifactRange = unique(cell2mat(arrayfun(@(x) (x-artifactTol:x+artifactTol)',...
    artifactIdx,'UniformOutput',false)));

artifactRange(artifactRange<1) = [];
artifactRange(artifactRange>ySize) = [];


%modified 12/4/2013 --> added a median filter
imgOutfilt = mat2gray(medfilt2(double(img), [3,3]))*255; %filt was 17x7
%now subtract mean from every row of original img
stdI = std(double(imgOutfilt),0,2);
meanI = mean(imgOutfilt,2);

%lowerMean = meanI - 1.2*mean(stdI);

%modified 2/27/2014
higherMean = mean(imgOutfilt,2) + 2*mean(stdI);
finalMean = repmat(mean(imgOutfilt(:)), size(higherMean));


finalMean(artifactRange) = higherMean(artifactRange);
artifactMask(artifactRange,:) = true;
imgOutfilt = uint8(double(imgOutfilt) - repmat(finalMean, [1, size(imgIn,2)]));
imgOut = imgOutfilt;

%State 4: Select Image ROI


%State 5: Create intensity Mask
%State 6: Create Derivative Mask

%State 7: Select Best Group


%--------------------------------------------------------------------------
%BIOMETRY FUNCTIONS
%--------------------------------------------------------------------------

function analysisOut = generateTableOutput(userData)
%This function creates a table for the Shadowgraph output
%TABLE OUTPUT
    %basicTable
        %SubjectID, ExperimentID, TotalShift, ImageShift, ImageTilt
    %coeffTable with variables:
        %Variable, Value, FitName, FitWindow        
    %fourierTable
        %see legacy table
    %biometryTable
        %legacy table + LensThickness, LensDiameter, LensMask
    %paramTable
        %Param_ImgROI, Param_ZoomRegion, Param_VariantFunction
        %Param_InitialRotation, Param_ImgRes, Param_nRecursions
        %Param_FitDivFactor% Param_StartFitTol, Param_FitWindow, Param_FitType

%..........................................................................
%First, the fits
%..........................................................................
S = struct([]);
Important_coeffs = struct([]);


%Cosine Fits

idx = 1;
cosFit_orig = userData.cosFit_orig;

S(idx).FitName = categorical(cellstr('Cosine'));
S(idx).InputDescription = cellstr('Original Segmentation Points');
S(idx).FitFcn = func2str(@cosine_tiltDecentration);
S(idx).FitWindow = NaN; %mm
S(idx).InputPoints = {cosFit_orig.segPoints};
S(idx).FitPoints = {cosFit_orig.xyFit};
S(idx).FitObject = [];
S(idx).Coefficients = force1D(cosFit_orig.coeff,1);

idx = idx+1;
cosFit_corr = userData.cosFit_corr;

S(idx).FitName = categorical(cellstr('Cosine'));
S(idx).InputDescription = cellstr('Corrected Segmentation Points');
S(idx).FitFcn = func2str(@cosine_tiltDecentration);
S(idx).FitWindow = NaN; %mm
S(idx).InputPoints = {cosFit_corr.segPoints};
S(idx).FitPoints = {cosFit_corr.xyFit};
S(idx).FitObject = [];
S(idx).Coefficients = [];

%Conic fits: only keep fits of original and corrected poitns
if isfield(userData, 'conicFit')
    
    conicFit = userData.conicFit;
    fits2keep = [1,3]; 
    fitsdesc = {'Original Segmentation Points','Corrected Segmentation Points'};
    for j = 1:numel(fits2keep)
        jConicFit = conicFit(fits2keep(j));
        
        idx = idx+1;
        S(idx).FitName = categorical(cellstr('Conic'));
        S(idx).InputDescription = fitsdesc(j);
        S(idx).FitFcn = {formula(jConicFit.fitObj{1}),... %anterior
            formula(jConicFit.fitObj{2})}; %posterior
        S(idx).FitWindow = userData.params.fitWindow; %mm
        S(idx).InputPoints = jConicFit.inputXY';
        S(idx).FitPoints = jConicFit.inlierFit';
        S(idx).FitObject = jConicFit.fitObj;
        S(idx).Coefficients = {coeffvalues(jConicFit.fitObj{1}),... %anterior
            coeffvalues(jConicFit.fitObj{2})}; %posterior
        
        Important_coeffs(end+1).Variable = categorical(cellstr('R'));
        Important_coeffs(end).Value = cell2mat(jConicFit.R');
        Important_coeffs(end).FitName = categorical(cellstr('Conic'));
        Important_coeffs(end).FitWindow = userData.params.fitWindow; %mm
        Important_coeffs(end).InputDescription = fitsdesc(j);
        
        Important_coeffs(end+1).Variable = categorical(cellstr('p'));
        Important_coeffs(end).Value = cell2mat(jConicFit.p');
        Important_coeffs(end).FitName = categorical(cellstr('Conic'));
        Important_coeffs(end).FitWindow = userData.params.fitWindow;
        Important_coeffs(end).InputDescription = fitsdesc(j);
        
    end
    
    
end

%Get the SPHERICAL Fit?

fitTable = struct2table(S);
fitTable.Properties.VariableUnits = {'','','','mm','mm','mm','',''};
fitTable.Properties.VariableDescriptions = ...
    {'Name of the fit', 'Description of the points used to create the fit',...
    'Function used to create the fit', 'Size of the window used for fit in mm',...
    'X and Y points used in the fit', 'X and Y points generated from fit',...
    'Fit object that can be used to recreate the fit','Fit coefficients'};


coeffTable = struct2table(Important_coeffs);
coeffTable.Properties.VariableUnits = {'','Base Unit: mm','','mm',''};
coeffTable.Properties.VariableDescriptions = {'Name of the variable',...
    'Coefficient Value [anterior, posterior]','Name of the fit',...
    'Window used to analyze the Fit','Description of Input that created the fit'};


%..........................................................................
%Next, Biometry using the corrected points
%..........................................................................

xyPoints = cosFit_corr.segPoints;
[biometryTable, fourierTable] = calculateLegacy_biometry(xyPoints);

%Calculate: CSA, Perimeter, Area, Thickness
imgRes = userData.xyRes(1); img = userData.selectedVariant{1};
cosFit_orig = userData.cosFit_orig;
segPoints_mm = cosFit_orig.segPoints;
shiftAmount = cosFit_orig.ImgShift;

[propsOut,lensMask] = estimateArea(img, segPoints_mm, imgRes, shiftAmount);

biometryTable.CSA =  propsOut.Area; %CSA %mm^2
biometryTable.Perimeter = propsOut.Perimeter; %mm
biometryTable.LensThickness = propsOut.Thickness; %mm ADDED
biometryTable.LensDiameter = propsOut.Diameter; %mm
%biometryTable.LensMask = {lensMask}; %logical image TOO MUCH SPACE

%..........................................................................
%Generate Parameter table
%..........................................................................
params = userData.params;
paramTable = table;
paramTable.Param_ImgROI = {find(userData.imgRoi)};

if isempty(userData.zoomRegion) %use lensMask to get zoom region
    %guess zoom region
    props = regionprops(lensMask,'BoundingBox'); tol = 100; %pxls
    bBox = props.BoundingBox;
    xlim = [bBox(1)-tol, bBox(1)+bBox(3)+tol];
    ylim = [bBox(2)-tol, bBox(2)+bBox(4)+tol];
    
    userData.zoomRegion = {xlim, ylim};
    
    paramTable.Param_ZoomRegion = userData.zoomRegion;
else
    paramTable.Param_ZoomRegion = userData.zoomRegion;
end


paramTable.Param_VariantFunction = ...
    cellstr(func2str(userData.ImgVariants(userData.variantIdx).fcn_cdata));

 
paramTable.Param_InitialRotation = userData.ImgRotated(1).nRotations;
paramTable.Param_ImgRes = params.xyRes(1:2);
paramTable.Param_nRecursions = params.nRecursions;
paramTable.Param_FitDivFactor = params.fitDivFactor;
paramTable.Param_StartFitTol = params.startFitTol;
paramTable.Param_FitWindow = params.fitWindow;
paramTable.Param_FitType = func2str(params.fitType);

paramTable.Properties.VariableUnits =  {'','pixels','','radians','mm',...
    '','mm','mm','mm',''};
paramTable.Properties.VariableDescriptions =  {'Image mask that delineates the ROI',...
    'X and Y Limits of zoom area in image',...
    'Function that manipulates original image to create the image used in analysis',...
    'Number of CCW rotations need to put the lens posterior down',...
    'Image resolution in pixels per mm', 'Number of recursions used to remove outliers',...
    'Divisor used in recursive fits', 'Initial tolerance of recursive fits',...
    'Size of window used for quad/conic fit',...
    'Type of fit used in recursive fits'};

%..........................................................................
%Start Creating the output
%..........................................................................

basicTable = table;
basicTable.SubjectID = cellstr('');
basicTable.ExperimentID = cellstr('');

basicTable.TotalShift = cosFit_orig.shiftAmount; 
basicTable.ImageShift = cosFit_orig.ImgShift;
basicTable.ImageTilt = cosFit_corr.tilt;
basicTable.ImageSize = size(userData.imgRoi);

basicTable.Fits = {fitTable};
basicTable.Coefficients = {coeffTable};


basicTable.Properties.VariableUnits = {'','','mm','mm','radians','pixels','',''};
basicTable.Properties.VariableDescriptions = {'','',...
    'Amount segmentation points were shifted in total',...
    'Initial shift guess (useful in displaying images)',...
    'Tilt value calculated by Cosine Fit',...
    'Size of the image in pixels','Table containing all fits',...
    'Table containing important coefficients'};

logTable = table;
logTable.CreatedBy = cellstr('');
logTable.CreatedDate = datetime('now','Format','M/d/yyyy');
% logTable.Notes = cellstr(''); DELETED

logTable.Properties.VariableDescriptions = {'First and last name of analyzer',...
    'Automatically generated timestamp in the format M/D/Y'};%'Notes relevant to analysis'

analysisOut = cat(2,basicTable, paramTable, biometryTable, fourierTable,...
    logTable);

function [biometryData,fourierData] = calculateLegacy_biometry(xyPointsIn)
%This function is entirely based on legacy code with the exception of the
%table output

% Parameters


%Switch X and Y again
yIn = xyPointsIn(:,1); xIn = xyPointsIn(:,2);

[T_ALL, R_ALL] = cart2pol(xIn, yIn);
T_ALL_AD = T_ALL;
R_ALL_AD = R_ALL;

% Fit to Fourier model  (10 coefficents)
% b = zeros(11,1);
a0 = [2.6; 0.2; -0.9; 0.01; 0.38; -0.03; -0.17; 0.03; 0.07; -0.01; -0.02];

lb = [1, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf];
ub = [5, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf];
data2 = T_ALL_AD;



fitDisp = 'off'; %'iter';


options = optimset('Algorithm', 'Trust-Region-Reflective', 'DerivativeCheck',...
    'on', 'Display', fitDisp); % 'LargeScale', 'off' removed 718
[b, resnorm, residual,exitflag,output, lambda, jacobian] = lsqcurvefit(@fourier_model2,...
    a0, data2, R_ALL_AD, lb, ub, options);

% Fit to Fourier model (20 coefficents)
% a = zeros(21,1);
a0 = [2.6; 0.2; -0.9; 0.01; 0.38; -0.03; -0.17; 0.03; 0.07; -0.01; -0.02; ...
    0.02; 0.01; 0.01; 0.01; 0.01; 0.01; 0.01; 0.01; 0.01; 0.01];

lb = [1, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, ...
    -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf];
ub = [5, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf,...
    Inf, Inf, Inf, Inf, Inf, Inf, Inf];
data2 = T_ALL_AD;

options = optimset('Algorithm', 'Trust-Region-Reflective', 'DerivativeCheck',...
    'on', 'Display', fitDisp);
[a, resnorm, residual,exitflag,output, lambda, jacobian] = lsqcurvefit(@fourier_model, ...
    a0, data2, R_ALL_AD, lb, ub, options);
% get lens contour from new coefficients
t =  -pi:2*pi/10000:pi;

a_rsme = a(1) + a(2).*cos(1.*t) + a(3).*cos(2.*t)+ a(4).*cos(3.*t)...
    + a(5).*cos(4.*t)+ a(6).*cos(5.*t)+ a(7).*cos(6.*t)+ ...
    a(8).*cos(7.*t)+ a(9).*cos(8.*t)+ a(10).*cos(9.*t)+ ...
    a(11).*cos(10.*t) + a(12).*cos(11.*t) + a(13).*cos(12.*t)+ a(14)...
    .*cos(13.*t)+ a(15).*cos(14.*t)+ a(16).*cos(15.*t)+ a(17).*cos(16.*t)+ ...
    a(18).*cos(17.*t)+ a(19).*cos(18.*t)+ a(20).*cos(19.*t)+ ...
    a(21).*cos(20.*t);

b_rsme = b(1) + b(2).*cos(1.*t) + b(3).*cos(2.*t)+ b(4).*cos(3.*t)...
    + b(5).*cos(4.*t)+ b(6).*cos(5.*t)+ b(7).*cos(6.*t)+ ...
    b(8).*cos(7.*t)+ b(9).*cos(8.*t)+ b(10).*cos(9.*t)+ ...
    b(11).*cos(10.*t);
ro = b_rsme;

% compute rmse 20
a_rsme_comp = sum(a(1) + a(2).*cos(1.*T_ALL) + a(3).*cos(2.*T_ALL)+ a(4).*cos(3.*T_ALL)...
    + a(5).*cos(4.*T_ALL)+ a(6).*cos(5.*T_ALL)+ a(7).*cos(6.*T_ALL)+ ...
    a(8).*cos(7.*T_ALL)+ a(9).*cos(8.*T_ALL)+ a(10).*cos(9.*T_ALL)+ ...
    a(11).*cos(10.*T_ALL)+ a(12).*cos(11.*T_ALL) + a(13).*cos(12.*T_ALL)+ ...
    a(14).*cos(13.*T_ALL)+ a(15).*cos(14.*T_ALL)+ a(16).*cos(15.*T_ALL)+ ...
    a(17).*cos(16.*T_ALL)+ a(18).*cos(17.*T_ALL)+ a(19).*cos(18.*T_ALL)+ ...
    a(20).*cos(19.*T_ALL)+ a(21).*cos(20.*T_ALL), 2);

% compute rmse 10
b_rsme_comp = sum(b(1) + b(2).*cos(1.*T_ALL) + b(3).*cos(2.*T_ALL)+ b(4).*cos(3.*T_ALL)...
    + b(5).*cos(4.*T_ALL)+ b(6).*cos(5.*T_ALL)+ b(7).*cos(6.*T_ALL)+ ...
    b(8).*cos(7.*T_ALL)+ b(9).*cos(8.*T_ALL)+ b(10).*cos(9.*T_ALL)+ ...
    b(11).*cos(10.*T_ALL),2);

rmse_10 = sqrt(sum((R_ALL-b_rsme_comp).^2)/length(R_ALL));
rmse_20 = sqrt(sum((R_ALL-a_rsme_comp).^2)/length(R_ALL));

rmse_10_corrected = sqrt(sum((R_ALL-b_rsme_comp).^2)/(length(R_ALL)-11));
rmse_20_corrected = sqrt(sum((R_ALL-a_rsme_comp).^2)/(length(R_ALL)-21));


%Choose structure or table layout
fourierData = table; %struct([]);

fourierData.Fourier10Coeff = force1D(b,1);
fourierData.Fourier20Coeff = force1D(a,1);
fourierData.Fourier10RMSE = [rmse_10, rmse_10_corrected];
fourierData.Fourier20RMSE = [rmse_20, rmse_20_corrected];

fourierData.Properties.VariableDescriptions =...
    {'10 Fourier fit coefficients with DC offset a0 (11 values)',...
    '20 Fourier fit coefficients with DC offset a0 (21 values)',...
    'Fourier 10 fit Root Mean Square Error: [RMSE10, RMSE10 - DOF]',...
    'Fourier 20 fit Root Mean Square Error: [RMSE10, RMSE10 - DOF]'};   

%%
%BIOMETRY HALF
%This function is entirely based on legacy code with the exception of the
%table output

% Parameters
theta = -pi:pi/2:pi;

%Remember that X and Y are switched
X_A = xIn; Y_A = yIn;

% half lens
% Find half lens for corrected original data
Whole_Lens = [X_A, Y_A];
[Y_A_I, Y_A_J] = ind2sub(size(Whole_Lens), find(Whole_Lens(:,2)>=0));

Half_Lens = [];

for ind_X = 1:length(Y_A_I)
    Half_Lens(ind_X, 1) = Whole_Lens(Y_A_I(ind_X),1);
    Half_Lens(ind_X, 2) = Whole_Lens(Y_A_I(ind_X),2);
end

sort_half = sortrows(Half_Lens, 1);
x_dist= abs(diff(sort_half(:,1)));

% Calculate volume based on equation v = h*pi*r^2 where h is the thickness
% of each slice... first slice is eliminated
r_squared = sort_half(:,2).^2;
h=[0; x_dist];
volume_data=pi.*sum(r_squared.*h);

%
clear Half_Lens r_squared
% Compute lens dimensions
dimen20 = a(1) + a(2).*cos(1.*theta) + a(3).*cos(2.*theta)+ a(4).*cos(3.*theta)...
    + a(5).*cos(4.*theta)+ a(6).*cos(5.*theta)+ a(7).*cos(6.*theta)+ ...
    a(8).*cos(7.*theta)+ a(9).*cos(8.*theta)+ a(10).*cos(9.*theta)+ ...
    a(11).*cos(10.*theta)+ a(12).*cos(11.*theta) + a(13).*cos(12.*theta)+ ...
    a(14).*cos(13.*theta)+ a(15).*cos(14.*theta)+ a(16).*cos(15.*theta)+ ...
    a(17).*cos(16.*theta)+ a(18).*cos(17.*theta)+ a(19).*cos(18.*theta)+ ...
    a(20).*cos(19.*theta)+ a(21).*cos(20.*theta);

dimen10 = b(1) + b(2).*cos(1.*theta) + b(3).*cos(2.*theta)+ b(4).*cos(3.*theta)...
    + b(5).*cos(4.*theta)+ b(6).*cos(5.*theta)+ b(7).*cos(6.*theta)+ ...
    b(8).*cos(7.*theta)+ b(9).*cos(8.*theta)+ b(10).*cos(9.*theta)+ ...
    b(11).*cos(10.*theta);

% VOLUME CALCULATIONS FROM COSINE FIT DATA
%10 COEFFICEINT

thetapos = [0:0.001:pi]';
r_pos_10 = b(1) + b(2)*cos(1*thetapos)+ b(3)*cos(2*thetapos)...
    + b(4)*cos(3*thetapos)+ b(5)*cos(4*thetapos)...
    + b(6)*cos(5*thetapos)+ b(7)*cos(6*thetapos)...
    + b(8)*cos(7*thetapos)+ b(9)*cos(8*thetapos)...
    + b(10)*cos(9*thetapos)+ b(11)*cos(10*thetapos);

[ZposCos_10, YposCos_10] = pol2cart(thetapos, r_pos_10);

total_10 = [ZposCos_10 YposCos_10];
sort_10 = sortrows(total_10, 1);
z_dist_10 = [0; abs(diff(sort_10(:,1)))];
y_10_squared = sort_10(:,2).^2;

Vol_10=pi*sum(y_10_squared.*z_dist_10);

%20 COEFFICENT
r_pos_20 = a(1) + a(2).*cos(1.*thetapos) + a(3).*cos(2.*thetapos)+ a(4).*cos(3.*thetapos)...
    + a(5).*cos(4.*thetapos)+ a(6).*cos(5.*thetapos)+ a(7).*cos(6.*thetapos)+ ...
    a(8).*cos(7.*thetapos)+ a(9).*cos(8.*thetapos)+ a(10).*cos(9.*thetapos)+ ...
    a(11).*cos(10.*thetapos)+ a(12).*cos(11.*thetapos) + a(13).*cos(12.*thetapos)+ ...
    a(14).*cos(13.*thetapos)+ a(15).*cos(14.*thetapos)+ a(16).*cos(15.*thetapos)+ ...
    a(17).*cos(16.*thetapos)+ a(18).*cos(17.*thetapos)+ a(19).*cos(18.*thetapos)+ ...
    a(20).*cos(19.*thetapos)+ a(21).*cos(20.*thetapos);

[ZposCos_20, YposCos_20] = pol2cart(thetapos, r_pos_20);

total_20 = [ZposCos_20 YposCos_20];
sort_20 = sortrows(total_20, 1);
z_dist_20 = [0; abs(diff(sort_20(:,1)))];
y_20_squared = sort_20(:,2).^2;

Vol_20=pi*sum(y_20_squared.*z_dist_20);
% SURFACE AREA FROM COSINE FIT
%10 coefficent
s_area10 = quadl(@fourier_model_SA, 0, pi, 1.e-6, 0, b)*2*pi;

%20 coefficent
s_area20 = quadl(@fourier_model_SA20, 0, pi, 1.e-6, 0, a)*2*pi;
% Biometry based on 10 coeffecient cosine fit
csa_10 = sum(b.^2.*[pi; pi/2; pi/2; pi/2; pi/2; pi/2; pi/2; pi/2; pi/2; ...
    pi/2; pi/2]);
csa_20 = sum(a.^2.*[pi; pi/2; pi/2; pi/2; pi/2; pi/2; pi/2; pi/2; pi/2; ...
    pi/2; pi/2; pi/2; pi/2; pi/2; pi/2; pi/2; pi/2; pi/2; pi/2;...
    pi/2; pi/2]);


%PUT IT IN A TABLE
biometryData = table();

%10 Coefficient Fit
biometryData.Fourier10_AnteriorThickness = dimen10(1); %bA_10; %mm 
biometryData.Fourier10_PosteriorThickness =  dimen10(3); %bP_10; %mm
biometryData.Fouier10_LensThickness = dimen10(1)+dimen10(3); %T_10 %mm (Ant + Pos Thickness)
biometryData.Fourier10_LensDiameter = dimen10(2)+dimen10(4); %D_10; %mm 
biometryData.Fourier10_CSA = csa_10; %mm^2
biometryData.Fouier10_SurfaceArea = s_area10; %mm^2
biometryData.Fourier10_Volume = Vol_10; %mm^3

%20 CoefficientFit
biometryData.Fourier20_AnteriorThickness = dimen20(1); %bA_20; %mm 
biometryData.Fourier20_PosteriorThickness = dimen20(3); %bP_20; %mm
biometryData.Fouier20_LensThickness = dimen20(1)+dimen20(3); %T_20 %mm (Ant + Pos Thickness)
biometryData.Fourier20_LensDiameter = dimen20(2)+dimen20(4); %D_20; %mm
biometryData.Fourier20_CSA = csa_20; %mm^2
biometryData.Fouier20_SurfaceArea = s_area20; %mm^2
biometryData.Fourier20_Volume = Vol_20; %mm^3


%Pixel Count on Raw Image
biometryData.CSA =  NaN; %CSA %mm^2 --> calculated later
biometryData.Perimeter = NaN; %mm  --> calculated later


%Volume from Centered and Aligned Data
biometryData.Volume = volume_data; %mm^3

%Add Metadata
biometryData.Properties.VariableUnits =...
    {'mm','mm','mm','mm','mm^2','mm^2','mm^3',...Fourier10
    'mm','mm','mm','mm','mm^2','mm^2','mm^3',...Fourier20
    'mm^2','mm','mm^3'}; %Pixel Count & Volume

biometryData.Properties.VariableDescriptions =...
    {'Fourier 10 Fit: Anterior Thickness',...
    'Fourier 10 Fit: Posterior Thickness',...
    'Fourier 10 Fit: Lens Thickness',...
    'Fourier 10 Fit: Lens Diameter',...
    'Fourier 10 Fit: Cross Sectional Area (CSA)',...
    'Fourier 10 Fit: Surface Area',...
    'Fourier 10 Fit: Volume',...
    'Fourier 20 Fit: Anterior Thickness',...
    'Fourier 20 Fit: Posterior Thickness',...
    'Fourier 20 Fit: Lens Thickness',...
    'Fourier 20 Fit: Lens Diameter',...
    'Fourier 20 Fit: Cross Sectional Area (CSA)',...
    'Fourier 20 Fit: Surface Area',...
    'Fourier 20 Fit: Volume',...
    'Cross Sectional Area (CSA) from Pixel Count',...
    'Perimter from Pixel Count',...
    'Volume from Pixel Count'};


function [propsOut,lensMask] = estimateArea(img, segPoints_mm, imgRes, shiftAmount)
%This function was created separately because I  no longer use the function
%that generated the image and area
%INPUT:
    %img - mxn array the same size of the image
    %segPoints_mm - mx2 array contianing the original segmentation points
    %imgRes - scalar containing the image resolution (pixels/mm)
    %shiftAmount - [xc_guess, yc_guess], posterior down
%OUTPUT:
    %propsOut - structure containing the following fields in mm:
        %Perimeter
        %Area
        %Lens Thickness
        %Diameter


%Sort the Segmentation Points
[theta, rho] = cart2pol(segPoints_mm(:,1), segPoints_mm(:,2));
sortedThetaRho = sortrows([theta, rho],1);
[xSorted,ySorted] = pol2cart(sortedThetaRho(:,1), sortedThetaRho(:,2));

%convert it back to pixels from mm
x = (xSorted + shiftAmount(1))*imgRes;
y = (ySorted + shiftAmount(2))*imgRes;

%Create a lens mask 
lensMask = roipoly(img, x, y);

%Get the region properties of the mask
propsAll = regionprops(lensMask, {'Area', 'Perimeter',...
    'MinorAxisLength','MajorAxisLength','Orientation'});
areaAll = propsAll.Area;

%use only the larges
propsLens = propsAll(areaAll == max(areaAll));


propsOut.Area = propsLens.Area/(imgRes*imgRes);
propsOut.Perimeter = propsLens.Perimeter/imgRes;
propsOut.Thickness = propsLens.MinorAxisLength/imgRes;
propsOut.Diameter = propsLens.MajorAxisLength/imgRes;
propsOut.Orientation = propsLens.Orientation;

%..........................................................................

%--------------------------------------------------------------------------
%NEW FUNCTIONS: SHADOWGRAPH
%--------------------------------------------------------------------------

function [intMask, gMask, iThresh, gThreshOut] = extractEdgeMask(img, imgType)
%input filtered image


%shift the image to be between

%%
%determine what an adequate intensity threshold is
[count, ~] = imhist(img); %count(1:50) = 0; %get rid of the first few values

isRidiculous = count> mean(count)+2*std(count);

if any(isRidiculous)
    count(isRidiculous) = max(count(~isRidiculous));
end

gmag = imgradient(img); gmag = uint8(mat2gray(gmag)*255);
gThresh = multithresh(gmag);

%find the mode of the image
[~,peakIdx] = max(count); imgMode = peakIdx - 1;%starts from zero

%data before the mode
beforePeak = count(1:peakIdx);
%first minimum before the peak
minIdx = find(beforePeak < 5,1,'last');

if isempty(minIdx)
    minIdx = find(beforePeak < min(beforePeak)+5,1,'last');
end

noiseMin = minIdx - 1;

%assume the noise is normally distributed
halfWidth = (peakIdx - minIdx);
noiseMax = imgMode+halfWidth;

%assume there is overlap between noise and data
iThresh = noiseMax - round(halfWidth/1.85);

switch imgType
    case 1
        filterSize = [20 10];
    case 2
        filterSize = [5 20];
    case 3
        filterSize = [5 20];
    case 4 %not verified
        filterSize = [5 20];
    otherwise
        filterSize = [5 20];
end

if ~isempty(halfWidth)
    gThreshOut = gThresh-halfWidth;
else
    gThreshOut = gThresh;
end
intMask = medfilt2(img > iThresh, filterSize);
gMask = medfilt2(gmag > gThreshOut);
% intMaskEdge = bwmorph(intMask, 'remove');

isPlot = false;
if isPlot
    figure(20); imagesc(img), colormap gray
    figure(21); imagesc(img>iThresh), colormap gray
    figure(22); imagesc(intMask), colormap gray
    figure(23); imagesc(gMask); colormap gray
    figure(24); imagesc(intMask & gMask), colormap gray
    
    figure(25); imagesc(intMaskEdge), colormap gray
    figure(26); imagesc(intMaskEdge & gMask), colormap gray
end

function userData = straightenLens(userData, singlePerimeterPoints)
%This function untilts and centers the lens using a cosine fit

imgRes = userData.xyRes(1);

if nargin < 2
    singlePerimeterPoints = userData.edgePoints;
end


%Switch the orientation so that it is NOT posterior up but posterior (messy) RIGHT
X_raw = singlePerimeterPoints(:,1);
Y_raw = singlePerimeterPoints(:,2);

offset_X = min(X_raw);
offset_y = min(Y_raw);

%find an approximate center for the lens

X_mm = X_raw./imgRes;
Y_mm = Y_raw./imgRes;


%estimate xc and yc
minX_val = min(X_mm); maxX_val = max(X_mm);
minX_ind = X_mm == minX_val; %logical indexing
maxX_ind = X_mm == maxX_val; %logical indexing

%find the x values corresponding to the min&max Y
minY_val = mean(Y_mm(minX_ind)); maxY_val = mean(Y_mm(maxX_ind));

xc_guess = mean([minX_val, maxX_val]); %+ offset_X;
yc_guess = mean([minY_val, maxY_val]); %+ offset_Y;

Xmm_orig = X_mm - xc_guess;
Ymm_orig = Y_mm - yc_guess;



%%
%tilt and decenter

%flip X and Y to make the equitorial axis vertical (this affects the cosine fit)
X_in = Ymm_orig; Y_in = Xmm_orig;


%first, make a guess of the orientation
[propsOut,lensMask] = estimateArea(userData.selectedVariant{1},...
    [X_in,Y_in], imgRes, [yc_guess, xc_guess]);


orient = propsOut.Orientation;

if orient > 0
    tilt_guess = 90 - orient; %positive theta: quadrant1
else
    tilt_guess = -(90 + orient); %negative theta: quadrant2
end

tilt_guess_radians = tilt_guess *(pi/180); %this is to be rotated CCW!!! Units --> rad
%convert to polar coordinates
[T_in, R_in] = cart2pol(X_in, Y_in);


%initialize fit
nCoeff = 20;
initialValues= zeros(nCoeff+3,1); %coeff + xc, yc, theta0
bound_params = [0.1, 0.1,10*(pi/180)]; %abs(xc, yc, theta0 [radians])

iniCoeff = inf(1,nCoeff-1);

lowerBounds = [2, -iniCoeff, -bound_params];
upperBounds = [5,  iniCoeff,  bound_params];

% iniVals = zeros(size(initialValues)); %coefficient input
data = [force1D(X_in); force1D(Y_in)];

if userData.runSilent
    fitDisp = 'off';
else
    fitDisp = 'iter';
end

options = optimset('Algorithm', 'Trust-Region-Reflective',...
    'DerivativeCheck','on', 'Display',fitDisp);
%%

%apply fit
[a, resnorm, residual, exitflag, output, lambda, jacobian] = lsqcurvefit(@cosine_tiltDecentration,...
    initialValues, data, R_in, lowerBounds, upperBounds, options);


%if my guess is too much different than the fit, use the guess SOOO RISKY!
if abs(tilt_guess_radians - a(end,1)) > 5 * (pi/180) %5 degrees
%     a(end,1) = tilt_guess_radians;
end

%apply the fit
fit_r = feval(@cosine_tiltDecentration, a, data); 

%sort it by theta
fit_polar_sorted = sortrows([T_in, fit_r]);

%convert fit to cartesian coordinates
[fit_x, fit_y] = pol2cart(fit_polar_sorted(:,1)', fit_polar_sorted(:,2)');

%%

%remove tilt using an affine tranformation matrix
Transform_cw = @(th, shift) [cos(th), sin(th), shift(1);...
                              -sin(th),  cos(th), shift(2);...
                              0,        0,       1      ]; %th --> rad          

%shift X and y untilt coorndinates based on fit result
xc = a(end-2,1); yc = a(end-1,1); theta_tilt = a(end,1);



T_mat = Transform_cw(theta_tilt, [xc,yc]); %create matrix
X_trans = [X_in,Y_in,ones(size(X_in))]'; %rearrange points for matrix input
Xp = T_mat * X_trans; %calculate new points X_prime

X_corr = Xp(1,:); Y_corr = Xp(2,:); %X and Y corrected

%..........................................................................
%correct the fit also
%..........................................................................
Xfit_trans = [fit_x; fit_y; ones(size(fit_x))]; %[x1,y2,z2;x2,y2,z2;...]
Xfit_p = T_mat * Xfit_trans;
Xfit_corr = Xfit_p(1,:); Yfit_corr = Xfit_p(2,:);
%..........................................................................


%%
%Make sure the center goes through the equator

%estimate xc and yc
minY_val2 = min(Yfit_corr); maxY_val2 = max(Yfit_corr);
minY_ind2 = Yfit_corr == minY_val2; %logical indexing
maxY_ind2 = Yfit_corr == maxY_val2; %logical indexing

%find the x values corresponding to the min&max Y
minX_val2 = mean(Yfit_corr(minY_ind2)); maxX_val2 = mean(Yfit_corr(maxY_ind2));

xc_final = mean([minX_val2, maxX_val2]); %+ offset_X;
yc_final = mean([minY_val2, maxY_val2]); %+ offset_Y;

X_out = X_corr - xc_final;
Y_out = Y_corr - yc_final;

Xfit_out = Xfit_corr - xc_final;
Yfit_out = Yfit_corr - yc_final;


%yc_guess in xc_total becuase it was calculated before coordinates were
%switched -xc because it is applied in the transform
xc_total = yc_guess - xc + xc_final; 
yc_total = xc_guess - yc + yc_final;

%%
%rotate corresponding image
% xyRes = userData.xyRes; imgRes = xyRes(1);
% 
imgVariant = userData.selectedVariant{1};
[ySize, xSize] = size(imgVariant(:,:,1));
xSample_mm = (1:xSize)/imgRes; ySample_mm = (1:ySize)/imgRes;
% 
% imgRoi = userData.imgRoi;
% if ~all(imgRoi(:))
%     roi_props = regionprops(imfill(imgRoi,'holes'),'BoundingBox','Area');
%     if numel(roi_props) > 1
%         [~, biggestROI_idx] = max([roi_props.Area]);
%         roi_props = roi_props(biggestROI_idx);
%     end
%     
%     bBox = roi_props.BoundingBox;
%     tol = 1; %mm
%     bBox_mm = bBox / imgRes;
%     bBox_mm(1) = bBox_mm(1) - xc_guess - tol;
%     bBox_mm(2) = bBox_mm(2) - yc_guess - tol;
%     bBox_mm(3:4) = bBox_mm(3:4) + tol;
%     
%     topLeftX = bBox_mm(1); topLeftY = bBox_mm(2);
%     bottomLeftY = topLeftY + bBox_mm(4);
% %     bBox_rect = [topLeftX, bottomLeftY, bBox_mm(3:4)];
% else
%     bBox_mm = [];
% end
% 
% 
% 
% 
% hFig = figure('Visible','on','Units','normalized','Position',[0 0 1 1]); 
% hImg = imagesc(xSample_mm-xc_guess, ySample_mm-yc_guess, imgVariant); 
% set(gca,'XTick',[],'YTick',[]);
% colormap gray; 
% hold on; hP(1) = plot(Y_in, X_in,'or');
% hP(2) = plot(fit_y, fit_x,'LineWidth',2);
% hold off;
% 
% if isempty(bBox_mm)
%     F = getframe(gca);
% else
%     F = getframe(gca);
% end
% delete(hFig);

%%
%OUTPUT
segPoints = [Y_out; X_out]';
xyFit = [Yfit_out; Xfit_out]';
shiftAmount = [yc_total, xc_total];

%unshifted
cosFit_orig.xScale = xSample_mm-xc_guess;
cosFit_orig.yScale = ySample_mm-yc_guess;
cosFit_orig.ImgShift = [xc_guess, yc_guess];
cosFit_orig.segPoints = [Y_in, X_in];
cosFit_orig.xyFit = [fit_y', fit_x'];
cosFit_orig.shiftAmount = [xc_guess, yc_guess];
cosFit_orig.coeff = a;

%corrected
cosFit_corr.segPoints = segPoints;
cosFit_corr.xyFit = xyFit;
cosFit_corr.shiftAmount = shiftAmount;
cosFit_corr.tilt = theta_tilt; %radians, clockwise
% cosFit_corr.fitSnapshot = F.cdata;
cosFit_corr.coeff = [];
cosFit_orig.fitObj = [];



userData.cosFit_orig = cosFit_orig;
userData.cosFit_corr = cosFit_corr;

function userData = updateWindow(userData)
%This function updates the points within the window and works on a an
%entire lens (not ant and pos surface)
%INPUT: the cosine fit must be ran before running this function
%OUTPUT: adds "isInWindow" field to costFit_orig & cosFit_corr

windowSize = userData.params.fitWindow;
dWindow = windowSize/2;

segPts_orig = userData.cosFit_orig.segPoints;
xyFit_orig = userData.cosFit_orig.xyFit;
segPts_corr = userData.cosFit_corr.segPoints;
xyFit_corr = userData.cosFit_corr.xyFit;

isInWindow_orig = {abs(segPts_orig(:,1)) < dWindow, abs(xyFit_orig(:,1)) < dWindow};
isInWindow_corr = {abs(segPts_corr(:,1)) < dWindow, abs(xyFit_corr(:,1)) < dWindow};

userData.cosFit_orig.isInWindow = isInWindow_orig; %{pts, fit}
userData.cosFit_corr.isInWindow = isInWindow_corr; %{pts, fit}

function output = splitPts(cosFit, isUseWindow)
%This function returns the anterior and posterior points within the window
%INPUT: 
    % cosFit_corr or cosFit_orig structure 
    % isUseWindow - logical scalar that determines if you want to return
                % only points within the valid window
%OUTPUT: {antPts, posPts; antFit, posFit}

if nargin < 2
    isUseWindow = true;
end

%corrected points
isInWindow = cosFit.isInWindow;

if isUseWindow
    segPoints = cosFit.segPoints(isInWindow{1},:);
    xyFit = cosFit.xyFit(isInWindow{2},:);
else
    segPoints = cosFit.segPoints;
    xyFit = cosFit.xyFit;
end

isAntPts = segPoints(:,2) <= 0;
isAntFit = xyFit(:,2) <= 0;

antPts = segPoints(isAntPts,:);
posPts = segPoints(~isAntPts,:);

antFit = xyFit(isAntFit,:);
posFit = xyFit(~isAntFit,:);

output = {antPts, antFit; posPts, posFit};

%--------------------------------------------------------------------------

%..........................................................................

function [aScanInt, aScanDerivative] = calculateThreshold(img, mode,...
    intMultFactor, derMultFactor, h)
%mode determines positive, negative or both
%aScanInt = {aScanImg, threshold}
%aScanDer = {aScanProcessed, derivativeThresh}


img = double(img);
avgInt = mean(img(:));
stdevI = std(img(:));
intThresh = max([avgInt + stdevI*intMultFactor,0]); %fail safe to make sure the min is 0


%Calculate Intensity
avgI = avgInt;
normalizedImg = img-avgI;
imgDer = getDerivative(normalizedImg, h);


switch mode %added 3/25/2013
    case 0 %positive derivative
        imgDer(imgDer < 0) = 0;
        imgD = imgDer;
    case 1 %negative derivative
        imgDer(imgDer > 0) = 0;
        imgDer = abs(imgDer);
        imgD = imgDer;
        
    case 2 %use both edges to create the threshold
        imgD = imgDer;
        imgD = abs(imgD);
end

%for calculations, remove all zero values
imgD(imgD==0)=[];
avgDiff = mean(imgD(:));
stdevDiff = std(imgD(:));
derivThresh = avgDiff + stdevDiff*derMultFactor;

aScanInt = {img, intThresh};
aScanDerivative = {imgDer, derivThresh};

function matrixOut = getDerivative(matrixIn, h)

A = matrixIn; B = circshift(matrixIn,[-h, 0]);
matrixOut = B-A; matrixOut(end-h+1:end,:) = [];
matrixOut = cat(1, zeros(h,size(matrixIn,2)), matrixOut);

function [peakInd, peakVals] = findMaximums(gradientCol,n, yRes, isFirstPeak, separationInMM)
%gradientCol = difference of image in the y dirrection for every h pixels
%for a particular column
%n = number of peaks to detect


minPixelSeparation = round(yRes*separationInMM);

nVals= numel(gradientCol);
%Make sure its a column array
gradientCol = force1D(gradientCol);

%Concatenate with column containing index number
gradientArray = cat(2,(1:nVals)', gradientCol);

%Sort it based on value
if isFirstPeak
    filteredArray = cat(2, flipud(gradientArray), [1;abs(diff(gradientArray(:,1)))]); %WARNING, DIDN'T SORT BY INTENSITY FIRST COME FIRST SERVE 5/29/2013
else
    sortedArray = sortrows(gradientArray, 2);
    sortedArray = cat(2, sortedArray, [1;abs(diff(sortedArray(:,1)))]);
    filteredArray = sortedArray;
end

%remove pixels with no intensity
filteredArray (filteredArray (:,2)==0,:) = [];

isTooClose = filteredArray(:,3) < minPixelSeparation;

%remove pixels separated by less than then minPixelSeparation
while any(isTooClose) && numel(isTooClose) > 1
    firstOne = find(isTooClose, 1,'last'); %counting from the bottom up
    stopInd = find(~isTooClose,1,'last');
    
    if isempty(stopInd) || stopInd > firstOne
        stopInd =  1;
    end
    
    filteredArray(stopInd:firstOne-1,:) = [];
    filteredArray(:,3) = [100; abs(diff(filteredArray(:,1)))];
    
    if stopInd ~= 1
        isTooClose = filteredArray(:,3) < minPixelSeparation;
    else
        isTooClose = false;
    end
end


if size(filteredArray,1) < n
    %     peakInd = [];
    %     peakVals = [];
    if ~isempty(filteredArray)
        peakInd = num2cell(filteredArray(:,1)); %peakInd{n} = []; peakInd = peakInd';
        peakVals = num2cell(filteredArray(:,2)); %peakVals{n} = []; peakVals = peakVals';
    else
        peakInd = [];
        peakVals = [];
    end
else
    peaksFound = sortrows(filteredArray(end-n+1:end,:),1);
    peakInd = num2cell(peaksFound(:,1)); %peaksFound(:,1);
    peakVals = num2cell(peaksFound(:,2)); %peaksFound(:,2);
end

if length(peakInd) < n
    peakInd{n,1} = [];
    peakVals{n,1} = [];
end

%..........................................................................


function [idxMaskAbove, idxMaskBelow,...
    idxAbove, idxBelow] = segmentMatrix(boundary, matrix)
%This function returns the indices of the matrix on either side of a given
%line/boundary (above and below).  The boundary is interpolated linearly for values of the
%boundary not explicitly given
%
%..........................................................................
%inputs:
%..........................................................................
%BOUNDARY: a 2D array containing the x and y values of the boundary
%           surface
%note - any deisred tolerance range must already be included
%MATRIX:  the matrix in which the bondary is drawn (for example an image)
%
%
%..........................................................................
%outputs
%..........................................................................
%idxMaskAbove - logical indexing mask of area of matrix above boundary
%idxMaskBelow - logical indexing maks of area of matrix below boundary
%
%
%idxAbove: a 1D matrix containing the indices in the matrix input about the boundary
%idxBelow: a 1D matrix containing the indices in the matrix input below the boundary
%
%==========================================================================

xIn = boundary(:,1); yIn = boundary(:,2);
[nRows, nCols] = size(matrix);

%make sure border spans from beginning to end of matrix by spreading the
%first and last values of the boundary to the edge of the matrix

if min(xIn) ~=1
    xIn = [1; xIn];
    yIn = [yIn(1); yIn];
end

if max(xIn) ~= nCols
    xIn(end+1) = nCols;
    yIn(end+1) = yIn(end);
end

[~, ia, ~] = unique(xIn, 'stable');
xUnique = xIn(ia); yUnique = yIn(ia);

%interpolate for misssing values of y
xFit = 1:nCols;
yFit = interp1(xUnique,yUnique,xFit,'nearest');

%make sure any y values out of range is brought to matrix boundary
yFit(yFit > nRows) = nRows;
yFit(yFit < 1) = 1;

% xInd = repmat(1:nCols, [nRows,1]);
yInd = repmat((1:nRows)', [1, nCols]);
% boundXInd = repmat(force1D(xFit,1), [nRows,1]);
boundYInd = repmat(force1D(yFit,1), [nRows,1]);

idxMaskAbove = yInd < boundYInd;
idxMaskBelow = yInd > boundYInd;

matIndices = reshape(1:numel(matrix), [nRows,nCols]);
idxAbove = force1D(matIndices(idxMaskAbove));
idxBelow = force1D(matIndices(idxMaskBelow));

%..........................................................................
% %This code takes longer and will give me actual indices
%..........................................................................

%%make sure x and y values don't come out of the matrix
% isInvalid = (xFit > nCols | xFit < 1) | (yFit > nRows | yFit < 1);
% xFit(isInvalid) = []; yFit(isInvalid) = [];
%
% %Convert to a cell array
% iSurfaceFit = num2cell([force1D(xFit), force1D(round(yFit))],2);
%
% %for each column, get the indices
% idxAbove = cellfun(@(idx) sub2ind(imgSize, 1:idx(2),repmat(idx(1),[1,idx(2)])),...
%     iSurfaceFit, 'UniformOutput', false);
% idxAbove = [idxAbove{:}];
%
%
% idxBelow = cellfun(@(idx) sub2ind(imgSize, idx(2):imgRows, repmat(idx(1),[1,idx(2)])),...
%     iSurfaceFit, 'UniformOutput', false);
% idxBelow = [idxBelow{:}];

drawnow; %useless line of code that keeps commented line of code in function


function [totalFit, centralFit, dataPoints,R,p] = fitDataPoints_mm(pointsIn, params)


xInlier = pointsIn(:,1); 
yInlier = pointsIn(:,2);

%in this case, we dont have additional possible candidates (yet...we can find them from the edges?)
allX = xInlier; 
allY = yInlier; 


if numel(xInlier) < 5
    totalFit = [];
    centralFit = [];
    dataPoints = [];
    R = [];
    p = [];
    return
end

xFit = min(xInlier):0.05:max(xInlier); %0.05mm increments

[yFit, p, pmm, inlierPoints] = runRecursion_mm(allX, allY, xInlier,yInlier,xFit, params);
R = 1./(2*pmm); %quick validity check
% disp(R)

totalFit = [force1D(xFit),force1D(yFit)];
centralFit = [force1D(xFitOut), force1D(yFitOut)];
dataPoints = inlierPoints;

function [inlierPoints, finalFitY, p_mm] = runRecursion_mm(allX, allY, rawX, rawY, xFit, params)
%This function fits points within a specified window recursively using a
%particular fit
%INPUT:
    %pointsIn - input the points within the window of interest
    %params - xyRes, nRecursions, startFitTol, fitDivTol
%NOTE: Units in MM!!

%INPUT:
    %rawX, raw Y --> known inliers; 
    %allX, allY --> possible inliers
    %xFit --> xPoints you'd like to fit with final curve
    %params: nRecursions, startFitTol, fitDivTol,  fitType

n = params.nRecursions;
fType = params.fitType;


startWidth = params.startFitTol; %0.1; %mm
divFactor = params.fitDivFactor; %1.5;

dw_mm = startWidth; %changed so the width gets smaller on every recursion

%PRELIMINARY ANALYSIS, use all points but use inliers as a guidance
p_prelim = feval(fType, rawX, rawY);
yFit_prelim = polyval(p_prelim,allX);
upperBound_prelim = yFit_prelim + dw_mm/2;
lowerBound_prelim = yFit_prelim - dw_mm/2;
prelimOutlier = bitor(allY > upperBound_prelim, allY < lowerBound_prelim);
xData = allX; yData = allY;
yData(prelimOutlier) = [];
xData(prelimOutlier) = [];

% xData = rawX; yData = rawY;

if n < 1
    n = 1;
end

isAnt = p_prelim(1) > 0;

isUseConic = false; %true;

if isAnt
    sNo = 3;
else
    sNo = 4;
end

for i = 1:n  
     
    if isUseConic
        [p_mm, iGOF] = newConicFit(xData, yData, sNo, true);
        yFit = feval(p_mm, xData);
    else
        p_mm = feval(fType, xData, yData);
        yFit = polyval(p_mm,xData);
    end
    
    upperBound = yFit + dw_mm/2;
    lowerBound = yFit - dw_mm/2;
    isOutlier = bitor(yData > upperBound, yData < lowerBound);
    yData(isOutlier) = [];
    xData(isOutlier) = [];
    
    if isUseConic
        finalFitY = feval(p_mm, xFit);
    else
        finalFitY = polyval(p_mm, xFit);
    end
    
    %PLOT RECURSIONS
    % upperBound(isOutlier) = []; lowerBound(isOutlier) = [];
    % figure(1225); plot(allX, allY, 'xr'), axis ij;
    % hold on; plot(xData, upperBound,'--m'); plot(xData, lowerBound,'--m');
    % plot(xData, yData, 'xg');  hold off;
    % pause;
     
    dw_mm = (startWidth/(divFactor*i));
end

if isUseConic
   yFinal = feval(p_mm, allX);
else
    yFinal = polyval(p_mm, allX);

end
upperBound_final = yFinal + dw_mm;
lowerBound_final = yFinal - dw_mm;
isFinalOutlier = bitor(allY > upperBound_final, allY < lowerBound_final);
finalX = force1D(allX(~isFinalOutlier)); finalY = force1D(allY(~isFinalOutlier));
inlierPoints  = [finalX, finalY];%[force1D(xData),force1D(yData)];


%..........................................................................

function [bestGroup, allGroups] = selectBestGroup(rawX, rawY, fType, isConvex, userData)
%This function selects the group of points that best fits the function
%..........................................................................
%INPUTS:
%..........................................................................
%rawX = 1D array containing x points
%rawY = 1D array containing y points
%fType = function handle to fitting function that outputs coefficients
%isConvex = should the polynomial be positive or negatively curved (valid
%           for quadratics
%userData = structure of user inputs containing y resolution (userData.xyRes(2))
%..........................................................................
%OUTPUTS
%..........................................................................
%bestGroup = structure containg statistics (x and y points, outliers etc)
%              of best fit
%allGroups = structure containting all groups created from original data


u_yRes = userData.xyRes(2);
dw_pxls = 0.15*u_yRes; %0.3mm

bestGroup = [];
allGroups = [];


%have a mincurvature
if userData.scanWidth > 10
    minA = 0.0005;
elseif userData.scanWidth > 7
    minA = 0.0005;
elseif userData.scanWidth < 7
    minA = 0.0001;
else
    minA = 0;
end


%break data points into groups
originalGroups = classifyData(rawX,rawY, userData);
groups = originalGroups;

if isempty(groups)
    groups(1).data = [force1D(rawX),force1D(rawY)];
    groups(1).n = length(rawX);
    originalGroups(1).n = groups(1).n;
end

nGroups  = length(groups);
allGData = cat(1,groups(:).data);



if nGroups == 1
    cond2Thresh = 1.3;
else
    cond2Thresh = 1.7;
end

cond2 = numel(rawX)/sum([groups.n])> cond2Thresh ;
nALines = userData.scanWidth*userData.xyRes(1);



if nGroups > 1 && userData.scanWidth < 10 && userData.imgType == 1
    %make a group out of all the valid groups
    
    groups(end+1).data = allGData;
    groups(end).n = length(allGData(:,1));
elseif  cond2 || (nALines < 600 && userData.scanWidth > 8)
    %make a group out of all raw points
    groups(end+1).data = [rawX, rawY];
    groups(end).n = length(rawX);
end


groupSizeThresh = max([originalGroups.n]);

if cond2
    groupSizeThresh = max([groupSizeThresh, round(length(rawX)/3)]);
end

allGData = cat(1,groups(:).data);
allX = allGData(:,1); allY = allGData(:,2);


nGroups = length(groups);
isGroupLegit = false(1,numel(groups));
%determine which group gives best fit
for i = 1:nGroups
    data = groups(i).data;
    
    if userData.imgType == 4
        %is exvivo lens?
        dataIncluded = data;
        uniqueData = uniqueByCol(dataIncluded,2); %changed this to 2
    elseif userData.imgType ~= 1 %is not cornea
        %use only the central data
        nPoints = size(data,1); cntrIdx = round(size(data,1)/2);
        halfPtsIncl = round((0.95*nPoints)/2);  %central 70 changed from .7
        dataIncluded = data(max([1,cntrIdx-halfPtsIncl]): min([size(data,1),cntrIdx+halfPtsIncl]),:);
        uniqueData = dataIncluded;
        
    else
        %is cornea
        dataIncluded = data;
        uniqueData = uniqueByCol(dataIncluded,1); %changed this to 2
    end
    
    groupX = uniqueData(:,1); groupY = uniqueData(:,2);
    
    p = ransac(uniqueData, dw_pxls, fType, 60, 8, isConvex, userData); %dataIn, dw_pxls, nIter, nPoints
    
    if all(~p) %if p has an error with RANSAC
        p = feval(fType,groupX,groupY);
    end
    
    yFit = polyval(p,allX);
    
    %Find the outliers of all points
    upperBound = yFit + dw_pxls/2;
    lowerBound = yFit - dw_pxls/2;
    isOutlier = bitor(allY > upperBound, allY < lowerBound);
    nOutliers = sum(isOutlier);
    inlierX = allX(~isOutlier); inlierY = allY(~isOutlier);
    nInliers = numel(inlierX);
    inlierFit = polyval(p,inlierX);
    sseInlier = mean((inlierFit - inlierY).^2); %sum((inlierFit - inlierY).^2);
    
    %Find the number of outliers within the group
    groupFit = polyval(p, groupX);
    groupUpperBound = groupFit + dw_pxls/2;
    groupLowerBound = groupFit - dw_pxls/2;
    isGroupOutlier = bitor(groupY > groupUpperBound, groupY < groupLowerBound);
    nGroupOutliers = sum(isGroupOutlier);
    groupSize(i) = groups(i).n;
    groupsseInlier = mean((groupFit - groupY).^2);
    
    %Save findings to a structure
    iGroup.global.inliers = [force1D(inlierX), force1D(inlierY)];
    iGroup.global.outliers = [force1D(allX(isOutlier)),force1D(allY(isOutlier))];
    iGroup.global.nOutliers = nOutliers;
    iGroup.global.nInliers = nInliers;
    iGroup.global.bounds = [force1D(allX), force1D(upperBound), force1D(lowerBound)];
    iGroup.global.sse = sseInlier;
    iGroup.global.fit = [allX, yFit];
    
    iGroup.local.inliers = [force1D(groupX(~isGroupOutlier)),...
        force1D(groupY(~isGroupOutlier))];
    iGroup.local.outliers = [force1D(groupX(isGroupOutlier)),...
        force1D(groupY(isGroupOutlier))];
    iGroup.local.nOutliers = nGroupOutliers;
    iGroup.local.bounds = [force1D(groupX),...
        force1D(groupUpperBound), force1D(groupLowerBound)];
    iGroup.local.nPoints = groupSize(i);
    iGroup.local.sse = groupsseInlier;
    iGroup.local.fit = [groupX, groupFit];
    
    iGroup.p = p;
    iGroup.yRange = [min(groupY), max(groupY)];
    iGroup.idx = i;
    
    iGroup.p = p;
    iGroup.yRange = [min(groupY), max(groupY)];
    iGroup.idx = i;
    
    if nInliers > 0 && ~any(isnan(yFit))
        isGroupLegit(i) = true;
    end
    
    groupAnalysis(i) = iGroup;
    
    
    
end %end for loop


%COST FUNCTION
%use: how it fits overall data (A) AND how it fits groups (B)
validGroupsIdx = 1:numel(groupAnalysis);
sseInlier = arrayfun(@(x) groupAnalysis(x).global.sse,validGroupsIdx);
groupsseInlier= arrayfun(@(x) groupAnalysis(x).local.sse,validGroupsIdx);
nOutliers = arrayfun(@(x) groupAnalysis(x).global.nOutliers,validGroupsIdx);
nGroupOutliers = arrayfun(@(x) groupAnalysis(x).local.nOutliers,validGroupsIdx);

overallCost = sseInlier/max(sseInlier) + (nOutliers/numel(allX));
groupCost = (groupsseInlier/max(groupsseInlier))./groupSize + (nGroupOutliers./groupSize);

if max(groupCost) < 1 %modified 2/27/2014
    costFcn = 0.15*overallCost + exp(3*groupCost);
else
    costFcn = 0.35*overallCost + 3*groupCost;
end

costFcn_cell = num2cell(costFcn);

[groupAnalysis(:).costFcn] = deal(costFcn_cell{:});

%remove all polynomials with the wrong curve
%add a constraint for minimum curvature 2/27/2014

polys = cat(1,groupAnalysis.p);
As = polys(:,1);

fcnStack = dbstack;
if userData.imgType ~= 5 || strcmpi(fcnStack(3).name, 'fit2surfaces')
    if isConvex
        isAValid = As>0 & abs(As) > minA;
    else
        isAValid = As< 0 & abs(As) > minA;
    end
else
    isAValid = true(size(As));
end

validAs = As(isAValid);

%if the scan width is greater than 4mm, use a minimium group size as a
%constraint

if userData.scanWidth > 5
    
    if userData.scanWidth > 12 && userData.imgType == 1
        groupSizeThresh = max([groupSizeThresh - 50,0]);
        
    elseif userData.scanWidth > 12 && userData.imgType == 2
        if isConvex
            groupSizeThresh = 0;
        else
            %do nothing
        end
    end
    isAboveMinSize = groupSize/min([groupSizeThresh, 250]) > 0.5;
    
else
    isAboveMinSize = groupSize/groupSizeThresh > 0; %all true
end



if ~any(validAs)
    allGroups  = groupAnalysis;
    return
else
    validGroups = groupAnalysis(isAValid & isAboveMinSize' & isGroupLegit');
    
end

if isempty(validGroups)
    validGroups = groupAnalysis(isAValid);
end

validCostFcn = [validGroups.costFcn];

%if the first group has comparably the lowest cost function use it
if (abs(min(validCostFcn)-validCostFcn(1)) < 0.3) && userData.imgType == 1
    gIdx = validGroups(1).idx;
else
    [minCost,validGIdx] = min(validCostFcn); %try using both outlier and sse
    gIdx = validGroups(validGIdx).idx; %get the overall group index
    
    if gIdx > numel(originalGroups) && nALines > 400 %if you're going to use the biggest group, make sure there is no close second
        [secondMinCost, otherIdx] = min(validCostFcn(1:end-1));
        
        if abs(minCost-secondMinCost) > std(validCostFcn);
            gIdx = otherIdx;
        end
    end
end


[groupAnalysis(:).isBest] = deal(false);
groupAnalysis(gIdx).isBest = true;


bestGroup = groupAnalysis(gIdx);

% if gIdx ~=3
%     groupAnalysis(end) = [];
% end

allGroups = groupAnalysis;

function finalGroups = classifyData(xData,yData, userData)

u_yRes = userData.xyRes(2);
isPlotAll = false;

makeColumnVector = @(array) reshape(array, [numel(array),1]);

if userData.scanWidth > 12 && userData.imgType == 1
    maxDy_mm = 0.3;
else
    %any pixel with greater than 0.3mm difference is considered an outlier
    maxDy_mm = 0.15; %mm
end

maxDy = maxDy_mm*u_yRes;

%make sure data is column
xData = makeColumnVector(xData);
yData = makeColumnVector(yData);

finalGroups = groupByHistogram(xData, yData, maxDy, true, userData);

if isempty(finalGroups)
    return
end

%if there is only one group of data, try the second method
% groups2 = groupByDerivative(xData, yData, maxDy);

%isPlot?
if isPlotAll
    plotGroups(finalGroups)
end

function finalGroups = groupByHistogram(xData, yData, maxDy, isAnalyze, userData)
isPlotAll = false;

isPlot = isPlotAll && isAnalyze;
bins = (min(yData):max(yData))';count =histc(yData,bins);
groupPercentage = 0.1; %min number of pixels to be considered a group
minGroupMembers = round(groupPercentage*numel(yData));

if isAnalyze
    artifactThresh = 0.05*max(xData); %min number of pixels to make an artifact
else
    artifactThresh = 0.5*max(xData); %min number of pixels to make an artifact
end

%if there aren't any groups increase the threshold a little
isGTMaxDy = false; divisor = 4;

while ~any(isGTMaxDy)
    histData = [bins,count]; avgCount = mean(count(count>0));
    dyThresh = max([avgCount/divisor,1]);
    hasCount = count>dyThresh; validValues = histData(hasCount,:);
    
    %Intermediate variables used to find the values that are greater than maxDy
    idx = (2:size(validValues,1))';
    isGTMaxDy = diff(validValues(:,1)) > maxDy;
    isGTIdx= idx(isGTMaxDy);
    
    if divisor ==1
        break
    else
        divisor = divisor -1;
    end
    
end

if isempty(isGTIdx) %even data distribution
    finalGroups = [];
    return
end

dataIdx = sort([isGTIdx-1; isGTIdx]);

%convert to y location index
dataLocation = validValues(dataIdx,1);

while any(diff(dataLocation)==0)
    idx = find(diff(dataLocation)==0)+1;
    dataLocation(idx) = dataLocation(idx)+1;
end

if dataLocation(1) == bins(1)
    isBlankStart = true;
else
    isBlankStart = false;
    dataLocation = [bins(1);dataLocation];
end

allRanges = [dataLocation, [dataLocation(2:end)-1; bins(end)]];

isEven = mod(1:size(allRanges,1),2)==0;

if isBlankStart
    isDataRange = isEven;
else
    isDataRange = ~isEven;
end

zeroRange = allRanges(~isDataRange,:);
dataRange = allRanges(isDataRange,:);

nGroups = size(dataRange,1);

if isempty(nGroups)
    finalGroups = [];
    return
end

%Find weighted means for every data range
%Group Data and find means
origData = [xData,yData]; origDataSorted = sortrows(origData,2);
ySorted = origDataSorted(:,2);
allGroups = struct([]);

for i = 1:nGroups
    iStart = dataRange(i,1); iEnd = dataRange(i,2);
    
    %convert to histogram idx
    [~,idx1] = min(abs(iStart-bins)); [~, idx2]= min(abs(iEnd-bins));
    
    iBins = bins(idx1:idx2); iCount = count(idx1:idx2);
    isArtifact = iCount > artifactThresh;
    %Recreate the original data set from the histogram
    %     zeroCount = iCount == 0;
    %     iCount2 = iCount; iCount2(zeroCount) = [];
    %     iBins2 = iBins; iBins2(zeroCount) = [];
    %     iA = cell2mat(cellfun(@(x) ones(x(1),1)+x(2),...
    %         num2cell(cat(2,iCount2, (0:numel(iCount2)-1)'),2), 'UniformOutput', false));
    %     recreatedData = iBins2(iA);
    %
    %     %Calculate weighted average
    %     iAvg = mean(recreatedData); %round(sum((iBins.*iCount))/sum(iCount));
    %     avgVals(i) = iAvg;
    %     stdVals(i) = std(recreatedData);
    
    
    if isPlot
        if i ==1
            figure;
            hold on;
        end
        plot(iBins,iCount,'Color', rand(1,3));
    end
    
    
    %Original Data %added 3/25/2013
    %find all data within the range
    idxStart = find(abs(ySorted-iStart)==min(abs(ySorted-iStart)),1,'first');
    idxEnd = find(abs(ySorted-iEnd)==min(abs(ySorted-iEnd)),1,'last');
    
    iData = origDataSorted(idxStart:idxEnd,:);
    iData = sortrows(iData,1); %sort by x values
    
    
    %PUT LINEAR ARTIFACT IN IT'S OWN SUBGROUP if artifact is on either end
    %of the group
    
    if isAnalyze && any(isArtifact)
        iY = iData(:,2);
        %remove artifact
        artifactIdx = iBins(isArtifact);
        artifactRange = unique(cell2mat(arrayfun(@(x) (x-2:x+2)', artifactIdx,'UniformOutput',false)));
        
        nArtifacts = numel(artifactRange);
        
        %make it a row vector
        artifactRange = reshape(artifactRange, [1,nArtifacts]);
        iYMatrix = repmat(iY, [1,nArtifacts]);
        artifactMatrix = repmat(artifactRange,[numel(iY),1]);
        
        %check for any artifacts in the bin
        isNotData = logical(sum(iYMatrix == artifactMatrix,2));
        
        iSubData = iData(isNotData,:);
        iSubX = iSubData(:,1); iSubY = iSubData(:,2);
        
        
        %if artifact is not in the middle of the group, separate it
        minPxDiff = 5;
        if all(abs(iSubX - mean(iData(:,1))) > minPxDiff)
            
            %         if any( bitor(abs(iSubX - iData(1,1)) < minPxDiff,...
            %                 abs(iSubX - iData(end,1)) < minPxDiff)) && all(abs(iSubX - mean(iData(:,1))) > minPxDiff)
            
            gIdx = numel(allGroups)+1;
            %Put all the artifacts in it's own group
            allGroups(gIdx).n = size(iSubData,1);
            allGroups(gIdx).avg = [mean(iSubX), mean(iSubY)];
            allGroups(gIdx).data = iSubData;
            allGroups(gIdx).stdev = std(iSubY);
            allGroups(gIdx).hasArtifact = any(isArtifact);
            allGroups(gIdx).artifactIdx = iBins(isArtifact);
            
            %connectivity of the layer %added 4/2/2013
            allGroups(gIdx).connectivity = std(abs(diff(iSubY)));
            
            %remove these horizontal artifacts from the group
            iData(isNotData,:) = [];
            
        else
            %continue
            
        end
        
        isArtifact = false;
        
    end
    
    
    iX = iData(:,1); iY = iData(:,2);
    gIdx = numel(allGroups)+1;
    allGroups(gIdx).n = size(iData,1);
    allGroups(gIdx).avg = [mean(iX),mean(iY)];
    allGroups(gIdx).data = iData;
    allGroups(gIdx).stdev = std(iY);
    allGroups(gIdx).hasArtifact = any(isArtifact);
    allGroups(gIdx).artifactIdx = iBins(isArtifact);
    %get the largest artifact idx %added 4/26/2013
    [~,id] = max(iCount(isArtifact));
    
    if ~isempty(id)
        id = id(1);
    end
    
    artifactBins = iBins(isArtifact);
    allGroups(gIdx).largestArtifactIdx = artifactBins(id);
    
    %connectivity of the layer %added 4/2/2013
    allGroups(gIdx).connectivity = std(abs(diff(iY)));
    
    
    %For plot legend
    stdVals(i) = std(diff(iY)); %std(iY);
    legendStr{i} = sprintf('sigma = %5.2f',stdVals(i));
    
    
    
end

if isPlot
    legend(legendStr);
    title('Histogram Distribution')
end


%If the distance between groups is less than 1 then run it again after
%removing points that are noise


%threshold is more than 5
% if dyThresh < 5
%     whoCalledMe = dbstack;
%
%     allGroups.connectivity
%     %and we are not in a nested function
%     if numel(whoCalledMe) >  2 && ~strcmpi(whoCalledMe(2).name)
%         groupByHistogram(xData, yData, maxDy, false, userData)
%
%     end
% end

finalGroups = allGroups;

if isAnalyze %only analyze if necessary
    
    switch userData.imgType
        case 1 %cornea (in-vivo)
            if userData.scanWidth > 12 && userData.imgType ==1
                minStd = 5;
                maxConn = 800; %should be 80
            elseif userData.scanWidth < 6
                minStd = 5;
                maxConn = 800; %should be 80
            elseif userData.scanWidth < 11
                minStd = 12;
                maxConn = 800; %should be 80
            else
                %remove anything with a stdev less than 6 or greater than 100
                minStd = 15; %25; modified: 5/12/2014
                maxConn = 800; %should be 200
            end
        case 2%lens (in-vivo)
            
            if userData.scanWidth > 12
                minStd = 0;
                maxConn = 40;
            elseif userData.scanWidth < 6
                minStd = 3;
                maxConn = 80; %should be 80
            elseif userData.scanWidth < 10
                minStd = 3;
                maxConn = 80; %should be 80
            else
                minStd = 4;
                maxConn = 80; %should be 80
            end
            
        otherwise
            %settings of lens for now (modified: 3/6/2014)
            if userData.scanWidth < 6
                minStd = 3;
                maxConn = 80; %should be 80
            else
                minStd = 10;
                maxConn = 80; %should be 80
            end
            
    end
    
    isArtifact = [finalGroups.stdev]  < minStd;
    
    finalGroups(isArtifact) = [];
    
    %also remove any groups that are KNOWN artifacts
    hasArtifact = [finalGroups.hasArtifact];
    finalGroups(hasArtifact) = [];
    
    gCount = [finalGroups(:).n];
    isEnoughMembers = gCount > minGroupMembers;
    if any(isEnoughMembers)
        finalGroups(~isEnoughMembers) = [];
    end
    
    
    
    isRandom = [finalGroups.connectivity]>maxConn;
    if any(~isRandom)
        finalGroups(isRandom) = [];
    end
    
end

function pOut = ransac(dataIn, dw_pxls, fType, maxIter, nPoints, isConvex, userData)
%This function runs n iterations on the x and y data using random x number
%of points returning the polynomial with the best fit
%
%Inputs: data - 2D arrray containing x and y data in columnar fashion [x1,y1;x2,y2],
%dw_pxls = number of pixels for tolerance range; scalar
%nIter = max number of iterations; scalar (optional) [default-60]
%nPoints = number of points randomly selected; scalar (optional) [default-8]
%
%Outputs: p - quadratic polynomial coefficients; returns false if there is
%an error

%Use only unique values in the second column
% xIn = dataIn(:,1); yIn = dataIn(:,2);

%Parse inputs
switch nargin
    case [0 1] %not enough inputs
        pOut = false;
        return
    case 2
        fType = @quadFit;
        maxIter = 60; %default max
        nPoints = 8; %default points
        isConvex = true;
    case 3
        maxIter = 60; %default max
        nPoints = 8; %default points
        isConvex = true;
    case 4
        nPoints = 8; %default points
        isConvex = true;
    case 5
        isConvex = true;
end

switch userData.imgType
    case 1 %in vivo cornea
        Rmin = 6.5; Rmax = 8.3; %mm
    case 2 %in vivo lens
        
        if isConvex
            Rmin = 6; Rmax = 20; %mm
        else
            Rmin = 6; Rmax = 50; %mm
        end
    otherwise
        Rmin = 5; Rmax = 30;
end

%cost threshold
costThresh = 20;

%take only unique y values to minimize the effect of horizontal artifacts
uniqueData = uniqueByCol(dataIn, 2);

%Find out if the percent discarded is too high
if size(uniqueData,1)/size(dataIn,1) < 0.5
    uniqueData = dataIn;
end

xIn = uniqueData(:,1); yIn = uniqueData(:,2);


dataSize = size(uniqueData,1);
xData = uniqueData(:,1); yData = uniqueData(:,2);

%if data size is smaller than number of points, redetermine number of
%points to be 1% of the data
if dataSize < nPoints
    nPoints = round(0.01*dataSize);
end

p={};
costFcn = deal([]);
isGo = true;
count = 0;

xRes = userData.xyRes(1); yRes = userData.xyRes(2);
Rmm = [];
while isGo
    %increment the count
    count = count+1;
    
    %if you haven't found anything, for the posterior lens, use only the
    %central datapoints
    if (count > maxIter/2) && userData.imgType == 2 && ~isConvex
        tol = 1*xRes;
        cntrIdx =  floor(dataSize/2);
        %use points closer to the center
        iIdx = randi(round([max([1,cntrIdx-tol]),...
            min([cntrIdx+tol,dataSize])]),[1,nPoints]);
    else
        
        %randomly select data points unique in the Y direction
        iIdx = randperm(dataSize);
        iIdx = iIdx(1:nPoints);
        
    end
    iX = xData(iIdx); iY = yData(iIdx);
    iXmm = iX/xRes; iYmm = iY/yRes;
    
    %fit the points using least squares
    iP = feval(fType,iX,iY);
    iPmm = feval(fType, iXmm, iYmm);
    iRmm = 1/(2*iPmm(1));
    iFit = polyval(iP,xIn);
    
    
    %create a tolerance range to classify outliers
    upperBound = iFit + dw_pxls/4;
    lowerBound = iFit - dw_pxls/4;
    isOutlier = bitor(yIn > upperBound, yIn < lowerBound);
    nOutliers = sum(logical(nonzeros(isOutlier)));
    
    %determine how well inliers fit to the curve
    inlierX = xIn(~isOutlier); inlierY = yIn(~isOutlier);
    inlierFit = polyval(iP, inlierX);
    sse = mean((inlierFit - inlierY).^2);
    %OTHER OPTIONS FOR SSE:
    %sum((inlierFit - inlierY).^2/max((inlierFit - inlierY).^2));
    %sum((inlierFit - inlierY).^2);
    %std((inlierFit - inlierY).^2);
    %mean((inlierFit - inlierY).^2);
    
    
    
    
    %if quadratic is in the right range, keep it
    isInRadiusRange = (abs(iRmm) < Rmax && abs(iRmm) > Rmin);
    isRightInflection = (iP(1) > 0 && isConvex) || (iP(1) < 0 && ~isConvex);
    if isRightInflection && isInRadiusRange
        
        p{end+1} = iP;
        
        %Cost Fucntion: root(error) + percent of outliers (makes outliers
        %dominating term)
        costFcn(end+1) = (sse)^0.25 + (nOutliers/dataSize)*100;
        Rmm(end+1) = iRmm;
    else %if quadratic is inflected (curved downward), skip it
        
        
        if count > maxIter
            isGo = false;
        end
        
        
        continue
    end
    
    if costFcn(end) < costThresh && count > maxIter %changed || to && (Modified: 1/30/2015)
        isGo = false;
    end
    
    
end %end while

if ~isempty(p)
    p = cat(1,p{:});
else
    pOut = false;
    return
end


[err,idxOut] = min(costFcn);
% disp(count)
% fprintf(1,'\nError is: %10.5f\n\n', err);
pOut = p(idxOut,:);

function output = uniqueByCol(data, col)
%This function returns that data with unique values in a particular column,
%keeps row with first occurence of the number
% input: data - columnar data, col - column number

if col > size(data,2) || col < 1
    output = false;
    return
end

%Get the column of interest
colData = data(:,col);
[~,uniqueIdx,~] = unique(colData,'stable');
output = data(uniqueIdx,:);


function grayMap(hAx)

colormap(hAx, repmat(linspace(0,1, 64)',[1,3]));


function userData = updateImgVariant(userData, variantIdx)
%This function updates the userData structure for the image variant.

if nargin < 2
    variantIdx = userData.variantIdx;
else 
    userData.variantIdx = variantIdx;
end

imgRotated = userData.ImgRotated(1).cdata;
ImgVariants = userData.ImgVariants;
variantLabel = ImgVariants.label;
fcn_variant = ImgVariants(variantIdx).fcn_cdata;
variantImg = fcn_variant(imgRotated);

userData.selectedVariant = {variantImg, variantLabel};


% --- Executes on button press in check_apply2all.
function check_apply2all_Callback(hObject, eventdata, handles)
% hObject    handle to check_apply2all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_apply2all

if get(hObject, 'Value')
    set(handles.check_overwrite, 'Enable', 'on');
else
    set(handles.check_overwrite, 'Enable', 'off');
end


% --- Executes on button press in check_overwrite.
function check_overwrite_Callback(hObject, eventdata, handles)
% hObject    handle to check_overwrite (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_overwrite


% --- Executes during object creation, after setting all properties.
function popup_variantSelect_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_variantSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function popup_RoiDisplay_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_RoiDisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
