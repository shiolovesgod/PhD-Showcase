function varargout = uiloadmultidir(varargin)
% UILOADMULTIDIR MATLAB code for uiloadmultidir.fig
%      UILOADMULTIDIR, by itself, creates a new UILOADMULTIDIR or raises the existing
%      singleton*.
%
%      H = UILOADMULTIDIR returns the handle to a new UILOADMULTIDIR or the handle to
%      the existing singleton*.
%
%      UILOADMULTIDIR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in UILOADMULTIDIR.M with the given input arguments.
%
%      UILOADMULTIDIR('Property','Value',...) creates a new UILOADMULTIDIR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before uiloadmultidir_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to uiloadmultidir_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help uiloadmultidir

% Last Modified by GUIDE v2.5 03-Feb-2017 13:49:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @uiloadmultidir_OpeningFcn, ...
                   'gui_OutputFcn',  @uiloadmultidir_OutputFcn, ...
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


% --- Executes just before uiloadmultidir is made visible.
function uiloadmultidir_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to uiloadmultidir (see VARARGIN)

% Choose default command line output for uiloadmultidir
handles.output = hObject;

% set(hObject, 'WindowStyle','modal');

%..........................................................................
%populate listbox with the current directory
%..........................................................................

prevDir = getappdata(hObject, 'prevDir');
if isempty(prevDir)
    dirName = cd;
else
    dirName = prevDir;
end

dirInfo = dir(dirName);

isDir = cat(1,dirInfo.isdir);
isDir(1:2) = false; %'.' and '..'

 set(handles.edit_imgDir, 'String', dirName);
 
if any(isDir)
    dirNames = sortrows({dirInfo(isDir).name}');
    hListbox = handles.listbox_imgFiles;
    nDirs = numel(dirNames);
    set(hListbox, 'String',dirNames,'Max', nDirs);
    
    setappdata(handles.b_loadDir, 'rootDir', dirName); 
    set(hListbox, 'Value',1);
    set(handles.b_ok, 'Enable', 'on');   
end

% Update handles structure
guidata(hObject, handles);

uiwait(hObject);



% UIWAIT makes uiloadmultidir wait for user response (see UIRESUME)
% uiwait(handles.f_uiloadimgdir);


% --- Outputs from this function are returned to the command line.
function varargout = uiloadmultidir_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


isOk = getappdata(hObject, 'isOk');

if isOk
    rootDir = getappdata(handles.b_loadDir, 'rootDir');
    hListbox = handles.listbox_imgFiles;
    fNames = get(hListbox, 'String');
    
    fullFiles = arrayfun(@(i) fullfile(rootDir, fNames{i}), 1:numel(fNames),...
        'UniformOutput', false);    
    fullFilesOut = fullFiles(get(hListbox,'Value'));
    
    isLoadSubjectInfo = get(handles.check_subjectInfo,'Value');
else
    fullFilesOut = [];
    isLoadSubjectInfo = [];
end


varargout{1} = fullFilesOut ;

varargout{2} = isLoadSubjectInfo;

%Delete the figure
delete(hObject);

% --- Executes on selection change in listbox_imgFiles.
function listbox_imgFiles_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_imgFiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_imgFiles contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_imgFiles


%==========================================================================
%Button Callback
%==========================================================================

% --- Executes on button press in b_loadDir.
function b_loadDir_Callback(hObject, eventdata, handles)
% hObject    handle to b_loadDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hFig = ancestor(hObject, 'figure');
prevDir = getappdata(hFig, 'prevDir');

rootDir = uigetdir(prevDir);

if ~rootDir
    return
end

setappdata(hFig,'prevDir', rootDir);


set(handles.edit_imgDir, 'String', rootDir);
setappdata(hObject, 'rootDir', rootDir);

populateListbox(handles);

function populateListbox(handles)

rootDir  = getappdata(handles.b_loadDir, 'rootDir');

if isempty(rootDir)
    rootDir = cd;
end

dirInfo = dir(rootDir);
isDir = cat(1,dirInfo.isdir);

dirNames = {dirInfo(isDir).name}';
dirNames(1:2) = []; %get rid of . and ...

%check to see if it has the excel file
if logical(get(handles.check_subjectInfo, 'Value'))
   
    isDirValid = false(numel(dirNames), 1);
    for j = 1:numel(dirNames)
        %keep it if it has the xls or xlsx file
        jD = dir(fullfile(rootDir, dirNames{j}, '*.xls'));
        jD2 = dir(fullfile(rootDir, dirNames{j}, '*.xlsx'));
        
        if ~isempty(jD) || ~isempty(jD2) || true
            isDirValid(j) = true;
            
        end
    end    
    dirNames(~isDirValid) = [];
end


nDirs = numel(dirNames);

hListbox = handles.listbox_imgFiles;
set(hListbox, 'String', sortrows(dirNames), 'Max', max([nDirs, 1]), 'Value',1);

if nDirs > 0
%     set(hListbox, 'Value', 1:nFiles);
    set(handles.b_ok, 'Enable','on');
else
    set(hListbox, 'Value', 0);    
    set(handles.b_ok, 'Enable','off');
end

% --- Executes on button press in b_ok.
function b_ok_Callback(hObject, eventdata, handles)
% hObject    handle to b_ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hFig = ancestor(hObject,'figure');
setappdata(hFig, 'isOk', true);

uiresume(hFig);


function checkDir(dirNames, ext)

%check the directory for the presence of a certain type of file

%==========================================================================
%CREATE FCN
%==========================================================================

% --- Executes during object creation, after setting all properties.
function edit_imgDir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_imgDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function listbox_imgFiles_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_imgFiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes when user attempts to close f_uiloadimgdir.
function f_uiloadimgdir_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to f_uiloadimgdir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


hObject = ancestor(hObject, 'figure');
% buttons = [handles.b_corneaFname, handles.b_lensFname, handles.b_retinaFname];
% 
% arrayfun(@(iHandle) setappdata(iHandle, 'img', []), buttons);
% arrayfun(@(iHandle) setappdata(iHandle, 'fname', ''), buttons);
setappdata(hObject, 'isOk', false);

if isequal(get(hObject, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    uiresume(hObject);
else
    % The GUI is no longer waiting, just close it
    delete(hObject);
end


% --- Executes during object creation, after setting all properties.
function popup_nImgs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_nImgs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in popup_nImgs.
function popup_nImgs_Callback(hObject, eventdata, handles)
% hObject    handle to popup_nImgs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_nImgs contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_nImgs


% --- Executes on button press in check_subjectInfo.
function check_subjectInfo_Callback(hObject, eventdata, handles)
% hObject    handle to check_subjectInfo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_subjectInfo
