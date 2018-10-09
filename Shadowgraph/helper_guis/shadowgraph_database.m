function varargout = shadowgraph_database(varargin)
% SHADOWGRAPH_DATABASE MATLAB code for shadowgraph_database.fig
%      SHADOWGRAPH_DATABASE, by itself, creates a new SHADOWGRAPH_DATABASE or raises the existing
%      singleton*.
%
%      H = SHADOWGRAPH_DATABASE returns the handle to a new SHADOWGRAPH_DATABASE or the handle to
%      the existing singleton*.
%
%      SHADOWGRAPH_DATABASE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SHADOWGRAPH_DATABASE.M with the given input arguments.
%
%      SHADOWGRAPH_DATABASE('Property','Value',...) creates a new SHADOWGRAPH_DATABASE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before shadowgraph_database_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to shadowgraph_database_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help shadowgraph_database

% Last Modified by GUIDE v2.5 27-Mar-2017 10:08:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @shadowgraph_database_OpeningFcn, ...
    'gui_OutputFcn',  @shadowgraph_database_OutputFcn, ...
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


% --- Executes just before shadowgraph_database is made visible.
function shadowgraph_database_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to shadowgraph_database (see VARARGIN)

% Choose default command line output for shadowgraph_database
handles.output = hObject;


%Intialize Filters
hFilters = [handles.popup_filterEye, handles.popup_filterImgType,...
    handles.popup_filterValidated, handles.popup_filterSegmented,...
    handles.popup_filterExcluded, handles.popup_filterNotes];

arrayfun(@(iHandle) popup_filterGeneric_Callback(iHandle, 'manualcall',...
    handles), hFilters);

%Set default tab
hAllTabs = [handles.tab_image, handles.tab_graph, handles.tab_table,...
    handles.tab_settings];
switchTab_Callback(hAllTabs(1), [], handles);

%Initialize table list (i.e. disable them until an image is loaded)
populateTableList(handles); %initialize graph plot settings in this function


%Initialize the subject menu delete
uimenu('Tag','menu_delete','Label','Remove Selected',...
    'Callback', @removeElement,'Parent', handles.cmenu_listbox);


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes shadowgraph_database wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = shadowgraph_database_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%=========================================================================
%MENUS
%=========================================================================

% ------------------------------------------------------------------------
function menu_file_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function menu_import_Callback(hObject, eventdata, handles)
% hObject    handle to menu_import (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function menu_importImage_Callback(hObject, eventdata, handles)
% hObject    handle to menu_importImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function menu_importDirectory_Callback(hObject, eventdata, handles)
% hObject    handle to menu_importDirectory (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


hFig = ancestor(hObject, 'figure');
baseDir = getappdata(hFig, 'baseDir');

[dirFileLocations, isLoadSubjectInfo] = uiloadmultidir(baseDir);

if isempty(dirFileLocations)
    return
end

%Check to see what is already loaded?


%parse the input and create a basis for database (do not keep file
%locations, just keep the relevant information:id, age, species etc.. If
%you don't know the id, use file information as an identifier)

[newTableContainer, outputFlag] = parseFileDirectory(dirFileLocations,...
    isLoadSubjectInfo);

oldTableContainer = getappdata(hFig, 'tableContainer');

%merge new and old containers
if ~isempty(oldTableContainer)
    tableContainer = mergeTableContainers(oldTableContainer,...
        newTableContainer, true); %isForceUnique = true
else
    tableContainer = newTableContainer;
end

tableContainer_clean = removeDatabaseDuplicates(tableContainer);

helpdlg('Load Completed Successfully');

setappdata(hFig, 'baseDir', baseDir);
setappdata(hFig, 'tableContainer', tableContainer_clean);

populateSubjectListbox(handles, tableContainer.Subject);
populateTableList(handles);
updateGraphPanel(handles);



function menu_importSubject_Callback(hObject, eventdata, handles)
% hObject    handle to menu_importSubject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function menu_load_Callback(hObject, eventdata, handles)
% hObject    handle to menu_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% tableContainer = getappdata(ancestor(hObject, 'figure'),'tableContainer');
% displaySagittalTable(tableContainer);
[fName, pName] = uigetfile('*.mat', 'Select the MATLAB file containing the subject data');


if ~fName
    return
end

%load the file into the workspace

load(fullfile(pName, fName));

if ~exist('tableContainer', 'var')
    errordlg('The file selected does not contain the subject information structure.',...
        'Invalid File.');
    return
end

newTableContainer = tableContainer;

if isempty(newTableContainer)
    return
end

%Overwriting for now, but need to merge with existing tableContainer
hFig = ancestor(hObject, 'figure');
oldTableContainer = getappdata(hFig, 'tableContainer');

%try to merge new and old in future: newTableContainer

%clean up the database
newTableContainer_clean = removeDatabaseDuplicates(newTableContainer);

%rewrite the table container
setappdata(hFig, 'tableContainer', newTableContainer_clean);
populateSubjectListbox(handles, newTableContainer_clean.Subject);
populateTableList(handles); %popups for adding fields
initializeTableSettings(handles); %default main table fields

updateGraphPanel(handles);
updateDataTableDisplay(handles); %update the main data table


function menu_save_Callback(hObject, eventdata, handles)
% hObject    handle to menu_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Determine if it's all or selected to be saved
% menuName = get(hObject,'label');

hFig = ancestor(hObject, 'figure');
tableContainer = getappdata(hFig, 'tableContainer');

if isempty(tableContainer)
    errordlg('There is no information to be saved');
    return
end

[fName, pName] = uiputfile('*.mat','Save Subject as:');

if ~pName
    return
end

set(hFig, 'Pointer','watch');
%Save the info to the MATLAB file
saveLoc = fullfile(pName, fName);
save(saveLoc, 'tableContainer','-v7.3');

helpdlg(sprintf('Database successfully saved in: %s', saveLoc),...
    'Save Successful');
set(hFig, 'Pointer','arrow');

function menu_analyze_Callback(hObject, eventdata, handles)
% hObject    handle to menu_analyze (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function menu_analyzeSagittal_Callback(hObject, eventdata, handles)
% hObject    handle to menu_analyzeSagittal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function menu_analyzeSagittalCurrent_Callback(hObject, eventdata, handles)
% hObject    handle to menu_analyzeSagittalCurrent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hFig = ancestor(hObject, 'figure');
hListboxE = handles.listbox_experiment;

%Determine which is experiment selected
tableContainer = getappdata(hFig, 'tableContainer');
selectedVal = get(hListboxE, 'Value');
expTable = getappdata(hListboxE, 'table');


if isempty(expTable) || isempty(selectedVal) || isempty(tableContainer)
    %nothing is exists
    return
end

expTableRow = expTable(selectedVal(1),:); %only analyze the first selection

if expTableRow{1,'ImgType'} == 'Coronal'
    %not ready for coronal images
    return
end

if isfield(tableContainer, 'Sagittal')
    sagittalTable = tableContainer.Sagittal;
else
    sagittalTable = [];
end

%Run the semiAutoSegment_shadowgraph program in verbose mode
thisAnalysisOut = sagittalSegment(tableContainer,...
    expTableRow, sagittalTable, false); %runSilent = false

if isempty(thisAnalysisOut)
    %user cancelled
    return
end

thisAnalysisOut{:,'CreatedBy'} = cellstr(strtrim(get(handles.edit_analyzedBy,'String')));

%append it to the sagittal table
if isempty(sagittalTable)
    Sagittal = thisAnalysisOut;
else
    Sagittal = cat(1,sagittalTable,thisAnalysisOut);
end


%Update the table container
tableContainer.Sagittal = Sagittal;
setappdata(hFig, 'tableContainer', tableContainer);
updateSubjectSummaryText(handles);

updateDataTableDisplay(handles); 
populateTableList(handles); %update popup for adding to data table

listbox_experiment_Callback(handles.listbox_experiment, [], handles);

function menu_analyzeMultiple_Callback(hObject, eventdata, handles)
% hObject    handle to menu_analyzeMultiple (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hFig = ancestor(hObject, 'figure');
tableContainer = getappdata(hFig, 'tableContainer');

if isempty(tableContainer)
    %No data found
    return
end

experimentTable = tableContainer.Experiment;
subjectTable = tableContainer.Subject;

%First determine which data to analyze


imgClass = {'Sagittal'}; %user can decide this later
imgType = 5; %alternative is 6 for Coronal

%Get all the experiments with this image type
qStr = buildQueryString('ImgType', imgClass);
query = {'ImgType', qStr};
[expTable_filtered, logicalIdx] = filterTable(experimentTable, query);

if isempty(expTable_filtered)
    return
end

%Determine if you would like to overwrite
isOverwrite = false; %for now

if isfield(tableContainer, 'Sagittal') && ~isempty(tableContainer.Sagittal) %Sagital analysis
    Sagittal = tableContainer.Sagittal;
else
    Sagittal = [];
    isOverwrite = true;
end

if ~isOverwrite
    %Remove those that already exist
    qStr = buildQueryString('ExperimentID', Sagittal.ExperimentID);
    query = {'ExperimentID', qStr};
    [tableOut, isExpDone] = filterTable(expTable_filtered, query);
    
    %remove it from the list if it's done already
    expTable_filtered(isExpDone,:) = [];
end

if isempty(expTable_filtered)
    %nothing left
    return
end


%Let users select from list dialog which they will use
[selection, ok] = listdlg('ListString',expTable_filtered.ExperimentID,...
    'SelectionMode','multiple', 'ListSize',[240 400],...
    'PromptString','Select files to analyze:','Name','Sagittal Analysis');

if ok && ~isempty(selection)
    expTable_filtered = expTable_filtered(selection, :);
else
    %user didn't select anything
    return
end


nRows = size(expTable_filtered,1);

%Get a waitbar going
hWaitbar = makeWaitbar('Analyze Sagittal Images', 'Batch Analysis');

analysisTable = Sagittal;
analysisFlag = false(nRows, 1);
for i = 1:nRows
    
    expTableRow = expTable_filtered(i,:);
    
    %Update waitbar
    if ishandle(hWaitbar)
        waitbar(i/nRows, hWaitbar,...
            sprintf('Analyzing Sagittal Image:\n%s (%d of %d)',...
            strrep(char(expTableRow.ExperimentID),'_',' '), i, nRows));
    else
        break
    end
    
    %Run the semiAutoSegment_shadowgraph program in silent mode
    thisAnalysisOut = sagittalSegment(tableContainer,...
        expTableRow, analysisTable, true); %runSilent = true;
    
    %If there is a problem
    if isempty(thisAnalysisOut)
        analysisFlag(i) = true;
        continue
    end
    
    if isempty(Sagittal)
        Sagittal = thisAnalysisOut;
    else
        Sagittal = cat(1,Sagittal,thisAnalysisOut);
    end
end

Sagittal{:,'CreatedBy'} = cellstr(strtrim(get(handles.edit_analyzedBy,...
    'String')));

if ~ishandle(hWaitbar)
    %user cancelled
    errordlg('User cancelled analysis.');
else
    %normal analaysis, give user report on failed and successful
    delete(hWaitbar)
    helpdlg('Analysis Complete');
end


%Update the table container
tableContainer.Sagittal = Sagittal;
setappdata(hFig, 'tableContainer', tableContainer);
updateSubjectSummaryText(handles);

function menu_edit_Callback(hObject, eventdata, handles)
% hObject    handle to menu_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function menu_analyzeFit_Callback(hObject, eventdata, handles)
% hObject    handle to menu_analyzeFit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_analyzeFit4order_Callback(hObject, eventdata, handles)
% hObject    handle to menu_analyzeFit4order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%.........................................................................
%Get subjects
%.........................................................................
hFig = ancestor(hObject, 'figure');
tableContainer = getappdata(hFig, 'tableContainer');

if isempty(tableContainer) || ~isfield(tableContainer, 'Sagittal') || isempty(tableContainer.Sagittal)
    %No data found
    errordlg('No analyses found.')
    return
end

%experimentTable = tableContainer.Experiment;
%subjectTable = tableContainer.Subject;

sagTable = tableContainer.Sagittal;

%Let users select from list dialog which they will use
[selection, ok] = listdlg('ListString',sagTable.ExperimentID,...
    'SelectionMode','multiple', 'ListSize',[240 400],...
    'PromptString','Select files to analyze:','Name','Sagittal Analysis');

if ok && ~isempty(selection)
    sagTable_selected = sagTable(selection, :);
else
    %user didn't select anything
    return
end

sagTable_selected = batchPoly4FitAppend(sagTable_selected);

sagTable(selection,:) = sagTable_selected;

%Update DB
tableContainer.Sagittal = sagTable;
setappdata(hFig, 'tableContainer', tableContainer);

%Make sure you can now use the fields in settings
populateTableList(handles, tableContainer); %popups for adding fields



function sagTable = batchPoly4FitAppend(sagTable)
%This function adds a polynomial fit to the analysis

%Get a waitbar going
fitTablesCell = sagTable.Fits;
coeffTablesCell = sagTable.Coefficients;

nRows = size(fitTablesCell,1);
if nRows > 3
    hWaitbar = makeWaitbar('Analyzing Polynomial Fits', 'Batch Analysis');
    hasWaitbar = true;
else 
    hasWaitbar = false;
end

%For each subject
for i = 1:numel(fitTablesCell)
    
    iFitTable = fitTablesCell{i};
    iCoeffTable = coeffTablesCell{i};
    validTableRows = iFitTable(iFitTable.FitName == 'Conic',:);
%     validTableRows = iTable(iTable.FitName == 'Conic' &...
%         cellfun(@isempty,strfind(lower(iTable.InputDescription),'corrected')),:);
    
    %Update waitbar
    if ~hasWaitbar 
        %do nothing
    elseif ishandle(hWaitbar)
        waitbar(i/nRows, hWaitbar,...
            sprintf('Processing Item: %d of %d', i, nRows));
    else
        break
    end
    
    %For original and corrected
    for m = 1:size(validTableRows, 1)
        mTableRow_fit = validTableRows(m,:);
        mPts = mTableRow_fit.InputPoints{1};
        mFit = mTableRow_fit.FitPoints{1}; %update after you fit the new points
        mInputDesc = mTableRow_fit.InputDescription;
        mFitWindow = mTableRow_fit.FitWindow;
        
        mFitObj = cell(1,length(mPts)); mGOF = mFitObj;
        mR = [NaN, NaN]; %newR
        %For anterior and posterior surfaces
        for j = 1:length(mPts)
            jXY = mPts{j};
            jFitX = mFit{j}(:,1);
            [jFitObj, mGOF{j}] = fit_poly4(jXY(:,1), jXY(:,2), j);
            
            mFit{j}(:,2) = feval(jFitObj, jFitX);
            mFitObj{j} = jFitObj;
            mR(j) = 1/(2*jFitObj.c);
        end
        
        mNewRow_fit = mTableRow_fit;
        mNewRow_fit{:,'FitName'} = categorical(cellstr('Poly4'));
        mNewRow_fit{:,'FitFcn'} = cellstr('a*(x-x0)^4 + b*(x-x0)^3 + c*(x-x0)^2 + d*(x-x0) + e');
        mNewRow_fit{:, 'FitPoints'} = {mFit};
        mNewRow_fit{:,'FitObject'} = {mFitObj};
        mNewRow_fit{:,'Coefficients'} = {{coeffvalues(mFitObj{1}),... %anterior
            coeffvalues(mFitObj{2})}}; %posterior
        
        %if the row for this fit already exists, replace it, otherwise add the
        isHere_fit = iFitTable.FitName == 'Poly4' & ...
            cellfun(@(rowName) strcmpi(mNewRow_fit.InputDescription, rowName), iFitTable.InputDescription);
       
        if any(isHere_fit)
            iFitTable(isHere_fit,:) = mNewRow_fit;
        else
            iFitTable(end+1,:) = mNewRow_fit;
        end
        
        
        
        mNewRow_coeff = iCoeffTable(1,:);
        mNewRow_coeff.Variable(1) = categorical(cellstr('R'));
        mNewRow_coeff.Value = mR;
        mNewRow_coeff.FitName(1) = categorical(cellstr('Poly4'));
        mNewRow_coeff.InputDescription(1) = mInputDesc;
         
        
        isHere_coeff = iCoeffTable.FitName == 'Poly4' &...
            cellfun(@(iStr) strcmpi(iStr, mInputDesc), iCoeffTable.InputDescription);
        
        if any(isHere_coeff)
            iCoeffTable(isHere_coeff,:) = mNewRow_coeff;
        else
            iCoeffTable(end+1,:) = mNewRow_coeff;
        end
        
    end
    
    %Put the table back into the dataset
    fitTablesCell{i} = iFitTable;
    coeffTablesCell{i} = iCoeffTable;
end

if ~hasWaitbar 
    %do nothing
elseif ishandle(hWaitbar)
    delete(hWaitbar)
else
    errordlg('User canceled analysis')
    return
end

sagTable.Fits = fitTablesCell;
sagTable.Coefficients = coeffTablesCell;

helpdlg('Batch fit complete','4th Order Fit')

function [fitObj, gof] = fit_poly4(xdata, ydata, sType)
%This funciton performs a fourth order fit on x and y data using
%Levenberg-Marquadt algorithm

[xData, yData] = prepareCurveData(xdata, ydata);

ft = fittype( 'a*(x-x0)^4 + b*(x-x0)^3 + c*(x-x0)^2 + d*(x-x0) + e', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.StartPoint = [0.001697 0.003 0.09 0.02 -2 0];

if sType == 2 %ant
    opts.StartPoint = -opts.StartPoint;
end

opts.Algorithm = 'Levenberg-Marquardt';

% Fit model to data.
[fitObj, gof] = fit( xData, yData, ft, opts );



function menu_fileConsolidateImg_Callback(hObject, eventdata, handles)
% hObject    handle to menu_fileConsolidateImg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

consolidateImages(handles);



%==========================================================================
%FITS: Forth order
%==========================================================================



%==========================================================================
%SEGMENT: Sagittal
%==========================================================================


function thisAnalysisOut = sagittalSegment(tableContainer,...
    expTableRow, analysisTable, runSilent)
%This function runs semiAutoSegment_shadowgraph on a single row of the
%experiment table and also checks the analysis table to see if 

subjectTable = tableContainer.Subject;

if nargin < 4
    runSilent = false;
end

if nargin < 3
    if isfield(tableContainer,'Sagittal')
        analysisTable = tableContainer.Sagittal;
    else
        analysisTable = [];
    end
end

%if more than one row, only analyze the first
if size(expTableRow, 1) > 1
    expTableRow = expTableRow(1,:);
end

thisAnalysisIn = [];

if ~isempty(analysisTable)
    %Find the corresponding analysis row
    isThisRow = strcmpi(analysisTable.ExperimentID, char(expTableRow.ExperimentID));
    if any(isThisRow)
        thisAnalysisIn = sortrows(analysisTable(isThisRow,:),'CreatedDate','descend');
        thisAnalysisIn = thisAnalysisIn(1,:);
    end
end

thisTableIn = innerjoin(subjectTable, expTableRow);
userData = createUserDataStructure(thisTableIn(1,:), thisAnalysisIn);

try
    thisAnalysisOut = semiAutoSegment_shadowgraph(userData, runSilent); %runSilent = true
catch ME
    %error with the automated fit
    disp(ME.message)
    thisAnalysisOut = [];
end

if isempty(thisAnalysisOut)
    %there was a problem
    return
else
    %add SubjectID, and ExperimentID, CreatedBy, CreatedDate, Notes
    thisAnalysisOut.SubjectID = expTableRow.SubjectID;
    thisAnalysisOut.ExperimentID = expTableRow.ExperimentID;
    
    %add the polynomial fit
    thisAnalysisOut = batchPoly4FitAppend(thisAnalysisOut);
end

function userData = createUserDataStructure(expTable, analysisTable)
%This function creates a userData structure for the semiAutoSegment program
%INPUT:
%expTable - a join of subjectID and experimentTable
%subjectTable is used to for the age field to determine window size
%analysisTable - previous output from the SemiAutoSegment program that contains
%settings you wish to copy (if it is for a different image, set
%the ROI to all true first)

%Initialize
windowSize = 6; %mm

if size(analysisTable,1) > 1
    analysisTable = analysisTable(1,:);
end
    

colNames = expTable.Properties.VariableNames;
if any(ismember(colNames, 'Age'))
    age = expTable.Age;
    
    %If the eye is less than 10 years old, make the window size 4mm
    if age < 10
        windowSize = 4; %mm
    end
end

%..........................................................................
%Read the Image
%..........................................................................
%Find the image file and read it
[img, imgLocation] = readImage(expTable);

xSize = size(img,2);
xyRes = repmat(1/0.0064102564,[1,2]); %pixels/mm
scanWidth = xSize/xyRes(1); %[width, height]

userData.fileLocation = imgLocation;
userData.scanWidth = scanWidth;
userData.imgSize = size(img);
userData.xyRes = xyRes;
userData.imgIn = img;
userData.ID = char(expTable.ExperimentID);

switch expTable.ImgType 
    case 'Sagittal'
        userData.imgType = 5; %shadowgraph sagittal
    case 'Coronal'
        userData.imgType = 6; %shadowgraph coronal
    otherwise 
        userData.imgType = 5; %default to sagittal for now
end


%..........................................................................
%Set parameters
%..........................................................................
if nargin > 1 && ~isempty(analysisTable)
    params.xyRes = xyRes; %should never change
    params.nRecursions = analysisTable.Param_nRecursions;
    params.fitType = str2func(analysisTable.Param_FitType); %should never change (it's not saved)
    params.fitDivFactor =analysisTable.Param_FitDivFactor;
    params.startFitTol = analysisTable.Param_StartFitTol;
    params.fitWindow = windowSize;
    
    nRotations = analysisTable.Param_InitialRotation;
    params.nRotations = nRotations;
    
    userData.params = params;
    userData.variantFcn = analysisTable.Param_VariantFunction;
    
    %convert linear indices to roi mask
    imgRoi_stored = cell2mat(analysisTable.Param_ImgROI);
    if any([size(imgRoi_stored,1),size(imgRoi_stored,2)]==1)
        imgRoi_idx = imgRoi_stored;
        imgRoi = false(size( imrotate(img(:,:,1), nRotations*90) ));
        imgRoi(imgRoi_idx) = true;
    else
        %may not be linear indicies
        imgRoi = imgRoi_stored;
    end
    userData.imgRoi = imgRoi;
    
    zoomRegion = analysisTable.Param_ZoomRegion;
    if ~isempty(zoomRegion)
        userData.zoomRegion = zoomRegion;
    else
        userData.zoomRegion = [];
    end
    
    userData.presetSettings = true;
else
    userData.params.fitWindow = windowSize;
    userData.presetSettings = false; %use defaults
end

function tableFinal = displaySagittalTable(tableContainer)
%Show Sagittal Table

if ~isfield(tableContainer, 'Sagittal') || isempty(tableContainer.Sagittal)
    return
end

%1. Show only the rows of the table with the most recent CreatedDate
%sort by experimentID then created date DESC, unique first occurence
Subject = tableContainer.Subject;
Experiment = tableContainer.Experiment;
Sagittal = tableContainer.Sagittal;

%Make sure experiment and subject rows are unique
[~,ia_E,~]=unique(Experiment.ExperimentID, 'stable');
Experiment = Experiment(ia_E,:);

[~,ia_S,~]=unique(Subject.SubjectID, 'stable');
Subject = Subject(ia_S,:);


tableData_all = sortrows(join(join(Sagittal, Experiment), Subject),...
    {'SubjectID', 'CreatedDate'},'descend');

[~, ia, ~] = unique(tableData_all.ExperimentID,'first');
tableData = tableData_all(ia,:);


%2. Merge the tables: Sagittal, Subject, Epxeriment,  -('Species','ImgType')
outputCols = {'SubjectID','ExperimentID' 'Eye','Age', 'PMT',...
    'LensThickness','Fourier20_LensThickness',...
    'LensDiameter','Fourier20_LensDiameter'};
tableOut = tableData(:, outputCols);

%3. Get the get coefficients from the coefficient table
%radius, asphericity
CoeffTable = tableData(:,'Coefficients');
CoeffVals = ...
    rowfun(@(tableIn_cell) unstack(tableIn_cell{1},'Value','Variable',...
    'AggregationFunction', @(x) x(1,:)),... %give me only the first/last
    CoeffTable,'OutputVariableNames',{'Coefficients'});
R_ant = CoeffVals.Coefficients.R(:,1);
R_pos = CoeffVals.Coefficients.R(:,2);
p_ant = CoeffVals.Coefficients.p(:,1);
p_pos = CoeffVals.Coefficients.p(:,2);
CoeffTableOut = table(R_ant,R_pos,p_ant,p_pos, 'VariableNames',...
    {'Conic_Ra', 'Conic_Rp', 'Conic_Pa', 'Conic_Pp'});


% mat2xl2(cat(1,tableFinal.Properties.VariableNames,table2cell(tableFinal)));
tableFinal = cat(2, tableOut, CoeffTableOut);



%==========================================================================
%LISTBOX CALLBACKS
%==========================================================================


function listbox_subject_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_subject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_subject contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_subject

hFig = ancestor(hObject, 'figure');
hListboxS = handles.listbox_subject; %also hObject
selectedVal = get(hListboxS, 'Value');
% hListboxE = handles.listbox_experiment;

tableContainer = getappdata(hFig, 'tableContainer');

if isempty(tableContainer) || ~isfield(tableContainer, 'Subject') || isempty(selectedVal)
    %this should never happen
    return
end

subjectTable = getappdata(hObject, 'table'); %tableContainer.Subject;

%Get Selected SubjectID
thisSubjectID = subjectTable.SubjectID(selectedVal);
%Alternative: = cellstr(get(hListboxS, 'String'));
%             = allSubjectIDs(selectedVal);

updateExperimentListbox(handles, thisSubjectID);

function listbox_experiment_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_experiment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_experiment contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_experiment

%show the relevnt analysis
updateAnalysisListbox(handles);

displaySingleSubject(handles);


%Update graph selection
selectPointByID(handles);

updateSummaryTable(handles); %TEMP --> THIS DOESNT GO HERE

function listbox_analysis_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_analysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_analysis contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_analysis

hListboxA = hObject; 
megaTable_Sag = getappdata(hListboxA, 'table');

if isempty(megaTable_Sag)
    return
end

%Display the analysis
selectedAnalysis = get(hListboxA,'Value');
displaySagittalAnalysis(handles.axes_img, megaTable_Sag(selectedAnalysis,:));

    
function removeElement(hMenu, ~)
%This function deletes an element from a listbox

hFig = ancestor(hMenu, 'figure');


hListbox = get(hFig, 'CurrentObject');
if ~strcmpi(get(hListbox,'Enable'),'on')
    errordlg('There is nothing to remove.','No Data Found');
end

tokens = regexp(get(hListbox, 'Tag'),'listbox_(.*)','tokens'); 

if isempty(tokens)
    return %i don't recognize the calling parent
end

listboxName = lower(char(tokens{1}));


tableContainer = getappdata(hFig, 'tableContainer');
thisTable = getappdata(hListbox, 'table');
selectedVal = get(hListbox, 'value');


if ~isfield(tableContainer, 'Experiment') || isempty(thisTable) || selectedVal < 1
    return
end

uVerify = questdlg(sprintf('Are you sure you would like to remove this %s from the database?',...
    listboxName),'Remove Element?','Yes','No','No');

if ~strcmpi(uVerify, 'Yes')
    return
end


subTable = tableContainer.Subject; 
allSubIDs_subTable = subTable.SubjectID;

expTable = tableContainer.Experiment;
allSubIDs_expTable = expTable.SubjectID;
allExpIDs_expTable = expTable.ExperimentID;

 
if strcmpi(listboxName, 'subject')
    %remove subject
    selected_SubID = force1D(thisTable.SubjectID(selectedVal),2);
    
    %find all subjects with this ID
    isRemoveRow_subTable =  ismember(allSubIDs_subTable, selected_SubID);
    isRemoveRow_expTable = ismember(allSubIDs_expTable, selected_SubID);
    selected_ExpID = allExpIDs_expTable(isRemoveRow_expTable);
    
    %Subject Table
    if any(isRemoveRow_subTable)
        subTable(isRemoveRow_subTable,:) = [];
        tableContainer.Subject = subTable;
    end
    
else
    selected_ExpID = force1D(thisTable.ExperimentID(selectedVal),2);
end
 
requestingTable = [upper(listboxName(1)),listboxName(2:end)];
tableContainer = removeExperimentFromTable(tableContainer, requestingTable,...
    selected_ExpID, 'ExperimentID');


%store updates
setappdata(hFig,'tableContainer', tableContainer);

%refresh displays
handles = guidata(hFig);
populateSubjectListbox(handles, tableContainer.Subject);
populateTableList(handles);
updateGraphPanel(handles);


%..........................................................................
%listbox manipulation functions
%..........................................................................

function populateSubjectListbox(handles, subjectTableIn)
%This function filters the subject table according to the selected settings

if nargin < 2
    tableContainer = getappdata(ancestor(handles.axes_img,'figure'),...
        'tableContainer');
    if isempty(tableContainer)
        subjectTableIn = [];
    else
        subjectTableIn = tableContainer.Subject;
    end
end

if ~isempty(subjectTableIn)
    [isOk, subjectTable] = filterSubjectTable(handles, subjectTableIn);
end

if ~isOk
    subjectTable = subjectTableIn;
end

hListboxS = handles.listbox_subject;

if isempty(subjectTable) 
    
    clearListbox();
    return
    
end

subjectID_S = subjectTable.SubjectID; %subjectID from the subject table

set(hListboxS, 'String', subjectID_S, 'Value',1, 'Enable','On'); %Select None
setappdata(hListboxS, 'table', subjectTable);

updateExperimentListbox(handles);
updateSubjectSummaryText(handles);

function updateExperimentListbox(handles, subjectID_In)
%This function updates the experiment listbox given a subjectID or multiple
%subjectIDs


if nargin < 2 %By default uses the subject IDs in
    hListboxS = handles.listbox_subject;
    subjectTable = getappdata(hListboxS, 'table');
    selectedVal = get(hListboxS, 'Value');
    subjectID_In = subjectTable.SubjectID(selectedVal);
end

hListboxE = handles.listbox_experiment;
hFig = ancestor(hListboxE, 'figure');

tableContainer = getappdata(hFig, 'tableContainer');
experimentTable = tableContainer.Experiment;


%Filter the table for the subjectID_In
queryStrOut = buildQueryString('SubjectID', subjectID_In);
query = {'SubjectID', queryStrOut};
valid_ExperimentRows = filterTable(experimentTable, query);
validExpID_filtered = getappdata(handles.panel_filter,'ExperimentID');

if ~isempty(validExpID_filtered) && iscell(validExpID_filtered) &&...
        ~all(strcmpi(validExpID_filtered, 'none')) && size(valid_ExperimentRows,1) > 0
    isValidRow = ismember(valid_ExperimentRows.ExperimentID,validExpID_filtered);
    
    valid_ExperimentRows = valid_ExperimentRows(isValidRow,:);
end

subjectID_E = valid_ExperimentRows.SubjectID;
if isempty(subjectID_E) %no valid experiments
    clearListbox([hListboxE, handles.listbox_analysis]);
    set([handles.check_flagAnt, handles.check_flagPos],'Enable','off');
    return
end

%Build string for listbox
listboxStr_table = rowfun(@(iID, iExpID) strrep(iExpID,iID,''),...
    valid_ExperimentRows, 'InputVariables', {'SubjectID', 'ExperimentID'});
listboxStr = cellstr(listboxStr_table.Var1);
%Alternative: listboxStr = valid_ExperimentRows.ExperimentID;

%Populate Listbox (need to update displays based on listbox)
set(hListboxE, 'String', listboxStr,'Value',1,'Enable','on');
setappdata(hListboxE, 'table', valid_ExperimentRows);

%Run the experiment listbox callback
funCall = dbstack;
if ~any(strcmpi({funCall.name},'pointCallback'))
    listbox_experiment_Callback(hListboxE, [], handles)
end

set([handles.check_flagAnt, handles.check_flagPos],'Enable','on');

function updateAnalysisListbox(handles, experimentID_In)

if nargin < 2 %By default uses the subject IDs in
    hListboxE = handles.listbox_experiment;
    experimentTable = getappdata(hListboxE, 'table');
    selectedVal = get(hListboxE, 'Value');
    experimentID_In = experimentTable.ExperimentID(selectedVal);
end

hListboxA = handles.listbox_analysis;
hFig = ancestor(hListboxA, 'figure');

tableContainer = getappdata(hFig, 'tableContainer');

if ~isfield(tableContainer, 'Sagittal')
    clearListbox(hListboxA);
    return
end

analysisTable = tableContainer.Sagittal;


%Filter the table for the experimentID_In
queryStrOut = buildQueryString('ExperimentID', experimentID_In);
query = {'ExperimentID', queryStrOut};
valid_AnalysisRows = filterTable(analysisTable, query);

if isempty(valid_AnalysisRows) %no valid analyses
    clearListbox(hListboxA);
    return
end


%Build string for listbox
listboxStr_table = rowfun(@(iDate) cellstr(datestr(iDate, 'dd|mmm|yy HH:MM PM')),...
    valid_AnalysisRows, 'InputVariables', {'CreatedDate'});
listboxStr = cellstr(listboxStr_table.Var1); 
%Alternative: listboxStr = valid_ExperimentRows.ExperimentID;

%Populate Listbox (need to update displays based on listbox)
set(hListboxA, 'String', flipud(listboxStr),'Value',1,'Enable','on');
setappdata(hListboxA, 'table', flipud(valid_AnalysisRows));

function clearListbox(hListbox)
%This function clears the contents of a listbox

if nargin < 1
    %clear all
    handles = guidata(ancestor(gcbf, 'figure'));
    hListbox = [handles.listbox_subject, handles.listbox_experiment,...
        handles.listbox_analysis];
end 

set(hListbox, 'String','Empty', 'Value',1, 'Max',1,'Enable','off');
arrayfun(@(iH) setappdata(iH, 'table',[]), hListbox);

function updateGraphPanel(handles)
%This function reads settings to update the graph plot

hGraphAx = handles.axes_graph;
hFig = ancestor(hGraphAx,'figure');
tableContainer = getappdata(hFig,'tableContainer');

h_xField = handles.popup_xField;
h_yField = handles.edit_yField;

xFieldData = getappdata(h_xField, 'data');
yFieldStr = get(h_yField, 'String');

hPanel = handles.panel_filter;
validExpIDs = getappdata(hPanel,'ExperimentID');


if isempty(tableContainer) || isempty(xFieldData) || isempty(yFieldStr)...
        || isempty(validExpIDs) || all(strcmpi(validExpIDs,'none'))
    %nothing to do
    cla(hGraphAx);
    setappdata(hGraphAx, 'hPlots',[]);
    setappdata(hGraphAx, 'hExcluded',[]);
    setappdata(hGraphAx, 'yData',[]);
    setappdata(hGraphAx, 'xData',[]);
    setappdata(hGraphAx, 'hLegend',[]);
    setappdata(hGraphAx, 'hHilight',[]);
    set(hGraphAx, 'UIcontextmenu','');
    return
end 

try
    yFieldData = evalFieldExpression(tableContainer, yFieldStr); %Value, Units, Label
catch
    return
end
 
%Plot the data & set callbacks
[hP, newXFieldData, newYFieldData] = plotGraph(hGraphAx, xFieldData, yFieldData);

%Show the outliers
if isempty(newXFieldData)
    return
end

plotExclusions(handles, hP, newXFieldData, newYFieldData);

%Update the context menu
updateGraphContextMenu(handles, newYFieldData);

%Hilight the current point
selectPointByID(handles)

%Turn off the settings flag
setappdata(hFig, 'settingsFlag', false);

function displaySettingChange(hObject, edata, handles)
%This function runs when the user displays ant or post
updateGraphPanel(handles);

function tableContainer = removeExperimentFromTable(tableContainer, requestingTable, expID, keyName)
%This function removes all rows related to a particular experimentID
%INPUT: expID --> cellstr array (multipl allowed)

expID = force1D(expID,2);
allTables = fieldnames(tableContainer);
switch requestingTable
    case {'Subject','Experiment'}
        tableNames = allTables;
    case 'Analysis'
        tableNames = intersect(allTables, 'Sagittal');
    otherwise %remove only from the requesting table
        tableNames = intersect(allTables, requestingTable);
end

if nargin < 4
    keyName = 'ExperimentID';
end

for i = 1:numel(tableNames)
    iTable = tableContainer.(tableNames{i});
    iColNames = iTable.Properties.VariableNames;
    
    if ~ismember(keyName,iColNames)
        continue %this field is not in this table
    end
    
    iKeyVals = iTable.(keyName);
    iIsRemove = ismember(iKeyVals, expID);
    
    if any(iIsRemove)
        iTable_mod = iTable;
        iTable_mod(iIsRemove,:) = [];  
        tableContainer.(tableNames{i}) = iTable_mod;
    end    
end


%..........................................................................
% Graph Plot Functions
%..........................................................................

%add or remove relevant fields of context menu

function [hP, newXFieldData, newYFieldData] = plotGraph(hAx, xFieldData, yFieldData)
%This function creates plots for each graph and associates the point with
%the image file/analysis

handles = guidata(ancestor(hAx, 'figure'));
xLabelStr = xFieldData(1).Label;
xUnits = xFieldData(1).Units;

nPlots = numel(yFieldData);
isAnt = get(handles.check_dispAnt, 'Value');
hP = zeros(1,nPlots)-1; %initialize handles


newXFieldData = [];
newYFieldData = [];



legendStr = {};

%Clear the axes and get ready to plot again
cla(hAx, 'reset');
hold(hAx, 'on');

%Filter the data, show only the desired data
hPanel = handles.panel_filter;
validExpIDs = getappdata(hPanel,'ExperimentID');

if isempty(validExpIDs) || all(strcmpi(validExpIDs,'none')) %nothing should be showing
    return
end

%also get the subject IDs
subTable = getappdata(handles.listbox_subject, 'table');

if isempty(subTable)
    return
end
validSubIDs = subTable.SubjectID;

dataFlag = false(nPlots, 1);

for i = 1:nPlots
    
    [iXData, iYData] = matchXandYData(handles, xFieldData(1), yFieldData(i));    
    
    
    if isempty(iXData)
        fprintf(1,'\n\nError with plot: %s vs %s\n\n',xFieldData(1).Label,...
            yFieldData(1).Label);
        continue %there is a problem
    end
    
    %get the points that correspond to valid IDs
    if strcmpi(iXData.KeyName, 'SubjectID')
        thisValidID = validSubIDs;
    else %keyName = ExperimentID
        thisValidID = validExpIDs;
    end
    
    %This separates outliers from points
    [iXData_new, iYData_new] = filterDataStructure(iXData, iYData,...
        thisValidID, isAnt);
    
    if isempty(iXData_new)
        dataFlag = false;
        continue %no data here
    end
    
    if isempty(newXFieldData)
        newXFieldData = iXData_new;
        newYFieldData = iYData_new;
    else
        newXFieldData(end+1) = iXData_new;
        newYFieldData(end+1) = iYData_new;
    end
    
    if isempty(iXData_new.Value)
        continue %no valid points to plot
    end
    
    iX = iXData_new.Value; iY = iYData_new.Value;
    
    
    hP(i) = plot(hAx, iX, iY,'.','MarkerSize',10);
    setappdata(hP(i),'xData', iXData_new);
    setappdata(hP(i),'yData', iYData_new);
    
    %Callback when the line is clicked
    set(hP(i),'ButtonDownFcn', @pointCallback);
    
    %LegendStr
    iUnits = yFieldData(i).Units;
    if ~isempty(iUnits)
        iLegendStr = sprintf('%s (%s)',yFieldData(i).Label, iUnits);
    else
        iLegendStr = yFieldData(i).Label;
    end
    
    legendStr{end+1} = strrep(iLegendStr,'_',' ');
end

grid(hAx, 'on');

xlabel(sprintf('%s (%s)', xLabelStr, xUnits));
hold(hAx, 'off');


if ~isempty(legendStr)
    hLegend = legend(hAx, legendStr);
    setappdata(hAx, 'hLegend', hLegend);
end

hP(dataFlag) = []; %this plot doesn't exist
setappdata(hAx, 'hPlots', hP);
setappdata(hAx, 'xData', newXFieldData); %Data structure in the order it is plotted
setappdata(hAx, 'yData', newYFieldData); %Data structure in the order it is plotted

function plotExclusions(handles, hP, xFieldData, yFieldData)
%This function plots the exclusions on the graph plot


hAx = ancestor(hP(1), 'axes');

if isempty(hAx)
    hAx = handles.axes_graph; %should not happen
end

outlierData = {xFieldData.Outliers};

if all(cellfun(@isempty,outlierData))
    %no outliers
    return
end

hCheck = handles.check_showExclusions;
if ~get(hCheck, 'Value')
    %no need to show outliers
    return
end
 
nPlots = numel(yFieldData);

hold(hAx, 'on');

hP_excluded = zeros(1,nPlots)-1; 

for i = 1:nPlots
    iHPlot = hP(i);
    
    thisXData = xFieldData(i);
    thisYData = yFieldData(i);
    iX_outlier = thisXData.Outliers; 
    iY_outlier = thisYData.Outliers; 
    
    if isempty(iX_outlier)
        continue
    end
    
    hP_excluded(i) = plot(hAx, iX_outlier, iY_outlier, 'xr','MarkerSize',8);
    
    
    %Switch values and outliers for callback
    xOutlierKey = thisXData.OutlierKeyValue;
    yOutlierKey = thisYData.OutlierKeyValue;
    
    thisXData.Outlier = thisXData.Value;
    thisXData.Value = iX_outlier; 
    thisXData.OutlierKeyValue = thisXData.KeyValue;
    thisXData.KeyValue = xOutlierKey;
    
    
    thisYData.Outlier = thisYData.Value;
    thisYData.Value = iY_outlier;
    thisYData.OutlierKeyValue = thisYData.KeyValue;
    thisYData.KeyValue = yOutlierKey;
    
    setappdata(hP_excluded(i), 'xData', thisXData);
    setappdata(hP_excluded(i), 'yData', thisYData);
    
    %Callback when the line is clicked
    set(hP_excluded(i),'ButtonDownFcn', @pointCallback);
    
    %Save the handle
    if ishandle(iHPlot)
        setappdata(iHPlot, 'hExcluded', hP_excluded);
    end
end

hold(hAx, 'on');

setappdata(hAx, 'hExcluded', hP_excluded)

%..........................................................................
%Point callbacks
%..........................................................................
function pointCallback(hObj, ~)
%This function is called when the user clicks a line

hAx = ancestor(hObj, 'axes');
hFig = ancestor(hAx, 'figure');
handles = guidata(hFig);

%make this point active
pointIdx = hilightPoint(hObj);

if isempty(pointIdx)
    return
end

%UPDATE the PANEL VIEWS, run the selection change function, select
thisXData = getappdata(hObj,'xData');
thisYData = getappdata(hObj,'yData'); 

tableContainer = getappdata(hFig, 'tableContainer');
expTable = tableContainer.Experiment;

if strcmpi(thisYData.KeyName,'SubjectID')
    %it won't be specific to experiment
    subjectID = thisXData.KeyValue{pointIdx};
    
    isThisRow = strcmpi(expTable.SubjectID, subjectID);
    expRow = expTable(isThisRow,:);
    experimentID = expRow.ExperimentID{1}; %char array
else %'ExperimentID'
    experimentID = thisXData.KeyValue{pointIdx}; %char array
    
    isThisRow = strcmpi(expTable.ExperimentID, experimentID);
    expRow = expTable(isThisRow,:);
    subjectID= expRow.SubjectID{1};
end


selectListboxID(handles, subjectID, experimentID);

%Update image displays: this function should update table too
displaySingleSubject(handles, expRow(1,:));

%show the relevnt analysis
updateAnalysisListbox(handles, experimentID);

%update summary table
updateSummaryTable(handles, experimentID)

%if the user double clicks, SWITCH TAB
x=1; %put the code if you want this functionality

function idx = hilightPoint(hLine, xyPos)
%This function makes a point on a given line active

%Initialize
idx = [];
if ~ishandle(hLine)
    %nothing to see here
    return
end

hAx = ancestor(hLine, 'axes');
hHilight = getappdata(hAx, 'hHilight');

if nargin < 2
    lastClickPos = get(hAx, 'CurrentPoint');
else
    lastClickPos = xyPos;
end
clickX = lastClickPos(1,1); clickY = lastClickPos(1,2);
x = get(hLine, 'XData'); y = get(hLine, 'YData');

%snap to vertex
[~, idx] = min(sqrt((x-clickX).^2 + (y-clickY).^2));
thisX = x(idx); thisY = y(idx);


if ( isempty(hHilight) || ~ishandle(hHilight) )
    
    hold(hAx, 'on');
    %Create it again
    hHilight = plot(hAx, thisX, thisY, 'o','Color',[153 0 120]./255,...
        'MarkerSize', 6);
    hold(hAx, 'off');
    %save it
    setappdata(hAx, 'hHilight', hHilight);
else
    %change the data only
    set(hHilight,'XData', thisX, 'YData',thisY,'Visible','on');
end

%..........................................................................

function [xOut, yOut] = matchXandYData(handles, xFieldData, yFieldData)
%This function pairs up the x and y field data, both x and y field
%structures must be scalars

%Initialize Output
xOut = []; yOut = [];

if isempty(yFieldData.Value) || isempty(xFieldData.Value)
    %There is NO Data
    return
end


%Get the information you need
tableContainer = getappdata(ancestor(handles.axes_img,'figure'), 'tableContainer');
expTable = tableContainer.Experiment;
lookUpTable = expTable(:,{'ExperimentID', 'SubjectID','OutlierFlag'});


fieldKey_X = xFieldData.KeyName;
fieldKey_Y = yFieldData.KeyName;

xTable = table(xFieldData.KeyValue, xFieldData.Value,'VariableNames',{'KeyX', 'XValue'});
yTable = table(yFieldData.KeyValue, yFieldData.Value, 'VariableNames',{'KeyY', 'YValue'});


if strcmpi(fieldKey_X, fieldKey_Y)
    if(size(xFieldData.Value,1) ==size(yFieldData.Value,1))
        
        mergedTable = innerjoin(yTable, xTable,'LeftKeys',...
            'KeyY', 'RightKeys','KeyX');
        
        if strcmpi(fieldKey_X,'SubjectID')
            outlierFlag = false(size(mergedTable,1),2); %no outliers possible
        else
            mergedTable2 = join(mergedTable, lookUpTable,'LeftKeys',...
                'KeyY', 'RightKeys', fieldKey_X);
            outlierFlag = mergedTable2{:,'OutlierFlag'};
        end
        
        xOut = xFieldData;
        xOut.Value = mergedTable{:, 'XValue'};
        xOut.OutlierFlag = outlierFlag;
        
        yOut = yFieldData;
        yOut.Value = mergedTable{:, 'YValue'};
        yOut.OutlierFlag = outlierFlag;
    end
    %if they arent the same size, we have a major issue
    return
end



if strcmpi(fieldKey_Y, 'ExperimentID')
    
    firstTable = yTable;
    secondTable = xTable;
    key1 = 'KeyY';
    key2 = 'KeyX';
    
elseif strcmpi(fieldKey_X, 'ExperimentID')
    firstTable = xTable;
    secondTable = yTable;
    key1 = 'KeyX';
    key2 = 'KeyY';
    
    
else
    %this should never happen.
    return
end

mergedTables = join(join(firstTable, lookUpTable,'LeftKeys',key1,...
        'RightKeys','ExperimentID'),secondTable,'LeftKeys','SubjectID',...
        'RightKeys',key2);

xOut = xFieldData; 
xOut.Value = mergedTables{:,'XValue'};
xOut.OutlierFlag = mergedTables{:,'OutlierFlag'};
xOut.KeyName = 'ExperimentID';
xOut.KeyValue = mergedTables{:,key1};


yOut = yFieldData; 
yOut.Value = mergedTables{:,'YValue'};
yOut.OutlierFlag = mergedTables{:,'OutlierFlag'};
yOut.KeyName = 'ExperimentID';
yOut.KeyValue = mergedTables{:,key1};

function [xOut, yOut] = filterDataStructure(xData, yData, validIDs, isAnt)
%This function filters and removes outliers of xFieldData and yFieldData
%INPUTS:
    %handles
    %xFieldDataIn--> x output of matchXandY function (SCALAR STRUCTURE)
    %yFieldDataIn--> y output of matchXandY function (SCALAR STRUCTURE)
    %validIDs --> validIDs filtered
    %isAnt --> whether we want the anterior or not


%Determine which rows are valid
isRowValid = ismember(xData.KeyValue,validIDs);

if ~any(isRowValid)
    xOut = [];
    yOut = [];
    return
end

%Determine the column index

xVal = xData.Value;
if ~isAnt && size(xVal,2) > 1
    colIdx_x = 2;
else
    colIdx_x = 1;
end

yVal = yData.Value;
if ~isAnt && size(yVal,2) > 1
    colIdx_y = 2;
else
    colIdx_y = 1;
end
    
%DetermineOutliers
allOutliers = yData.OutlierFlag;

if size(yVal,2) == 1
    isOutlier = any(allOutliers,2);
else 
    isOutlier = allOutliers(:, colIdx_y);
end
 
isKeepPoint = bitand(isRowValid, ~isOutlier);
isKeepOutlier = bitand(isRowValid, isOutlier);

xOut = rmfield(xData, {'OutlierFlag'});
xOut.Value = xVal(isKeepPoint, colIdx_x);
xOut.KeyValue = xData.KeyValue(isKeepPoint);
xOut.Outliers = xVal(isKeepOutlier,colIdx_x);
xOut.OutlierKeyValue = xData.KeyValue(isKeepOutlier);

yOut = rmfield(yData, {'OutlierFlag'});
yOut.Value = yVal(isKeepPoint, colIdx_y);
yOut.KeyValue = yData.KeyValue(isKeepPoint);
yOut.Outliers = yVal(isKeepOutlier,colIdx_y);
yOut.OutlierKeyValue = yData.KeyValue(isKeepOutlier);


function selectPointByID(handles, experimentID)
%This function highlights a point on the graph given an ID

hAx = handles.axes_graph;
hPlots = getappdata(hAx, 'hPlots');
hListboxE = handles.listbox_experiment;

if nargin < 2
    %display selected subject in experiment 
    
    selectedVal = get(hListboxE, 'Value');
    expTable = getappdata(hListboxE, 'table');
    
    if isempty(selectedVal) || isempty(expTable)
        return
    end
    
    
    expRow = expTable(selectedVal, :);
    experimentID= expRow.ExperimentID;
end


if isempty(hPlots) || any(~ishandle(hPlots))
    return
end

nPlots = numel(hPlots);
xFieldData = getappdata(hAx, 'xData');
yFieldData = getappdata(hAx, 'yData');

hasExpID_Y = strcmpi({yFieldData.KeyName},'ExperimentID');
hasExpID_X = strcmpi({xFieldData.KeyName},'ExperimentID');

if any(hasExpID_Y)
    expIDs = yFieldData(find(hasExpID_Y,1,'first')).KeyValue;
    
elseif any(hasExpID_X)
    expIDs = xFieldData(find(hasExpID_X,1,'first')).KeyValue;
else
    %no expIDs, must use subjectIDs, very rare case
    subjectIDs = xFieldData(1).KeyValue;
    tableContainer = getappdata(ancestor(hAx, 'figure'),'tableContainer');
    if isempty(tableContainer)
        return
    else 
        return %no good solution for now
    end
end

isThisID = strcmpi(expIDs, experimentID);

if ~any(isThisID)
    hHilight = getappdata(hAx, 'hHilight');
    if ishandle(hHilight);
        set(hHilight, 'Visible','off');
    end
    return
end


for j = 1:nPlots
    jX = xFieldData(j).Value(isThisID);
    jY = yFieldData(j).Value(isThisID);
    jY = jY(1);%just make sure it's only one point
    jX = jX(1);
    
    hilightPoint(hPlots(j), [jX, jY]);
end

function updateGraphContextMenu(handles, yFieldData)
%This function updates the context menu for the graph plot given the data
%structure for y fields. The structure has fields: Label, Units, and Data

hCmenu = handles.cmenu_graph;
hKids = get(hCmenu, 'Children');
hGraphAx = handles.axes_graph;
iHPlots = getappdata(hGraphAx, 'hPlots');

if nargin < 2
    yFieldData = getappdata(hGraphAx, 'yData');
end


if isempty(yFieldData)
    delete(hKids);
    set(hGraphAx, 'UIcontextmenu','');
    return
end


existingMenuLabels = get(hKids, 'Label');
newLabels = force1D({yFieldData.Label});

%check to see if the menus are exactly the same
if numel(newLabels) == numel(hKids) &&...
        all(strcmpi(sort(existingMenuLabels), sort(newLabels)))
   %if they are the same, make sure the data exists in the same order 
   [~, menuIdx] = sort(existingMenuLabels);
   [~, dataIdx] = sort(newLabels);
   
   newKids = hKids(menuIdx);
   yFieldData = yFieldData(dataIdx);
   
   %Update the containers for  data & menus
   setappdata(hGraphAx, 'yData', yFieldData);
   set(hCmenu, 'Children',newKids);
else
    %recreate it
    delete(hKids);
    newKids = zeros(numel(newLabels),1) - 1; %initialize handles
    
    for j = 1:numel(newLabels)
        newKids(j) = uimenu('Parent', hCmenu, 'Label', newLabels{j},...
            'Callback', @(hObj, edata) togglePlotVisibility(hObj),...
            'Checked','on','Enable','on'); %Visible by default
        
    end
    
end

if ~isempty(iHPlots)
    isDisable = ~ishandle(iHPlots);
    set(newKids(isDisable),'Enable','off','Checked','off');
end


set(hGraphAx, 'UIcontextmenu',hCmenu);

function togglePlotVisibility(hObject, newState)
%This funciton turns a plot on or off depending on the current state

menuIdx = get(hObject, 'Position');
state = get(hObject, 'Checked'); %you have to manually toggle the state

handles = guidata(ancestor(hObject, 'figure'));
hGraphAx = handles.axes_graph;
hP = getappdata(hGraphAx, 'hPlots');


if isempty(hP) || numel(hP) <  menuIdx
    %Cant find the plot handle or it is disabled
    return
end

if strcmpi(get(hObject, 'Enable'),'off')
    state = 'on'; %turn it off
end

if nargin < 2
    switch state
        case 'on' %currently on, turn it off
            newState = 'off';
        case 'off' %currently off, turn it on
            newState = 'on';
    end
end

%uncheck the menu, update the: menu check mark, plot visibility, 
set(hP(menuIdx), 'Visible',newState);
set(hObject, 'Checked',newState);


hExcluded = getappdata(hP(menuIdx), 'hExcluded'); %handle passed to parent plot

if ~isempty(hExcluded)
    set(hExcluded(1), 'Visible', newState);
end

function updateGraphLegend(hGraphAx)
%This function is not finished or necessary
handles = guidata(ancestor(hGraphAx, 'figure'));
yFieldData = getappdata(hGraphAx,'yData');
hKids = get(handles.cmenu_graph, 'Children');
isChecked = strcmp(get(hKids, 'Checked'),'on');

legendStr = yFieldData(isChecked);
currentStr = cellstr(get(hLegend,'String'));


hPlots = getappdata(hGraphAx, 'hPlots');
isVisible = strcmpi(get(hPlot, 'Visible'),'on');
cellfun(@(iP) getappdata(iP, ''))

if isempty(currentSr)
    newStr = yFieldData(menuIdx.Label);
else
    switch menuIdx
        case 1 %its the first
        case numel(currentString)
    end
    currentStr = newStr
end

legend(hGraphAx,'String',newStr)

function isOk = selectListboxID(handles, subjectID, experimentID)
%This function selects a subject in the listbox given an ID

%Initialize Variable
isOk = true;

hListboxS = handles.listbox_subject;

subTable = getappdata(hListboxS, 'table');

if isempty(subTable)
    %nothing exists
    isOk = false;
    return
end

%Search the subject listbox, because it may be filtered
isThisSubject = strcmpi(subTable.SubjectID, subjectID);

%See if you found the subject in the current listbox view
if ~any(isThisSubject)
    isOk = false;
    return
end

subjectIdx = find(isThisSubject);

%Run the listbox callback
set(hListboxS, 'Value', subjectIdx(1));
listbox_subject_Callback(hListboxS, 'manualcall', handles);


%..........................................................................
%Check ExperimentID
%..........................................................................
if nargin < 3
    return
end

hListboxE = handles.listbox_experiment;
expTable = getappdata(hListboxE, 'table');

if isempty(expTable)
    isOk = false;
    return %nothing to select
end

isThisExpID = strcmpi(expTable.ExperimentID, experimentID);

if ~any(isThisExpID)
    isOk = false;
    return
end

expIdx = find(isThisExpID);
set(hListboxE, 'Value', expIdx(1));




%-------------------------------------------------------------------------
%FILE MANIPULATION
%-------------------------------------------------------------------------

function [tableContainerOut, outputFlag] = parseFileDirectory(dirLocations,...
    isLoadSubjectInfo)
%This function takes parent directory names as input and extracts file
%information basics.  It then extracts the Excel Data, and converts the
%information to database table format
%INPUT
%dirLocations - cell containing full directory locations for files
%being loaded
%isLoadSubjectInfo - logical scalar: determines whether or not to read
%Excel File

nDirs = numel(dirLocations);
outputFlag = false(nDirs, 1);

%Create a waitbar
hWaitbar = makeWaitbar('Load Image Directory', 'Batch Load');

tableContainerOut = struct('Subject',[], 'Experiment',[],...
    'LegacyCoronal',[],'LegacySagittal',[]);


for i=1:nDirs
    
    iFullDir = dirLocations{i};
    
    %Update waitbar
    if ishandle(hWaitbar)
        thisFolder = ancestorDir(iFullDir,0);
        waitbar(i/nDirs, hWaitbar,...
            sprintf('Loading Directory: %s (%d of %d)',...
            thisFolder, i, nDirs));
    else
        break
    end
    
    
    %----------------------------------------------------------------------
    %Excel File
    %----------------------------------------------------------------------
    %Check to see if it has an excelfile
    [hasXLFile, XL_fullpath, nXLFiles] = getXLfilename(iFullDir, true);
    
    
    %use the excel file to read subject information
    if isLoadSubjectInfo && hasXLFile %state: 11
        
        if nXLFiles > 1
            [legacyTableCell_raw, legacyFlags_cell] = cellfun(@(iFile) readXLfile(iFile), XL_fullpath,...
                'UniformOutput', false);
            legacyTable_xl = cat(1, legacyTableCell_raw{:}); %two rows
            legacyFlags_xl = cat(1, legacyFlags_cell{:});
        else
            [legacyTable_xl, legacyFlags_xl] = readXLfile(XL_fullpath);
        end
        
    elseif hasXLFile %state: 10
        %THIS CODE IS NOT WRITTEN YET, SO CONTINUE
        
        
        if nXLFiles > 1
            [legacyTable_xl, legacyFlags_xl] = readXLfile(XL_fullpath{1},'Eye Data');
        else
            [legacyTable_xl, legacyFlags_xl] = readXLfile(XL_fullpath,'Eye Data');
        end
        
    else %state 00, 01
        outputFlag(i) = true;
        legacyTable_xl = readXLfile('',{'Eye Data'});
        [subjectID, eyePos] = parseEyeDataFromDir(thisFolder);
        legacyTable_xl.SubjectID = cellstr(upper(subjectID));
        
        if ~isempty(eyePos)
            legacyTable_xl.Eye = deal(categorical(eyePos));
        end
        
    end
    
    %----------------------------------------------------------------------
    %Get the image
    %----------------------------------------------------------------------
    imgType = {'Sagittal', 'Coronal'};
    imgTable = parseImgsFromDir(iFullDir, imgType, true); %isReturnAll = true
    
    if isempty(imgTable)
        continue
    end
    
    %This currently works ONLY with pre-analyzed images
    iTableContainer = convertLegacyTable_xl(legacyTable_xl, imgTable, isLoadSubjectInfo);
    
    if i == 1
        tableContainerOut = iTableContainer;
    else
        tableContainerOut = mergeTableContainers(tableContainerOut,...
            iTableContainer, true); %isForceUnique= true
    end
    
end


if ishandle(hWaitbar)
    delete(hWaitbar);
else
    errordlg('User Cancelled the Directory Load.');
    outputFlag = true;
end

function [subjectID, eyePos] = parseEyeDataFromDir(parentFolder)
%This function tries to search the parent folder for eye data


subjectID = regexp(parentFolder, '\w*-\d*','match');
eyePos = regexp(upper(parentFolder), 'O[S/D]','match');

if isempty(eyePos)
    eyePos = {'OX'}; %unknown Eye
end

if ~isempty(subjectID) && ~isempty(eyePos)
    subjectID = strrep(parentFolder, eyePos,'');
else
    subjectID = parentFolder;
end


%..........................................................................
% Extract Images from Directory
%..........................................................................

function [hasImg, img_fullpath, nImgs] = getImageFromDir(parentDir, imgType, isReturnAll)
%This function is used to check a directory for a Sagittal or Coronal
%Image If more than one file is found, it returns the first. The image must
% be jpg
%INPUT:
%parentDir - string: 'rootDir\parentFolder'
%imgType - string: 'Sagittal' or 'Coronal', not case sensitive
%isReturnAll - logical scalar: number of files to return
%OUTPUT:
%hasImg - logical
%img_fullpath - string: 'rootDir\parentFolder\imgFileName.jpg'


if nargin < 2
    isReturnAll = false; %%
end

%Initialize variables
hasImg = true;
img_fullpath = [];
nImgs = 0;

%first check for an excel file
iImgFile = dir(fullfile(parentDir, sprintf('*%s*.jpg',imgType)));

if isempty(iImgFile) %check for another image file using a variation of image type
    
    %check for the image with different case
    allJPGs = dir(fullfile(parentDir, '*.jpg'));
    
    if isempty(allJPGs)
        hasImg = false;
        return
    end
    
    isThisImgValid = cellfun(@(iJPG) ~isempty(strfind(lower(iJPG),...
        lower(imgType(1:3)))), {allJPGs.name}');
    
    if ~any(isThisImgValid)
        hasImg = false;
        return
    end
    
    iImgFile = allJPGs(isThisImgValid);
    
end

%print a warning letting the user know that multiple files were
%found and you are using the first
nImgs = numel(iImgFile);
if nImgs > 1 && ~isReturnAll
    delimStr = '------------------------------------------------';
    fprintf(1,['\n\n%s\nWARNING: Mulitple Images files found\n%s\n\n',...
        'Folder: %s\n\nFiles:\n\n%s']', delimStr, delimStr, ancestorDir(parentDir,0));
    
    cellfun(@disp,{iImgFile.name}');
    fprintf(1,'\nThe first file will be used: \n%s\n\n%s\n', iImgFile(1).name, delimStr);
end

%what if there are multiple XL files??? just take the first...
if isReturnAll && nImgs > 1
    img_fullpath = arrayfun(@(i) fullfile(parentDir, iImgFile(i).name),...
        1:nImgs, 'UniformOutput', false);
else
    img_fullpath = fullfile(parentDir, iImgFile(1).name);
end

function imgTableOut = parseImgsFromDir(parentDir, imgType, isReturnAll)
%This function reads images from the parent directory and returns them as a
%table
%INPUT:
%parentDir: string - directory to search images for imgType
%string or cell: containing the keywords to search the filename
%'Example: Sagittal' or {'Sagittal', 'Coronal'}
%OUTPUT:
%imgTableOut: table - m rows with the columns (ImgType, Img_originalLocation)


imgTableOut = [];

if nargin < 3
    isReturnAll = true;
end


if iscell(imgType)
    nImgTypes = numel(imgType);
else
    imgType = {imgType};
    nImgTypes = 1;
end

for i = 1:nImgTypes
    
    iImgType = imgType{i};
    
    [iHasImg, iImg_fullpath, nImgs] = getImageFromDir(parentDir, iImgType, isReturnAll);
    
    
    if ~iHasImg
        continue
    end
    
    if nImgs > 1
        iTable = table(categorical(repmat({iImgType}, [nImgs,1])),force1D(iImg_fullpath,2),...
            'VariableNames',{'ImgType','Img_originalLocation'});
    else
        
        iTable = table(categorical({iImgType}),{iImg_fullpath},...
            'VariableNames',{'ImgType','Img_originalLocation'});
    end
    
    iTable.ImgNo = (1:nImgs)';
    
    if isempty(imgTableOut)
        imgTableOut = iTable;
    else
        imgTableOut = cat(1, imgTableOut, iTable);
    end
    
end




%-------------------------------------------------------------------------
% TABLE MANIPULATION
%-------------------------------------------------------------------------

function containerOut = mergeTableContainers(tableContainer1,...
    tableContainer2, isForceUnique)
%This function merges the tables in container2 with container1.  A container is a
%structure with tables.
%INPUT
%tableContainer(1/2) - a structure containing tables.  The fieldnames in
%the containers must match
%isForceUnique - boolean (whether or not the output can have duplicates)
%OUTPUT
%containerOut - tableConteiner1 with all of the information in
%tableContainer2 added

if nargin < 3
    isForceUnique = true;
end

names2 = fieldnames(tableContainer2);
containerOut = tableContainer1;

for i = 1:numel(names2)
    iField = names2{i};
    
    if ~isfield(tableContainer1, iField) || isempty(tableContainer1.(iField))
        %It's not there, add it
        containerOut.(iField) = tableContainer2.(iField);
    else
        
        oldTable = tableContainer1.(iField);
        newTable = tableContainer2.(iField);
        
        if isForceUnique
            try
                mergedTable = union(oldTable,newTable); %try union first
            catch
                mergedTable = joinTableRows({oldTable,newTable}); %if error, try something else
            end
        else
            mergedTable = cat(1, oldTable, newTable); %union = alternative to cat
        end
        
        containerOut.(iField) = mergedTable;
    end
end

function mergedTableOut = joinTableRows(tablesIn)
%This function returns only unique rows in the tables, but the tables must
%have the same number of rows and this should only be used in a catch
%statement if the union function returns an error

%join all the tables together
mergedTable = cat(1, tablesIn{:});

%find the table with the most columsn
colNames = mergedTable.Properties.VariableNames;
[nRows, nCols] = size(mergedTable);

isThisColProblematic = false(1,nCols);
isThisRowDuplicate = nan(nRows,nCols);


for i = 1:nCols
    thisCol = colNames{i};
    
    try
        unique(mergedTable(:, thisCol));
        
    catch ME
        isThisColProblematic(i) = true;
        
        if strcmpi(ME.identifier, 'MATLAB:table:VarUniqueMethodFailed')
            %find out if there are any duplicates
            %if you have xy points in a cell,
            [isDuplicate, Xout, iX] = findDuplicates1D(mergedTable{:,thisCol});
            isThisRowDuplicate(:, i) = force1D(isDuplicate,2);
            
        else %if we don't know what the problem is, keep rolling
            rethrow(ME);
        end
        
    end
end

%For each column check to see if the values are unique
%Remove problematic Columns

if any(isThisRowDuplicate(:)) %it's highly unlikely that there will be duplicates
    %do something
    mergedTable_noProblems = mergedTable(:, ~isThisColProblematic);
    [~, ia, ~] = unique(mergedTable_noProblems,'stable');
    
    isRowDuplicate_unique = true(nRows,1);
    isRowDuplicate_unique(ia) = false;
    
    isKeepRow = ~any(cat(2,isRowDuplicate_unique,isThisRowDuplicate),2);
    
    mergedTableOut = mergedTable(isKeepRow,:); %this may be wrong though...
    %sanityCheck = numel(ic) == sum(isKeepRow);
else
    %output is plain cat
    mergedTableOut= mergedTable;
end

function [isDuplicate, Xout, iX] = findDuplicates1D(X)
%X is a cell array of whatever you'd like to compare, it will ONLY work
%properly on 1D arrays
%INPUT:
%X - 1D cell array of anything that can be compared with @isequaln
%OUTPUT
%isDuplicate - logical vector the size = numel(X)
%Xout - cell output of unique rows
%iX - index number of unique elements


nElements = numel(X);


Xout = {}; iX = [];
isDuplicate = true(nElements, 1);
for i = 1:nElements
    thisElement = X{i};
    
    iIsDuplicate = false;
    
    for j = 1:numel(Xout)
        jIsEqual = isequaln(Xout{j}, thisElement);
        
        if jIsEqual
            iIsDuplicate = true;
            break
        end
    end
    
    isDuplicate(i) = iIsDuplicate;
    
    if ~iIsDuplicate
        Xout = cat(1, Xout, X{i});
        iX = cat(1, iX, i);
    end
    
end

function tableContainer = convertLegacyTable_xl(legacyTable_xl, imgTable, isLoadSubjectInfo)
%This function converts the raw excel legacy data to separate tables needed
%for the database.  The input is the output of the @readXLfile function.
%INPUT:
%legacyTable_xl: table with multiple rows if there are multiple files
%imgTable: table with multiple rows for sagittal and coronal images
%imgType, subjectID, experimentDate, eyePos
%isLoadSubjectInfo: logical scalar, whether or not to load the legacy
%data
%OUTPUT
%tableContiner: structure with fields:
%Subject
%SubjectID, Age, Species
%Experiment
%ExperimentID, SubjectID, PMT, ExperimentDate, ImgType,
%Img_originalLocation, Img_filename, Img_serverDirectory,
%Img_localDirectory, AcquiredBy
%LegacyCoronal
%SubjectID, ExperimentID, Coronal section of Legacy Table
%LegacySagittal
%SubjectID, ExperimentID, All other sections of Legcay Table


if nargin < 3
    isLoadSubjectInfo = true;
end

%Create the Subject Table
subjectTable = legacyTable_xl(:,{'SubjectID', 'Age','Species'});
subjectTable{:,'Gender'} = categorical({''}); %Initialize Gender Field
subjectTable = writeTableMetadata(subjectTable, 'Gender','Gender of subject', '');


%Create Experiment Table
imgTableHybrid = imgTable;
imgTableHybrid(:,'SubjectID') = legacyTable_xl(1,'SubjectID');
imgTableHybrid(:,'ExperimentDate') = legacyTable_xl(1,'ExperimentDate');
imgTableHybrid(:,'Eye') = legacyTable_xl(1,'Eye');
imgTableHybrid(:,'PMT') = legacyTable_xl(1,'PMT');


imgTableHybrid = copyTableMetadata(legacyTable_xl, imgTableHybrid);
experimentTable = generateExperimentTable(imgTableHybrid);


%Create Legacy Table
coronal_expID = experimentTable(experimentTable{:, {'ImgType'}}=='Coronal',...
    {'ExperimentID'});
sagittal_expID = experimentTable(experimentTable{:, {'ImgType'}}=='Sagittal',...
    {'ExperimentID'});

tableContainer.Subject = subjectTable;
tableContainer.Experiment = experimentTable;

%assume all legacy files come from the same image
if ~isLoadSubjectInfo
    return %don't append legacy information
    
elseif ~isfield(table2struct(legacyTable_xl),'Fourier10Coeff')
    legacyTable_coronal = [];
    legacyTable_sagittal = [];
else
    coronal_colNames = {'SubjectID','Coronal_AvgR', 'Coronal_MaxR',...
        'Coronal_MinR', 'Coronal_StdR', 'Coronal_CSA','Coronal_Circumference'};
    legacyTable_coronal = legacyTable_xl(:, coronal_colNames);
    
    if isempty(coronal_expID) %There is no image for this data
        legacyTable_coronal{:, 'ExperimentID'} = strrep(sagittal_expID.ExperimentID(1),...
            'Sag','CorX');
    else
        legacyTable_coronal(:, 'ExperimentID') = coronal_expID(1, 'ExperimentID');
    end
    
    
    sagittal_colNames = {'SubjectID', 'OriginalPoints','CorrectedPoints',...
        'Fourier10Coeff', 'Fourier20Coeff', 'Fourier10RMSE', 'Fourier20RMSE',...
        'Fourier10_AnteriorThickness', 'Fourier10_PosteriorThickness',...
        'Fouier10_LensThickness', 'Fourier10_LensDiameter', 'Fourier10_CSA',...
        'Fouier10_SurfaceArea','Fourier10_Volume',...End Fourier10
        'Fourier20_AnteriorThickness', 'Fourier20_PosteriorThickness',...
        'Fourier20_LensThickness','Fourier20_LensDiameter','Fourier20_CSA',...
        'Fouier20_SurfaceArea','Fourier20_Volume',...End Fourier20
        'CSA','Perimeter', 'Volume',...End Biometry
        'conic_window','conic_R','conic_p','conic_z0',...
        'conic_rsquare','conic_rmse','conic_sse','conic_dfe',...
        'conic_adjrsquare',...End Conic Fit
        'sph_window','sph_R','sph_z0',...
        'sph_rsquare','sph_rmse','sph_sse',...
        'sph_dfe','sph_adjrsquare'}; %End Spherical Fit
    
    legacyTable_sagittal = legacyTable_xl(:, sagittal_colNames);
    
    if isempty(sagittal_expID) %The image is missing
        legacyTable_sagittal{:, 'ExperimentID'} = strrep(coronal_expID.ExperimentID(1),...
            'Cor','SagX');
    else
        legacyTable_sagittal(:, 'ExperimentID') = sagittal_expID(1, 'ExperimentID');
    end
end


tableContainer.LegacyCoronal = legacyTable_coronal;
tableContainer.LegacySagittal = legacyTable_sagittal;

function tableOut = generateExperimentTable(tableIn)
%This function generates experiment ID given subjectID, eye position, and
%date.
%INPUT
%tableIn - a table that has m rows and at least the following columns:
%SubjectID - string
%Eye - string
%ExperimentDate - string (M/D/Y)
%ImgType - string Sagittal/Coronal
%Img_originalLocation
%PMT
%OUTPUT
%experimentID - table experiment ID in table format with the columns:
%SubjectID
%ExperimentID
%
%OBSOLETE: tableOut - a table with two columns:
%ExperimentID - string
%ID Flag - logical: false if there ws an error generating the ID


%..........................................................................
%First, Generate Experiment ID
%..........................................................................
subjectID = tableIn{:, 'SubjectID'};
dateStr = tableIn{:,'ExperimentDate'};
eyePos = cellstr(char(tableIn{:, 'Eye'}));
imgType = cellstr(char(tableIn{:, 'ImgType'}));
imgNo = tableIn{:, 'ImgNo'};


%Change the format of the date object to: YYYYmmmDD Example: 2016Jun05
try
    dateStr_formatted = cellfun(@(iStr) datestr(datetime(iStr,'InputFormat','M/d/uuuu'),...
        'yyyymmmdd'), dateStr,'UniformOutput', false);
catch
    try
        dateStr_formatted = cellfun(@(iStr) datestr(datetime(iStr,'InputFormat','M/d/uuuu'),...
            'yyyymmmdd'), dateStr,'UniformOutput', false);
    catch
        dateStr_formatted = repmat({'UnknownDt'},[numel(subjectID),1]);
    end
end

experimentID = arrayfun(@(iRow) strjoin({[subjectID{iRow}, eyePos{iRow,:}],dateStr_formatted{iRow},...
    [imgType{iRow}(1:3),sprintf('%02d',imgNo(iRow))]},'_'), 1:numel(subjectID),'UniformOutput', false)';

%..........................................................................


expIDTable = table(force1D(experimentID,2), 'VariableNames', {'ExperimentID'});
expIDTable.Properties.VariableDescriptions = {'The Experiment ID is: SubjectID, Eye Position, Experiment Date, Image Type, and Image Number'};

tableOut = cat(2, tableIn(:,{'SubjectID'}), expIDTable,...
    tableIn(:,{'Eye','PMT','ExperimentDate','ImgType','Img_originalLocation'}));

%Additional Fields
tableOut{:,{'Img_filename'}} = {''};
tableOut{:,{'Img_serverDirectory'}} = {''};
tableOut{:,{'Img_localDirectory'}} = {''};
tableOut{:,{'AcquiredBy'}} = {'Ashik Mohamed'};
tableOut{:,{'OutlierFlag'}} = false(1,2); %outliers: eye, ant, pos
tableOut{:,{'ValidatedFlag'}} = false;
tableOut{:,{'Notes'}} = {''};

tableMeta = {'Img_filename', 'Filename of consolidated image in local directory: filename.jpg';...
    'Img_serverDirectory', 'Location of the image file on the cloud server';...
    'Img_localDirectory', 'Location of the image file on the local computer';...
    'AcquiredBy','Research team member that acquired the image.';...
    'OutlierFlag', 'Whether to exclude the lens surface';...
    'ValidatedFlag','Whether the image has been validated';...
    'Notes','Notes about the image'};

for j = 1:size(tableMeta, 1)
    tableOut = writeTableMetadata(tableOut, tableMeta{j,1}, tableMeta{j,2},'');
end

function tableContainer = removeDatabaseDuplicates(tableContainer)
%This function removes duplicate ID keys from Subject and Experiment tables

subjectTable = tableContainer.Subject; 
subjectNew = removeTableDuplicates(subjectTable, 'SubjectID');
%Delete empty subjectIDs
subjectNew = subjectNew(~cellfun(@isempty, subjectNew.SubjectID),:);

expTable = tableContainer.Experiment;
[expNew, oldIDs] = removeTableDuplicates(expTable, 'ExperimentID');


tableContainer.Subject = subjectNew;
tableContainer.Experiment = expNew; 


%Update the Legacy tables based on the new IDs
if ~isfield(tableContainer, 'LegacySagittal')
    return 
end

legacySag = tableContainer.LegacySagittal;
legacyCor = tableContainer.LegacyCoronal;

newIDs = expNew.ExperimentID;
isChanged = ~arrayfun(@(idx) strcmpi(newIDs{idx}, oldIDs{idx}), 1:numel(oldIDs));
changedIdx = find(isChanged);

for i = 1:numel(changedIdx)
    iOldID = oldIDs{changedIdx(i)};
    iNewID = newIDs{changedIdx(i)};
    
    isThisID_Sag = strcmpi(legacySag.ExperimentID,iOldID);
    isThisID_Cor = strcmpi(legacyCor.ExperimentID,iOldID);
    
    if any(isThisID_Sag)
        legacySag{isThisID_Sag,'ExperimentID'} = cellstr(iNewID);
    elseif any(isThisID_Cor)
        legacyCor{isThisID_Cor,'ExperimentID'} = cellstr(iNewID);
    else
        continue
    end

end



if isfield(tableContainer, 'Sagittal') 
    sagTable = tableContainer.Sagittal;
    
    if isempty(getMetadata(sagTable, 'CSA')) 
        %Add metadata to the Sagittal Table
        sagTable = writeTableMetadata(sagTable,'CSA',...
            'Cross Sectional Area: count of pixels within the lens ROI','mm^2');
        sagTable = writeTableMetadata(sagTable,'Perimeter',...
            'Count of pixels around the lens ROI','mm');
        sagTable = writeTableMetadata(sagTable,'Thickness',...
            'Maximum height of the lens','mm');
        sagTable = writeTableMetadata(sagTable,'Diameter',...
            'Maximum width of the lens','mm');
        
        tableContainer.Sagittal = sagTable;
    end
end




tableContainer.LegacySagittal = legacySag;
tableContainer.LegacyCoronal = legacyCor;

function [tableOut, oldID] = removeTableDuplicates(tableIn, keyField)
%This function removes duplicates of a key field and copies information
%from duplicate rows if it is not empty.  They keyField must be a string.
%KeyFields: {'SubjectID', 'ExperimentID','AnalysisID'}


nRows = size(tableIn,1);

%make sure the keyField is unique
[~, ia, ~] = unique(tableIn.(keyField));
isUnique = false(nRows,1);
isUnique(ia) = true;
duplicateRows = tableIn(~isUnique, :);
uniqueRows = tableIn(isUnique,:);

oldID = tableIn.(keyField);

%..........................................................................
%For Experiment Table
%..........................................................................
if any(strcmpi(tableIn.Properties.VariableNames, 'Img_originalLocation'))...
        && strcmpi(keyField,'ExperimentID')
    
    %these are from different images, and shouldn't be considerd
    %duplicates, I need to change the Experiment IDs
    
    dupExpIDs_trunk = cellfun(@(iID) iID(1:end-2),...
        duplicateRows.(keyField),'UniformOutput',false);
    [~,ia2,~] = unique(dupExpIDs_trunk);
    duplicateRows2 = duplicateRows(ia2,:); %make sure you have unique duplicates
    dupExpIDs_trunk2 = dupExpIDs_trunk(ia2);
    for j = 1:size(duplicateRows2, 1)
        %rename all of the Experiment IDs with new numbers
        %Find all rows with the experiment ID without last two digits
        
        
        jExpID_trunk = dupExpIDs_trunk2{j};
        jIsThisRow_all = ~cellfun(@isempty,...
            strfind(tableIn.(keyField), jExpID_trunk(1:end-2)));
        
        %find all rows with the same ID
        jRowsWithThisID = tableIn(jIsThisRow_all, :);
        
        %make sure it's a different file name
        [~, ja, ~] = unique(jRowsWithThisID.Img_originalLocation);
        jIsNotDuplicateRow = jRowsWithThisID(ja,:);
        
        %change the name of the id to a new ID
        jNewIDs = arrayfun(@(iNo) sprintf('%s%02d',jExpID_trunk,iNo),...
            1:size(jIsNotDuplicateRow,1),'UniformOutput', false);
        
        tableIn.(keyField)(jIsThisRow_all) = jNewIDs';
    end
    
    
    %now you can remove all the duplicates from tableIn
    [~, ia3, ~] = unique(tableIn.(keyField));
    isUnique = false(nRows,1);
    isUnique(ia3) = true;
    tableOut = tableIn(isUnique,:);
    oldID = oldID(ia3);
    return
end


%..........................................................................
%For Everything Else
%..........................................................................
oldID = []; %ID won't change
uniqueCells = table2cell(uniqueRows);
duplicateCells = table2cell(duplicateRows);

for i = 1:size(duplicateRows,1)
    %Copy the data if it's missing from the uniqueRows
    iRow_duplicate = duplicateRows(i,:);
    isThisRow = strcmpi(uniqueRows.(keyField),...
        iRow_duplicate.(keyField));
    
    iRow_mainCell = uniqueCells(isThisRow,:);
    iRow_dupCell = duplicateCells(i,:);
    
    isPresent_main =  cellfun(@(iInput) strcmpi(class(iInput),'categorical')...
        || (~isempty(iInput) && ~any(isnan(iInput))) ,...
        iRow_mainCell);
    
    if all(isPresent_main)
        continue
    end
    
    isPresent_duplicate = cellfun(@(iInput) strcmpi(class(iInput),'categorical')...
        || (~isempty(iInput) && ~any(isnan(iInput))) ,...
        iRow_dupCell);
    
    isOverwrite = ~isPresent_main & isPresent_duplicate;
    
    if any(isOverwrite)
        iRow_mainCell(isOverwrite) = iRow_dupCell(isOverwrite);
        uniqueRows(isThisRow,:) = cell2table(iRow_mainCell);
    else
        continue
    end
end

tableOut = uniqueRows;


%..........................................................................
%table metadata functions
%..........................................................................
function tableOut = copyTableMetadata(tableOld, tableNew)
%This function copies the descriptions and variable units from table Old to
%tableNew

oldDescrips = tableOld.Properties.VariableDescriptions;
oldUnits = tableOld.Properties.VariableUnits;
tableOut = tableNew;

isOldUnitEmpty = isempty(oldUnits);
isOldDescripsEmpty = isempty(oldDescrips);

%If there is no information in the old table, return
if isempty(oldDescrips) && isempty(oldUnits)
    return
end

newColNames = tableNew.Properties.VariableNames;
oldColNames = tableOld.Properties.VariableNames;

%if it is completely empty, initialize it before the loop
newDescrips = tableNew.Properties.VariableDescriptions;
if isempty(newDescrips)
    newDescrips = repmat({''},[1, numel(newColNames)]);
end

newUnits = tableNew.Properties.VariableUnits;
if isempty(newUnits)
    newUnits = repmat({''},[1, numel(newColNames)]);
end


for i = 1:numel(newColNames)
    
    thisName = newColNames{i};
    
    %Check to see if the column in the new table is in the old table
    if isfield(table2struct(tableOld), thisName)
        colIdx = find(cellfun(@(iOldName) strcmpi(iOldName, thisName),...
            oldColNames));
        
        %Copy the descriptions if the old one isn't empty
        if ~isOldDescripsEmpty && ~isempty(oldDescrips{colIdx})
            this_oldDescrip = oldDescrips{colIdx};
            newDescrips{i} = this_oldDescrip;
        end
        
        %Copy the units if the old one isn't empty
        if ~isOldUnitEmpty && ~isempty(oldUnits{colIdx})
            this_oldUnit = oldUnits{colIdx};
            newUnits{i} = this_oldUnit;
        end
    end
    
end

tableOut.Properties.VariableDescriptions = newDescrips;
tableOut.Properties.VariableUnits = newUnits;

function tableOut = writeTableMetadata(tableIn, colName, description, units)
%This function writes metadata to tables
%INPUT
    %tableIn
    %colName - string
    %description - string
    %units - string
%OUTPUT
    %tableOut - same table with metadata


if nargin < 3
    description = '';
end

if nargin < 4
    units = '';
end


colNames = categorical(tableIn.Properties.VariableNames);
tableOut = tableIn;

colIdx = find(colNames == colName,1);
if isempty(colIdx)
    return
end

descr_new = tableIn.Properties.VariableDescriptions;
units_new = tableIn.Properties.VariableUnits;
nRows = size(tableIn, 1);

if isempty(descr_new)
    descr_new = repmat({''},[1,nRows]);
end

if isempty(units_new)
    units_new = repmat({''},[1,nRows]);
end

descr_new{colIdx} = description;
units_new{colIdx} = units;

tableOut.Properties.VariableDescriptions = descr_new;
tableOut.Properties.VariableUnits = units_new;

function [colUnits, colDescription] = getMetadata(tableIn, colName)
%This function gets the column description and units of a given column name

%Initialize Output
colUnits = [];
colDescription = [];

colNames = categorical(tableIn.Properties.VariableNames); 

colIdx = find(colNames == colName,1);
if isempty(colIdx)
    %Column not Found
    return
end

colUnits = tableIn.Properties.VariableUnits{colIdx};
colDescription = tableIn.Properties.VariableDescriptions{colIdx};

%..........................................................................
%table filter functions
%..........................................................................

function [tableOut, logicalIdx] = filterTable(tableIn, query, logCombFcn)
%This function filters the tabel based on column names given the query for
%each column
%INPUT
%tableIn
%query - cell containing the column name and the query for that column
%{colName1, queryStr1; colName2, queryStr2;...colNameN, queryStrN};
%colNames - names of the columns being queried
%queryStr -
%cell array with queries, each row of the cell array refers to a
%different colName query, nRows = numel(colName)
%the name of the column must be used in the correct MATLAB syntex
%in the query string:
%Example: 'strcmpi(SubjectID, 'RIE023') || strcmpi(SubjectID,'RIE001')'
%         'Age <= 5 && Age >=3'
%combFcn - an anonymous function that determines how the queries are
%combined, OR, AND, a mixture.
%Example: @(logArray) logArray(:,1) & (logArray(:,2) | logArray(:,3))

if nargin < 3
    logCombFcn = @(logArray) any(logArray,2); %uses OR
end

nRows = size(tableIn,1);
nQueries = size(query,1);
isRowValid = false(nRows, nQueries);

for i = 1:nQueries
    
    iQuery = query{i,2};
    iColName = query{i,1};
    thisCol = tableIn.(iColName);
    
    
    if ~isempty(strfind(iColName,'Date'))
        %this is a string that should be a date
        thisCol = datetime(thisCol);
    end
    
    %if its categorical and the query is using strcmp, make it a string
    if strcmpi(class(thisCol),'categorical') && ~isempty(strfind(iQuery,'strcmp'))
        thisCol = cellstr(thisCol);
    end
    
    %make sure the column is described as the column name
    
    %Put a variable into memory with the same name as the column name
    eval(sprintf('%s = thisCol;',iColName));
    
    %Run the query
    eval(sprintf('iIsRowValid = %s;',iQuery));
    
    isRowValid(:,i) = iIsRowValid;
end

if nQueries > 1
    logicalIdx = logCombFcn(isRowValid);
else
    logicalIdx = isRowValid;
end
tableOut = tableIn(logicalIdx,:);

function qStrOut = buildQueryString(colName, qString, logOp)
%this function builds queries for string classes using strcmpi
%INPUT:
%colName - character array containing the column name
%qStr - cell array containing the name of each string you are looking
%for. Example: qString = {'Subject1', 'Subject2'};
%logOp - character array containing a logical operator: logOp = '||'
%OUPTUT
%qStrOut - character array containing the query string
%Example: 'strcmpi(SubjectID, 'Subject1') || strcmpi(SubjectID,'Subject2')'

if nargin < 3
    logOp = '|'; %default the logical operator to OR
end

qString = cellstr(qString);

%build the string elements
qStr_split = cellfun(@(iStr) sprintf('strcmpi(%s,''%s'')', colName, iStr),...
    qString, 'UniformOutput', false);

%concatenate it with the logical operator
qStrOut = strjoin(qStr_split, logOp);


%--------------------------------------------------------------------------
%Extract Data from Excel Files
%--------------------------------------------------------------------------

function [hasXLFile, iXL_fullpath, nFiles] = getXLfilename(parentDir, isReturnAll)
%This function is used to check a directory for an Excel File (*.xls then
%*.xlsx) If more than one file is found, it returns the first

if nargin < 2
    isReturnAll = false;
end

%Initialize variables
hasXLFile = true;
iXL_fullpath = [];
nFiles = [];

%first check for an excel file
iXLfile = dir(fullfile(parentDir, '*.xls'));

if isempty(iXLfile) %check for an advanced excel file
    iXLfile = dir(fullfile(parentDir, '*.xlsx'));
    
    if isempty(iXLfile)
        hasXLFile = false;
        return
    end
end

%print a warning letting the user know that multiple files were
%found and you are using the first
nFiles = numel(iXLfile);
if nFiles > 1 && ~isReturnAll
    delimStr = '================================================';
    fprintf(1,['\n\n%s\nWARNING: Mulitple Excel files found\n%s\n\n',...
        'Folder: %s\n\nFiles:\n\n%s']', delimStr, delimStr, ancestorDir(parentDir,0));
    
    cellfun(@disp,{iXLfile.name}');
    fprintf(1,'\nThe first file will be used: \n%s\n\n%s\n', iXLfile(1).name, delimStr);
end

%what if there are multiple XL files??? just take the first...
if nFiles > 1 && isReturnAll
    iXL_fullpath = arrayfun(@(i) fullfile(parentDir, iXLfile(i).name), 1:nFiles,...
        'UniformOutput', false);
else
    iXL_fullpath = fullfile(parentDir, iXLfile(1).name);
end

function [tableOut, xlFlags] = readXLfile(iXL_fullpath, xl_sheetNames)
%This function searches the parent directory for an excel file and returns
%the data if it is present. If there is no XL file or the sheet cannot be
%read, it returns a blank/NaN legacy table
%INPUT:
%iXL_fullpath - string: 'rootDir\parentDir\XLfile.xls'
%xl_sheetNames - cell: {'sheetName1', 'sheetName2'...}
%OUTPUT:
%tableOut - legacy table containing ALL columns from the Excel sheet.
%Table is initialized if the sheet is not found or cannot be read
%xlFlags - logical array (same number of elements as xl_sheetNames).
%Value = true if sheet could not be read


%initialize function variables
tableOut = [];



if nargin < 2
    xl_sheetNames = {'Eye Data', 'Fourier Coefficents',... %NOTE: Coefficents is mispelled
        'Original', 'Centered And Aligned', 'Biometry','Curvature'};
end



%..................................................................
%READ XL SHEET
%..................................................................


nSheets = numel(xl_sheetNames);
xlFlags = false(1,nSheets);

xl_output = cell(1,nSheets);
tableOut = table([]); %Table containing information from Excel file

for j = 1:nSheets
    jSheetName = xl_sheetNames{j};
    
    %read each sheet
    [iTable, errFlag, errMsg] = readLegacyXLsheet(iXL_fullpath, jSheetName);
    
    
    xl_output{j} = iTable; %individualized output
    xlFlags(j) = errFlag;
    
    if isempty(tableOut)
        tableOut = iTable;
    elseif iscell(iTable)
        tableOut = cat(2, tableOut, iTable{:});
    else
        tableOut = cat(2, tableOut, iTable);
    end
    
end

%..........................................................................
%Read XL File: Sub-Functions
%..........................................................................

function [tableOut, errFlag, errMsg] = readLegacyXLsheet(iXL_fullpath, jSheetName)
%This function reads a shadowgraph Excel sheet given the full filepath and
%sheetname (as it appears in the Excel sheet).  If it cannot read the
%sheet, it initializes the row to NaNs.
%INPUT:
%iXL_fullpath: full file location of excel file: 'rootDir\parentFolder\filename.xls'
%jSheetName: name of the sheet you'd like to read. Each "switch case"
%reflects a sheet in the Excel file verbatim
%OUTPUT:
%tableOut - depending on the sheet, it will return a cell {table1,
%table2} or a table
%errFlag - the error flag becomes true if the catch case is activated
%errMsg - error structure of catch statement


%Initialize outputs
tableOut = [];
errFlag = false;
errMsg = [];


switch jSheetName
    case 'Eye Data'
        try
            [~,~,jData] = xlsread(iXL_fullpath, jSheetName,'A1:B6');
            %Choose structure or table layout
            eyeData = table; %struct([]);
            
            if ischar(jData{2,2}) %Age, verify datatype
                jData{2,2} = str2double(jData{2,2});
            end
            
            if ischar(jData{5,2}) %PMT, verify datatype
                jData{5,2} = str2double(jData{5,2});
            end
            
            
            eyeData.SubjectID = jData(1,2); %make sure this is a cell
            eyeData.Age = jData{2,2}; %years
            eyeData.Species = categorical(jData(3,2)); %Human/Cyno cell string
            eyeData.Eye = categorical(jData(4,2)); %OS/OD
            eyeData.PMT = jData{5,2}; %Post Mortem Time (hours)
            eyeData.ExperimentDate = jData(6,2); %cell string
            
            tableOut = addLegacyMetadata(eyeData, 'subject');
        catch errMsg
            errFlag = true;
            tableOut = initializeXL_Legacy('subject');
        end
        
        
        
    case 'Fourier Coefficents' %coefficients is MISPELLED in the Excel sheet
        
        try
            [~,~,jData1] = xlsread(iXL_fullpath, jSheetName,'A1:M2');
            [~,~,jData2] = xlsread(iXL_fullpath, jSheetName,'A4:W5');
            
            fourier10Labels = jData1(1,:); fourier10Data = jData1(2,:);
            fourier20Labels = jData2(1,:); fourier20Data = jData2(2,:);
            
            %Choose structure or table layout
            fourierData = table; %struct([]);
            
            fourierData.Fourier10Coeff = force1D(cell2mat(fourier10Data(1:end-2)),1);
            fourierData.Fourier20Coeff = force1D(cell2mat(fourier20Data(1:end-2)),1);
            fourierData.Fourier10RMSE = [fourier10Data{end-1:end}];
            fourierData.Fourier20RMSE = [fourier20Data{end-1:end}];
            
            tableOut = addLegacyMetadata(fourierData, 'fourier');
        catch errMsg
            errFlag = true;
            tableOut = initializeXL_Legacy('fourier');
        end
        
        
    case 'Original'
        
        try
            [~,~,jData] = xlsread(iXL_fullpath, jSheetName);
            
            originalPts = cell2mat(jData);
            origPtsData = table({originalPts}, 'VariableNames',{'OriginalPoints'});
            tableOut = addLegacyMetadata(origPtsData, 'original_points');
        catch errMsg
            errFlag = true;
            tableOut = initializeXL_Legacy('original_points');
        end
        
        
    case 'Centered And Aligned'
        
        try
            [~,~,jData] = xlsread(iXL_fullpath, jSheetName);
            centeredPts = cell2mat(jData);
            corrPtsData = table({centeredPts}, 'VariableNames',{'CorrectedPoints'});
            tableOut = addLegacyMetadata(corrPtsData, 'corrected_points');
        catch errMsg
            errFlag = true;
            tableOut = initializeXL_Legacy('corrected_points');
        end
        
    case 'Biometry'
        
        try
            biometryData = table; %struct([]);
            
            %10 Coefficient Fit
            [~,~,jData1] = xlsread(iXL_fullpath, jSheetName,'A2:G3');
            
            biometryData.Fourier10_AnteriorThickness = jData1{2,1}; %mm
            biometryData.Fourier10_PosteriorThickness = jData1{2,2}; %mm
            biometryData.Fouier10_LensThickness = jData1{2,3}; %mm (Ant + Pos Thickness)
            biometryData.Fourier10_LensDiameter = jData1{2,4}; %mm
            biometryData.Fourier10_CSA = jData1{2,5}; %mm^2
            biometryData.Fouier10_SurfaceArea = jData1{2,6}; %mm^2
            biometryData.Fourier10_Volume = jData1{2,7}; %mm^3
            
            %20 CoefficientFit
            [~,~,jData2] = xlsread(iXL_fullpath, jSheetName,'A6:G7');
            
            biometryData.Fourier20_AnteriorThickness = jData2{2,1}; %mm
            biometryData.Fourier20_PosteriorThickness = jData2{2,2}; %mm
            biometryData.Fourier20_LensThickness = jData2{2,3}; %mm
            biometryData.Fourier20_LensDiameter = jData2{2,4}; %mm
            biometryData.Fourier20_CSA = jData2{2,5}; %mm^2
            biometryData.Fouier20_SurfaceArea = jData2{2,6}; %mm^2
            biometryData.Fourier20_Volume = jData2{2,7}; %mm^3
            
            
            %Pixel Count on Raw Image
            [~,~,jData3] = xlsread(iXL_fullpath, jSheetName,'A10:B12');
            
            biometryData.CSA = jData3{2,1}; %mm^2
            biometryData.Perimeter = jData3{2,2}; %mm
            
            
            %Volume from Centered and Aligned Data
            [~,~,jData4] = xlsread(iXL_fullpath, jSheetName,'A14:A15');
            
            biometryData.Volume = jData4{2,1}; %mm^3
            
            %Coronal Data
            [~,~,jData5] = xlsread(iXL_fullpath, jSheetName,'A18:F19');
            
            coronalData = table;
            
            coronalData.Coronal_AvgR = jData5{2,1}; %mm
            coronalData.Coronal_MaxR = jData5{2,2}; %mm
            coronalData.Coronal_MinR = jData5{2,3}; %mm
            coronalData.Coronal_StdR = jData5{2,4}; %mm
            coronalData.Coronal_CSA = jData5{2,5}; %mm^2
            coronalData.Coronal_Circumference = jData5{2,6}; %mm
            
            
            tableOut{1} = addLegacyMetadata(biometryData, 'biometry');
            tableOut{2} = addLegacyMetadata(coronalData, 'coronal');
            
        catch errMsg
            errFlag = true;
            tableOut{1} = initializeXL_Legacy('biometry');
            tableOut{2} = initializeXL_Legacy('coronal');
            
        end
        
        
    case 'Curvature'
        try
            
            %Curvature: 4mm Conic Fit
            [~,~,jData1] = xlsread(iXL_fullpath, jSheetName,'A2:I4');
            
            conicCurveData = table;
            conicCurveData.conic_window = 4; %4mm window
            conicCurveData.conic_R = [jData1{2,2},jData1{3,2}];
            conicCurveData.conic_p = [jData1{2,3},jData1{3,3}];
            conicCurveData.conic_z0 = [jData1{2,4},jData1{3,4}];
            conicCurveData.conic_rsquare = [jData1{2,5},jData1{3,5}];
            conicCurveData.conic_rmse = [jData1{2,6},jData1{3,6}];
            conicCurveData.conic_sse = [jData1{2,7},jData1{3,7}];
            conicCurveData.conic_dfe = [jData1{2,8},jData1{3,8}];
            conicCurveData.conic_adjrsquare = [jData1{2,8},jData1{3,8}];
            
            
            %Curvature: 3mm Spherical Fit
            [~,~,jData2] = xlsread(iXL_fullpath, jSheetName,'A7:H9');
            
            sphCurveData = table;
            sphCurveData.sph_window = 3; %3mm window
            sphCurveData.sph_R = [jData2{2,2}, jData2{3,2}];
            sphCurveData.sph_z0 = [jData2{2,3}, jData2{3,3}];
            sphCurveData.sph_rsquare = [jData2{2,4}, jData2{3,4}];
            sphCurveData.sph_rmse = [jData2{2,5}, jData2{3,5}];
            sphCurveData.sph_sse = [jData2{2,6}, jData2{3,6}];
            sphCurveData.sph_dfe = [jData2{2,7}, jData2{3,7}];
            sphCurveData.sph_adjrsquare = [jData2{2,8}, jData2{3,8}];
            
            
            tableOut{1} = addLegacyMetadata(conicCurveData, 'conCurve');
            tableOut{2} = addLegacyMetadata(sphCurveData, 'sphCurve');
            
        catch errMsg
            errFlag = true;
            tableOut{1} = initializeXL_Legacy('conCurve');
            tableOut{2} = initializeXL_Legacy('sphCurve');
            
        end
end

function tableOut = initializeXL_Legacy(tableName)
%This function initializes sub-table in legacy megatable for Excel Data

switch lower(tableName)
    
    case 'subject'
        %SHEET 1: Eye Data
        eyeData = table({''},nan,categorical({''}),categorical({''}),nan,{''},...
            'VariableNames',{'SubjectID','Age', 'Species', 'Eye','PMT','ExperimentDate'});
        
        tableOut = addLegacyMetadata(eyeData,'subject');
        
    case 'fourier'
        %SHEET 2: Fourier 'Coefficents'
        fourierData =  table(nan(1, 11),nan(1,21),nan(1,2),nan(1,2),...
            'VariableNames',...
            {'Fourier10Coeff', 'Fourier20Coeff', 'Fourier10RMSE','Fourier20RMSE'});
        tableOut = addLegacyMetadata(fourierData, 'fourier');
        
    case 'original_points'
        %SHEET 3: Original
        origPtsData = table({[]}, 'VariableNames',{'OriginalPoints'});
        tableOut = addLegacyMetadata(origPtsData, 'original_points');
        
    case 'corrected_points'
        %SHEET 4: Centered And Aligned
        corrPtsData = table({[]}, 'VariableNames',{'CorrectedPoints'});
        tableOut = addLegacyMetadata(corrPtsData, 'corrected_points');
        
    case 'biometry'
        %SHEET 5: Biometry & Coronal
        
        biometryData = table(nan,nan,nan,nan,nan,nan,nan,...Fourier10
            nan,nan,nan,nan,nan,nan,nan,...Fourier20
            nan, nan,nan,...CSA, Perimeter, Volume
            'VariableNames',{'Fourier10_AnteriorThickness', 'Fourier10_PosteriorThickness',...
            'Fouier10_LensThickness', 'Fourier10_LensDiameter', 'Fourier10_CSA',...
            'Fouier10_SurfaceArea','Fourier10_Volume',...End Fourier10
            'Fourier20_AnteriorThickness', 'Fourier20_PosteriorThickness',...
            'Fourier20_LensThickness','Fourier20_LensDiameter','Fourier20_CSA',...
            'Fouier20_SurfaceArea','Fourier20_Volume',...End Fourier20
            'CSA','Perimeter', 'Volume'});
        
        tableOut = addLegacyMetadata(biometryData, 'biometry');
    case 'coronal'
        coronalData = table(nan, nan, nan, nan, nan, nan,...
            'VariableNames',{'Coronal_AvgR', 'Coronal_MaxR','Coronal_MinR'...
            'Coronal_StdR','Coronal_CSA','Coronal_Circumference'});
        
        tableOut = addLegacyMetadata(coronalData, 'coronal');
        
    case 'concurve'
        %SHEET 6: Curvature (Conic & Spherical)
        conicCurveData = table(nan, nan(1,2),nan(1,2),nan(1,2),...
            nan(1,2),nan(1,2), nan(1,2), nan(1,2),nan(1,2),...
            'VariableNames',{'conic_window','conic_R','conic_p','conic_z0',...
            'conic_rsquare','conic_rmse','conic_sse','conic_dfe','conic_adjrsquare'});
        
        tableOut = addLegacyMetadata(conicCurveData, 'conCurve');
    case 'sphcurve'
        
        sphCurveData = table(nan, nan(1,2),nan(1,2),...
            nan(1,2),nan(1,2), nan(1,2), nan(1,2),nan(1,2),...
            'VariableNames',{'sph_window','sph_R','sph_z0',...
            'sph_rsquare','sph_rmse','sph_sse','sph_dfe','sph_adjrsquare'});
        
        tableOut = addLegacyMetadata(sphCurveData, 'sphCurve');
    case 'original_points'
        
        origPtsData = table({}, 'VariableNames',{'OriginalPoints'});
        tableOut = addLegacyMetadata(origPtsData, 'original_points');
    case 'corrected_points'
        
        corrPtsData = table({}, 'VariableNames',{'CorrectedPoints'});
        tableOut = addLegacyMetadata(corrPtsData, 'corrected_points');
end

function tableOut = addLegacyMetadata(tableIn, tableName)
%This is a very specific sub-function designed to work with data of a
%particular order.  It adds descriptions and units to tables created from
%reading excel data

switch lower(tableName)
    case 'subject'
        tableIn.Properties.VariableUnits = {'', 'years','',...
            '', 'hours', 'M/D/Y'};
        tableIn.Properties.VariableDescriptions = {'Eye Source, Year, and ID number',...
            'Age of the eye in years', 'Species of the eye',...
            'Left (OS) or Right (OD) eye',...
            'Post Mortem Time before imaging', 'Date of the Experiment'};
    case 'fourier'
        %no need to add units, just variable descriptions
        tableIn.Properties.VariableDescriptions =...
            {'10 Fourier fit coefficients with DC offset a0 (11 values)',...
            '20 Fourier fit coefficients with DC offset a0 (21 values)',...
            'Fourier 10 fit Root Mean Square Error: [RMSE10, RMSE10 - DOF]',...
            'Fourier 20 fit Root Mean Square Error: [RMSE10, RMSE10 - DOF]'};
    case 'biometry'
        tableIn.Properties.VariableUnits =...
            {'mm','mm','mm','mm','mm^2','mm^2','mm^3',...Fourier10
            'mm','mm','mm','mm','mm^2','mm^2','mm^3',...Fourier20
            'mm^2','mm','mm^3'}; %Pixel Count & Volume
        
        tableIn.Properties.VariableDescriptions =...
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
    case 'coronal'
        tableIn.Properties.VariableUnits =...
            {'mm','mm','mm','mm','mm^2','mm'}; %Coronal Data
        tableIn.Properties.VariableDescriptions =...
            {'Coronal: Average Radius',...
            'Coronal: Max Radius',...
            'Coronal: Min Radius',...
            'Coronal: Std Radius',...
            'Coronal: Cross Sectional Area',...
            'Coronal: Circumference'};
    case 'concurve'
        tableIn.Properties.VariableUnits =...
            {'mm', 'mm', 'mm','mm','','mm','mm','',''};
        
        tableIn.Properties.VariableDescriptions =...
            {'Window size used to fit segmentation points',...
            'Radius of curvature: [anterior, posterior]',...
            'Shape factor: p',...
            'Offset: z0',...%note: there is no x0  because it should be 0
            'Goodness of Fit: R Square',...
            'Goodness of Fit: Root Mean Square Error (RMSE)',...
            'Goodness of Fit: Sum Square of the Error (SSE)',...
            'Goodness of Fit: Degrees of freedom (DFE)',...
            'Goodness of Fit: Adjusted R Square'};
        
    case 'sphcurve'
        
        tableIn.Properties.VariableUnits =...
            {'mm', 'mm', 'mm','','mm','mm','',''};
        
        tableIn.Properties.VariableDescriptions=...
            {'Window size used to fit segmentation points',...
            'Radius of curvature: [anterior, posterior]',...
            'Offset: z0',...%note: there is no x0  because it should be 0
            'Goodness of Fit: R Square',...
            'Goodness of Fit: Root Mean Square Error (RMSE)',...
            'Goodness of Fit: Sum Square of the Error (SSE)',...
            'Goodness of Fit: Degrees of freedom (DFE)',...
            'Goodness of Fit: Adjusted R Square'};
    case 'original_points'
        tableIn.Properties.VariableUnits ={'mm'};
        tableIn.Properties.VariableDescriptions =...
            {'Points from the original segmentation in mm with images posterior up.'};
    case 'corrected_points'
        tableIn.Properties.VariableUnits ={'mm'};
        tableIn.Properties.VariableDescriptions =...
            {'These are the points from the untilted and centered segmentation with images posterior up.'};
end

tableOut = tableIn;

%--------------------------------------------------------------------------

%..........................................................................
%SWITCHING Tabs
%..........................................................................
function switchTab_Callback(hObject, eventdata, handles)
%This function is called on a tab change

hFig = ancestor(handles.axes_img, 'figure');
hAllTabs = [handles.tab_image, handles.tab_graph, handles.tab_table,...
    handles.tab_settings];

%the one with an odd color is the last tab
allColors = get(hAllTabs,'BackgroundColor'); allColors = cat(1,allColors{:});
hLastTab = hAllTabs(~all(allColors == repmat(mode(allColors), [size(allColors,1),1]),2));

if strcmpi(get(hLastTab,'String'), 'settings')
    updateGraphPanel(handles);    
end

%if it was just settings, check for the settings flag

%Turn off all tabs
set(hAllTabs, 'Value', 0,'BackgroundColor',[0.94, 0.94, 0.94],...
    'ForegroundColor','k','FontWeight','normal');
set(hObject, 'Value', 1, 'BackgroundColor', [0.071, 0.212, 0.141],...
    'ForegroundColor','w', 'FontWeight','bold');

%Turn off all panels
hAllPanels = [handles.panel_single, handles.panel_table, handles.panel_settings,...
    handles.panel_axesImg, handles.panel_axesGraph];
set(hAllPanels, 'Visible','off');

switch lower(get(hObject,'String'))
    case 'image'
        %Turn on the panel
        set([handles.panel_single, handles.panel_axesImg], 'Visible','on');
        
    case 'graph'
        %Turn on the panel
        set([handles.panel_single, handles.panel_axesGraph], 'Visible','on');
        
        settingsFlag = getappdata(hFig, 'settingsFlag');
        
        if ~isempty(settingsFlag) && settingsFlag
            %re-run plots
            updateGraphPanel(handles)
        end
        
    case 'table'
        
        if strcmpi(get(hLastTab,'String'), 'settings')
            updateDataTableDisplay(handles);  %if any filters are changed
        end
        
        %Turn on the panel
        set(handles.panel_table, 'Visible','on');
        
    case 'settings'
        %Turn on the panel
        set(handles.panel_settings, 'Visible','on');
    
    otherwise
        %should never happen
end

function updateSubjectSummaryText(handles)
%This function updates the text_subjectCount control with information about
%database

hFig = ancestor(handles.axes_img,'figure');
tableContainer = getappdata(hFig, 'tableContainer');

if isempty(tableContainer) || ~isfield(tableContainer, 'Experiment')
    nLoaded = 0;
    nDisplayed = 0;
    nExcluded = 0;
    nAnalyzed = 0;
    totalImgs = 0;
else
    subTable = tableContainer.Subject;
    expTable = tableContainer.Experiment;
    
    nLoaded = size(subTable, 1);
    totalImgs = size(expTable,1);
    
    %Get number of subjects displayed
    subTable_lbox = getappdata(handles.listbox_subject,'table');
    if isempty(subTable_lbox)
        nDisplayed = 0;
    else
        nDisplayed = size(subTable_lbox,1);
    end
    
    %Get nExcluded
    outlierFlag = expTable.OutlierFlag;
    nExcluded = max([0,numel(find(all(outlierFlag,2)))]);
    
    %Get nAnalyzed 
    if ~isfield(tableContainer, 'Sagittal')
        nAnalyzed = 0; 
    else
        nAnalyzed = size(tableContainer.Sagittal,1);
    end
    
end

textboxStr = sprintf(['Subjects Visible: %d of %d\n',...
    'Images Analyzed: %d of %d\n',...
    'Images Excluded: %d of %d\n'],...
    nDisplayed, nLoaded,...
    nAnalyzed, totalImgs,...
    nExcluded, totalImgs);

set(handles.text_dbSummary,'String', textboxStr);


function b_send2xl_Callback(hObject, eventdata, handles)
% hObject    handle to b_send2xl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hTable = handles.table_data;
tabledata = get(hTable, 'Data');

if isempty(tabledata)
    errordlg('No Data Found')
end

colNames = get(hTable, 'ColumnName');

data4XL = cat(1, colNames', tabledata);

mat2xl(data4XL);

%==========================================================================
%SETTINGS PANEL: 
%==========================================================================

function b_localDir_Callback(hObject, eventdata, handles)
% hObject    handle to b_localDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Change the local directory of the images       
hFig = ancestor(hObject, 'figure');


tableContainer = getappdata(hFig, 'tableContainer');

if isempty(tableContainer) || ~isfield(tableContainer, 'Experiment')
    errordlg('No subjects have been found.  Please load data first',...
        'No Database Loaded');
    return
end

prevDir = getappdata(hFig, 'prevDir');

localDir = uigetdir(prevDir);

if ~localDir
    return
end



% isOk = setTableValue(handles, 'Experiment',...
%     'Img_localDirectory', cellstr(localDir),[]);
%keyValue=[] means apply to all images


setappdata(hFig,'prevDir', localDir);
set(handles.text_localDir, 'String', localDir,'ToolTipString', localDir);

set(hObject, 'ToolTipString', localDir);
setappdata(hFig, 'localDir', localDir);

updateLocalDir(handles, localDir);

function updateLocalDir(handles, localDir)
%This function updates the local directory location 

hFig = ancestor(handles.axes_img,'figure');

%This function updates the local directory of all images
tableContainer = getappdata(hFig, 'tableContainer');

if isempty(tableContainer) || ~isfield(tableContainer,'Experiment')
    return
end

if nargin < 2 || isempty(localDir)
    localDir = get(handles.text_localDir, 'String');
end

if isempty(localDir) || ~exist(localDir, 'dir')
    errordlg('Please set the local directory and try again.',...
        'Directory Not Set');
    return
end

expTable = tableContainer.Experiment;

imgFname = expTable{:,'Img_filename'};


%Create location
fileLoc = cellfun(@(iFname) fullfile(localDir, iFname), imgFname,...
    'UniformOutput', false);

%make sure it's there
isFileHere = cellfun(@(iLoc) logical(exist(iLoc, 'file')), fileLoc);

if ~any(isFileHere)
    errordlg('Image files are not found in the local directory.',...
        'Files Not Found');
    return
end

%update structure
expTable{isFileHere, 'Img_localDirectory'} = cellstr(localDir);

tableContainer.Experiment = expTable;

setappdata(hFig, 'tableContainer', tableContainer);

helpdlg('Local directory update successful.  Please remember to save changes.',...
    'Update Complete');


%..........................................................................

function popup_xField_Callback(hObject, eventdata, handles)
% hObject    handle to popup_xField (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_xField contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_xField

%Find the seletected value field and add the data to the popup appdata


thisVal = get(hObject, 'Value');

if strcmpi(get(hObject, 'Enable'),'off') || thisVal(1) == 0
    setappdata(hObject, 'data',[]); %data structure
    return
end

popupStr = get(hObject, 'String');
thisPopupStr = popupStr(thisVal); % can only be one

%convert the selected fields to data structure
[parentTableName, fieldName] = parseFieldString(thisPopupStr);
dataStructure = analyzeFieldString(handles, parentTableName, fieldName);

setappdata(hObject, 'data', dataStructure);
setappdata(ancestor(hObject,'figure'),'settingsFlag',true);

function b_addYField_Callback(hObject, eventdata, handles)
% hObject    handle to b_addYField (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hPopup = handles.popup_tableList1;
hListbox = handles.listbox_fieldList1;
hEdit = handles.edit_yField;
currentStr = get(hEdit, 'Str');


%Get selected field name
thisVal_left = get(hListbox, 'Value');
if strcmpi(get(hListbox, 'Enable'),'off') || thisVal_left(1) == 0
    return
end

listboxStr = cellstr(get(hListbox, 'String'));
thisField = listboxStr(thisVal_left); % can only be one

parentTableName = getappdata(hPopup,'parentTableName'); %char array

if isempty(parentTableName)
   parentTableName = updateFieldNames(hPopup);
end

%convert the selected fields to data structure
newStr = cellfun(@(iStr) sprintf('<<%s::%s>>', parentTableName, iStr),...
    thisField,'UniformOutput', false);

strOut = strjoin(cat(2,{currentStr}, force1D(newStr,1)),' ');

set(hEdit, 'String', strOut);

setappdata(ancestor(hObject,'figure'),'settingsFlag',true);

function setSettingsFlag(hObject, eventdata, handles)

setappdata(ancestor(hObject,'figure'),'settingsFlag',true);

%..........................................................................


function updateDataTableDisplay(handles, dataTable)
%This function updates the main table data display

if nargin < 2
    dataTable = readTableSettings(handles);
end

if isempty(dataTable)
    return
end

colNames = dataTable.Properties.VariableNames;

if ~any(ismember(colNames, 'ExperimentID'))
    
    %check valid SubjectIDs
    subTable = getappdata(handles.listbox_subject, 'table');
    validIDs = subTable.SubjectID;
    allIDs = dataTable.SubjectID;
else
    validIDs = getappdata(handles.panel_filter,'ExperimentID');
    allIDs = dataTable.ExperimentID; 
end

isValid = cellfun(@(iID) any(strcmpi(validIDs, iID)), allIDs);

if ~any(isValid)
    return
end

dataTable_valid = dataTable(isValid,:);

%check to see if anything is filtered

%update table
set(handles.table_data, 'Data', table2cell(dataTable_valid),...
    'ColumnName', colNames); 
setappdata(handles.table_data, 'table', dataTable);

function initializeTableSettings(handles)
%This function initializes the table output if there is Sagittal Data

hListbox_right = handles.listbox_tableColumns;
tableContainer = getappdata(ancestor(handles.axes_img, 'figure'), 'tableContainer');

if numel(get(hListbox_right, 'String')) > 0 || isempty(tableContainer)
    %don't overwrite existing fields
    return
end


%set right string
str = {'Subject::SubjectID', 'Experiment::ExperimentID',...   
    'Subject::Age', 'Subject::Gender'}';


if isfield(tableContainer,'Sagittal')
    %add the Sagittal Fields
    str = cat(1, str,...
        {'Sagittal::LensDiameter','Sagittal::LensThickness',...
        'Sagittal::Fouier10_LensThickness', 'Sagittal::Fourier10_LensDiameter',...
        'Sagittal::Fouier20_LensThickness', 'Sagittal::Fourier20_LensDiameter',...
        'Sagittal::Fit: Conic R', 'Sagittal::Fit: Conic p',...
        'LegacySagittal::Fouier10_LensThickness', 'LegacySagittal::Fourier10_LensDiameter',...
        'LegacySagittal::Fourier20_LensThickness','LegacySagittal::Fourier20_LensDiameter',...
        'LegacySagittal::conic_p','LegacySagittal::conic_R'}');
end

set(hListbox_right,'String', str); 

function b_addField_Callback(hObject, eventdata, handles)
% hObject    handle to b_addField (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hPopup = handles.popup_tableList2;
hListbox_left = handles.listbox_fieldList2;
hListbox_right = handles.listbox_tableColumns;

%Get selected field name
thisVal_left = get(hListbox_left, 'Value');
if strcmpi(get(hListbox_left, 'Enable'),'off') || thisVal_left(1) == 0
    return
end

listboxStr_left = cellstr(get(hListbox_left, 'String'));
thisField_left = listboxStr_left(thisVal_left); % can only be one

parentTableName = getappdata(hPopup,'parentTableName'); %char array

if isempty(parentTableName)
    parentTableName = updateFieldNames(hPopup);
end

newStr = cellfun(@(iStr) sprintf('%s::%s', parentTableName, iStr),...
    thisField_left,'UniformOutput', false);


%Append it to the table on the right
currentStr = get(hListbox_right,'String');

if isempty(currentStr)
    listboxStrOut = newStr;
else
    listboxStrOut = cat(1,force1D(currentStr), force1D(cellstr(newStr)));
end

%remove duplicates
listboxStrOut = unique(listboxStrOut,'stable');

%Update the string
set(hListbox_right, 'String', listboxStrOut,'Enable','on', 'Max', numel(listboxStrOut));
setappdata(ancestor(hObject,'figure'),'settingsFlag',true);

function b_removeField_Callback(hObject, eventdata, handles)
% hObject    handle to b_removeField (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hListbox_right = handles.listbox_tableColumns;

%Get selected field name
thisVal_right = get(hListbox_right, 'Value');
if strcmpi(get(hListbox_right, 'Enable'),'off') || isempty(thisVal_right) || thisVal_right(1) == 0
    return
end

listboxStr_right = get(hListbox_right, 'String');

if isempty(listboxStr_right)
    %there's nothing to delete
    set(hListbox_right,'Enable','off');
    return
end


currentStr = get(hListbox_right, 'String'); newStr = currentStr;

newStr(thisVal_right) = [];

if isempty(newStr)
     enableStr = 'off';
else
    enableStr = 'on';
end

%Update handles, listbox, and structure
set(hListbox_right, 'String',newStr, 'Enable', enableStr,'Value',1,...
    'Max', numel(newStr));
setappdata(ancestor(hObject,'figure'),'settingsFlag',true);

function moveColumnName(hObject, eventdata, handles)
%This function moves a column up or down in the tableColumns listbox

hListbox_right = handles.listbox_tableColumns;

%Get selected field name
thisVal_right = sort(get(hListbox_right, 'Value'));
if strcmpi(get(hListbox_right, 'Enable'),'off') || thisVal_right(1) == 0
    return
end

listboxStr_right = cellstr(get(hListbox_right, 'String'));

if isempty(listboxStr_right)
    %there's nothing to move
    set(hListbox_right,'Enable','off');
    return
end

currentStr = get(hListbox_right, 'String'); newStr = currentStr;
nElements = numel(listboxStr_right);

dir = regexp(lower(get(hObject, 'Tag')),'(up|down)','match');
dir = dir{1};

newValSelected = ones(numel(thisVal_right),1);

for j = 1:numel(thisVal_right)
    jVal = thisVal_right(j);
    switch dir
        case 'up'
            jNewPosIdx = jVal - 1;
        case 'down'
            jNewPosIdx = jVal+1;
    end
    bothIdx = sort([jVal, jNewPosIdx]);
    
    
    if jNewPosIdx > nElements || jNewPosIdx < 1
        newValSelected(j) = jVal;
        %no room to move
        continue
    else
        newValSelected(j) = jNewPosIdx;
    end
    
    newStr(bothIdx) = flipud(newStr(bothIdx));
    
    
end

%Set the new data settings
set(hListbox_right, 'String',newStr,'Value', newValSelected);

setappdata(ancestor(hObject,'figure'),'settingsFlag',true);

%..........................................................................
%Filter Subjects
%..........................................................................

function popup_filterGeneric_Callback(hObject, eventdata, handles)
%This function is the callback for a filter field callback

hPanel = handles.panel_filter;
filters = getappdata(hPanel, 'filters');
if ~isempty(hObject) && ishandle(hObject)
    thisTag = regexp(get(hObject,'Tag'),'popup_filter(.*)','tokens');
    callbackName = lower(char(thisTag{1}));
    
    
    if strcmpi(get(hObject, 'Style'),'popupmenu');
        allStr = get(hObject, 'String');
        selectedVal = get(hObject, 'Value');
        selectedStr = allStr{selectedVal}; %as a char array
        
    else
        selectedStr = get(hObject,'String');
    end
    
    allFlag = strcmpi(selectedStr, 'All') || strcmpi(selectedStr,'None')...
        || isempty(selectedStr);
else
    return
end

if allFlag
    filters.(callbackName) = [];
else
    
    switch callbackName
        case {'eye','imgtype'}
            filters.(callbackName) = selectedStr;
        case {'segmented','validated','excluded'}
            filters.(callbackName) = isempty(strfind(lower(selectedStr), 'not'));
        case {'notes'}
            filters.(callbackName) = selectedStr;
        otherwise
            %do nothing
    end
end

setappdata(hPanel, 'filters', filters);

if ~strcmpi(eventdata, 'manualcall')
    populateSubjectListbox(handles);
    updateGraphPanel(handles);
end

function [isOk, subjectTableOut] = filterSubjectTable(handles, subjectTableIn, filters)

hPanel = handles.panel_filter;
hFig = ancestor(handles.axes_img,'figure');
tableContainer = getappdata(hFig, 'tableContainer');

%Initialize Variable
isOk = true;
subjectTableOut = [];

if isempty(tableContainer) || ~isfield(tableContainer, 'Experiment')
    %nothing loaded
    isOk = false;
    return
end


if nargin < 3
    filters = getappdata(hPanel, 'filters');
end

if nargin < 2 || isempty(subjectTableIn)
    subjectTableIn = tabelContainer.Subject;
end

if ~isempty(filters)
   
    expTable = tableContainer.Experiment;    
    isThisRow = true(size(expTable,1),1);
    
    allExpIDs = expTable.ExperimentID;
    
    
    %Segmentation Filter
    if ~isempty(filters.segmented)
        hasSegmentation = isfield(tableContainer, 'Sagittal');
        

        if filters.segmented %we want those that are segmented
            if hasSegmentation
                validExpIDs = intersect(allExpIDs, tableContainer.Sagittal.ExperimentID);
            else
                %return nothing, but isOk = true;
                subjectTableOut = [];
                setappdata(hPanel,'ExperimentID','none');
                return
            end
            
        else %we want those that arent segmented
            if ~hasSegmentation
                %all subjects valid
                validExpIDs = allExpIDs;
            else 
                %only those that arent in the sagittal table are valid
                validExpIDs = setdiff(allExpIDs, tableContainer.Sagittal.ExperimentID);
            end
            
        end
        
        isThisRow = bitand(isThisRow,...
            cellfun(@(iID) any(strcmpi(validExpIDs, iID)), allExpIDs));
    end
    

    %Filter Eye
    if ~isempty(filters.eye)
        isThisRow = bitand(isThisRow, expTable.Eye == filters.eye);
    end
    
    %Filter ImgType
    if ~isempty(filters.imgtype)
        isThisRow = bitand(isThisRow, expTable.ImgType == filters.imgtype);
    end

    %Filter Validation
    if ~isempty(filters.validated)
        isThisRow = bitand(isThisRow, expTable.ValidatedFlag == filters.validated);
    end
    
    
    %Filter Excluded
    if ~isempty(filters.excluded)
        isThisRow = bitand(isThisRow, all(expTable.OutlierFlag,2) == filters.excluded);
    end
    
    %Filter Notes
    if ~isempty(filters.notes)
        isThisRow = bitand(isThisRow, cellfun(@isempty,...
            strfind(lower(expTable.Notes), lower(strtrim(filter.notes)))));
    end
    
    %None of the images valid
    if ~any(isThisRow)
        subjectTableOut = [];
        setappdata(hPanel,'ExperimentID','none');
        return
    end
    
    
    setappdata(hPanel,'ExperimentID',expTable.ExperimentID(isThisRow));
    validSubjectIDs = expTable.SubjectID(isThisRow);
    uniqueSubjectIDs = unique(validSubjectIDs,'stable');
    
    isSubjectValid = cellfun(@(iID) any(strcmpi(iID,uniqueSubjectIDs)),...
        subjectTableIn.SubjectID);
    subjectTableOut = subjectTableIn(isSubjectValid,:);
    
    
else
    %nothing to filter
    subjectTableOut = subjectTableIn;
    return
    
end

function [dataTableOut,dataStructure] = readTableSettings(handles)
%This function extracts the values and data from the tableList

hListbox_right = handles.listbox_tableColumns;

%Initialize Output
dataStructure = [];
dataTableOut = [];

if numel(get(hListbox_right, 'String')) < 1
    return
end

%get the table_field pairs and make sure there are no duplicates
parentsANDfields_string = unique(get(hListbox_right, 'String'),'stable');

tableContainer = getappdata(ancestor(handles.axes_img, 'figure'),'tableContainer');
lookUpTable = tableContainer.Experiment(:, {'SubjectID','ExperimentID'});

finalTableOrder = [];

for i = 1:numel(parentsANDfields_string)
    iPair = parentsANDfields_string{i};
    
    iToken = regexp(iPair,'(.*)::(.*)','tokens');
    iParent = iToken{1}{1}; iField = iToken{1}{2};
    
    %convert the selected fields to data structure
    iData = analyzeFieldString(handles, iParent, iField);
    iData.Parent = iParent; iData.Field = iField;
    iKeyName = iData.KeyName; iKeyValue = iData.KeyValue;

    
    %make sure the field type is compatable with matlab
    switch lower(class(iData.Value))
        case 'categorical'
            iData.Value = cellstr(iData.Value);
        case 'datetime'
            iData.Value = cellstr(iData.Value);
        otherwise
            %do nothing
    end

            
    if strcmpi(iField, iData.KeyName)
        iTable = table(iData.Value, 'VariableNames', {iField});
        finalTableOrder = cat(2, finalTableOrder, cellstr(iField));
    else
        iTableLabel = strtrim(strrep(iField,'Fit:',''));
        iTableLabel = strrep(iTableLabel,' ','_');
        
        
        if iscell(iData.Value) && istable(iData.Value{1})
            continue %do not parse tables
        elseif size(iData.Value, 2) == 2
            iTable = table(iKeyValue, iData.Value(:,1), iData.Value(:,2),...
                'VariableNames', {iKeyName,...
                [iParent(1:3), '_Ant_', iTableLabel], [iParent(1:3),'_Pos_', iTableLabel]});
        else
            iTable = table(iKeyValue, iData.Value,...
                'VariableNames', {iKeyName, [iParent(1:3),'_', iTableLabel]});
        end
        finalTableOrder = cat(2, finalTableOrder, iTable.Properties.VariableNames(2:end));
    end
    
    %keep only unique rows of iTable
    
    switch iParent            
        case 'Sagittal'
            %keep only the most recent
            iParentTable = tableContainer.(iParent);
            [sortedTable, iSortIdx] = sortrows(iParentTable, 'CreatedDate','descend');
            [~, iKeepIdx] = unique(sortedTable.ExperimentID, 'first');
            iTableRowIdx = iSortIdx(iKeepIdx);
            iTable = iTable(sort(iTableRowIdx),:);
        otherwise
            [~, iKeepIdx] = unique(iTable(:,1),'stable');
            iTable = iTable(iKeepIdx,:);
    end
    
    if i == 1
        dataStructure = iData;
        dataTable = iTable;
    else
        dataStructure = cat(1,force1D(dataStructure), force1D(iData)); 
        
        %cat the table with 
        if any(ismember(dataTable.Properties.VariableNames, iKeyName))
            dataTable = outerjoin(dataTable, iTable,'MergeKeys', true);
        else
            %use lookup table
            if strcmpi(iKeyName, 'ExperimentID')
                
                firstTable = iTable;
                secondTable = dataTable;
           
            elseif ismember(dataTable.Properties.VariableNames, 'ExperimentID')
                firstTable = dataTable;
                secondTable = iTable;
            else
                %this should never happen.
                return
            end
            
            dataTable = join(join(firstTable, lookUpTable,'LeftKeys','ExperimentID',...
                'RightKeys','ExperimentID'),secondTable,'LeftKeys','SubjectID',...
                'RightKeys','SubjectID');

        end
        
    end
    
end

%reorganize table fields
dataTableOut = dataTable(:,finalTableOrder);

%Store the information
setappdata(hListbox_right,'data', dataStructure);
setappdata(hListbox_right,'dataTable', dataTableOut)


%..........................................................................
function populateTableList(handles, tableContainer)
%This function poupulates the table list popup

hFig = ancestor(handles.axes_img, 'figure');

if nargin < 2
    tableContainer = getappdata(hFig, 'tableContainer');
end
hPopups = [handles.popup_tableList1, handles.popup_tableList2];
hListboxes = [handles.listbox_fieldList1, handles.listbox_fieldList2];

%Make sure the graph updates
setappdata(hFig, 'settingsFlag',true);

if isempty(tableContainer)
    %Clear fieldX
    set(handles.popup_xField,'String','No fields found','Enable','off',...
        'Value',1);
    
    %Clear all the data from the table fields
    arrayfun(@(iH) set(iH, 'String','No tables found.'), hPopups);
    set(hPopups, 'Enable','off','Value',1);
    
    %Clear listbox string and appdata
    arrayfun(@(iH) set(iH, 'String','No fields found.'), hListboxes);
    set(hListboxes, 'Enable','off','Value',1, 'Max',1);
    arrayfun(@(iH) setappdata(iH, 'parentTable',''), hListboxes);
    return
end

tableNames = fieldnames(tableContainer);
arrayfun(@(iH) set(iH, 'String', tableNames,'Value',1,...
    'Enable','on'), hPopups);

updateFieldNames(hPopups(1));
updateFieldNames(hPopups(2));

%This also populates the x field popup list
set(handles.popup_xField, 'Enable','on','Value',1);



fieldNameContainer = updateFieldNamesContainer(handles, tableContainer);

parentTableList = fieldnames(fieldNameContainer);
nElements = structfun(@numel, fieldNameContainer);
parentTables = arrayfun(@(idx) repmat(parentTableList(idx),...
    [nElements(idx),1]),1:numel(parentTableList), 'UniformOutput', false);
parentTables = cat(1, parentTables{:});

groupedFields = struct2cell(fieldNameContainer);
groupedFields = cat(2,groupedFields{:});
parent_child = [parentTables,groupedFields'];

xFieldStr = arrayfun(@(idx) strjoin(parent_child(idx,:),'::'),1:size(parent_child, 1),...
    'UniformOutput', false)';

setappdata(handles.popup_xField, 'parentTables', parentTables);
setappdata(handles.popup_xField, 'fields', groupedFields');

%Set the string
isAgeField = strcmpi(groupedFields, 'Age'); %make this the default X

if any(isAgeField)
    initialVal = find(isAgeField);
else
    initialVal = 1;
end

set(handles.popup_xField, 'String',xFieldStr, 'Value',initialVal,'Enable','on');
%Run the callback
popup_xField_Callback(handles.popup_xField, [], handles);

%Initialize y-field
if isempty(get(handles.edit_yField, 'String'))
    set(handles.edit_yField, 'String','<<LegacySagittal::Fourier20_Volume>>');
end

%Inialize table Columns?
x=1;


%..........................................................................

function [parentTableName, fieldName, match] = parseFieldString(stringIn)
%This function parses the user input/field string <<parentTable::fieldName>>


if iscell(stringIn)
    stringIn = char(stringIn);
end

% tokens = regexp('<<Sagittal:Fit R_ant>>','<<(\w.*):(\w.*)>>','tokens')
% regexp('Sagittal:Volume','<?<?(\w*):(\w*)>?>?','tokens')
% regexp('<<Sagittal::Volume>> + <<Car::Dog>>','<?<?(\w*)::(\w*)>?>?','tokens')

[tokens, match] = regexp(strtrim(stringIn),'<?<?(.*)::(.*[^>])>?>?','tokens','match');

if isempty(tokens)
    parentTableName = [];
    fieldName = [];
    match = [];
    return
end


tokens_out = cat(1,tokens{:});

parentTableName = tokens_out(:,1); fieldName = tokens_out(:,2);

function output = evalFieldExpression(tableContainer, yFieldStrIn)
%Parse the yFieldStr.  Different expressions are separated using a
%semi-colon (;)

yFieldStr = strsplit(yFieldStrIn, ';');
nExpressions = numel(yFieldStr);

output = struct('Value',[], 'Label','','Units','');
for j = 1:nExpressions
    jExp = yFieldStr{j};
    [jParentName, jFieldName, jMatch] = parseFieldString(jExp);
    
    nTables = numel(jParentName);
    
    if ~isfield(tableContainer, jParentName)
        output(j).Label = sprintf('%s::%s',jParentName{1},jExp);
        continue
    end
    
    kVals = cell(1,nTables); kUnits= cell(1,nTables);
    jLabel = jExp;
    for k = 1:nTables
        kTable = tableContainer.(jParentName{k});
        [kVals{k}, kUnits{k}] = extractTableValue(kTable,...
            jFieldName{k});
        kValStr = strrep(jMatch{k}, jMatch{k},sprintf('kVals{%d}',k));
        jLabel = strrep(jLabel, jMatch{k}, jFieldName{k});
    end
    
    try
        jOutput = eval(kValStr);
    catch
        fprintf(1,'\n\nError evaluating the expression: \n%s\n\n', kValStr);
    end
    
    output(j).Value = jOutput;
    output(j).Label = sprintf('%s::%s',jParentName{1},jLabel);
    output(j).Units = kUnits{1}; %This may not be accurate
    output(j).OutlierFlag = [];
    
    %Add the key value and name
    try
        output(j).KeyName = 'ExperimentID';
        output(j).KeyValue = kTable.ExperimentID;
    catch
        output(j).KeyName = 'SubjectID';
        output(j).KeyValue = kTable.SubjectID;
    end
    
end


%..........................................................................

function dataStructure = analyzeFieldString(handles, parentTableName, fieldNames)
%This function gets data from tables, if no tableName is given, it looks
%through all tables/
%OUTPUT
    %dataStructure - structure with fields: Value, Units, Label


%Initialize Output
dataStructure = struct([]);

tableContainer = getappdata(ancestor(handles.axes_img,'figure'),'tableContainer');

if isempty(tableContainer)
    return
end

if ~iscell(fieldNames)
    fieldNames = {fieldNames};
end

nFields = numel(fieldNames);

%Make sure the number of parent tables match number of fields
if ~iscell(parentTableName) || (numel(parentTableName) == 1 && nFields > 1)
    parentTableName = cellstr(parentTableName);
    parentTableName = repmat(cellstr(parentTableName), nFields, 1);
end


for i = 1:nFields
    iField = fieldNames{i};   
    iParentName = parentTableName{i};
    
    iTable = tableContainer.(iParentName);
    
    isCoeffField = ~isempty(strfind(iField, 'Fit'));
    
    if ~isCoeffField
        [iValOut, iUnitsOut] = extractTableValue(iTable, iField);
    else
        %if it's a coeffcient
        iFieldC = regexp(iField, '[pPrR]','match');
        [iValOut, iUnitsOut] = extractTableValue(iTable, ['Fit: ',iFieldC{1}]);
    end
    
    dataStructure(i).Value = iValOut; %This should work for everything else
    dataStructure(i).Units = iUnitsOut;
    dataStructure(i).Label = sprintf('%s::%s',iParentName,iField);
    dataStructure(i).OutlierFlag = [];
    
    try
        dataStructure(i).KeyName = 'ExperimentID';
        dataStructure(i).KeyValue = iTable.ExperimentID; %try this first
    catch
        dataStructure(i).KeyName = 'SubjectID';
        dataStructure(i).KeyValue = iTable.SubjectID; %all tables have subjectID
    end

    
end

%..........................................................................
function selectedTableName = updateFieldNames(hObject, eventdata, handles)
%This function is for updating the field names of the child listbox
%associated with hObject (popup parent). Inputs don't really matter.  It
%is the callback for popup_tableList selection change
%
%hObject --> handles.popup_tableList#(1/2)
%
%NOTE: A table should never be empty

hFig = ancestor(hObject, 'figure');

if nargin < 3
    handles = guidata(hFig);
end

%Get the name of the selected table
allTables = get(hObject,'String');
selectedTableName = allTables{get(hObject, 'Value')};

fieldNameContainer = getappdata(hFig,'fieldNameContainer');

if isempty(fieldNameContainer)
    fieldNameContainer = updateFieldNamesContainer(handles);
end

%Get the actual table
newFieldNames = fieldNameContainer.(selectedTableName);

hTag = get(hObject, 'Tag');
listboxNo = hTag(end);
hPopup = eval(['handles.listbox_fieldList',listboxNo]);

setappdata(hPopup, 'parentTableName',selectedTableName);
set(hPopup,'String', newFieldNames, 'Value',1,'Enable','on',...
    'Max',numel(newFieldNames));

function fieldNameContainer = updateFieldNamesContainer(handles, tableContainer)
%This function creates a structure with the same fields as table container
%to hold the fieldNames and puts in the appdata of the figure

hFig = ancestor(handles.axes_img, 'figure');

if nargin < 2
    tableContainer = getappdata(hFig, 'tableContainer');
end

if isempty(tableContainer)
    setappdata(hFig, 'fieldNameContainer',[]);
    return
end

fieldNameContainer = [];
allTableNames = fieldnames(tableContainer);

for i = 1:numel(allTableNames);
    iTableName = allTableNames{i};
    iTable = tableContainer.(iTableName);
    
    if isempty(iTable)
        continue
    end
    
    %Get the names of the fields
    iFieldNames = iTable.Properties.VariableNames;
    
    %parse fieldnames    
    switch iTableName
        case 'Subject'
            %all fields are ok to display
            %SubjectID, Age, Species
        case 'Experiment'
            %should be ok, not all helpful, but ok
            %SubjectID, ExperimentID, PMT, ExperimentDated, ImgType,
            %Img_originalLocation, Img_filename, Img_serverDirectory,
            %Img_localDirectory, AcquiredBy
            
            iFields2remove = {'OutlierFlag'};
            
            iFieldNames = setdiff(iFieldNames, iFields2remove);
            
        case 'LegacyCoronal'
            %All Ok
        case 'LegacySagittal'
            %remove the points & coefficients
            iFields2remove = {'OriginalPoints', 'CorrectedPoints',...
                'Fourier10Coeff','Fourier20Coeff', 'Fourier10RMSE', 'Fourier20RMSE'};
            
            iFieldNames = setdiff(iFieldNames, iFields2remove);
        case 'Sagittal'
            
            %remove all the params
            isRemoveField = ~cellfun(@isempty, strfind(iFieldNames, 'Param_'));
            
            iFieldNames(isRemoveField) = [];
            
            %Remove additional fields
            iFields2remove = {'TotalShift', 'ImageShift', 'ImageSize',...
                'Fourier10Coeff', 'Fourier20Coeff', 'Fourier10RMSE', 'Fourier20RMSE',...
                };
            
            iFieldNames = setdiff(iFieldNames, iFields2remove);
            
            %add the fit coeffiecients
            coeffDescr = cellfun(@(coeffTable) cat(2, coeffTable.FitName,...
                coeffTable.Variable), iTable.Coefficients(1),...
                'UniformOutput', false); coeffDescr = coeffDescr{1};
            
            coeffStr = arrayfun(@(iRow) strjoin([{'Fit:'},...
                cellstr(coeffDescr(iRow,:))]), 1:size(coeffDescr,1),...
                'UniformOutput', false); 
            
            iFields2add = unique(coeffStr);
            
            iFieldNames = union(iFieldNames, iFields2add);
    end
    
    %Put the fieldNames
    fieldNameContainer.(iTableName) = iFieldNames;
    
end

%Put it in the figure appdata
setappdata(hFig,'fieldNameContainer', fieldNameContainer);

function fieldNameSelectionChange

hTag = get(hObject,'Tag');
listboxNo = get(hTag(end));

hListbox = ['handles.listbox_fieldList',listboxNo];

if ~ishandle(hListbox)
    %something is wrong
    return
end

switch str2double(get(hTag(end)))
    case 1
    case 2
end


%--------------------------------------------------------------------------
%WORKING WITH SUMMARY TABLE
%--------------------------------------------------------------------------

function updateSummaryTable(handles, experimentID)

hFig = ancestor(handles.axes_img,'figure');
hAx = handles.axes_graph;
tableContainer  = getappdata(hFig, 'tableContainer');
hTable = handles.table_single;

if isempty(tableContainer) || ~isfield(tableContainer, 'Experiment')
    clearTable(hTable);
    return
end


%Get the currently selected experimentID if user doesn't input one
if nargin < 2
    hListboxE = handles.listbox_experiment;
    selectedVal = get(hListboxE, 'Value');
    expTable = getappdata(hListboxE,'table');
    
    if ~isempty(selectedVal) && selectedVal > 0 && ~isempty(expTable)
        experimentID = expTable.ExperimentID(selectedVal);
    else
        experimentID = '';
    end
    expRow = expTable(selectedVal, :);
else
    
    expTable = tableContainer.Experiment;
    isThisID = strcmpi(expTable.ExperimentID,experimentID);
    
    expRow = expTable(isThisID, :);
    
end

%Get the Sagittal Data or Coronal
if isempty(experimentID) || isempty(expRow)
    clearTable(hTable);
    return
end



if expRow.ImgType(1) == 'Sagittal'
    legacyTableFields = {'Experiment::PMT',...'<<Subject::Age>>','<<Experiment::Eye>>',...
        '<<LegacySagittal::Fourier20_LensThickness>>',...
        '<<LegacySagittal::Fourier20_LensDiameter>>',...
        '<<LegacySagittal::conic_R>>','LegacySagittal::conic_p'};
    
    currentTableFields = {'Experiment::PMT','<<Sagittal::LensThickness>>',...'<<Subject::Age>>','<<Experiment::Eye>>',...
        '<<Sagittal::LensDiameter>>',...
        '<<Sagittal::Fit_Conic_R>>','Sagittal::Fit_Conic_p'};
    
    summaryTableColumnNames = {'PMT', 'Lens T', 'Lens D',...'Age', 'Eye', 
        'Ant R', 'Pos R','Ant p', 'Pos p'};
    
    
elseif expRow.ImgType(1)  == 'Coronal'
    legacyTableFields = {'Experiment::PMT','LegacyCoronal::Coronal_AvgR',...'<<Subject::Age>>','<<Experiment::Eye>>',...
        'LegacyCoronal::Coronal_MaxR','LegacyCoronal::Coronal_MinR'...
        'LegacyCoronal::Coronal_StdR', 'LegacyCoronal::Coronal_CSA',...
        'LegacyCoronal::Coronal_Circumference'};
    
    currentTableFields = {}; %nothing for now
    
    summaryTableColumnNames = {'PMT', 'Avg R',...'Age', 'Eye', 
        'Max R', 'Min R', 'Std R','CSA','Circum.'};
    
else
     clearTable(hTable);
    return
end


if ~isempty(legacyTableFields)
    [finalLegacyTable, finalLegacyCell] = generateTableFromFields(tableContainer,...
    legacyTableFields, experimentID);
else
    finalLegacyTable = [];
    finalLegacyCell = [];
end


if ~isempty(currentTableFields)
    [finalCurrentTable, finalCurrentCell] = generateTableFromFields(tableContainer,...
        currentTableFields, experimentID);
else
    finalCurrentTable = [];
    finalCurrentCell = [];
end 

rowNames = {'New','Legacy'};
rowNames([isempty(finalCurrentCell), isempty(finalLegacyCell)]) = [];

try
    summaryTableData = cat(1, finalCurrentCell, finalLegacyCell);
    
    %Set the table data
    set(hTable, 'ColumnName',summaryTableColumnNames,'Data', summaryTableData,...
        'RowName',rowNames);
    
    
    setappdata(hTable, 'legacyTable', finalLegacyTable);
    setappdata(hTable, 'currentTable', finalCurrentTable);
    setappdata(hTable, 'imgType', expRow.ImgType(1));
    setappdata(hTable, 'expRow', expRow);
catch
    set(hTable, 'Data',[]);
    
end

%Set the data for subject info
subTable = tableContainer.Subject;
thisSubjectTable = join(expRow, subTable);
thisSubjectTable = thisSubjectTable(1,:); %just encase there is more than one


%Set the edit boxes
set(handles.edit_subjectID, 'String',thisSubjectTable.SubjectID);
setappdata(handles.edit_subjectID,'prevValue',thisSubjectTable.SubjectID);

set(handles.edit_experimentID, 'String', experimentID);
setappdata(handles.edit_experimentID,'prevValue',experimentID);

set(handles.edit_age, 'String', thisSubjectTable.Age);
setappdata(handles.edit_age,'prevValue',thisSubjectTable.Age);

set(handles.edit_notes, 'String', expRow.Notes{1});
setappdata(handles.edit_notes, 'prevValue', expRow.Notes{1});

%set the outlierFlags
outlierFlag = expRow.OutlierFlag;
set(handles.check_flagAnt,'Value', outlierFlag(1));
set(handles.check_flagPos,'Value', outlierFlag(2));

function clearTable(hTable)

set(hTable, 'Data', []);
setappdata(hTable, 'legacyTable', []);
setappdata(hTable, 'currentTable', []);
setappdata(hTable, 'imgType', []);
setappdata(hTable, 'expRow', []);


function [finalTable, finalCell] = generateTableFromFields(tableContainer, stringIn, experimentID)
%This function generates a table given the string

if nargin < 2
    experimentID = '';
end


if ~iscell(stringIn)
    stringIn = cellstr(stringIn);
end

%Get the names of the parents and fields in original order
[allParents_cell, allFields_cell] = cellfun(@(iStr) parseFieldString(iStr),...
    stringIn, 'UniformOutput', false);

allParents = cat(1,allParents_cell{:});
allFields = cat(1, allFields_cell{:});

%remove parentTable if it's not in the table container
existingTables = fieldnames(tableContainer);
invalidTables = setdiff(allParents, existingTables);


if ~isempty(invalidTables)
    isParentValid = force1D(~cellfun(@(iParent) any(strcmpi(invalidTables, iParent)),...
        allParents));
else
    isParentValid = true(size(allParents));
end



if ~all(isParentValid)
    %nothing exists, exit
    finalTable = [];
    finalCell = [];
    return
end

uniqueValidParents = unique(allParents(isParentValid));

% nTables = numel(uniqueValidParents);


%Make sure the sagittal table only has unique values

%Create a mega table (should be alphabetized, Experiment should always be
%the first table
tableContainer = makeKeyValuesUnique(tableContainer);

tableOut = cellfun(@(iName) tableContainer.(iName), uniqueValidParents,...
    'UniformOutput', false);

megaTable = tableOut{1}; %This will contain duplicates for various fields
for i = 2:numel(tableOut)
    
    if isempty(tableOut{i})
        continue
    end
    megaTable = outerjoin(megaTable, tableOut{i},'Type','left','MergeKeys',true);
end

megaColNames = megaTable.Properties.VariableNames';
isFieldValid = cellfun(@(iField) any(strcmpi(megaColNames, iField)),...
    allFields);


%keep only the necessary rows
if ~isempty(experimentID)
    %only keep the rows that match the experimentID
    isThisID = strcmpi(megaTable.ExperimentID,experimentID);
    
    if any(isThisID)
        megaTable = megaTable(isThisID,:);
    end 
end

%add coefficient fields to the mega table
isCoeffField = isParentValid & ~isFieldValid;
if any(isCoeffField) %coeffients
    
    %must use extractTableValue function
    coeffField = allFields(isCoeffField);
    coeffColIdx = find(isCoeffField);
    for i = 1:numel(coeffField)
        iField = coeffField(i);
        iValue = extractTableValue(megaTable, iField);
        
        
        if ~isempty(iValue)
            megaTable{:,iField} = iValue;
            isFieldValid(coeffColIdx(i)) = true;
        end
    end
    
    
end


%Extract all the fields needed
validFields = allFields(isFieldValid);
finalTable = megaTable(:, validFields);


% for i = 1:nTables
%     iTableName = uniqueValidParents(i);
%     isThisFields = strcmpi(allParents, iTableName);
%     iFields = allFields(isThisFields);
%     
%     iTableOut = megaTable(:, iFields);
%     output{i} = iTableOut;
%     colIdx{i} = find(isThisFields);
% end
% 
% %Concatenate each table
% finalTable = cat(2,output{:});
% 




% finalColumnOrder = cell2mat(colIdx);
nRows = size(finalTable,1);
nColumns = numel(stringIn);
cellOutput = cell(nRows,nColumns);
% finalCell(:, finalColumnOrder) = table2cell(finalTable);
cellOutput(:, isParentValid) = table2cell(finalTable);

finalCell = {};
for j = 1:size(cellOutput,2)
    jCell = cellOutput(:,j);
    
    if strcmpi(class(jCell{1}),'categorical')
        jCell = force1D(cellfun(@char, jCell,'UniformOutput', false),2);
        finalCell = cat(2,finalCell,vertcat(jCell{:}));
    else
        finalCell = cat(2,finalCell, num2cell(vertcat(jCell{:})));
    end
    
end

function isOK = setTableValue(handles, tableName, fieldName, fieldValue, keyValue)

hFig = ancestor(handles.axes_img,'figure');
tableContainer = getappdata(hFig, 'tableContainer');

%Initilalize Output
isOK = false;

if isempty(tableContainer)
    return
end

if strcmpi(tableName, 'Subject')
    keyName = 'SubjectID';
else
    keyName = 'ExperimentID';
end
    
if nargin < 5
   if strcmpi(keyName,'SubjectID')
       hListboxS = handles.listbox_subject;
       subTable = getappdata(hListboxS, 'table');
       selectedVal = get(hListboxS, 'Value');
       if isempty(subTable) || isempty(selectedVal) || selectedVal < 1
           return
       end
       
       keyValue = subTable.SubjectID(selectedVal(1));
       
       if iscell(keyValue)
           keyValue = char(keyValue);
       end
       
       
   else
       hListboxE = handles.listbox_experiment;
       expTable = getappdata(hListboxE, 'table');
       selectedVal = get(hListboxE, 'Value');
       if isempty(expTable) || isempty(selectedVal) || selectedVal < 1
           return
       end
       
       keyValue = expTable.ExperimentID(selectedVal(1));
       
       if iscell(keyValue)
           keyValue = char(keyValue);
       end
       
   end
end


thisTable = tableContainer.(tableName);

if isempty(keyValue) && nargin == 5
    %Set all rows of table
    isThisRow = true(size(thisTable,1),1);
else
    isThisRow = strcmpi(thisTable.(keyName), keyValue);
end

if ~any(isThisRow)
    isOK = [];
    return
end

thisTable{isThisRow,fieldName} = fieldValue;

tableContainer.(tableName) = thisTable;

setappdata(hFig, 'tableContainer', tableContainer);
isOK = true;



function [valOut, unitsOut] = extractTableValue(tableIn, fieldName)
%This function extracts specific values from a table specified by
%fieldName.  This does not accept multiple tables

if iscell(fieldName)
    fieldName = fieldName{1};
end

colNames = tableIn.Properties.VariableNames;
colDescription = tableIn.Properties.VariableDescriptions;
colUnits = tableIn.Properties.VariableUnits;


if ~isempty(strfind(fieldName,'Fit')) 
    if ~isempty(strfind(lower(fieldName), 'r'))
        %retrieve it from the coefficient table
        if ~isempty(strfind(lower(fieldName), 'poly'))
            ftype = 'Poly4';
        else
            ftype = 'Conic';
        end
        
        R_cell = cellfun(@(iTable) getCoeff(iTable, 'R','Corrected',ftype),...
            tableIn.Coefficients, 'UniformOutput', false);
        R_mat = cat(1,R_cell{:});
        
        unitsOut = 'mm';
        valOut = R_mat;
    elseif ~isempty(strfind(lower(fieldName), 'p'))
        %retrieve it from the coefficient table
        p_cell = cellfun(@(iTable) getCoeff(iTable, 'p','Corrected','Conic'),...
            tableIn.Coefficients, 'UniformOutput', false);
        p_mat = cat(1,p_cell{:});
        
        unitsOut = '';
        valOut = p_mat;
    else 
        %I have no idea what it is
        unitsOut = '';
        valOut = {};
        disp(['Field Not Found ',fieldName]);
    end
else
        
        isThisCol = strcmpi(colNames, fieldName);
        if any(isThisCol)
            
            %Get the Units
            thisUnits = colUnits(isThisCol);
            unitsOut = thisUnits{1};
            valOut = tableIn.(fieldName);
            
        else
            %what if I can't find it
            unitsOut = '';
            valOut = {};
            disp(['Field Not Found ',fieldName]);
        end
        
end

function val = getCoeff(fitsTable, varName, varDescription, fitType)
%This function returns the coefficient value with name and description
%input.  If it is empty, it returns an NaN
%INPUT
    %tableIn - fits table (not as a cell)
    %varName - categorical string 'R', 'p'
    %varDescription - character a ray {not case sensitive}
    %fitType - categorical string 'Cosine', 'Conic'
%OUTPUT
    %val - value of the coefficients [ant, pos] where it exists

    if isempty(fitsTable)
        val = NaN;
        return
    end
    if nargin > 3 %all filters
        isThisVar = fitsTable.Variable == varName & fitsTable.FitName == fitType &... %case sensitive
    ~cellfun(@isempty,strfind(lower(fitsTable.InputDescription), lower(varDescription))); %case insensitive
    elseif nargin > 2 %varName, varDescrip
        isThisVar = fitsTable.Variable == varName &... %case sensitive
    ~cellfun(@isempty,strfind(lower(fitsTable.InputDescription), lower(varDescription))); %case insensitive
    elseif nargin > 1 %varName
        isThisVar = fitsTable.Variable == varName; %case sensitive
    else %all Rows
        isThisVar = true(size(fitsTable), 1);
    end


if any(isThisVar)
    val = fitsTable(isThisVar,'Value');
    val = val.Value(1,:); %only the first match
else
    firstVal = fitsTable(1,'Value');
    firstVal  = firstVal.Value(1,:);
    val = nan(size(firstVal));
end

function tableContainer = makeKeyValuesUnique(tableContainer)
%This function removes duplicates of Sagittal table and keeps only the
%latest so that the join function can be used

if isfield(tableContainer,'Sagittal')
    sagTable = sortrows(tableContainer.Sagittal,{'ExperimentID','CreatedDate'},'descend');
    [~, ia, ~] = unique(sagTable.ExperimentID,'first');
    tableContainer.Sagittal = sagTable(ia,:);
end


%--------------------------------------------------------------------------
%WORKING WITH IMAGES
%--------------------------------------------------------------------------

function displaySingleSubject(handles, expRow)
%The input expRow MUST be a row from the Experiment table

%Get handles
hListboxE = handles.listbox_experiment;
hFig = ancestor(hListboxE, 'figure');
hAx = handles.axes_img;

if nargin < 2
    %display selected subject in experiment 
    
    selectedVal = get(hListboxE, 'Value');
    
    if isempty(selectedVal)
        return
    end
    
    expTable = getappdata(hListboxE, 'table');
    expRow = expTable(selectedVal, :);
end

tableContainer = getappdata(hFig, 'tableContainer');

if isempty(tableContainer)
    return
end

%Display plot

megaTable_Sag = getExperimentRow(tableContainer, expRow.ExperimentID(1));
hListboxA = handles.listbox_analysis;
setappdata(hListboxA, 'table', megaTable_Sag);

if isempty(megaTable_Sag)
    %Display the image only
    hImg = displayOriginalImage(hAx, expRow);
    
else
    %Display the analysis
    selectedAnalysis = get(hListboxA,'Value');
    hImg = displaySagittalAnalysis(hAx, megaTable_Sag(selectedAnalysis,:));
end

%__________________________________________________________________________

function hImg = displayOriginalImage(hAx, thisRow)
%This function displays an image on an image given the row from the
%experiment table hImg is empty if it was unsuccessful

handles = guidata(ancestor(hAx, 'figure'));

if size(thisRow, 1) > 1
    thisRow = thisRow(1,:);
end

[img, ~] = readImage(thisRow);

cla(hAx,'reset'); setappdata(hAx, 'hPoints',[]);
if isempty(img)
    title('Image Not Found');
    hImg = [];
    return
end

if strcmpi(get(handles.tool_stretch, 'State'),'on')
    hImg = imagesc(img,'Parent', hAx);
else
    hImg = imshow(img,'Parent', hAx);
end


set(hAx, 'XTick',[], 'YTick',[]);
 
set(handles.label_imgTitle, 'String',...
    strrep(char(thisRow.ExperimentID),'_','  '));

setappdata(hAx, 'hImg', hImg);
setappdata(hAx, 'img', img);

function hImg = displaySagittalAnalysis(hAx, megaTable)
%megaTable --> subject, experiment, and analysis data

if size(megaTable, 1) > 1
    megaTable = megaTable(1,:);
end

%1. Read the image
[img, ~] = readImage(megaTable);

cla(hAx,'reset'); setappdata(hAx, 'hPoints',[]);
if isempty(img)
    title('Image Not Found');
    hImg = [];
    return
end

%2. Rotate the image as necessary
nRotations = megaTable.Param_InitialRotation;

if nRotations > 0
    %rotate the image CCW
    imgR = imrotate(img, nRotations*90);
else
    imgR = img;
end

%3. Create Image Display Scale
[ySize, xSize] = size(imgR(:,:,1));
imgRes = megaTable.Param_ImgRes;
shiftAmount = megaTable.ImageShift;
xScale = (1:xSize)/imgRes(1) - shiftAmount(1); 
yScale = (1:ySize)/imgRes(1) - shiftAmount(2);


%4. Display Centered Image
hImg = imagesc(xScale , yScale ,...
    imgR,'Parent', hAx);


%5. Show Segmentation and Fits
fitCell = megaTable.Fits; fitTable = fitCell{1};

%Determine Fit Type
hFig = ancestor(hAx, 'figure');
handles = guidata(hFig);
allStr = get(handles.popup_plotFitType,'String');
fitType = allStr{get(handles.popup_plotFitType, 'Value')};

if ~any(fitTable.FitName == fitType)
    fitType = 'Conic'; %default back to conic
end

% plotFit(hAx, fitTable, 'Cosine', 'Original', 'ry'); %fitType, inputDesc, colors
plotFit(hAx, fitTable, fitType, 'Original', 'bm'); %fitType, inputDesc, colors

%6. Zoom to Roi
zoomToRoi(hAx, megaTable);
 

%7. Tilt axis
tilt_deg = megaTable.ImageTilt * (180/pi);
[az, el] = view(hAx);
view(hAx, az+tilt_deg,el);

setappdata(hAx,'hImg', hImg);
setappdata(hAx, 'img', imgR);

%8. Title String
handles = guidata(ancestor(hAx, 'figure'));
set(handles.label_imgTitle, 'String',...
    sprintf('%s {%s Fit}',strrep(char(megaTable.ExperimentID),'_','  '), fitType));

%..........................................................................

function [img, imgLocation] = readImage(imgTable,iRow)
%This function tries to find the image and load it
%INPUT
%Single row of the Experiment table

if nargin < 2
    iRow = 1;
end

origLoc = char(imgTable.Img_originalLocation(iRow));
img = []; %initialize
imgLocation =[];

if exist(origLoc, 'file')
    img = imread(origLoc);
    imgLocation = origLoc;
else
    localDir = char(imgTable.Img_localDirectory(iRow));
    
    if isempty(imgTable.Img_filename{iRow})
        [~, imgFname, ext] = fileparts(origLoc);
        imgFname = [imgFname, ext];
    else
        imgFname = char(imgTable.Img_filename(iRow));
    end
    localLocation = fullfile(localDir, imgFname);
    if exist(localLocation, 'file')
       img = imread(localLocation);
       imgLocation = localLocation;
    else
        warning('Cannot find the image file!');
    end
    
end

function megaTable = getExperimentRow(tableContainer, ExperimentID)
%MegaTable --> subject, experiment, and sagittal table

%Initialize Output
megaTable = [];

if ~isfield(tableContainer, 'Sagittal') 
    return
end

if iscell(ExperimentID)
    ExperimentID = cell2mat(ExperimentID);
end

subjectTable = tableContainer.Subject;
experimentTable = tableContainer.Experiment;
sagittalTable = tableContainer.Sagittal;


%Get the subject row from each table
isThisRow = strcmpi(sagittalTable.ExperimentID, ExperimentID);

if ~any(isThisRow)
    return
end

thisSagID = sagittalTable(isThisRow, :);

thisRow = join(join(thisSagID, subjectTable),experimentTable);

%by default, display the latest (you can also have them select from listbox)
thisRow = sortrows(thisRow, {'CreatedDate'},'descend');

megaTable = thisRow;

function hP = plotFit(hAx, fitTable, fitType, inputDescription, colors)
%This function plots the points denoted by fitType and inputDescription
%INPUT:
    %hAx - handle of the axes to plot on
    %fitTable - Fit table from semiAutoSegment_shadowgraph output (field: fits)
    %fitType - categorical string {'Cosine', 'Conic'}, case sensitive
    %inputDescription - string {'Original Segmentation Points','Corrected
        %Segmentation Points'} can be partial string and is case
        %insensitive


if nargin < 5
    colors = ('rb')';
else
    if ischar(colors)
        colors = force1D(colors,2);
    end
end


if nargin < 4
    inputDescription = 'original'; %should always be present
end

if nargin < 3 || isempty(fitType)
    fitType = 'Conic'; %'Cosine';
end


%make sure fit type is title case
fitType = [upper(fitType(1)), lower(fitType(2:end))];

%find the row for conic fit of corrected segmentation points

isThisFit = fitTable.FitName == fitType &...
    ~cellfun(@isempty,(strfind(lower(fitTable.InputDescription), lower(inputDescription))));

if ~any(isThisFit)
    return
end

fitRow = fitTable(isThisFit,:);
fitRow = fitRow(1,:); %only take the first row

xyPoints = fitRow.InputPoints; xyPoints = xyPoints{1};
xyFit = fitRow.FitPoints; xyFit = xyFit{1};


hold(hAx, 'on');
hP_new = [];

for j = 1:numel(xyPoints)
    %For each surface
    jPoints = xyPoints{j}; jFit = xyFit{j};
    jX_points = jPoints(:,1); jY_points = jPoints(:,2);
    jX_fit = jFit(:,1); jY_fit = jFit(:,2);
    
    hP_new(end+1) = plot(hAx, jX_points, jY_points,'x','Color', colors(1,:));
    hP_new(end+1) = plot(hAx, jX_fit, jY_fit,'LineWidth',2,'Color', colors(2,:));
    
    
end

hold(hAx, 'off');

hP = getappdata(hAx, 'hPoints');

if ~isempty(hP)
    hP(~ishandle(hP)) = [];
    hP = cat(1,hP_new);
else 
    hP = hP_new;
end

setappdata(hAx, 'hPoints', hP);

function zoomToRoi(hAx, megaTable)
%The table must have the fields: zoomRegion, imgRoi, ImgShift, ImgRes
%imgR is an image the same size as the rotated image

imgRes = megaTable.Param_ImgRes; imgRes = imgRes(1);
shiftAmount = megaTable.ImageShift;


zoomRegion = megaTable.Param_ZoomRegion;


%zoom as necessary
if isempty(zoomRegion) 
    
    imgRoi = megaTable.Param_ImgROI; imgRoi = imgRoi{1};
    imgSize = megaTable.ImageSize;
    
    %make it a mask
    if size(imgRoi, 2) == 1
        imgRoi_Idx = imgRoi;
        imgRoi = false(imgSize);
        imgRoi(imgRoi_Idx) = true;
    end
    
    if ~all(imgRoi(:))
    
        %use the imgROI to determine zoom region
        props = regionprops(imgRoi,'BoundingBox','Area'); tol = 100; %pxls
        [~, idx_sorted] = sort(props.Area,'descend');
        props = props(idx_sorted);
        bBox = props(1).BoundingBox;
        xlim = [bBox(1)-tol, bBox(1)+bBox(3)+tol];
        ylim = [bBox(2)-tol, bBox(2)+bBox(4)+tol];
        
        zoomRegion = {xlim, ylim};
        
    end
end

if isempty(zoomRegion)
    %if it's still empty, return
    return
end

zoomRegion_mm = zoomRegion;
zoomRegion_mm{1} = (zoomRegion{1}./imgRes) - shiftAmount(1);
zoomRegion_mm{2} = (zoomRegion{2}./imgRes) - shiftAmount(2);

zoom(ancestor(hAx,'figure'), 'reset');
set(hAx, {'XLim', 'YLim'}, zoomRegion_mm);
axis(hAx,'manual'); %to make sure it stays zoomed
drawnow;

function consolidateImages(handles)
%this function consolidates all image files to a single folder using
%filename: experimentID

hFig = ancestor(handles.axes_img,'figure');


%.........................................................................
%make sure filenames are unique
%.........................................................................
tableContainer = getappdata(hFig, 'tableContainer');

if isempty(tableContainer) || ~isfield(tableContainer,'Experiment')
    return
end

expTable = tableContainer.Experiment;

[~,ia,~] = unique(expTable.Img_originalLocation, 'stable');
expTable = expTable(ia,:);

%.........................................................................
%Get the new file directory
%.........................................................................

localDir = get(handles.text_localDir, 'String');

%see if any of the files have a localDirectory
if isempty(localDir)
    priorLocalDirs = unique(expTable.Img_localDirectory);
    
    isValidDir = cellfun(@(iDir) exist(iDir, 'dir'), priorLocalDirs);
    
    if ~any(isValidDir)
        errordlg('Please set the local image directory in the settings panel.',...
            'No Image Directory');
        return
    end
    
    priorLocalDirs = priorLocalDirs(isValidDir);
    
    [selection, isOk] = listdlg('PromptString','Please select the local directory.',...
        'SelectionMode','single', 'Name','Select Directory',...
        'ListSize', [450 100], 'ListString',priorLocalDirs);
    
    if ~isOk
        return
    end
    
    localDir = priorLocalDirs(selection);
    
end


if isempty(localDir) || ~exist(localDir, 'dir')
    errordlg('Please set the local image directory in the settings panel.',...
        'No Directory Found');   
    return
end



%.........................................................................
%start copy process
%.........................................................................


isFileHere = logical(cellfun(@exist, expTable.Img_originalLocation));

if ~any(isFileHere)
    errordlg('Image files not found on local computer');
    return
end

srcLocation = expTable.Img_originalLocation(isFileHere);


newExt = '.jpg';
newFname = cellfun(@(iName) [iName, newExt],...
    expTable.ExperimentID(isFileHere), 'UniformOutput',false);

newLocation = cellfun(@(iName) fullfile(localDir, iName),...
    newFname, 'UniformOutput', false);

nImgs = numel(newLocation);

set(hFig, 'Pointer','watch');
%Copy the files to the new location
[status,msg, msgID] = arrayfun(@(idx) copyfile(srcLocation{idx}, newLocation{idx}),...
    1:nImgs,'UniformOutput', false);
set(hFig, 'Pointer','arrow');

if iscell(status)
    status = cell2mat(status);
end

isFileCopied = isFileHere;
isFileCopied(isFileHere) = force1D(status,2); 

%Prompt User

if any(isFileCopied)
    fprintf(1,'\n\n%d Files Copied Successfully\n',numel(find(status)));
    cellfun(@disp, newLocation(status));
    disp(' ');
else
    disp(msg{1});
    return
end

%.........................................................................
%update database
%.........................................................................

expTable{isFileCopied, 'Img_localDirectory'} = cellstr(localDir);
expTable{isFileCopied, 'Img_filename'} = force1D(newFname(isFileCopied),2);

tableContainer.Experiment = expTable;

setappdata(hFig, 'tableContainer', tableContainer);

helpdlg(sprintf('%d images copied to folder:\n%s',...
    numel(find(status)),localDir));

% --------------------------------------------------------------------
% CREATE FUNCTIONS
% --------------------------------------------------------------------

function listbox_subject_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_subject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function listbox_experiment_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_experiment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function listbox_analysis_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_analysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_subjectID_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_subjectID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_notes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_notes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_experimentID_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_experimentID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_age_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_age (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%..........................................................................
%Settings Panel
%..........................................................................

function popup_xField_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_xField (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_yField_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_yField (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function listbox_tableColumns_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_tableColumns (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function listbox_fieldList1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_fieldList1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function listbox_fieldList2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_fieldList2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popup_tableList2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_tableList2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popup_tableList1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_tableList1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popup_plotFitType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_plotFitType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%..........................................................................
%Filter Panel
%..........................................................................

function popup_filterEye_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_filterEye (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popup_filterImgType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_filterImgType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popup_filterValidated_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_filterValidated (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popup_filterSegmented_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_filterSegmented (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popup_filterExcluded_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_filterExcluded (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popup_filterNotes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_filterNotes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%..........................................................................
%Additional Settings
%..........................................................................


function edit_analyzedBy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_analyzedBy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function text_localDir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_localDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --------------------------------------------------------------------


% --------------------------------------------------------------------



%..........................................................................
%Single Subject: Edit boxes callback
%..........................................................................

function fieldEdit_Callback(hObject, eventdata, handles)


prevValue = getappdata(hObject, 'prevValue');
changeFlag = true;
propertyName = [];
graphFlag = false;

switch get(hObject,'Tag')
    case {'edit_experimentID','edit_subjectID'}
        %don't allow editing of these properties
    case 'edit_notes'
        tableName = 'Experiment';
        fieldName = 'Notes';
        fieldValue = cellstr(get(hObject, 'String'));
 
        propertyName = 'String';
        
    case 'edit_age'
        
        tableName = 'Experiment';
        fieldName = 'Age';
        fieldValue = str2double(get(hObject, 'String'));
        
        if isempty(fieldValue)
            fieldValue = prevValue;
            changeFlag = false;
        end
        
        propertyName = 'String';
        
    case {'check_flagAnt','check_flagPos'}
        tableName = 'Experiment';
        fieldName = 'OutlierFlag';
        fieldValue = logical([get(handles.check_flagAnt, 'Value'),...
            get(handles.check_flagPos,'Value')]);
        
        valOut = get(hObject, 'Value');
        
        propertyName = 'Value';
        
        graphFlag = true;
end

if changeFlag
    isOK = setTableValue(handles, tableName, fieldName, fieldValue);
else
    isOK = false;
end

if ~isOK
    valOut = prevValue;
elseif ~exist('valOut', 'var')
    valOut = fieldValue; %This is for the checkboxes because both must be set even though one is edited
end


set(hObject, propertyName, valOut);

if graphFlag
    updateGraphPanel(handles);
    updateSubjectSummaryText(handles);
end

% --------------------------------------------------------------------
function tool_stretch_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to tool_stretch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hTool = hObject; 
hAx = handles.axes_img;
hImg = getappdata(hAx, 'hImg');

if isempty(hImg) || ~ishandle(hImg)
    return
end


switch get(hTool,'State')
    case 'on'
        axis(hAx, 'normal');
    case 'off'
        axis(hAx, 'image');
end


% --- Executes on button press in check_flagPos.
function check_flagPos_Callback(hObject, eventdata, handles)
% hObject    handle to check_flagPos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_flagPos


% --------------------------------------------------------------------




% --- Executes on button press in b_localDir.


% --------------------------------------------------------------------
function cmenu_listbox_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------


% --------------------------------------------------------------------
