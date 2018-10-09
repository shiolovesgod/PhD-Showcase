function imgOut = generateGenericImg(imgType)

imgOut = [];

if ischar(imgType)
    imgType = lower(imgType);
end

switch imgType 
    case {1,'cornea'}
        fname = 'corneaImg.tiff';
        fStr = 'Cornea';
    case {2,'lens'}
        fname = 'lensImg.tiff';
        fStr = 'Lens';
    case {3, 'retina'}
        fname = 'retinaImg.tiff';
        fStr = 'Retina';
end

[mFilePath, ~, ~]= fileparts(mfilename('fullpath'));

if isempty(mFilePath)
   mFilePath = strtok((cd),'\');
end

imgDir = fullfile(mFilePath, 'SegInit');

dirFlag = ~isdir(imgDir);

if dirFlag
    fileFlag = true;
else
    fileFlag = ~exist(fullfile(imgDir, fname),'file');
end


if fileFlag
    uiwait(msgbox(sprintf(['A sample %s Image was not found. ',...
        'Please select an image that will serve as a prototype for all %s Images'], fStr, fStr),...
        'Select Prototype Image','help','modal'));
    
    %create it and store it
    [fnameLoaded, pnameLoaded] = uigetfile('*.tiff', 'Select Sample Image','MultiSelect','off');
    
    if ~fnameLoaded
        return
    end
    
    imgLoaded = imread(fullfile(pnameLoaded, fnameLoaded));
    hFig = figure; imagesc(imgLoaded); colormap gray
    title(sprintf('Prototype for %s Image', fStr));
    
    uInput = questdlg('Is this the prototype you would like to use for segmentation?',...
        'Accept Prototype?','Yes','No','Yes');
    
    if ishandle(hFig) 
        delete(hFig); 
    end
    
    if ~strcmpi(uInput, 'yes')
        return
    end
    
    if dirFlag
        mkdir(imgDir)
    end

    [status, msg] = copyfile(fullfile(pnameLoaded, fnameLoaded), fullfile(imgDir, fname));
    
    if ~status
        errordlg(sprintf('Problem copying file to the new folder:\n %s', msg),...
            'Write Error','modal');
        return
    end
    
    imgOut = imgLoaded;
    
else
    %just load it
    imgOut = imread(fullfile(imgDir, fname));
end

    