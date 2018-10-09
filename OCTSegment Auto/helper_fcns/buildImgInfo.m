function [imgInfo, outputFlag] = buildImgInfo(imgFilename)
%This function returns a structure with the image information given the
%full file path of an image:
%
%INPUT:
    %imgFilename (string) - contains full file path of image file
%OUTPUT:
    %imgInfo (structure) - information of the image containing the
    %following fields:
        %imgSize - [ySize, xSize] or [nRows, nCols]
        %xyRes - [xResolution, yResolution] in pixels per mm
        %scanWidth - (scalar) scanWidth in mm
        %scanDepth - (scalar) scanDepth in mm
	%
    %outputFlag (logical) - false if the imgInfo is accurate, true if it's
        %not

img = imread(imgFilename); 
[ySize, xSize] = size(img);

info = imfinfo(imgFilename);
outputFlag = false; %if something goes wrong give the user a flag
if isfield(info, 'ImageDescription')
    description = info.ImageDescription;
    
    resStr = regexp(description, '\d+\.\d+', 'match');
    res =str2double(resStr);
else
    %Get user to input scan width
    uInput = inputdlg('Please input the image scan width (mm).',...
        'Scan Width', 1, {'8'});
    
    %If it's empty, throw a flag
    if isempty(uInput)
        outputFlag = true;
        return
    end
    
    scanWidth = str2double(regexp(uInput, '\d+', 'match'));
    
    %If you can't understand it, throw a flag
    if isempty(scanWidth)
        outputFlag = true;
        errordlg(sprintf('User entered invalid scan width: %s\n Please input numbers in mm.',...
            [char(39),uInput,char(39)]),'Error','modal');
        return
    end
    
    %calculate the resolution
    res = [scanWidth/xSize, 10.43/ySize];
    
end

%Calculate xyResolution (circular, I know)
scanWidth = xSize/res(1); %[width, height]
scanDepth = 10.43;
xyRes = [xSize/scanWidth, ySize/scanDepth];

%Prepare Output
imgInfo = [];
imgInfo.imgSize = [ySize, xSize];
imgInfo.xyRes = xyRes;
imgInfo.scanWidth = scanWidth; 
imgInfo.scanDepth = scanDepth;
