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
%     warndlg('Enlarge your ROI');
    validPoints =[];
    return
end

validPoints = double(sortrows(cat(1,allData{:}),1));


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

