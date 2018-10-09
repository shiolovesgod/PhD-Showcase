
function [allSegmentationPts, areaImg] = segmentFromCenter(img, roiBounds, boundary)
%This function finds the boundaries of an round shape given the center of
%the shape and the bounds of the ROI
%INPUTS
%img: image to be segmented (2D)
%roiBounds: the coordinates of the boundaries:
%   [leftCntrX, leftCntrY;
%   rightCntrX, rightCntrY;
%   topCntrX, topCntrY;
%   bottomCntrX, bottomCntrY]
%boundary: 'inner' or 'outer'


%Parse inputs
if nargin < 3
    boundary = 'outer';
end

if strcmpi(boundary, 'inner')
    occurence = 'first';
else
    occurence = 'last';
end

%We will use the x-axis only to find the center coordinate and to divide
%top and bottom halves
xAxis = roiBounds(1:2,:); yAxis = roiBounds(3:4,:);
lensCntrX = round(sum(roiBounds(1:2,1))/2); lensCntrY = round(sum(roiBounds(3:4,2))/2);
imgSize = size(img);
nCols = size(img,2);

%separate into thirds
nSegments = 7;
leftmostCoords = round(sum(xAxis,1)/nSegments); leftSegmentX = leftmostCoords(1);
rightSegmentX= leftSegmentX*(nSegments-1);

%..........................................................................
%DIVIDE THE IMAGE INTO L,R, TOP & BOTTOM HALVES ACCORDING TO THE LENS CNTR
%..........................................................................
idxMat = reshape(1:numel(img), size(img));

%Get the horizontal and vertical limits
leftBoundX = xAxis(1,1); rightBoundX = xAxis(2,1);
topBoundY = yAxis(1,2); bottomBoundY = yAxis(2,2);

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


