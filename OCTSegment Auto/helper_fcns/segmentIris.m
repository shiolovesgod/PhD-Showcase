function [lensCoords, lensFlag, dataRangeMask, lensMaskF_thick] = segmentIris(imgIn, imgType, isPlot)
%This function locates the y indices of data ranges
%INPUTS: imgIn - plain filtered image (without mean removed from each row)

%Initialize
lensCoords = [];
lensFlag = true;
imgIn = double(imgIn); %just encase class is uint8

img = imgIn - repmat(mean(imgIn,2),[1, size(imgIn,2)]);
diff_colSum = sum(imgIn - img,2);


%Do a  median filter on iDiff and subtact the mean
diffF = medfilt1(diff_colSum,100); diffF = diffF-mean(diffF);

%Everything below the baseline becomes zero
diffF_baseline = diffF;
diffF_baseline(diffF < 0) = 0;%mean(iDiffF); %change the baseline

%Use the new mean as a threshold to create a step function
stepDiffF = zeros(size(diffF_baseline));
meanThresh = mean(diffF_baseline(diffF_baseline > 0));

stepBaseline = mean(diffF_baseline);

if imgType == 2 && false
    stepBaseline = 0.6*stepBaseline; %MODIFIED 12/3/2015
end

stepDiffF(diffF_baseline > stepBaseline) = 1;
stepEdges = diff(stepDiffF);
positiveEdgeIdx = find(stepEdges>0);
negativeEdgeIdx = find(stepEdges<0);
[peakVal,dataPeakIdx, peakWidth] = findpeaks(diffF,'MinPeakHeight',...
    mean(diffF)-0.5*std(diff_colSum),...+std(iDiff_colSum)
    'MinPeakDistance',100,...
    'Annotate','extents'); %


%use the positive and negative edges to find the dataRange of the image
%if only one edge is found, threshold is too  high. Use find peaks to get
%the other data limit
if numel(positiveEdgeIdx) < 2
    %Always take the first and last peaks for cornea, first and biggest peaks
    %for lens 
    
    %which edge did I find? lens or cornea? does the top or bottom half
    %have more data? if there is more data in the top half, i have the
    %bottom
    
    %use the edge furthest away from the current peak
    [~,edgeIdxNo] = max(abs(positiveEdgeIdx-dataPeakIdx));
    
%     if imgType == 2 || sum(sum(imgIn(1:positiveEdgeIdx,:))) < sum(sum(imgIn(negativeEdgeIdx:end, :)))
%         edgeIdxNo = 1;
%     else
%         edgeIdxNo = numel(dataPeakIdx);
%     end
    
    if abs(dataPeakIdx(edgeIdxNo) - positiveEdgeIdx) > 100
        %exaggerate peak width a little
        positiveEdgeIdx = [round(dataPeakIdx(edgeIdxNo) - peakWidth(edgeIdxNo));...
            positiveEdgeIdx];
        negativeEdgeIdx = [round(dataPeakIdx(edgeIdxNo) + peakWidth(edgeIdxNo));...
            negativeEdgeIdx];
    else
        %Couldn't find the first/surface surface!
        %TELL ME THAT THERE'S A PROBLEM
        disp('Linda, listen. Honey, there is a problem with an edge')
    end
elseif numel(positiveEdgeIdx) > 2
    
    positiveEdgeIdx = positiveEdgeIdx([1, end]);
    negativeEdgeIdx = negativeEdgeIdx([1, end]);
    
    if peakVal(1) > meanThresh*.25 && imgType == 2 %MODIFIED 12/3/2015
        positiveEdgeIdx(1) = dataPeakIdx(1);    
    end
else
     if peakVal(1) > meanThresh*.25 && imgType == 2 %MODIFIED 12/3/2015
        positiveEdgeIdx(1) = dataPeakIdx(1);    
    end
end


%replace edges below zero or above image size
positiveEdgeIdx(positiveEdgeIdx < 1) = 1;
negativeEdgeIdx(negativeEdgeIdx > size(imgIn,1)) = size(imgIn, 1);

%where the data of the image is located %THIS IS NOT EXACT (an estimate)
dataRange = sortrows([positiveEdgeIdx, negativeEdgeIdx],1);


%create a mask of where data is located (RGB image, B channel is empty)
[surf1DataRangeMask, surf2DataRangeMask, surf3DataRangeMask,safeLensMask] = deal(zeros(size(img)));
surf1DataRangeMask(dataRange(1,1):dataRange(1,2),:) = 1;
surf2DataRangeMask(dataRange(end,1):dataRange(end,2),:) = 1;
dataRangeMask = cat(3, surf1DataRangeMask,surf2DataRangeMask,surf3DataRangeMask);
safeLensMask(dataRange(1,end):end,:) = 1; %SURELY includes the lens

% figure(22); imagesc(imgIn), colormap gray;
% hold on; plot(repmat([0 400],[numel(dataRange),1])',...
%     repmat(dataRange(:),[1,2])','Color','r');
% return;
 

%use data range mask to find lens ROI

%ESTIMATE LENS WIDTH
lensFlag = false;


%The lens is characterized by a bright reflectivity
lensThresh = multithresh(img, 3);
lensThresh = lensThresh(end) + 10;%lensThresh = max([150, lensThresh(end)]);
lensMask = img > lensThresh & safeLensMask; %dataRangeMask(:,:,2); %use the lensMask
lensMaskF = medfilt2(lensMask,[5,5]);
lensMaskF = imfill(lensMaskF, 'holes'); %fill holes
lensMaskF = bwmorph(bwmorph(lensMaskF, 'thicken'),'thicken'); %thicken 
lensMaskF = bwmorph(lensMaskF, 'bridge');
lensMaskF_thick = bwdist(lensMaskF) < 8;
lensMaskF_thick = bwmorph(lensMaskF_thick,'bridge');


%Remove small connected components
cc = bwconncomp(lensMaskF_thick);
ccSize = cellfun(@numel,cc.PixelIdxList);

if cc.NumObjects > 2
    isTooSmall  = ccSize < 0.6*mean(ccSize);
else 
    isTooSmall = repmat(false, [cc.NumObjects,1]);
end
pixelList= regionprops(cc,'PixelList'); 
ccArea = regionprops(cc,'Area'); ccArea = cell2mat(struct2cell(ccArea));
pixelList = pixelList(~isTooSmall); 



if numel(pixelList) > 1
    pixelList = {pixelList.PixelList};
    
    %go back to the original mask...no padding
    for i = 1:numel(pixelList)
        iMask = zeros(size(imgIn));
        iPixels = pixelList{i};
        iMaskIdx = sub2ind(size(iMask), iPixels(:,2), iPixels(:,1));
        iMask(iMaskIdx) = 1;
        [iI, iJ]= find(iMask & lensMask);
        pixelList{i} = [iJ, iI];
    end


    
    %leftmost, topmost, rightmost, and bottom-most
    boundCoords = cellfun(@getGroupBounds, pixelList,'UniformOutput', false);
    boundCoords = cat(1, boundCoords{:});   
    
    %get leftmost pixels
    leftBounds = cat(1,boundCoords.left); %both xy coords
    rightBounds = cat(1,boundCoords.right);
    
    if imgType==2 
        %use the posterior lens to help narrow right answer
        antRoiImg = surf1DataRangeMask .* img;
        antRoiImg(antRoiImg < 1)= nan;
        antRoiThresh = mean(antRoiImg(:),1,'omitnan') ;
        antRoi = surf1DataRangeMask & img > antRoiThresh;% lensThresh-5;
        [~, antRoiJ] = find(antRoi);
        leftX = min(antRoiJ); rightX = max(antRoiJ);
        [~,rightGroupIdx] = min(abs(leftBounds(:,1)-rightX)); 
        [~,leftGroupIdx] = min(abs(rightBounds(:,1)-leftX));
    else
        [~,leftGroupIdx] = min(leftBounds(:,1)); %only x coords
        [~,rightGroupIdx] = max(rightBounds(:,1));
    end
    
    %xy and y value for left and right lens coords
    lensCoords = [boundCoords(leftGroupIdx).right;... 
        boundCoords(rightGroupIdx).left];
else
    lensCoords = [];
    
end

% 
% figure(23); plot(sum(img(dataRange(1,2):dataRange(2,1),:),1));
% figure(24); imagesc(lensMaskF_thick); %lensMaskF
% title(sprintf('Image: %d/%d', idxNo, numel(randImgNos)));
% pause(1); continue % %keyboard;


%     
% figure(23);clf; imagesc(imgIn), colormap gray;
% figure(24); clf; imagesc(lensRoi), colormap gray; 
% pause(0.1); continue


function boundOut = getGroupBounds(iGroup)
%leftmost, topmost, rightmost, and bottom-most

%get leftmost point
xGroupSort = sortrows(iGroup,[1,2]);
yGroupSort = sortrows(iGroup , [2,1]);
boundOut.left = xGroupSort(1,:);
boundOut.top = yGroupSort(1,:);
boundOut.right = xGroupSort(end,:); 
boundOut.bottom = yGroupSort(end,:);
