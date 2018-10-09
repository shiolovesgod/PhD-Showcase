function [xyOut, ptWeight] = getMagnetWeights(xyVals, imgIn, xyRes, xFitIn)
%
%INPUTS:  xyVals:   [x1, y1; x2,y2;...xn,yn]
        % imgIn:   image loaded from file without filtering
        % xyRes:    [xRes, yRes]
%OUTPUT:  
        % ptWeight: [x1y1_weight....xnyn_weight];


%%
%Initialize

[ySize, xSize] = size(imgIn); %UINPUT

isPlot = false;


isNormalize = true; isByALine = true; %normalize, if yes...by ALine?
isUseBrightest = false; 
isUseScaling = true;  scaleMin = 2; scaleMax = 6;
%%


%%
%==========================================================================
%Pick the center of the median filtered image
%==========================================================================
xPoints = xyVals(:,1); yPoints = xyVals(:,2);
xRes = xyRes(1); yRes = xyRes(2);

%do a quad fit, do a 1mm fit create a tolerance range
tolRange = 0.25*yRes;  %250um or 0.25mm (generous)
pFit = polyfit(xPoints, yPoints,2);
xFit = 1:xSize;
yFit = polyval(pFit, xFit);
yLow = yFit - tolRange/2;
yHigh = yFit + tolRange/2;

%get only the x with the first y
xyPoints_firstX = sortrows(xyVals,[1,2]); [~,ia,~] = unique(xyPoints_firstX(:,1), 'first');
xyPoints_firstX = xyVals(ia,:);

if isPlot
    figure(289); hAxMain = gca; imagesc(imgIn); hold on; plot(xyPoints_firstX(:,1), xyPoints_firstX(:,2),'xr'); colormap gray
    hold on; plot(xFit, yHigh, 'b'); plot(xFit, yLow, 'b');
    yIdx = 1:size(imgIn,1);
end


%%
%find bright a-line and min peak
lowXY = [xFit', yLow'];
highXY = [xFit', yHigh'];
roiPoints = [lowXY; flipud(highXY)];
roiMask = poly2mask(roiPoints(:,1), roiPoints(:,2), ySize, xSize);
roiMask(:,1) = circshift(roiMask(:,2),[1,0]);


%%
%NORMALIZE IMAGE
%Get sample from top and bottom regions
normWidth = 10; %10 pixels at top and bottom
yNormTop = yLow +normWidth;
yNormBottom = yHigh - normWidth;

topXY = [xFit', yNormTop'];
roiTOP = [topXY; flipud(lowXY)];
roiMaskTOP = poly2mask(roiTOP(:,1), roiTOP(:,2), ySize, xSize);
roiMaskTOP(:,1) = circshift(roiMaskTOP(:,2),[1,0]);

bottomXY = [xFit', yNormBottom'];
roiBOTTOM = [bottomXY; flipud(highXY)];
roiMaskBOTTOM = poly2mask(roiBOTTOM(:,1), roiBOTTOM(:,2), ySize, xSize);
roiMaskBOTTOM(:,1) = circshift(roiMaskBOTTOM(:,2),[1,0]);

%divide roi into above and below line
roiAbove = [lowXY; flipud([xFit', yFit'-1])];%flipud([xFit',ones(size(yFit))'])];% 
roiMaskTopHalf = poly2mask(roiAbove(:,1), roiAbove(:,2), ySize, xSize);
roiMaskTopHalf(:,1) = circshift(roiMaskTopHalf(:,2),[1,0]);

roiBelow = [[xFit', yFit'+1]; flipud(highXY)];%flipud([xFit', ones(size(yFit))'+ySize])];%
roiMaskBottomHalf = poly2mask(roiBelow(:,1), roiBelow(:,2), ySize, xSize);
roiMaskBottomHalf(:,1) = circshift(roiMaskBottomHalf(:,2),[1,0]);

%get intensity of sample region
topSample = double(roiMaskTOP).*double(imgIn);
% topSample(~roiMaskTOP) = nan;

bottomSample = double(roiMaskBOTTOM).*double(imgIn);
% bottomSample(~roiMaskBOTTOM) = nan;
%%
meanImgFilter = ones(size(imgIn))+mean(imgIn(:));

if isByALine
    topMeanVal = sum(topSample,1)./sum(roiMaskTOP,1);
    bottomMeanVal = sum(bottomSample,1)./sum(roiMaskBOTTOM,1);
    
    topMean = repmat(topMeanVal,[ySize,1]);
    
    if any(isnan(topMean(:)))
        %replace nans with a mean value 
        topMean(isnan(topMean)) = mean(topMean(:),1,'omitnan');
    end
    
    bottomMean = repmat(bottomMeanVal,[ySize,1]);
    if any(isnan(bottomMean(:)))
        bottomMean(isnan(bottomMean)) = mean(bottomMean(:),1,'omitnan');
    end
    
else
    topMeanVal = mean(topSample(:));
    bottomMeanVal = mean(bottomSample(:));
    
    topMean = zeros(size(imgIn)) + topMeanVal;
    bottomMean =  zeros(size(imgIn)) + bottomMeanVal;
end

meanImgFilter(roiMaskTopHalf) = topMean(roiMaskTopHalf);
meanImgFilter(roiMaskBottomHalf) = bottomMean(roiMaskBottomHalf);

%%
%CREATE ROI
%determine window size
dWindow = round(range(xFitIn)/xRes); %6; %mm should i continue to use 6mm?
dWindow_px = dWindow*xyRes(1);
cntr = round(mean(xFitIn)); %xSize/2;  [~, cntr] = min(abs(diff(yFit) - 0)); %
leftWin = round(cntr-dWindow_px/2); 
rightWin = round(cntr+dWindow_px/2);
[~,leftWindIdx]=min(abs(xFit - leftWin));
[~,rightWindIdx]=min(abs(xFit - rightWin));

% roiMask2 = roiMask; roiMask2(:,180:190)=0;
roiMask2 = roiMask; roiMask2(:,[1:leftWindIdx,rightWindIdx:xSize]) = 0;

%..........................................................................
%USE BRIGHTEST PIXELS from raw OCU as anchors?
%..........................................................................

if isUseBrightest
    bases =  [10 2 1.45];
    isBrightest = false(size(imgIn));
    for iBase = bases
        Iraw = iBase.^(double(imgIn)/10);
        minThresh = mean(Iraw(:)) + 4*std(Iraw(:));
        minThreshMat = mean(Iraw,1) + 4*std(Iraw,1);
        Iraw_mask = Iraw > minThresh;
        
        %Determine whether or not that A-Line should be considered
        %Only if that A-Line has no ROI already
        isInclude = ~repmat(logical(sum(isBrightest,1)),[ySize,1]);
        
        isBrightest = isInclude & Iraw_mask | isBrightest;
        %     keyboard;
        
    end
    
    brightestPixelsMask = isBrightest & roiMask2;
    
    if any(brightestPixelsMask(:))
        %should I be using this method?
        
        
    end
    
    
    hasNoPixelsMask = false(size(imgIn)); %all alines that dont'have any pixels
    hasNoPixels = ~logical(sum(isBrightest & roiMask,1));
    hasNoPixelsMask(:,hasNoPixels)= true;
    brightestPixelsValidArea = isBrightest | hasNoPixelsMask;
    
    roiMask2 = roiMask2 & brightestPixelsValidArea;
else
    brightestPixelsMask = false;    
end
%..........................................................................

imgCrop = (double(imgIn).*roiMask2);

if isNormalize
    imgCrop = imgCrop./meanImgFilter;
end

%%
if isUseScaling    
    roiImg = double(imgIn).*roiMask;
    roiImgStat = roiImg; roiImgStat(roiImg==0) = nan;
    iStd = std(roiImgStat,[],1,'omitnan');
    iMean = mean(roiImgStat,1,'omitnan');
    iMax = max(roiImgStat, [],1);
    iMedian = median(roiImgStat,1,'omitnan');
   
    iMaxSmoothed = medfilt1(iMax,10);
    saturationRatio = numel(find(iMax == 255))/400;
    
    if(saturationRatio > 0.03)
        fitFcn = @(b,x) b(1)+b(2)*exp(-(abs(x-b(3))/b(4)).^(b(5)/2));
        startVals = [mean(iMax), max(iMax),200, 20,3];
        weights = ones(numel(iMax),1)-0.5; weights(150:250) = 10;
        
        try
            mdl = fitnlm(1:numel(iMax),iMaxSmoothed', fitFcn, startVals,...
                'Weights',weights);
        catch
            %no weights
            mdl = fitnlm(1:numel(iMax),iMaxSmoothed', fitFcn, startVals);
        end
        
        %rescale model
        coeff = mdl.Coefficients{:,1};
        coeff(1) = 1; coeff(2)=2; %(between 1 & 2)?
        normalizedFit = fitFcn(coeff,1:numel(iMax));
        normalizedFitMatrix = repmat(normalizedFit,[size(imgIn,1),1]);
        
      else
        %use a constant scale throught the image
        normalizedFitMatrix = repmat(mean([scaleMin, scaleMax]), size(imgIn));
    end
    
    
end


%create weights from the raw img
%turn everything in that region into a point;
%create points
cropIdx = find(logical(imgCrop));
[cropYOrig, cropXOrig] = ind2sub(size(imgCrop),cropIdx);


%weight --> combination of intensity and distance from orig fit?
cropInt = mat2gray(imgCrop);% imgCrop;%
% minInt = min(cropInt(:)); meanInt = mean(cropInt(:)); stdInt = std(cropInt(:));
% intWeight = mat2gray(cropInt, [0, mean(cropInt(:)) + 11*std(cropInt(:))]);
% figure; imagesc(cropInt); colormap gray;

cropGrad =  mat2gray(conv2(double(imgIn), [1;-1],'same').*roiMask2); %mat2gray(
% figure; imagesc(cropGrad); colormap hot;
gradEffect = mat2gray(cropInt.*cropGrad);

%Linear indices of orignal fit
fitIdx = sub2ind(size(imgIn), round(yFit), xFit);

%find the difference between the new possible pixels and the original fit
%line 
[cropIdxGrid, fitIdxGrid] = meshgrid(cropIdx',fitIdx);

distFromFit = min(abs(cropIdxGrid-fitIdxGrid),[],1)';

distFromFit_norm = mat2gray(distFromFit);

%Weight Function


cropIntMat = cropInt(cropIdx);
% gradEffectMat = gradEffect(cropIdx);
% ptWeightAll = exp(4*cropIntMat)+ exp(1*gradEffectMat) - exp(4*distFromFit_norm);
% ptWeightAll = exp(3*gradEffectMat); %- exp(2*distFromFit_norm);
% ptWeightAll = exp(3*cropIntMat) - exp(2*distFromFit_norm);


if isUseScaling
    ptWeightAll = exp(normalizedFitMatrix(cropIdx).*cropIntMat);
else
    ptWeightAll = exp(4*cropIntMat);%10.^(cropIntMat);%cropIntMat;%cropIntMat;%gradEffectMat;%cropIntMat% %10.^cropIntMat;
end
ptWeightAll = ptWeightAll + abs(min(ptWeightAll))+1e-10;
% ptWeight = distFromFit_norm;

if any(brightestPixelsMask(:))
    increaseAmount = mean(ptWeightAll) + 3*std(ptWeightAll); %"anchor" pixels are increased by this amount
    brightestPxlIdx = find(logical(brightestPixelsMask));
    [bPxlMat, allIdxMat] = meshgrid(brightestPxlIdx, cropIdx);
    isHere = bPxlMat == allIdxMat;
    brightestPxl_relativeIdxLogical = logical(sum(isHere,2)); %where its located in cropIdx
    ptWeightAll(brightestPxl_relativeIdxLogical) = ptWeightAll(brightestPxl_relativeIdxLogical) + increaseAmount;
end


% brightestPxlIdx = find(logical(brightestPixelsMask));
%%
%REMOVE THOSE BELOW 2*STD? 
%=========================================================================
cropX = cropXOrig; cropY = cropYOrig;
ptWeight = ptWeightAll;

%=========================================================================

if isPlot
    cropIdx2 = cropIdx;
    imgWeight = zeros(size(imgCrop));
    imgWeight(cropIdx2) = ptWeight;
    figure(35); imagesc(imgWeight); %hold on;plot(xyVals(:,1), xyVals(:,2),'.g');
    colorbar; colormap hot
end

% cropX_mm = cropX/xyRes(1); cropY_mm = cropY/xyRes(2);

xyOut = [cropX, cropY];