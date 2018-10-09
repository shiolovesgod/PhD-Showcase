function [lensMask, surroundMask] = estimateShadowgraphROIs(img_raw, isPlot)
%This function uses the raw shadowgraph image (loaded using imread) as
%input and outputs the lensROI logical mask along with the surrounding
%region containing the sutures
%INPUT:
%   img_raw = imread(fileLocation); you can also input a variantImg (mxnx1)
%OUTPUT:
%   lensMask - logical image size of img_raw containing the lensMask
%               estimate
%   surroundMask - logical image size of img_raw containing the estimate of
%               surrounding region

if nargin < 2
    isPlot = false;
end


if size(img_raw,3) > 1
    imgR = img_raw(:,:,2); %default to the blue channel
else
    imgR = img_raw;
end


intThresh = round(graythresh(imgR)*255);
isBrightPixel = imgR > intThresh*1.5;

se = strel('disk',90); %se should be big enough to close the holes

%filter the bright pixels using something similiar to a median filter
%(but faster)
isBrightPixel_mod = imclose(isBrightPixel, se);

%get the connected components
cc = bwconncomp(isBrightPixel_mod);

%use only the cc with the most pixels
nPixels = cellfun(@numel,cc.PixelIdxList);
[~,idx_biggest] = max(nPixels);
mask_biggestCC = zeros(cc.ImageSize);
mask_biggestCC(cc.PixelIdxList{idx_biggest}) = true;

%find the hole within the mask
filled_biggestCC = imfill(mask_biggestCC);
hole_biggestCC = ~mask_biggestCC & filled_biggestCC;

surroundMask = mask_biggestCC; %surrounding region minus the lens
lensMask = hole_biggestCC;

if ~any(lensMask(:)) %it couldn't be found
    lensMask = [];
    return
end
%check to see if lesnMask has more than one cc if it does, use the largest
%one

cc_lens = bwconncomp(lensMask);

%cancel that (was trying to find sutures) 
% darkImg = imgR(:,:,2) < intThresh*1.2;
% iEdge = darkImg & surroundMask;
% iEdge2 = imopen(iEdge, strel('disk',3));
% props = regionprops(iEdge,'Eccentricity','PixelIdxList');
% isValid = [props.Eccentricity] > 0.8;



%keep only the largest roi for the lens
if cc_lens.NumObjects > 1
    nPixels_lens = cellfun(@numel,cc_lens.PixelIdxList);
    [~,idx] = max(nPixels_lens);
    lensMask = zeros(size(lensMask));
    lensMask(cc_lens.PixelIdxList{idx}) = true;
end


if isPlot
    %show it on the original image
    edge_lens = edge(lensMask); [lensY, lensX] = find(edge_lens);
    edge_suture = edge(filled_biggestCC); [sutureY, sutureX] = find(edge_suture);
    
    %PLOT
    figure(30); imagesc(imgR);
    hold on; plot(lensX, lensY, '.r'); plot(sutureX, sutureY, '.b');
    
end