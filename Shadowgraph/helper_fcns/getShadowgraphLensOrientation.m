function nRotations = getShadowgraphLensOrientation(lensMask)
%This function returns the orientation of the shadowgraph lens given the
%lens ROI mask
%INPUT: lensMask - lensROI output from "estimateShadowGraphROIs" function
%OUTPUT: nRotations - number of counterclockwise (ccw) rotations needed to
                %be posterior down (0 --> posterior down, 1 --> posterior
                %left, 2--> posterior up, 3 --> posterior right);

% cc = bwconncomp(lensMask);
props = regionprops(lensMask, 'Centroid','Extrema', 'Eccentricity','Orientation');

eccThresh = 0.2;

if props.Eccentricity < eccThresh
    nRotations = [];
    return %won't be conclusive
end

%otherwise use the orientation 
orientation = props.Orientation; 

[~,orient_idx] = min(abs(abs(orientation) - [0 90]));
isVertical = orient_idx == 2;

%determine the location of the centroid w.r.t the equator


extrema = props.Extrema;
[~, leftIdx] = min(extrema(:,1)); [~, rightIdx] = max(extrema(:,1)); 
[~, bottomIdx] = min(extrema(:,2)); [~, topIdx] = max(extrema(:,2));

if isVertical %vertical major axis (SEE PAGE 6 of SW notebook 3)
    idxCol = 1;
    equator = mean([extrema(bottomIdx, idxCol), extrema(topIdx, idxCol)]);
    centroid = round(props.Centroid(idxCol));
    bit1 = (centroid - equator) > 0;
else %isHorizontal
    idxCol = 2;
    equator = mean([extrema(leftIdx, idxCol), extrema(rightIdx, idxCol)]);
    centroid = round(props.Centroid(idxCol));
    
    bit1 = (equator - centroid) > 0;
end

bit2 = isVertical;
nRotations = bin2dec(sprintf('%d%d', [bit1, bit2]));


    