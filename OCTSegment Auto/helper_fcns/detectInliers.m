function xyInliers = detectInliers(validPixelMask, xyFit, fitTol)
%fitTol in pixels


[yPts, xPts] = find(validPixelMask);

xFit = xyFit(:,1); yFit = xyFit(:,2);
yLow = yFit - fitTol/2;
yHigh = yFit + fitTol/2;

lowXY = [xFit, yLow];
highXY = [xFit, yHigh];
roiPoints = [lowXY; flipud(highXY)];
isInlier = inpolygon(xPts, yPts, roiPoints(:,1), roiPoints(:,2));

xInliers = force1D(xPts(isInlier)); yInliers = force1D(yPts(isInlier));

%for every inlier, go a-line by a-line and pick the closest one to the fit
xyInliers = [];
for iX = xFit'
    isThisX = xInliers == iX;
    iXIdx = xFit == iX;
    thisYFit = yFit(iXIdx);
    
    if ~any(isThisX)
        continue
    end
    
    thisX = xInliers(isThisX); thisY = yInliers(isThisX);
    [~,closestYIdx] = min(abs(thisY - thisYFit));
    xyInliers(end+1,:) = [iX, thisY(closestYIdx)];
    
end