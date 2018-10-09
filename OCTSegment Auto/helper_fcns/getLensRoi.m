function lensRoiMask = getLensRoi(img, lensCoords, roiTol)
%creates an ROI around lens points +-roiTol on top and bottom and 1/8
%roiTol left and right
%INPUTS: 
        %img - only needed for the size
        %lensCoords - [x1,y1;x2;y2]
        %roiTol - scalar with tolerance in pixels

mc = polyfit(lensCoords(:,1), lensCoords(:,2),1);
xFit = 1:size(img,2);%min(lensX)-5:max(lensX)+5;
yFit = polyval(mc, xFit);
yFit_top = yFit - roiTol; yFit_bottom = yFit + roiTol;

xLeft = max([1, lensCoords(1,1)-round(roiTol/8)]);
xRight = min([lensCoords(2,1)+round(roiTol/8), size(img,2)]);
lensRoiMask = poly2mask([xLeft, xRight, xRight, xLeft, xLeft],...
    [yFit_top(1), yFit_top(end), yFit_bottom(end), yFit_bottom(1), yFit_top(1)],...
    size(img,1), size(img,2));
