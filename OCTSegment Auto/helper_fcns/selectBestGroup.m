
function [bestGroup, allGroups] = selectBestGroup(rawX, rawY, fType, isConvex, userData)
%This function selects the group of points that best fits the function
%..........................................................................
%INPUTS:
%..........................................................................
%rawX = 1D array containing x points
%rawY = 1D array containing y points
%fType = function handle to fitting function that outputs coefficients
%isConvex = should the polynomial be positive or negatively curved (valid
%           for quadratics
%userData = structure of user inputs containing y resolution (userData.xyRes(2))
%..........................................................................
%OUTPUTS
%..........................................................................
%bestGroup = structure containg statistics (x and y points, outliers etc)
%              of best fit
%allGroups = structure containting all groups created from original data


u_yRes = userData.xyRes(2);
dw_pxls = 0.15*u_yRes; %0.3mm

bestGroup = [];
allGroups = [];


%have a mincurvature
if userData.scanWidth > 10
    minA = 0.0005;
elseif userData.scanWidth > 7
    minA = 0.0005;
elseif userData.scanWidth < 7
    minA = 0.0001;
else
    minA = 0;
end


if isfield(userData, 'surfNo')
    switch userData.surfNo
        case 1
            idealR = 7.9;
        case 2
            idealR = 6.4; %ADJUST THIS ONCE YOU KNOW
        case 3
            idealR = 15;
        case 4
            idealR = 8;
        otherwise
            idealR = [];
    end
else
    switch userData.imgType
        case 1
            idealR = 7.9;
        case 2
            idealR = []; %ADJUST THIS ONCE YOU KNOW
        case 3
            idealR = [];
    end
end


%break data points into groups
originalGroups = classifyData(rawX,rawY, userData);
groups = originalGroups;

if isempty(groups)
    groups(1).data = [force1D(rawX),force1D(rawY)];
    groups(1).n = length(rawX);
    originalGroups(1).n = groups(1).n;
end

nGroups  = length(groups);
allGData = cat(1,groups(:).data);



if nGroups == 1
    cond2Thresh = 1.3;
else
    cond2Thresh = 1.7;
end

cond2 = numel(rawX)/sum([groups.n])> cond2Thresh ;
nALines = userData.scanWidth*userData.xyRes(1);



if nGroups > 1 && userData.scanWidth < 10 && userData.imgType == 1
    %make a group out of all the valid groups
    
    groups(end+1).data = allGData;
    groups(end).n = length(allGData(:,1));
elseif  cond2 || (nALines < 600 && userData.scanWidth > 8)
    %make a group out of all raw points
    groups(end+1).data = [rawX, rawY];
    groups(end).n = length(rawX);
end


groupSizeThresh = max([originalGroups.n]);

if cond2
    groupSizeThresh = max([groupSizeThresh, round(length(rawX)/3)]);
end

allGData = cat(1,groups(:).data);
allX = allGData(:,1); allY = allGData(:,2);


nGroups = length(groups);
isGroupLegit = false(1,numel(groups));
%determine which group gives best fit
for i = 1:nGroups
    data = groups(i).data;
    
    if userData.imgType == 4
        %is cornea
        dataIncluded = data;
        uniqueData = uniqueByCol(dataIncluded,2); %changed this to 2
    elseif userData.imgType ~= 1 %is not cornea
        %use only the central data
        nPoints = size(data,1); cntrIdx = round(size(data,1)/2);
        halfPtsIncl = round((0.95*nPoints)/2);  %central 70 changed from .7
        dataIncluded = data(max([1,cntrIdx-halfPtsIncl]): min([size(data,1),cntrIdx+halfPtsIncl]),:);
        uniqueData = dataIncluded;
        
    else
        %is cornea
        dataIncluded = data;
        uniqueData = uniqueByCol(dataIncluded,1); %changed this to 2
    end
    
    groupX = uniqueData(:,1); groupY = uniqueData(:,2);
    
    try
        [p, R, radiusFlag, inflectionFlag, outlierFlag] = robustFit(uniqueData,...
            dw_pxls, fType, isConvex, userData);
        p = p';
    catch ME
        [radiusFlag, inflectionFlag, outlierFlag] = deal(true);
        p = [];
    end
    
    if (radiusFlag||inflectionFlag||outlierFlag) %&& (nGroups == 1 || i==1)
        %size(uniqueData,1)*100;
        [p2,R2] = ransacNew(uniqueData, dw_pxls, fType,6000, 10, isConvex, userData); %dataIn, dw_pxls, nIter, nPoints
        if ~isempty(p2) && p2(1)
            p = p2;
            R = R2;
        else
            R = inf; %disqualify this Group
        end
        
        
    end
    
    if all(~p) %if p has an error with RANSAC
        p2 = feval(fType,groupX,groupY);
        
        if ~isempty(p2)
            p = p2;
        end
    end
    
    yFit = polyval(p,allX);
    
    %Find the outliers of all points
    upperBound = yFit + dw_pxls/2;
    lowerBound = yFit - dw_pxls/2;
    isOutlier = bitor(allY > upperBound, allY < lowerBound);
    nOutliers = sum(isOutlier);
    inlierX = allX(~isOutlier); inlierY = allY(~isOutlier);
    nInliers = numel(inlierX);
    inlierFit = polyval(p,inlierX);
    sseInlier = mean((inlierFit - inlierY).^2); %sum((inlierFit - inlierY).^2);
    
    %Find the number of outliers within the group
    groupFit = polyval(p, groupX);
    groupUpperBound = groupFit + dw_pxls/2;
    groupLowerBound = groupFit - dw_pxls/2;
    isGroupOutlier = bitor(groupY > groupUpperBound, groupY < groupLowerBound);
    nGroupOutliers = sum(isGroupOutlier);
    groupSize(i) = groups(i).n;
    groupsseInlier = mean((groupFit - groupY).^2);
    
    %Save findings to a structure
    iGroup.global.inliers = sortrows([force1D(inlierX), force1D(inlierY)]);
    iGroup.global.outliers = sortrows([force1D(allX(isOutlier)),force1D(allY(isOutlier))]);
    iGroup.global.nOutliers = nOutliers;
    iGroup.global.nInliers = nInliers;
    iGroup.global.bounds = [force1D(allX), force1D(upperBound), force1D(lowerBound)];
    iGroup.global.sse = sseInlier;
    iGroup.global.fit = [allX, yFit];
    
    iGroup.local.inliers = [force1D(groupX(~isGroupOutlier)),...
        force1D(groupY(~isGroupOutlier))];
    iGroup.local.outliers = [force1D(groupX(isGroupOutlier)),...
        force1D(groupY(isGroupOutlier))];
    iGroup.local.nOutliers = nGroupOutliers;
    iGroup.local.bounds = [force1D(groupX),...
        force1D(groupUpperBound), force1D(groupLowerBound)];
    iGroup.local.nPoints = groupSize(i);
    iGroup.local.sse = groupsseInlier;
    iGroup.local.fit = [groupX, groupFit];
    
    iGroup.p = p;
    iGroup.R = R;
    iGroup.yRange = [min(groupY), max(groupY)];
    iGroup.idx = i;
    
%     iGroup.p = p;
%     iGroup.yRange = [min(groupY), max(groupY)];
%     iGroup.idx = i;
    
    if nInliers > 0 && ~any(isnan(yFit))
        isGroupLegit(i) = true;
    end
    
    groupAnalysis(i) = iGroup;
    
    
    
end %end for loop


%COST FUNCTION
%use: how it fits overall data (A) AND how it fits groups (B)
validGroupsIdx = 1:numel(groupAnalysis);
sseInlier = arrayfun(@(x) groupAnalysis(x).global.sse,validGroupsIdx);
groupsseInlier= arrayfun(@(x) groupAnalysis(x).local.sse,validGroupsIdx);
nOutliers = arrayfun(@(x) groupAnalysis(x).global.nOutliers,validGroupsIdx);
nGroupOutliers = arrayfun(@(x) groupAnalysis(x).local.nOutliers,validGroupsIdx);
groupR = [groupAnalysis.R];

overallCost = sseInlier/max(sseInlier) + (nOutliers/numel(allX));
groupCost = (groupsseInlier/max(groupsseInlier))./groupSize + (nGroupOutliers./groupSize);

if max(groupCost) < 1 %modified 2/27/2014
    costFcn = 0.15*overallCost + exp(3*groupCost);
else
    costFcn = 0.35*overallCost + 3*groupCost;
end

if ~isempty(idealR)
    costFcn = costFcn + abs((abs(groupR) - idealR))/idealR;
end

costFcn_cell = num2cell(costFcn);

[groupAnalysis(:).costFcn] = deal(costFcn_cell{:});

%remove all polynomials with the wrong curve
%add a constraint for minimum curvature 2/27/2014

polys = cat(1,groupAnalysis.p);
As = polys(:,1);

fcnStack = dbstack;
if userData.imgType ~= 4 || strcmpi(fcnStack(3).name, 'fit2surfaces')
    if isConvex
        isAValid = As>0 & abs(As) > minA;
    else
        isAValid = As< 0 & abs(As) > minA;
    end
else
    isAValid = true(size(As));
end

validAs = As(isAValid);

%if the scan width is greater than 4mm, use a minimium group size as a
%constraint
if userData.imgType == 1
    minGroupSize = 0.5;
else
    minGroupSize = 0.2;
end

if userData.scanWidth > 5
    
    if userData.scanWidth > 12 && userData.imgType == 1
        groupSizeThresh = max([groupSizeThresh - 50,0]);
        
    elseif userData.scanWidth > 12 && userData.imgType == 2
        if isConvex
            groupSizeThresh = 0;
        else
            %do nothing
        end
    end
    isAboveMinSize = groupSize/min([groupSizeThresh, 250]) > minGroupSize;
    
else
    isAboveMinSize = groupSize/groupSizeThresh > 0; %all true
end



if ~any(validAs)
    allGroups  = groupAnalysis;
    return
else
    validGroups = groupAnalysis(isAValid & isAboveMinSize' & isGroupLegit');
    
end

if isempty(validGroups)
    validGroups = groupAnalysis(isAValid);
end

validCostFcn = [validGroups.costFcn];

%if the first group has comparably the lowest cost function use it
[~,idxMin] = min(validCostFcn);

if ~isempty(idealR)
    Rdiff = abs([validGroups.R] - idealR);
end

if (abs(min(validCostFcn)-validCostFcn(1)) < 0.3) && userData.imgType == 1
    gIdx = validGroups(1).idx;
elseif (abs(min(validCostFcn)-validCostFcn(1)) < 0.5) && (mean(validGroups(idxMin).yRange)> 1800 || Rdiff(1) < Rdiff(idxMin))
    gIdx = validGroups(1).idx;
else
    [minCost,validGIdx] = min(validCostFcn); %try using both outlier and sse
    gIdx = validGroups(validGIdx).idx; %get the overall group index
    
    if gIdx > numel(originalGroups) && nALines > 400 %if you're going to use the biggest group, make sure there is no close second
        [secondMinCost, otherIdx] = min(validCostFcn(1:end-1));
        
        if abs(minCost-secondMinCost) > std(validCostFcn);
            gIdx = otherIdx;
        end
    end
end


[groupAnalysis(:).isBest] = deal(false);
groupAnalysis(gIdx).isBest = true;


bestGroup = groupAnalysis(gIdx);

% if gIdx ~=3
%     groupAnalysis(end) = [];
% end

allGroups = groupAnalysis;

function pOut = ransac(dataIn, dw_pxls, fType, maxIter, nPoints, isConvex, userData)
%This function runs n iterations on the x and y data using random x number
%of points returning the polynomial with the best fit
%
%Inputs: data - 2D arrray containing x and y data in columnar fashion [x1,y1;x2,y2],
%dw_pxls = number of pixels for tolerance range; scalar
%nIter = max number of iterations; scalar (optional) [default-60]
%nPoints = number of points randomly selected; scalar (optional) [default-8]
%
%Outputs: p - quadratic polynomial coefficients; returns false if there is
%an error

%Use only unique values in the second column
% xIn = dataIn(:,1); yIn = dataIn(:,2);

%Parse inputs
switch nargin
    case [0 1] %not enough inputs
        pOut = false;
        return
    case 2
        fType = @quadFit;
        maxIter = 60; %default max
        nPoints = 8; %default points
        isConvex = true;
    case 3
        maxIter = 60; %default max
        nPoints = 8; %default points
        isConvex = true;
    case 4
        nPoints = 8; %default points
        isConvex = true;
    case 5
        isConvex = true;
end

switch userData.imgType
    case 1 %in vivo cornea
        Rmin = 6.5; Rmax = 8.3; %mm
    case 2 %in vivo lens
        
        if isConvex
            Rmin = 6; Rmax = 20; %mm
        else
            Rmin = 6; Rmax = 50; %mm
        end
    otherwise
        Rmin = 5; Rmax = 30;
end

%cost threshold
costThresh = 20;

%take only unique y values to minimize the effect of horizontal artifacts
uniqueData = uniqueByCol(dataIn, 2);

%Find out if the percent discarded is too high
if size(uniqueData,1)/size(dataIn,1) < 0.5
    uniqueData = dataIn;
end

xIn = uniqueData(:,1); yIn = uniqueData(:,2);


dataSize = size(uniqueData,1);
xData = uniqueData(:,1); yData = uniqueData(:,2);

%if data size is smaller than number of points, redetermine number of
%points to be 1% of the data
if dataSize < nPoints
    nPoints = round(0.01*dataSize);
end

p={};
costFcn = deal([]);
isGo = true;
count = 0;

xRes = userData.xyRes(1); yRes = userData.xyRes(2);
Rmm = [];
while isGo
    %increment the count
    count = count+1;
    
    %if you haven't found anything, for the posterior lens, use only the
    %central datapoints
    if (count > maxIter/2) && userData.imgType == 2 && ~isConvex
        tol = 1*xRes;
        cntrIdx =  floor(dataSize/2);
        %use points closer to the center
        iIdx = randi([max([1,cntrIdx-tol]),...
            min([cntrIdx+tol,dataSize])],[1,nPoints]);
    else
        
        %randomly select data points unique in the Y direction
        iIdx = randperm(dataSize);
        iIdx = iIdx(1:nPoints);
        
    end
    iX = xData(iIdx); iY = yData(iIdx);
    iXmm = iX/xRes; iYmm = iY/yRes;
    
    %fit the points using least squares
    iP = feval(fType,iX,iY);
    iPmm = feval(fType, iXmm, iYmm);
    iRmm = 1/(2*iPmm(1));
    iFit = polyval(iP,xIn);
    
    
    %create a tolerance range to classify outliers
    upperBound = iFit + dw_pxls/4;
    lowerBound = iFit - dw_pxls/4;
    isOutlier = bitor(yIn > upperBound, yIn < lowerBound);
    nOutliers = sum(logical(nonzeros(isOutlier)));
    
    %determine how well inliers fit to the curve
    inlierX = xIn(~isOutlier); inlierY = yIn(~isOutlier);
    inlierFit = polyval(iP, inlierX);
    sse = mean((inlierFit - inlierY).^2);
    %OTHER OPTIONS FOR SSE:
    %sum((inlierFit - inlierY).^2/max((inlierFit - inlierY).^2));
    %sum((inlierFit - inlierY).^2);
    %std((inlierFit - inlierY).^2);
    %mean((inlierFit - inlierY).^2);
    
    
    
    
    %if quadratic is in the right range, keep it
    isInRadiusRange = (abs(iRmm) < Rmax && abs(iRmm) > Rmin);
    isRightInflection = (iP(1) > 0 && isConvex) || (iP(1) < 0 && ~isConvex);
    if isRightInflection && isInRadiusRange
        
        p{end+1} = iP;
        
        %Cost Fucntion: root(error) + percent of outliers (makes outliers
        %dominating term)
        costFcn(end+1) = (sse)^0.25 + (nOutliers/dataSize)*100;
        Rmm(end+1) = iRmm;
    else %if quadratic is inflected (curved downward), skip it
        
        
        if count > maxIter
            isGo = false;
        end
        
        
        continue
    end
    
    if costFcn(end) < costThresh && count > maxIter %changed || to && (Modified: 1/30/2015)
        isGo = false;
    end
    
    
end %end while

if ~isempty(p)
    p = cat(1,p{:});
else
    pOut = false;
    return
end


[err,idxOut] = min(costFcn);
% disp(count)
% fprintf(1,'\nError is: %10.5f\n\n', err);
pOut = p(idxOut,:);


function [pOut, Rout] = ransacNew(dataIn, dw_pxls, fType, maxIter, nPoints, isConvex, userData)
%This function runs n iterations on the x and y data using random x number
%of points returning the polynomial with the best fit
%
%Inputs: data - 2D arrray containing x and y data in columnar fashion [x1,y1;x2,y2],
%dw_pxls = number of pixels for tolerance range; scalar
%nIter = max number of iterations; scalar (optional) [default-60]
%nPoints = number of points randomly selected; scalar (optional) [default-8]
%
%Outputs: p - quadratic polynomial coefficients; returns false if there is
%an error

%Use only unique values in the second column
% xIn = dataIn(:,1); yIn = dataIn(:,2);

%Parse inputs
switch nargin
    case [0 1] %not enough inputs
        pOut = false;
        return
    case 2
        fType = @quadFit;
        maxIter = 60; %default max
        nPoints = 8; %default points
        isConvex = true;
    case 3
        maxIter = 60; %default max
        nPoints = 8; %default points
        isConvex = true;
    case 4
        nPoints = 8; %default points
        isConvex = true;
    case 5
        isConvex = true;
end


if isfield(userData, 'surfNo')
    switch userData.surfNo
        case 1
            Rmin = 5; Rmax = 9;
        case 2
            Rmin = 5.5; Rmax= 10;
        case 3
            Rmin = 7; Rmax = 20;
        case 4
            Rmin = 7; Rmax = 50;
    end
else
    switch userData.imgType
        case 1 %in vivo cornea
            Rmin = 5; Rmax = 9; %mm
        case 2 %in vivo lens
            
            if isConvex
                Rmin = 3; Rmax = 20; %mm
            else
                Rmin = 3; Rmax = 50; %mm
            end
        otherwise
            Rmin = 5; Rmax = 30;
    end
end

%cost threshold
costThresh = 20;

%take only unique y values to minimize the effect of horizontal artifacts
uniqueData = uniqueByCol(dataIn, 2);


%Find out if the percent discarded is too high
if size(uniqueData,1)/size(dataIn,1) < 0.5
    uniqueData = dataIn;
end

xIn = uniqueData(:,1); yIn = uniqueData(:,2);


dataSize = size(uniqueData,1);
xData = uniqueData(:,1); yData = uniqueData(:,2);

%if data size is smaller than number of points, redetermine number of
%points to be 1% of the data
if dataSize <= nPoints
    nPoints = max([round(0.01*dataSize),4]);
end

if dataSize < 4
    pOut = [];
    Rout = [];
    return
end

maxUniqueSubsets = factorial(dataSize)/(factorial(dataSize - nPoints)*factorial(nPoints));

maxIter = min([maxIter, maxUniqueSubsets]);

p={}; Rmm = [];
costFcn = deal([]);
isGo = true;
count = 0;

xRes = userData.xyRes(1); yRes = userData.xyRes(2);
Rmm = [];
while isGo
    %increment the count
    count = count+1;
    
    %if you haven't found anything, for the posterior lens, use only the
    %central datapoints
    if (count > maxIter/2) && userData.imgType == 2 && ~isConvex
        tol = 1*xRes;
        cntrIdx =  floor(dataSize/2);
        %use points closer to the center
        iIdx = randi([max([1,cntrIdx-tol]),...
            min([cntrIdx+tol,dataSize])],[1,nPoints]);
    else
        
        %randomly select data points unique in the Y direction
        iIdx = randperm(dataSize);
        iIdx = iIdx(1:nPoints);
        
    end
    iX = xData(iIdx); iY = yData(iIdx);
    iXmm = iX/xRes; iYmm = iY/yRes;
    
    %fit the points using least squares
    iP = feval(fType,iX,iY);
    iPmm = feval(fType, iXmm, iYmm); 
    iRmm = 1/(2*iPmm(1));
    iFit = polyval(iP,xIn);
    
    
    %create a tolerance range to classify outliers
    upperBound = iFit + dw_pxls/4;
    lowerBound = iFit - dw_pxls/4;
    isOutlier = bitor(yIn > upperBound, yIn < lowerBound);
    nOutliers = sum(logical(nonzeros(isOutlier)));
    
    %determine how well inliers fit to the curve
    inlierX = xIn(~isOutlier); inlierY = yIn(~isOutlier);
    inlierFit = polyval(iP, inlierX);
    sse = mean((inlierFit - inlierY).^2);
    %OTHER OPTIONS FOR SSE:
    %sum((inlierFit - inlierY).^2/max((inlierFit - inlierY).^2));
    %sum((inlierFit - inlierY).^2);
    %std((inlierFit - inlierY).^2);
    %mean((inlierFit - inlierY).^2);
    
    
    
    
    %if quadratic is in the right range, keep it
    isInRadiusRange = (abs(iRmm) < Rmax && abs(iRmm) > Rmin);
    isRightInflection = (iP(1) > 0 && isConvex) || (iP(1) < 0 && ~isConvex);
    if isRightInflection && isInRadiusRange
        
        p{end+1} = iP;
        Rmm(end+1) = iRmm;
        %Cost Fucntion: root(error) + percent of outliers (makes outliers
        %dominating term)
        costFcn(end+1) = (sse)^0.25 + (nOutliers/dataSize)*100;
        Rmm(end+1) = iRmm;
    else %if quadratic is inflected (curved downward), skip it
        
        
        if count > maxIter
            isGo = false;
        end
        
        
        continue
    end
    
    if count > maxIter %changed || to && (Modified: 1/30/2015) modified 1/4
        isGo = false;
    end
    
    
end %end while

if min(costFcn) > costThresh
    pFlag = true;
end

if ~isempty(p)
    p = cat(1,p{:});
else
    pOut = false;
    Rout = [];
    return
end


[err,idxOut] = min(costFcn);
% disp(count)
% fprintf(1,'\nError is: %10.5f\n\n', err);
pOut = p(idxOut,:);
Rout = Rmm(idxOut);

if false
    
    iP = feval(fType,iX,iY);
    iPmm = feval(fType, iXmm, iYmm);
    iRmm = 1/(2*iPmm(1));
    iFit = polyval(pOut,xIn);
    figure(288);clf; plot(uniqueData(:,1), uniqueData(:,2), '.r'); axis ij;
    hold on; plot(xIn, iFit,'g');
    title(sprintf('This is good: R = %3.2f (Err = %3.2f)\nCount = %d',...
        Rmm(idxOut), err,count));
    pause(.3);
end


function [pOut, Rout, radiusFlag, inflectionFlag, outlierFlag] = robustFit(dataIn,...
    dw_pxls, fType, isConvex, userData)

if isfield(userData, 'surfNo')
    switch userData.surfNo
        case 1
            Rmin = 5; Rmax = 9;
        case 2
            Rmin = 5.5; Rmax= 10;
        case 3
            Rmin = 7; Rmax = 20;
        case 4
            Rmin = 7; Rmax = 50;
    end
else
    switch userData.imgType
        case 1 %in vivo cornea
            Rmin = 5; Rmax = 9; %mm
        case 2 %in vivo lens
            
            if isConvex
                Rmin = 3; Rmax = 20; %mm
            else
                Rmin = 3; Rmax = 50; %mm
            end
        otherwise
            Rmin = 5; Rmax = 30;
    end
end


xRes = userData.xyRes(1); yRes = userData.xyRes(2);

%take only unique y values to minimize the effect of horizontal artifacts
uniqueData = uniqueByCol(dataIn, 2);

%Find out if the percent discarded is too high
if size(uniqueData,1)/size(dataIn,1) < 0.5
    uniqueData = dataIn;
end

xIn = uniqueData(:,1); yIn = uniqueData(:,2);


dataSize = size(uniqueData,1);
xData = uniqueData(:,1); yData = uniqueData(:,2);


x = xData; y = yData;
x_mm = x/xRes; y_mm = y/yRes;

%fit the points using least squares
fitModel = fitlm(x,y,'quadratic','RobustOpts', 'on');
iP = fitModel.Coefficients; iP = flipud(iP{:,1});

fitModel_mm= fitlm(x_mm,y_mm,'quadratic','RobustOpts', 'on');
iPmm = fitModel_mm.Coefficients; iPmm = flipud(iPmm{:,1});
iRmm = 1/(2*iPmm(1));
iFit = polyval(iP,xIn);


%create a tolerance range to classify outliers
upperBound = iFit + dw_pxls/4;
lowerBound = iFit - dw_pxls/4;
isOutlier = bitor(yIn > upperBound, yIn < lowerBound);
nOutliers = sum(logical(nonzeros(isOutlier)));

%determine how well inliers fit to the curve
inlierX = xIn(~isOutlier); inlierY = yIn(~isOutlier);
inlierFit = polyval(iP, inlierX);
sse = mean((inlierFit - inlierY).^2);

%if quadratic is in the right range, keep it
isInRadiusRange = (abs(iRmm) < Rmax && abs(iRmm) > Rmin);
isRightInflection = (iP(1) > 0 && isConvex) || (iP(1) < 0 && ~isConvex);

radiusFlag = ~isInRadiusRange;
inflectionFlag = ~isRightInflection;
outlierFlag = nOutliers/numel(yData) > 0.8;
pOut = iP;
Rout = iRmm;

function output = uniqueByCol(data, col)
%This function returns that data with unique values in a particular column,
%keeps row with first occurence of the number
% input: data - columnar data, col - column number

if col > size(data,2) || col < 1
    output = false;
    return
end

%Get the column of interest
colData = data(:,col);
[~,uniqueIdx,~] = unique(colData,'stable');
output = data(uniqueIdx,:);