
function finalGroups = classifyData(xData,yData, userData, maxDy, isAnalyze)

isPlotAll = false;

if nargin < 5
    isAnalyze = true;
end

if nargin < 4
    
    u_yRes = userData.xyRes(2);
    if userData.scanWidth > 12 && userData.imgType == 1
        maxDy_mm = 0.3;
    else
        %any pixel with greater than 0.3mm difference is considered an outlier
        maxDy_mm = 0.15; %mm
    end
    
    maxDy = maxDy_mm*u_yRes;
end

%make sure data is column
xData = force1D(xData);
yData = force1D(yData);

finalGroups = groupByHistogram(xData, yData, maxDy, isAnalyze, userData);

if isempty(finalGroups)
    return
end

%if there is only one group of data, try the second method
% groups2 = groupByDerivative(xData, yData, maxDy);

%isPlot?
if isPlotAll %|| true
    plotGroups(finalGroups)
end

function plotGroups(groups)

figure(212);clf; hAx = gca;
hold on
for i = 1:numel(groups)
    iData = groups(i).data; iColor = rand(1,3);
    xMin = min(iData(:,1)); xMax= max(iData(:,1));
    iAvg = groups(i).avg(2);
    iPlot = plot(iData(:,1), iData(:,2), 'x','Color', iColor);
    plot([xMin, xMax],[iAvg, iAvg], 'Color', iColor);
end

graphLabels = cellfun(@(n) {sprintf('Group %d',n); sprintf('Group %d',n)}, num2cell(1:numel(groups)),...
    'UniformOutput',false);
legend(cat(1,graphLabels{:}), 'Location', 'EastOutside');
axis ij;

hold off

function finalGroups = groupByHistogram(xData, yData, maxDy, isAnalyze, userData)
%userData has scanWidth and imgType

isPlotAll = false;

isPlot = isPlotAll && isAnalyze;
bins = (min(yData):max(yData))';count =histc(yData,bins);
groupPercentage = 0.1; %min number of pixels to be considered a group
minGroupMembers = round(groupPercentage*numel(yData));

if isAnalyze
    artifactThresh = 0.05*max(xData); %min number of pixels to make an artifact
else
    artifactThresh = 0.5*max(xData); %min number of pixels to make an artifact
end

%if there aren't any groups increase the threshold a little
isGTMaxDy = false; divisor = 4;

while ~any(isGTMaxDy)
    histData = [bins,count]; avgCount = mean(count(count>0));
    dyThresh = max([avgCount/divisor,1]);
    hasCount = count>dyThresh; validValues = histData(hasCount,:);
    
    %Intermediate variables used to find the values that are greater than maxDy
    idx = (2:size(validValues,1))';
    isGTMaxDy = diff(validValues(:,1)) > maxDy;
    isGTIdx= idx(isGTMaxDy);
    
    if divisor ==1
        break
    else
        divisor = divisor -1;
    end
    
end

if isempty(isGTIdx) %even data distribution
    finalGroups = [];
    return
end

dataIdx = sort([isGTIdx-1; isGTIdx]);

%convert to y location index
dataLocation = validValues(dataIdx,1);

while any(diff(dataLocation)==0)
    idx = find(diff(dataLocation)==0)+1;
    dataLocation(idx) = dataLocation(idx)+1;
end

if dataLocation(1) == bins(1)
    isBlankStart = true;
else
    isBlankStart = false;
    dataLocation = [bins(1);dataLocation];
end

allRanges = [dataLocation, [dataLocation(2:end)-1; bins(end)]];

isEven = mod(1:size(allRanges,1),2)==0;

if isBlankStart
    isDataRange = isEven;
else
    isDataRange = ~isEven;
end

zeroRange = allRanges(~isDataRange,:);
dataRange = allRanges(isDataRange,:);

nGroups = size(dataRange,1);

if isempty(nGroups)
    finalGroups = [];
    return
end

%Find weighted means for every data range
%Group Data and find means
origData = [xData,yData]; origDataSorted = sortrows(origData,2);
ySorted = origDataSorted(:,2);
allGroups = struct([]);

for i = 1:nGroups
    iStart = dataRange(i,1); iEnd = dataRange(i,2);
    
    %convert to histogram idx
    [~,idx1] = min(abs(iStart-bins)); [~, idx2]= min(abs(iEnd-bins));
    
    iBins = bins(idx1:idx2); iCount = count(idx1:idx2);
    isArtifact = iCount > artifactThresh;
    %Recreate the original data set from the histogram
    %     zeroCount = iCount == 0;
    %     iCount2 = iCount; iCount2(zeroCount) = [];
    %     iBins2 = iBins; iBins2(zeroCount) = [];
    %     iA = cell2mat(cellfun(@(x) ones(x(1),1)+x(2),...
    %         num2cell(cat(2,iCount2, (0:numel(iCount2)-1)'),2), 'UniformOutput', false));
    %     recreatedData = iBins2(iA);
    %
    %     %Calculate weighted average
    %     iAvg = mean(recreatedData); %round(sum((iBins.*iCount))/sum(iCount));
    %     avgVals(i) = iAvg;
    %     stdVals(i) = std(recreatedData);
    
    
    if isPlot
        if i ==1
            figure;
            hold on;
        end
        plot(iBins,iCount,'Color', rand(1,3));
    end
    
    
    %Original Data %added 3/25/2013
    %find all data within the range
    idxStart = find(abs(ySorted-iStart)==min(abs(ySorted-iStart)),1,'first');
    idxEnd = find(abs(ySorted-iEnd)==min(abs(ySorted-iEnd)),1,'last');
    
    iData = origDataSorted(idxStart:idxEnd,:);
    iData = sortrows(iData,1); %sort by x values
    
    
    %PUT LINEAR ARTIFACT IN IT'S OWN SUBGROUP if artifact is on either end
    %of the group
    
    if isAnalyze && any(isArtifact)
        iY = iData(:,2);
        %remove artifact
        artifactIdx = iBins(isArtifact);
        artifactRange = unique(cell2mat(arrayfun(@(x) (x-2:x+2)', artifactIdx,'UniformOutput',false)));
        
        nArtifacts = numel(artifactRange);
        
        %make it a row vector
        artifactRange = reshape(artifactRange, [1,nArtifacts]);
        iYMatrix = repmat(iY, [1,nArtifacts]);
        artifactMatrix = repmat(artifactRange,[numel(iY),1]);
        
        %check for any artifacts in the bin
        isNotData = logical(sum(iYMatrix == artifactMatrix,2));
        
        iSubData = iData(isNotData,:);
        iSubX = iSubData(:,1); iSubY = iSubData(:,2);
        
        
        %if artifact is not in the middle of the group, separate it
        minPxDiff = 5;
        if all(abs(iSubX - mean(iData(:,1))) > minPxDiff)
            
            %         if any( bitor(abs(iSubX - iData(1,1)) < minPxDiff,...
            %                 abs(iSubX - iData(end,1)) < minPxDiff)) && all(abs(iSubX - mean(iData(:,1))) > minPxDiff)
            
            gIdx = numel(allGroups)+1;
            %Put all the artifacts in it's own group
            allGroups(gIdx).n = size(iSubData,1);
            allGroups(gIdx).avg = [mean(iSubX), mean(iSubY)];
            allGroups(gIdx).data = iSubData;
            allGroups(gIdx).stdev = std(iSubY);
            allGroups(gIdx).hasArtifact = any(isArtifact);
            allGroups(gIdx).artifactIdx = iBins(isArtifact);
            
            %connectivity of the layer %added 4/2/2013
            allGroups(gIdx).connectivity = std(abs(diff(iSubY)));
            
            %remove these horizontal artifacts from the group
            iData(isNotData,:) = [];
            
        else
            %continue
            
        end
        
        isArtifact = false;
        
    end
    
    
    iX = iData(:,1); iY = iData(:,2);
    gIdx = numel(allGroups)+1;
    allGroups(gIdx).n = size(iData,1);
    allGroups(gIdx).avg = [mean(iX),mean(iY)];
    allGroups(gIdx).data = iData;
    allGroups(gIdx).stdev = std(iY);
    allGroups(gIdx).hasArtifact = any(isArtifact);
    allGroups(gIdx).artifactIdx = iBins(isArtifact);
    %get the largest artifact idx %added 4/26/2013
    [~,id] = max(iCount(isArtifact));
    
    if ~isempty(id)
        id = id(1);
    end
    
    artifactBins = iBins(isArtifact);
    allGroups(gIdx).largestArtifactIdx = artifactBins(id);
    
    %connectivity of the layer %added 4/2/2013
    allGroups(gIdx).connectivity = std(abs(diff(iY)));
    
    
    %For plot legend
    stdVals(i) = std(diff(iY)); %std(iY);
    legendStr{i} = sprintf('sigma = %5.2f',stdVals(i));
    
    
    
end

if isPlot
    legend(legendStr);
    title('Histogram Distribution')
end


%If the distance between groups is less than 1 then run it again after
%removing points that are noise


%threshold is more than 5
% if dyThresh < 5
%     whoCalledMe = dbstack;
%
%     allGroups.connectivity
%     %and we are not in a nested function
%     if numel(whoCalledMe) >  2 && ~strcmpi(whoCalledMe(2).name)
%         groupByHistogram(xData, yData, maxDy, false, userData)
%
%     end
% end

finalGroups = allGroups;

if isAnalyze %only analyze if necessary
    
    if isfield(userData, 'surfNo')
        switch userData.surfNo
            case 1
                minStd = 12;
                maxConn = 800; %should be 80
            case 2
                minStd = 5;
                maxConn = 800; %should be 80
            case 3
                minStd = 2;
                maxConn = 80; %should be 80
            case 4
                minStd = 2;
                maxConn = 80; %should be 80
        end
    else
        switch userData.imgType
            case 1 %cornea (in-vivo)
                if userData.scanWidth > 12 && userData.imgType ==1
                    minStd = 5;
                    maxConn = 800; %should be 80
                elseif userData.scanWidth < 6
                    minStd = 5;
                    maxConn = 800; %should be 80
                elseif userData.scanWidth < 11
                    minStd = 5;%12; MODIFIED: 10/20/2015
                    maxConn = 800; %should be 80
                else
                    %remove anything with a stdev less than 6 or greater than 100
                    minStd = 15; %25; modified: 5/12/2014
                    maxConn = 800; %should be 200
                end
            case 2%lens (in-vivo)
                
                if userData.scanWidth > 12
                    minStd = 0;
                    maxConn = 40;
                elseif userData.scanWidth < 6
                    minStd = 3;
                    maxConn = 80; %should be 80
                elseif userData.scanWidth < 10
                    minStd = 2; %3? in semiAutoSegment (7-29-2015)
                    maxConn = 80; %should be 80
                else
                    minStd = 4;
                    maxConn = 80; %should be 80
                end
                
            otherwise
                %settings of lens for now (modified: 3/6/2014)
                if userData.scanWidth < 6
                    minStd = 3;
                    maxConn = 80; %should be 80
                else
                    minStd = 10;
                    maxConn = 80; %should be 80
                end
                
        end
    end
    
    isArtifact = [finalGroups.stdev]  < minStd;
    
    finalGroups(isArtifact) = [];
    
    %also remove any groups that are KNOWN artifacts
    hasArtifact = [finalGroups.hasArtifact];
    finalGroups(hasArtifact) = [];
    
    gCount = [finalGroups(:).n];
    isEnoughMembers = gCount > minGroupMembers;
    if any(isEnoughMembers)
        finalGroups(~isEnoughMembers) = [];
    end
    
    
    
    isRandom = [finalGroups.connectivity]>maxConn;
    if any(~isRandom)
        finalGroups(isRandom) = [];
    end
    
end
