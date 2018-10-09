%Read file list
[fnameIn, pname] = uigetfile('*.tiff','MultiSelect', 'on');

if ~pname
    return
end

if ~iscell(fnameIn)
    fnameIn = {fnameIn};
end

nFiles = numel(fnameIn);
%%

listStr = {'Cornea', 'Lens', 'Retina', 'Cornea+Lens','Cornea+Lens+Retina'};
[uInput,isSlcnValid] = listdlg('ListString',listStr, 'Name','Select Image',...
    'PromptString','Select image type:',...
    'SelectionMode','Single', 'ListSize',[160 130]);

if ~isSlcnValid
    return
end

fnameCell = cell(3,nFiles);

switch uInput 
    case {1,2,3}
        imgTypeLoaded = uInput;
        fnameCell(uInput,:) = fnameIn; %only that type loaded
    case 4 %will give an error if too many or too little files are loaded
        imgTypeLoaded = [1 2];
        fnameCell = reshape(fnameIn, [2,nFiles/2]);
    case 5 %cautious of error
        imgTypeLoaded = [1 2 3];
        fnameCell = reshape(fnameIn, [3, nFiles/3]);
end



%%
%Which image type would you like to analyze
imgTypeStr = {'Cornea', 'Lens','Retina'}; 
imgTypeStr = imgTypeStr(imgTypeLoaded);
imgTypeStr = cat(2,imgTypeStr, 'All');

if numel(imgTypeStr) > 2
    [uInput,isSlcnValid] = listdlg('ListString', imgTypeStr, 'Name',...
        'Select Layer', 'PromptString', 'What layers would you like to segment?',...
        'SelectionMode','Single', 'ListSize',[160 80]);
    
    if ~isSlcnValid
        return
    end
    
    imgTypeIn = uInput;
    
    
else
    imgTypeIn = imgTypeLoaded;
    fnameIn = fnameCell(:,imgTypeIn);
end

if uInput == numel(imgTypeStr)
    %analyze all
    fnameIn = fnameCell;
    processAllFlag = true;
else
    fnameIn = fnameCell(uInput,:);
    processAllFlag = false;
end


%%
%Generate prototypes
clear prototype;
prototype = generateImgPrototype();
prototypeImgs = cat(3,prototype.img);
prototypeImgsFD = cat(3,prototype.imgFD);
allHist = {prototype.histogram};
%%

%Get the imgInfo structure (imgSize, xyRes, scanWidth, scanDepth)
imgInfo = buildImgInfo(fullfile(pname, fnameIn{1}));


settings = struct([]); 
settings(1).xyRes = imgInfo.xyRes; 
settings.imgSize = imgInfo.imgSize;
settings.scanWidth = imgInfo.scanWidth;
settings.scanDepth = imgInfo.scanDepth;
settings.imgType = imgTypeIn; %NEED TO SET THE IMAGE TYPE
%%
%default params
settings.sigma = 6;
settings.fitDivFactor = 1.5;
settings.intFilterSize = [10 5]; %row col? yes
settings.derFilterSize = [5 10];
settings.distBetweenPeaks = 1+imgTypeIn; %mm 2,3
settings.nRecursions = 5;
settings.fitTol = 0.1; %start fit tolerance width in mm
settings.lensROITol = 100; 
settings.artifactTol = 10;
settings.fitType = @quadfit;  
settings.prototype = prototype; %fields: img, histogram, imgFD

xyRes = settings.xyRes;
moduleNo = 7;%:8; %[1 2 3 4 5 6 7 8];
plotModule = 7; %max(mod); %which modules to plot

%%
maxNoOfImgs = 200;

if processAllFlag
    %take max number of images from a random location, but always start
    %with the cornea
    imgProcessEnd = min([nImgs, maxNoOfImgs+mod(maxNoOfImgs,numel(imgTypeLoaded))]);
    imgProcessStart = randi(nImgs-imgProcessEnd, 1,1);
    imgProcessStart = imgProcessStart-mod(imgProcessStart,numel(imgTypeLoaded))+1;
    imgs2process = (1:imgProcessEnd)+ imgProcessStart-1;
else
    %select random subset of images (max 200)
    nImgs = numel(fnameIn);
    randImgNos = randi(nImgs, 1, min([nImgs, maxNoOfImgs]));
    imgNo = randImgNos(1);
    imgs2process = randImgNos;
end


%%
fitOut = cell(6,size(fnameCell,2)); %each layer gets its own row
Rout = nan(6,size(fnameCell,2));

for imgNo = imgs2process
    %%
    clear iImg iImgdb iImgF
    idxNo = find(randImgNos == imgNo);
    iImgdb = imread(fullfile(pname, fnameIn{imgNo})); %how far i am along
    [settings.imgType, imgSetNo] = ind2sub(size(fnameCell), imgNo);
    
    iOutput = segmentImage(iImgdb, settings); %struct: points, curve, R, flag, confidence
%     keyboard;
    
    iPoints = force1D({iOutput.points});
    iR = force1D([iOutput.R]);
    
    switch settings.imgType
        case 1 %check for the third layer and add it 
            fitOut(1:numel(iPoints),imgSetNo) = iPoints;
            Rout(1:numel(iR),imgSetNo) = iR;
        case 2 %replace anterior lens if necessary
            
            nPoints_AL = size(iOutput(2).points,1);
            R_AL = iOutput(2).R;
            
            if nPoints_AL > size(fitOut{3,imgSetNo},1)
                %rewrite AL
                fitOut{3,imgSetNo} = iOutput(2).points;
                Rout(3,imgSetNo) = iOutput(2).R;
            end
            
            %add posterior lens
            fitOut{4,imgSetNo} = iOutput(1).points;
            Rout(4,imgSetNo) = iOutput(1).R;
            
        case 3
    end
end


%%
hFig1 = figure; hFig2 = figure; hFig3 = figure; hFig4 = figure;
for imgNo = randImgNos
%%
clear iImg iImgdb iImgF
idxNo = find(randImgNos == imgNo);
iImgdb = imread(fullfile(pname, fnameIn{imgNo}));
iHist = imhist(iImgdb);

iImgFD = (abs(fftshift(fft2(iImgdb)))); %log
fprintf(1,'\nIdx No: %d; Img No: %d\n',idxNo, imgNo);



%works only for cornea and retina. cornea is sometimes classified as lens
%(false positive)
[~,iImgTypeCalc] = min(cellfun(@(refHist) pdist2(iHist',refHist'), allHist));

if iImgTypeCalc == 2
    %need further analysis to make sure it's not a cornea
    %THIS IS NOT FOOL PROOF
    [~,iImgTypeCalc] = max([corr2(cImgFD, iImgFD), corr2(lImgFD, iImgFD)]);
%     figure(hFig1); imagesc(iImgdb); title(sprintf('ImgType = %d', iImgType2));
%     pause(1);
%     continue
end

   
if iImgTypeCalc ~= imgTypeIn
    figure(hFig1); imagesc(iImgdb); title(sprintf('ImgType = %d', iImgTypeCalc));
    pause(1);
    keyboard;
    %     continue
end

%CHECKS IMAGE TYPE
% figure(hFig1); imagesc(iImgdb); colormap gray;
% title(sprintf('Image Type: %d', iImgTypeCalc)); pause(0.5); continue

if iImgTypeCalc == 3
    continue
end

settings.imgType = iImgTypeCalc;
imgType = iImgTypeCalc;

%%
%==========================================================================
%MODULE 1: Filter the Image
%==========================================================================

if (moduleNo>=1)
    iImgF = double(medfilt2(iImgdb,[10 10]));    
end



% figure(hFig1); imagesc(iImgF), colormap gray; 
% hold on; plot(([1,imgInfo.imgSize(2)]), [50 50]); 
% title(sprintf('Img No: %d (ImgType: %d)', imgNo, iImgTypeCalc));
% pause(1); continue



%%
%==========================================================================
%MODULE 2: Remove artifacts (only in corneal images)
%==========================================================================

if (moduleNo>=2)
    
    %Remove the mean of each row 
    iImg = iImgF - repmat(mean(iImgF,2),[1, size(iImgF,2)]);
    iImg_colSum = sum(iImg, 2);
    iDiff_colSum = sum(iImgF - iImg,2);
    
    %if the mean is too high, it is considered as an artifact
    artifactThresh = mean(iDiff_colSum)+3*std(iDiff_colSum);
    artifactIdx = [];
    
    %Get the mask for the artifacts
    if any(iDiff_colSum>artifactThresh) %you have artifacts
        [~,artifactPeakIdx, peakWidth] = findpeaks(iDiff_colSum, 'MinPeakHeight',...
            artifactThresh-1,'Width', 'halfprom', 'Annotate','extents');
        
        %if the peak width is too wide...it's data
        isData = peakWidth > 15;
        artifactPeakIdx(isData) = [];  peakWidth(isData) = [];
        if any(artifactPeakIdx)
            artifactBounds = round([artifactPeakIdx-peakWidth/2,artifactPeakIdx + peakWidth/2]);
            artifactBounds(artifactBounds<1)= 1;
            artifactBounds(artifactBounds > size(iImgF,1)) = size(iImgF,2);
            artifactRange = arrayfun(@(idx) round(artifactBounds(idx,1)):round(artifactBounds(idx,2)),...
                1:size(artifactBounds,1),'UniformOutput', false);
            artifactIdx = [artifactRange{:}];
        end
    end
    
    if false
        figure(hFig1); clf; imagesc(iImg), colormap gray;
        figure(hFig2); clf; plot(iDiff_colSum,1:numel(iDiff_colSum),'r');
        axis ij; axis tight; pause(0.1); continue
    end
    
    artifactMask = zeros(size(iImgF));
    artifactMask(artifactIdx, :) = 1;
    artifactMask(1:50,:) = 1; %dump the noise at the top of the image
    
    
     if any(plotModule == 2)
        figure(hFig1); imagesc(iImgF); colormap gray; 
        figure(hFig2); imagesc(iImg); colormap gray;
        figure(hFig3); imagesc(artifactMask); colormap gray;       
        pause(0.1);
        continue
     end
    
    
end


%%

%==========================================================================
%MODULE 3: Categorize Image/Find Broad ROIs 
%==========================================================================


if (moduleNo >= 3) && iImgTypeCalc ~=3 %no lens in retinal image
    
    isPlotMod3 = any(plotModule==3);    
    [lensCoords, lensFlag, dataRangeMask] = segmentIris(iImgF, isPlotMod3);
else
    lensCoords = [];
end


if ~isempty(lensCoords)
    lensROITol = settings.lensROITol;
    lensRoi = getLensRoi(iImgF, lensCoords, lensROITol);
else
    lensRoi = true(size(iImgF));
    lensCoords = [1,1;1,size(iImg,2)]; lensFlag = true;
end


if isPlotMod3
    figure(22); imagesc(iImg.*lensRoi); colormap gray;
    if ~isempty(lensCoords)
        hold on;
        plot(lensCoords(:,1), lensCoords(:,2),'.r');
        hold off
        pause(.1);
    end
    continue
end


%%
%==========================================================================
%MODULE 4: Select an ROI
%==========================================================================

%Get AntLensRoi just encase  you have problems with the lens
if ~lensFlag
    switch iImgTypeCalc
        case 1
            [iHLensIdx, jHLensIdx] = find(dataRangeMask(:,:,2));
            lensHIdxStart = min(iHLensIdx);
            
            %everything before the lens starts
            antRoi = false(size(iImg));
            antRoi(1:lensHIdxStart,:) = true;
            posRoi = antRoi;
        case 2
            
            [iLensIdx, jLensIdx]=find(lensRoi);
            lensIdxStart = min(jLensIdx); lensIdxEnd = max(jLensIdx);
            
            [iHLensIdx, jHLensIdx] = find(dataRangeMask(:,:,1));
            
            if range(iHLensIdx) > 100
                lensHIdxStart = max([min(iHLensIdx)-50,1]); lensHIdxEnd = max(iHLensIdx)-100;
                lensHStart = min([lensHIdxStart, lensHIdxEnd]);
                lensHEnd = max([lensHIdxStart,lensHIdxEnd]);
            elseif range(iHLensIdx) > 50
                lensHIdxStart = max([min(iHLensIdx),1]); lensHIdxEnd = max(iHLensIdx);
                lensHStart = min([lensHIdxStart, lensHIdxEnd]);
                lensHEnd = max([lensHIdxStart,lensHIdxEnd]);
            else
                %everything before last roi
                lensHStart = max([min(iHLensIdx)-50,1]);
                [iHLensIdx2, jHLensIdx2] = find(dataRangeMask(:,:,2));
                lensHEnd = min(iHLensIdx2)-100;
                

            end
            
            antRoi= false(size(iImg));
            antRoi(lensHStart:lensHEnd, lensIdxStart:lensIdxEnd) = true;

            posRoi = lensRoi;
        case 3
    end
else
    [antRoi, posRoi] = deal(true(size(iImg)));
end





if iImgTypeCalc == 2 && any(plotModule == 4)
    figure(hFig1); imagesc(antRoi.*iImg); colormap gray
    pause(.1); continue
end



%%
%==========================================================================
%MODULE 5: Calculate Intensity Mask
%==========================================================================

if (moduleNo >=5)
    
    if all(antRoi(:) == posRoi(:))
        iImgFMask = iImgF; iImgFMask(~antRoi) = nan;
        imgMean = mean(iImgFMask(:),1,'omitnan');
        imgStd = std(iImgFMask(:),1,'omitnan');
    else
        imgMean = mean(iImgF(:));
        imgStd = std(iImgF(:));
    end
    intThresh = imgMean + 0.5*imgStd;
    intMask = iImgF > intThresh; %don't use iImg, use iImgF  
    
end


%%
%==========================================================================
%MODULE 6: Calculate Derivative Mask
%==========================================================================

if (moduleNo >=6)
    %ANTERIOR SURFACE
    h = 3; %(x+h) - x
    intMultFactor = 0.5; derMultFactor = 1.5;
    
    [~, diffData_ant] = calculateThreshold(iImg, 0,...
        intMultFactor, derMultFactor, h); %NOT iIMGF
    aScanD_ant = diffData_ant{1}; threshD = diffData_ant{2};
    
    if all(antRoi(:) == posRoi(:)) && ~all(posRoi(:))
         [~, diffData_ant2] = calculateThreshold(iImg.*antRoi, 0,...
        intMultFactor, derMultFactor, h); %NOT iIMGF
        threshD = diffData_ant2{2};
    end
    
    derMask_ant = aScanD_ant > threshD;
    
    
    %POSTERIOR SURFACE
    [~, diffData_pos] = calculateThreshold(iImg, 1,...
        intMultFactor, derMultFactor, h);
    aScanD_pos = diffData_pos{1}; threshD_pos = diffData_pos{2};

    if all(antRoi(:) == posRoi(:)) && ~all(posRoi(:))
         [~, diffData_pos2] = calculateThreshold(iImg.*posRoi, 1,...
        intMultFactor, derMultFactor, h); %NOT iIMGF
        threshD_pos = diffData_pos2{2};
    end
    
    derMask_pos = aScanD_pos > threshD_pos;
    
end


%%
%==========================================================================
%MODULE 7a: Calculate Derivative Mask
%==========================================================================

if any(moduleNo >=7)
    %SET the surface number
    if imgType == 2
        sNo_ant = 4;
    else
        sNo_ant = 1 + 2*(imgType-1);
    end
    
    %discard noise at the top of the image
    top50Mask = false(size(iImg)); 
    if iImgTypeCalc ~=1
        top50Mask(1:50,:) = true;
    else
        
        [iHCornea, jHCornea] = find(dataRangeMask(:,:,1));
        if min(iHCornea) > 50
             top50Mask(1:50,:) = true;
        end
        
    end
    
    %7a: Valid Pixels
    ant_validPixelMask = intMask & derMask_ant & ~top50Mask & antRoi; %& ~artifactMask;
    
    aScanDMasked = aScanD_ant; aScanDMasked(~ant_validPixelMask) = 0;
    
    %use the brightest peak on the first go
    validPoints_ant = findValidPoints(aScanDMasked, false, settings);
    
    rawAntX = validPoints_ant(:,1); rawAntY = validPoints_ant(:,2);
    settings.surfNo = sNo_ant;
    
    %Select points that pertain to a single group
    [bestGroup_ant, ant_allGroups] = selectBestGroup(rawAntX, rawAntY, @quadfit,...
        true, settings);
    
    
    if isempty(bestGroup_ant)
        antData = validPoints_ant;
        iR_ant = -1;
        antFlag = true;
    else
        antData = bestGroup_ant.global.inliers;
        antData_mm = [antData(:,1)./imgInfo.xyRes(1),...
            antData(:,2)./imgInfo.xyRes(:,2)];
        p_mm = polyfit(antData_mm(:,1), antData_mm(:,2),2);
        iR_ant = 1/(2*p_mm(1));
        antFlag = false;
    end 
    
    if iR_ant > 15 || iR_ant < 3 %something's up
        keyboard;
    end
    
    
    switch settings.imgType
        case 1
            fitLimits = [1, size(iImg,2)];
        case 2
            fitLimits = lensCoords(:,1)';
    end
    
    
    
    %Do a magnet fit instead of a recursive fit
    isCorrected = false; isMagnet = false;  fType = 'quad';%'conic'; 
    

    iXFit = fitLimits(1):fitLimits(2);
    
    [antFit, ant_R, antFit_mm, ant_curve, ant_curve_mdl] = fitDataPoints(antData, xyRes,...
        sNo_ant, fType, iImgdb, iXFit, isMagnet, isCorrected);
   
    
    %Include all the inliers from the valid pixel mask with a very small
    %tolerance window
    fitTol = round(0.03 * imgInfo.xyRes(2)); %0.01mm
    antPoints = detectInliers(ant_validPixelMask, antFit, fitTol);
    
    
    
%      if any(plotModule == 7)
        
%         figure(hFig2); imagesc(ant_validPixelMask); colormap gray;
        figure(hFig1); imagesc(iImg); colormap gray;
        hold on; plot(validPoints_ant(:,1), validPoints_ant(:,2),'.r');
        plot(antFit(:,1), antFit(:,2),'b');
        plot(antPoints(:,1), antPoints(:,2),'om'); 
        pause(0.5);
%         continue
%      end
%%    
%==========================================================================
%Module 7b: Posterior Surface
%==========================================================================
    %SET surface Number
    if imgType == 2
        sNo_pos = 3;%3;
    else
        sNo_pos = 2 + 2*(imgType-1);
    end
    
    if imgType == 1
        xFit = antFit(:,1); yFit = antFit(:,2);
        %Discard everything above anterior surface
        idx2Exclude = segmentMatrix([force1D(xFit), force1D(round(yFit+0.4*xyRes(2)))], iImg);
        tempImgMask = ones(size(iImg)); tempImgMask(idx2Exclude) = 0;
        
        %discard everything a below maximum distance away
        bottomIdx2Exclude = segmentMatrix([force1D(xFit), force1D(round(yFit+0.4*xyRes(2)*6))], iImg);
        tempImgMask(~bottomIdx2Exclude) = 0;
        posRoi = tempImgMask;
    end


    validPixelMask_pos = intMask & derMask_pos & posRoi;

    aScanDMasked_pos = aScanD_pos; aScanDMasked_pos(~validPixelMask_pos)=0;
    
    validPoints_pos = findValidPoints(aScanDMasked_pos, false, settings);
    
    if ~any(validPoints_pos)
        posPoints = [];
        return
    end
    
    rawPosX = validPoints_pos(:,1); rawPosY = validPoints_pos(:,2);
    
    %BestGroup
    settings.surfNo = sNo_pos;
    isConvex = settings.imgType == 1;
    [bestGroup_pos, pos_allGroups] = selectBestGroup(rawPosX, rawPosY, @quadfit,...
        isConvex, settings);
    
    if isempty(bestGroup_pos)
        posData = validPoints_pos;
        posFlag = true;
    else
        posData = bestGroup_pos.global.inliers;
        posFlag = false;
    end
    
      
    %Do a magnet fit instead of a recursive fit
    isCorrected = false; isMagnet = false;  fType = 'quad';%'conic'; 
    

    %iXFit = fitLimits(1):fitLimits(2);
    
    [posFit, pos_R, pos_xyFit_mm, pos_curve, pos_curve_mdl] = fitDataPoints(posData, xyRes,...
        sNo_pos, fType, iImgdb, iXFit, isMagnet, isCorrected);
   
    
    %Include all the inliers from the valid pixel mask with a very small
    %tolerance window
    %fitTol = round(0.01 * imgInfo.xyRes(2)); %0.01mm
    posPoints = detectInliers(validPixelMask_pos, posFit, fitTol);
    

    if isempty(posPoints) && iImgTypeCalc == 2
        disp('Listen Linda, I cant find the anterior lens'); 
        posFlag = true;
    end
        
    
    
    if any(plotModule == 7)
        
        figure(hFig2); imagesc(ant_validPixelMask); colormap gray;
        figure(hFig1); imagesc(iImg); colormap gray;
        hold on; plot(validPoints_ant(:,1), validPoints_ant(:,2),'.r');
        plot(antFit(:,1), antFit(:,2),'b');
        plot(antPoints(:,1), antPoints(:,2),'og')
        
        
%           figure(hFig3); imagesc(iImg); colormap gray;
        hold on; plot(validPoints_pos(:,1), validPoints_pos(:,2),'.m');
        plot(posFit(:,1), posFit(:,2),'b');
        if ~isempty(posPoints)
            plot(posPoints(:,1), posPoints(:,2),'og')
        else
            plot(posData(:,1), posData(:,2),'.g');
        end
        
        
%         plot(antData(:,1), antData(:,2), '.g'); hold off;
        title(sprintf('R_a_n_t = %3.2fmm (%d points)    R_p_o_s = %3.2fmm (%d points)',...
            ant_R, size(antPoints,1),  pos_R, size(posPoints,1)));
%         title(sprintf('Image: %d/%d', idxNo, numel(randImgNos)));
        pause(.5);
        
        continue
    end
    
    
end





%%
%==========================================================================
%MODULE 8: Calculate Derivative Mask
%==========================================================================

%%
end

