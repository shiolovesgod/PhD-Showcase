function [output, lensCoords] = segmentImage(imgIn, settings)
%INPUTS
    %imgIn - raw unfiltered db Image
    %Settings used:
        %xyRes = [xRes, yRes]
        %prototype = 3x1 struct with the following fields [fcn generateImgPrototype]
            %img, imgFD, histogram
            %prototype(1,2,3) --> [cornea info, lens info, retina info]
        %surfNo (usually set within this function before use)
        %imgType (usually changed within this function if necessary)
        %scanWidth (default: 8mm)
            
%OUTPUT
    %output - mx1 structure with the following fields
        %points - mx2
        %curve - cfit object
        %R - scalar double
        %flag - scalar booleean
        %confidence - scalar double (unassigned right now)
      



if ~isfield(settings, 'scanWidth')
    settings.scanWidth = 8; %default scan width = 8mm
end


if ~isfield(settings, 'imgType')
    settings.imgType = [];
end

xyRes = settings.xyRes;
prototype = settings.prototype;
imgTypeIn = settings.imgType;
lensROITol = settings.lensROITol;




%Initilalize Outputs
output = struct('points', [],'curve', [],'R', [],'flag', [],'confidence',[]);


%Initial settings
moduleNo = 7;%:8; %[1 2 3 4 5 6 7 8];%
plotModule = [];%7; %max(mod); %which modules to plot


figNo = 268;

if ~isempty(plotModule)
    hFig1 = figure(figNo); clf; hFig2 = figure(figNo+1); clf;
%     hFig3 = figure(figNo+3); hFig4 = figure(figNo+4);
end
%%

try
    
    imgTypeCalc = getImgType(imgIn, prototype);
catch
    %protype not same size
    imgTypeCalc = -1;
end
   
if ~isempty(imgTypeIn) && imgTypeCalc ~= imgTypeIn
%     figure(88);clf; imagesc(imgIn); title(sprintf('ImgType = %d', imgTypeCalc));
%     pause(1);
%     keyboard;
    %     return
    
    if true %imgTypeCalc ~= 3 %MODIFIED 12-2-2015
        imgTypeCalc = imgTypeIn;
    end
%     imgTypeCalc = 3;%imgTypeIn;
end

%CHECKS IMAGE TYPE
% figure(hFig1); imagesc(iImgdb); colormap gray;
% title(sprintf('Image Type: %d', imgTypeCalc)); pause(0.5); return

if imgTypeCalc == 3
    lensCoords = [];
    return
end

imgType = imgTypeCalc;
settings.imgType= imgType;

if ~isfield(settings, 'distBetweenPeaks')
    settings.distBetweenPeaks = imgTypeIn+1;
end

%%
%==========================================================================
%MODULE 1: Filter the Image
%==========================================================================

if (moduleNo>=1)
    imgF = double(medfilt2(imgIn,[2 2]));    
end


% figure(hFig1); imagesc(iImgF), colormap gray; 
% hold on; plot(([1,imgSize(2)]), [50 50]); 
% title(sprintf('Img No: %d (ImgType: %d)', imgNo, imgTypeCalc));
% pause(1); return



%%
%==========================================================================
%MODULE 2: Remove artifacts (only in corneal images)
%==========================================================================

if (moduleNo>=2)
    
    %Remove the mean of each row 
    img = imgF - repmat(mean(imgF,2),[1, size(imgF,2)]);
    iImg_colSum = sum(img, 2);
    iDiff_colSum = sum(imgF - img,2);
    
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
            artifactBounds(artifactBounds > size(imgF,1)) = size(imgF,2);
            artifactRange = arrayfun(@(idx) round(artifactBounds(idx,1)):round(artifactBounds(idx,2)),...
                1:size(artifactBounds,1),'UniformOutput', false);
            artifactIdx = [artifactRange{:}];
        end
    end
    
    if false
        figure(hFig1); clf; imagesc(iImg), colormap gray;
        figure(hFig2); clf; plot(iDiff_colSum,1:numel(iDiff_colSum),'r');
        axis ij; axis tight; pause(0.1); return
    end
    
    artifactMask = zeros(size(imgF));
    artifactMask(artifactIdx, :) = 1;
    artifactMask(1:50,:) = 1; %dump the noise at the top of the image
    
    
     if any(plotModule == 2)
        figure(hFig1); imagesc(imgF); colormap gray; 
        figure(hFig2); imagesc(img); colormap gray;
        figure(hFig3); imagesc(artifactMask); colormap gray;       
        pause(0.1);
        return
     end
    
    
end


%%

%==========================================================================
%MODULE 3: Categorize Image/Find Broad ROIs 
%==========================================================================


if (moduleNo >= 3) && imgTypeCalc ~=3 %no lens in retinal image
    
    isPlotMod3 = any(plotModule==3);    
    [lensCoords, lensFlag, dataRangeMask, lensMaskF] = segmentIris(imgF, imgType,isPlotMod3);
else
    lensCoords = [];
end


if ~isempty(lensCoords)
    lensRoi = getLensRoi(imgF, lensCoords, lensROITol);
else
    lensRoi = true(size(imgF));
    lensCoords = [1,1;1,size(img,2)]; lensFlag = true;
end


if isPlotMod3
    figure(22); imagesc(img.*lensRoi); colormap gray;
    if ~isempty(lensCoords)
        hold on;
        plot(lensCoords(:,1), lensCoords(:,2),'.r');
        hold off
%         pause(.1);
    end
    return
end


%%
%==========================================================================
%MODULE 4: Select an ROI
%==========================================================================

%Get AntLensRoi just encase  you have problems with the lens
if ~lensFlag
    switch imgTypeCalc
        case 1
            [iHLensIdx, jHLensIdx] = find(dataRangeMask(:,:,2));
            lensHIdxStart = min(iHLensIdx);
            
            %everything before the lens starts
            antRoi = false(size(img));
            antRoi(1:lensHIdxStart,:) = true;
            
            %use only central region for poorly defined posterior cornea
            posRoi = antRoi;
            posRoi(1) = 0;
        case 2
            
            [iLensIdx, jLensIdx]=find(lensRoi);
            lensIdxStart = min(jLensIdx); lensIdxEnd = max(jLensIdx);
            
            [iHLensIdx, jHLensIdx] = find(dataRangeMask(:,:,1));
            iHCenter = unique(mean(iHLensIdx));
            %get the center of that range and offset it by x pixels
            
            lensAntTol = 200; %pixels
            if range(iHLensIdx) > 50 %really wide
                %take the brightest peak in the ROI and get a tol range
                [~,brightestPeak] = max(sum(img.*dataRangeMask(:,:,1),2));
                lensHStart =  max([brightestPeak - lensAntTol,1]);
                lensHEnd = brightestPeak + lensAntTol;

            else
                %everything before last roi
                lensHStart = max([min(iHLensIdx)-lensAntTol,1]);
                [iHLensIdx2, jHLensIdx2] = find(dataRangeMask(:,:,2));
                lensHEnd = min(iHLensIdx2)-100;
                
            end
%             if range(iHLensIdx) > 100
%                 lensHStart = min([lensHIdxStart, lensHIdxEnd]);
%                 lensHEnd = max([lensHIdxStart,lensHIdxEnd]);
%             elseif range(iHLensIdx) > 50
%                 lensHIdxStart = max([min(iHLensIdx),1]); lensHIdxEnd = max(iHLensIdx);
%                 lensHStart = min([lensHIdxStart, lensHIdxEnd]);
%                 lensHEnd = max([lensHIdxStart,lensHIdxEnd]);
%             else
%                 %everything before last roi
%                 lensHStart = max([min(iHLensIdx)-50,1]);
%                 [iHLensIdx2, jHLensIdx2] = find(dataRangeMask(:,:,2));
%                 lensHEnd = min(iHLensIdx2)-100;
%                 
% 
%             end
            
            antRoi= false(size(img));
            antRoi(lensHStart:lensHEnd, lensIdxStart:lensIdxEnd) = true;

            posRoi = lensRoi;
        case 3
    end
else
    [antRoi, posRoi] = deal(true(size(img)));
end





if imgTypeCalc == 2 && any(plotModule == 4)
    figure(hFig1); imagesc(antRoi.*img); colormap gray
    pause(.1); return
end



%%
%==========================================================================
%MODULE 5: Calculate Intensity Mask
%==========================================================================

if (moduleNo >=5)
    
    
    if numel(find(antRoi)) > numel(find(posRoi)) 
        biggerMask = antRoi;
        smallerMask = posRoi;
    else
        biggerMask = posRoi;
        smallerMask = antRoi;
    end 
    
    %if pos mask is a subset of anterior roi.
    if all(all((bitor(antRoi, posRoi) == biggerMask)));
        iImgFMask = imgF; iImgFMask(~smallerMask) = nan;
        imgMean = mean(iImgFMask(:),1,'omitnan');
        imgStd = std(iImgFMask(:),1,'omitnan');
    else
        imgMean = mean(imgF(:));
        imgStd = std(imgF(:));
    end
    intThresh = imgMean + 0.5*imgStd;
    intMask = imgF > intThresh; %don't use iImg, use iImgF  
    
end


%%
%==========================================================================
%MODULE 6: Calculate Derivative Mask
%==========================================================================

if (moduleNo >=6)
    %ANTERIOR SURFACE
    h = 3; %(x+h) - x
    intMultFactor = 0.5; derMultFactor = 1.5;
    
    [~, diffData_ant] = calculateThreshold(img, 0,...
        intMultFactor, derMultFactor, h); %NOT iIMGF
    aScanD_ant = diffData_ant{1}; threshD = diffData_ant{2};
    
    if all(antRoi(:) == posRoi(:)) && ~all(posRoi(:))
         [~, diffData_ant2] = calculateThreshold(img.*antRoi, 0,...
        intMultFactor, derMultFactor, h); %NOT iIMGF
        threshD = diffData_ant2{2};
    end
    
    derMask_ant = aScanD_ant > threshD;
    
    
    %POSTERIOR SURFACE
    [~, diffData_pos] = calculateThreshold(img, 1,...
        intMultFactor, derMultFactor, h);
    aScanD_pos = diffData_pos{1}; threshD_pos = diffData_pos{2};

    if all(antRoi(:) == posRoi(:)) && ~all(posRoi(:))
         [~, diffData_pos2] = calculateThreshold(img.*posRoi, 1,...
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
    top50Mask = false(size(img)); 
    if imgTypeCalc ~=1
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
        antData_mm = [antData(:,1)./xyRes(1),...
            antData(:,2)./xyRes(:,2)];
        p_mm = polyfit(antData_mm(:,1), antData_mm(:,2),2);
        iR_ant = 1/(2*p_mm(1));
        antFlag = false;
    end 
    
    if iR_ant > 15 || iR_ant < 3 %something's up
        antFlag = true;
    end
    
    
    switch imgType
        case 1
            fitLimits = [1, size(img,2)];
        case 2
            if lensFlag
                fitLimits = [1, size(img,2)];
            else
                fitLimits = lensCoords(:,1)';
            end
    end
    
    
    
    %Do a magnet fit instead of a recursive fit
    isCorrected = false; isMagnet = false;  fType = 'quad';%'conic'; 
    

    iXFit = fitLimits(1):fitLimits(2);
    
    [antFit, ant_R, antFit_mm, ant_curve, ant_curve_mdl] = fitDataPoints(antData, xyRes,...
        sNo_ant, fType, imgIn, iXFit, isMagnet, isCorrected);
   
    
    %Include all the inliers from the valid pixel mask with a very small
    %tolerance window
    fitTol = round(0.03 * xyRes(2)); %0.01mm
    antPoints = detectInliers(ant_validPixelMask, antFit, fitTol);
    
    output(1).points = antPoints; 
    output(1).R = ant_R;
   
     if any(plotModule == 7)
        
%         figure(hFig2); imagesc(ant_validPixelMask); colormap gray;
        figure(hFig1); imagesc(img); colormap gray;
        hold on; plot(validPoints_ant(:,1), validPoints_ant(:,2),'.r');
        plot(antFit(:,1), antFit(:,2),'b');
        plot(antPoints(:,1), antPoints(:,2),'om'); 
%         pause(0.5);
%         return
     end
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
        idx2Exclude = segmentMatrix([force1D(xFit), force1D(round(yFit+0.4*xyRes(2)))], img);
        tempImgMask = ones(size(img)); tempImgMask(idx2Exclude) = 0;
        
        %discard everything a below maximum distance away
        bottomIdx2Exclude = segmentMatrix([force1D(xFit), force1D(round(yFit+0.4*xyRes(2)*6))], img);
        tempImgMask(~bottomIdx2Exclude) = 0;
        posRoi = tempImgMask;
        
        %ALSO....rescale the intensity mask for the posterior cornea
        iImgFMask_pos = imgF; iImgFMask_pos(~posRoi) = nan;
        imgMean_pos = mean(iImgFMask_pos(:),1,'omitnan');
        imgStd_pos = std(iImgFMask_pos(:),1,'omitnan');
        
        intThresh_pos = imgMean_pos + 0.5*imgStd_pos;
        intMask_pos = imgF > intThresh_pos; %don't use iImg, use iImgF
    else
        intMask_pos = intMask;
    
    end
    
    
    validPixelMask_pos = intMask_pos & derMask_pos & posRoi;

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
        sNo_pos, fType, imgIn, iXFit, isMagnet, isCorrected);
   
    
    %Include all the inliers from the valid pixel mask with a very small
    %tolerance window
    %fitTol = round(0.01 * xyRes(2)); %0.01mm
    posPoints = detectInliers(validPixelMask_pos, posFit, fitTol);
    

    if isempty(posPoints) && imgTypeCalc == 2
        disp('Listen Linda, I cant find the anterior lens'); 
        posFlag = true;
    end
    
    output(2).points = posPoints;
    output(2).R = pos_R;
    ouptut(2).flag = posFlag;
    
%==========================================================================
%Module 7c: Anterior Lens on Cornea
%==========================================================================
    if imgTypeCalc == 1 && ~all(lensRoi(:))
        %Look for the lens
        sNo_pos2 = 3;
        
        
        validPixelMask_pos2 = intMask & derMask_ant & lensRoi;
        
        aScanDMasked_pos2 = aScanD_ant; aScanDMasked_pos2(~validPixelMask_pos2)=0;
        
        validPoints_pos2 = findValidPoints(aScanDMasked_pos2, false, settings);
        
        if ~any(validPoints_pos2)
            posPoints = [];
            return
        end
        
        rawPosX2 = validPoints_pos2(:,1); rawPosY2 = validPoints_pos2(:,2);
        
        %BestGroup
        settings.surfNo = sNo_pos2;
        isConvex2 = settings.imgType == 1;
        [bestGroup_pos2, pos_allGroups2] = selectBestGroup(rawPosX2, rawPosY2, @quadfit,...
            isConvex2, settings);
        
        if isempty(bestGroup_pos2)
            posData2 = validPoints_pos;
            posFlag2 = true;
        else
            posData2 = bestGroup_pos2.global.inliers;
            posFlag2 = false;
        end
        
        
        %Do a magnet fit instead of a recursive fit
        isCorrected = false; isMagnet = false;  fType = 'quad';%'conic';
        
        
        %iXFit = fitLimits(1):fitLimits(2);
        
        [posFit2, pos_R2, pos_xyFit_mm2, pos_curve2, pos_curve_mdl2] = fitDataPoints(posData2, xyRes,...
            sNo_pos2, fType, imgIn, iXFit, isMagnet, isCorrected);
        
        
        %Include all the inliers from the valid pixel mask with a very small
        %tolerance window
        %fitTol = round(0.01 * xyRes(2)); %0.01mm
        posPoints2 = detectInliers(validPixelMask_pos2, posFit2, fitTol);
        output(3).points = posPoints2;
        output(3).R = pos_R2;
        
    else
        posPoints2 = [];
        pos_R2 = [];
        lensFlag2= false;
    end
        
   
    
    
    if any(plotModule == 7)
        
        figure(hFig2); imagesc(ant_validPixelMask); colormap gray;
        figure(hFig1); imagesc(img); colormap gray;
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
        
        
        if ~isempty(posPoints2)
            plot(validPoints_pos2(:,1), validPoints_pos2(:,2),'.m')
            plot(posFit2(:,1), posFit2(:,2),'b');
            plot(posPoints2(:,1), posPoints2(:,2),'og');   
            plot(posData2(:,1), posData2(:,2),'xc');
            
        end
        
%         plot(antData(:,1), antData(:,2), '.g'); hold off;
        title(sprintf('R_a_n_t = %3.2fmm (%d points)    R_p_o_s = %3.2fmm (%d points)',...
            ant_R, size(antPoints,1),  pos_R, size(posPoints,1)));
%         title(sprintf('Image: %d/%d', idxNo, numel(randImgNos)));
%         pause(.5);
        
        return
    end
    
    
    
end

