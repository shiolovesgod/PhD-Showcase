function imgTypeCalc = getImgType(imgIn, prototype)
%This function classifies an image as a cornea, lens or retina. Blinks are
%classified as a retina.

if nargin < 2 %no prototype in
    
    cImg = generateGenericImg(1);
    lImg = generateGenericImg(2);
    rImg = generateGenericImg(3);
    
    prototype(1).img = cImg;
    prototype(2).img = lImg;
    prototype(3).img = rImg;
    %imgType: 1-cornea, 2-lens, 3-retina
    prototypeImgs = cat(3,cImg,lImg,rImg);
    
    cImgFD = (abs(fftshift(fft2(cImg))));
    lImgFD = (abs(fftshift(fft2(lImg))));
    rImgFD = (abs(fftshift(fft2(rImg))));
    
    prototype(1).imgFD = cImgFD;
    prototype(2).imgFD = lImgFD;
    prototype(3).imgFD = rImgFD;
    
    prototypeImgsFD = cat(3,cImgFD,lImgFD,rImgFD);
    
    %%
    %Get each histogram
    for j = 1:numel(prototype)
        jImg = prototype(j).img;
        if isempty(jImg)
            continue
        end
        prototype(j).histogram = imhist(jImg);
    end
end

allFD = {prototype.imgFD};
allHist = {prototype.histogram};


imgHist = imhist(imgIn);

imgFD = (abs(fftshift(fft2(imgIn))));
%works only for cornea and retina. cornea is sometimes classified as lens
%(false positive)

[~,imgTypeCalc] = min(cellfun(@(refHist) pdist2(imgHist',refHist'), allHist));

if imgTypeCalc == 2
    %need further analysis to make sure it's not a cornea
    [~,imgTypeCalc] = max(cellfun(@(refImgFD) corr2(refImgFD, imgFD), allFD(1:2)));
%     [~,imgTypeCalc] = max([corr2(cImgFD, imgFD), corr2(lImgFD, imgFD)]);
%     figure(hFig1); imagesc(iImgdb); title(sprintf('ImgType = %d', iImgType2));
%     pause(1);
end

   