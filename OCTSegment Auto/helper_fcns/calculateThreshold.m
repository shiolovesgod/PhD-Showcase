
function [aScanInt, aScanDerivative] = calculateThreshold(img, mode,...
    intMultFactor, derMultFactor, h)
%mode determines positive, negative or both
%aScanInt = {aScanImg, threshold}
%aScanDer = {aScanProcessed, derivativeThresh}


img = double(img);
avgInt = mean(img(:));
stdevI = std(img(:));
intThresh = max([avgInt + stdevI*intMultFactor,0]); %fail safe to make sure the min is 0


%Calculate Intensity
avgI = avgInt;
normalizedImg = img-avgI;
imgDer = getDerivative(normalizedImg, h);


switch mode %added 3/25/2013
    case 0 %positive derivative
        imgDer(imgDer < 0) = 0;
        imgD = imgDer;
    case 1 %negative derivative
        imgDer(imgDer > 0) = 0;
        imgDer = abs(imgDer);
        imgD = imgDer;
        
    case 2 %use both edges to create the threshold
        imgD = imgDer;
        imgD = abs(imgD);
end

%for calculations, remove all zero values
imgD(imgD==0)=[];
avgDiff = mean(imgD(:));
stdevDiff = std(imgD(:));
derivThresh = avgDiff + stdevDiff*derMultFactor;

aScanInt = {img, intThresh};
aScanDerivative = {imgDer, derivThresh};

function matrixOut = getDerivative(matrixIn, h)

A = matrixIn; B = circshift(matrixIn,[-h, 0]);
matrixOut = B-A; matrixOut(end-h+1:end,:) = [];
matrixOut = cat(1, zeros(h,size(matrixIn,2)), matrixOut);
