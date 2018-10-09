function img = setWhiteBack(img, alpha, isPlot)
%sets image background to white
%img = setWhiteBack(img, alpha, isPlot)

if nargin == 2
    isPlot = false;
end

allZeros = [];

%modify anything that already has all zeros
for k = 1:size(img,3)
    if isempty(allZeros)
        allZeros = sum(img,3)== 0;
    end
    
    iImg = img(:,:,k);
    iImg(allZeros) = 1;
    img(:,:,k) = iImg;
end


alpha(alpha>0) = 1;

switch class(img)
    case 'uint16'
        maxValue = (2^16)-1;
    case 'unit8'
        maxValue = (2^8)-1;
    case 'unit32'
        maxValue = (2^64)-1;
    case 'unit64'
        maxValue = (2^64)-1;
    otherwise
        maxValue = max(img(:));
end


for k = 1:size(img, 3)
    img(:,:,k) = img(:,:,k).*alpha;
end

allZeros = [];

for k = 1:size(img,3)
    if isempty(allZeros)
        allZeros = sum(img,3)== 0;
    end
    
    iImg = img(:,:,k);
    iImg(allZeros) = maxValue;
    img(:,:,k) = iImg;
end

if isPlot
    figure; imagesc(img);
end