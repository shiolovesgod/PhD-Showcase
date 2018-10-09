function prototype = generateImgPrototype()
%generates a "protoype" of each image type: cornea, lens, retina
%OUTPUT
    %prototype = 3x1 struct with the following fields
        %img, imgFD, histogram
        %prototype(1,2,3) --> [cornea info, lens info, retina info]


cImg = generateGenericImg(1);
lImg = generateGenericImg(2);
rImg = generateGenericImg(3);

prototype(1).img = cImg;
prototype(2).img = lImg;
prototype(3).img = rImg;
%imgType: 1-cornea, 2-lens, 3-retina

cImgFD = (abs(fftshift(fft2(cImg))));
lImgFD = (abs(fftshift(fft2(lImg))));
rImgFD = (abs(fftshift(fft2(rImg))));

prototype(1).imgFD = cImgFD;
prototype(2).imgFD = lImgFD;
prototype(3).imgFD = rImgFD;


%%
%Get each histogram
for j = 1:numel(prototype)
    jImg = prototype(j).img;
    if isempty(jImg)
        continue
    end
    prototype(j).histogram = imhist(jImg);
end
