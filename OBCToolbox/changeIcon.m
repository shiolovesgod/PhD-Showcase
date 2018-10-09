function changeIcon(hFig, img, alpha)

if ~strcmp(get(hFig, 'Type'), 'figure')
    return
end

%temporarily create image filename
tempImgFileName = 'x123nertempppImg.png';

if nargin > 2
    imwrite(img, tempImgFileName, 'png','Alpha', alpha);
else
    imwrite(img, tempImgFileName, 'png');
end

jframe  = get(hFig, 'javaframe');
jIcon = javax.swing.ImageIcon(tempImgFileName);
jframe.setFigureIcon(jIcon)

delete(tempImgFileName)