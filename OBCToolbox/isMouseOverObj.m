function output = isMouseOverObj(hFig,hObj)

set([hFig hObj],'Units','Pixels')
imgPos = get(hObj,'Position');

figPos=get(hFig,'Position');
absImgPos = figPos(1:2)+imgPos(1:2);
absImgPos = [absImgPos, absImgPos+imgPos(3:4)];
xBound = absImgPos([1,3]); yBound = absImgPos([2,4]);

mousePointer = get(0,'PointerLocation');
xPointCoord = mousePointer(1); yPointCoord = mousePointer(2);

xCondition = (xPointCoord >= xBound(1)) && xPointCoord <=xBound(2)+10; %+10 means they can click to the right of the box
yCondition = (yPointCoord >= yBound(1)) && yPointCoord <=yBound(2);

output = all([xCondition, yCondition]);