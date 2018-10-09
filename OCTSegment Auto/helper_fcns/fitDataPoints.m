

function [xyFit, R, xyFit_mm, curve, curve_mdl] = fitDataPoints(xyVals, xyRes, sNo, fType, imgIn,xFit, isMagnet, isCorrected)
%INPUTS
    %xyVals - mx2, input in pixels of the surface segmentaiton points
    %sNo - scalar, surface number
    %fType - string, fit type {'quad', 'conic', 'spherical', 'line'}
    %xFit

isPlot = false; %???? let user decide


if isempty(xFit)
    xFit = 1:xSize;
end

if isMagnet %NOT THE LATEST VERSION OF MAGNET FIT
    [xyOut, ptWeight] = getMagnetWeights(xyVals, imgIn, xyRes, xFit);
    cropX = xyOut(:,1);s
    cropY = xyOut(:,2);

else
    xyValsFilt = filterNoOverlap(xyVals, 5); %maybe i need smaller filters for surfaces with less points
    cropX = xyValsFilt(:,1); 
    cropY = xyValsFilt(:,2);
    ptWeight = ones(size(cropX));
end

%%
%..........................................................................
%RUN WEIGHTED FIT
%..........................................................................

cropX_mm = cropX/xyRes(1); cropY_mm = cropY/xyRes(2);
xFit_mm = xFit/xyRes(1);


switch fType
    case 'quad' %b = a, b, c of ax^2+bx+c
        fitFcn = @(b,x) b(1)*x.^2 + b(2)*x + b(3);
        curve_mdl = fitlm(cropX_mm, cropY_mm,'quadratic',...
            'Weight', ptWeight);%,'CoefficientNames',{'a','b','c'});

        fStr = fittype('a*x^2+b*x+c');
        coeffVals = flipud(curve_mdl.Coefficients{:,1}); %a, b, c
        R = 1/(2*coeffVals(1));
        
    case 'conic' %turn this back into a curve fit?
        [curve_mdl,~] = newConicFit(cropX_mm, cropY_mm, sNo, isCorrected, ptWeight);
        
        fStr = fittype(curve_mdl.Formula.Expression);
        coeffVals = curve_mdl.Coefficients{:,1};
        R = 1/coeffVals(1);
         
    case 'spherical'
        [curve_mdl, ~] = circFit([cropX_mm, cropY_mm], sNo, ptWeight);
               
        fStr = fittype(curve_mdl.Formula.Expression);
        coeffVals = curve_mdl.Coefficients{:,1};
        R = coeffVals(1);
        
    case 'line'
        curve_mdl = fitlm(cropX_mm, cropY_mm,'constant',...
            'Weight', ptWeight);%,'CoefficientNames',{'a','b','c'});
        
        fStr = fittype('0*x+b'); %y = mx + b
        coeffVals = curve_mdl.Coefficients{1,1};
        R = inf;
end

%TURN THE ANSWER INTO A cFit
coeffStr = sprintf('coeffVals(%d),',1:numel(coeffVals));
coeffStr(end) =[]; 
makeCFitStr = sprintf('cfit(fStr, %s);', coeffStr);
curve = eval(makeCFitStr);

%Prepare output: NOTE --> THIS DOES NOT GIVE A SMOOTH SURFACE
yFit_mm = feval(curve_mdl, force1D(xFit_mm,2));
yFit_pxls = yFit_mm*xyRes(2);
xyFit = [force1D(xFit,2), force1D(yFit_pxls,2)];
xyFit_mm = [force1D(xFit_mm, 2), force1D(yFit_mm,2)];
