
function [idxMaskAbove, idxMaskBelow,...
    idxAbove, idxBelow] = segmentMatrix(boundary, matrix)
%This function returns the indices of the matrix on either side of a given
%line/boundary (above and below).  The boundary is interpolated linearly for values of the
%boundary not explicitly given
%
%..........................................................................
%inputs:
%..........................................................................
%BOUNDARY: a 2D array containing the x and y values of the boundary
%           surface
%note - any deisred tolerance range must already be included
%MATRIX:  the matrix in which the bondary is drawn (for example an image)
%
%
%..........................................................................
%outputs
%..........................................................................
%idxMaskAbove - logical indexing mask of area of matrix above boundary
%idxMaskBelow - logical indexing maks of area of matrix below boundary
%
%
%idxAbove: a 1D matrix containing the indices in the matrix input about the boundary
%idxBelow: a 1D matrix containing the indices in the matrix input below the boundary
%
%==========================================================================

xIn = boundary(:,1); yIn = boundary(:,2);
[nRows, nCols] = size(matrix);

%make sure border spans from beginning to end of matrix by spreading the
%first and last values of the boundary to the edge of the matrix

if min(xIn) ~=1
    xIn = [1; xIn];
    yIn = [yIn(1); yIn];
end

if max(xIn) ~= nCols
    xIn(end+1) = nCols;
    yIn(end+1) = yIn(end);
end

[~, ia, ~] = unique(xIn, 'stable');
xUnique = xIn(ia); yUnique = yIn(ia);

%interpolate for misssing values of y
xFit = 1:nCols;
yFit = interp1(xUnique,yUnique,xFit,'nearest');

%make sure any y values out of range is brought to matrix boundary
yFit(yFit > nRows) = nRows;
yFit(yFit < 1) = 1;

% xInd = repmat(1:nCols, [nRows,1]);
yInd = repmat((1:nRows)', [1, nCols]);
% boundXInd = repmat(force1D(xFit,1), [nRows,1]);
boundYInd = repmat(force1D(yFit,1), [nRows,1]);

idxMaskAbove = yInd < boundYInd;
idxMaskBelow = yInd > boundYInd;

matIndices = reshape(1:numel(matrix), [nRows,nCols]);
idxAbove = force1D(matIndices(idxMaskAbove));
idxBelow = force1D(matIndices(idxMaskBelow));

%..........................................................................
% %This code takes longer and will give me actual indices
%..........................................................................

%%make sure x and y values don't come out of the matrix
% isInvalid = (xFit > nCols | xFit < 1) | (yFit > nRows | yFit < 1);
% xFit(isInvalid) = []; yFit(isInvalid) = [];
%
% %Convert to a cell array
% iSurfaceFit = num2cell([force1D(xFit), force1D(round(yFit))],2);
%
% %for each column, get the indices
% idxAbove = cellfun(@(idx) sub2ind(imgSize, 1:idx(2),repmat(idx(1),[1,idx(2)])),...
%     iSurfaceFit, 'UniformOutput', false);
% idxAbove = [idxAbove{:}];
%
%
% idxBelow = cellfun(@(idx) sub2ind(imgSize, idx(2):imgRows, repmat(idx(1),[1,idx(2)])),...
%     iSurfaceFit, 'UniformOutput', false);
% idxBelow = [idxBelow{:}];

drawnow; %useless line of code that keeps commented line of code in function

