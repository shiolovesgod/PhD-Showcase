%%This function is a very basic function to be used when pasting a given
%%matrix to Microsoft Excel.  Select the first cell desired for matrix
%%content and simply type the formula into MATLAB.
%NOTE: If there are multiple instances of Excel open, the first instance
%opened will be considered as the active sheet.
%%

function A = mat2xl(Ain, excelFN, isVisible)

switch nargin
    case 1
        isVisible = true; %must be T/F (or 1/0)
        excelFN = '';
    case 2
        isVisible = true; % must be fullfile(pathname, filename) with xls extension
end

%This is to be used with an already open version of Excel

try
    hExcel = actxGetRunningServer('Excel.Application');
catch
    hExcel = actxserver('Excel.Application');
    set(hExcel, 'Visible', isVisible);
    invoke(hExcel.Workbooks,'Add');
end

%Find the cell range
alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';

%Get numeric values for row and column
row = get(hExcel.ActiveCell, 'Row');
col = get(hExcel.ActiveCell, 'Column');

if col <= 26
    colname = alphabet(col);
    
else
    ind1 = floor(col/26);
    ind2 = mod(col,26);
    colname = [alphabet(ind1), alphabet(ind2)];
end

startcellname = [colname,int2str(row)];

%Find the end of the pasting region
A = processInput(Ain);
[len wid] = size(A);

endrow = len+row-1;
endcol = wid+col-1;

if endcol <= 26
    endcolname = alphabet(endcol);
else
    ind3 = floor(endcol/26);
    ind4 = mod(endcol,26);
    endcolname = [alphabet(ind3), alphabet(ind4)];
end

endcellname = [endcolname, int2str(endrow)];

range = [startcellname,':',endcellname];

%Select the Range
Select(Range(hExcel, sprintf('%s',range)));
DataRange = get(hExcel,'Selection');

%Check Contents
warn = 'Yes'; 
cellcontents = get(DataRange, 'Value');

if iscell(cellcontents)
    addit = sum(sum(cellfun(@(x) all(isnan(x)), cellcontents)));
    [a, b] = size(cellcontents);
    isWarnUser = addit ~= a*b;
else
    isWarnUser = ~isnan(cellcontents);
end
    
    if isWarnUser
        warn = questdlg('Are you sure you want to overwrite current cells?','Content','Yes','No','No');
    end

    
    if strcmpi (warn,'Yes')
        %Write the matrix to the range
        set(hExcel.Selection,'Value', A);
    end

 
if ~isempty(excelFN)
    hExcel.ActiveWorkbook.SaveAs(excelFN);
end

if ~isVisible
    hExcel.Quit;
end

function Aout = processInput(Ain)

Aout = {};

if iscell(Ain)
    %for each column, append the matrix to a column
    for iCol = 1:size(Ain,2)
        nRowsInThisCol = 0;
        colStart = size(Aout,2) +1;
        
        for iRow = 1:size(Ain,1)
            thisCell = Ain{iRow,iCol};
            
            isString = ischar(thisCell);
            
            if isempty(thisCell) %added 06/22/2015
                [nRows, nCols] = deal(1);
                thisCell = [];
                isSingular = true; 
            elseif ~isString || size(thisCell,1) > 1
                [nRows, nCols] = size(thisCell);
                isSingular = false;
            else
                [nRows, nCols] = deal(1);
                isSingular = true; 
            end
            
            rowStart = nRowsInThisCol+1;

            if isSingular
                Aout(rowStart:rowStart+nRows-1,colStart:colStart+nCols-1) = {thisCell};
            else
                Aout(rowStart:rowStart+nRows-1,colStart:colStart+nCols-1) = num2cell(thisCell);
            end
            nRowsInThisCol = nRowsInThisCol + nRows;
        end
    end
    
    
else
    Aout = Ain;
end

%%
%CONVERT FROM MATLAB 2 XL 
%col_# = qColNo; col_letter --> qAddress
%    rem = mod(qColNo,nLetters); quo = floor(qColNo/nLetters); 
%     
%      if rem==0
%             rem = nLetters;
%      end
%         
%     qAddress = alphabet(rem);
%     
%     while quo > 0 % do..while loop
%        rem = mod(quo, nLetters);  quo = floor(quo/nLetters); 
%         if rem==0
%             rem = nLetters;
%         end
%         qAddress = [alphabet(rem), qAddress];        
%     end
    