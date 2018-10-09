function [uniqueStr, duplicateInd] = uniqueStringList(inputString)
%[uniqueStr, duplicateInd] = uniqueStringList(inputString)
%returns the string without duplicataes and the indices of removed
%duplicates
%input should be a cell string
%case insensitive
%modify this to allow the user to select whether they want unique columns,
%rows or all

nElements = numel(inputStr);
indArray = 1:nElements;
duplicateInd = [];
stringList = inputString;

for iInd = 1:nElements
    %Remove the current element that we are comparing
    currentStr = stringList(iInd);
    remainingList = stringList(~(iInd == indArray));
    
    %Compare the current string to the remaining list
    isDuplicate = cellfun(@(strElement) strcmpi(strElement, stringList),...
        remainingList, 'UniformOutput', false);
    
   duplicateInd = [duplicateInd, indArray(is)]
    
end
