function cellOutput = ypStructOut(pgNo)

category = 'Grocery Stores';
zipCode = 33186;
sortorder = 'distance';

%the expression for each field
expBname ='<h3 class="business-name fn org">.*?<a.*?>(?<bname>.)*?</a>.*?</h3>';
expAdr = '<span class="street-address">.*?(?<adr>[\S ]+).*?</span>';
expCity = '<span class="locality">(?<city>.)*?</span>';
expState = '<span class="region">(?<state>.)*?</span>';
expZip = '<span class="postal-code">(?<zipCode>.)*?</span>';
expMiles = '<div class="distance">.*?(?<distance>.)*?</div>';
expPhone = '<span class="business-phone phone">(?<phone>.)*?</span>';
allExp = {expBname, expAdr, expCity, expState, expZip, expPhone, expMiles};
expTitle = '<title>(?<pgTitle>.)*?</title>';

%Initialize page number if it is not provided
if nargin < 1
    isPgNoDefined = false;
    pgNo = 1;
else
    isPgNoDefined = true;
end

pgNoIndx = 1;
thisPg = pgNo(pgNoIndx);
allStoresStruct = struct([]);
fprintf(1,'\nYellow Pages Results: %s near %d\n\n', category, zipCode);

%Get the file from online (for a given page
pgResult = urlread(createYPurl(zipCode, category, thisPg, sortorder));
[~,S] = regexp(pgResult, expTitle,'tokens','names');

checkPage = @(S) isempty(strfind(S.pgTitle,'No Matches Found'));


while checkPage(S)
    
    %This gets the information for each store and puts it in a cell array
    allStores = regexp(pgResult, '<div class="info">(?<tagvalue>.*?)</div>','match');

    for i = 1: numel(allStores)
        
        iStore = allStores{i};
        
        %This generates the outputs for each store
        [~, structOut] = regexp(iStore, allExp,'tokens','names');
        
        
        if ~any(cellfun(@isempty, structOut))
            
            %A structure with the information for each store
            thisStore = cell2struct(cellfun(@struct2cell,structOut)',cellfun(@fieldnames,structOut));
            
            
            %FORMATTING OF PARTICULAR FIELDS
            %Replace &amp; with an ampersand only (&) in the businessname
            thisStore.bname = regexprep(thisStore.bname, '&amp;', '&');
            
            %Make sure the phone number is numbers only
            thisStore.phone = str2num(char(regexp(thisStore.phone,'\d','match'))');
            
            %Make sure the distance is miles only (no units)
            newDistance = cell2mat(regexp(thisStore.distance, ('\d*\.\d*|\d*'),'match'));
            
            if ~isempty(newDistance)
                thisStore.distance = newDistance;
            end
            
            %Create a full address
            thisStore.fullAddress = [thisStore.adr, ' ',thisStore.city,...
                ' ',thisStore.state, ' ',thisStore.zipCode];
            
            if isempty(allStoresStruct)
                allStoresStruct = thisStore;
            else
                allStoresStruct(end+1) = thisStore; 
            end
        end
    end
    
    fprintf(1,'Page %d complete.\n', thisPg);
    
    %get the next page
    pgNoIndx = pgNoIndx + 1;
    
    if isPgNoDefined
        try
            thisPg = pgNo(pgNoIndx);
        catch
            break
        end
    else
        thisPg = pgNoIndx;
    end
    
    pgResult = urlread(createYPurl(33186, 'Grocery Stores', thisPg, sortorder));
    [~,S] = regexp(pgResult, expTitle,'tokens','names');
    
end


%remove duplicate fields
phoneNumbers = [allStoresStruct.phone];
isDuplicate = [false, diff(phoneNumbers) == 0];
allStoresStruct(isDuplicate) = [];


%Put it in a cell so that you can paste to excel
nEntries = numel(allStoresStruct);
nFields = numel(fieldnames(allStoresStruct));
cellOutput = cell(nEntries + 1, nFields);
colHeaders = fieldnames(allStoresStruct);


for fieldIndx = 1:nFields
   iField =  colHeaders{fieldIndx};
   cellOutput(1,fieldIndx) = {iField};
   cellOutput(2:end, fieldIndx) = cellfun(@(x) allStoresStruct(x).(iField),...
       num2cell(1:nEntries),'UniformOutput', false);
end


%newLineChars = strfind(fileX, char(10));