function output = isThisDuplicate()

%Isolate Process Name
currentDir=mfilename('fullpath');
all = strfind(currentDir,'\'); firstIndex = all(end-1)+1; secondIndex = all(end)-1;
foldername= currentDir(firstIndex:secondIndex);
processName = [foldername,'.exe'];

%Check it against current tasks
[status, task]=dos('tasklist /fo "csv" /nh');

output = isempty(strfind(status,processName));

