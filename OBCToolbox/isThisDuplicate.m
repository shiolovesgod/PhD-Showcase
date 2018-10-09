function output = isThisDuplicate(currentDir)

%Isolate Process Name
all = strfind(currentDir,'\'); firstIndex = all(end-1)+1; secondIndex = all(end)-1;
foldername= currentDir(firstIndex:secondIndex);
processName = [foldername,'.exe'];

%Check it against current tasks
[status, task]=dos('tasklist /fo "csv" /nh');

output = numel(strfind(task,processName))>1;

if output
    h = actxserver('WScript.Shell');    
    h.AppActivate('Camera Select');
    h.AppActivate('New');
    %     hDlg = helpdlg('Duplicate Found, program terminating');
    %     waitfor(hDlg);
end


