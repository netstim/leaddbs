function createReadme(session, groups)
% function createReadme
%
% Creates a Readme.txt text file that describes a scanning session.
%
% djh, 5/2005
% based on Ben Backus' version for mrLoadRet 2.*

if nargin < 2
  help createReadme
  return
end

% Check if a Readme file already exists
if exist('Readme.txt','file')
    questionString = 'Readme.txt already exists. Do you want to continue, which will create a new file called Readme.txt?';
    buttonName = questdlg(questionString, 'Warning', 'Yes', 'No', 'No');
    pause(.1);  % Prevent hanging
    if strcmp(buttonName, 'No')
        return
    end
end

% Open Readme.txt
[fid, message] = fopen('Readme.txt', 'w');
    if (~exist('message','var'))
        message='default message';
    end
if fid == -1
    warndlg(message);
    return
end

% Session identifiers
fprintf(fid, '%s\n', ['Description: ',session.description]);
fprintf(fid, '%s\n', ['Subject: ',session.subject]);
fprintf(fid, '%s\n', ['Operator: ',session.operator]);
fprintf(fid, '%s\n', ['Magnet: ',session.magnet]);
fprintf(fid, '%s\n', ['Coil: ',session.coil]);
fprintf(fid, '%s\n', ['Protocol name: ',session.protocol]);

% Info about each scan in each group
for g = 1:length(groups)
    groupname = groups(g).name;
    scanParams = groups(g).scanParams;

    % Write group name
    fprintf(fid,'\n%s\n',groupname);
    
    % Write cquisition parameters:
    fprintf(fid,'\n%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\n',...
        'scan','description','fileName','voxelSize','dataSize',...
        'totalFrames','junkFrames','nFrames','framePeriod');
    nScans = length(scanParams);
    for iScan = 1:nScans
        fprintf(fid,'%d\t',iScan);
        description = scanParams(iScan).description;
        fprintf(fid,'%s\t',description);
        fileName = scanParams(iScan).fileName;
        fprintf(fid,'%s\t',fileName);
        vsize = scanParams(iScan).voxelSize;
        fprintf(fid,'%s%g%s%g%s%g%s\t','[',vsize(1),' ',vsize(2),' ',vsize(3),']');
        dsize = scanParams(iScan).dataSize;
        fprintf(fid,'%s%d%s%d%s%d%s\t','[',dsize(1),' ',dsize(2),' ',dsize(3),']');
        totalFrames = scanParams(iScan).totalFrames;
        fprintf(fid,'%d\t',totalFrames);
        junkFrames = scanParams(iScan).junkFrames;
        fprintf(fid,'%d\t',junkFrames);
        nFrames = scanParams(iScan).nFrames;
        fprintf(fid,'%d\t',nFrames);
        framePeriod = scanParams(iScan).framePeriod;
        fprintf(fid,'%g\t\n',framePeriod);
    end
    fprintf(fid,'\n\n');
end

% Close Readme.txt
status = fclose(fid);
if status == -1
    warndlg(messsage);
    return
end

% Copy notes from an existing file
createReadmeAppendNotes;




function createReadmeAppendNotes
% function createReadmeAppendNotes
%
% Appends notes file(s) to Readme file. Loops, prompting 
% user for file pathnames until Cancel is selected.
% Called by mrCreateReadme.m
%
% djh, 9/4/01

for i=1:5
    
    % Dialog box to get notes file
    pathStr = mlrGetPathStrDialog(pwd,'Select notes to append to Readme.txt','*.*');
    if isempty(pathStr)
        return
    end
    
    % Read contents of notes file
    [fid, message] = fopen(pathStr,'r');
    if fid == -1
        warndlg(messsage);
        return
    end
    [A,rcount] = fread(fid,inf);
    status = fclose(fid);
    if status == -1
        warndlg(messsage);
        return
    end
    
    % Append to Readme.txt
    [fid, message] = fopen('Readme.txt','a');
    if fid == -1
        warndlg(messsage);
        return
    end
    wcount = fwrite(fid,A);
    if wcount ~= rcount
        warning(['createReadmeAppendNotes failed to append entire notes file, ',pathStr]);
    end
    status = fclose(fid);
    if status == -1
        warndlg(messsage);
        return
    end
end
