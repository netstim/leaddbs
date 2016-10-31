% mrReadme - make a readme file
%
%      usage: [  ] = mrReadme( session, groups )
%         by: denis schluppeck
%       date: 2008-06-12
%        $Id$:
%     inputs: session, groups (structures usually stored in MLR variable)
%    outputs: 
%
%    purpose: similar functionality to createReadme, but no GUI buttons to press
%             % based on Ben Backus' version for mrLoadRet 2.* and djh's improvements
%
%        e.g: mrReadme(session, groups)
%
function [  ]=mrReadme( session, groups )

if nargin < 2
  help mrReadme
  return
end

% Check if a Readme file already exists
if exist('Readme.txt','file')
    if ~askuser('(mrReadme) Readme.txt already exists. Do you want to create a new file called Readme.txt?')
      % user cancelled, return, but Readme file maybe out of date, so add a note:
      [fid, message] = fopen('Readme.txt', 'a');
      if (~exist('message','var'))
        message='default message';
      end
      
      if fid == -1
	% problems opening readme file for writing?
	mrWarnDlg(message);
	return
      else
	% write the note at the end of the current file:
	fseek(fid, 0, 1);
	fprintf(fid, '%s\n', ['** NB!']);
	fprintf(fid, '%s\n', ['** this readme file may be out of sync with the current mrSession.mat']);
	fclose(fid);
	return
      end
    end
end

% either Readme.txt didn't exist of user wants to create a new one:
% Open Readme.txt
[fid, message] = fopen('Readme.txt', 'w');
if (~exist('message','var'))
  message='default message';
end
if fid == -1
  mrWarnDlg(message);
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

fprintf(fid,'** Readme created: %s **\n\n', datestr(now));

% Close Readme.txt
status = fclose(fid);
if status == -1
    mrWarnDlg(messsage);
    return
end

% Copy notes from an existing file
if askuser('Do you want to append any notes?')
  createReadmeAppendNotes;
end
return

function createReadmeAppendNotes
% function createReadmeAppendNotes
%
% Appends notes file(s) to Readme file.
% ask user if he/she wants to append something?
% Called by mrReadMe.m
%
% ds 2008/06/12 - based on code by djh


% Dialog box to get notes file
pathStr = mlrGetPathStrDialog(pwd,'Select notes to append to Readme.txt','*.*');
if isempty(pathStr)
  return
end

% Read contents of notes file
[fid, message] = fopen(pathStr,'r');
if fid == -1
  mrWarnDlg(messsage);
  return
end
[A,rcount] = fread(fid,inf);
status = fclose(fid);
if status == -1
  mrWarnDlg(messsage);
  return
end

% Append to Readme.txt
[fid, message] = fopen('Readme.txt','a');
if fid == -1
  mrWarnDlg(messsage);
  return
end
wcount = fwrite(fid,A);
if wcount ~= rcount
  warning(['createReadmeAppendNotes failed to append entire notes file, ',pathStr]);
end
status = fclose(fid);
if status == -1
  mrWarnDlg(messsage);
  return
end



return