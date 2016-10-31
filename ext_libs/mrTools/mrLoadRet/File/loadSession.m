function [session, groups, version] = loadSession(dirPathStr)
% function [session, groups, mrLoadRetVersion] = loadSession([dirPathStr])
%
% Loads the mrSESSION and groups structures from mrSESSION.mat.
%
% dirPathStr: directory that contains the mrSESSION.mat file
% defaults to current directory (pwd)
%
% djh, 2/17/2001
% 7/16/02 djh, update mrSESSION to version 3.01
% 11/2004 djh, update to version 4.0

% current version of mrLoadRet
version = mrLoadRetVersion;

% Default dirPathStr = pwd
curDir = pwd;
if ieNotDefined('dirPathStr')
    dirPathStr = curDir;
end
pathStr = fullfile(dirPathStr,'mrSession.mat');

% Load mrSESSION.mat to get session and groups
if exist(pathStr,'file')
    load(pathStr)
    if ieNotDefined('session') 
      disp(sprintf('(loadSession) No session variable in mrSession.mat (probably wrong version file)'));
      session = [];groups = [];version = [];
      return
    end
    
    % Check version numbers for consistency
    if (session.mrLoadRetVersion > version)
      disp('(loadSession) This session was created with mrLoadRet version %i the current mrLoadRet version is %i.',session.mrLoadRetVersion,version);
    end
else
%    mrErrorDlg(['No mrSESSION.mat file in directory: ',dirPathStr]);
    session = [];
    groups = [];
end

% check that all scanParams are valid. This is useful if the
% scanParams optional fields have changed, so that old session
% files can be compatible with changes
for g = 1:length(groups)
  clear newScanParams;
  for s = 1:length(groups(g).scanParams)
    % go through and validate scanParams
    [tf newScanParams(s)] = isscan(groups(g).scanParams(s));
    % check for invalid scna
    if ~tf
      mrWarnDlg(sprintf('(loadSession) Scan %i in group %i is invalid',s,g));
    end
  end
  if exist('newScanParams','var')
    % now set the scan params to this newly validated scanParams. Normally this
    % won't change anything, but if a new field has been added, then the new
    % scanParams will include a default value on this field
    groups(g).scanParams = newScanParams;
  end
end
