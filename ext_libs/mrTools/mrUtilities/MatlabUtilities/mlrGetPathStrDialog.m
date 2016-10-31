function pathstr = getPathStrDialog(startPathStr,title,filterspec,multiselect)
%
% fullPath = getPathStrDialog([startDir],[title],[filterspec],[multiselect])
%
% Calls uigetfile to display a dialog box that allows the user to pick a file.
%
% startDir: starting point for the dialogbox
% title: window title of dialog box, passed on to uigetfile.
% filterspec: passed onto uigetfile. default: '*'
% multiselect: 'on' or 'off', allows user to pick multiple files.
%
% Returns full path string if multiselct is 'off'
% Returns cell array of path strings if multiselect is 'on'
%
% djh, 2/5/2001

if ieNotDefined('startPathStr')
	startPathStr = pwd;
end
if ~exist(startPathStr,'dir')
	startPathStr = pwd;
end
if ieNotDefined('title')
	title = 'Pick a file';
end
if ieNotDefined('filterspec')
	filterspec = {'*.mat','MAT files'; '*.*','All files'};
end
if ieNotDefined('multiselect')
	multiselect = 'off';
end

% save current file position
cDir=pwd;
% move to appropriate directory
chdir(startPathStr);
% call uigetfile
[filename, dirPath] = uigetfile(filterspec,title,'multiselect',multiselect);
% move back to original directory
chdir(cDir);
	
if isequal(filename,0)
	% cancel pressed
	pathstr = [];
elseif iscell(filename)
	% cell array of path strings
	pathstr = cell(size(filename));
	for p = 1:length(filename)
		pathstr{p} = fullfile(dirPath,filename{p});
	end
else
	% return full path string
	pathstr = fullfile(dirPath,filename);
end

% make sure that if it is a multiselect, we always return a cell array
if isstr(pathstr) && strcmp(multiselect,'on')
  tmp = pathstr;
  pathstr = {};
  pathstr{1} = tmp;
end
