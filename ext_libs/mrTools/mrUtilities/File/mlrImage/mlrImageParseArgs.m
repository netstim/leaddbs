% mlrImageParseArgs.m
%
%      usage: [imageArgs otherArgs] = mlrImageParseArgs(args)
%         by: justin gardner
%       date: 09/04/11
%    purpose: Used by mlrImageLoad, mlrImageHeaderLoad and mlrVol
%             for parsing input arguments into images and other
%             arguments. Lets user pass in arguments for these
%             functions in the following ways (example is for
%             mlrVol but mlrImageLoad, mlrImageHeaderLoad and
%             any other function that uses this accepts the same
%             types of arguments
%       
%             Filename:
%             mlrVol mprage01.nii
%             mlrVol mprage01.hdr
%             mlrVol mprage01.fid
%
%             ScanNum / groupNum
%
%             v = newView;
%             mlrVol(v,'scanNum=2','groupNum=1')
%
%             Data structure
%
%             d = rand(100,100,200);
%             mlrVol(d)
%
%             Data structure and header
%
%             [d h] = cbiReadNifti('filename.hdr');
%             mlrVol(d,h)
%
%             Pick from canonical volume directory
%
%             mlrVol canonical
%
%             Pick from current directory
%
%             mlrVol .
%
function [imageArgs otherArgs] = mlrImageParseArgs(args)

% init return arguments
imageArgs = {};
otherArgs = {};

% check input arguments
if ~any(nargin == [1])
  help mlrImageParseArgs
  return
end

% filterspecs - this isn't used right now, since to put up
% a dialog where you can select both files and directories
% required going to java and there it was a bit complicated
% how to use the filterspec
filterspec = {'*.hdr;*.nii', 'Nifti Files (*.hdr, *.nii)';'*.sdt;*.edt;','SDT/SPR or EDT/EPR Files (*.sdt, *.spr)'};

% arguments for different filetypes
altArgs = {'orient','xMin','xMax','yMin','yMax','zMin','zMax','swapXY','swapXZ','swapYZ','flipX','flipY','flipZ','shiftX','shiftY','shiftZ','rotateXY','rotateXZ','rotateYZ','interpMethod','applyToHeader','applyToData','rescale'};
loadArgs = {'kspace','movepro','movepss','receiverNum'};

% go through the list looking for image arguments
iArg = 1;nArgs = length(args);
while iArg <= nArgs
  % string arguments that have = signs means they are an argument not a filename
  if isstr(args{iArg}) && ~isempty(strfind(args{iArg},'='))
    % assume that after we see a = sign then all the rest
    % of the arguments are not image filenames
    otherArgs = {args{iArg:end}};
    break;
  elseif isstr(args{iArg})
    % if it is a string that has no extension then see if 
    % it is a special word
    ext = getext(args{iArg});
    if isempty(ext)
      % put on the default extension
      filenameWithDefaultExt = setext(args{iArg},mrGetPref('niftiFileExtension'));
      % and see if that exists
      if isfile(filenameWithDefaultExt) 
	imageArgs{end+1} = filenameWithDefaultExt;
	iArg = iArg+1;
      % or the filename w/out any extension
      elseif isfile(args{iArg}) || (~strcmp(args{iArg},'.') && isdir(args{iArg}))
	disp(sprintf('(mlrImageParseArgs) Unknown filetype for %s (missing extension)',args{iArg}));
	iArg = iArg+1;
      % get from canonical directory, if the name is any one of the following
      elseif any(strcmp({'canonical','volume','volumedirectory','volumedir','voldir'},lower(stripext(args{iArg}))))
	% put up the path/str dialog box
	filename = getPathStrDialog(mrGetPref('volumeDirectory'),'Choose a volume',filterspec,'off');
	% if user hit cancel, then clear out all arguments and return
	if isempty(filename)
	  imageArgs = {};
	  return
	end
	% otherwise use name that the user has supplied a filename
	imageArgs{end+1} = filename;
	iArg = iArg+1;
      % if any of the following then bring up a dialog box
      elseif strcmp('.',args{iArg}) || any(strcmp({'dialog','ask'},lower(stripext(args{iArg}))))
	filename = getPathStrDialog('.','Choose a volume',filterspec,'off');
	% if user hit cancel, then clear out all arguments and return
	if isempty(filename)
	  imageArgs = {};
	  return
	end
	% otherwise use name that the user has supplied
	imageArgs{end+1} = filename;
	iArg = iArg+1;
      else
	% unknown, so just assume all other arguments are parameters
	otherArgs = {args{iArg:end}};
	break;
      end
    else
      % nonempty extension, assume that this is a filename
      imageArgs{end+1}.filename = args{iArg};
      iArg = iArg+1;
      % check for arguments associated with this load filename
      % like movepro or xMin, etc. see above loadArgs and altArgs list
      % keep doing check until we don't pick up any new
      % arguments. This allows user to put the laod and alt
      % arguments in any order they want to
      iArgStart = nan;
      imageArgs{end}.altArgs = {};
      imageArgs{end}.loadArgs = {};
      while (iArgStart ~= iArg)
	iArgStart = iArg;
	[iArg imageArgs{end}.loadArgs] = getOtherArgs(args,iArg,loadArgs,imageArgs{end}.loadArgs);
	[iArg imageArgs{end}.altArgs] = getOtherArgs(args,iArg,altArgs,imageArgs{end}.altArgs);
      end
    end
  elseif isview(args{iArg})
    % if we have a view then collect any additional qualifying
    % arguments-look for scanNum and groupNum args that can get passed to getArgs
    jArg = iArg+1;
    while jArg <= nArgs
      if isstr(args{jArg})
	% see if they contain the string scanNum or groupNum
	if ~isempty(strfind(lower(args{jArg}),'scannum')) || ~isempty(strfind(lower(args{jArg}),'groupnum'))
	  % if they don't have an = sign then the argumnet
	  % comes after so increment jArg by 1
	  if isempty(strfind(args{jArg},'='))
	    jArg = jArg+1;
	  end
	  jArg = jArg+1;
	else
	  break;
	end
      else
	break;
      end
    end
    % now we have all the potential scanNum and groupNum args
    % that can get passed to getArgs, so let getArgs handle it
    scanNum = [];groupNum = [];
    getArgs({args{iArg+1:jArg-1}},{'scanNum=[]','groupNum=[]'});
    % get the correct filename, from the correct scan and group
    filename = viewGet(args{iArg},'tseriespathstr',scanNum,groupNum);
    if isempty(filename)
      disp(sprintf('(mlrImageParseArgs) Could not find tSeries in view'));
    else
      imageArgs{end+1} = filename;
    end
    % set iArg to the next argument
    iArg = jArg;
  % if it is a structure, see if it has a data field
  elseif isstruct(args{iArg})
    % check structures
    if isfield(args{iArg},'filename')
      imageArgs{end+1} = args{iArg};
    elseif isfield(args{iArg},'data')
      imageArgs{end+1} = args{iArg};
    elseif isfield(args{iArg},'d')
      % d is also acceptable, but switch it to a data filed
      imageArgs{end+1} = rmfield(args{iArg},'d');
      imageArgs{end}.data = d;
    elseif mlrImageIsHeader(args{iArg})
      imageArgs{end+1} = args{iArg};
    end
    iArg = iArg+1;
  elseif isnumeric(args{iArg})
    % if it is a data structure, then collect its header
    % if there is one in the next argument
    imageArgs{end+1}.data = double(args{iArg});
    iArg = iArg+1;
    if iArg <= nArgs
      if isstruct(args{iArg})
	imageArgs{end}.h = args{iArg};
	iArg = iArg+1;
      end
    end
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   getPathStrDialog   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%
function pathName = getPathStrDialog(startPath,dialogTitle,filterSpec,multiSelect)

% import the java class
import javax.swing.JFileChooser;

% use current path if not passed in
if (nargin == 0) || isempty(startPath) || strcmp(startPath,'.')
  startPath = pwd;
else
  % snippet of code to handle the ~
  curpwd = pwd;
  cd(startPath);
  startPath = pwd;
  cd(curpwd);
end

% get a java chooser object
jchooser = javaObjectEDT('javax.swing.JFileChooser', startPath);

% set to choose directories
jchooser.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);

% set dialog title
if nargin > 1
  jchooser.setDialogTitle(dialogTitle);
end

% open the dialog
status = jchooser.showOpenDialog([]);

% get status and return
if status == JFileChooser.APPROVE_OPTION
  jFile = jchooser.getSelectedFile();
  pathName = char(jFile.getPath());
else
  pathName = [];
end

%%%%%%%%%%%%%%%%%%%%%%
%%   getOtherArgs   %%
%%%%%%%%%%%%%%%%%%%%%%
function  [iArg loadArgs] = getOtherArgs(args,iArg,argNames,loadArgs)

if nargin < 4,loadArgs = {};end

% go search in args for getArgs style arguments
% that match the argNames list and put those into
% loadArgs
while (iArg<=length(args))
  % check for a string with one of the argument names in it
  if isstr(args{iArg}) && any(strcmp(argNames,strtok(args{iArg},'=')))
    % ok, this is an argument keep it in the list
    loadArgs{end+1} = args{iArg};
    % if there is no equal sign then the next argument has
    % to be passed as well
    if isempty(strfind(args{iArg},'='))    
      iArg = iArg+1;
      % and keep the argument in our list
      if iArg <= length(args)
	loadArgs{end+1} = args{iArg};
      end
    end
    iArg = iArg+1;
  else
    % no argument name, stop
    break;
  end
end

