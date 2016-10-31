% mr3to4.m
%
%      usage: mr3to4()
%        $Id$
%         by: justin gardner
%       date: 03/12/07
%    purpose: take an mrLoadRet3 directory structure and create
%             a mrLoadRet4 compatibale one. This is done w/out
%             changing anything in the mrLoadRet3 directories and
%             only making appropriate links to files. It will
%             not write over any old files.
%             once this is run once, call with the argument 3 or 4
%             to switch between mrLoadRet3 and mrLoadRet4


function retval = mr3to4(vernum)

% check arguments
if ~any(nargin == [0 1])
  help mr4to3
  return
end

% install mrLoadRet 3.1 paths, if none exist
if ~exist('mlrImageReadNifti')
  mrPaths(3.1);
end

% use for debugging
% if set to one, no actions will be taken
global debugflag;
debugflag = 0;

% file names
sessionFilename = 'mrSESSION.mat';
oldSessionFilename = 'mrSession4.mat';

% first make sure we have a session file
if ~isfile(sessionFilename)
  disp('(mr3to4) Could not find mrSession.mat file');
  return
end

% now load it, and make sure that it is the correct version
m = load(sessionFilename);
if isfield(m,'session') && isfield(m,'groups') && (m.session.mrLoadRetVersion >= 4)
  currentMrSessionVersion = 4;
elseif isfield(m,'dataTYPES') && isfield(m,'mrSESSION') && (m.mrSESSION.mrLoadRetVersion >= 3)
  currentMrSessionVersion = 3;
else
  disp('(mr3to4) Unknown mrLoadRet version');
  return
end

% if called with one argument, all we need to do is swap
% the mrSession file
if (nargin == 1)
  if isfile('mrSession4.mat') && isfile('mrSESSION3.mat')
    if (currentMrSessionVersion == 3) && (vernum == 4)
	disp(sprintf('Swapping in mrSession 4 file'));
	mysystem('mv -f mrSession.mat mrSESSION3.mat');
	mysystem('cp mrSession4.mat mrSession.mat');
    elseif (currentMrSessionVersion == 4) && (vernum == 3)
	disp(sprintf('Swapping in mrSession 3.1 file'));
	mysystem('mv -f mrSession.mat mrSESSION4.mat');
	mysystem('cp mrSession3.mat mrSESSION.mat');
    end
  else
    disp('(mr3to4) Run mr4to3 with no arguments first');
  end
  return
end

% make sure we have version 3
if currentMrSessionVersion ~= 3
  disp('(mr3to4) This is not a valid mrLoadRet 3 directory');
  return
end

% print out information on raw group
printBlockBegin('PFiles');
scanParams = m.dataTYPES(1).scanParams;
functionals = m.mrSESSION.functionals;
for i = 1:length(functionals)
  disp(sprintf('%i: %s [%s] tr=%0.3f',i,functionals(i).PfileName,num2str(functionals(i).effectiveResolution),functionals(i).framePeriod));
end
printBlockEnd;

% look in anatomy directory for anatomy files
if isdir('Raw/Anatomy/Inplane')
  anatdir = dir('Raw/Anatomy/Inplane/*.img');
  if length(anatdir)
    printBlockBegin('Anatomy files');
    for i = 1:length(anatdir)
      disp(sprintf('%i: %s',i,anatdir(i).name));
    end
  else
    disp(sprintf('(mr3to4) No anatomy found'));
  end
  printBlockEnd;
end

% check for an old mrSESSION (4) file
moveOldSessionFile = 0;
if isfile(oldSessionFilename)
  moveOldSessionFile = 1;
  printBlockBegin('Old mrSession.mat');
  newMrSessionName = sprintf('%s_%s.mat',stripext(oldSessionFilename),datestr(now,'yymmdd_HH:MM:SS'));
  disp(sprintf('%s will be moved to -> %s',oldSessionFilename,newMrSessionName));
  printBlockEnd;
end

% make sure everything is ok to go
if ~askuser('Are these the correct data files?')
  return
end

% ok, make the appropriate directories
printBlockBegin('Making directories');
dirs = {'MotionComp' 'MotionComp/TSeries' 'Raw/TSeries'};
for j = 1:length(dirs)
  if (~isdir(sprintf('%s',dirs{j}))) || isempty(dir(dirs{j}))
    disp(sprintf('(mr3to4) Making directory: %s',dirs{j}));
    mysystem(sprintf('mkdir %s',dirs{j}))
  end
end
printBlockEnd;

% link anatomy directory
if ~isdir('Anatomy')
  disp(sprintf('(mr3to4) Linking anatomy directory'));
  if ~debugflag
    mysystem('ln -s Raw/Anatomy/Inplane Anatomy');
  end
end

% now make links in Raw/Pfiles to the motion comp data
filenames = makeLinks(m.mrSESSION.functionals,'../../Raw/Pfiles/','MotionComp/TSeries','../../Raw/Pfiles_preMC');

preMCScanParams = m.mrSESSION.functionals;
for i = 1:length(filenames)
  % strip off _mcf.img and replace
  preMCScanParams(i).PfileName = sprintf('%s.img',filenames{i}(1:end-8));
end

% now make links to Raw/Pfiles_preMC
filenamesRaw = makeLinks(preMCScanParams,'../../Raw/Pfiles_preMC','Raw/TSeries');
printBlockEnd;

% move the old mr session
printBlockBegin(sprintf('Making new %s (4) file',oldSessionFilename));

% deal with old mrSession file
if moveOldSessionFile
  disp(sprintf('(mr3to4) Moving %s to %s',oldSessionFilename,newMrSessionName));
  mysystem(sprintf('mv %s %s',oldSessionFilename,newMrSessionName));
end

% now create new mrSESSION variable
m.session.mrLoadRetVersion = 4.5;
m.session.description = m.mrSESSION.sessionCode;
m.session.subject = m.mrSESSION.subject;
m.session.operator = 'XX';
m.session.magnet = 'Allegra 3T';
m.session.coil = 'XXX';
m.session.protocol = 'XXX';

% make raw group
m.groups(1).name = 'Raw';
for i = 1:length(preMCScanParams)
  newScanParams(i).dataSize = [preMCScanParams(i).fullSize length(preMCScanParams(i).slices)];
  newScanParams(i).description = m.dataTYPES(1).scanParams(i).annotation;
  newScanParams(i).fileName = preMCScanParams(i).PfileName;
  newScanParams(i).fileType = 'Nifti';
  newScanParams(i).framePeriod = preMCScanParams(i).framePeriod;
  newScanParams(i).junkFrames = preMCScanParams(i).junkFirstFrames;
  newScanParams(i).nFrames = preMCScanParams(i).nFrames;
  newScanParams(i).niftiHdr = preMCScanParams(i).analyzeInfo;
  newScanParams(i).totalFrames = preMCScanParams(i).totalFrames;
  newScanParams(i).voxelSize = preMCScanParams(i).effectiveResolution;
  newScanParams(i).originalFileName = [];
  newScanParams(i).originalGroupName = [];
  [tf thisScanParams] = isscan(newScanParams(i));
  if isscan(thisScanParams);
    m.groups(1).scanParams(i) = orderfields(thisScanParams);
  else
    keyboard
  end
  % get the stimfilename
  if isfield(m.mrSESSION,'doer')
    expnum = [];
    % search for matching file name, since the numbering is different for doer/mrloadret
    for j = 1:length(m.mrSESSION.doer.exp)
      if isfield(m.mrSESSION.doer.exp{j},'filename') && strcmp(stripext(m.mrSESSION.doer.exp{j}.filename),stripext(functionals(i).PfileName))
	expnum = j;
      end
    end
    if ~isempty(expnum)
      [pathstr m.groups(1).auxParams(i).stimFileName] = fileparts(m.mrSESSION.doer.exp{expnum}.stimfilename);
      m.groups(1).auxParams(i).stimFileName = sprintf('%s.mat',stripext(m.groups(1).auxParams(i).stimFileName));
    end
  end
end

% make motion comp group
m.groups(2).name = 'MotionComp';
for i = 1:length(filenames)
  newScanParams(i).originalFileName{1} = m.groups(1).scanParams(i).fileName;
  newScanParams(i).junkFrames = 0;
  newScanParams(i).originalGroupName{1} = 'Raw';
  newScanParams(i).fileName = filenames{i};
  [tf thisScanParams] = isscan(newScanParams(i));
  if isscan(thisScanParams)
    m.groups(2).scanParams(i) = orderfields(thisScanParams);
  else
    keyboard
  end
  m.groups(2).auxParams(i).auxParams = 1;
end

if ~debugflag
  disp(sprintf('(mr3to4) Saving mrSESSION.mat'));
  session = m.session;
  groups = m.groups;
  eval(sprintf('save %s session groups -V6',oldSessionFilename));
  mysystem('cp mrSESSION.mat mrSESSION3.mat');
end
printBlockEnd
mr3to4(4);

%%%%%%%%%%%%%%%%%%%%%%%%
% make links to files
%%%%%%%%%%%%%%%%%%%%%%%%
function destFilenames = makeLinks(scanParams,sourcedir,destdir,sourceHeaderDir)

if ~exist('sourceHeaderDir'),sourceHeaderDir = '';,end

origDir = pwd;
fprintf(sprintf('Current directory is %s \n', origDir));
fprintf(sprintf('Entering %s \n', destdir));
cd(destdir);

scanParams = scanParams;
printBlockBegin(sprintf('Linking %s files',sourcedir));

for i = 1:length(scanParams)
  % get sourcename,
  sourceFilename = fullfile(sourcedir,scanParams(i).PfileName);
  destFilenames{i} = scanParams(i).PfileName;
  destFilename = fullfile(destdir,scanParams(i).PfileName);

  % check to see if the img file is already there
  if isfile(destFilename)
    disp(sprintf('%s already exists',destFilename));
  elseif ~isfile(sourceFilename)
    disp(sprintf('(mr3to4) Could not find original file %s',sourceFilename));
  else
    disp(sprintf('Linking %s to %s',sourceFilename,destFilename));
    mysystem(sprintf('ln -s %s .',sourceFilename));
  end
  % now link the hdr file
  destFilename = sprintf('%s.hdr',stripext(destFilename));
  sourceFilename = sprintf('%s.hdr',stripext(sourceFilename));
  % if we have a source header dir that mean to copy the
  % matching header from the source dir instead if it exists
  if ~isempty(sourceHeaderDir)
    revertSourceFilename = sourceFilename;
    % the filename will have the trailing _mcf.hdr taken off
    if ~isempty(findstr('_mcf',sourceFilename))
      sourceFilename = sprintf('%s.hdr',stripext(stripext(getLastDir(sourceFilename)),'_'));
    else
      sourceFilename = sprintf('%s.hdr',stripext(getLastDir(sourceFilename)));
    end
    sourceFilename = fullfile(sourceHeaderDir,sourceFilename);
    % if the hdr doesn't exist, then go back to the one we were
    % going to use
    if ~isfile(sourceFilename)
      sourceFilename = revertSourceFilename;
    end
  end
  if isfile(destFilename)
    disp(sprintf('%s already exists',destFilename));
  elseif ~isfile(sourceFilename)
    disp(sprintf('(mr3to4) Could not find original file %s',sourceFilename));
  else
    disp(sprintf('Linking %s to %s',sourceFilename,destFilename));
    mysystem(sprintf('ln -s %s %s',sourceFilename,getLastDir(destFilename)));
  end
  % now link the dicom file
  destFilename = sprintf('%s-header.txt',stripext(destFilename));
  sourceFilename = sprintf('%s-header.txt',stripext(sourceFilename));
  if isfile(destFilename)
    disp(sprintf('%s already exists',destFilename));
  else
    disp(sprintf('Linking %s to %s',sourceFilename,destFilename));
    mysystem(sprintf('ln -s %s .',sourceFilename));
  end
end
printBlockEnd;

fprintf(sprintf('Returning to %s \n', origDir));
cd(origDir)

% swap session to 4
mr3to4(4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% draw a line separator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function printBlockBegin(blockname)

disp(sprintf(blockname));
disp('===================================================================');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% draw a line separator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function printBlockEnd

mrDisp(sprintf('\n'));
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function used for debugging
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mysystem(str)

global debugflag;

% for debugging just print out commands
% and don't actually run them
if (debugflag)
  disp(str);
else
  system(str);
end




function Xform = myInplane2VolXform(rot,trans,scaleFac)
%
% Xform = inplane2VolXform(rot,trans,scaleFac)
%
% Returns 4x4 homogeneous tranform that tranforms from inplane to
% volume.
%
% djh/gmb, '97
%
% Modification:
% - Flip first 2 rows and cols so that it deals with
% (y,x,z) coords instead of (x,y,z).  DJH, 7/98.

A=diag(scaleFac(2,:))*rot*diag(1./scaleFac(1,:));
b = (scaleFac(2,:).*trans)';

Xform = zeros(4,4);
Xform(1:3,1:3)=A;
Xform(1:3,4)=b;
Xform(4,4)=1;

Xform([1 2],:) = Xform([2 1],:);
Xform(:,[1 2]) = Xform(:,[2 1]);
    
    
