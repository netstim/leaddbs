% mr4to3.m
%
%      usage: mr4to3()
%        $Id$
%         by: justin gardner
%       date: 03/12/07
%    purpose: take an mrLoadRet4 directory structure and create
%             a mrLoadRet3.1 compatibale one. This is done w/out
%             changing anything in the mrLoadRet4 directories and
%             only making appropriate links to files. It will
%             not write over any old files.
%             once this is run once, call with the argument 3 or 4
%             to switch between mrLoadRet3 and mrLoadRet4


function retval = mr4to3(vernum)

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
sessionFilename = 'mrSession.mat';
bestRotFilename = 'bestRotVol.mat';
oldSessionFilename = 'mrSESSION3.mat';

% first make sure we have a session file
if ~isfile(sessionFilename)
  disp('(mr4to3) Could not find mrSession.mat file');
  return
end

% now load it, and make sure that it is the correct version
m = load(sessionFilename);
if isfield(m,'session') && isfield(m,'groups') && (m.session.mrLoadRetVersion >= 4)
  currentMrSessionVersion = 4;
elseif isfield(m,'dataTYPES') && isfield(m,'mrSESSION') && (m.mrSESSION.mrLoadRetVersion >= 3)
  currentMrSessionVersion = 3;
else
  disp('(mr4to3) Unknown mrLoadRet version');
  keyboard
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
    disp('(mr4to3) Run mr4to3 with no arguments first');
  end
  return
end

% make sure we have version 4
if currentMrSessionVersion ~= 4
  disp('(mr4to3) This is not a valid mrLoadRet 4 directory');
  return
end
% Check that we have necessary groups
rawGroup = 0;motionCompGroup = 0;averagesGroup = 0;
for i = 1:length(m.groups)
    if strcmp(m.groups(i).name,'Raw')
        rawGroup = i;
    elseif strcmp(m.groups(i).name,'MotionComp')
        motionCompGroup = i;
    elseif strcmp(m.groups(i).name, 'Averages')
        averagesGroup = i;
    end
end

% if there is no raw group, then quit
if ~rawGroup
  disp('(mr4to3) Could not find raw data group');
  return
end

% print out information on raw group
printBlockBegin('Raw files');
scanParams = m.groups(rawGroup).scanParams;
for i = 1:length(scanParams)
  disp(sprintf('%i: %s (%s) [%s] tr=%0.3f',i,scanParams(i).fileName,scanParams(i).description,num2str(scanParams(i).dataSize),scanParams(i).framePeriod));
end
printBlockEnd;

% print out information on MotionComp group
if motionCompGroup
  printBlockBegin('Motion corrected files');
  scanParams = m.groups(motionCompGroup).scanParams;
  for i = 1:length(scanParams)
    disp(sprintf('%i: %s (%s) [%s] tr=%0.3f',i,scanParams(i).fileName,scanParams(i).description,num2str(scanParams(i).dataSize),scanParams(i).framePeriod));
  end
  printBlockEnd;
end

% print out information on Averages group
if averagesGroup
  printBlockBegin('Averaged files');
  scanParams = m.groups(averagesGroup).scanParams;
  for i = 1:length(scanParams)
    disp(sprintf('%i: %s (%s) [%s] tr=%0.3f',i,scanParams(i).fileName,scanParams(i).description,num2str(scanParams(i).dataSize),scanParams(i).framePeriod));
  end
  printBlockEnd;
end


% look in anatomy directory for anatomy files
if isdir('Anatomy')
  anatdir = dir('Anatomy/*.img');
  if length(anatdir)
    printBlockBegin('Anatomy files');
    for i = 1:length(anatdir)
      disp(sprintf('%i: %s',i,anatdir(i).name));
    end
  else
    disp(sprintf('(mr4to3) No anatomy found'));
    return
  end
  % choose an anatomy
  if length(anatdir)>1
    anatnum = getnum(sprintf('Choose which anatomy to use (1:%i)',length(anatdir)),1:length(anatdir));
  else
    anatnum = 1;
  end
  printBlockEnd;
end

% load the alignment file
alignment = [];
if isfile(bestRotFilename)
  printBlockBegin(sprintf('Alignment file (%s)',bestRotFilename));
  load(bestRotFilename);
  % make the alignment
  if exist('xform','var')
    alignment = xform;
    disp(num2str(alignment));
  elseif exist('rot','var') & exist('trans','var') & exist('scaleFac','var')
    alignment = myInplane2VolXform(rot,trans,scaleFac);
    disp(num2str(alignment));
  else
    disp('(mr4to3) Invalid alignment file.')
  end
  printBlockEnd;
end

% check for an old mrSESSION (3.1) file
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
%dirs = {'Raw/Anatomy' 'Raw/Anatomy/Inplane' 'Raw/Pfiles' 'Raw/Pfiles_preMC' 'Inplane'};
dirs = {'Raw/Anatomy' 'Raw/Anatomy/Inplane' 'Raw/Pfiles' 'Raw/Pfiles_preMC' 'Inplane','Inplane/Original','Inplane/Original/TSeries','Inplane/Original/TSeries/Analyze'};

if averagesGroup
    dirs{end+1} = 'Inplane/Averages';
    dirs{end+1} = 'Inplane/Averages/TSeries';
    dirs{end+1} = 'Inplane/Averages/TSeries/Analyze';
end

for j = 1:length(dirs)
  if (~isdir(sprintf('%s',dirs{j})))
    disp(sprintf('(mr4to3) Making directory: %s',dirs{j}));
    mysystem(sprintf('mkdir %s',dirs{j}))
  end
end
printBlockEnd;

% now make links in Raw/Pfiles_preMC for the raw data
filenames{rawGroup} = makeLinks(m.groups(rawGroup).scanParams,'../../Raw','Raw/Pfiles_preMC');

% now make links in Raw/Pfiles to the motion comp data
if motionCompGroup
  filenames{motionCompGroup} = makeLinks(m.groups(motionCompGroup).scanParams,'../../MotionComp','Raw/Pfiles');
  % we will make the originals from motionCompGroup
  originalGroup = motionCompGroup;
  % also make link into tseries
  makeLinks(m.groups(motionCompGroup).scanParams,'../../../../MotionComp','Inplane/Original/TSeries/Analyze');
else
  % use the un motion compd data if motion comp is not done
  makeLinks(m.groups(rawGroup).scanParams,'../../Raw','Raw/Pfiles');
  % also make link into tseries
  makeLinks(m.groups(rawGroup).scanParams,'../../../../Raw','Inplane/Original/TSeries/Analyze');
  originalGroup = rawGroup;
end

if averagesGroup
    %filenames{averagesGroup} = makeLinks(m.groups(averagesGroup).scanParams,'../../Averages','Raw/Averages');
    filenames{averagesGroup} = makeLinks(m.groups(averagesGroup).scanParams,'../../../../Averages','Inplane/Averages/TSeries/Analyze');
end    
printBlockEnd;

if ~isfile('Inplane/Inplane.hdr') && anatnum
  printBlockBegin('Saving anatomy file');
  [anatimg anathdr] = mlrImageReadNifti(fullfile('Anatomy',anatdir(anatnum).name));
  cbiWriteNifti('Inplane/Inplane.hdr',anatimg,anathdr);
  printBlockEnd;
else
  [anatimg anathdr] = mlrImageReadNifti('Inplane/Inplane.hdr');
end

% move the old mr session
printBlockBegin(sprintf('Making new %s (3.1) file',oldSessionFilename));

% deal with old mrSession file
if moveOldSessionFile
  disp(sprintf('(mr4to3) Moving %s to %s',oldSessionFilename,newMrSessionName));
  mysystem(sprintf('mv %s %s',oldSessionFilename,newMrSessionName));
end


% now create new mrSESSION variable
m.mrSESSION.mrLoadRetVersion = 3.1;
m.mrSESSION.sessionCode = m.session.description;
m.mrSESSION.description = m.session.description;
m.mrSESSION.subject = m.session.subject;
m.mrSESSION.examNum = [];

% set inplane info
m.mrSESSION.inplanes.analyzeInfo = anathdr;
m.mrSESSION.inplanes.fullSize = anathdr.dim(2:3)';
m.mrSESSION.inplanes.voxelSize = anathdr.pixdim(2:4)';
m.mrSESSION.inplanes.FOV = anathdr.dim(2:3)' .* anathdr.pixdim(2:3)';
m.mrSESSION.inplanes.nSlices = anathdr.dim(4);
m.mrSESSION.inplanes.origin = [0 0 0]; % make up this one
m.mrSESSION.inplanes.spacing = 0; % make up this one.
m.mrSESSION.inplanes.examNum = [];
m.mrSESSION.inplanes.crop = [];
m.mrSESSION.inplanes.cropSize = anathdr.dim(2:3)';

m.mrSESSION.alignment = alignment;

% make the data types structure and fill out 
% the functionals field in mrSESSION
m.dataTYPES(1).name = 'Original';
scanParams = m.groups(originalGroup).scanParams;
for i = 1:length(scanParams)
  % scan params
  m.dataTYPES(1).scanParams(i).annotation = scanParams(i).description;
  m.dataTYPES(1).scanParams(i).nFrames = scanParams(i).nFrames;
  m.dataTYPES(1).scanParams(i).framePeriod = scanParams(i).framePeriod;
  m.dataTYPES(1).scanParams(i).slices = [1:scanParams(i).dataSize(3)];
  m.dataTYPES(1).scanParams(i).cropSize = scanParams(i).dataSize(1:2);

  % blocked analysis (just make up)
  m.dataTYPES(1).blockedAnalysisParams(i).blockedAnalysis = 1;
  m.dataTYPES(1).blockedAnalysisParams(i).detrend = 1;
  m.dataTYPES(1).blockedAnalysisParams(i).inhomoCorrect = 1;
  m.dataTYPES(1).blockedAnalysisParams(i).temporalNormalization = 0;
  m.dataTYPES(1).blockedAnalysisParams(i).nCycles = 10;

  % event analysis set to 0
  m.dataTYPES(1).eventAnalysisParams(i).eventAnalysis = 0;

  % set functional info in mrSESSION
  m.mrSESSION.functionals(i).analyzeInfo = scanParams(i).niftiHdr;
  m.mrSESSION.functionals(i).PfileName = filenames{originalGroup}{i};
  m.mrSESSION.functionals(i).voxelSize = scanParams(i).voxelSize;
  m.mrSESSION.functionals(i).totalFrames = scanParams(i).totalFrames;
  m.mrSESSION.functionals(i).junkFirstFrames = scanParams(i).junkFrames;
  m.mrSESSION.functionals(i).nFrames = scanParams(i).nFrames;
  m.mrSESSION.functionals(i).slices = 1:scanParams(i).dataSize(3);
  m.mrSESSION.functionals(i).fullSize = scanParams(i).dataSize(1:2);
  m.mrSESSION.functionals(i).cropSize = scanParams(i).dataSize(1:2);
  m.mrSESSION.functionals(i).crop = [];
  m.mrSESSION.functionals(i).effectiveResolution = scanParams(i).voxelSize;
  m.mrSESSION.functionals(i).framePeriod = scanParams(i).framePeriod;
  m.mrSESSION.functionals(i).reconParams = [];
  m.mrSESSION.functionals(i).rawFileType = 'Analyze';
end

% make the data types structure and fill out 
% the functionals field in mrSESSION
if averagesGroup
  m.dataTYPES(2).name = 'Averages';
  scanParams = m.groups(averagesGroup).scanParams;
  for i = 1:length(scanParams)
    % scan params
    m.dataTYPES(2).scanParams(i).annotation = scanParams(i).description;
    m.dataTYPES(2).scanParams(i).nFrames = scanParams(i).nFrames;
    m.dataTYPES(2).scanParams(i).framePeriod = scanParams(i).framePeriod;
    m.dataTYPES(2).scanParams(i).slices = [1:scanParams(i).dataSize(3)];
    m.dataTYPES(2).scanParams(i).cropSize = scanParams(i).dataSize(1:2);

    % blocked analysis (just make up)
    m.dataTYPES(2).blockedAnalysisParams(i).blockedAnalysis = 1;
    m.dataTYPES(2).blockedAnalysisParams(i).detrend = 1;
    m.dataTYPES(2).blockedAnalysisParams(i).inhomoCorrect = 1;
    m.dataTYPES(2).blockedAnalysisParams(i).temporalNormalization = 0;
    m.dataTYPES(2).blockedAnalysisParams(i).nCycles = 10;

    % event analysis set to 0
    m.dataTYPES(2).eventAnalysisParams(i).eventAnalysis = 0;

  end
end




if ~debugflag
  disp(sprintf('(mr4to3) Saving mrSESSION.mat'));
  mrSESSION = m.mrSESSION;
  dataTYPES = m.dataTYPES;
  eval(sprintf('save %s mrSESSION dataTYPES -V6',oldSessionFilename));
end
printBlockEnd

% switch the mrSessions
mr4to3(3);

%%%%%%%%%%%%%%%%%%%%%%%%
% make links to files
%%%%%%%%%%%%%%%%%%%%%%%%
function destFilenames = makeLinks(scanParams,sourcedir,destdir)

origDir = pwd;
fprintf(sprintf('Current directory is %s \n', origDir));
fprintf(sprintf('Entering %s \n', destdir));
cd(destdir);

scanParams = scanParams;
printBlockBegin(sprintf('Linking %s files',sourcedir));

for i = 1:length(scanParams)
  % get sourcename,
  sourceFilename = fullfile(sourcedir,'TSeries',scanParams(i).fileName);
  % if there is an original filename, then use that
  if isfield(scanParams(i),'originalFileName') && ...
	(length(scanParams(i).originalFileName) == 1)
    % append the groupName
    destFilenames{i} = sprintf('%s_%s.img',stripext(scanParams(i).originalFileName{1}),sourcedir);
    destFilename = fullfile(destdir,destFilenames{i});
    destFilename = strcat('Scan', num2str(i), '.img');
  % otherwise just use the filename
  else
      destFilenames{i} = scanParams(i).fileName;
      %destFilename = fullfile(destdir,scanParams(i).fileName);
      destFilename = strcat('Scan', num2str(i), '.img');
  end
  % check to see if the img file is already there
  if isfile(destFilename)
    disp(sprintf('%s already exists',destFilename));
  elseif ~isfile(sourceFilename)
    disp(sprintf('(mr4to3) Could not find original file %s',sourceFilename));
  else
    disp(sprintf('Linking %s to %s',sourceFilename,destFilename));
    mysystem(sprintf('ln -s %s %s',sourceFilename, destFilename));
  end
  % now move the hdr file
  destFilename = sprintf('%s.hdr',stripext(destFilename));
  sourceFilename = sprintf('%s.hdr',stripext(sourceFilename));
  if isfile(destFilename)
    disp(sprintf('%s already exists',destFilename));
  elseif ~isfile(sourceFilename)
    disp(sprintf('(mr4to3) Could not find original file %s',sourceFilename));
  else
    disp(sprintf('Linking %s to %s',sourceFilename,destFilename));
    mysystem(sprintf('ln -s %s %s',sourceFilename, destFilename));
  end
end
printBlockEnd;

fprintf(sprintf('Returning to %s \n', origDir));
cd(origDir)

if ~isfile('mrSession4.mat')
  mysystem('cp mrSession.mat mrSession4.mat');
end

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
    
    
