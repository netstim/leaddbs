% mlrAnatDBPut.m
%
%        $Id:$ 
%      usage: mlrAnatDBPut(subjectID,filePath,fileType,<verbose=0>,<comments=[]>)
%         by: justin gardner
%       date: 06/22/15
%    purpose: Puts files into repo. The filePath should point to a directory outside of mlrAnatDB and
%             that file or directory will be copied into an appropriate place into mlrAnatDBm added to the
%             repo, committed and pushed (if the user accepts)
%
%             Typically you specify what kind of file you have (roi, freesurfer, baseAnat, localizer) and
%             this will put it into the correct directory
% 
%             e.g. To put mlr surfaces from a view
%             mlrAnatDBPut(25,v,'mlrBaseAnat');
%
%             e.g. To put a directory of rois
%             mlrAnatDBPut(25,'~/data/rois','roi');
%
%             e.g. To put a freesurfer directory
%             mlrAnatDBPut(25,'~/data/freesurfer','freesurfer');
%
%             Valid types can be found in the code below: roi, freesurfer, baseAnat, freesurfer, canonical, localizer
%
function tf = mlrAnatDBPut(subjectID,filePath,fileType,varargin)

tf = false;
% check arguments
if nargin < 2
  help mlrAnatDBPut;
  return;
end

% format the subject id
subjectID = mlrAnatDBSubjectID(subjectID);

% and get arguments
getArgs(varargin,{'fileDir=[]','largefiles=[]','verbose=0','comments=[]','freesurfer=[]','useHardLinks=1'});

% check file
if ~iscell(filePath) && ~isview(filePath) && (~isfile(filePath) && ~isdir(filePath))
  disp(sprintf('(mlrAnatDBPut) Could not find: %s',filePath));
  return
end

% if there are no comments then ask user for comments
if isempty(comments)
  options.Resize='on';
  comments = inputdlg('Enter commit comments','Commit comments',[1 100],{''},options);
  if ~isempty(comments)
    comments = comments{1};
  else
    comments = '';
  end
end

v = [];
% figure out what path to put it under and whether it is a large file or not
if isempty(fileDir)
  switch lower(fileType)
   % Add ROIs to DB
   case {'roi','rois'}
    fileType = 'rois';
    fileDir = 'mlrROIs';
    if isempty(largefiles),largefiles = false;end
    if isview(filePath)
      v = filePath;
      % if a view then have user select which bases to get
      [filePath filePathNiftiROIs filePathScreenShots]= getROIFilenames(v,subjectID);
      % put the niftis and screen shots into database
      mlrAnatDBPut(subjectID,filePathNiftiROIs,'niftiROIs','comments',comments);
      mlrAnatDBPut(subjectID,filePathScreenShots,'screenShots','comments',comments);
    elseif isdir(filePath)
      filePath = getFilenames(filePath,'*.mat');
    end
   % Add Nifti ROIs to DB
   case {'niftiroi','niftirois'}
    fileType = 'niftirois';
    fileDir = 'niftiROIs';
    if isempty(largefiles),largefiles = false;end
    if ~iscell(filePath) && isdir(filePath)
      filePath = getFilenames(filePath,'*.mat');
    end
   % Add Screen shots to DB
   case {'screenshot','screenshots'}
    fileType = 'screenshots';
    fileDir = 'screenShots';
    if isempty(largefiles),largefiles = false;end
    if ~iscell(filePath) && isdir(filePath)
      filePath = getFilenames(filePath,'*.png');
    end
   % Add freesurfer directory to DB
   case 'freesurfer'
    fileType = 'freesurfer';
    fileDir = 'anatomy';
    if isempty(largefiles),largefiles = true;end
   % Add canonical to DB
   case {'3d','canonical','anatomy','anat'}
    fileDir = 'anatomy';
    if isempty(largefiles),largefiles = false;end
   % Add mlrBase Anatomies to DB
   case {'baseanat','mlrbaseanat','mlrbaseanatomy','mlrbaseanatomies'}
    fileType = 'mlrbaseanat';
    fileDir = 'mlrBaseAnatomies';
    if isempty(largefiles),largefiles = false;end
    if isview(filePath)
      % if a view then have user select which bases to get
      filePath = getBaseFilenames(filePath);
    elseif isdir(filePath)
      filePath = getFilenames(filePath,'*.mat');
    end
   % Add surfaces to DB
   case {'surfaces','surface'}
    fileDir = 'surfaces';
    if isempty(largefiles),largefiles = false;end
    if isdir(filePath),filePath = getFilenames(filePath,{'*.off','*.vff','*.hdr','*.img','*.nii'});end
   % Add localizer sessions to DB
   case {'localizer','localizers'}
    fileType = 'localizers';
    fileDir = 'localizers';
    if isempty(largefiles),largefiles = true;end
  end    
end
if isempty(filePath),return,end
if isempty(fileDir)
  disp(sprintf('(mlrAnatDBPut) Must specify a vaild fileType or fileDir'));
  return
end

% get local repo
if ~largefiles
  [localRepo] = mlrAnatDBGetRepo(subjectID);
  if isempty(localRepo),return,end
else
  [localRepo localRepoLargeFiles] = mlrAnatDBGetRepo(subjectID);
  if isempty(localRepo)||isempty(localRepoLargeFiles),return,end
end

% get location to copy to 
if largefiles
  destDir = fullfile(localRepoLargeFiles,fileDir);
else
  destDir = fullfile(localRepo,fileDir);
end

% make directory if necessary
if ~isdir(destDir), mkdir(destDir); end

% make sure filePath is a cell array (i.e. one cell for each file)
filePath = cellArray(filePath);

% remember current directory
curpwd = pwd;

% hard link or copy each file using force
for iFile = 1:length(filePath)
  % check if we are trying to copy onto ourself, then don't bother
  if ~strcmp(fileparts(filePath{iFile}),destDir)
    if useHardLinks
      disp(sprintf('(mlrAnatDBPut) Hard linking %s to %s',filePath{iFile},destDir));
      if isdir(filePath{iFile})
	% rsync if a directory
	[status,result] = system(sprintf('rsync -a --link-dest=%s %s/ %s',filePath{iFile},filePath{iFile},fullfile(destDir,getLastDir(filePath{iFile}))));
      else
	% link if a file
	[status,result] = system(sprintf('ln -f %s %s',filePath{iFile},fullfile(destDir,getLastDir(filePath{iFile}))));
	% if failed, then..
	if status
	  % try copy
	  success = copyfile(filePath{iFile},fullfile(destDir,getLastDir(filePath{iFile})),'f');
	  status = ~success;
	end
      end
      % check success
      if status
	disp(sprintf('(mlrAnatDBPut) Could not link/copy file to repo: %s',result));
	return
      end
    else
      disp(sprintf('(mlrAnatDBPut) Copying %s to %s',filePath{iFile},destDir));
      success = copyfile(filePath{iFile},fullfile(destDir,getLastDir(filePath{iFile})),'f');
      % check success
      if ~success,return,end
    end
  end
  % get where it was copied to
  toPath{iFile} = fullfile(fileDir,getLastDir(filePath{iFile}));
end

% check base anatomies for flat off files that should
% get put in the repo as well.
if strcmp(lower(fileType),'mlrbaseanat')
  % look for base matfiles
  for iFile = 1:length(filePath)
    % now look for any with .mat which contain the base structures
    if strcmp(getext(filePath{iFile}),'mat');
      % look for any flats with a flatFileName  
      b = load(filePath{iFile});
      if isfield(b,'base') && isfield(b.base,'coordMap') && isfield(b.base.coordMap,'flatFileName');
	% get the flat file
	flatFile = fullfile(b.base.coordMap.path,b.base.coordMap.flatFileName);
	% if it exists then put it into the list, and make a new name for it
	if isfile(flatFile)
	  % get the original name and check that it is a file
	  flatFrom = fullfile(b.base.coordMap.path,b.base.coordMap.flatFileName);
	  if isfile(flatFrom)
	    % make the new name
	    flatFileSaveName = setext(b.base.name,'off');
	    % and also change the name in this base
	    b.base.coordMap.flatFileName = flatFileSaveName;
	    % and save the base
	    cd(localRepo);
	    save(toPath{iFile},'-struct','b','base');
	    % now link the flat off file into surfaces
	    flatTo = fullfile('surfaces',flatFileSaveName);
	    [status,result] = system(sprintf('ln -f %s %s',flatFrom,flatTo));
            % if failed, then..
	    if status
	      % try copy
	      success = copyfile(flatFrom,flatTo,'f');
	      status = ~success;
	    end
	    cd(curpwd);
	    if status
	      disp(sprintf('(mlrAnatDBPut) Could not save off file for flat %s',flatFrom));
	    else
	      % put it in the list of stuff to commit
	      toPath{end+1} = flatTo;
	    end
	  else
	    disp(sprintf('(mlrAnatDBPut) Could not find off file for flat: %s',flatFrom));
	  end
	end
      end
    end
  end
end

% make links to canonical for surfaces
if strcmp(lower(fileType),'surfaces')
  % get list of nifti files in directory - these should be canonical
  canonicals = getFilenames(destDir,{'*.hdr','*.img','*.nii'});
  cd(localRepo);
  % make links to canonical
  for iFile = 1:length(canonicals)
    linkFrom = fullfile('surfaces',getLastDir(canonicals{iFile}));
    toPath{end+1} = setext(subjectID,getext(canonicals{iFile}));
    mysystem(sprintf('ln -sfh %s %s',linkFrom,toPath{end}));
  end
  % make a link to freesurfer dir(this will not be tracked by 
  % hg because thereis an .hgignore file already there to hide it)
  cd('surfaces');
  linkFrom = fullfile('..','anatomy',freesurfer);
  mysystem(sprintf('ln -sfh %s freesurfer',linkFrom));
  % and also make a text file that contains the freesurfer link
  if isfile('.freesurfer')
    mysystem(sprintf('rm -f .freesurfer'));
  end
  % store what the file links to in the .freesurfer file
  % this is so the mlrAnatDBGetRepo can remake the link as needed
  % since hg will not store a link to outside the repo
  mysystem(sprintf('echo %s > .freesurfer',fullfile('anatomy',freesurfer)));
  toPath{end+1} = fullfile('surfaces','.freesurfer');
end

% if this is a freesurfer directory, then go look for surfRelax
% so that we can remove it and then add it back later to the surfaces dir
if strcmp(lower(fileType),'freesurfer')
  surfRelaxDir = fullfile(localRepoLargeFiles,toPath{1},'surfRelax');
  if isdir(surfRelaxDir)
    % remove it (we are going to add it back later to surfaces)
    rmdir(surfRelaxDir,'s');
    surfRelaxDir = fullfile(filePath{1},'surfRelax');
  else
    surfRelaxDir = '';
  end
end

if any(strcmp(lower(fileType),{'freesurfer','localizers'}))
  % update the branch numbers for both the small and large files
  branchNum = mlrAnatDBGetBranchNum(localRepo);
  mlrAnatDBSetBranchNum(localRepo,branchNum+1);
  mlrAnatDBSetBranchNum(localRepoLargeFiles,branchNum+1);
elseif any(strcmp(lower(fileType),{'rois','mlrbaseanat'}))
  % update the branch numbers for small files
  branchNum = mlrAnatDBGetBranchNum(localRepo);
  mlrAnatDBSetBranchNum(localRepo,branchNum+1);
end

% now do the commit part
if largefiles
  % add and commit to large files repo
  cd(localRepoLargeFiles);
  comments = addCommit(toPath,largefiles,comments,verbose);
else
  % add and commit to the repo
  cd(localRepo);
  addCommit(toPath,false,comments,verbose);
end

cd(curpwd);

% ask if the user wants to push
if ~verbose || isequal(questdlg(sprintf('Do you want to push to the central repo: %s? This can take several minutes depending on your connection. If you choose no now, you will need to push using mlrAnatDBPush later.',mrGetPref('mlrAnatDBCentralRepo')),'Do push?','Yes','No','Yes'),'Yes')
  mlrAnatDBPush(subjectID);
end

% for freesurfer, go and add the surfRelax directory as surfaces
if strcmp(lower(fileType),'freesurfer') && ~isempty(surfRelaxDir)
  tf = mlrAnatDBPut(subjectID,surfRelaxDir,'surfaces','comments',comments,'freesurfer',filePath{1});
  if ~tf,return,end
end

% success.
tf = true;

% the code after this is just to put stuff on to the wiki page...

% wiki info
mlrAnatDBWikiHostname = mrGetPref('mlrAnatDBWikiHostname');
mlrAnatDBWikiDirname = mrGetPref('mlrAnatDBWikiDirname');
if isempty(mlrAnatDBWikiHostname) || isempty(mlrAnatDBWikiDirname),return,end

% get the index page
indexPageLoc = sprintf('%s/pages/mlranatdb/index.txt',mlrAnatDBWikiDirname);
[status,indexPage] = system(sprintf('ssh %s ''cat %s''',mlrAnatDBWikiHostname,indexPageLoc));
if status ~= 0
  disp(sprintf('(mlrAnatDBPut) Could not update wiki. Unable to get index file from ssh://%s/%s',mlrAnatDBWikiHostname,indexPageLoc));
  return
end
% parse entries
indexPageScan = textscan(indexPage,'%s%s%s%s%s%s%s%s','Headerlines',1,'Delimiter','|');
% go look for this subjectID
matchLine = [];
for iSubject = 1:length(indexPageScan{1})
  if ~isempty(findstr(subjectID,indexPageScan{2}{iSubject}))
    matchLine = iSubject;
  end
end
% if it is empty then need to add a line
if isempty(matchLine)
  indexPageScan{1}{end+1} = '';
  indexPageScan{2}{end+1} = sprintf('[[%s',subjectID);
  indexPageScan{3}{end+1} = sprintf('%s]]',subjectID);
  indexPageScan{4}{end+1} = '';
  indexPageScan{5}{end+1} = '';
  indexPageScan{6}{end+1} = '';
  indexPageScan{7}{end+1} = '';
  indexPageScan{8}{end+1} = '';
  matchLine = length(indexPageScan{1});
end
% if this is a segmentation or an roi then update that
if strcmp(fileType,'freesurfer')
  indexPageScan{4}{matchLine} = mlrAnatDBGetUsername;
  indexPageScan{5}{matchLine} = datestr(now);
elseif strcmp(fileType,'rois')
  indexPageScan{6}{matchLine} = mlrAnatDBGetUsername;
  indexPageScan{7}{matchLine} = datestr(now);
end
% set the last update field
indexPageScan{8}{matchLine} = datestr(now);

% write temporary file
indexTempFilename = fullfile(mrGetPref('mlrAnatDBLocalRepo'),'index.txt');
f = fopen(indexTempFilename,'w');
if f == -1
  disp(sprintf('(mlrAnatDBPut) Could not open temporary index file %s',indexTempFilename));
  return
end
fprintf(f,'^SubjectID^Segmentation by^Segmentation date^ROIs by^ROIs date^Last updated^\n');
for iLine = 1:length(indexPageScan{1})
  fprintf(f,'|%s |%s |%s |%s |%s |%s |%s |\n',indexPageScan{2}{iLine},indexPageScan{3}{iLine},indexPageScan{4}{iLine},indexPageScan{5}{iLine},indexPageScan{6}{iLine},indexPageScan{7}{iLine},indexPageScan{8}{iLine});
end
fclose(f);
% write it back
[status,result] = system(sprintf('scp %s %s:%s',indexTempFilename,mlrAnatDBWikiHostname,fullfile(mlrAnatDBWikiDirname,'pages','mlranatdb','index.txt')));
if (status~=0)
  disp(sprintf('(mlrAnatDBPut) Could not write wiki index file: %s',result));
  return
end
% remove temporary file
mysystem(sprintf('rm -f %s',indexTempFilename));

% get branchNumber
branchNum = mlrAnatDBGetBranchNum(localRepo);

% now that we have written the index page to the wiki, read the individual subject page
subjectPageLoc = sprintf('%s/pages/mlranatdb/%s.txt',mlrAnatDBWikiDirname,subjectID);
[status,subjectPage] = system(sprintf('ssh %s ''cat %s''',mlrAnatDBWikiHostname,subjectPageLoc));
% parse entries
if ~isequal(status,0)
  % start a new page
  for i = 1:6
    subjectPageScan{i} = {};
  end
  restOfPage = '';
else
  % find out where log ends
  endLogPos = findstr('end log',subjectPage);
  if isempty(endLogPos)
    endLogPos = length(subjectPage);
  else
    endLogPos = endLogPos-1;
  end
  restOfPage = subjectPage(endLogPos+9:end);
  % convert [[ and ]] to quotes
  replaceLoc = findstr('[[',subjectPage(1:endLogPos));
  for iReplaceLoc = 1:length(replaceLoc)
    subjectPage(replaceLoc(iReplaceLoc))='"';
  end
  replaceLoc = findstr(']]',subjectPage(1:endLogPos));
  for iReplaceLoc = 1:length(replaceLoc)
    subjectPage(replaceLoc(iReplaceLoc)+1)='"';
  end
  % read the page in
  subjectPageScan = textscan(subjectPage(1:endLogPos),'%q%q%q%q%q%d','Headerlines',2,'Delimiter','|');
end
% write temporary file
subjectTempFilename = fullfile(mrGetPref('mlrAnatDBLocalRepo'),sprintf('%s.txt',subjectID));
f = fopen(subjectTempFilename,'w');
if f == -1
  disp(sprintf('(mlrAnatDBPut) Could not open temporary subject file %s',subjectTempFilename));
  return
end
% write header
fprintf(f,sprintf('====== %s log ======\n',subjectID));
% write header of table
fprintf(f,'^Updated by^Update type^Update time^Comments^branchNum^\n');
% write new entries
dateStrNow = datestr(now);
if strcmp(fileType,'rois')
  fprintf(f,'|%s |[[#%s_%s|%s]] |%s |%s |%i |\n',mlrAnatDBGetUsername,'ROIs',fixBadChars(dateStrNow,{' ','_'}),fileType,dateStrNow,comments,branchNum);
elseif strcmp(fileType,'screenshots')
  fprintf(f,'|%s |[[#%s_%s|%s]] |%s |%s |%i |\n',mlrAnatDBGetUsername,'ScreenShots',fixBadChars(dateStrNow,{' ','_'}),fileType,dateStrNow,comments,branchNum);
else
  fprintf(f,'|%s |%s |%s |%s |%i |\n',mlrAnatDBGetUsername,fileType,datestr(now),comments,branchNum);
end
% write old entries
for iLine = 1:length(subjectPageScan{1})
  if subjectPageScan{3}{iLine}(1) == '['
    % print with a link: [[ ]]
    fprintf(f,'|%s |[%s] |%s |%s |%i |\n',subjectPageScan{2}{iLine},strtrim(subjectPageScan{3}{iLine}),subjectPageScan{4}{iLine},subjectPageScan{5}{iLine},subjectPageScan{6}(iLine));
  else
    % print a normal field
    fprintf(f,'|%s |%s |%s |%s |%i |\n',subjectPageScan{2}{iLine},subjectPageScan{3}{iLine},subjectPageScan{4}{iLine},subjectPageScan{5}{iLine},subjectPageScan{6}(iLine));
  end
end
fprintf(f,'end log\n');
% write out rois into log
if strcmp(fileType,'rois')
  fprintf(f,'\n====== ROIs %s ======\n',dateStrNow);
  % cycle through ROIs and put information into table
  if isview(v)
    fprintf(f,'Created by: %s\n\n',mlrAnatDBGetUsername);
    fprintf(f,'^roiName ^date ^createdFromSession ^createdOnBase ^voxelSize ^numVoxels ^volume (mmxmmxmm) ^\n');
    for iROI = 1:length(filePath)
      % get roi name
      [roiPath roiName] = fileparts(filePath{iROI});
      roiName = stripext(roiName);
      % get roi number
      roiNum = viewGet(v,'roiNum',roiName);
      if ~isempty(roiNum)
	roi = viewGet(v,'roi',roiNum);
	fprintf(f,'|%s |%s |%s |%s |%s |%i |%0.1f |\n',roi.name,roi.date,roi.createdFromSession,roi.createdOnBase,mlrnum2str(roi.voxelSize),size(roi.coords,2),size(roi.coords,2)*prod(roi.voxelSize));
      end
    end
  end
% write out screenshots
elseif strcmp(fileType,'screenshots')
  fprintf(f,'\n====== ScreenShots %s ======\n',dateStrNow);
  fprintf(f,'Created by: %s\n\n',mlrAnatDBGetUsername);
  for iFile = 1:length(filePath)
    filename = getLastDir(filePath{iFile});
    filename = sprintf('%s_%s_%s.%s',subjectID,stripext(filename),fixBadChars(dateStrNow),getext(filename));
    fprintf(f,'<html><div></html>\n');
    % get image size if we have imread
    if exist('imread') > 0
      im = imread(filePath{iFile});
      imSize = size(im);
      imageRatio = imSize(1)/imSize(2);
      fprintf(f,'{{%s?600x%i|}}\n',filename,round(600*imageRatio));
    else
      % just write it without any scale factor
      fprintf(f,'{{%s|}}',filename);
    end
    fprintf(f,'<html></div></html>\n');
    % write the png file
    [status,result] = system(sprintf('scp %s %s:%s',filePath{iFile},mlrAnatDBWikiHostname,fullfile(mlrAnatDBWikiDirname,'media','mlranatdb',filename)));
    if (status~=0)
      disp(sprintf('(mlrAnatDBPut) !!! Could not write wiki image file: %s !!!',filename));
    end
  end
end
% print rest of file
fprintf(f,'%s',restOfPage);
fclose(f);
% write it back
[status,result] = system(sprintf('scp %s %s:%s',subjectTempFilename,mlrAnatDBWikiHostname,fullfile(mlrAnatDBWikiDirname,'pages','mlranatdb',sprintf('%s.txt',subjectID))));
if (status~=0)
  disp(sprintf('(mlrAnatDBPut) Could not write wiki subject file: %s',result));
  return
end
% remove temporary file
mysystem(sprintf('rm -f %s',subjectTempFilename));
%%%%%%%%%%%%%%%%%%%
%    addCommit    %
%%%%%%%%%%%%%%%%%%%
function comments = addCommit(toPath,largefiles,comments,verbose)

% commit to repo
if largefiles
  disppercent(-inf,sprintf('(mlrAnatDBAddCommitPush) Committing files to repo %s. This may take a minute or two...',pwd));
else
  disp(sprintf('(mlrAnatDBAddCommitPush) Committing files to repo %s',pwd));
end

% add file to repo
for iFile = 1:length(toPath)
  if largefiles
    [status,result] = mysystem(sprintf('hg add %s --large',toPath{iFile}));
  else
    [status,result] = mysystem(sprintf('hg add %s',toPath{iFile}));
  end

  if status ~= 0
    mrWarnDlg(sprintf('(mlrAnatDBPlugin) Could not add files to local repo. Have you setup your config file for hg?'));
    return
  end
end

% commit
[status,result] = mysystem(sprintf('hg commit -m ''%s''',comments));
if largefiles,disppercent(inf);,end

%%%%%%%%%%%%%%%%%%%%%%
%    getFilenames    %
%%%%%%%%%%%%%%%%%%%%%%
function filePathOut = getFilenames(filePath,matchStr)

%getfilenames that match string
matchStr = cellArray(matchStr);

% go through and find matching files
filePathOut = {};
for iMatch = 1:length(matchStr)
  dirPath = dir(fullfile(filePath,matchStr{iMatch}));
  for iFile = 1:length(dirPath)
    filePathOut{end+1} = fullfile(filePath,dirPath(iFile).name);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   getBaseFilenames   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%
function filePath = getBaseFilenames(v)

% get list of ROIs to save
baseList = selectInList(v,'bases','Select MLR Base Anatomies to save');

% list of extensions that might be generated
extList = {'hdr','img','nii','mat'};

% default to no files
filePath = {};

% save each anat locally, getting the filenames as the ones to put into the repo
for iBase = baseList
  thisAnat = saveAnat(v,iBase,false,false);
  [thisAnatPath thisAnatFilename] = fileparts(thisAnat);
  % get all the files associated with thisAnat
  for iExt = 1:length(extList)
    filename = fullfile(thisAnatPath,setext(thisAnatFilename,extList{iExt}));
    if isfile(filename)
      filePath{end+1} = filename;
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%   getROIFilenames   %%
%%%%%%%%%%%%%%%%%%%%%%%%%
function [filePath niftiFilePath screenShotsFilePath] = getROIFilenames(v,subjectID)

% default to no files
filePath = {};
niftiFilePath = {};
screenShotsFilePath = {};

% get list of ROIs to save
roiList = selectInList(v,'rois','Select ROI(s) to save');
if isempty(roiList),return,end

% get info about the repo
localRepo = mlrAnatDBGetRepo(subjectID);
branchNum = mlrAnatDBGetBranchNum(localRepo);
localRepoNiftiROI = fullfile(localRepo,'niftiROIs');
localRepoScreenShots = fullfile(localRepo,'screenShots');

% session name
sessionName = getLastDir(viewGet(v,'homeDir'));
baseName = viewGet(v,'baseName');

% get username
userName = mlrAnatDBGetUsername;

% current settings
curROI = viewGet(v,'curROI');
curBase = viewGet(v,'curBase');
curROIGroup = viewGet(v,'roiGroup');
showROI = viewGet(v,'showROIs');

% start some lists
displayOnBase = {};

% go through list of chosen ROIs and make sure that
% see if they all have valid createdOnBase and displayOnBase
paramsInfo = {};
baseNames = putOnTopOfList('None',viewGet(v,'baseNames'));
missingCreatedOn = {};
missingDisplayOn = {};
for iROI = roiList
  % get the roi
  roi = viewGet(v,'roi',iROI);
  % check createdOnBase
  if isempty(viewGet(v,'baseNum',roi.createdOnBase)) 
    % mark as missing
    missingCreatedOn{end+1} = roi.name;
    paramsInfo{end+1} = {sprintf('%s_createdOnBase',fixBadChars(roi.name)),baseNames,sprintf('Choose a base to use for createdOn for %s. If you set to none then will not be able to export a Nifit ROI but otherwise will be ok',roi.name)};
  end
  % check displayOnBase
  if isempty(viewGet(v,'baseNum',roi.displayOnBase))
    % mark as missing
    missingDisplayOn{end+1} = roi.name;
    paramsInfo{end+1} = {sprintf('%s_displayOnBase',fixBadChars(roi.name)),baseNames,sprintf('Choose a base to use for displayOn for %s. If you set to none then will not be able to make a snapshot for uploading to the wiki but otherwise will be ok',roi.name)};
  end
end
% if there are bases that are missing then put up the dialog to choose
if ~isempty(paramsInfo)
  mrWarnDlg('Some ROIs for export to mlrAnatDB are missing createdOnBase or displayOnBase - if createdOnBase is missing then the ROIs will not be exported as Nifti. If displayOnBase is missing will not be able to export a snapshot to the wiki. Otherwise export will be fine. If you want to fix the bases you can set them in the dialog');
  params = mrParamsDialog(paramsInfo,'Fix ROI created/displayOnBase');
  if isempty(params),return,end
  % now fix them
  for iROI = 1:length(missingCreatedOn)
    val = params.(sprintf('%s_createdOnBase',missingCreatedOn{iROI}));
    % fix createdOnBase if not set to none
    if ~isequal(val,'None')
      v = viewSet(v,'roiCreatedOnBase',val,viewGet(v,'roiNum',missingCreatedOn{iROI}));
    end
  end
  % now fix displayOnBase
  for iROI = 1:length(missingDisplayOn)
    val = params.(sprintf('%s_displayOnBase',missingDisplayOn{iROI}));
    % fix displayOnBase if not set to None
    if ~isequal(val,'None')
      v = viewSet(v,'roiDisplayOnBase',val,viewGet(v,'roiNum',missingDisplayOn{iROI}));
    end
  end
end

% go through the list of chosen ROIS
for iRoi = roiList
  % get the roi
  roi = viewGet(v,'roi',iRoi);
  % set user who is saving this
  v = viewSet(v,'roiCreatedBy',userName,iRoi);
  % set subjectID
  v = viewSet(v,'roiSubjectID',subjectID,iRoi);
  % set session if not already
  if isempty(roi.createdFromSession)
    % set that the ROI is created from this session
    v = viewSet(v,'roiCreatedFromSession',sessionName,iRoi);
  end
  % if createdOnBase is not set
  if isempty(roi.createdOnBase)
    % set to current base name
    v = viewSet(v,'roiCreatedOnBase',baseName,iRoi);
  end
  % if displayOnBase is not set
  if isempty(roi.displayOnBase)
    % set to current base name
    v = viewSet(v,'roiDisplayOnBase',baseName,iRoi);
  end
  % set the branch num of repo
  if isempty(roi.branchNum)
    % set that the ROI is created from this session
    v = viewSet(v,'roiBranchNum',branchNum+1,iRoi);
  end
  % save it  
  filePath{end+1} = saveROI(v,iRoi,false);
  % set the .mat extension
  filePath{end} = fullfile(fileparts(filePath{end}),setext(getLastDir(filePath{end}),'mat'));
  
  % now convert into a nifti file and save that into the repo
  roiName = viewGet(v,'roiName',iRoi);
  createdOnBase = viewGet(v,'roiCreatedOnBase',iRoi);
  if ~isempty(createdOnBase)
    baseNum = first(viewGet(v,'baseNum',createdOnBase));
    % if the base exists
    if isempty(baseNum)
      % does not exist
      disp(sprintf('(mlrAnatDBPut) Could not find base %s which ROI %s was created on. Not exporting to Nifti',createdOnBase,viewGet(v,'roiName',iRoi)));
    else
      % does exist, so try to export as nifti. Get the header
      hdr = viewGet(v,'baseHdr',baseNum);
      % make it have the dimensions of the anatomy it was created from
      coordMap = viewGet(v,'baseCoordMap',baseNum);
      if ~isempty(coordMap)
	hdr.dim(1:4) = [3 coordMap.dims];
      end
      % set the current roi (to export)
      v = viewSet(v,'curROI',iRoi);
      niftiFilePath{end+1} = fullfile(localRepoNiftiROI,setext(roiName,'nii'));
      % export the ROI to nifti
      mlrExportROI(v,niftiFilePath{end},'hdr',hdr);
    end
  end

  % remember which displayOnBases have been seen
  displayOnBase{end+1} = viewGet(v,'roiDisplayOnBase',iRoi);
  if isempty(displayOnBase{end})
    displayOnBase = displayOnBase{1:end-1};
  end
end

% now try to big up base images with the rois drawn on them for saving
uniqueDisplayOnBase = unique(displayOnBase);
% for each base, see if it is loaded
for iBase = 1:length(uniqueDisplayOnBase)
  % get base num and make sure that it exists
  baseNum = viewGet(v,'baseNum',uniqueDisplayOnBase{iBase});
  if isempty(baseNum)
    disp(sprintf('(mlrAnatDBPut) Base %s is not loaded, so cannot make snapshot on base',uniqueDisplayOnBase{iBase}));
    continue;
  end
  % ok. Now load the base into the view
  v = viewSet(v,'curBase',baseNum);
  % set the roi group to all rois that have this base as their displayOnBase
  v = viewSet(v,'roiGroup',viewGet(v,'roiName',roiList(strcmp(uniqueDisplayOnBase{iBase},displayOnBase))));
  v = viewSet(v,'showROIs','group perimeter');
  % make screen shot
  mrPrint(v,'useDefault=1','roiSmooth=0');
  % ask whether user wants to save screen shot
  paramsInfo = {{'saveScreenShot',uniqueDisplayOnBase{iBase},'Save this image as a screen shot to DB'}};
  params = mrParamsDialog(paramsInfo,'Do you want to save this screen shot?');
  if ~isempty(params)
    
    screenShotsFilePath{end+1} = fullfile(localRepoScreenShots,setext(params.saveScreenShot,'png'));
    print(selectGraphWin(true),'-dpng',screenShotsFilePath{end});
  end
  closeGraphWin;
end

% set back
v = viewSet(v,'curROI',curROI);
v = viewSet(v,'curBase',curBase);
v = viewSet(v,'roiGroup',curROIGroup);
v = viewSet(v,'showROIs',showROI);
refreshMLRDisplay(v);


%%%%%%%%%%%%%%%%%%
%    mysystem    %
%%%%%%%%%%%%%%%%%%
function [status,result] = mysystem(command)

disp(sprintf('(mlrAnatDBPut): %s',command));
[status,result] = system(command,'-echo');

