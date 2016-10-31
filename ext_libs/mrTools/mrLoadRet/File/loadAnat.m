function [view anatFilePath] = loadAnat(view,anatFileName,anatFilePath)
%
%        $Id$
% view = loadAnat(view,[anatFileName],[anatFilePath])
%
% Loads an anatomy array and sets view.baseVolumes to include.
%
% anatFileName: Anatomy file. If not specified or empty, prompts user to
% select the file. Anatomy file must be nifti format. anatFileName must be
% in the <homedir>/Anatomy subdirectory. Use symbolic links to point to a
% shared directory for volume anatomy files.
%
% anatFilePath is the path of where to open up the dialog. This
% will bey returned so that the GUI can open up in the same
% place each time.
%
% djh, 1/9/98
% 5/2005, djh, update to mrLoadRet-4.0

%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the anatomy file %
%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0,help loadAnat,return,end

if ieNotDefined('anatFilePath')
    anatFilePath = '';
end

% Open dialog box to have user choose the file
if ieNotDefined('anatFileName')
    if ieNotDefined('anatFilePath')
        startPathStr = viewGet(view,'anatomyDir');
    else
        startPathStr = anatFilePath;
    end
    filterspec = {'*.img;*.nii;*.nii.gz','Nifti/Analyze files'};
    title = 'Choose anatomy file';
    pathStr = mlrGetPathStrDialog(startPathStr,title,filterspec,'on');
elseif ieNotDefined('anatFilePath')
    pathStr = fullfile(viewGet(view,'anatomyDir'),anatFileName);
else
    pathStr = fullfile(anatFilePath,anatFileName);
end

% Aborted
if ieNotDefined('pathStr')
  return
end

% make sure we have a cell array of paths to load
pathStr = cellArray(pathStr);

% cycle through all paths and load each one
for pathNum = 1:length(pathStr)
  % Check whether extension of the file is img or nii
  % matlab function "fileparts" takes the last .foo as extension!
  [path,name,extension] = fileparts(pathStr{pathNum});
  % set default extension, if extension not specified
  if isempty(extension)
    pathStr{pathNum} = setext(pathStr{pathNum},mrGetPref('niftiFileExtension'),0);
 % extension is not .nii or .img
  elseif ~any(strcmp(extension,{'.nii', '.img', '.hdr','.gz'}))
    mrWarnDlg(['(loadAnat) File type ',extension,' is not a valid anatomy file format']);
    return
  end

  % File does not exist
  if ~exist(pathStr{pathNum},'file')
    mrWarnDlg(['(loadAnat) File ',pathStr{pathNum},' not found']);
    return
  end

  anatFilePath = path;

  % Load nifti file
  h = mrMsgBox(['(loadAnat) Loading volume: ',pathStr{pathNum},'. Please wait']);
  if ishandle(h), close(h), end
  hdr = mlrImageHeaderLoad(pathStr{pathNum});
  if isempty(hdr),return,end
  hdr = mlrImageGetNiftiHeader(hdr);

  % get volume dimension...
  % hdr.dim(1) should probably be the number of dimensions
  % but some files don't seem to have this set properly
  % so this check seems to work
  volumeDimension = sum((hdr.dim(2:end)~=1)&(hdr.dim(2:end)~=0));

  % Error if it dimension is greater than 4D.
  if (volumeDimension > 4)
    mrErrorDlg(['(loadAnat) Volume must be 3D or 4D. This file contains a ',num2str(volumeDimension),'D array.']);
  end

  % Handle 4D file
  if (volumeDimension == 4)
    paramsInfo = {{'frameNum',0,'incdec=[-1 1]',sprintf('minmax=[0 %i]',hdr.dim(5)),'This volume is a 4D file, to display it as an anatomy you need to choose a particular time point or take the mean over all time points. Setting this value to 0 will compute the mean, otherwise you can select a particular timepoint to display'}};
    params = mrParamsDialog(paramsInfo,'Choose which frame of 4D file. 0 for mean');
    drawnow
    if isempty(params)
      return
    end
    % if frameNum is set to 0, take the mean
    if params.frameNum == 0
      % load the whole thing and average across volumes
      [vol hdr] = mlrImageLoad(pathStr{pathNum});
      if isempty(vol),return,end
      vol = nanmean(vol,4);
    else
      % load a single volume
      [vol hdr] = mlrImageLoad(pathStr{pathNum},'volNum',params.frameNum);
      if isempty(vol),return,end
    end
  else
    % if 3D file just load
    [vol hdr] = mlrImageLoad(pathStr{pathNum});
    if isempty(vol),return,end
  end

  % add in a missing qform based on voxel dimensios. This
  % won't have correct info, but will allow the rest of the
  % code to run and can be fixed by a correct alignment
  if isempty(hdr.qform)
    mrWarnDlg(sprintf('(loadAnat) !!!! Missing qform, making one based only on voxel dimensiosn !!!!'));
    hdr.qform = diag(hdr.pixdim);
    hdr.qform(4,4) = 1;
  end

  % now check for an assoicated .mat file, this is created by
  % mrLoadRet and contains parameters for the base anatomy.
  % (it is the base structure minus the data and hdr) This is
  % essential for flat files
  base = [];
  matFilename = sprintf('%s.mat',stripext(pathStr{pathNum}));
  if isfile(matFilename)
    l = load(matFilename);
    if isfield(l,'base')
      base = l.base;
    end
  end

  % now put volume into an orientation that mrLoadRet viewer
  % will show as correct L/R and have the axial/sagittal/coronal
  % buttons work correctly. Only do this for anatomy files
  % and not surfaces or flat maps.
  if ~isfield(base,'type') || (base.type == 0)
    % get what the current orientation is and save it (for export like mlrExportROI to be able
    % to export to the original orientation
    axisLabels = mlrImageGetAxisLabels(hdr.qform);
    base.originalOrient = axisLabels.orient;
    % Convert volume into LPI orientation
    [vol hdr base.xformFromOriginal] = mlrImageOrient('LPI',vol,hdr);
    % get the nifti header from the mlrImage header
    hdr = mlrImageGetNiftiHeader(hdr);
    % permutation and sliceIndex and rotate take on default values
    % now that we load and fix the volume above
    permutationMatrix = eye(3);
    sliceIndex = [1 2 3];
    if ~isfield(base,'rotate') || isempty(base.rotate)
      base.rotate = 90;
    end
    if ~isfield(base,'curSlice') || isempty(base.curSlice)
      sliceOrientation = viewGet(view,'sliceOrientation');
      if ~isempty(sliceOrientation) && any(sliceOrientation == [1 2 3])
	% set current slice to the slice half way through the volume
	voldim = size(vol);
	base.curCoords = ceil(voldim/2);
      end
    end
  else
    % get the nifti header from the mlrImage header
    hdr = mlrImageGetNiftiHeader(hdr);
    % and set permutation and sliceIndex the old way, based on the header
    [permutationMatrix sliceIndex] = getPermutationMatrix(hdr);
  end

  % Warning if no alignment information in the header.
  if ~hdr.sform_code
    mrWarnDlg('(loadAnat) No base coordinate frame in the volume header. (i.e. sform is not set). To fix this, run mrAlign and do an alignment or Set Base Coordinate Frame. If this anatomy was taken in the same session, then overlay alignment will use the qform and should look fine.');
  end

  %%%%%%%%%%%%%%%%%%%%%%
  % Add it to the view %
  %%%%%%%%%%%%%%%%%%%%%%

  % Set required structure fields (additional fields are set to default
  % values when viewSet calls isbase).
  base.name = name;
  base.data = vol;
  base.hdr = hdr;
  base.permutationMatrix = permutationMatrix;

  % Add it to the list of base volumes and select it
  view = viewSet(view,'newBase',base);
end
