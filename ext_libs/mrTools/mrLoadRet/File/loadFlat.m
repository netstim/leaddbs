% loadFlat.m
%
%        $Id$
%      usage: loadFlat(v,flatFileName)
%         by: justin gardner
%       date: 07/20/07
%    purpose: 
%
function v = loadFlat(v,flatFileName)

% check arguments
if ~any(nargin == [0 1 2])
  help loadFlat
  return
end

% get mrGlobals and view
mrGlobals;

% Open dialog box to have user choose the file
if ieNotDefined('flatFileName')
  if ieNotDefined('flatFilePath')
    startPathStr = mrGetPref('volumeDirectory');
  else
    startPathStr = flatFilePath;
  end
  filterspec = {'*.mat','Matlab flat file'};
  title = 'Choose flat file';
  pathStr = mlrGetPathStrDialog(startPathStr,title,filterspec,'on');
else
  pathStr = flatFileName;
end

% Aborted
if ieNotDefined('pathStr')
  return
end

% make into a cell array
pathStr = cellArray(pathStr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% now load the base anatomy file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filterspec = {'*.hdr','Nifti file header'};
title = 'Choose volume anatomy file header for this flat image';
anatPathStr = mlrGetPathStrDialog(fileparts(pathStr{1}),title,filterspec);

% Strip extension to make sure it is .mat
anatPathStr = [stripext(anatPathStr),'.hdr'];

% File does not exist
if ~exist(anatPathStr,'file')
  mrWarnDlg(['File ',anatPathStr,' not found']);
  return
end

% load the anatomy file header
hdr = mlrImageReadNiftiHeader(anatPathStr);

% Extract permutation matrix to keep track of slice orientation.
% This logic which is admittedly arcane is duplicated in mrAlignGUI. If you
% make changes here, please update that function as well.
permutationMatrix = getPermutationMatrix(hdr);

% get parameters for renderinf flat maps
paramsInfo = {{'resolution',1,'minmax=[1 inf]','incdec=[-1 1]','Resolution of flat map. Make number larger for finer resolution. The finer the resolution the longer it will take to render the flat map'},...
	      {'threshold',1,'type=checkbox','Threshold the flat map to two levels of gray (one for low and one for high curvature)'}};
params = mrParamsDialog(paramsInfo,'Flat map rendering options');
if isempty(params),return,end
drawnow
% now go through all flat file names, load and convert
% into images, then into base anatomies
for fileNum = 1:length(pathStr)
  % Strip extension to make sure it is .mat
  filename = [stripext(pathStr{fileNum}),'.mat'];

  % File does not exist
  if ~exist(filename,'file')
    mrWarnDlg(['File ',filename,' not found']);
    continue
  end

  % load the flat file
  flat = load(filename);

  % check its fields
  if ~isfield(flat,'curvature') || ~isfield(flat,'gLocs2d') || ~isfield(flat,'gLocs3d')
    mrWarnDlg(sprintf('(loadFlat) %s is not a flat file generated using SurfRelax',pathStr{fileNum}));
    continue
  end


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % generate the flat image
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % first get coordinates
  xmin = min(flat.gLocs2d(:,1));
  xmax = max(flat.gLocs2d(:,1));
  ymin = min(flat.gLocs2d(:,2));
  ymax = max(flat.gLocs2d(:,2));
  x = xmin:(1/params.resolution):xmax;
  y = ymin:(1/params.resolution):ymax;

  % now we need to interp the curvature on to a regular
  % grid. we do this by finding the x,y point that is closest
  % to the one on the grid.
  h = mrWaitBar(0,sprintf('Creating flat image for %s',filename));
  for i = 1:length(x)
    mrWaitBar(i/length(x),h);
    for j = 1:length(y)
      % find nearest point in curvature
      dist = (flat.gLocs2d(:,1)-x(i)).^2+(flat.gLocs2d(:,2)-y(j)).^2;
      flat.pos(i,j) = first(find(min(dist)==dist));
      % points that are greater than a distance of 5 away are
      % probably not in the flat patch so mask them out
      if (min(dist) < 5)
	flat.mask(i,j) = 1;
	flat.baseCoords(i,j,:) = flat.gLocs3d(flat.pos(i,j),:);
	flat.map(i,j) = flat.curvature(flat.pos(i,j));
      else
	flat.mask(i,j) = 0;
	flat.baseCoords(i,j,:) = [0 0 0];
	flat.map(i,j) = 0;
      end
    end
  end
  mrCloseDlg(h);

  % now blur image
  flat.blurMap(:,:) = blur(flat.map(:,:));
  % threshold image
  flat.median = median(flat.blurMap(:));
  flat.thresholdMap(:,:) = (flat.blurMap(:,:)>median(flat.blurMap(:)))*0.5+0.25;
  flat.thresholdMap = blur(flat.thresholdMap);
  % mask out points not on map
  flat.thresholdMap(~flat.mask(:)) = 0;

  % now generate a base structure
  clear base;
  base.hdr = hdr;
  base.name = getLastDir(pathStr{fileNum});
  base.permutationMatrix = permutationMatrix;

  % load all the flat maps into the base. We
  % need to make all the flat images have
  % the same width and height.
  if params.threshold
    base.data(:,:,1) = flat.thresholdMap;
    base.range = [0 1];
    base.clip = [0 1];
  else
    base.data(:,:,1) = flat.map;
  end    
  base.coordMap.coords(:,:,1,1) = flat.baseCoords(:,:,3);
  base.coordMap.coords(:,:,1,2) = hdr.dim(3)-flat.baseCoords(:,:,1)+1;
  base.coordMap.coords(:,:,1,3) = hdr.dim(4)-flat.baseCoords(:,:,2)+1;
  base.coordMap.dims = hdr.dim([2 3 4])';

  v = viewSet(v,'newBase',base);
end



