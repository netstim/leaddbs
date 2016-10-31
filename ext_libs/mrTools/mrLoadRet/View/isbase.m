function [tf base] =  isbase(base)
% function [tf base] =  isbase(base)
%
% Checks to see if it is a valid base structure. Can be called with
% either one or two output arguments:
%
% tf =  isbase(base)
% [tf base] =  isbase(base)
%
% tf is logical 1 (true) if base is a valid base structure.
% tf is logical 0 (false) if it is not.
% 
% If called with two output arguments then an attempt is made to make it
% into a valid base structure by setting optional fields to default
% values.
% 
% djh, 2007

if (nargout == 2)
  % Add optional fields and return true if the overlay with optional
  % fields is valid.
  requiredFields = {'data','hdr','name','permutationMatrix'};
  optionalFields = {'range',[min(base.data(:)) max(base.data(:))];
		    'clip',defaultClip(base.data);
		    'coordMap',[];
		    'rotate',0;
		    'surfaceRotate',0;
		    'tilt',0;
		    'curCoords',[];
		    'sliceOrientation',[],;
		    'type',[],;
		    'gamma',1,;
		    'vol2tal',[];
		    'vol2mag',[];
                    'talInfo',[];
		    'originalOrient',[];
		    'xformFromOriginal',[];
		    'alpha',1;
		    'displayOverlay',true;
		    'multiDisplay',false;
		    'multiAxis',0;
		    'overlay',[];
		    'overlayAlpha',0.5;
		    'overlays',{};
		    'curCorticalDepth',[];
		    'h',[];
		    'plane',[];
		    'fascicles',[]};
else
  % Return 0 if the overlay structure is missing any fields required or
  % optional (since w/out changing the base structure it is invalid).
  requiredFields = {'clip','coordMap','curCoords','data','hdr','name','permutationMatrix','range','rotate','surfaceRotate','sliceOrientation','type','gamma','tilt','vol2tal','vol2mag','talInfo','originalOrient','xformFromOriginal','alpha','displayOverlay','multiDisplay','multiAxis','overlay','overlayAlpha','overlays','h','plane','fascicles'};
  optionalFields = {'curCorticalDepth',0};
end

% Initialize return value
tf = true;
if ieNotDefined('base')
  tf = false;
  return
end
if ~isstruct(base)
  tf = false;
  return
end

% hack to change talPoints field to talinfo
if isfield(base,'talPoints')
  disp(sprintf('(isbase) Changing talPoints field to talInfo'));
  if ~isfield(base,'talInfo') || isempty(base.talInfo)
    base.talInfo = base.talPoints;
  end
  % remove the legacy field name
  base = rmfield(base,'talPoints');
end

% Check required fields
for f = 1:length(requiredFields)
  fieldName = requiredFields{f};
  if ~isfield(base,fieldName)
    % mrWarnDlg(['Invalid base, missing field: ',fieldName]);
    tf = false;
    return
  end
end

%convert curSlice to curCoords for backward compatibility pre-multiaxis (up to commit 7266921 on Nov 9th 2014)
if isfield(base,'type') && ~base.type && isfield(base,'curSlice') && isfield(base,'sliceOrientation') ...
  && ~isempty(base.type) && ~isempty(base.curSlice) && ~isempty(base.sliceOrientation) 
    base.curCoords = ceil(size(base.data)/2);
    %this is taken from the old version of viewGet(..,'baseSliceIndex')
    switch base.sliceOrientation %find the saved view orientation and get slice index
      case 3   % This used to be saggital
        [m,index] = max(base.permutationMatrix * [1 0 0]');
        base.sliceOrientation = 1; %saggital view is now 1
      case 2   % Coronal
        [m,index] = max(base.permutationMatrix * [0 1 0]');
      case 1   % This used to be axial
        [m,index] = max(base.permutationMatrix * [0 0 1]');
        base.sliceOrientation = 3; %axial view is now 3
    end
    base.curCoords(index)=base.curSlice; %set the current Slice in the current view orientation to whatever was saved in curSlice
end

% Optional fields and defaults
for f = 1:size(optionalFields,1)
  fieldName = optionalFields{f,1};
  default = optionalFields{f,2};
  if ~isfield(base,fieldName)  
    base.(fieldName) = default;
  end
end

% remove any fields that are not required or optional
if nargout == 2
  baseFieldNames = fieldnames(base);
  for f = 1:length(baseFieldNames)
    % check required fields
    if ~any(strcmp(baseFieldNames{f},requiredFields))
      % check optional fields, (only check first column of
      % field names, not the default values...
      match = strcmp(baseFieldNames{f},optionalFields);
      if ~any(match(:,1))
	disp(sprintf('(isbase) Removing unnecessary field %s from base',baseFieldNames{f}));
	base = rmfield(base,baseFieldNames{f});
      end
    end
  end
end

% order the fields
base = orderfields(base);

% auto set the base type if it is set to empty
if isempty(base.type)
  base.type = ~isempty(base.coordMap);
end

% validate coordMap field
if ~isempty(base.coordMap) && any(base.type==[1 2]);
  [tf base.coordMap] = validateCoordMap(base.coordMap,base.type);
  if ~tf,disp(sprintf('(isbase) Invalid baseCoordMap'));end
end

% make sure that name does not include full path
base.name = getLastDir(base.name);

%%%%%%%%%%%%%%%%%%%%%
%%   defaultClip   %%
%%%%%%%%%%%%%%%%%%%%%
function clip = defaultClip(image)
% Choose default clipping based on histogram
histThresh = length(image(~isnan(image)))/1000;
[cnt, val] = hist(image(~isnan(image)),100);
goodVals = find(cnt>histThresh);
clipMin = val(min(goodVals));
clipMax = val(max(goodVals));
%if that doesn't work, take the 5th and 95th percentiles
if clipMin==clipMax
  im = sort(image(~isnan(image)));
  clipMin = im(round(0.05*length(im)));
  clipMax = im(round(0.95*length(im)));
end
clip = [clipMin,clipMax];

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   validateCoordMap   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tf coordMap] = validateCoordMap(coordMap,baseType)

% First, we fix old coordMaps to have correctly named fields
if isfield(coordMap,'flatDir')
  disp(sprintf('(isbase) Updating format of coordMap'));
  newCoordMap.path = coordMap.flatDir;
  newCoordMap.dims = coordMap.dims;
  newCoordMap.flatFileName = coordMap.flatFileName;
  newCoordMap.innerCoordsFileName = coordMap.innerFileName;
  newCoordMap.outerCoordsFileName = coordMap.outerFileName;
  newCoordMap.curvFileName = coordMap.curvFileName;
  newCoordMap.anatFileName = coordMap.anatFileName;
  if isfield(coordMap,'coords')
    newCoordMap.coords = coordMap.coords;
  end
  newCoordMap.innerCoords = coordMap.innerCoords;
  newCoordMap.outerCoords = coordMap.outerCoords;
  coordMap = newCoordMap;
elseif isfield(coordMap,'innerFileName')
  disp(sprintf('(isbase) Updating format of coordMap'));
  if isfield(coordMap,'path')
    newCoordMap.path = coordMap.path;
  else
    newCoordMap.path = '';
  end
  newCoordMap.dims = coordMap.dims;
  newCoordMap.innerSurfaceFileName = coordMap.innerFileName;
  newCoordMap.innerCoordsFileName = coordMap.innerCoordsFileName;
  newCoordMap.outerSurfaceFileName = coordMap.outerFileName;
  newCoordMap.outerCoordsFileName = coordMap.outerCoordsFileName;
  newCoordMap.curvFileName = coordMap.curvFileName;
  newCoordMap.anatFileName = coordMap.anatFileName;
  if isfield(coordMap,'coords')
    newCoordMap.coords = coordMap.coords;
  end
  newCoordMap.innerCoords = coordMap.innerCoords;
  newCoordMap.outerCoords = coordMap.outerCoords;
  newCoordMap.innerVtcs = coordMap.innerVtcs;
  newCoordMap.outerVtcs = coordMap.outerVtcs;
  newCoordMap.tris = coordMap.tris;
  coordMap = newCoordMap;
elseif isfield(coordMap,'inner')
  disp(sprintf('(isbase) Updating format of coordMap'));
  if isfield(coordMap,'path')
    newCoordMap.path = coordMap.path;
  else
    newCoordMap.path = '';
  end
  newCoordMap.innerSurfaceFileName = coordMap.inner;
  newCoordMap.innerCoordsFileName = coordMap.inner;
  newCoordMap.outerSurfaceFileName = coordMap.outer;
  newCoordMap.outerCoordsFileName = coordMap.outer;
  newCoordMap.curvFileName = coordMap.curv;
  newCoordMap.anatFileName = coordMap.anatomy;
  if isfield(coordMap,'coords')
    newCoordMap.coords = coordMap.coords;
  end
  newCoordMap.innerCoords = coordMap.innerCoords;
  newCoordMap.outerCoords = coordMap.outerCoords;
  newCoordMap.innerVtcs = coordMap.innerVtcs;
  newCoordMap.outerVtcs = coordMap.outerVtcs;
  newCoordMap.tris = coordMap.tris;
  coordMap = newCoordMap;
end

if baseType == 1
  requiredFields = {'path','dims','flatFileName','innerCoordsFileName','outerCoordsFileName','curvFileName','anatFileName','innerCoords','outerCoords'};
  optionalFields = {'coords',[]};
elseif baseType == 2
  requiredFields = {'path','dims','innerSurfaceFileName','innerCoordsFileName','outerSurfaceFileName','outerCoordsFileName','curvFileName','anatFileName','innerCoords','outerCoords','innerVtcs','outerVtcs','tris'};
  optionalFields = {'coords',[]};
end  

% Check required fields
tf = true;
for f = 1:length(requiredFields)
  fieldName = requiredFields{f};
  if ~isfield(coordMap,fieldName)
    disp(sprintf('(isbase) Missing coordMap field: %s',fieldName));
    tf = false;
  end
end

% Optional fields and defaults
for f = 1:size(optionalFields,1)
  fieldName = optionalFields{f,1};
  default = optionalFields{f,2};
  if ~isfield(coordMap,fieldName)  
    coordMap.(fieldName) = default;
  end
end

% make sure the filenames are set correctly
if baseType == 2
  if strcmp(coordMap.innerCoordsFileName,'Same as surface')
    coordMap.innerCoordsFileName = coordMap.innerSurfaceFileName;
  end
  if strcmp(coordMap.outerCoordsFileName,'Same as surface')
    coordMap.outerCoordsFileName = coordMap.outerSurfaceFileName;
  end
end


  
