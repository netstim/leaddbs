% baseInfo.m
%
%      usage: baseInfo(view)
%         by: justin gardner
%       date: 09/28/07
%    purpose: 
%
function baseInfo(v)

% check arguments
if ~any(nargin == [1])
  help baseInfo
  return
end

% get base info
scanNum = viewGet(v,'curScan');
groupNum = viewGet(v,'curGroup');
baseDims = viewGet(v,'baseDims');
baseQform = viewGet(v,'baseqform');
baseSform = viewGet(v,'baseSform');
baseSformCode = viewGet(v,'baseSformCode');
baseVolPermutation = viewGet(v,'baseVolPermutation');
baseVoxelSize = viewGet(v,'baseVoxelSize');
baseName = viewGet(v,'baseName');
baseCoordMap = viewGet(v,'baseCoordMap');
baseGamma = viewGet(v,'baseGamma');
baseRange = viewGet(v,'baseRange');
baseClip = viewGet(v,'baseClip');
baseType = viewGet(v,'baseType');
base2tal = viewGet(v,'base2tal');
vol2mag = viewGet(v,'baseVol2mag');
vol2tal = viewGet(v,'baseVol2tal');
baseAlpha = viewGet(v,'baseAlpha');
baseMultiDisplay = viewGet(v,'baseMultiDisplay');
baseDisplayOverlay = viewGet(v,'baseDisplayOverlay');

% set parameters
paramsInfo = {};
paramsInfo{end+1} = {'baseName',baseName,'editable=0','The name of the base anatomy'};
paramsInfo{end+1} = {'voxelSize',baseVoxelSize,'editable=0','Voxel dimensions in mm'};
paramsInfo{end+1} = {'baseDim',baseDims,'editable=0','Dimensions of base anatomy'};
paramsInfo{end+1} = {'qform',baseQform,'editable=0','Qform matrix specifies the transformation to the scanner coordinate frame'};
paramsInfo{end+1} = {'sform',baseSform,'editable=0','Sform matrix is set by mrAlign and specifies the transformation to the canonical base volume in the magnet coordinates if sform_code = 1, and to talairach if sform_code = 3'};
paramsInfo{end+1} = {'sform_code',baseSformCode,'editable=0','Sform code. This is 0 if the sform hase never been set. 1 if the sform is the xformation to the canonical base volume in magnet coordinates, and 3 if it is to talairach coordinates'};
paramsInfo{end+1} = {'vol2mag',vol2mag,'editable=0','This is the xformation that takes the canonical base that this base anatomy was aligned to into magnet coordinates. To set this field, you can use export from mrAlign.'};
paramsInfo{end+1} = {'vol2tal',vol2tal,'editable=0','This is the xformation that takes the canonical base that this base anatomy was aligned to into talairach coordinates. To set this field, you can export talairach coordinates from mrAlign.'};
paramsInfo{end+1} = {'clip',baseClip,'editable=0','Clip values for display'};
paramsInfo{end+1} = {'range',baseRange,'editable=0','Range of values in anatomy image'};
paramsInfo{end+1} = {'gamma',baseGamma,'editable=0','Gamma for display'};
paramsInfo{end+1} = {'baseType',baseType,'editable=0','Type of base. 0 = inplane. 1 = flat. 2 = surface'};
paramsInfo{end+1} = {'baseAlpha',baseAlpha,'editable=0','Alpha value with which this base will be drawn'};
paramsInfo{end+1} = {'baseDisplayOverlay',baseDisplayOverlay,'editable=0','Whether this base will display overlays - this is useful sometimes for fascicles and other bases that may not need to be shown with an overlay'};
paramsInfo{end+1} = {'baseMultiDisplay',baseMultiDisplay,'editable=0','Toggles whether this base will be shown with other bases at the same time - like when you show fascicles with a cortical surface at the same time'};

% add baseCoordMap info for flat files
if baseType == 1
  paramsInfo{end+1} = {'path',baseCoordMap.path,'editable=0','Directory from which this flat map was originally created'};
  paramsInfo{end+1} = {'flatFileName',baseCoordMap.flatFileName,'editable=0','Name of original off file from which this flat map was created'};
  paramsInfo{end+1} = {'innerCoordsFileName',baseCoordMap.innerCoordsFileName,'editable=0','Name of inner mesh (aka gray matter mesh) from which this flat map was created'};
  paramsInfo{end+1} = {'outerCoordsFileName',baseCoordMap.outerCoordsFileName,'editable=0','Name of outer mesh (aka white matter mesh) from which this flat map was created'};
  paramsInfo{end+1} = {'curvFileName',baseCoordMap.curvFileName,'editable=0','Name of curvature file from which this flat map was created'};
  paramsInfo{end+1} = {'anatFileName',baseCoordMap.anatFileName,'editable=0','Name of anatomy file from which the xform for this flat map was taken'};
  paramsInfo{end+1} = {'viewFlatOnSurface',[],'type=pushbutton','buttonString=View flat on surface','callback',@viewFlatOnSurface,'passParams=1','Click to view flat on the surface meshes'};
end

% add info for surfaces
if baseType == 2
  paramsInfo{end+1} = {'path',baseCoordMap.path,'editable=0','Directory from which this surface map was originally created'};
  paramsInfo{end+1} = {'innerSurfaceFileName',baseCoordMap.innerSurfaceFileName,'editable=0','Inner surface'};
  paramsInfo{end+1} = {'innerCoordsFileName',baseCoordMap.innerCoordsFileName,'editable=0','Inner surface coordinates'};
  paramsInfo{end+1} = {'outerSurfaceFileName',baseCoordMap.outerSurfaceFileName,'editable=0','Outer surface'};
  paramsInfo{end+1} = {'outerCoordsFileName',baseCoordMap.outerCoordsFileName,'editable=0','Outer surface coordinates'};
  paramsInfo{end+1} = {'curv',baseCoordMap.curvFileName,'editable=0','Name of curvature file from which this surface was created'};
  paramsInfo{end+1} = {'anatomy',baseCoordMap.anatFileName,'editable=0','Name of anatomy file from which the xform for this surface was taken'};
end

% bring up dialog
mrParamsDialog(paramsInfo,'Base anatomy information');

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   viewFlatOnSurface   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function params = viewFlatOnSurface(params,varargin);

thispwd = pwd;
if isdir(params.path)
  cd(params.path);
else
  mrWarnDlg(sprintf('Directory %s does not exist, please find the anatomy folder',params.path));
  pathStr = uigetdir(mrGetPref('volumeDirectory'),'Find anatomy directory');
  if pathStr == 0,return,end
  cd(pathStr);
end

mrFlatViewer(params.flatFileName,params.outerCoordsFileName,params.innerCoordsFileName,params.curvFileName,params.anatFileName);
cd(thispwd);