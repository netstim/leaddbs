function [view  userCancel] = newROI(view,name,select,color,xform,voxelSize,coords,xformCode,vol2tal,vol2mag)

% function view = newROI(view,[name],[select],[color],[xform],[voxelSize],[coords],[xformCode],[vol2tal],[vol2mag])
%
% Makes new empty ROI, adds it to view.ROIs, and selects it.
%
% name: name (string) for the ROI.
%   Default: 'ROI<number>' where <number> is chosen (to avoid conflicts with
%   other ROI names) as one plus the number of previously existing ROIs.
% select: if non-zero, chooses the new ROI as the selectedROI
%    Default: 1.
% color: sets color for drawing the ROI.
%    Default: 'b'
% xform: 4x4 transform from functionals to ROI. 
%    Default:  uses the xform from the current base anatomy
% voxelSize: 3-vector specifying the voxel size of the ROI
%    Default: uses the voxel size of the current base anatomy
% coords: 4xN array of ROI coordinates (bottom row filled with ones) in the
%    reference frame of the xform.
%    Default: []
%
% djh, 7/2005 (modified from mrLoadRet-3.1)

userCancel = 1;
if isempty(viewGet(view,'curBase')) & ieNotDefined('xform') & ieNotDefined('voxelSize')
  mrErrorDlg('You must load a base anatomy before creating an ROI.');
end

if ieNotDefined('name')
  % go through roi names and get the largest numbered
  % roi name, i.e. ROI4 then make then new name ROI5
  maxnum = 0;
  for i = 1:length(view.ROIs)
    if regexp(view.ROIs(i).name,'^ROI\d+$')
      maxnum = max(maxnum,str2num(view.ROIs(i).name(4:end)));
    end
  end
  name=sprintf('ROI%.0f',maxnum+1);
end
if ieNotDefined('select')
  select = 1;
end
if ieNotDefined('color')
  color = 'black';
end
if ieNotDefined('sformCode')
  baseNum = viewGet(view,'currentBase');
  sformCode = viewGet(view,'baseSformCode',baseNum);
end
if ieNotDefined('xform')
  baseNum = viewGet(view,'currentBase');
  % if the baseSformCode == 0 then use the bases Qform matrix
  % since it means that the alignment has not been set.
  if viewGet(view,'baseSformCode',baseNum) == 0
    xform = viewGet(view,'baseQform',baseNum);
    sformCode = 0;
  else
    xform = viewGet(view,'baseSform',baseNum);
  end
end
if ieNotDefined('voxelSize')
  baseNum = viewGet(view,'currentBase');
  voxelSize = viewGet(view,'baseVoxelSize',baseNum);
end
if ieNotDefined('coords')
  coords = [];
end
if ieNotDefined('vol2mag')
  baseNum = viewGet(view,'currentBase');
  vol2mag = viewGet(view,'baseVol2mag',baseNum);
end
if ieNotDefined('vol2tal')
  baseNum = viewGet(view,'currentBase');
  vol2tal = viewGet(view,'baseVol2tal',baseNum);
end

colors = putOnTopOfList(color,color2RGB);
roiParams{1} = {'name',name,'Name of roi, avoid using punctuation and space'};
roiParams{2} = {'color',colors,'type=popupmenu','The color that the roi will display in'};
roiParams{3} = {'notes','','Brief notes about the ROI'};
params = mrParamsDialog(roiParams,'Create a new ROI');
if isempty(params),return,end

% Set required fields. Additional (optional) optional fields are set by
% isroi which is called by viewSet newROI.
ROI.name = params.name;
ROI.viewType = view.viewType;
ROI.color = params.color;
ROI.xform = xform;
ROI.sformCode = sformCode;
ROI.voxelSize = voxelSize;
ROI.coords = coords;
ROI.notes = params.notes;
ROI.vol2mag = vol2mag;
ROI.vol2tal = vol2tal;

% add some additional information about where this was created
baseName = viewGet(view,'basename');
ROI.createdOnBase = baseName;
ROI.displayOnBase = baseName;
ROI.createdFromSession = getLastDir(viewGet(view,'homeDir'));

% Add it to the view
[view tf]= viewSet(view,'newROI',ROI);
 
% The user could still have canceled (when there is a name conflict)
% so check for that
if ~tf,return,end

% also get the ROI back, because the name may have changed (if
% there was a conflict)
ROInum = viewGet(view,'nROIs');

% Select it and reset view.prevCoords
if select
  if (ROInum > 0)
    view = viewSet(view,'currentROI',ROInum);
    view = viewSet(view,'prevROIcoords',[]);
  end
end
userCancel = 0;
return;
