function view = modifyROI(view,coords,xform,voxelSize,sgn)
%
%        $Id$
%
% view = modifyROI(view,coords,xform,voxelSize,[sgn])
%
% Adds/deletes coordinates to/from the current ROI.
%
% coords: 4xN array of ROI coordinates (bottom row filled with ones) in the
% reference frame of the xform.
%
% xform: 4x4 homogeneous transform matrix that transforms the coords to the
% current rois coordinate reference frame. 
%
% voxelSize: 3-vector providing voxel size of new coordinates in mm. This
% is used along with the voxel size of the ROI when transforming the
% coordinates.
%
% sgn: If sgn = 1 [default], adds user-specified coordinates to selected
% ROI in current slice. If sgn = 0, removes those coordinates from the
% ROI. if sgn=-1, replaces coordinates
% djh, 7/2005 (modified from mrLoadRet-3.1)


% Error if no current ROI
ROInum = viewGet(view,'currentROI');
if isempty(ROInum)
	mrErrorDlg('No ROI currently selected.');
end

% Get current ROI coords
curCoords = viewGet(view,'roiCoords',ROInum);

% Save prevCoords for undo
view = viewSet(view,'prevROIcoords',curCoords);

% Transform to ROI reference frame
roiVoxelSize = viewGet(view,'roiVoxelSize',ROInum);
newCoords = xformROIcoords(coords,xform,voxelSize,roiVoxelSize);

% Merge/remove coordinates
if sgn<0
  coords=newCoords;
elseif sgn==0
  coords = removeCoords(newCoords,curCoords);
elseif sgn>0
  coords = mergeCoords(curCoords,newCoords);
end
view = viewSet(view,'ROIcoords',coords);

return;


function coords = mergeCoords(coords1,coords2)
%
% coords = mergeCoords(coords1,coords2)
%
% Merges coords1 and coords2, removing duplicates, used for
% example to merge ROI coordinates into one big ROI.
%
% coords, coords1, and coords2: 3xN arrays of (x,y,z) coordinates
% dims is size of volume
%
% djh, 7/98
% djh, 2/2001, dumped coords2Indices & replaced with union(coords1',coords2','rows')

if ~isempty(coords1) & ~isempty(coords2)
    coords = union(coords1',coords2','rows');
    coords = coords';
else
    coords = [coords1 coords2];
end

return


function coords = removeCoords(coords1,coords2)
%
% coords = removeCoords(coords1,coords2)
%
% Removes coords1 from coords2, used for example to remove
% coordinates from an ROI.
%
% coords, coords1, and coords2: 3xN arrays of (x,y,z) coordinates
% dims is size of volume
%
% djh, 7/98
% djh, 2/2001, dumped coords2Indices & replaced with setdiff(coords1',coords2','rows')

if ~isempty(coords1) && ~isempty(coords2)
    coords = setdiff(coords2',coords1','rows');
    coords = coords';
else
    coords=coords2;
end
return



% Test
view = MLR.views{1};
coords = [round(viewGet(view,'baseDims')/2)';1];
xform = viewGet(view,'baseXform');
%voxelSize = viewGet(view,'baseVoxelSize');
voxelSize = [4 4 4];
view = modifyROI(view,coords,xform,voxelSize,0);
viewGet(view,'roiCoords')


