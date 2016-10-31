function view = combineROIs(view,roi1,roi2,action,newName)

% function view = combineROIs(view,roi1,roi2,action,[newName])
%
% Logical combination (union, intersection, xor, set difference) of ROIs.
% Modifies roi1 by combining it with roi2
%
% roi1 and roi2 can be either ROI names, ROI numbers, or empty (for current ROI).
% Default: current ROI
%
% if newName is specified, will create a new ROI with that name
% otherwise roi1 is modified
%
% action must be empty or one of the following strings: 
%     'Intersection', 'Union', 'XOR', 'A not B'
% If empty, then 'Intersection' is used as the default.
%
% djh 8/2007

nROIs = viewGet(view,'numberofROIs');
if ieNotDefined('roi1')
  roi1 = viewGet(view,'currentROI');
end
if isstr(roi1)
  roi1 = viewGet(view,'roiNum',roi1);
end
if ~isnumeric(roi1) | (roi1 < 1) | (roi1 > nROIs)
  mrErrorDlg('Invalid ROI');
end
if ieNotDefined('roi2')
  roi2 = viewGet(view,'currentROI');
end
if isstr(roi2)
  roi2 = viewGet(view,'roiNum',roi2);
end
if ~isnumeric(roi2) | (roi2 < 1) | (roi2 > nROIs)
  mrErrorDlg('Invalid ROI');
end
if ieNotDefined('action')
  action = 'Intersection';
end

% Get coordinates
roiVoxelSize1 = viewGet(view,'roiVoxelSize',roi1);
roiCoords1 = viewGet(view,'roiCoords',roi1);
roiVoxelSize2 = viewGet(view,'roiVoxelSize',roi2);
roiCoords2 = viewGet(view,'roiCoords',roi2);
roi2roi = viewGet(view,'roi2roi',roi2,roi1);

% Transform coords or roi2, using xformROIcoords to supersample the
% coordinates
if (size(roiCoords1,1)==3)
  roiCoords1(4,:) = 1;
end
if (size(roiCoords2,1)==3)
  roiCoords2(4,:) = 1;
end

roiCoords2 = round(xformROIcoords(roiCoords2,roi2roi,roiVoxelSize2,roiVoxelSize1));


% Transpose because matlab functions work on rows, not cols
coords1 = roiCoords1';
coords2 = roiCoords2';

action = lower(action);

% Combine
switch action
  case 'intersection'
     if isempty(coords1)
       newCoords = coords1;
     elseif isempty(coords2)
       newCoords = coords2;
     else       
       newCoords = intersect(coords1,coords2,'rows');
     end
  case 'union'
     if isempty(coords1)
       newCoords = coords2;
     elseif isempty(coords2)
       newCoords = coords1;
     else       
       newCoords = union(coords1,coords2,'rows');
     end
  case 'xor'
     if isempty(coords1)
       newCoords = coords2;
     elseif isempty(coords2)
       newCoords = coords1;
     else       
       newCoords = setxor(coords1,coords2,'rows');
     end
  case 'a not b'
     if isempty(coords1)
       newCoords = coords1;
     elseif isempty(coords2)
       newCoords = coords1;
     else       
       newCoords = setdiff(coords1,coords2,'rows');
     end
   otherwise
     error('unknown action: %s',action);
end

% Transpose back 
newCoords = newCoords';

% Select ROI and modify it
if ieNotDefined('newName')
  view = viewSet(view,'currentROI',roi1);
  % add coordinates, they are all roi1 coordinates, so we pass
  % the identity as the xform
  view = modifyROI(view,newCoords,eye(4),roiVoxelSize1,-1); %JB: added new option -1 to replace coordinates, so that Undo can be used
% make a new one
else
  roi = viewGet(view,'ROI',roi1);
  roi.name = newName;
  roi.coords = newCoords;
  view = viewSet(view,'newROI',roi);
  view = viewSet(view,'currentROI',viewGet(view,'roiNum',newName));
end
