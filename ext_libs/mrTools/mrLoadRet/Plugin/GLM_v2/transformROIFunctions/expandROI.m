% expandROI.m
%
%        $Id: expandROI.m 1969 2010-12-19 19:14:32Z julien $ 
%      usage: transformedRoi = expandROI(roi,margin,<kernelType>)
%         by: julien besle
%       date: 11/01/2011
%
%    purpose: expands roi coordinates in 3D by margin voxels
%      input:   - margin: number of voxels that will be added (removed if negative) around the ROI.
% kernelType:   - (optional) type of kernel to convolve the ROI with: 
%                     sphere (default), cube, disc or square.
%                     the margin is the radius of the sphere/disc or the half-side of the square/cube 
%                     discs and squares are in the X-Y plane  

function roi = expandROI(roi,margin,kernelType)

if ~ismember(nargin,[2 3])
  help expandROI;
  return
end

if ieNotDefined('kernelType')
  kernelType = 'sphere';
end
if ~ismember(kernelType,{'sphere','cube','disc','square'})
  mrWarnDlg(['(expandROI) unknown kernel type ' kernelType]);
  roi=[];
  return
end

if numel(margin)==1
  switch(kernelType)
    case {'sphere','cube'}
      margin = [margin margin margin];
    case {'disc','square'}
      margin = [margin margin 0];
  end
elseif numel(margin)==2
  if ismember(kernelType,{'disc','square'})
    margin = [margin 0];
  else
    mrWarnDlg('(expandROI) margin must be a 1 or 3 element vector');
    roi=[];
    return
  end
elseif numel(margin)==3
  if ismember(kernelType,{'disc','square'})
    margin = margin(1:2);
    mrWarnDlg('(expandROI) Ignoring third element of margin');
  end
else
  mrWarnDlg('(expandROI) margin must be a 1, 2 or 3 element vector');
  roi=[];
  return
end
  
if diff(margin>0,2)
  mrWarnDlg('(expandROI) All elements of margin must have the same sign');
  roi=[];
  return
end
  

trim = any(margin<0);
margin = abs(margin);

boxCoords = [min(roi.coords(1:3,:),[],2)-margin'  max(roi.coords(1:3,:),[],2)+margin'];

%shift coordinates so that the boxes starts at 1 on all dimensions
voxelShift = -boxCoords(:,1)+1;

boxCoords = boxCoords+repmat(voxelShift,1,2);
roiCoords = roi.coords(1:3,:)+repmat(voxelShift,1,size(roi.coords,2));

volume = zeros(boxCoords(:,2)');
volume(sub2ind(boxCoords(:,2)',roiCoords(1,:),roiCoords(2,:),roiCoords(3,:)))=1;

if trim
  volume = 1-volume;
end

switch(kernelType)
  case 'sphere'
    [sphereX,sphereY,sphereZ] = ndgrid(-margin(1):margin(1),-margin(2):margin(2),-margin(3):margin(3));
    kernel = zeros(2*margin(1)+1,2*margin(2)+1,2*margin(3)+1);
    kernel(sqrt((sphereX./margin(1)).^2 + (sphereY./margin(2)).^2 + (sphereZ./margin(3)).^2)<1)=1;
  case 'cube'
%     kernel = zeros(2*margin(1)+2,2*margin(2)+2,2*margin(3)+2);
%     kernel(2:end-1,2:end-1,2:end-1)=1;
    kernel = ones(2*margin(1),2*margin(2),2*margin(3));
  case 'square'
    kernel = ones(2*margin(1),2*margin(2));
  case 'disc'
    [discX,discY] = ndgrid(-margin(1):margin(1),-margin(2):margin(2));
    kernel = zeros(2*margin(1)+1,2*margin(2)+1);
    kernel(sqrt((discX./margin(1)).^2 + (discY./margin(2)).^2)<1)=1;
end

volume = logical(convn(volume,kernel,'same'));
if trim
  volume = ~volume;
end
[newCoordsX,newCoordsY,newCoordsZ] = ind2sub(boxCoords(:,2)',find(volume));
roi.coords = [newCoordsX-voxelShift(1) newCoordsY-voxelShift(2) newCoordsZ-voxelShift(3)]';


