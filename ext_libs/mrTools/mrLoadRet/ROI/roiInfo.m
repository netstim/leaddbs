% roiInfo.m
%
%        $Id$
%      usage: roiInfo()
%         by: justin gardner
%       date: 01/16/08
%    purpose: 
%
function retval = roiInfo(v,roinum)

% check arguments
if ~any(nargin == [1 2])
  help roiInfo
  return
end

if ieNotDefined('roiNum')
  roiNum = viewGet(v,'currentROI');
end

if isempty(roiNum)
  mrWarnDlg('(roiInfo) No ROI currently selected.');
  return
elseif length(roiNum)>1
  mrWarnDlg('(roiInfo) This function cannot be used when several ROIs are selected');
  return
end

roiName = viewGet(v,'roiName',roiNum);
roiDate = viewGet(v,'roidate',roiNum);
roiColor = viewGet(v,'roicolor',roiNum);
roiVoxelSize = viewGet(v,'roivoxelsize',roiNum);
roiVolume = viewGet(v,'roivolume',roiNum);
roiSformCode = viewGet(v,'roiSformCode',roiNum);
roiXform = viewGet(v,'roixform',roiNum);
roiNotes = viewGet(v,'roiNotes',roiNum);
vol2mag = viewGet(v,'roiVol2mag',roiNum);
vol2tal = viewGet(v,'roiVol2tal',roiNum);
roiSubjectID = viewGet(v,'roiSubjectID',roiNum);
roiCreatedBy = viewGet(v,'roiCreatedBy',roiNum);
roiCreatedOnBase = viewGet(v,'roiCreatedOnBase',roiNum);
roiDisplayOnBase = viewGet(v,'roiDisplayOnBase',roiNum);
roiCreatedFromSession = viewGet(v,'roiCreatedFromSession',roiNum);


% check to see which base anatomy this roi aligns with
baseMatch = {};
for bnum = 1:viewGet(v,'numberOfBaseVolumes')
  % get the base voxelSize and xfrom
  baseVoxelSize = viewGet(v,'baseVoxelSize',bnum);
  baseXform = viewGet(v,'baseXform',bnum);
  % if it matches, then put it in thee list of matching base names
  if isequal(baseXform,roiXform) && isequal(baseVoxelSize,roiVoxelSize)
    baseMatch{end+1} = viewGet(v,'baseName',bnum);
  end
end
if isempty(baseMatch),baseMatch = 'No matching base anatomy';,end
if length(baseMatch)==1,baseMatch = baseMatch{1};end

paramsInfo = {{'name',roiName,'editable=0','The name of the ROI'},...
  {'notes',roiNotes,'editable=0','Notes associated with ROI'},...
  {'date',roiDate,'editable=0','The date of creation'},...
  {'color',roiColor,'editable=0','ROI color'},...
  {'voxelsize',roiVoxelSize,'editable=0','Voxel dimensions in mm'},...
  {'volume',roiVolume,'editable=0','Volume of ROI in cubic mm'},...
  {'sformCode',roiSformCode,'editable=0','Sform code for ROI. 1 means vol2mag 2 means vol2talarich'},...
  {'xform',roiXform,'editable=0','xform matrix specifies the transformation to the canoncial volume in magnet coordinates'},...
  {'vol2mag',vol2mag,'editable=0','xform matrix specifies the transformation of the canonical base volume to magnet coordinates'},...
  {'vol2tal',vol2tal,'editable=0','xform matrix specifies the transformation of the canonical base volume to talairach coordinates'},...
  {'baseMatch',baseMatch,'editable=0','The base volume that has the same voxel size and xform as this ROI. This is the base volume on which the ROI was originally defined. If there is no matching base anatomy, it means that the ROI was defined on a different base volume than the one you have loaded.'},...
  {'ROICoords',[],'type=pushbutton','buttonString=Show ROI coordinates','callback',@showCurrentROICoords,'callbackArg',v,'Print the coordinates for this ROI into the matlab window. Note that these will be the actual ROI coordinates not transformed into the scan coordinates. If you want the variable ROICoords set to the coordinates in your matlab workspace, you can hold the shift key down as you press this button.'},...
  {'ROIScanCoords',[],'type=pushbutton','buttonString=Show scan coordinates','callback',@showCurrentROIScanCoords,'callbackArg',v,'Print the coordinates transformed into the scan coordinates for thie ROI to the matlab window. If you want the variable ROICoords set to the coordinates in your matlab workspace, you can hold the shift key down as you press this button.'},...
  {'ROIBaseCoords',[],'type=pushbutton','buttonString=Show base coordinates','callback',@showCurrentROIBaseCoords,'callbackArg',v,'Print the coordinates transformed into the base coordinates for thie ROI to the matlab window. If you want the variable ROICoords set to the coordinates in your matlab workspace, you can hold the shift key down as you press this button.'},...
  {'subjectID',roiSubjectID,'SubjectID for this ROI. This is usually set by mlrAnatDB when you add the ROI to the DB','editable=0'},...
  {'createdBy',roiCreatedBy,'Person who created this ROI. Usually set by mlrAnatDB when you add the ROI to the DB','editable=0'},...
  {'createdFromSession',roiCreatedFromSession,'Session this was created from','editable=0'},...
  {'createdOnBase',roiCreatedOnBase,'The base anatomy this roi was created on','editable=0'},...
  {'displayOnBase',roiDisplayOnBase,'The base anatomy that this roi is best displayed on','editable=0'}};
% give ability to findROI for non baseCoordMapped ROIs
if viewGet(v,'baseType') == 0
  paramsInfo{end+1} = {'findROI',[],'type=pushbutton','buttonString=Find ROI','callback',@findROI,'callbackArg',v,'Go to the closest slice for which this ROI has some coordinates.'};
end

mrParamsDialog(paramsInfo,'ROI information',1.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper function, called by ROI Info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function retval = showCurrentROIScanCoords(view)

retval = [];

% get the roi
roiNum = viewGet(view,'currentROI');
if isempty(roiNum),return,end

% get the current scan number
scanNum = viewGet(view,'curScan');
% get the coordinates
coords = getROICoordinates(view,roiNum,scanNum);

% and display them to the buffer
disp(sprintf('ROI %s: n=%i',viewGet(view,'roiName',roiNum),size(coords,2)));
dispCoords(coords);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper function, called by ROI Info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function retval = showCurrentROIBaseCoords(view)

retval = [];

% get the roi
roiNum = viewGet(view,'currentROI');
if isempty(roiNum),return,end

% get the current scan number
scanNum = viewGet(view,'curScan');
% get the coordinates
coords = getROICoordinates(view,roiNum,0);

% and display them to the buffer
disp(sprintf('ROI %s: n=%i',viewGet(view,'roiName',roiNum),size(coords,2)));
dispCoords(coords);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper function, called by ROI Info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function retval = showCurrentROICoords(view)

retval = [];

% get the roi
roiNum = viewGet(view,'currentROI');
if isempty(roiNum),return,end

% just get roi coordinates
coords = viewGet(view,'ROICoords',roiNum);

% display to buffer
disp(sprintf('ROI %s: n=%i',viewGet(view,'roiName',roiNum),size(coords,2)));
dispCoords(coords);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper function, also called by ROI Info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dispCoords(coords)

% if shift is held down then just dump as an array
% that can be used
if any(strcmp(get(gcf,'CurrentModifier'),'shift'))
  % just dump as an array
  disp(sprintf('Setting variable ROICoords'));
  evalin('base',sprintf('ROICoords = [%s;%s;%s];',num2str(coords(1,:)),num2str(coords(2,:)),num2str(coords(3,:))));
  return
end
% and display them to the buffer
numCols = 20;
xline = 'x:';yline = 'y:';sline = 's:';colnum = 0;
for i = 1:size(coords,2)
  xline = sprintf('%s%4.0i',xline,coords(1,i));
  yline = sprintf('%s%4.0i',yline,coords(2,i));
  sline = sprintf('%s%4.0i',sline,coords(3,i));
  colnum = colnum + 1;
  if (colnum == numCols)
    disp(sprintf('Coordinates %i:%i',i-numCols+1,i));
    disp(xline);disp(yline);disp(sline);
    colnum = 0;
    xline = 'x:';yline = 'y:';sline = 's:';
  end

end
if colnum
    disp(sprintf('Coordinates %i:%i',i-colnum+1,size(coords,2)));
  disp(xline);disp(yline);disp(sline);
  colnum = 0;
end

