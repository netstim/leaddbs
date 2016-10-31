function view = loadROI(view,filename,select,startPathStr)
%
% view = loadROI(view,[filename],[select],[startPath])
%
% Loads an ROI and adds it to view.ROIs.
%
% If filename is not specified, prompts user to select a file. If filename
% is specified, it loads from view.subdir/ROIs/name.mat. Filename can
% be a string specifying an ROI file or it can be a cell array of
% filenames to load multiple ROIs at once.
%
% select: If non-zero, selects the loaded ROI as the current (selected)
% ROI. Default: 1.
%
% The file must contain a structure or structures, each with the following
% fields:
% - name: string
% - viewType: 'Volume', 'Surface', or 'Flat'
% - color:
% - coords: 4xN array of ROI coordinates
% - xform: 4x4 homogeneous transform matrix that maps roi coordinates to base coordinate frame.
% - voxelSize: 3-vector specifying the ROI voxel size.
%
% In addition the 'name' field is set to the variable name to ensure
% consistency.
%
% djh, 1/9/98
% 8/2005, djh, updated to mrLoadRet-4.0

mrGlobals

if ieNotDefined('select')
    select = 1;
end

% Path to roi
if ieNotDefined('startPathStr')
  startPathStr = viewGet(view,'roiDir');
end

% Complete pathStr
if ieNotDefined('filename')
    pathStr = mlrGetPathStrDialog(startPathStr,'Choose one or more ROIs','*.mat','on');
else
    if iscell(filename)
        pathStr = cell(size(filename));
        for p=1:length(pathStr)
          pathStr{p} = fullfile(startPathStr,filename{p});
          pathStr{p} = [stripext(pathStr{p}),'.mat'];
        end
    else
        pathStr = {[stripext(fullfile(startPathStr,filename)),'.mat']};
    end
end
if isempty(pathStr),return,end
if ~iscell(pathStr)
    pathStr = {pathStr};
end

% check for any name collision here
roiNames = viewGet(view,'roiNames');
nameCollision = 0;nameCollisionStr = '';
for p = 1:length(pathStr)
  loadNames{p} = stripext(getLastDir(pathStr{p}));
  if any(strcmp(loadNames{p},roiNames))
    nameCollision(p) = 1;
    nameCollisionStr = sprintf('%s ''%s''',nameCollisionStr,loadNames{p});
  else
    nameCollision(p) = 0;
  end
end
replaceDuplicates = 0;
if any(nameCollision)
  paramsInfo{1} = {'nameCollision',{'Replace all','Deal with individually'},'type=popupmenu',sprintf('ROIs %s are already loaded. For these ROIs, you can either just replace with the new ones from disk, or be asked individually what to do about them.',nameCollisionStr)};
  params = mrParamsDialog(paramsInfo,sprintf('Some ROIs are already loaded',nameCollisionStr));
  if isempty(params),return,end
  if strcmp(params.nameCollision,'Replace all')
    replaceDuplicates = 1;
  end
end
% Load the file. Loop through the variables that were loaded and add
% each of them as a new ROI, setting roi.fieldnames as we go.
for p = 1:length(pathStr)
    if exist(pathStr{p},'file')
        s = load(pathStr{p});
        varNames = fieldnames(s);
        roi = eval(['s.',varNames{1}]);
        roi.name = varNames{1};
        % Add it to the view
        [view tf] = viewSet(view,'newROI',roi,replaceDuplicates);
	if tf
	  % Select it and reset view.prevCoords
	  if select
            ROInum = viewGet(view,'numberofROIs');
            if (ROInum > 0)
	      view = viewSet(view,'currentROI',ROInum);
            end
	  end
        end
    else
        mrWarnDlg(['ROI ',pathStr{p},' not found.']);
    end
end

return;

% Test/debug
view = loadROI(MLR.views{1},'ROI1');
view = loadROI(MLR.views{1},{'ROI1','ROI2'});
view = loadROI(MLR.views{1});

