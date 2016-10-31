function pathStr = saveROI(view,roiName,confirm,pathStr)
%
% saveROI(view,[roiName],[confirm],[pathStr])
%
% Saves an ROI to a file. The filename = ROI.name.
%
% roiName: Either the ROI name or number.
%          Default: current overlay.
% confirm: If filename already exists, prompt user to over-write.
%          Default: uses mrLoadRet 'verbose' preference or 0 (if preference
%          not defined.
%
%          If pathStr is specified, will store ROI in that 
%          directory instead of the default session one:
%          e.g.
%          saveROI(v,'V1',0,'~/data/mlrAnatDB/s0001/ROIs')
%
% djh, 9/2005

pathStr = [];
if ieNotDefined('roiName')
	roiNum = viewGet(view,'currentROI');
	roiName = viewGet(view,'roiName',roiNum);
end

if ieNotDefined('confirm')
    pref = mrGetPref('overwritePolicy');
    confirm = ~strcmp(pref,'Overwrite');
end

% support saving multiple ROIs at once
if iscell(roiName)
    for r = 1:length(roiName)
        % recurse
        saveROI(view,roiName{r},confirm);
    end
    % and finish
    return
end

if isstr(roiName)
    roiNum = viewGet(view,'roiNum',roiName);
elseif isnumeric(roiName)
	roiNum = roiName;
	roiName = viewGet(view,'roiName',roiNum);
else
	myErrorDlg(['Bad ROI name: ',roiName]);
end

% Assign local variable with roiName = roi
roi = viewGet(view,'roi',roiNum);
% fix characters that are not allowed in variable names
roiName = fixBadChars(roiName);
eval([roiName,'=roi;']);

% path to file
filename = roiName;
if ieNotDefined('pathStr')
  pathStr = fullfile(viewGet(view,'roiDir'),filename);
else
  pathStr = fullfile(pathStr,filename);
end

% Write, though check for over-writing
%
saveFlag = 'Yes';
if exist([pathStr,'.mat'],'file')
	if confirm
		saveFlag = questdlg([pathStr,' already exists. Overwrite?'],...
			'Save ROI?','Yes','No','No');
	end
end
if strcmp(saveFlag,'Yes')
	fprintf('Saving %s...',pathStr);
	saveString = ['save(pathStr,','''',roiName,'''',');'];
	eval(saveString);
	fprintf('done\n');
else
	fprintf('ROI not saved...');
end
return;

% Test/debug
saveROI(MLR.views{1},'ROI1');
saveROI(MLR.views{1});
