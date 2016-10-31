function saveOverlay(view,overlayName,analysisName,confirm)
%
%   saveOverlay(view,[overlayName],[analysisName],[confirm])
%
% Saves an overlay to a file. The filename = overlay.name.
%
% overlayName: Can be either the name or the number.
%          Default: current overlay.
% analysisname: Can be either the name or the number.
%          Default: current analysis
% confirm: If nonzero, prompts user to determine what to do if filename
%          already exists. Otherwise, relies of 'overwritePolicy'
%          preference to choose what to do. Default: 0.
%
% djh, 7/2004
% djh, 7/2006 updated with analysisName
% djh, 5/2007 added stuff to handle overwrite (copied from saveAnalysis)


if ieNotDefined('overlayName')
  overlayNum = viewGet(view,'currentOverlay');
  overlayName = viewGet(view,'overlayName',overlayNum);
end
if isstr(overlayName)
  overlayNum = viewGet(view,'overlayNum',overlayName);
elseif isnumeric(overlayName)
  overlayNum = overlayName;
  overlayName = viewGet(view,'overlayName',overlayNum);
else
  myErrorDlg(['Bad overlay name: ',overlayName]);
end

if ieNotDefined('analysisName')
  analysisNum = viewGet(view,'currentAnalysis');
  analysisName = viewGet(view,'analysisName',analysisNum);
end
if isstr(analysisName)
  analysisNum = viewGet(view,'analysisNum',analysisName);
elseif isnumeric(analysisName)
  analysisNum = analysisName;
  analysisName = viewGet(view,'analysisName',analysisNum);
else
  myErrorDlg(['Bad analysis name: ',analysisName]);
end

if ieNotDefined('confirm')
  confirm = 0;
end

% Path
filename = [analysisName,'-',overlayName,'.mat'];
pathStr = viewGet(view,'overlayDir');

% Check for over-writing
if isfile(fullfile(pathStr,filename))

  % get overwrite preference or ask user
  saveMethod = mrGetPref('overwritePolicy');
  saveMethodTypes = {'Merge','Rename','Overwrite'};
  saveMethod = find(strcmp(saveMethod,saveMethodTypes));
  if confirm || isempty(saveMethod) || (saveMethod==0)
    % ask the user what to do
    paramsInfo = {{'saveMethod',saveMethodTypes,'type=popupmenu','Choose how you want to save the variable'}};
    params = mrParamsDialog(paramsInfo,sprintf('%s already exists',filename));
    if isempty(params),return,end
    saveMethod = find(strcmp(params.saveMethod,saveMethodTypes));
  end

  if saveMethod == 1
    disp('(saveOverlay) Merging with old analysis');

    % load the old overlay
    view = loadOverlay(view,filename,pathStr,'tmpOverlay');
    oldOverlayNum = viewGet(view,'overlayNum','tmpOverlay',analysisNum);
    oldOverlay = viewGet(view,'overlay',oldOverlayNum,analysisNum);
    oldOverlayParams = viewGet(view,'overlayParams',oldOverlayNum,analysisNum);
    oldOverlayData = viewGet(view,'overlayData',[],oldOverlayNum,analysisNum);
    oldOverlayType = viewGet(view,'overlayType',oldOverlayNum,analysisNum);
    oldOverlayGroupName = viewGet(view,'overlayGroupName',oldOverlayNum,analysisNum);
    oldOverlayReconcileFunction = viewGet(view,'overlayReconcileFunction',oldOverlayNum,analysisNum);

    % get the new overlay
    newOverlay = viewGet(view,'overlay',overlayNum,analysisNum);
    newOverlayParams = viewGet(view,'overlayParams',overlayNum,analysisNum);
    newOverlayData = viewGet(view,'overlayData',[],overlayNum,analysisNum);
    newOverlayType = viewGet(view,'overlayType',overlayNum,analysisNum);
    newOverlayGroupName = viewGet(view,'overlayGroupName',overlayNum,analysisNum);
    newOverlayReconcileFunction = viewGet(view,'overlayReconcileFunction',overlayNum,analysisNum);
    overlayMergeFunction = viewGet(view,'overlayMergeFunction',overlayNum,analysisNum);

    % check if they have the same Type and merge them
    if strcmp(oldOverlayType,newOverlayType) && strcmp(oldOverlayGroupName,newOverlayGroupName)
      % reconcile them
      [oldOverlayParams oldOverlayData] = feval(oldOverlayReconcileFunction,oldOverlayGroupName,...
        oldOverlayParams,oldOverlayData);
      [newOverlayParams newOverlayData] = feval(newOverlayReconcileFunction,newOverlayGroupName,...
        newOverlayParams,newOverlayData);
      % and combine them
      [mergedParams,mergedData] = feval(overlayMergeFunction,newOverlayGroupName,...
        oldOverlayParams,newOverlayParams,oldOverlayData,newOverlayData);
      newOverlay.params = mergedParams;
      newOverlay.data = mergedData;
      
      % replace overlay with the newly merged one
      view = viewSet(view,'deleteoverlay',oldOverlayNum,analysisNum);
      view = viewSet(view,'deleteoverlay',overlayNum,analysisNum);
      view = viewSet(view,'newoverlay',newOverlay,analysisNum);
    else
      mrWarnDlg('(saveOverlay) Merge failed. Save overlay aborted.');
      return
    end

  elseif saveMethod == 2
    % put up a dialog to get new filename
    [filename pathStr] = uiputfile({'*.mat'},'Enter new name to save overlay',fullfile(pathStr,filename));
    % user hit cancel
    if isequal(filename,0)
      return
    end
    % otherwise accept the filename, make sure it has a .mat extension
    filename = [stripext(filename),'.mat'];

  elseif saveMethod == 3
    % this is the easiest, just overwrite
    disp('(saveAnalysis) Overwriting old analysis');
  end
end

% Assign local variable with overlayName = overlay for later saving
overlayNum = viewGet(view,'overlayNum',overlayName,analysisNum);
overlay = viewGet(view,'overlay',overlayNum,analysisNum);
overlayName = fixBadChars(overlayName);
eval([overlayName,'=overlay;']);

% Finally, write the file
pathStr = fullfile(pathStr,filename);
fprintf('Saving %s...',pathStr);
saveString = ['save(pathStr,','''',overlayName,'''',');'];
eval(saveString);
fprintf('done\n');

return
