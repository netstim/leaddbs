function saveAnalysis(view,analysisName,confirm)
%
%   saveAnalysis(view,[analysisName],[confirm])
%
% Saves an analysis to a file. The filename = analysis.name.
%
% analysisname: Can be either the name or the number.
%          Default: current analysis
% confirm: If nonzero, prompts user to determine what to do if filename
%          already exists. Otherwise, relies of 'overwritePolicy'
%          preference to choose what to do. Default: 0.
%
%
% djh, 7/2006
% jlg, 4/2007 added stuff to handle overwrite

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
analysisdir = viewGet(view,'analysisdir',[],analysisNum);
filename = [analysisName,'.mat'];
pathStr = analysisdir;

% Write, though check for over-writing
if isfile(fullfile(pathStr,filename))
  % get the preference how to deal with what to do with over-writing
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
    disp(sprintf('(saveAnalysis) Merging with old analysis'));

    % load the old analysis
    view = loadAnalysis(view,filename,pathStr,'tmpAnalysis');
    oldAnalNum = viewGet(view,'analysisNum','tmpAnalysis');
    oldAnal = viewGet(view,'analysis',oldAnalNum);
    oldType = viewGet(view,'analysisType',oldAnalNum);
    oldGroupName = viewGet(view,'analysisGroupName',oldAnalNum);
    oldParams = viewGet(view,'analysisParams',oldAnalNum);
    oldData = viewGet(view,'analysisData',oldAnalNum);
    numOldOverlays = viewGet(view,'numberofOverlays',oldAnalNum);
    oldReconcileFunction = viewGet(view,'reconcileFunction',oldAnalNum);

    % get the new analysis
    newAnal = viewGet(view,'analysis',analysisNum);
    newType = viewGet(view,'analysisType',analysisNum);
    newGroupName = viewGet(view,'analysisGroupName',analysisNum);
    newParams = viewGet(view,'analysisParams',analysisNum);
    newData = viewGet(view,'analysisData',analysisNum);
    numNewOverlays = viewGet(view,'numberofOverlays',analysisNum);
    newReconcileFunction = viewGet(view,'reconcileFunction',analysisNum);
    mergeFunction = viewGet(view,'mergeFunction',analysisNum);

    % Check if they have the same Type and merge them
    if strcmp(oldType,newType) && strcmp(oldGroupName,newGroupName)
      [oldParams oldData] = feval(oldReconcileFunction,oldGroupName,oldParams,oldData);
      [newParams newData] = feval(newReconcileFunction,newGroupName,newParams,newData);
      [mergedParams,mergedData] = feval(mergeFunction,newGroupName,oldParams,newParams,oldData,newData);
      newAnal.params = mergedParams;
      newAnal.d = mergedData;

      % loop through overlays and merge the params and data
      for oldOverlayNum = 1:numOldOverlays
        oldOverlay = viewGet(view,'overlay',oldOverlayNum,oldAnalNum);
        oldOverlayParams = viewGet(view,'overlayParams',oldOverlayNum,oldAnalNum);
        oldOverlayData = viewGet(view,'overlayData',[],oldOverlayNum,oldAnalNum);
        oldOverlayType = viewGet(view,'overlayType',oldOverlayNum,oldAnalNum);
        oldOverlayName = viewGet(view,'overlayName',oldOverlayNum,oldAnalNum);
        oldOverlayGroupName = viewGet(view,'overlayGroupName',oldOverlayNum,oldAnalNum);
        oldOverlayReconcileFunction = viewGet(view,'overlayReconcileFunction',oldOverlayNum,oldAnalNum);
        matchedOverlay = 0;

        for newOverlayNum = 1:numNewOverlays
          newOverlay = viewGet(view,'overlay',newOverlayNum,analysisNum);
          newOverlayParams = viewGet(view,'overlayParams',newOverlayNum,analysisNum);
          newOverlayData = viewGet(view,'overlayData',[],newOverlayNum,analysisNum);
          newOverlayType = viewGet(view,'overlayType',newOverlayNum,analysisNum);
          newOverlayGroupName = viewGet(view,'overlayGroupName',newOverlayNum,analysisNum);
          newOverlayReconcileFunction = viewGet(view,'overlayReconcileFunction',newOverlayNum,analysisNum);
          overlayMergeFunction = viewGet(view,'overlayMergeFunction',newOverlayNum,analysisNum);

          % see if there is a matching overlay
          if strcmp(oldOverlayType,newOverlayType) && strcmp(oldOverlayGroupName,newOverlayGroupName)
            matchedOverlay = 1;
            % reconcile them
            [oldOverlayParams oldOverlayData] = feval(oldOverlayReconcileFunction,oldOverlayGroupName,...
              oldOverlayParams,oldOverlayData);
            [newOverlayParams newOverlayData] = feval(newOverlayReconcileFunction,newOverlayGroupName,newOverlayParams,newOverlayData);
            % and combine them
            [mergedParams,mergedData] = feval(overlayMergeFunction,newOverlayGroupName,...
              oldOverlayParams,newOverlayParams,oldOverlayData,newOverlayData);
            newOverlay.params = mergedParams;
            newOverlay.data = mergedData;
            newAnal.overlays(newOverlayNum) = newOverlay;
            disp(sprintf('(saveAnalysis) Merged overlay %s from old analysis',oldOverlayName));
          end
        end
        % if there was no match, then simply add the overlay to the analysis struct
        if ~matchedOverlay
          %view = viewSet(view,'newoverlay',oldOverlay,analysisNum);
	  % we have added the overlay to the one in the mlr structure,
	  % but below we are using newAnal to set the new analysis,
	  % so we retrieve the overlays and stick it on to our
          % newAnal structure
	  %updatedAnalysis = viewGet(view,'analysis',analysisNum);
	  %newAnal.overlays = updatedAnalysis.overlays;
	  
	  % this used to do the above to tack on an old overlay, but that
	  % didn't work if there was a combination of old overlays and
	  % new overlays (since it overwrote newAnal.overlays, it undid
	  % anything that was done above when you found a matchedOverlay).
	  % so now we just tack on the overlay. Note that it should
	  % already be reconciled, since the old analysis was installed
	  % into MLR. So this should be the same as above -jg.
	  newAnal.overlays(end+1) = oldOverlay;
          disp(sprintf('(saveAnalysis) Added overlay %s from old analysis',oldOverlayName));
        end
      end

      % replace analysis with the newly merged one
      view = viewSet(view,'deleteAnalysis',oldAnalNum);
      view = viewSet(view,'deleteAnalysis',analysisNum);
      view = viewSet(view,'newAnalysis',newAnal);
    else
      mrWarnDlg('(saveAnalysis) Merge failed. Save analysis aborted.');
      return
    end

  elseif saveMethod == 2
    % put up a dialog to get new save name
    [filename pathStr] = uiputfile({'*.mat'},'Enter new name to save analysis as',fullfile(pathStr,filename));
    % user hit cancel
    if isequal(filename,0)
      return
    end
    % otherwise accept the filename, make sure it has a .mat extension
    filename = sprintf('%s.mat',stripext(filename));
    % change the analysis name
    view = viewSet(view,'analysisName',stripext(filename),analysisNum);
    analysisName = viewGet(view,'analysisName',analysisNum);
  elseif saveMethod == 3
    % this is the easiest, just overwrite
    disp(sprintf('(saveAnalysis) Overwriting old analysis'));
  end
end

% Assign local variable with analysisName = analysis for later saving
analysisNum = viewGet(view,'analysisNum',analysisName);
analysis = viewGet(view,'analysis',analysisNum);
analysisName = fixBadChars(analysisName);
eval([analysisName,'=analysis;']);

% Finally, write the file
pathStr = fullfile(pathStr,filename);
fprintf('Saving %s...',pathStr);
if getfield(whos('analysis'),'bytes')<2e9
  saveString = ['save(pathStr,','''',analysisName,'''',');'];
else %if the structure is more than 2Gb
  %check if matlab/the computer is 64bit
  architecture = computer('arch');
  if str2num(architecture(end-1:end))~=64
    mrErrorDlg('(saveAnalysis) Analysis is more than 2Gb, but Matlab is 32bit. Cannot save analysis.');
  else
    mrWarnDlg('(saveAnalysis) Analysis is more than 2Gb, using option -v7.3 to save');
    saveString = ['save(pathStr,','''',analysisName,'''',',''-v7.3'');'];
  end
end
  
eval(saveString);
fprintf('done\n');

return;
