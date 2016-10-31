% exportTal2mrLR.m
%
%        $Id: exportTal2mrLR.m,v 1.0 2008/01/10 
%      usage: exportTal2mrLR(vol2tal, vol2mag)
%         by: shani offen, based on saveSform by justin gardner
%       date: 2008-Jan-10
%    purpose: saves vol2tal and vol2mag xforms to mrLoadRet-4.5 mrSession variable 
%             and to the base, scanparams, and ROIs. Is based on JG's saveSform code
%
% as of 2008Jan25, it doesn't do anything
% it still needs to export to scan, roi, and base

function retval = exportTal2mrLR(vol2tal, vol2mag, talInfo)

% check arguments
if ~(nargin == 3)
  help exportTal2mrLR
  return
end

% first check this directory for mrSession.m
path = '';
if isfile('mrSession.mat')
  % if there is a session then ask if the user wants to export to this directory
  answer = questdlg(sprintf('Save Talairach information to %s?',getLastDir(pwd)),'Export');
  if strcmp(answer,'Cancel'),return,end % Cancel~=No; if they return No, the next 'if' loop takes over.
  if strcmp(answer,'Yes') 
    path = pwd;
  end
end

% if path is not set then ask the user to choose a path
if isempty(path)
  [filename path] = uigetfile('*.mat','Choose a mrSession.mat file to save Talairach information to');
  if filename == 0,return,end
  if ~strcmp(filename,'mrSession.mat')
    mrWarnDlg(sprintf('(exportTal2mrLR) %s is not an mrSession file',filename));
    return
  end
end

% start a view in the corresponding location
cd(path);
v = newView('Volume');
mrGlobals
% just make sure that the home dir matches
if ~strcmp(MLR.homeDir,path)
  answer = questdlg(sprintf('mrLoadRet is open on session %s? Ok to close and open on %s',...
                            getLastDir(MLR.homeDir),getLastDir(path)));
  if ~strcmp(answer,'Yes')
    mrWarnDlg(sprintf('(exportTal2mrLR) Could not open a view to %s',getLastDir(path)));
    deleteView(v);
    return
  end
  % clear MLR and start over
  deleteView(v);
  clear global MLR;
  v = newView('Volume');
end

% ask the user what to export to (bases, rois, and/or groups)
listOps = {'Bases_','ROIs_','EPIs_'};
for i = 1:3
  paramsInfo{i}{1} = listOps{i};
  paramsInfo{i}{2} = 1;
  paramsInfo{i}{3} = 'type=checkbox';
  paramsInfo{i}{4} = sprintf('Export Talairach information to %s',listOps{i});
end
paramsInfo{3}{4} = 'You will have a chance to choose groups later';
talParams = mrParamsDialog(paramsInfo,'Choose types of files to export Talairach information to');
if isempty(talParams),return,end
clear paramsInfo

% Go through the options one at a time




%***************************
%     EPIs
%***************************

% save to ***SCANS***
if(talParams.EPIs_)
  % ask the user which groups to export to
  groupNames = viewGet(v,'groupNames');
  for i = 1:length(groupNames)
    paramsInfo{i}{1} = groupNames{i};
    paramsInfo{i}{2} = 1;
    paramsInfo{i}{3} = 'type=checkbox';
    paramsInfo{i}{4} = sprintf('Export Talairach information to group %s',groupNames{i});
  end
  groupParams = mrParamsDialog(paramsInfo,'Choose groups to export Talairach information to');
  if isempty(groupParams),return,end
  clear paramsInfo

  % now go through and update the session parameters
  for iGroup = 1:viewGet(v, 'numberofGroups')
    % see if we are supposed to update the group
    if groupParams.(viewGet(v,'groupName',iGroup))
      for iScan = 1:viewGet(v, 'nScans', iGroup);

        % check if anything is being replaced:
        sform_code = viewGet(v,'scansformcode');        
        curScanVol2Tal = viewGet(v, 'scanVol2Tal', iScan, iGroup);
        curScanVol2Mag = viewGet(v, 'scanVol2Mag', iScan, iGroup);
          
        % * First skip scans with major problems *
        
        % if the sform_code == 3, then this scan was aligned to a base with
        % a tal transform, so changing it is complicated. Right now we're not
        % going to deal with that situation, but we might later
        if sform_code == 3 % if was originally aligned to the Vol with old talXform
          mrWarnDlg(sprintf(['(exportTal2mrLR) Scan %i in group %s  has been '...
                             'directly aligned to a base'... 
                    ' with a Talairach transformation (e.g., its sform_code is 3),'...
                    ' so changing its TalInfo will create too many problems. We may'...
                    ' be able to support this type of change later, but for now, you'...
                    ' need to stick with the transform youve got.'],iScan, groupNames{iGroup}));
          % I had actually already written some code to deal with this, which I'll
          % stick in a subfunction below so we can try to implement it later.
          skipScan = 1;
          
        else % The other thing to check is if the scan was aligned to some other volume
          
          if and(~isequal(curScanVol2Mag,vol2mag),~isempty(curScanVol2Mag))
            sameVol = questdlg(sprintf(['Scan %i in %s group is aligned to a different volume. '...
                               'This is almost certainly an error. Please indicate whether to continue '...
                               'anyway, or skip saving the Talairach transform to this scan (recommended).'],...
                               iScan,groupNames{iGroup}),...
                               'Skip this scan?','Skip','Continue','Abort','Skip');
            if any(strcmp(sameVol,'Skip'),strcmp(sameVol,''))
              skipScan = 1;
            elseif strcmp(sameVol,'Continue')
              skipScan = 0;
              disp(sprintf('vol2mag for scan %i in %s group has been changed',iScan,groupNames{iGroup}))              
            elseif strcmp(sameVol,'Abort')
              mrWarnDlg('(exportTal2mrLR) Exiting and not saving the Talairach information.')
              return
            end
          else % if sform_code = 1 and no vol2mag error, don't skip the scan
            skipScan = 0;
          end % checking for a vol2mag mismatch
        end % checking for problems that mean we should skip the scan
        
        if ~skipScan % if we don't need to skip the scan, do some more checks and add the talInfo
            
          % update user on what (if anything)  is changing:
          % first check if there already is a vol2tal
          if ~isequal(curScanVol2Tal,vol2tal) %check if there's a change
            if isempty(curScanVol2Tal) % if old vol2tal was empty, then it's easy
              disp(sprintf('vol2tal for scan %i in %s group has been added',iScan,groupNames{iGroup}));
            else  
              disp(sprintf('vol2tal for scan %i in %s group has been changed',iScan,groupNames{iGroup}));
            end
          end
          % also check if adding a new vol2mag
          if isempty(curScanVol2Mag)
            disp(sprintf('vol2mag for scan %i in %s group has been added',iScan,groupNames{iGroup}))
          end
          
          % finally, make the changes!
          v = viewSet(v,'scanvol2tal',vol2tal,iScan,iGroup); % set the fields
          v = viewSet(v,'scanvol2mag',vol2mag,iScan,iGroup);
          v = viewSet(v,'scanTalInfo',talInfo,iScan,iGroup)' %******* I think I need to write this! ********
          
        end % checking if need to skip the scan
        clear skipScan curScanVol2Tal curScanVol2Mag sform_code
      end % for iScan
    end % if this group gets updated 
  end % for iGroup
end % checking if user wanted to save to scans 
% save the session
saveSession

%***************************
%     ROIs 
%***************************


if(talParams.ROIs_)
  % let user choose ROIs
  
  v = loadROI(v);

  for iROI = 1:viewGet(v, 'nrois');

    % check if anything is being replaced:
    ROIsform_code = viewGet(v, 'roisformcode', iROI);
    curRoiVol2Tal = viewGet(v, 'roiVol2Tal', iROI);
    curRoiVol2Mag = viewGet(v, 'roiVol2Mag', iROI);
          
    % Check if ROI was defined relative to a different base and needs to be skipped
    if ROIsform_code == 3 % if was originally aligned to the Vol with old talXform
      mrWarnDlg(sprintf(['(exportTal2mrLR) ROI %s has been defined on a base '...
                         ' with a Talairach transformation (e.g., its sform_code is 3),'...
                         ' so changing its TalInfo will create too many problems. We may'...
                         ' be able to support this type of change later, but for now, you'...
                         ' need to stick with the transform youve got.'],viewGet(v,'roiname',iROI)));
      % I had actually already written some code to deal with this, which I'll
      % stick in a subfunction below so we can try to implement it later.
      skipRoi = 1;
          
    else % The other thing to check is if the ROI was defined on other volume
          
      if and(~isequal(curRoiVol2Mag,vol2mag),~isempty(curRoiVol2Mag))
        sameVol = questdlg(sprintf(['ROI %s is defined relative to a different volume. '...
                            'This is almost certainly an error. Please indicate whether to continue '...
                            'anyway, or skip saving the Talairach transformation to this this ROI (recommended).'],...
                            viewGet(v,'roiname',iROI)),...
                           'Skip this ROI?','Skip','Continue','Abort','Skip');
        if any(strcmp(sameVol,'Skip'),strcmp(sameVol,''))
          skipRoi = 1;
        elseif strcmp(sameVol,'Continue')
          skipRoi = 0;
          disp(sprintf('vol2mag for ROI %s has been changed',viewGet(v,'roiname',iROI)))              
        elseif strcmp(sameVol,'Abort')
          mrWarnDlg('(exportTal2mrLR) Exiting and not saving the Talairach information to ROIs or bases.')
          return
        end
      else % if sform_code = 1 and no vol2mag error, don't skip the ROI
        skipRoi = 0;
      end % checking for a vol2mag mismatch
    end % checking for problems that mean we should skip the ROI
    
    if ~skipRoi % if we don't need to skip the scan, do some more checks and add the talInfo
      
      % update user on what (if anything)  is changing:
      % first check if there already is a vol2tal
      if ~isequal(curRoiVol2Tal,vol2tal) %check if there's a change
        if isempty(curRoiVol2Tal) % if old vol2tal was empty, then it's easy
          disp(sprintf('vol2tal for ROI %s has been added',viewGet(v,'roiname',iROI)));
        else  
          disp(sprintf('vol2tal for ROI %s has been modified',viewGet(v,'roiname',iROI)));
        end
      end
      % also check if adding a new vol2mag
      if isempty(curRoiVol2Mag)
        disp(sprintf('vol2mag for ROI %s has been added',viewGet(v,'roiname',iROI)));
      end
      
      % finally, make the changes!
      v = viewSet(v,'roivol2tal',vol2tal,iROI); % set the fields
      v = viewSet(v,'roivol2mag',vol2mag,iROI);
      saveROI(v,iROI,0);
      
    end % checking if need to skip the ROI
    clear skipRoi curRoiVol2Tal curRoiVol2Mag ROIsform_code    
  end % for iROI = 1:viewGet(v, 'nrois');
end % checking if user wanted to save to ROIs 

%***************************
%     Bases
%***************************


if(talParams.Bases_)
  % let user choose bases - but bases seem to have to be loaded one at a time
  
  disp('You will now choose base anatomies one at a time to save the Talaiarach information to.')
  basesOver = 0;
  while ~basesOver
    askBase = questdlg('Would you like to choose another base to add the Talaiarch transform to?',...
                       'Save Talairach transform to another base?');
    if ~strcmp(askBase,'Yes'), basesOver = 1; end
    
    v = loadAnat(v);
  end
  
  if viewGet(v,'numbase') % make sure at least one base is loaded
    for iBase = 1:viewGet(v, 'numbase');

      % check if anything is being replaced:
      basesform_code = viewGet(v, 'basesformcode', iBase);
      curBaseVol2Tal = viewGet(v, 'baseVol2Tal', iBase);
      curBaseVol2Mag = viewGet(v, 'baseVol2Mag', iBase);
      
      
      if basesform_code == 3 % if was originally aligned to the Vol with old talXform
        mrWarnDlg(sprintf(['(exportTal2mrLR) Base %s has been defined from a base '...
                           ' with a Talairach transformation (e.g., its sform_code is 3),'...
                           ' so changing its TalInfo will create too many problems. We may'...
                           ' be able to support this type of change later, but for now, you'...
                           ' need to stick with the transform youve got.'],viewGet(v,'basename',iBase)));
        % I had actually already written some code to deal with this, which I'll
        % stick in a subfunction below so we can try to implement it later.
        skipBase = 1;
        
      else % Check if base was aligned relative to a different base and needs to be skipped
          
        if and(~isequal(curBaseVol2Mag,vol2mag),~isempty(curBaseVol2Mag))
          sameVol = questdlg(sprintf(['Base %s is defined relative to a different volume. '...
                              'This is almost certainly an error. Please indicate whether to continue '...
                              'anyway, or skip saving the Talairach transformation to this this base (recommended).'],...
                                     viewGet(v,'basename',iBase)),...
                             'Skip this Base?','Skip','Continue','Abort','Skip');
          if any(strcmp(sameVol,'Skip'),strcmp(sameVol,''))
            skipBase = 1;
          elseif strcmp(sameVol,'Continue')
            skipBase = 0;
            disp(sprintf('vol2mag for Base %s has been changed',viewGet(v,'basename',iBase)))              
          elseif strcmp(sameVol,'Abort')
            mrWarnDlg('(exportTal2mrLR) Exiting and not saving the Talairach information to base anatomies.')
            return
          end
        else % if sform_code = 1 and no vol2mag error, don't skip the ROI
          skipBase = 0;
        end % checking for a vol2mag mismatch
      end % checking for problems that mean we should skip the base
      
      if ~skipBase % if we don't need to skip the base, do some more checks and add the talInfo  
        
        % update user on what (if anything)  is changing:
        % first check if there already is a vol2tal
        if ~isequal(curBaseVol2Tal,vol2tal) %check if there's a change
          if isempty(curBaseVol2Tal) % if old vol2tal was empty, then it's easy
            disp(sprintf('vol2tal for Base %s has been added',viewGet(v,'basename',iBase)));
          else  
            disp(sprintf('vol2tal for Base %s has been modified',viewGet(v,'basename',iBase)));
          end
        end
        % also check if adding a new vol2mag
        if isempty(curBaseVol2Mag)
          disp(sprintf('vol2mag for Base %s has been added',viewGet(v,'basename',iBase)));
        end
        
        % finally, make the changes!
        v = viewSet(v,'basevol2tal',vol2tal,iBase); % set the fields
        v = viewSet(v,'basevol2mag',vol2mag,iBase);
        saveAnat(v,iBase,0,0);
        
      end % checking if need to skip the Base
      clear skipBase curBaseVol2Tal curBaseVol2Mag basesform_code    
    end % for iROI = 1:viewGet(v, 'nrois');
  end % make sure at least one base is loaded
end % checking if user wanted to save to bases





%%%%%%%%%%%%%%
function saveTalWithCode3
% I just copied and pasted my old code here so I can use it later if we decide to 
% support changes when sform_code = 3

        curScanVol2Tal = viewGet(v, 'scanVol2Tal', iScan, iGroup);
        curScanVol2Mag = viewGet(v, 'scanVol2Mag', iScan, iGroup);
        
        % first check vol2tal
        if ~isequal(curScanVol2Tal,vol2tal) %check if there's a change
          if isempty(curScanVol2Tal) % if old vol2tal was empty, then it's easy
            disp(sprintf('vol2tal for scan %i in %s group has been added',iScan,groupNames{iGroup}))
          else  % but if replacing a vol2tal, might need to adjust the alignments
            disp(sprintf(['vol2tal for scan %i in %s group has been modified, and sform has been'...
                          ' adjusted if necessary to the new Talairach xform. No need to realign.']...
                         ,iScan, groupNames{iGroup}))
            % need to fix the sform if it was aligned to the vol with a dif tal xform
            sform_code = viewGet(v,'scansformcode');
            if sform_code == 3 % if was originally aligned to the Vol with old talXform
              newSform = vol2tal*inv(curScanVol2Tal)*viewGet(v,'scansform');
              v = viewSet(v,'scansform',newSform,iScan,iGroup); % change in session parameters
                                                                % also need to save it to the nifti header
              filename = viewGet(v,'tseriespath',iScan,iGroup);
              if isfile(filename)
                hdr = mlrImageReadNiftiHeader(filename);
                hdr = cbiSetNiftiSform(hdr,newSform);
                % unfortunately that code sets sform_code to 1 so need to undo that!
                hdr.sform_code = 3;
                hdr = cbiWriteNiftiHeader(hdr,filename);
                hdr = mlrImageReadNiftiHeader(filename); %have to do this again for some reason;
                v = viewSet(v,'niftiHdr',hdr,iScan,iGroup); %also save to session params
              else
                disp(sprintf('(exportTal2Scans) Could not open file %s',filename));
              end % check if file exists
            end % if not previously aligned to Talairach, don't change sform
          end % check if overwriting vol2tal
          v = viewSet(v,'scanvol2tal',vol2tal,iScan,iGroup); % set the field
        end % if they're equal,no need to change anything
        
        % now do vol2mag, easier because won't need to reset s-form
        if ~isequal(curScanVol2Mag,vol2mag)
          if isempty(curScanVol2Mag)
            disp(sprintf('vol2mag for scan %i in %s group has been added',iScan,groupNames{iGroup}))
          else % very weird to overwrite the vol2mag, might be a mistake
            disp(sprintf(['vol2mag for scan %i in %s group is being changed. This is '...
                          'surprising and might be an error.'], iScan, groupNames{iGroup}))
          end
          v = viewSet(v,'scanvol2mag',vol2mag,iScan,iGroup); % set the field
        end % if they're equal,no need to change anything
