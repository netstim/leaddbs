% mlrAnatDBPlugin
%
%        $Id:$ 
%      usage: mlrAnatDBPlugin(action,<v>)
%         by: justin gardner
%       date: 12/28/2014
%    purpose: Plugin function for mercurial based anatomy database
%
function retval = mlrAnatDBPlugin(action,v)

% check arguments
if ~any(nargin == [1 2])
  help DefaultPlugin
  return
end

switch action
 case {'install','i'}
  % check for a valid view
  if (nargin ~= 2) || ~isview(v)
     disp(sprintf('(mlrAnatDBPlugin) Need a valid view to install plugin'));
  else
    % add the Add for mlrAnatDB menu
    mlrAdjustGUI(v,'add','menu','Anat DB','/File/ROI','Callback',@mlrAnatDB,'Separator','on');
    mlrAdjustGUI(v,'add','menu','Anat DB Preferences','/File/Anat DB/','Callback',@mlrAnatDBPreferences);
    mlrAdjustGUI(v,'add','menu','Load ROIs from Anat DB','/File/Anat DB/Anat DB Preferences','Callback',@mlrAnatDBLoadROIs,'Separator','on');
    mlrAdjustGUI(v,'add','menu','Load Base Anatomies from Anat DB','/File/Anat DB/Load ROIs from Anat DB','Callback',@mlrAnatDBLoadBaseAnatomies);
    mlrAdjustGUI(v,'add','menu','Import Surface from Anat DB','/File/Anat DB/Load Base Anatomies from Anat DB','Callback',@mlrAnatDBImportSurface,'Separator','on');
    mlrAdjustGUI(v,'add','menu','Add Session to Anat DB','/File/Anat DB/Import Surface from Anat DB','Callback',@mlrAnatDBAddSession,'Separator','on');
    mlrAdjustGUI(v,'add','menu','Add ROIs to Anat DB','/File/Anat DB/Add Session to Anat DB','Callback',@mlrAnatDBAddROIs);
    mlrAdjustGUI(v,'add','menu','Add Base Anatomies to Anat DB','/File/Anat DB/Add ROIs to Anat DB','Callback',@mlrAnatDBAddBaseAnatomies);
    mlrAdjustGUI(v,'add','menu','Examine ROI in Anat DB','/File/Anat DB/Add Base Anatomies to Anat DB','Callback',@mlrAnatDBExamineROI,'Separator','on');
    mlrAdjustGUI(v,'add','menu','Merge and Check ROIs for Anat DB','/File/Anat DB/Examine ROI in Anat DB','Callback',@mlrAnatDBMergeCheck);
    mlrAdjustGUI(v,'add','menu','Reverse pRF','/File/Anat DB/Merge and Check ROIs for Anat DB','Callback',@mlrAnatDBReversepRF,'Separator','on');
    mlrAdjustGUI(v,'add','menu','Get Reversed pRF','/File/Anat DB/Reverse pRF','Callback',@mlrAnatDBgetReversedpRF);


    % add the callback that will tell mlrAnatDB that a base has been added
    % this is so that we can update the fields that point to the 
    % surfaces from which the anatomy was built (so that flatViewer and
    % makeFlat work)
    v = viewSet(v,'callback','newBase',@mlrAnatDBBaseChange);
    % return true to indicate successful plugin
    retval = true;
   end
 % return a help string
 case {'help','h','?'}
   retval = 'This plugin support exporting sessions and ROIs to a git managed repository';
 otherwise
   disp(sprintf('(mlrAnatDBPlugin) Unknown command %s',action));
end

%%%%%%%%%%%%%%%%%%%
%    mlrAnatDB    %
%%%%%%%%%%%%%%%%%%%
function mlrAnatDB(hObject,eventdata)

% code-snippet to get the view from the hObject variable.
v = viewGet(getfield(guidata(hObject),'viewNum'),'view');

% get repo locations
centralRepo = mrGetPref('mlrAnatDBCentralRepo');
localRepoTop = mrGetPref('mlrAnatDBLocalRepo');

% see if the preference is set
if isempty(centralRepo) || isempty(localRepoTop) 
  % do not enable any thing, because we don't have correct
  % preferences set
  mlrAdjustGUI(v,'set','Load ROIs from Anat DB','Enable','off');
  mlrAdjustGUI(v,'set','Load Base Anatomies from Anat DB','Enable','off');
  mlrAdjustGUI(v,'set','Import Surface from Anat DB','Enable','off');
  mlrAdjustGUI(v,'set','Add Session to Anat DB','Enable','off');
  mlrAdjustGUI(v,'set','Add ROIs to Anat DB','Enable','off');
  mlrAdjustGUI(v,'set','Add Base Anatomies to Anat DB','Enable','off');
  mlrAdjustGUI(v,'set','Examine ROI in Anat DB','Enable','off');
  return
end

% see if we are in an Anat DB session
if ~mlrAnatDBInLocalRepo(v)
  % if not, then only allow add session and examine ROI
  mlrAdjustGUI(v,'set','Load ROIs from Anat DB','Enable','on');
  mlrAdjustGUI(v,'set','Load Base Anatomies from Anat DB','Enable','on');
  mlrAdjustGUI(v,'set','Import Surface from Anat DB','Enable','on');
  mlrAdjustGUI(v,'set','Add Session to Anat DB','Enable','on');
  mlrAdjustGUI(v,'set','Add ROIs to Anat DB','Enable','off');
  mlrAdjustGUI(v,'set','Add Base Anatomies to Anat DB','Enable','on');
  if viewGet(v,'nROIs')
    mlrAdjustGUI(v,'set','Examine ROI in Anat DB','Enable','on');
  else
    mlrAdjustGUI(v,'set','Examine ROI in Anat DB','Enable','off');
  end
else
  % otherwise don't offer add seesion, but add everything else
  % contingent on whether there are ROIs loaded and so forth
  mlrAdjustGUI(v,'set','Load ROIs from Anat DB','Enable','on');
  mlrAdjustGUI(v,'set','Load Base Anatomies from Anat DB','Enable','on');
  mlrAdjustGUI(v,'set','Import Surface from Anat DB','Enable','on');
  mlrAdjustGUI(v,'set','Add Session to Anat DB','Enable','off');

  % see if any bases are loaded
  if viewGet(v,'numBase')
    mlrAdjustGUI(v,'set','Add Base Anatomies to Anat DB','Enable','on');
  else
    mlrAdjustGUI(v,'set','Add Base Anatomies to Anat DB','Enable','off');
  end    
  
  % see if we have any rois loaded, and gray out Add/AnatDB/ROIs menu accordingly
  if viewGet(v,'nROIs')
    mlrAdjustGUI(v,'set','Add ROIs to Anat DB','Enable','on');
    mlrAdjustGUI(v,'set','Examine ROI in Anat DB','Enable','on');
  else
    mlrAdjustGUI(v,'set','Add ROIs to Anat DB','Enable','off');
    mlrAdjustGUI(v,'set','Examine ROI in Anat DB','Enable','off');
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    mlrAnatDBAddSession    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mlrAnatDBAddSession(hObject,eventdata)

% code-snippet to get the view from the hObject variable.
v = viewGet(getfield(guidata(hObject),'viewNum'),'view');

% remember what directory we started in
curpwd = pwd;

% Warn user that we are about to close session and will need to open it up again
% in the mlrAnatDB
if strcmp(questdlg(sprintf('(mlrAnatDBPlugin) Will now copy (using hard links) your current session into mlrAnatDB. To do so, will need to temporarily close the current session and then reopen in the mlrAnatDB session. Your current work will be saved as usual through the mrLastView mechanism which stores all your current settings. Also, this will not take any more hard disk space, since the files will be copied as hard links. Click OK to continue, or cancel to cancel this operation. If you hit cancel, you will be able to run File/Anat DB/Add Session at a later time.'),'mlrAnatDBPlugin','Ok','Cancel','Cancel'),'Cancel')
  cd(curpwd);
  return
end
  
% get subjectID
subjectID = mlrAnatDBSubjectID(v);

% get home directory
homeDir = viewGet(v,'homeDir');

% user said we could close, so do it
mrQuit;

% put the session into the repo
if ~mlrAnatDBPut(subjectID,homeDir,'localizer')
  % failure, so go back to the other directory
  cd(homeDir);
  mrLoadRet;
  cd(curpwd);
  return
end

% everything went ok, switch directories and start up over there
[localRepo localRepoLargeFiles] = mlrAnatDBGetRepo(subjectID);
cd(fullfile(localRepoLargeFiles,'localizers',getLastDir(homeDir)));
mrLoadRet;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    mlrAnatDBAddROIs    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mlrAnatDBAddROIs(hObject,eventdata)

% code-snippet to get the view from the hObject variable.
v = viewGet(getfield(guidata(hObject),'viewNum'),'view');

% put base anatomies into repo
mlrAnatDBPut(v,v,'rois');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    mlrAnatDBAddBaseAnatomies    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mlrAnatDBAddBaseAnatomies(hObject,eventdata)

% code-snippet to get the view from the hObject variable.
v = viewGet(getfield(guidata(hObject),'viewNum'),'view');

% put base anatomies into repo
mlrAnatDBPut(v,v,'mlrBaseAnat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    mlrAnatDBgetReversedpRF    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mlrAnatDBgetReversedpRF(hObject,eventdata)
%%
v = viewGet(getfield(guidata(hObject),'viewNum'),'view');

curpwd = pwd;

% get subjectID
subjectID = mlrAnatDBSubjectID(v);

% pull repos, note we need the large file repo (unfortunately), because
% otherwise we don't have access to the pRF analysis
[~, localRepoLargeFiles] = mlrAnatDBGetRepo(subjectID);

% try to find the pRF analysis automatically
pRFloc = '';
locdir = fullfile(localRepoLargeFiles,'localizers',sprintf('%s*',subjectID));
localizers = dir(locdir);
if length(localizers)>1
    warning('Implement choosing between multiple localizers');
    keyboard
elseif isempty(localizers)
    warning('Did not find any localizers for this subject.');
    return
else
    pRFloc = fullfile(localRepoLargeFiles,'localizers',localizers(1).name);
end

if strcmp(questdlg(sprintf('(mlrAnatDBPlugin) Before you run this make sure you selected the analysis and scan you want the new overlays to be added to!!!'),'mlrAnatDBPlugin','Ok','Cancel','Cancel'),'Cancel')
  return
end
if strcmp(questdlg(sprintf('(mlrAnatDBPlugin) I need to close your current session while I load the overlays from the pRF session. Is that okay?'),'mlrAnatDBPlugin','Ok','Cancel','Cancel'),'Cancel')
  return
end

mrQuit;

cd(pRFloc);

v = newView();
if strcmp(questdlg(sprintf('(mlrAnatDBPlugin) Please pick the analysis you want to load overlays from.'),'mlrAnatDBPlugin','Ok','Cancel','Cancel'),'Cancel')
    mrQuit
    mrQuit(0);
    cd(curpwd)
    mrLoadRet;
    return
end
v = loadAnalysis(v);
overlays = v.analyses{1}.overlays;

opts = [];
chosen = [];
while isempty(chosen)
    disp(sprintf('\nAll overlays available'));
    for i = 1:length(overlays)
        cOver = overlays(i);
        disp(sprintf('Overlay %i: %s',i,cOver.name));
        if strfind(cOver.name,'Overlap')
            chosen = [chosen i];
        end
        opts = [opts i];
    end
    disp(sprintf('\nOverlays currently selected'));
    for ci = chosen
        disp(sprintf('Overlay %i: %s',ci,overlays(ci).name));
    end
    nchosen = input('Input a new array (e.g. [1 3]) or hit enter to continue: ');
    if ~isempty(nchosen)
        chosen = nchosen;
    end
end

scan2magFrom = viewGet(v,'scan2mag');

cOverlays = overlays(chosen);

% reopen the original mlr directory
mrQuit(0);
cd(curpwd)
v = mrLoadRet;

% transform cOverlays into the current scan space
scan2magTo = viewGet(v,'scan2mag');

scan2scan = inv(scan2magFrom) * scan2magTo;
analysis = viewGet(v,'analysis');
toDims = viewGet(v,'scandims');
[x,y,z] = ind2sub(toDims,1:(toDims(1)*toDims(2)*toDims(3)));
out = scan2scan*[x;y;z;ones(size(x))];

for ci = 1:length(cOverlays)
    if length(cOverlays(ci).data)>1
        %
        disp('(mlrAnatDBPlugin) I don''t know which scan to pull from. Please choose...');
        warning('NOT IMPLEMENTED');
        keyboard
    else
        dat = cOverlays(ci).data{1};
    end
    lOverlay = interp3(dat,out(2,:),out(1,:),out(3,:),'nearest');
    lOverlaydat = reshape(lOverlay,toDims);
    [v,analysis] = mrDispOverlay(lOverlaydat,1,analysis,v,sprintf('overlayName=%s',cOverlays(ci).name));
    refreshMLRDisplay(viewGet(v,'viewNum'));
    mglWaitSecs(1); % I don't know what's wrong here but sometimes this code doesn't work, this might fix it??
end

saveAnalysis(v,analysis.name);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    mlrAnatDBReversepRF    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mlrAnatDBReversepRF(hObject,eventdata)

v = viewGet(getfield(guidata(hObject),'viewNum'),'view');

% find pRF analysis
opts = [];
chosen = -1;
while chosen == -1 || isempty(chosen)
    for i = 1:length(v.analyses)
        cAnal = v.analyses{i};
        disp(sprintf('Analysis %i: %s',i,cAnal.name));
        if strfind(cAnal.name,'pRF')
            chosen = i;
        end
        opts = [opts i];
    end
    chosen = input('Which analysis do you want to pick? Choose one by typing in the number: ');
end

analysis = v.analyses{i};

% define stimulus range
polar = [];
while isempty(polar)
    resp = input('Is your input in polar [p] or cartesian [c] coordinates? ','s');
    if strcmpi(resp,'p')
        polar = 1;
        warning('Polar angle coordinates are not implemented');
        return
    elseif strcmpi(resp,'c')
        polar = 0;
    end
end

knownAnalyses = {'cohcon'};
knownCoords = {{[3.5 12 -7 7],[-12 -3.5 -7 7]}};
knownNames = {{'right','left'}};
knownStr = knownAnalyses{1};
for ki = 2:length(knownAnalyses)
    knownStr = [knownStr,', ',sprintf('''%s''',knownAnalyses{ki})];
end

anName = '';
names = '';
if ~polar
    resp = [];
    while ~iscell(resp)
        resp = input(sprintf('Input your stimulus as a cell.\nEach cell is an array with x0,x1,y0,y1 defining a rectangular space.\nYou can also input %s.\n',knownStr));
        
        for ki = 1:length(knownAnalyses)
            if strcmp(resp,knownAnalyses{ki})
                disp(sprintf('You requested the %s stimulus.',knownAnalyses{ki}));
                resp = knownCoords{ki};
                anName = knownAnalyses{ki};
                names = knownNames{ki};
            end
        end
        if any(cellfun(@(x) length(x)~=4,resp))
            disp('Badly conditioned responses, each array should be length 4');
            resp = [];
        end
    end
    stimulusCoords = resp;
else
    warning('Implement me!');
end

if isempty(anName)
    anName = input('Please input the name of your experiment (e.g. cohcon)\n','s');
end

dims = viewGet(v,'scandims');
total = dims(1)*dims(2)*dims(3);

% figureo out where in analysis.overlays the data is stored
whichDim = 0;
for i = 1:length(analysis.overlays(1).data)
    if isequal(dims,size(analysis.overlays(1).data{i}))
        whichDim = i;
    end
end

data.r2 = analysis.overlays(1).data{whichDim};
data.pa = analysis.overlays(2).data{whichDim};
data.ecc = analysis.overlays(3).data{whichDim};
data.std = analysis.overlays(4).data{whichDim};

% f = figure;
stimulusOverlap = zeros(dims(1),dims(2),dims(3),length(stimulusCoords));
disppercent(-inf,sprintf('Computing %i voxels...',total));
count = 0;
for x = 1:dims(1)
    for y = 1:dims(2)
        for z = 1:dims(3)
            stimulusOverlap(x,y,z,:) = computeOverlap(data,stimulusCoords,x,y,z,0);
            count = count+1;
            disppercent(count/total);
        end
    end
end
disppercent(inf);

if isempty(names)
    warning('implement me!');
    keyboard
end

[v, analysis] = mrDispOverlay(stimulusOverlap(:,:,:,1),viewGet(v,'curscan'),analysis,v,sprintf('overlayName=%sOverlap%s',names{1},anName));
[v, analysis] = mrDispOverlay(stimulusOverlap(:,:,:,2),viewGet(v,'curscan'),analysis,v,sprintf('overlayName=%sOverlap%s',names{2},anName));

refreshMLRDisplay(viewGet(v,'viewNum'));

saveAnalysis(v,analysis.name);

function overlap = computeOverlap(data,stimulusCoords,x,y,z,f)

% helper function for mlrAnatDBReversepRF

pa = data.pa(x,y,z);
ecc = data.ecc(x,y,z);
std = data.std(x,y,z);

if isnan(pa) || isnan(ecc) || (std <= 0)
    overlap = nan;
    return
else
    [cx,cy] = pol2cart(pa,ecc);
end

p = zeros(1,length(stimulusCoords));
for i = 1:length(stimulusCoords)
    sc = stimulusCoords{i};
    p(i) = mvncdf([sc(1) sc(3)],[sc(2) sc(4)],[cx,cy],[std 0;0 std]);
end

%%
if f>0
    figure(f)
    clf
    [sfx, sfy] = meshgrid(-20:20,-20:20);

    % out = gauss([1 cx cy std std 0 0],mfx,mfy);
    outs = gauss([1 cx cy std std 0 0],sfx,sfy);

    figure
    imagesc(sfx(1,:),sfy(:,1),outs);
    hold on
    rectangle('Position',[stimulusCoords{1}(1) stimulusCoords{1}(3) stimulusCoords{1}(2)-stimulusCoords{1}(1) stimulusCoords{1}(4)-stimulusCoords{1}(3)]);
    rectangle('Position',[stimulusCoords{2}(1) stimulusCoords{2}(3) stimulusCoords{2}(2)-stimulusCoords{2}(1) stimulusCoords{1}(4)-stimulusCoords{1}(3)]);

    title(sprintf('Left %0.3f Right %0.3f',p(2),p(1)));
end

overlap = p;

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    mlrAnatDBMergeCheck    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mlrAnatDBMergeCheck(hObject,eventdata)

% code-snippet to get the view from the hObject variable. 
v = viewGet(getfield(guidata(hObject),'viewNum'),'view');

disp(sprintf('\n\n\n\n\n\n\n\n(mlrAnatDBMergeCheck) Searching for ROIs that fit standards'));
pfxs = {'l','r'};
standards = {'V1','V2','V3','V4','V3a','V3b','V7','LO1','LO2','MT'};
% alternates = {'hV4','hMT+'};

% Check for all l/r ROIs and ask to rename
roiNames = viewGet(v,'roiNames');
if isempty(roiNames)
    disp('(mlrAnatDBMergeCheck) No ROIs to check'); return
end

nofind = {};
for pi = 1:length(pfxs)
    for si = 1:length(standards)
        searchfor = sprintf('%s%s',pfxs{pi},standards{si});
        [found, idx] = checkROIs(roiNames,searchfor,pfxs{pi},standards{si});
        if found==0
            disp(sprintf('(mlrAnatDBMergeCheck) You have no ROI: %s',searchfor));
            nofind{end+1} = searchfor;
        elseif found==2
            disp(sprintf('(mlrAnatDBMergeCheck) ROI %s appears to be mislabeled as %s.',searchfor,roiNames{idx}));
                        nofind{end+1} = searchfor;

        elseif found==3               
            disp(sprintf('(mlrAnatDBMergeCheck) You have no ROI: %s, you have %s which is similar...',searchfor,roiNames{idx}));
                        nofind{end+1} = searchfor;

        end
    end
end

if ~isempty(nofind)
    disp(sprintf('\n(mlrAnatDBMergeCheck) Please define the missing ROIs.\n\t\t\tYour ROIs may be mis-named.'));
    return
end
disp(sprintf('(mlrAnatDBMergeCheck) Found all standards'));


% Merge lV1+rV1 into V1, etc...
for si = 1:length(standards)
    cur = standards{si};
    curL = sprintf('%s%s',pfxs{1},standards{si});
    curR = sprintf('%s%s',pfxs{2},standards{si});
    [found, ~] = checkROIs(roiNames,cur,'','');
    [foundL, ~] = checkROIs(roiNames,curL,'','');
    [foundR, ~] = checkROIs(roiNames,curR,'','');
    if ~(found==1)
        if foundL&&foundR
            disp(sprintf('(mlrAnatDBMergeCheck) No ROI %s found. Computing the union of %s and %s',cur,curL,curR));
            v = combineROIs(v,curL,curR,'Union',cur);
        else
            disp('Something went wrong...');
            keyboard
        end
    else
        disp(sprintf('(mlrAnatDBMergeCheck) ROI %s already exists, skipping',cur));
    end
end
%%%%%
% HELPER FUNCTION TO CHECK FOR ROIS
%%%%%

function [found, idx] = checkROIs(rois,searchfor,prefix,suffix)
found = 0;
idx = -1;
for ri = 1:length(rois)
    if strcmp(rois{ri},searchfor)
        idx = ri; found = 1; return
    end
end
for ri = 1:length(rois)
    rois{ri} = lower(rois{ri});
end
suffix = lower(suffix);
% didn't find, try to find one that is similar?
for ri = 1:length(rois)
    if strcmp(rois{ri},searchfor)
        found = 3; idx = ri;
        return
    end
    if ~isempty(strfind(rois{ri},suffix))
        idx = ri; found = 2; 
        if isempty(strfind(rois{ri},prefix)) || length(rois{ri})~=(length(prefix)+length(suffix))
            found = 3;
        end
        return
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    mlrAnatDBExamineROI    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mlrAnatDBExamineROI(hObject,eventdata)

% code-snippet to get the view from the hObject variable. 
v = viewGet(getfield(guidata(hObject),'viewNum'),'view');

% select an ROI to examine
roiNames = viewGet(v,'roiNames');
curROI = viewGet(v,'curROI');
if (curROI >1) && (length(roiNames)>=curROI)
  roiNames = putOnTopOfList(roiNames{curROI},roiNames);
end
paramsInfo = {{'roiToExamine',roiNames,'Choose which ROI to examine in its original localizer session'}};
params = mrParamsDialog(paramsInfo,'Choose ROI');
if isempty(params),return,end

% get the original session that the ROI was defined in and make
% sure we have it in the repo
roi = viewGet(v,'roi',params.roiToExamine);

% check for createdFromSession
if isempty(roi.createdFromSession)
  mrWarnDlg(sprintf('(mlrAnatDBPlugin:mlrAnatDBExamineROI) The ROI %s does not have the field createdFromSession set. Not sure which session it was created from',roi.name));
  return
end

% Check that we have mercurial installed correctly
if ~mlrAnatDBCheckHg,return,end

% get subject ID
subjectID = mlrAnatDBSubjectID(v);
if isempty(subjectID),return,end

% get the local repos
[localRepoSubject localRepoSubjectLargeFiles] = mlrAnatDBGetRepo(subjectID);
if isempty(localRepoSubject) || isempty(localRepoSubjectLargeFiles),return,end

% check to see if session exists in repo
createdFromSession = fullfile(localRepoSubjectLargeFiles,roi.createdFromSession);
if ~isdir(createdFromSession)
  mrWarnDlg(sprintf('(mlrAnatDBPlugin:mlrAnatDBExamineROI) Could not find session %s in mlrAnatDB',createdFromSession));
  return
end
  
% now confirm that this is what the user really wants to do
if strcmp(questdlg(sprintf('(mlrAnatDBPlugin:mlrAnatDBExamineROI) Will now close this current session (saving work as always) and will load up the localizer session where ROI: %s was defined',params.roiToExamine),'Switch to Localizer session','Ok','Cancel','Cancel'),'Cancel')
  return
end

% ok, user said we could close, so do it
mrQuit;
cd(createdFromSession);
mrLoadRet
v = getMLRView;

% check to see if ROI is loaded
if ~any(strcmp(viewGet(v,'roiNames'),params.roiToExamine))
  % then load it
  v = loadROI(v,params.roiToExamine,0,fullfile(localRepoSubject,'mlrROIs'));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    mlrAnatDBLoadROIs    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mlrAnatDBLoadROIs(hObject,eventdata)

% code-snippet to get the view from the hObject variable.
v = viewGet(getfield(guidata(hObject),'viewNum'),'view');

% Check that we have mercurial installed correctly
if ~mlrAnatDBCheckHg,return,end

% get subject ID
subjectID = mlrAnatDBSubjectID(v);
if isempty(subjectID),return,end

% get the local repos
localRepoSubject = mlrAnatDBGetRepo(subjectID);
if isempty(localRepoSubject),return,end

% load the rois
v = loadROI(v,[],[],fullfile(localRepoSubject,'mlrROIs'));

% and refresh
refreshMLRDisplay(v);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    mlrAnatDBLoadBaseAnatomies    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mlrAnatDBLoadBaseAnatomies(hObject,eventdata)

% code-snippet to get the view from the hObject variable. 
v = viewGet(getfield(guidata(hObject),'viewNum'),'view');

% Check that we have mercurial installed correctly
if ~mlrAnatDBCheckHg,return,end

% get subject ID
subjectID = mlrAnatDBSubjectID(v);
if isempty(subjectID),return,end

% get the local repos
localRepoSubject = mlrAnatDBGetRepo(subjectID);
if isempty(localRepoSubject),return,end

% load the rois
v = loadAnat(v,[],fullfile(localRepoSubject,'mlrBaseAnatomies'));

% and refresh
refreshMLRDisplay(v);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    mlrAnatDBImportSurface    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mlrAnatDBImportSurface(hObject,eventdata)

% code-snippet to get the view from the hObject variable. 
v = viewGet(getfield(guidata(hObject),'viewNum'),'view');

% Check that we have mercurial installed correctly
if ~mlrAnatDBCheckHg,return,end

% get subject ID
subjectID = mlrAnatDBSubjectID(v);
if isempty(subjectID),return,end

% get the local repos
localRepoSubject = mlrAnatDBGetRepo(subjectID);
if isempty(localRepoSubject),return,end

% get the surface
filterspec = {'*.off','Off Surface file (*.off)'};
title = 'Choose outer surface file';
pathStr = mlrGetPathStrDialog(fullfile(localRepoSubject,'surfaces'),title,filterspec,'off');

% open the path up
if ~isempty(pathStr)
  base = importSurfaceOFF(pathStr);
  if ~isempty(base)
    v = viewSet(v, 'newbase', base);
    refreshMLRDisplay(v);
  end
end



%%%%%%%%%%%%%%%%%%
%    mysystem    %
%%%%%%%%%%%%%%%%%%
function [status,result] = mysystem(command)

disp(sprintf('(mlrAnatDBPlugin): %s',command));
[status,result] = system(command,'-echo');


%%%%%%%%%%%%%%%%%%%%
%    baseChange    %
%%%%%%%%%%%%%%%%%%%%
function v = mlrAnatDBBaseChange(v)

% get the base coord map for the new base
baseCoordMap = viewGet(v,'baseCoordMap');
% if not empty, we will check to see if its file
% pointers point to mlrAnatDB correctly
if ~isempty(baseCoordMap)
  % check path
  subjectID = getLastDir(fileparts(baseCoordMap.path));
  % update the path if it was already pointing to a subjectID directory
  if ~isempty(subjectID) && (subjectID(1) == 's')
    % change the path to point to the local mlrAnatDB
    baseCoordMap.path = mlrReplaceTilde(fullfile(mrGetPref('mlrAnatDBLocalRepo'),subjectID,getLastDir(baseCoordMap.path)));
    % and update
    v = viewSet(v,'baseCoordMap',baseCoordMap);
    % and see if the repo exists
    if isempty(mlrAnatDBGetRepo(subjectID,'noPull=1'))
      mlrAnatDBGetRepo(subjectID);
    end
  end
end
