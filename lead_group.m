function varargout = lead_group(varargin)
% LEAD_GROUP MATLAB code for lead_group.fig
%      LEAD_GROUP, by itself, creates a new LEAD_GROUP or raises the existing
%      singleton*.
%
%      H = LEAD_GROUP returns the handle to a new LEAD_GROUP or the handle to
%      the existing singleton*.
%
%      LEAD_GROUP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LEAD_GROUP.M with the given input arguments.
%
%      LEAD_GROUP('Property','Value',...) creates a new LEAD_GROUP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before lead_group_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to lead_group_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help lead_group

% Last Modified by GUIDE v2.5 23-Feb-2023 11:08:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @lead_group_OpeningFcn, ...
    'gui_OutputFcn',  @lead_group_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before lead_group is made visible.
function lead_group_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to lead_group (see VARARGIN)

handles.prod='group';
handles.callingfunction='lead_group';

% add recentgroups groups...
ea_initrecent(handles, 'groups');

% Choose default command line output for lead_group
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes lead_group wait for user response (see UIRESUME)
% uiwait(handles.leadfigure);

options.earoot = ea_getearoot;
options.prefs = ea_prefs;
setappdata(handles.leadfigure,'earoot',options.earoot);

% Build popup tables:

% atlassets:
atlases = dir(ea_space(options,'atlases'));
atlases = {atlases(cell2mat({atlases.isdir})).name};
atlases = atlases(cellfun(@(x) ~strcmp(x(1),'.'), atlases));
atlases{end+1} = 'Use none';

set(handles.atlassetpopup,'String', atlases);
[~, defix]=ismember(options.prefs.machine.defaultatlas, atlases);
if defix
    set(handles.atlassetpopup,'Value',defix);
end

% setup vat functions
funcs = ea_regexpdir(ea_getearoot, 'ea_genvat_.*\.m$', 0);
funcs = regexp(funcs, '(ea_genvat_.*)(?=\.m)', 'match', 'once');
names = cellfun(@(x) eval([x, '(''prompt'');']), funcs, 'Uni', 0);

if ~options.prefs.env.dev
    ossdbsInd = find(contains(names, 'OSS-DBS'));
    funcs(ossdbsInd) = [];
    names(ossdbsInd) = [];
end

setappdata(handles.leadfigure,'genvatfunctions',funcs);
setappdata(handles.leadfigure,'vatfunctionnames',names);

% set background image
set(gcf,'color','w');
im=imread([options.earoot,'icons',filesep,'logo_lead_group.png']);
image(im);
axis off;
axis equal;

try
    priorselection=find(ismember(fiberscell,stimparams.usefiberset)); % retrieve prior selection of fiberset.
    set(handles.fiberspopup,'Value',priorselection);
catch    % reinitialize using third entry.
    set(handles.fiberspopup,'Value',1);
end

if get(handles.fiberspopup,'Value')>length(get(handles.fiberspopup,'String'))
    set(handles.fiberspopup,'Value',length(get(handles.fiberspopup,'String')));
end

% Labels:
labeling = dir([ea_space(options,'labeling'),'*.nii']);
if ~isempty(labeling)
    labeling = cellfun(@(x) {strrep(x, '.nii', '')}, {labeling.name});
    set(handles.labelpopup,'String', labeling);

    % Initialize parcellation popupmenu
    parc = find(ismember(labeling, options.prefs.lg.defaultParcellation));
    if ~isempty(parc)
        set(handles.labelpopup,'Value',parc);
    end
end

% Set connectome popup
modlist = ea_genmodlist([],[],options,'dmri');
modlist{end+1}='Do not calculate connectivity stats';
set(handles.fiberspopup,'String',modlist);
set(handles.fiberspopup,'Value',length(modlist));

% set version text:
set(handles.versiontxt,'String',['v',ea_getvsn('local')]);

% make listboxes multiselectable:
set(handles.patientlist,'Max',100,'Min',0);
set(handles.grouplist,'Max',100,'Min',0);
set(handles.clinicallist,'Max',100,'Min',0);

if options.prefs.env.dev
    disp('Running in Developer Mode...')
    set(handles.mercheck,'Visible','on')
end

if ~isempty(varargin) && isfile(GetFullPath(varargin{1})) % Path to group analysis file provided as input
    groupFilePath = GetFullPath(varargin{1});
    load(groupFilePath, 'M');
    M.root = [fileparts(groupFilePath), filesep];
    set(handles.groupdir_choosebox,'String',M.root);
    set(handles.groupdir_choosebox,'TooltipString', M.root);
    setappdata(handles.leadfigure, 'M', M);
    try
        setappdata(handles.leadfigure, 'S', M.S);
        setappdata(handles.leadfigure, 'vatmodel', M.S(1).model);
    end
    ea_addrecent(handles,{M.root},'groups');
else
    M=getappdata(gcf,'M');
    if isempty(M)
        % initialize Model variable M
        M=ea_initializeM;
    end
end

setappdata(gcf,'M',M);
setappdata(handles.leadfigure, 'options', options);

ea_refresh_lg(handles);

ea_firstrun(handles,options);

ea_menu_initmenu(handles,{'prefs','transfer','group'},options.prefs);

ea_processguiargs(handles,varargin)

ea_bind_dragndrop(handles.leadfigure, ...
    @(obj,evt) DropFcn(obj,evt,handles), ...
    @(obj,evt) DropFcn(obj,evt,handles));

ea_ListBoxRenderer(handles.recentgroups);
ea_ListBoxRenderer(handles.atlassetpopup);
ea_ListBoxRenderer(handles.labelpopup);
ea_ListBoxRenderer(handles.fiberspopup);


% --- Drag and drop callback to load patdirs.
function DropFcn(~, event, handles)

% check if dropping area is in patient listbox
if event.Location.getX < 325 && event.Location.getX > 24 && ...
   event.Location.getY < 322 && event.Location.getY > 137
    target = 'patientList';
else
    target = 'groupDir';
end

switch event.DropType
    case 'file'
        folders = event.Data;
    case 'string'
        folders = {event.Data};
end

if strcmp(target, 'groupDir')
    % Save data for previous selected group folder
    if ~strcmp(handles.groupdir_choosebox.String,'Choose Dataset Directory') % group dir still not chosen
        ea_busyaction('on',handles.leadfigure,'group');
        disp('Saving data...');
        % save M
        ea_refresh_lg(handles);
        M = getappdata(handles.leadfigure,'M');

        disp('Saving data to disk...');
        try
            save(ea_getGroupAnalysisFile(handles.groupdir_choosebox.String),'M','-v7.3');
        catch
            warning('Data could not be saved.');
            keyboard
        end
        disp('Done.');
        ea_busyaction('off',handles.leadfigure,'group');
    end

    % Multiple folder dragged or not a proper BIDS folder
    if length(folders) > 1
        ea_error('Please drag either a dataset root folder or a group analysis folder into Lead Group!', simpleStack = 1);
    end

    if isfile(folders{1}) % Group analysis mat file dragged
        if ~isempty(regexp(folders{1}, ['derivatives\', filesep, 'leadgroup\', filesep, '.+\', filesep, 'dataset-.+_analysis-.+\.mat$'], 'match', 'once'))
            % Group analysis file within dataset folder
            groupdir = [fileparts(folders{1}), filesep];
            load(folders{1}, 'M');

            datasetFolder = regexp(groupdir, ['(.*)(?=\', filesep, 'derivatives\', filesep, 'leadgroup)'], 'match', 'once');
            if isfile(fullfile(datasetFolder, 'miniset.json'))
                for p = 1:size(M.patient.list,1)
                    [~, patient_tag] = fileparts(M.patient.list{p});
                    M.patient.list{p} = fullfile(datasetFolder, 'derivatives', 'leaddbs', patient_tag);
                end
                M.root = groupdir;
                save(folders{1}, 'M')
            end
        elseif ~isempty(regexp(folders{1}, ['\', filesep, 'dataset-.+_analysis-.+\.mat$'], 'match', 'once'))
            % Orphan group analysis file, will create proper dataset folder
            [groupdir, analysisFile] = ea_genDatasetFromGroupAnalysis(folders{1});
            load(analysisFile, 'M');
        else
            ea_error('Not a Lead Group Analysis file!', simpleStack = 1);
        end
    else % Dataset root folder or group analysis folder dragged
        if ~contains(folders{1}, ['derivatives', filesep, 'leadgroup', filesep]) && ~isfolder(fullfile(folders{1}, 'derivatives'))
            analysisFile = ea_regexpdir(folders{1}, '^dataset-[^\W_]+_analysis-[^\W_]+\.mat$', 0);
            if ~isempty(analysisFile)
               folders{1} = ea_genDatasetFromGroupAnalysis(analysisFile{1});
            end
        end
        analysisFile = ea_getGroupAnalysisFile(folders{1});
        if isempty(analysisFile) % Create new analysis file in case not found
            analysisFile = ea_genGroupAnalysisFile(folders{1});
        end
        groupdir = [fileparts(analysisFile), filesep];
        load(analysisFile, 'M');

        datasetFolder = regexp(groupdir, ['(.*)(?=\', filesep, 'derivatives\', filesep, 'leadgroup)'], 'match', 'once');
        if isfile(fullfile(datasetFolder, 'miniset.json'))
            for p = 1:size(M.patient.list,1)
                [~, patient_tag] = fileparts(M.patient.list{p});
                M.patient.list{p} = fullfile(datasetFolder, 'derivatives', 'leaddbs', patient_tag);
            end
            M.root = groupdir;
            save(analysisFile, 'M')
        end

    end

    set(handles.groupdir_choosebox, 'String', groupdir);
    set(handles.groupdir_choosebox, 'TooltipString', groupdir);

    ea_busyaction('on',handles.leadfigure,'group');

    setappdata(handles.leadfigure,'M',M);

    try
        setappdata(handles.leadfigure,'S',M.S);
        setappdata(handles.leadfigure,'vatmodel',M.S(1).model);
    end

    ea_busyaction('off',handles.leadfigure,'group');
    ea_refresh_lg(handles);
else
    if strcmp(handles.groupdir_choosebox.String, 'Choose Dataset Directory')
        ea_error('Please choose a group directory first to store the group analysis!', simpleStack = 1);
    end

    nonexist = cellfun(@(x) ~exist(x, 'dir'), folders);
    if any(nonexist)
        fprintf('\nExcluded non-existent/invalid folder:\n');
        cellfun(@disp, folders(nonexist));
        fprintf('\n');
        folders(nonexist) = [];
    end

    if ~isempty(folders)
        M=getappdata(handles.leadfigure,'M');

        M.patient.list=[M.patient.list; folders];
        M.patient.group=[M.patient.group; ones(length(folders),1)];
        options=ea_setopts_local(handles);

        tS=ea_initializeS(['gs_',M.guid],options,handles);
        if isempty(M.S)
            M=rmfield(M,'S');
            M.S(1:length(folders))=tS;
        else
            M.S(end+1:end+length(folders))=tS;
        end
        setappdata(handles.leadfigure, 'M', M);
        ea_refresh_lg(handles);
        % save M
        M=getappdata(handles.leadfigure,'M');
        save(ea_getGroupAnalysisFile(handles.groupdir_choosebox.String),'M','-v7.3');
    end
end


% --- Outputs from this function are returned to the command line.
function varargout = lead_group_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in groupdir_choosebox.
function groupdir_choosebox_Callback(hObject, eventdata, handles)
% hObject    handle to groupdir_choosebox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Save data for previous selected group folder
if ~strcmp(handles.groupdir_choosebox.String,'Choose Dataset Directory') % group dir still not chosen
    ea_busyaction('on',handles.leadfigure,'group');
    disp('Saving data...');
    % save M
    ea_refresh_lg(handles);
    M=getappdata(handles.leadfigure,'M');
    disp('Saving data to disk...');
    try
        save(ea_getGroupAnalysisFile(handles.groupdir_choosebox.String),'M','-v7.3');
    catch
        warning('Data could not be saved.');
        % keyboard
    end
    disp('Done.');
    ea_busyaction('off',handles.leadfigure,'group');
end

% groupdir=ea_uigetdir(ea_startpath,'Choose Dataset Directory');
groupdir = uigetdir;

if ~groupdir % user pressed cancel
    return
else
    if ~contains(groupdir, ['derivatives', filesep, 'leadgroup', filesep]) && ~isfolder(fullfile(groupdir, 'derivatives'))
        analysisFile = ea_regexpdir(groupdir, '^dataset-[^\W_]+_analysis-[^\W_]+\.mat$', 0);
        if ~isempty(analysisFile) % Load group analysis outside of dataset
            if length(analysisFile) > 1
                [~, guid] = fileparts(fileparts(analysisFile));
                index = listdlg('PromptString', 'Select Group Analysis', 'ListString', guid, 'SelectionMode', 'single', 'CancelString', 'Cancel');
                if ~isempty(index)
                    analysisFile = analysisFile(index);
                else
                    return;
                end
            end
            groupdir = ea_genDatasetFromGroupAnalysis(analysisFile{1});
        else % initialize group analysis in an empty folder.
            ea_cprintf('CmdWinWarnings', 'Initialize new dataset folder: %s\n', groupdir);
            ea_mkdir(fullfile(groupdir, 'derivatives', 'leaddbs'));
            ea_mkdir(fullfile(groupdir, 'derivatives', 'leadgroup'));
        end
    end
    analysisFile = ea_getGroupAnalysisFile(groupdir);
    if isempty(analysisFile) % Create new analysis file in case not found
        analysisFile = ea_genGroupAnalysisFile(groupdir);
    end
    groupdir = fileparts(analysisFile);
end

ea_load_group(handles,groupdir);


% --- Executes on selection change in patientlist.
function patientlist_Callback(hObject, eventdata, handles)
% hObject    handle to patientlist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns patientlist contents as cell array
%        contents{get(hObject,'Value')} returns selected item from patientlist
M=getappdata(gcf,'M');

M.ui.listselect=get(handles.patientlist,'Value');
set(handles.grouplist,'Value',M.ui.listselect);

setappdata(gcf,'M',M);
ea_refresh_lg(handles);


% --- Executes during object creation, after setting all properties.
function patientlist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to patientlist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in addptbutton.
function addptbutton_Callback(hObject, eventdata, handles)
% hObject    handle to addptbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(handles.groupdir_choosebox.String, 'Choose Dataset Directory')
    ea_error('Please choose a group directory first to store the group analysis!', simpleStack = 1);
end

M=getappdata(handles.leadfigure,'M');

folders=ea_uigetdir(ea_startpath,'Select Patient folders..');
M.patient.list=[M.patient.list;folders'];
M.patient.group=[M.patient.group;ones(length(folders),1)];
options=ea_setopts_local(handles);

tS=ea_initializeS(['gs_',M.guid],options,handles);

if isempty(M.S)
    M=rmfield(M,'S');
    M.S(1:length(folders))=tS;
else
    try
        M.S(end+1:end+length(folders))=tS;
    catch
        tS.volume=[0,0];
        tS.sources=[1:4];
        M.S(end+1:end+length(folders))=tS;
    end
end

setappdata(handles.leadfigure,'M',M);
setappdata(handles.leadfigure,'S',M.S);
ea_refresh_lg(handles);
% save M
M=getappdata(handles.leadfigure,'M');
save(ea_getGroupAnalysisFile(handles.groupdir_choosebox.String),'M','-v7.3');


% --- Executes on button press in removeptbutton.
function removeptbutton_Callback(hObject, eventdata, handles)
% hObject    handle to removeptbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
M=getappdata(handles.leadfigure,'M');

deleteentry=get(handles.patientlist,'Value');

M.patient.list(deleteentry)=[];

M.patient.group(deleteentry)=[];

try M.elstruct(deleteentry)=[]; end

for cvar=1:length(M.clinical.vars)
    try
        M.clinical.vars{cvar}(deleteentry,:)=[];
    end
end

if isfield(M,'S')
    try
    M.S(deleteentry)=[];
    end
    setappdata(handles.leadfigure, 'S', M.S);
end

try
    M.stats(deleteentry)=[];
end
setappdata(handles.leadfigure,'M',M);
ea_refresh_lg(handles);


% --- Executes on button press in vizbutton.
function vizbutton_Callback(hObject, eventdata, handles)
% hObject    handle to vizbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc;
M=getappdata(gcf,'M');
ea_busyaction('on',handles.leadfigure,'group');
% set options
options=ea_setopts_local(handles);
options.leadprod = 'group';
options.groupdir = M.root;

% set pt specific options
[options.root, options.patientname] = fileparts(handles.groupdir_choosebox.String);
options.root = [options.root, filesep];

options.expstatvat.do=M.ui.statvat;
options.native=0;

try
    options.numcontacts=size(M.elstruct(1).coords_mm{1},1);
catch
    warning('Localizations seem not properly defined.');
end

options.elmodel=M.elstruct(1).elmodel;
options=ea_resolve_elspec(options);
options.prefs=ea_prefs(options.patientname);
options.d3.verbose='on';

options.d3.hlactivecontacts=get(handles.highlightactivecontcheck,'Value');
options.d3.showactivecontacts=get(handles.showactivecontcheck,'Value');
options.d3.showpassivecontacts=get(handles.showpassivecontcheck,'Value');
options.d3.elrendering=M.ui.elrendering;
options.d3.mirrorsides=M.ui.mirrorsides;
try options.d3.regressorcolormap = M.ui.regressorcolormap; end
try options.d3.isomatrix=M.isomatrix; end
try options.d3.isomatrix_name=M.isomatrix_name; end

options.d2.write=0;

options.d2.atlasopacity=0.15;

options.d3.isovscloud=M.ui.isovscloudpopup;
options.d3.showisovolume=M.ui.showisovolumecheck;

options.d3.colorpointcloud=M.ui.colorpointcloudcheck;
options.d3.exportBB=0;	% don't export brainbrowser struct by default

options.normregressor=M.ui.normregpopup;

% Prepare isomatrix (includes a normalization step if M.ui.normregpopup
% says so:

for reg=1:length(options.d3.isomatrix)
    try
        options.d3.isomatrix{reg}=ea_reformat_isomatrix(options.d3.isomatrix{reg},M,options);
    end
end

if ~strcmp(handles.groupdir_choosebox.String,'Choose Dataset Directory') % group dir still not chosen
    disp('Saving data...');
    % save M
    save(ea_getGroupAnalysisFile(handles.groupdir_choosebox.String),'M','-v7.3');
    disp('Done.');
end

% export VAT-mapping
if options.expstatvat.do % export to nifti volume
    ea_exportvatmapping(M,options,handles);
end
options.groupmode=1;

% overwrite active contacts information with new one from S (if present).
ptidx=get(handles.patientlist,'Value');
try
    for pt=1:length(M.elstruct)
        M.elstruct(pt).activecontacts=M.S(pt).activecontacts;
    end
end

try
    for pt=1:length(M.elstruct)
        M.elstruct(pt).groupcolors=M.groups.color;
    end
end
options.groupmode=1;
options.modality=3; % use template image
options.patient_list=M.patient.list;

% mer development
vizstruct.elstruct=M.elstruct(ptidx);
uipatdirs=handles.patientlist.String(ptidx);
npts=length(uipatdirs);
if options.prefs.env.dev && get(handles.mercheck,'Value')
    filename=fullfile(options.root,options.patientname,'ea_groupvisdata.mat');
    if exist(filename,'file')
       choice = questdlg(sprintf('Group Data Found. Would you like to load %s now?',filename),...
           'Yes','No');
    end

    % Get vizstruct
    if ~exist('choice','var') || strcmpi(choice,'No')
        for pt=1:length(M.elstruct)
            options.uipatdirs{1}=uipatdirs{pt};
            M.merstruct(pt)=ea_getmerstruct(options);
        end

        for pt=1:length(M.elstruct)
            ea_progress(pt/npts, 'Loading microelectrode recordings from patient %d of %d\n', pt, npts);
            M.merstruct(pt).group=handles.grouplist.String(pt);
            [M.merstruct(pt).root,M.merstruct(pt).name]=fileparts(uipatdirs{pt});
            M.merstruct(pt).root(end+1)=filesep;
            M.merstruct(pt).ptdir=uipatdirs{pt};
            mua=load(fullfile(uipatdirs{pt},'ea_recordings.mat'));

            try
                mua.right=rmfield(mua.right,'CSPK');
                mua.left=rmfield(mua.left,'CSPK');
            catch
                mua.right=rmfield(mua.right,'CElectrode');
                mua.left=rmfield(mua.left,'CElectrode');
            end

            M.merstruct(pt).mua=mua;
        end
        disp('**Done loading')
        options = rmfield(options,'uipatdirs');
        vizstruct.merstruct=M.merstruct(ptidx);

        save(fullfile(options.root,options.patientname,'ea_groupelvisdata.mat'),...
            'options','vizstruct');
    elseif strcmpi(choice,'Yes')
        load(filename,'vizstruct')
    end
end

% amend .pt to identify which patient is selected (needed for isomatrix).
for pt=1:length(ptidx)
    M.elstruct(ptidx(pt)).pt=ptidx(pt);
end

elmodels = [{'Patient specified'};ea_resolve_elspec];
whichelmodel = elmodels{M.ui.elmodelselect};
% account for electrode model specified in lead group
if ~strcmp(whichelmodel,'Patient specified')
    arcell=repmat({whichelmodel},length(ptidx),1);
    [M.elstruct(ptidx).elmodel]=arcell{:};
end

resultfig=ea_elvis(options,M.elstruct(ptidx));

try % zoom on coordinates.
    coords={M.elstruct(:).coords_mm};
    for c=1:length(coords)
        call(c,:)=nanmean([coords{c}{1};coords{c}{2}]);
    end
    ea_zoomcenter(resultfig.CurrentAxes, mean(call), 5);
catch
    zoom(3);
end

% show VAT-mapping
if options.expstatvat.do % export to nifti volume
    pobj.plotFigureH=resultfig;
    pobj.color=[0.9,0.2,0.3];

    pobj.openedit=1;
    hshid=ea_datahash(M.ui.listselect);
    ea_roi([options.root,options.patientname,filesep,'statvat_results',filesep,'models',filesep,'statvat_',M.clinical.labels{M.ui.clinicallist},'_T_nthresh_',hshid,'.nii'],pobj);
end

ea_busyaction('off',handles.leadfigure,'group');


% --- Executes on selection change in clinicallist.
function clinicallist_Callback(hObject, eventdata, handles)
% hObject    handle to clinicallist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns clinicallist contents as cell array
%        contents{get(hObject,'Value')} returns selected item from clinicallist
M=getappdata(gcf,'M');

M.ui.clinicallist=get(handles.clinicallist,'Value');
setappdata(gcf,'M',M);
ea_refresh_lg(handles);


% --- Executes during object creation, after setting all properties.
function clinicallist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to clinicallist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in addvarbutton.
function addvarbutton_Callback(hObject, eventdata, handles)
% hObject    handle to addvarbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
M=getappdata(gcf,'M');
M.ui.clinicallist=length(M.clinical.labels)+1;
[numat,nuvar]=ea_get_clinical(M);
if ~isempty(numat) % user did not press cancel
    M.clinical.vars{end+1}=numat;
    M.clinical.labels{end+1}=nuvar;
end
set(handles.clinicallist,'Value',M.ui.clinicallist);
% store model and refresh UI
setappdata(gcf,'M',M);

ea_refresh_lg(handles);


function [mat,matname]=ea_get_clinical(M)
try
    mat=M.clinical.vars{M.ui.clinicallist};
catch % new variable
    mat=[];
end
try
    matname=M.clinical.labels{M.ui.clinicallist};
catch
    matname='New variable';
end
[numat,nuname]=ea_edit_regressor(M);

if ~isempty(numat) % user did not press cancel
    mat=numat;
    matname=nuname;
end


% --- Executes on button press in removevarbutton.
function removevarbutton_Callback(hObject, eventdata, handles)
% hObject    handle to removevarbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
M=getappdata(gcf,'M');

% delete data
M.clinical.vars(get(handles.clinicallist,'Value'))=[];
M.clinical.labels(get(handles.clinicallist,'Value'))=[];

% store model and refresh UI
setappdata(gcf,'M',M);
ea_refresh_lg(handles);


function [pathname] = ea_uigetdir(start_path, dialog_title)
% Pick a directory with the Java widgets instead of uigetdir

import javax.swing.JFileChooser;

if nargin == 0 || strcmp(start_path,'') % || start_path == 0 % Allow a null argument.
    start_path = pwd;
end

jchooser = javaObjectEDT('javax.swing.JFileChooser', start_path);

jchooser.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
if nargin > 1
    jchooser.setDialogTitle(dialog_title);
end

jchooser.setMultiSelectionEnabled(true);

status = jchooser.showOpenDialog([]);

if status == JFileChooser.APPROVE_OPTION
    jFile = jchooser.getSelectedFiles();
    pathname{size(jFile, 1)}=[];
    for i=1:size(jFile, 1)
        pathname{i} = char(jFile(i).getAbsolutePath);
    end

elseif status == JFileChooser.CANCEL_OPTION
    pathname = [];
else
    error('Error occured while picking file.');
end


% --- Executes on selection change in grouplist.
function grouplist_Callback(hObject, eventdata, handles)
% hObject    handle to grouplist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns grouplist contents as cell array
%        contents{get(hObject,'Value')} returns selected item from grouplist
M=getappdata(gcf,'M');

M.ui.listselect=get(handles.grouplist,'Value');

set(handles.patientlist,'Value',M.ui.listselect);

setappdata(gcf,'M',M);
ea_refresh_lg(handles);


% --- Executes during object creation, after setting all properties.
function grouplist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to grouplist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in plusgroupbutton.
function plusgroupbutton_Callback(hObject, eventdata, handles)
% hObject    handle to plusgroupbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
M=getappdata(gcf,'M');
M.patient.group(get(handles.patientlist,'Value'))=M.patient.group(get(handles.patientlist,'Value'))+1;
setappdata(gcf,'M',M);
ea_refresh_lg(handles);


% --- Executes on button press in minusgroupbutton.
function minusgroupbutton_Callback(hObject, eventdata, handles)
% hObject    handle to minusgroupbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
M=getappdata(gcf,'M');
if M.patient.group(get(handles.patientlist,'Value'))>1
    M.patient.group(get(handles.patientlist,'Value'))=M.patient.group(get(handles.patientlist,'Value'))-1;
end
setappdata(gcf,'M',M);
ea_refresh_lg(handles);


% --- Executes on button press in reviewvarbutton.
function reviewvarbutton_Callback(hObject, eventdata, handles)
% hObject    handle to reviewvarbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
M=getappdata(gcf,'M');

% store in model as variables

%M.clinical.vars{get(handles.clinicallist,'Value')}(isnan(M.clinical.vars{get(handles.clinicallist,'Value')}))=0;
[M.clinical.vars{get(handles.clinicallist,'Value')},M.clinical.labels{get(handles.clinicallist,'Value')}]=ea_get_clinical(M);


% store model and refresh UI
setappdata(gcf,'M',M);

ea_refresh_lg(handles);


% --- Executes on button press in moveptupbutton.
function moveptupbutton_Callback(hObject, eventdata, handles)
% hObject    handle to moveptupbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
M=getappdata(gcf,'M');
whichmoved=get(handles.patientlist,'Value');

if whichmoved(1)==1 % first entry anyways
    return
end

ix=1:length(M.patient.list);
ix(whichmoved)=ix(whichmoved)-1;
ix(whichmoved-1)=ix(whichmoved-1)+1;

M.patient.list=M.patient.list(ix);
M.patient.group=M.patient.group(ix);
M.ui.listselect=whichmoved-1;
for c=1:length(M.clinical.vars)
    M.clinical.vars{c} = M.clinical.vars{c}(ix,:);
end
try
    M.S = M.S(ix);
end
try
    M=rmfield(M,'elstruct');
end
try
    M=rmfield(M,'stats');
end
setappdata(gcf,'M',M);

set(handles.patientlist,'Value',whichmoved-1);
ea_refresh_lg(handles);


% --- Executes on button press in moveptdownbutton.
function moveptdownbutton_Callback(hObject, eventdata, handles)
% hObject    handle to moveptdownbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

M=getappdata(gcf,'M');
whichmoved=get(handles.patientlist,'Value');

if whichmoved(end)==length(M.patient.list) % last entry anyways
    return
end

ix=1:length(M.patient.list);
ix(whichmoved)=ix(whichmoved)+1;
ix(whichmoved+1)=ix(whichmoved+1)-1;

M.patient.list=M.patient.list(ix);
M.patient.group=M.patient.group(ix);
M.ui.listselect=whichmoved+1;
for c=1:length(M.clinical.vars)
    M.clinical.vars{c} = M.clinical.vars{c}(ix,:);
end
try
    M.S = M.S(ix);
end
try
    M=rmfield(M,'elstruct');
end
try
    M=rmfield(M,'stats');
end
setappdata(gcf,'M',M);

set(handles.patientlist,'Value',whichmoved+1);
ea_refresh_lg(handles);


% --- Executes on button press in calculatebutton.
function calculatebutton_Callback(hObject, eventdata, handles)
% hObject    handle to calculatebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ea_refresh_lg(handles);

M=getappdata(gcf,'M');

% set options
options=ea_setopts_local(handles);
%stimname=ea_detstimname(options);

options.groupmode = 1;
options.groupid = M.guid;
options.groupdir = M.root;

if isfield(M.ui, 'stimSetMode') && M.ui.stimSetMode
    options.stimSetMode = 1;
else
    options.stimSetMode = 0;
end

% determine if fMRI or dMRI
mods=get(handles.fiberspopup,'String');
mod=mods{get(handles.fiberspopup,'Value')};
switch mod
    case {'Patient''s fiber tracts', 'Patient''s fMRI time courses'}
        fibersfile=mod;
    case 'Do not calculate connectivity stats'
    otherwise % load fibertracts once and for all subs here.
        [fibersfile.fibers,fibersfile.fibersidx]=ea_loadfibertracts([ea_getconnectomebase('dmri'),mod,filesep,'data.mat']);
end

[selection]=ea_groupselectorwholelist(M.ui.listselect,M.patient.list);

for pt=selection
    % set pt specific options
    [options.root, options.patientname] = fileparts(M.patient.list{pt});
    options.root = [options.root, filesep];

    options = ea_getptopts(fullfile(options.root, options.patientname), options);

    fprintf('\nProcessing %s...\n\n', options.patientname);
    try
        options.numcontacts=size(M.elstruct(pt).coords_mm{1},1);
    catch % no localization present or in wrong format.
        ea_error(['Please localize ',options.patientname,' first.']);
    end
    options.d3.verbose='off';
    options.d3.elrendering=1;	% hard code to viz electrodes in this setting.
    options.d3.exportBB=0;	% don't export brainbrowser struct by default
    options.d3.colorpointcloud=0;

    options.d3.hlactivecontacts=get(handles.highlightactivecontcheck,'Value');
    options.d3.showactivecontacts=get(handles.showactivecontcheck,'Value');
    options.d3.showpassivecontacts=get(handles.showpassivecontcheck,'Value');
    try
        options.d3.isomatrix=M.isomatrix;
    catch
        options.d3.isomatrix={};
    end
    try
        options.d3.isomatrix_name=M.isomatrix_name;
    catch
        options.d3.isomatrix_name={};
    end
    options.normregressor=M.ui.normregpopup;
    for reg=1:length(options.d3.isomatrix)
        try
            options.d3.isomatrix{reg}=ea_reformat_isomatrix(options.d3.isomatrix{reg},M,options);
        end
    end

    options.d3.isovscloud=M.ui.isovscloudpopup;
    options.d3.showisovolume=M.ui.showisovolumecheck;
    options.d3.exportBB=0;
    options.expstatvat.do=0;
    try
        options.expstatvat.vars=M.clinical.vars(M.ui.clinicallist);
        options.expstatvat.labels=M.clinical.labels(M.ui.clinicallist);
        options.expstatvat.pt=pt;
    end
    options.expstatvat.dir=M.root;

    % Step 1: Re-calculate closeness to subcortical atlases.
    options.leadprod = 'group';
    options.patient_list=M.patient.list;
    options.d3.mirrorsides=0;

    resultfig=ea_elvis(options,M.elstruct(pt));

    if ~isfield(options.subj, 'norm')
        ea_cprintf('CmdWinWarnings', 'Running in Miniset mode: %s...\n', options.subj.subjId);
        volumespresent=0;
    elseif isempty(dir([options.subj.norm.transform.inverseBaseName, '*']))
        ea_cprintf('CmdWinWarnings', 'Tranformation not found for %s...\n', options.subj.subjId);
        volumespresent=0;
    else
        volumespresent=1;
    end

    % Step 2: Re-calculate VAT
    if isfield(M,'S')
        try
            setappdata(resultfig,'curS',M.S(pt));
        catch
            ea_error(['Stimulation parameters for ',M.patient.list{pt},' are missing.']);
        end
        vfnames=getappdata(handles.leadfigure,'vatfunctionnames');

        [~,ix]=ismember(M.vatmodel,vfnames);
        vfs=getappdata(handles.leadfigure,'genvatfunctions');
        try
            ea_genvat=eval(['@',vfs{ix}]);
        catch
            keyboard
        end

        options=getappdata(resultfig,'options'); % selected atlas could have refreshed.

        options.orignative=options.native; % backup
        options.native=~ea_getprefs('vatsettings.estimateInTemplate'); % see whether VTAs should be directly estimated in template space or not
        if options.native && ~volumespresent
            ea_cprintf('CmdWinWarnings', 'Calculating VTA in template space since patient folder %s is incomplete.\n', options.subj.subjId);
            options.native=0;
        end

        setappdata(handles.leadfigure,'resultfig',resultfig);
        setappdata(resultfig,'elstruct',M.elstruct(pt));
        setappdata(resultfig,'elspec',options.elspec);

        if options.native % Reload native space coordinates
            coords = ea_load_reconstruction(options);
        else
            coords = M.elstruct(pt).coords_mm;
        end

        vatCalcPassed = [0 0];
        stimparams = struct();
        if strcmp(M.S(pt).model, 'OSS-DBS (Butenko 2020)')
            if options.prefs.machine.vatsettings.butenko_calcAxonActivation
                feval(ea_genvat,M.S(pt),options);
                fprintf('\n');
                warning('off', 'backtrace');
                warning('OSS-DBS axon activation mode detect, skipping calc stats for %s!\n', options.patientname);
                warning('on', 'backtrace');
                continue;
            else
                [vatCalcPassed, stimparams] = feval(ea_genvat,M.S(pt),options);
            end
        else
            for side=1:2
                try
                    [vtafv,vtavolume] = feval(ea_genvat,coords,M.S(pt),side,options,['gs_',M.guid],handles.leadfigure);
                    vatCalcPassed(side) = 1;
                catch
                    vtafv=[];
                    vtavolume=0;
                    vatCalcPassed(side) = 0;
                end
                stimparams(1,side).VAT(1).VAT = vtafv;
                stimparams(1,side).volume = vtavolume;
            end
        end

        options.native=options.orignative; % restore
        setappdata(resultfig,'stimparams',stimparams(1,:));
    end
    % Calc VAT stats (atlas intersection and volume)
    if all(vatCalcPassed)
        ea_calc_vatstats(resultfig,options);
    else
        try
            ea_error(sprintf(['An error occured when building the VTA mesh/headmodel for %s.\n',...
                'Try re-calculating this VTA with a different atlas or with no atlas.'],...
                options.patientname));
        catch
            continue;
        end
    end

    % Step 3: Re-calculate connectivity from VAT to rest of the brain.
    if all(vatCalcPassed) && ~strcmp(mod,'Do not calculate connectivity stats')
        % Convis part:
        parcs=get(handles.labelpopup,'String');
        selectedparc=parcs{get(handles.labelpopup,'Value')};
        directory=[options.root,options.patientname,filesep];
        if ischar(fibersfile)
            switch mod
                case 'Patient''s fMRI time courses'
                    ea_error('Group statistics for fMRI are not yet supported. Sorry, check back later!');
                    pV=spm_vol([ea_space(options,'labeling'),selectedparc,'.nii']);
                    pX=spm_read_vols(pV);
                    ea_cvshowvatfmri(resultfig,pX,directory,filesare,handles,pV,selectedparc,mod,options);
                otherwise
                    ea_cvshowvatdmri(resultfig,directory,{mod,'gs'},selectedparc,options);
            end
        else
            ea_cvshowvatdmri(resultfig,directory,{fibersfile,'gs'},selectedparc,options);
        end
    end

    close(resultfig);
end
%% processing done here.
ea_refresh_lg(handles);

set(handles.calculatebutton, 'BackgroundColor', [0.93,0.93,0.93]);
set(handles.explorestats, 'Enable', 'on');
set(handles.exportstats, 'Enable', 'on');


% --- Executes on selection change in fiberspopup.
function fiberspopup_Callback(hObject, eventdata, handles)
% hObject    handle to fiberspopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns fiberspopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from fiberspopup
M=getappdata(gcf,'M');

if isfield(M.ui, 'connectomename') && ...
   strcmp(eventdata.Source.String{eventdata.Source.Value}, M.ui.connectomename)
    connChanged = 0;
else
    connChanged = 1;
end

M.ui.connectomename = eventdata.Source.String{eventdata.Source.Value};
setappdata(gcf,'M',M);

ea_refresh_lg(handles);

if ~isempty(M.patient.list) && connChanged
    set(handles.calculatebutton, 'BackgroundColor', [0.1;0.8;0.1]);
    set(handles.explorestats, 'Enable', 'off');
    set(handles.exportstats, 'Enable', 'off');
end


% --- Executes during object creation, after setting all properties.
function fiberspopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fiberspopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in labelpopup.
function labelpopup_Callback(hObject, eventdata, handles)
% hObject    handle to labelpopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns labelpopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from labelpopup
M=getappdata(gcf,'M');

labelChanged = 1;
if isfield(M.ui, 'labelpopup')
    if isnumeric(M.ui.labelpopup) % Old format, labelpopup is numeric
        if eventdata.Source.Value == M.ui.labelpopup
            labelChanged = 0;
        end
    else % New format, labelpopup is labeling parcellation name
        if strcmp(eventdata.Source.String{eventdata.Source.Value}, M.ui.labelpopup)
            labelChanged = 0;
        end
    end
end

M.ui.labelpopup = eventdata.Source.String{eventdata.Source.Value};
setappdata(gcf,'M',M);
ea_refresh_lg(handles);

if ~isempty(M.patient.list) && labelChanged
    set(handles.calculatebutton, 'BackgroundColor', [0.1;0.8;0.1]);
    set(handles.explorestats, 'Enable', 'off');
    set(handles.exportstats, 'Enable', 'off');
end


% --- Executes during object creation, after setting all properties.
function labelpopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to labelpopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in atlassetpopup.
function atlassetpopup_Callback(hObject, eventdata, handles)
% hObject    handle to atlassetpopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns atlassetpopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from atlassetpopup
M=getappdata(gcf,'M');

atlasChanged = 1;
if isfield(M.ui, 'atlassetpopup')
    if isnumeric(M.ui.atlassetpopup) % Old format, atlassetpopup is numeric
        if eventdata.Source.Value == M.ui.atlassetpopup
            atlasChanged = 0;
        end
    else % New format, atlassetpopup is atlas name
        if strcmp(eventdata.Source.String{eventdata.Source.Value}, M.ui.atlassetpopup)
            atlasChanged = 0;
        end
    end
end

M.ui.atlassetpopup = eventdata.Source.String{eventdata.Source.Value};
setappdata(gcf,'M',M);
ea_refresh_lg(handles);

if ~isempty(M.patient.list) && atlasChanged
    set(handles.calculatebutton, 'BackgroundColor', [0.1;0.8;0.1]);
    set(handles.explorestats, 'Enable', 'off');
    set(handles.exportstats, 'Enable', 'off');
end


% --- Executes during object creation, after setting all properties.
function atlassetpopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to atlassetpopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function options=ea_setopts_local(handles)

options.earoot=ea_getearoot;
options.verbose=3;
options.sides=1:2; % re-check this later..
options.atlasset=get(handles.atlassetpopup,'String');
try
    options.atlasset=options.atlasset{get(handles.atlassetpopup,'Value')};
catch % too many entries..
    set(handles.atlassetpopup,'Value',1);
    options.atlasset=1;
end
options.fiberthresh=1;
options.writeoutstats=1;
options.labelatlas=get(handles.labelpopup,'String');
try
    options.labelatlas=options.labelatlas{get(handles.labelpopup,'Value')};
catch % too many entries..
    set(handles.labelpopup,'Value',1);
    options.labelatlas=1;
end
options.writeoutpm = 0;
options.colormap=parula(64);
options.d3.write=1;
options.d3.prolong_electrode=2;
options.d3.writeatlases=1;
options.macaquemodus=0;


% --- Executes on button press in opensubgui.
function opensubgui_Callback(hObject, eventdata, handles)
% hObject    handle to opensubgui (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
M=getappdata(gcf,'M');
[selection]=ea_groupselectorwholelist(M.ui.listselect,M.patient.list);

lead_dbs('loadsubs',M.patient.list(selection));


% --- Executes on button press in choosegroupcolors.
function choosegroupcolors_Callback(hObject, eventdata, handles)
% hObject    handle to choosegroupcolors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

M=getappdata(gcf,'M');

for g=unique(M.patient.group)'
    M.groups.color(ismember(M.groups.group,g),:)=...
        ea_uisetcolor(M.groups.color(ismember(M.groups.group,g),:),['Group ',num2str(g),':']);
end
M.groups.colorschosen=1;

setappdata(gcf,'M',M);
ea_refresh_lg(handles);


% --- Executes on button press in setstimparamsbutton.
function setstimparamsbutton_Callback(hObject, eventdata, handles)
% hObject    handle to setstimparamsbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
M=getappdata(gcf,'M');

% try
%     uicell=inputdlg('Enter Variable name for Voltage-Parameters','Enter Stimulation Settings...',1);
%     uidata.U=evalin('base',uicell{1});
% catch
%     warning('Stim-Params could not be evaluated. Please Try again.');
%     return
% end
% try
%     uicell=inputdlg('Enter Variable name for Impedance-Parameters','Enter Stimulation Settings...',1);
%     uidata.Im=evalin('base',uicell{1});
% catch
%     warning('Stim-Params could not be evaluated. Please Try again.');
%     return
% end

options = ea_setopts_local(handles);
options.leadprod = 'group';
options.groupid = M.guid;
options.native = 0;
ea_refresh_lg(handles);

ea_stimparams(M.elstruct, handles.leadfigure, options);


% --- Executes on button press in highlightactivecontcheck.
function highlightactivecontcheck_Callback(hObject, eventdata, handles)
% hObject    handle to highlightactivecontcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of highlightactivecontcheck
M=getappdata(gcf,'M');
M.ui.hlactivecontcheck=get(handles.highlightactivecontcheck,'Value');

setappdata(gcf,'M',M);
ea_refresh_lg(handles);


% --- Executes on button press in showpassivecontcheck.
function showpassivecontcheck_Callback(hObject, eventdata, handles)
% hObject    handle to showpassivecontcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showpassivecontcheck

M=getappdata(gcf,'M');
M.ui.showpassivecontcheck=get(handles.showpassivecontcheck,'Value');


setappdata(gcf,'M',M);
ea_refresh_lg(handles);


% --- Executes on button press in showactivecontcheck.
function showactivecontcheck_Callback(hObject, eventdata, handles)
% hObject    handle to showactivecontcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showactivecontcheck

M=getappdata(gcf,'M');
M.ui.showactivecontcheck=get(handles.showactivecontcheck,'Value');

setappdata(gcf,'M',M);
ea_refresh_lg(handles);


% --- Executes on button press in showisovolumecheck.
function showisovolumecheck_Callback(hObject, eventdata, handles)
% hObject    handle to showisovolumecheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showisovolumecheck
M=getappdata(gcf,'M');
M.ui.showisovolumecheck=get(handles.showisovolumecheck,'Value');

setappdata(gcf,'M',M);
ea_refresh_lg(handles);


% --- Executes on button press in statvatcheck.
function statvatcheck_Callback(hObject, eventdata, handles)
% hObject    handle to statvatcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of statvatcheck
M=getappdata(gcf,'M');

M.ui.statvat=get(handles.statvatcheck,'Value');
setappdata(gcf,'M',M);


% --- Executes on button press in mercheck.
function mercheck_Callback(hObject, eventdata, handles)
% hObject    handle to mercheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of mercheck

M=getappdata(gcf,'M');

M.ui.mer=get(handles.mercheck,'Value');
setappdata(gcf,'M',M);


% --- Executes on selection change in normregpopup.
function normregpopup_Callback(hObject, eventdata, handles)
% hObject    handle to normregpopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns normregpopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from normregpopup
M=getappdata(gcf,'M');
M.ui.normregpopup=get(handles.normregpopup,'Value');
setappdata(gcf,'M',M);
ea_refresh_lg(handles);


% --- Executes during object creation, after setting all properties.
function normregpopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to normregpopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when user attempts to close leadfigure.
function leadfigure_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to leadfigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure

ea_busyaction('on',gcf,'group');
if ~strcmp(handles.groupdir_choosebox.String,'Choose Dataset Directory') % group dir still not chosen
    disp('Saving data...');
    % save M
    ea_refresh_lg(handles);
    M=getappdata(hObject,'M');
    disp('Saving data to disk...');
    try
        save(ea_getGroupAnalysisFile(handles.groupdir_choosebox.String),'M','-v7.3');
    catch
        warning('Data could not be saved.');
        keyboard
    end
    disp('Done.');
    disp('Bye for now.');
end
ea_busyaction('off',gcf,'group');
delete(hObject);


% --- Executes on button press in explorestats.
function explorestats_Callback(hObject, eventdata, handles)
% hObject    handle to explorestats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ea_refresh_lg(handles);
ea_lg_stats(handles.leadfigure);


% --- Executes on button press in viz2dbutton.
function viz2dbutton_Callback(hObject, eventdata, handles)
% hObject    handle to viz2dbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc;
M=getappdata(gcf,'M');

ea_busyaction('on',gcf,'group');

% set options
options=ea_setopts_local(handles);

% set pt specific options
[options.root, options.patientname] = fileparts(handles.groupdir_choosebox.String);
options.root = [options.root, filesep];

options.numcontacts=size(M.elstruct(1).coords_mm{1},1);
options.elmodel=M.elstruct(1).elmodel;
options=ea_resolve_elspec(options);
options.prefs=ea_prefs(options.patientname);
options.d3.verbose='on';

options.d3.elrendering=M.ui.elrendering;
options.d3.hlactivecontacts=get(handles.highlightactivecontcheck,'Value');
options.d3.showactivecontacts=get(handles.showactivecontcheck,'Value');
options.d3.showpassivecontacts=get(handles.showpassivecontcheck,'Value');
try
    options.d3.isomatrix=M.isomatrix;
catch
    options.d3.isomatrix={};
end
try
    options.d3.isomatrix_name=M.isomatrix_name;
catch
    options.d3.isomatrix_name={};
end

options.expstatvat.do=M.ui.statvat;

options.d2.showlegend=0;

options.d3.isovscloud=M.ui.isovscloudpopup;
options.d3.showisovolume=M.ui.showisovolumecheck;
options.d3.colorpointcloud=M.ui.colorpointcloudcheck;
options.normregressor=M.ui.normregpopup;

options.d2.write=1;

options.d2.atlasopacity=0.15;
options.groupmode=1;
options.groupid=M.guid;
options.modality=3; % use template image
options=ea_amendtoolboxoptions(options);

if strcmp(options.atlasset,'Use none')
    options.d2.writeatlases=0;
else
    options.d2.writeatlases=1;
end

% Prior Results are loaded here inside the function (this way, function
% can be called just by giving the patient directory.

% Prepare isomatrix (includes a normalization step if M.ui.normregpopup
% says so:

if options.d3.showisovolume || options.expstatvat.do % regressors be used - iterate through all
    allisomatrices=options.d3.isomatrix;
    allisonames=options.d3.isomatrix_name;
    for reg=1:length(allisomatrices)
        options.d3.isomatrix=allisomatrices{reg};
        options.d3.isomatrix_name=allisonames{reg};
        M.isomatrix=allisomatrices{reg};
        M.isomatrix_name=allisonames{reg};
        options.shifthalfup=0;
        try
            options.d3.isomatrix=ea_reformat_isomatrix(options.d3.isomatrix,M,options);
            if size(options.d3.isomatrix{1},2)==3 % pairs
                options.shifthalfup=1;
            end
        end

        if ~strcmp(handles.groupdir_choosebox.String,'Choose Dataset Directory') % group dir still not chosen
            ea_refresh_lg(handles);
            disp('Saving data...');
            % save M
            save(ea_getGroupAnalysisFile(handles.groupdir_choosebox.String),'M','-v7.3');
            disp('Done.');
        end

        % export coordinate-mapping
        if options.d3.showisovolume % export to nifti volume
            ea_exportisovolume(M.elstruct(get(handles.patientlist,'Value')),options);
        end
        % export VAT-mapping
        if options.expstatvat.do % export to nifti volume
            ea_exportvatmapping(M,options,handles);
        end

        ea_out2d(M,options,handles);
    end
else
    ea_out2d(M,options,handles);
end
ea_busyaction('off',gcf,'group');


function ea_out2d(M,options,handles)

for pt=1:length(M.patient.list)
    for side=1:2
        try
            M.elstruct(pt).activecontacts{side}=M.S(pt).activecontacts{side};
        end
    end
end
cuts=ea_writeplanes(options,M.elstruct(get(handles.patientlist,'Value')));


% --- Executes on button press in specify2doptions.
function specify2doptions_Callback(hObject, eventdata, handles)
% hObject    handle to specify2doptions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
options.prefs=ea_prefs('');
options.groupmode=1;
options.native=0;
ea_spec2dwrite(options);


% --- Executes on selection change in recentgroups.
function recentgroups_Callback(hObject, eventdata, handles)
% hObject    handle to recentgroups (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns recentgroups contents as cell array
%        contents{get(hObject,'Value')} returns selected item from recentgroups
ea_busyaction('on',handles.leadfigure,'group');
ea_recentcallback(handles, 'groups');
ea_busyaction('off',handles.leadfigure,'group');


% --- Executes during object creation, after setting all properties.
function recentgroups_CreateFcn(hObject, eventdata, handles)
% hObject    handle to recentgroups (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in specify3doptions.
function specify3doptions_Callback(hObject, eventdata, handles)
% hObject    handle to specify3doptions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ea_lg_3dsetting(handles.leadfigure)


% --- Executes on button press in exportstats.
function exportstats_Callback(hObject, eventdata, handles)
% hObject    handle to exportstats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
M = getappdata(gcf,'M');
exportFile = strrep(ea_getGroupAnalysisFile(M.root), '.mat', '_desc-stats_export.mat');
[file, path] = uiputfile('*.mat','Export DBS Stats as...', exportFile);
if file % make sure user didnt press cancel
    ea_lg_exportstats(M, [path, file]);
    fprintf('\nDBS Stats exported to:\n%s\n\n', [path, file]);
end



% --- Executes on button press in minisetbutton.
function minisetbutton_Callback(hObject, eventdata, handles)
% hObject    handle to minisetbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
M = getappdata(gcf,'M');
ea_generate_min_dataset(M);
