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

% Last Modified by GUIDE v2.5 18-Jun-2016 08:34:44

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

% Choose default command line output for lead_group
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes lead_group wait for user response (see UIRESUME)
% uiwait(handles.lg_figure);


% Build popup tables:

% atlassets:
options.earoot=[ea_getearoot];
setappdata(handles.lg_figure,'earoot',options.earoot);
as=dir([options.earoot,'atlases',filesep]);
asc=cell(0);
cnt=1;
for i=1:length(as)
    if as(i).isdir
        asc{cnt}=as(i).name;
        cnt=cnt+1;
    end
end
excludes={'.','..'};
asc(ismember(asc,excludes))=[];

asc{end+1}='Use none';

set(handles.atlassetpopup,'String',asc);



% setup vat functions

cnt=1;
earoot=[ea_getearoot];
ndir=dir([earoot,'ea_genvat_*.m']);
for nd=length(ndir):-1:1
    [~,methodf]=fileparts(ndir(nd).name);
    try
        [thisndc]=eval([methodf,'(','''prompt''',')']);
        ndc{cnt}=thisndc;
        genvatfunctions{cnt}=methodf;
        cnt=cnt+1;
    end
end
setappdata(gcf,'genvatfunctions',genvatfunctions);
setappdata(gcf,'vatfunctionnames',ndc);




% get electrode model specs and place in popup
set(handles.elmodelselect,'String',[{'Patient specified'},ea_resolve_elspec]);

% set background image
set(gcf,'color','w');
im=imread('ea_logo.png');
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


ll=dir([options.earoot,'templates',filesep,'labeling',filesep,'*.nii']);
for lab=1:length(ll)
    [~,n]=fileparts(ll(lab).name);
    labelcell{lab}=n;
end

set(handles.labelpopup,'String',labelcell);

try
    priorselection=find(ismember(labelcell,stimparams.labelatlas)); % retrieve prior selection of fiberset.
    if length(priorselection)==1
        set(handles.labelpopup,'Value',priorselection); % set to prior selection
    else % if priorselection was a cell array with more than one entry, set to use all
        set(handles.labelpopup,'Value',lab+1); % set to use all

    end
catch    % reinitialize using third entry.
    set(handles.labelpopup,'Value',1);


end


% set version text:
set(handles.versiontxt,'String',['v',ea_getvsn('local')]);


% make listboxes multiselectable:

set(handles.patientlist,'Max',100,'Min',0);
set(handles.grouplist,'Max',100,'Min',0);
set(handles.vilist,'Max',100,'Min',0);
set(handles.fclist,'Max',100,'Min',0);
set(handles.clinicallist,'Max',100,'Min',0);


M=getappdata(gcf,'M');
if isempty(M)
    % initialize Model variable M
    M=initializeM;
end
setappdata(gcf,'M',M);
refreshvifc(handles);




% --- Outputs from this function are returned to the command line.
function varargout = lead_group_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


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
refreshvifc(handles);



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
M=getappdata(gcf,'M');

folders=ea_uigetdir('/','Select Patient folders..');
M.patient.list=[M.patient.list;folders'];
M.patient.group=[M.patient.group;ones(length(folders),1)];

setappdata(gcf,'M',M);
refreshvifc(handles);
% save M
M=getappdata(gcf,'M');
save([get(handles.groupdir_choosebox,'String'),'LEAD_groupanalysis.mat'],'M','-v7.3');





% --- Executes on button press in removeptbutton.
function removeptbutton_Callback(hObject, eventdata, handles)
% hObject    handle to removeptbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
M=getappdata(gcf,'M');

deleteentry=get(handles.patientlist,'Value');

M.patient.list(deleteentry)=[];

M.patient.group(deleteentry)=[];

try M.elstruct(deleteentry)=[]; end
try    M.stimparams(deleteentry)=[]; end

for cvar=1:length(M.clinical.vars)
    M.clinical.vars{cvar}(deleteentry,:)=[];
end

try
    M.stats(deleteentry)=[];
end


setappdata(gcf,'M',M);
refreshvifc(handles);




% --- Executes on button press in vizbutton.
function vizbutton_Callback(hObject, eventdata, handles)
% hObject    handle to vizbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc;
M=getappdata(gcf,'M');
ea_busyaction('on',handles.lg_figure,'group');
% set options
options=ea_setopts_local(handles);
% set pt specific options
options.root=[fileparts(fileparts(get(handles.groupdir_choosebox,'String'))),filesep];
[~,options.patientname]=fileparts(fileparts(get(handles.groupdir_choosebox,'String')));


options.expstatvat.do=M.ui.statvat;

try

options.numcontacts=size(M.elstruct(1).coords_mm{1},1);
catch
    warning('Localizations seem not properly defined.');
end
options.elmodel=M.elstruct(1).elmodel;
options=ea_resolve_elspec(options);
options.prefs=ea_prefs(options.patientname);
options.d3.verbose='on';

options.d3.elrendering=M.ui.elrendering;
options.d3.hlactivecontacts=get(handles.highlightactivecontcheck,'Value');
options.d3.showactivecontacts=get(handles.showactivecontcheck,'Value');
options.d3.showpassivecontacts=get(handles.showpassivecontcheck,'Value');
try options.d3.isomatrix=M.isomatrix; end
try options.d3.isomatrix_name=M.isomatrix_name; end


options.d2.write=0;

options.d2.atlasopacity=0.15;

options.d3.isovscloud=M.ui.isovscloudpopup;
options.d3.showisovolume=M.ui.showisovolumecheck;

options.d3.colorpointcloud=M.ui.colorpointcloudcheck;
options.normregressor=M.ui.normregpopup;

% Prepare isomatrix (includes a normalization step if M.ui.normregpopup
% says so:

try options.d3.isomatrix=ea_reformat_isomatrix(options.d3.isomatrix,M,options); end

if ~strcmp(get(handles.groupdir_choosebox,'String'),'Choose Group Directory') % group dir still not chosen
    disp('Saving data...');
    % save M
    save([get(handles.groupdir_choosebox,'String'),'LEAD_groupanalysis.mat'],'M','-v7.3');
    disp('Done.');
end

% export VAT-mapping
if options.expstatvat.do % export to nifti volume
    ea_exportvatmapping(M,options,handles);
end

% overwrite active contacts information with new one from S (if present).
try
    for pt=1:length(M.elstruct)
        M.elstruct(pt).activecontacts=M.S(pt).activecontacts;
    end
end
resultfig=ea_elvis(options,M.elstruct(get(handles.patientlist,'Value')));

ea_busyaction('off',handles.lg_figure,'group');


% --- Executes on button press in corrbutton.
function corrbutton_Callback(hObject, eventdata, handles)
% hObject    handle to corrbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ea_busyaction('on',gcf,'group');

stats=preparedataanalysis(handles);


assignin('base','stats',stats);


% perform correlations:


if size(stats.corrcl,2)==1 % one value per patient
    try stats.vicorr.nboth=(stats.vicorr.nboth/2)*100; end
    try stats.vicorr.nright=(stats.vicorr.nright/2)*100; end
    try stats.vicorr.nleft=(stats.vicorr.nleft/2)*100; end

    if ~isempty(stats.vicorr.both)
        %ea_corrplot([stats.corrcl,stats.vicorr.both],'Volume Intersections, both hemispheres',stats.vc_labels);
        ea_corrplot([stats.corrcl,stats.vicorr.nboth],'VI_BH',stats.vc_labels,handles);
    end
    %     if ~isempty(stats.vicorr.right)
    %         %ea_corrplot([stats.corrcl,stats.vicorr.right],'Volume Intersections, right hemisphere',stats.vc_labels);
    %         ea_corrplot([stats.corrcl,stats.vicorr.nright],'VI_RH',stats.vc_labels,handles);
    %     end
    %     if ~isempty(stats.vicorr.left)
    %         %ea_corrplot([stats.corrcl,stats.vicorr.left],'Volume Intersections, left hemisphere',stats.vc_labels);
    %         ea_corrplot([stats.corrcl,stats.vicorr.nleft],'VI_LH',stats.vc_labels,handles);
    %     end
    %     if ~isempty(stats.fccorr.both)
    %         %ea_corrplot([stats.corrcl,stats.fccorr.both],'Fibercounts, both hemispheres',stats.fc_labels);
    %         ea_corrplot([stats.corrcl,stats.fccorr.nboth],'FC_BH',stats.fc_labels,handles);
    %     end
    %     if ~isempty(stats.fccorr.right)
    %         %ea_corrplot([stats.corrcl,stats.fccorr.right],'Fibercounts, right hemisphere',stats.fc_labels);
    %         ea_corrplot([stats.corrcl,stats.fccorr.nright],'FC_RH',stats.fc_labels,handles);
    %     end
    %     if ~isempty(stats.fccorr.left)
    %         %ea_corrplot([stats.corrcl,stats.fccorr.left],'Fibercounts, left hemisphere',stats.fc_labels);
    %         ea_corrplot([stats.corrcl,stats.fccorr.nleft],'FC_LH',stats.fc_labels,handles);
    %     end

elseif size(stats.corrcl,2)==2 % one value per hemisphere
    try stats.vicorr.nboth=(stats.vicorr.nboth)*100; end
    try stats.vicorr.nright=(stats.vicorr.nright)*100; end
    try stats.vicorr.nleft=(stats.vicorr.nleft)*100; end
    if ~isempty(stats.vicorr.both)

        %ea_corrplot([stats.corrcl(:),[stats.vicorr.right;stats.vicorr.left]],'Volume Intersections, both hemispheres',stats.vc_labels);
        ea_corrplot([stats.corrcl(:),[stats.vicorr.nright;stats.vicorr.nleft]],'VI_BH',stats.vc_labels,handles);
    end
    %     if ~isempty(stats.vicorr.right)
    %         %ea_corrplot([stats.corrcl(:,1),stats.vicorr.right],'Volume Intersections, right hemisphere',stats.vc_labels);
    %         ea_corrplot([stats.corrcl(:,1),stats.vicorr.nright],'VI_RH',stats.vc_labels,handles);
    %     end
    %     if ~isempty(stats.vicorr.left)
    %         %ea_corrplot([stats.corrcl(:,2),stats.vicorr.left],'Volume Intersections, left hemisphere',stats.vc_labels);
    %         ea_corrplot([stats.corrcl(:,2),stats.vicorr.nleft],'VI_LH',stats.vc_labels,handles);
    %     end
    %     if ~isempty(stats.fccorr.both)
    %         %ea_corrplot([stats.corrcl(:),[stats.fccorr.right;stats.fccorr.left]],'Fibercounts, both hemispheres',stats.fc_labels);
    %         ea_corrplot([stats.corrcl(:),[stats.fccorr.right;stats.fccorr.left]],'FC_BH',stats.fc_labels,handles);
    %     end
    %     if ~isempty(stats.fccorr.right)
    %         %ea_corrplot([stats.corrcl(:,1),stats.fccorr.right],'Fibercounts, right hemisphere',stats.fc_labels);
    %         ea_corrplot([stats.corrcl(:,1),stats.fccorr.nright],'FC_RH',stats.fc_labels,handles);
    %     end
    %     if ~isempty(stats.fccorr.left)
    %         %ea_corrplot([stats.corrcl(:,2),stats.fccorr.left],'Fibercounts, left hemisphere',stats.fc_labels);
    %         ea_corrplot([stats.corrcl(:,2),stats.fccorr.nleft],'FC_LH',stats.fc_labels,handles);
    %     end
    %
else
    ea_error('Please select a regressor with one value per patient or per hemisphere to perform this correlation.');
end
ea_busyaction('off',gcf,'group');



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
refreshvifc(handles);

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

refreshvifc(handles);


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

if ~isempty(numat); % user did not press cancel
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
refreshvifc(handles);






% --- Executes on selection change in vilist.
function vilist_Callback(hObject, eventdata, handles)
% hObject    handle to vilist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns vilist contents as cell array
%        contents{get(hObject,'Value')} returns selected item from vilist
M=getappdata(gcf,'M');


M.ui.volumeintersections=get(handles.vilist,'Value');

% store model and refresh UI
setappdata(gcf,'M',M);
refreshvifc(handles);



% --- Executes during object creation, after setting all properties.
function vilist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vilist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in fclist.
function fclist_Callback(hObject, eventdata, handles)
% hObject    handle to fclist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns fclist contents as cell array
%        contents{get(hObject,'Value')} returns selected item from fclist
M=getappdata(gcf,'M');


M.ui.fibercounts=get(handles.fclist,'Value');

% store model and refresh UI
setappdata(gcf,'M',M);
refreshvifc(handles);

% --- Executes during object creation, after setting all properties.
function fclist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fclist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function [pathname] = ea_uigetdir(start_path, dialog_title)
% Pick a directory with the Java widgets instead of uigetdir

import javax.swing.JFileChooser;

if nargin == 0 || strcmp(start_path,'') || start_path == 0 % Allow a null argument.
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



function refreshvifc(handles)
ea_busyaction('on',handles.lg_figure,'group');



% get model data

M=getappdata(handles.lg_figure,'M');

if strcmp(get(handles.groupdir_choosebox,'String'),'Choose Group Directory') % not set yet.
    ea_busyaction('off',handles.lg_figure,'group');
    return
end


% refresh group list
set(handles.grouplist,'String',M.patient.group);
if length(get(handles.patientlist,'String'))<max(M.ui.listselect)

    M.ui.listselect=1;
end
try set(handles.grouplist,'Value',M.ui.listselect);  end



% refresh patient list
set(handles.patientlist,'String',M.patient.list);
try set(handles.patientlist,'Value',M.ui.listselect); end


% refresh clinical list
set(handles.clinicallist,'String',M.clinical.labels);
try set(handles.clinicallist,'Value',M.ui.clinicallist); end


if get(handles.clinicallist,'Value')>length(get(handles.clinicallist,'String'))
    set(handles.clinicallist,'Value',length(get(handles.clinicallist,'String')));
end


% set isomatrix from variable in clinical list
try
    M.isomatrix=M.clinical.vars{get(handles.clinicallist,'Value')};
    M.isomatrix_name=M.clinical.labels{get(handles.clinicallist,'Value')};
end

% refresh selections on VI and FC Lists:
try

    if max(M.ui.volumeintersections)>length(get(handles.vilist,'String'))
        set(handles.vilist,'Value',1);
    else
        set(handles.vilist,'Value',M.ui.volumeintersections);
    end
end

try
    if M.ui.fibercounts>length(get(handles.fclist,'String'))
        set(handles.fclist,'Value',1);
    else
        set(handles.fclist,'Value',M.ui.fibercounts);
    end
end

M.groups.group=unique(M.patient.group); % STN, GPi, Thalamus, cZi
%groupcolors=squeeze(ind2rgb(round([1:9]*(64/9)),jet));


if ~isfield(M.groups,'color')
    M.groups.color=repmat(0.2,length(M.groups.group),3);
end

if ~isequal(size(M.groups.color),[length(M.groups.group),3])
    M.groups.color=repmat(0.2,length(M.groups.group),3);
end


% add graph metrics to connectome graph-metrics popup:


thisparc=get(handles.labelpopup,'String');
thisparc=thisparc{get(handles.labelpopup,'Value')};
try
    gmdir=dir([M.patient.list{1},filesep,'connectomics',filesep,thisparc,filesep,'graph',filesep,'*.nii']);

    gms{1}='';
    for gm=1:length(gmdir)
        gms{gm}=gmdir(gm).name;
    end
    set(handles.lc_graphmetric,'String',gms);
end




%% modalities for VAT metrics:

% dMRI:
cnt=1;
options.prefs=ea_prefs('');
options.earoot=[ea_getearoot];
try
    directory=[M.patient.list{1},filesep];
    modlist=ea_genmodlist(directory,thisparc,options);
    modlist{end+1}='Patient-specific fiber tracts';
    modlist{end+1}='Do not calculate connectivity stats';
    set(handles.fiberspopup,'String',modlist);
    if get(handles.fiberspopup,'Value')>length(modlist);
        set(handles.fiberspopup,'Value',1);
    end
end


% update UI



% % make setstimparams button green if set.
% if isfield(M,'stimparams')
%     if isfield(M.stimparams,'U')
%         set(handles.setstimparamsbutton,'BackgroundColor',[0.1;0.8;0.1]);
%     else
%         set(handles.setstimparamsbutton,'BackgroundColor',[0.93,0.93,0.93]);
%     end
% else
%     set(handles.setstimparamsbutton,'BackgroundColor',[0.93,0.93,0.93]);
% end

S=getappdata(handles.lg_figure,'S');

if ~isempty(S)
    set(handles.setstimparamsbutton,'BackgroundColor',[0.1;0.8;0.1]);
    M.S=S;
    M.S=ea_activecontacts(M.S);

    M.vatmodel=getappdata(handles.lg_figure,'vatmodel');
else
    set(handles.setstimparamsbutton,'BackgroundColor',[0.93,0.93,0.93]);
end


% make choosecolors button green if chosen.
if isfield(M.groups,'colorschosen')
    set(handles.choosegroupcolors,'BackgroundColor',[0.1;0.8;0.1]);
else
    set(handles.choosegroupcolors,'BackgroundColor',[0.93,0.93,0.93]);
end

% update checkboxes:

try set(handles.showactivecontcheck,'Value',M.ui.showactivecontcheck); end
try set(handles.showpassivecontcheck,'Value',M.ui.showpassivecontcheck); end
try set(handles.highlightactivecontcheck,'Value',M.ui.hlactivecontcheck); end
try set(handles.showisovolumecheck,'Value',M.ui.showisovolumecheck); end
try set(handles.statvatcheck,'Value',M.ui.statvat); end
try set(handles.colorpointcloudcheck,'Value',M.ui.colorpointcloudcheck); end
try set(handles.lc_smooth,'Value',M.ui.lc.smooth); end



% update selectboxes:
try set(handles.elrenderingpopup,'Value',M.ui.elrendering); end
try set(handles.atlassetpopup,'Value',M.ui.atlassetpopup); end
try set(handles.fiberspopup,'Value',M.ui.fiberspopup); end
try set(handles.labelpopup,'Value',M.ui.labelpopup); end
try set(handles.elmodelselect,'Value',M.ui.elmodelselect); end
try set(handles.normregpopup,'Value',M.ui.normregpopup); end
try set(handles.labelpopup,'Value',M.ui.lc.parcellation); end
try set(handles.lc_normalization,'Value',M.ui.lc.normalization); end
try set(handles.lc_graphmetric,'Value',M.ui.lc.graphmetric); end

% update enable-disable-dependencies:
try
    if M.ui.elrendering==3
        try set(handles.colorpointcloudcheck,'Enable','on'); end
    else
        try set(handles.colorpointcloudcheck,'Enable','off'); end
    end
end
% hide detachbutton if already detached:
try
    if M.ui.detached
        set(handles.detachbutton,'Visible','off');
    end
end

%% patient specific part:
if ~isempty(M.patient.list)
    
    % add modalities to NBS stats metric popup:
    
    tryparcs=dir([M.patient.list{1},filesep,'connectomics',filesep,thisparc,filesep,'*_CM.mat']);
if isempty(tryparcs)
    set(handles.lc_metric,'String','No data found.');
else
    avparcs=ones(length(tryparcs),1);
    for sub=1:length(M.patient.list)
        for parc=1:length(tryparcs)
            if ~exist([M.patient.list{1},filesep,'connectomics',filesep,thisparc,filesep,tryparcs(parc).name],'file');
                avparcs(parc)=0;
            end
        end
    end
    tryparcs=tryparcs(avparcs);
    pcell=cell(length(tryparcs),1);
    for p=1:length(pcell)
        [~,pcell{p}]=fileparts(tryparcs(p).name);
    end
    set(handles.lc_metric,'String',pcell);
end
    
    
    for pt=1:length(M.patient.list)
        % set stimparams based on values provided by user
        for side=1:2

            M.stimparams(pt,side).usefiberset=get(handles.fiberspopup,'String');
            try
                M.stimparams(pt,side).usefiberset=M.stimparams(pt,side).usefiberset{M.ui.fiberspopup};
            catch
                M.stimparams(pt,side).usefiberset=length(get(handles.fiberspopup,'String'));
                M.ui.fiberspopup=length(get(handles.fiberspopup,'String'));
            end
            M.stimparams(pt,side).labelatlas=get(handles.labelpopup,'String');
            M.stimparams(pt,side).labelatlas=M.stimparams(pt,side).labelatlas(M.ui.labelpopup);
            M.stimparams(pt,side).showfibers=1;
            M.stimparams(pt,side).fiberthresh=1;

            M.stimparams(pt,side).showconnectivities=1;
        end
        % load localization

        [~,pats{pt}]=fileparts(M.patient.list{pt});

        M.elstruct(pt).group=M.patient.group(pt);
        M.elstruct(pt).groupcolors=M.groups.color;
        M.elstruct(pt).groups=M.groups.group;

        options.sides=1:2;
        options.native=0;
        try

            [options.root,options.patientname]=fileparts(M.patient.list{pt});
            options.root=[options.root,filesep];

            [coords_mm,trajectory,markers,elmodel,manually_corrected]=ea_load_reconstruction(options);
            if M.ui.elmodelselect==1 % use patient specific elmodel
                if exist('elmodel','var')
                    M.elstruct(pt).elmodel=elmodel;
                else
                    M.elstruct(pt).elmodel='Medtronic 3389'; % use default for older reconstructions that did not store elmodel.
                end
            else
                elmodels=get(handles.elmodelselect,'String');
                M.elstruct(pt).elmodel=elmodels{get(handles.elmodelselect,'Value')};

            end
            M.elstruct(pt).coords_mm=coords_mm;
            M.elstruct(pt).trajectory=trajectory;

            M.elstruct(pt).name=[pats{pt}];
            if ~exist('markers','var') % backward compatibility to old recon format

                for side=1:2
                    markers(side).head=coords_mm{side}(1,:);
                    markers(side).tail=coords_mm{side}(4,:);
                    normtrajvector=(markers(side).tail-markers(side).head)./norm(markers(side).tail-markers(side).head);
                    orth=null(normtrajvector)*(options.elspec.lead_diameter/2);
                    markers(side).x=coords_mm{side}(1,:)+orth(:,1)';
                    markers(side).y=coords_mm{side}(1,:)+orth(:,2)'; % corresponding points in reality
                end
            end
            M.elstruct(pt).markers=markers;

        catch
            if pt>1 % first patient has worked but some other patient seems not to have worked.
                try
                    M.elstruct(1).coords_mm; % probe if error happens in pt. 1 ? if not show warning
                    warning(['No reconstruction present for ',pats{pt},'. Please check.']);
                end
            end
        end

    end



    % load stats for group

    for pt=1:length(M.patient.list)


        % (re-)load stats
        try
            load([M.patient.list{pt},filesep,'ea_stats']);
            M.stats(pt).ea_stats=ea_stats;
        end

        if ~isfield(M,'stats')
            % if no stats  present yet, return.
            setappdata(handles.lg_figure,'M',M);
            set(handles.lg_figure,'name','LEAD-DBS Group Analysis');
            break
        end

        priorvilist=M.vilist;
        try % try using stats from patient folder.
            M.vilist=ea_stats.atlases.names;
        catch
            try % try using stats from M-file.
                M.vilist=M.stats(pt).ea_stats.atlases.names;
            catch
                M.vilist={};
            end
        end
        % check and compare with prior atlas intersection list.

        if ~isempty(priorvilist) && ~isequal(priorvilist,M.vilist)

            warning('Patient stats are inhomogeneous. Please re-run group analysis (Section Prepare DBS stats).');
        end



        priorfclist=M.fclist;
        try % try using stats from patient folder.
            M.fclist=ea_stats.stimulation(1).ft(1).labels{1};
            fcdone=1;
        catch
            try % try using stats from M-file.
                M.fclist=M.stats(pt).ea_stats.stimulation(1).ft(1).labels{1};
                fcdone=1;
            catch
                M.fclist={};
                fcdone=0;
            end
        end

        % check and compare with prior fibertracking list.

        if fcdone

            if ~isempty(priorfclist) && ~isequal(priorfclist,M.fclist)
                warning('Trying to analyse inhomogeneous patient group. Please re-run single subject lead analysis with patients using always the same labeling atlas.');
            end

        end
    end

    try
        setappdata(handles.lg_figure,'elstruct',elstruct);
    end
else
    M.vilist={};
    M.fclist={};
end


% store everything in Model
setappdata(handles.lg_figure,'M',M);

% refresh UI

set(handles.vilist,'String',M.vilist);
set(handles.fclist,'String',M.fclist);



ea_busyaction('off',handles.lg_figure,'group');


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
refreshvifc(handles);




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
refreshvifc(handles);



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
refreshvifc(handles);


% --- Executes on button press in ttestbutton.
function ttestbutton_Callback(hObject, eventdata, handles)
% hObject    handle to ttestbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


stats=preparedataanalysis(handles);

assignin('base','stats',stats);

% perform t-tests:

if ~isempty(stats.vicorr.both)
    ea_ttest(stats.vicorr.both(repmat(logical(stats.corrcl),1,size(stats.vicorr.both,2))),stats.vicorr.both(~repmat(logical(stats.corrcl),1,size(stats.vicorr.both,2))),'Volume Intersections, both hemispheres',stats.vc_labels);
end
if ~isempty(stats.vicorr.right)
    ea_ttest(stats.vicorr.right(repmat(logical(stats.corrcl),1,size(stats.vicorr.right,2))),stats.vicorr.right(~repmat(logical(stats.corrcl),1,size(stats.vicorr.right,2))),'Volume Intersections, right hemisphere',stats.vc_labels);
end
if ~isempty(stats.vicorr.left)
    ea_ttest(stats.vicorr.left(repmat(logical(stats.corrcl),1,size(stats.vicorr.left,2))),stats.vicorr.left(~repmat(logical(stats.corrcl),1,size(stats.vicorr.left,2))),'Volume Intersections, left hemisphere',stats.vc_labels);
end


if ~isempty(stats.vicorr.nboth)
    ea_ttest(stats.vicorr.nboth(repmat(logical(stats.corrcl),1,size(stats.vicorr.both,2))),stats.vicorr.both(~repmat(logical(stats.corrcl),1,size(stats.vicorr.both,2))),'Normalized Volume Intersections, both hemispheres',stats.vc_labels);
end
if ~isempty(stats.vicorr.nright)
    ea_ttest(stats.vicorr.nright(repmat(logical(stats.corrcl),1,size(stats.vicorr.right,2))),stats.vicorr.right(~repmat(logical(stats.corrcl),1,size(stats.vicorr.right,2))),'Normalized Volume Intersections, right hemisphere',stats.vc_labels);
end
if ~isempty(stats.vicorr.nleft)
    ea_ttest(stats.vicorr.nleft(repmat(logical(stats.corrcl),1,size(stats.vicorr.left,2))),stats.vicorr.left(~repmat(logical(stats.corrcl),1,size(stats.vicorr.left,2))),'Normalized Volume Intersections, left hemisphere',stats.vc_labels);
end

if ~isempty(stats.fc.fccorr)
    ea_ttest(stats.fc.fccorr(repmat(logical(stats.corrcl),1,size(stats.fc.fccorr,2))),stats.fc.fccorr(~repmat(logical(stats.corrcl),1,size(stats.fc.fccorr,2))),'Fibercounts',stats.vc_labels);
end

if ~isempty(stats.fc.nfccorr)
    ea_ttest(stats.fc.nfccorr(repmat(logical(stats.corrcl),1,size(stats.fc.nfccorr,2))),stats.fc.nfccorr(~repmat(logical(stats.corrcl),1,size(stats.fc.nfccorr,2))),'Normalized Fibercounts',stats.vc_labels);
end





function [stats]=preparedataanalysis(handles)

M=getappdata(gcf,'M');


%M.stats(get(handles.vilist,'Value'))

% Get volume intersections:
vicnt=1; ptcnt=1;

howmanyvis=length(get(handles.vilist,'Value'));
howmanypts=length(get(handles.patientlist,'Value'));

vicorr_right=zeros(howmanypts,howmanyvis); vicorr_left=zeros(howmanypts,howmanyvis); vicorr_both=zeros(howmanypts,howmanyvis);
nvicorr_right=zeros(howmanypts,howmanyvis); nvicorr_left=zeros(howmanypts,howmanyvis); nvicorr_both=zeros(howmanypts,howmanyvis);
vc_labels={};

for vi=get(handles.vilist,'Value') % get volume interactions for each patient from stats
    for pt=get(handles.patientlist,'Value')
        usewhichstim=length(M.stats(pt).ea_stats.stimulation); % always use last analysis!
        for side=1:size(M.stats(pt).ea_stats.stimulation(usewhichstim).ft,2)
            for vat=1;
                if side==1 % right hemisphere
                    vicorr_right(ptcnt,vicnt)=vicorr_right(ptcnt,vicnt)+M.stats(pt).ea_stats.stimulation(usewhichstim).vat(side,vat).AtlasIntersection(vi);
                    nvicorr_right(ptcnt,vicnt)=nvicorr_right(ptcnt,vicnt)+M.stats(pt).ea_stats.stimulation(usewhichstim).vat(side,vat).nAtlasIntersection(vi);
                elseif side==2 % left hemisphere
                    vicorr_left(ptcnt,vicnt)=vicorr_left(ptcnt,vicnt)+M.stats(pt).ea_stats.stimulation(usewhichstim).vat(side,vat).AtlasIntersection(vi);
                    nvicorr_left(ptcnt,vicnt)=nvicorr_left(ptcnt,vicnt)+M.stats(pt).ea_stats.stimulation(usewhichstim).vat(side,vat).nAtlasIntersection(vi);
                end
                vicorr_both(ptcnt,vicnt)=vicorr_both(ptcnt,vicnt)+M.stats(pt).ea_stats.stimulation(usewhichstim).vat(side,vat).AtlasIntersection(vi);
                nvicorr_both(ptcnt,vicnt)=nvicorr_both(ptcnt,vicnt)+M.stats(pt).ea_stats.stimulation(usewhichstim).vat(side,vat).nAtlasIntersection(vi);

            end
        end


        % check if all three values have been served. if not, set to zero
        % (e.g. if there was no stimulation at all on one hemisphere, this
        % could happen.

        ptcnt=ptcnt+1;

    end
    vc_labels{end+1}=M.stats(pt).ea_stats.atlases.names{vi};

    ptcnt=1;
    vicnt=vicnt+1;
end


% Get fibercounts (here first ft is always right hemispheric, second always left hemispheric). There will always be two fts used.:
howmanyfcs=length(get(handles.fclist,'Value'));

fccnt=1; ptcnt=1;
fccorr_right=zeros(howmanypts,howmanyfcs);
nfccorr_right=zeros(howmanypts,howmanyfcs);
fccorr_left=zeros(howmanypts,howmanyfcs);
nfccorr_left=zeros(howmanypts,howmanyfcs);
fccorr_both=zeros(howmanypts,howmanyfcs);
nfccorr_both=zeros(howmanypts,howmanyfcs);
fc_labels={};
for fc=get(handles.fclist,'Value') % get volume interactions for each patient from stats
    for pt=get(handles.patientlist,'Value')
        usewhichstim=length(M.stats(pt).ea_stats.stimulation); % always use last analysis!
        fccorr_right(ptcnt,fccnt)=M.stats(pt).ea_stats.stimulation(usewhichstim).ft(1).fibercounts{1}(fc);
        nfccorr_right(ptcnt,fccnt)=M.stats(pt).ea_stats.stimulation(usewhichstim).ft(1).nfibercounts{1}(fc);
        fccorr_left(ptcnt,fccnt)=M.stats(pt).ea_stats.stimulation(usewhichstim).ft(2).fibercounts{1}(fc);
        nfccorr_left(ptcnt,fccnt)=M.stats(pt).ea_stats.stimulation(usewhichstim).ft(2).nfibercounts{1}(fc);
        fccorr_both(ptcnt,fccnt)=M.stats(pt).ea_stats.stimulation(usewhichstim).ft(1).fibercounts{1}(fc)+M.stats(pt).ea_stats.stimulation(usewhichstim).ft(2).fibercounts{1}(fc);
        nfccorr_both(ptcnt,fccnt)=M.stats(pt).ea_stats.stimulation(usewhichstim).ft(1).nfibercounts{1}(fc)+M.stats(pt).ea_stats.stimulation(usewhichstim).ft(2).nfibercounts{1}(fc);
        ptcnt=ptcnt+1;

    end
    ptcnt=1;
    fccnt=fccnt+1;
    fc_labels{end+1}=M.stats(pt).ea_stats.stimulation(usewhichstim).ft(1).labels{1}{fc};

end


% prepare outputs:

vicorr.both=vicorr_both;
vicorr.left=vicorr_left;
vicorr.right=vicorr_right;
vicorr.nboth=nvicorr_both;
vicorr.nleft=nvicorr_left;
vicorr.nright=nvicorr_right;
fccorr.both=fccorr_both;
fccorr.nboth=nfccorr_both;
fccorr.right=fccorr_right;
fccorr.nright=nfccorr_right;
fccorr.left=fccorr_left;
fccorr.nleft=nfccorr_left;

% clinical vector:
corrcl=M.clinical.vars{get(handles.clinicallist,'Value')};

corrcl=corrcl(get(handles.patientlist,'Value'),:);

clinstrs=get(handles.clinicallist,'String');
vc_labels=[clinstrs(get(handles.clinicallist,'Value')),vc_labels]; % add name of clinical vector to labels
fc_labels=[clinstrs(get(handles.clinicallist,'Value')),fc_labels]; % add name of clinical vector to labels


stats.corrcl=corrcl;
stats.vicorr=vicorr;
stats.fccorr=fccorr;
stats.vc_labels=vc_labels;
stats.fc_labels=fc_labels;









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

refreshvifc(handles);


% --- Executes on button press in moveptdownbutton.
function moveptdownbutton_Callback(hObject, eventdata, handles)
% hObject    handle to moveptdownbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
M=getappdata(gcf,'M');
whichmoved=get(handles.patientlist,'Value');
if length(whichmoved)>1; return; end % more than one selected..
if whichmoved==1 % first entry anyways
    return
end
ix=1:length(M.patient.list);
ix(whichmoved)=ix(whichmoved)-1;
ix(whichmoved-1)=ix(whichmoved-1)+1;

M.patient.list=M.patient.list(ix);
M.patient.group=M.patient.group(ix);
try
    M=rmfield(M,'stats');
end
try
    M=rmfield(M,'elstruct');
end
setappdata(gcf,'M',M);
set(handles.patientlist,'Value',whichmoved-1);

refreshvifc(handles);

% --- Executes on button press in moveptupbutton.
function moveptupbutton_Callback(hObject, eventdata, handles)
% hObject    handle to moveptupbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

M=getappdata(gcf,'M');
whichmoved=get(handles.patientlist,'Value');
if length(whichmoved)>1; return; end % more than one selected..
if whichmoved==length(M.patient.list) % last entry anyways
    return
end
ix=1:length(M.patient.list);
ix(whichmoved)=ix(whichmoved)+1;
ix(whichmoved+1)=ix(whichmoved+1)-1;

M.patient.list=M.patient.list(ix);
M.patient.group=M.patient.group(ix);
try
    M=rmfield(M,'stats');
end
try
    M=rmfield(M,'elstruct');
end
setappdata(gcf,'M',M);
set(handles.patientlist,'Value',whichmoved+1);

refreshvifc(handles);

% --- Executes on button press in calculatebutton.
function calculatebutton_Callback(hObject, eventdata, handles)
% hObject    handle to calculatebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

refreshvifc(handles);

M=getappdata(gcf,'M');

% set options
options=ea_setopts_local(handles);
stimname=ea_detstimname();

    % determine if fMRI or dMRI
    mods=get(handles.fiberspopup,'String');
    mod=mods{get(handles.fiberspopup,'Value')};
    switch mod
        case {'Patient-specific fiber tracts','rest_tc'}
            fibersfile=mod;
        case 'Do not calculate connectivity stats'
        otherwise % load fibertracts once and for all subs here.
                [fibersfile.fibers,fibersfile.fibersidx]=ea_loadfibertracts([ea_getearoot,'fibers',filesep,mod,'.mat']);
    end


for pt=M.ui.listselect

    % set pt specific options

    % own fileparts to support windows/mac/linux slashes even if they come
    % from a different OS.
    if isempty(strfind(M.patient.list{pt},'/'))
        lookfor='\';
    else
        lookfor='/';
    end

    slashes=strfind(M.patient.list{pt},lookfor);
    if ~isempty(slashes)
        options.patientname=M.patient.list{pt}(slashes(end)+1:end);
        options.root=M.patient.list{pt}(1:slashes(end));

    else
        options.patientname=M.patient.list{pt};
        options.root='';
    end

    disp(['Processing ',options.patientname,'.']);
try
    options.numcontacts=size(M.elstruct(pt).coords_mm{1},1);
catch % no localization present or in wrong format.
    ea_error(['Please localize ',options.patientname,' first.']);
end
    options.elmodel=M.elstruct(pt).elmodel;
    options=ea_resolve_elspec(options);
    options.prefs=ea_prefs(options.patientname);
    options.d3.verbose='off';

    options.native=0;


    options.d3.elrendering=M.ui.elrendering;
    options.d3.hlactivecontacts=get(handles.highlightactivecontcheck,'Value');
    options.d3.showactivecontacts=get(handles.showactivecontcheck,'Value');
    options.d3.showpassivecontacts=get(handles.showpassivecontcheck,'Value');
    try options.d3.isomatrix=M.isomatrix; end

    options.d3.isovscloud=M.ui.isovscloudpopup;
    options.d3.showisovolume=M.ui.showisovolumecheck;

    options.expstatvat.do=0;
    try
        options.expstatvat.vars=M.clinical.vars(M.ui.clinicallist);
        options.expstatvat.labels=M.clinical.labels(M.ui.clinicallist);
        options.expstatvat.pt=pt;
    end
    options.expstatvat.dir=M.ui.groupdir;
    processlocal=0;

        if M.ui.detached
            processlocal=1;
            mkdir([M.ui.groupdir,options.patientname]);
            options.root=M.ui.groupdir;
            %    options.patientname='tmp';

            ea_stats=M.stats(pt).ea_stats;
            coords_mm=M.elstruct(pt).coords_mm;
            trajectory=M.elstruct(pt).trajectory;
            save([M.ui.groupdir,options.patientname,filesep,'ea_stats'],'ea_stats');
            save([M.ui.groupdir,options.patientname,filesep,'ea_reconstruction'],'coords_mm','trajectory');
        end

        if ~exist(options.root,'file') % data is not there. Act as if detached. Process in tmp-dir.
            processlocal=1;
            warning('on');
            warning('Data has been detached from group-directory. Will process locally. Please be aware that you might loose this newly-processed data once you re-attach the single-patient data to the analysis!');
            warning('off');
            mkdir([M.ui.groupdir,options.patientname]);
            options.root=M.ui.groupdir;
            % options.patientname='tmp';

            ea_stats=M.stats(pt).ea_stats;
            coords_mm=M.elstruct(pt).coords_mm;
            trajectory=M.elstruct(pt).trajectory;
            save([M.ui.groupdir,options.patientname,filesep,'ea_stats'],'ea_stats');
            save([M.ui.groupdir,options.patientname,filesep,'ea_reconstruction'],'coords_mm','trajectory');
        end


    %delete([options.root,options.patientname,filesep,'ea_stats.mat']);

    % Step 1: Re-calculate closeness to subcortical atlases.

    resultfig=ea_elvis(options);
    % save scene as matlab figure


    options.modality=ea_checkctmrpresent(M.patient.list{pt});
    if options.modality(1) % prefer MR
        options.modality=1;
    else
        if options.modality(2)
            options.modality=2;
        else
            options.modality=1;
            warning(['No MR or CT volumes found in ',M.patient.list{pt},'.']);
        end
    end




    % Step 2: Re-calculate VAT
    if isfield(M,'S')
        try
            setappdata(resultfig,'S',M.S(pt));
        catch
            ea_error(['Stimulation parameters for ',M.patient.list{pt},' are missing.']);
        end
        vfnames=getappdata(handles.lg_figure,'vatfunctionnames');

        [~,ix]=ismember(M.vatmodel,vfnames);
        vfs=getappdata(handles.lg_figure,'genvatfunctions');

        ea_genvat=eval(['@',vfs{ix}]);

        for side=1:2
            setappdata(resultfig,'elstruct',M.elstruct(pt));
            setappdata(resultfig,'elspec',options.elspec);
            [stimparams(1,side).VAT(1).VAT,volume]=feval(ea_genvat,M.elstruct(pt).coords_mm,M.S(pt),side,options,stimname);
            stimparams(1,side).volume=volume;
        end

        setappdata(resultfig,'stimparams',stimparams(1,:));
    end
    % Step 3: Re-calculate connectivity from VAT to rest of the brain.
    if ~strcmp(mod,'Do not calculate connectivity stats')

        % Convis part:




        parcs=get(handles.labelpopup,'String');
        selectedparc=parcs{get(handles.labelpopup,'Value')};
        directory=[options.root,options.patientname,filesep];
        if ischar(fibersfile)
            switch mod
                case 'rest_tc'
                    ea_error('Group statistics for fMRI are not yet supported. Sorry, check back later!');
                    pV=spm_vol([options.earoot,'templates',filesep,'labeling',filesep,selectedparc,'.nii']);
                    pX=spm_read_vols(pV);
                    ea_cvshowvatfmri(resultfig,pX,directory,filesare,handles,pV,selectedparc,options);
                otherwise
                    ea_cvshowvatdmri(resultfig,directory,{mod,stimname},selectedparc,options);
            end
        else
            ea_cvshowvatdmri(resultfig,directory,{fibersfile,stimname},selectedparc,options);
        end
    end
    close(resultfig);

    if processlocal % gather stats and recos to M
        load([M.ui.groupdir,options.patientname,filesep,'ea_stats']);
        load([M.ui.groupdir,options.patientname,filesep,'ea_reconstruction']);

        M.stats(pt).ea_stats=ea_stats;
        M.elstruct(pt).coords_mm=coords_mm;
        M.elstruct(pt).trajectory=trajectory;
        setappdata(gcf,'M',M);

        save([M.ui.groupdir,'LEAD_groupanalysis.mat'],'M','-v7.3');
        try      movefile([options.root,options.patientname,filesep,'LEAD_scene.fig'],[M.ui.groupdir,'LEAD_scene_',num2str(pt),'.fig']); end
        %rmdir([M.ui.groupdir,'tmp'],'s');


    end
end
%% processing done here.

% calculate group-results for expstatvat if required:
% if options.expstatvat.do
%     disp('Averaging VAT-Stat files to a group result.');
%     for side=1:2
%         switch side
%             case 1
%                 si='rh';
%             case 2
%                 si='lh';
%         end
%         fis=dir([M.ui.groupdir,'statvat_results',filesep,'*_',si,'.nii']);
%         ea_dispercent(0,['Reading in files, ',si]);
%         for fi=1:length(fis);
%             ea_dispercent(fi/length(fis));
%             vV=spm_vol([M.ui.groupdir,'statvat_results',filesep,fis(fi).name]);
%             if ~exist('thisvat','var')
%                 tmp=spm_read_vols(vV);
%                 thisvat=zeros([size(tmp),length(fis)]);
%                 thisvat(:,:,:,1)=tmp;
%                 clear tmp
%             else
%                 try
%                 thisvat(:,:,:,fi)=spm_read_vols(vV);
%                 catch
%                     keyboard
%                 end
%             end
%         end
%         ea_dispercent(100,'end');
%         disp('Averaging...');
%         thisvat(thisvat==0)=nan;
%
%         thisvat=nanmean(thisvat,4);
%         disp('Done.');
%         disp('Saving file...');
%         %thisvat=thisvat/fi; % simple mean here.
%         vV.fname=[M.ui.groupdir,'statvat_results',filesep,si,'_mean.nii'];
%         vV.dt=[64,0];
%         spm_write_vol(vV,thisvat);
%         disp('Done.');
%     end
% end

refreshvifc(handles);





% --- Executes on selection change in fiberspopup.
function fiberspopup_Callback(hObject, eventdata, handles)
% hObject    handle to fiberspopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns fiberspopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from fiberspopup
M=getappdata(gcf,'M');
M.ui.fiberspopup=get(handles.fiberspopup,'Value');
setappdata(gcf,'M',M);
refreshvifc(handles);

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
M.ui.labelpopup=get(handles.labelpopup,'Value');
setappdata(gcf,'M',M);
refreshvifc(handles);

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
M.ui.atlassetpopup=get(handles.atlassetpopup,'Value');
setappdata(gcf,'M',M);
refreshvifc(handles);

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

options.earoot=[ea_getearoot];
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
options.writeoutpm=1;
options.colormap=jet;
options.d3.write=1;
options.d3.prolong_electrode=2;
options.d3.writeatlases=1;
options.macaquemodus=0;


% --- Executes on button press in groupdir_choosebox.
function groupdir_choosebox_Callback(hObject, eventdata, handles)
% hObject    handle to groupdir_choosebox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

nudir=[uigetdir];

if ~nudir % user pressed cancel
    return
end
ea_busyaction('on',handles.lg_figure,'group');

nudir=[nudir,filesep];
M=initializeM;


set(handles.groupdir_choosebox,'String',nudir);

try % if file already exists, load it (and overwrite M).
    load([nudir,'LEAD_groupanalysis.mat']);
catch % if not, store it saving M.
    save([nudir,'LEAD_groupanalysis.mat'],'M','-v7.3');
end

M.ui.groupdir=nudir;
setappdata(handles.lg_figure,'M',M);
try
setappdata(handles.lg_figure,'S',M.S);
setappdata(handles.lg_figure,'vatmodel',M.S(1).model);
end

ea_busyaction('off',handles.lg_figure,'group');

refreshvifc(handles);



% --- Executes on button press in opensubgui.
function opensubgui_Callback(hObject, eventdata, handles)
% hObject    handle to opensubgui (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
M=getappdata(gcf,'M');
lead('loadsubs',M.patient.list(M.ui.listselect));


% --- Executes on button press in choosegroupcolors.
function choosegroupcolors_Callback(hObject, eventdata, handles)
% hObject    handle to choosegroupcolors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


M=getappdata(gcf,'M');

for g=unique(M.patient.group)'
    M.groups.color(ismember(M.groups.group,g),:)=...
        uisetcolor(M.groups.color(ismember(M.groups.group,g),:),['Group ',num2str(g),':']);
end
M.groups.colorschosen=1;

setappdata(gcf,'M',M);
refreshvifc(handles);

function M=initializeM
M=struct;
M.patient.list={};
M.patient.group=[];

M.clinical.vars={};
M.clinical.labels={};
M.vilist={};
M.fclist={};
M.ui=struct;
M.ui.listselect=1;
M.ui.elrendeting=1;
M.ui.clinicallist=1;
M.ui.hlactivecontcheck=0;
M.ui.showpassivecontcheck=1;
M.ui.showactivecontcheck=1;
M.ui.showisovolumecheck=0;
M.ui.isovscloudpopup=1;
M.ui.atlassetpopup=1;
M.ui.fiberspopup=3;
M.ui.labelpopup=1;
M.ui.volumeintersections=1;
M.ui.fibercounts=1;
M.ui.elrendering=1;
M.ui.statvat=0;
M.ui.elmodelselect=1;
M.ui.detached=0;
M.ui.normregpopup=1;
M.ui.colorpointcloudcheck=0;
M.ui.lc.parcellation=1;
M.ui.lc.graphmetric=1;
M.ui.lc.normalization=1;
M.ui.lc.smooth=1;

M.S=[];
M.vatmodel=[];


% --- Executes on button press in setstimparamsbutton.
function setstimparamsbutton_Callback(hObject, eventdata, handles)
% hObject    handle to setstimparamsbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
M=getappdata(gcf,'M');

% try
%     uicell=inputdlg('Enter Variable name for Voltage-Parameters','Enter Stimulation Settings...',1);
% uidata.U=evalin('base',uicell{1});
% catch
%     warning('Stim-Params could not be evaluated. Please Try again.');
%     return
% end
% try
%         uicell=inputdlg('Enter Variable name for Impedance-Parameters','Enter Stimulation Settings...',1);
% uidata.Im=evalin('base',uicell{1});
% catch
%     warning('Stim-Params could not be evaluated. Please Try again.');
%     return
% end

options=ea_setopts_local(handles);
refreshvifc(handles);

ea_stimparams(M.elstruct,handles.lg_figure,options);



% --- Executes on button press in highlightactivecontcheck.
function highlightactivecontcheck_Callback(hObject, eventdata, handles)
% hObject    handle to highlightactivecontcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of highlightactivecontcheck
M=getappdata(gcf,'M');
M.ui.hlactivecontcheck=get(handles.highlightactivecontcheck,'Value');


setappdata(gcf,'M',M);
refreshvifc(handles);



% --- Executes on selection change in elrenderingpopup.
function elrenderingpopup_Callback(hObject, eventdata, handles)
% hObject    handle to elrenderingpopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns elrenderingpopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from elrenderingpopup

M=getappdata(gcf,'M');
M.ui.elrendering=get(handles.elrenderingpopup,'Value');


setappdata(gcf,'M',M);
refreshvifc(handles);


% --- Executes during object creation, after setting all properties.
function elrenderingpopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to elrenderingpopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in showpassivecontcheck.
function showpassivecontcheck_Callback(hObject, eventdata, handles)
% hObject    handle to showpassivecontcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showpassivecontcheck

M=getappdata(gcf,'M');
M.ui.showpassivecontcheck=get(handles.showpassivecontcheck,'Value');


setappdata(gcf,'M',M);
refreshvifc(handles);


% --- Executes on button press in showactivecontcheck.
function showactivecontcheck_Callback(hObject, eventdata, handles)
% hObject    handle to showactivecontcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showactivecontcheck

M=getappdata(gcf,'M');
M.ui.showactivecontcheck=get(handles.showactivecontcheck,'Value');


setappdata(gcf,'M',M);
refreshvifc(handles);


% --- Executes on button press in showisovolumecheck.
function showisovolumecheck_Callback(hObject, eventdata, handles)
% hObject    handle to showisovolumecheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showisovolumecheck
M=getappdata(gcf,'M');
M.ui.showisovolumecheck=get(handles.showisovolumecheck,'Value');

setappdata(gcf,'M',M);
refreshvifc(handles);




% --- Executes on selection change in isovscloudpopup.
function isovscloudpopup_Callback(hObject, eventdata, handles)
% hObject    handle to isovscloudpopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns isovscloudpopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from isovscloudpopup
M=getappdata(gcf,'M');
M.ui.isovscloudpopup=get(handles.isovscloudpopup,'Value');
setappdata(gcf,'M',M);
refreshvifc(handles);

% --- Executes during object creation, after setting all properties.
function isovscloudpopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to isovscloudpopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function [pU,pIm]=getstimparams(M)
try
    for pt=1:length(M.patient.list)
        for side=1:2
            pU{side}(pt,:)=M.stimparams(pt,side).U;
            pIm{side}(pt,:)=M.stimparams(pt,side).Im;
        end
    end
catch
    pU{1}=zeros(length(M.patient.list),size(M.elstruct(1).coords_mm{1},1));
    pIm{1}=repmat(1000,length(M.patient.list),size(M.elstruct(1).coords_mm{1},1));
    pU{2}=zeros(length(M.patient.list),size(M.elstruct(1).coords_mm{1},1));
    pIm{2}=repmat(1000,length(M.patient.list),size(M.elstruct(1).coords_mm{1},1));
end





% --- Executes on button press in statvatcheck.
function statvatcheck_Callback(hObject, eventdata, handles)
% hObject    handle to statvatcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of statvatcheck
M=getappdata(gcf,'M');

M.ui.statvat=get(handles.statvatcheck,'Value');
setappdata(gcf,'M',M);








% --- Executes on selection change in elmodelselect.
function elmodelselect_Callback(hObject, eventdata, handles)
% hObject    handle to elmodelselect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns elmodelselect contents as cell array
%        contents{get(hObject,'Value')} returns selected item from elmodelselect
M=getappdata(gcf,'M');
M.ui.elmodelselect=get(handles.elmodelselect,'Value');
setappdata(gcf,'M',M);
refreshvifc(handles);

% --- Executes during object creation, after setting all properties.
function elmodelselect_CreateFcn(hObject, eventdata, handles)
% hObject    handle to elmodelselect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in detachbutton.
function detachbutton_Callback(hObject, eventdata, handles)
% hObject    handle to detachbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
choice = questdlg('Would you really like to detach the group data from the single-patient data? This means that changes to single-patient reconstructions will not be updated into the group analysis anymore. This should only be done once all patients have been finally localized and an analysis needs to be fixed (e.g. after publication or when working in collaborations). Please be aware that this step cannot be undone!', ...
    'Detach Group data from single patient data...', ...
    'No, abort.','Yes, sure!','No, abort.');
% Handle response
switch choice
    case 'No, abort.'
        return
    case 'Yes, sure!'


        M=getappdata(gcf,'M');

        for pt=1:length(M.patient.list)

            slashes=findstr('/',M.patient.list{pt});
            if isempty(slashes)
                slashes=findstr('\',M.patient.list{pt});
            end
            ptname=M.patient.list{pt}(max(slashes)+1:end);


            M.patient.list{pt}=ptname;

        end
        M.ui.detached=1;

end

setappdata(gcf,'M',M);
refreshvifc(handles);


% --- Executes on button press in colorpointcloudcheck.
function colorpointcloudcheck_Callback(hObject, eventdata, handles)
% hObject    handle to colorpointcloudcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of colorpointcloudcheck
M=getappdata(gcf,'M');
M.ui.colorpointcloudcheck=get(handles.colorpointcloudcheck,'Value');
setappdata(gcf,'M',M);
refreshvifc(handles);

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
refreshvifc(handles);

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


% --- Executes when user attempts to close lg_figure.
function lg_figure_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to lg_figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure

ea_busyaction('on',gcf,'group');
if ~strcmp(get(handles.groupdir_choosebox,'String'),'Choose Group Directory') % group dir still not chosen
    disp('Saving data...');
    % save M
    refreshvifc(handles);
    M=getappdata(hObject,'M');
    try
        save([get(handles.groupdir_choosebox,'String'),'LEAD_groupanalysis.mat'],'M','-v7.3');
    catch
        warning('Data could not be saved.');
        keyboard
    end
    disp('Bye for now.');
end
ea_busyaction('off',gcf,'group');
delete(hObject);


% --- Executes on button press in targetreport.
function targetreport_Callback(hObject, eventdata, handles)
% hObject    handle to targetreport (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
refreshvifc(handles);
M=getappdata(gcf,'M');
ea_gentargetreport(M);


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
options.root=[fileparts(fileparts(get(handles.groupdir_choosebox,'String'))),filesep];
[~,options.patientname]=fileparts(fileparts(get(handles.groupdir_choosebox,'String')));


options.numcontacts=size(M.elstruct(1).coords_mm{1},1);
options.elmodel=M.elstruct(1).elmodel;
options=ea_resolve_elspec(options);
options.prefs=ea_prefs(options.patientname);
options.d3.verbose='on';

options.d3.elrendering=M.ui.elrendering;
options.d3.hlactivecontacts=get(handles.highlightactivecontcheck,'Value');
options.d3.showactivecontacts=get(handles.showactivecontcheck,'Value');
options.d3.showpassivecontacts=get(handles.showpassivecontcheck,'Value');
try options.d3.isomatrix=M.isomatrix; end
try options.d3.isomatrix_name=M.isomatrix_name; end

options.expstatvat.do=M.ui.statvat;

options.d2.showlegend=0;


options.d3.isovscloud=M.ui.isovscloudpopup;
options.d3.showisovolume=M.ui.showisovolumecheck;
options.d3.colorpointcloud=M.ui.colorpointcloudcheck;
options.normregressor=M.ui.normregpopup;

options.d2.write=1;

options.d2.atlasopacity=0.15;

options.modality=3; % use template image


if strcmp(options.atlasset,'Use none');
    options.d2.writeatlases=1;
else
    options.d2.writeatlases=1;
end

% Prior Results are loaded here inside the function (this way, function
% can be called just by giving the patient directory.

% Prepare isomatrix (includes a normalization step if M.ui.normregpopup
% says so:

try options.d3.isomatrix=ea_reformat_isomatrix(options.d3.isomatrix,M,options); end

if ~strcmp(get(handles.groupdir_choosebox,'String'),'Choose Group Directory') % group dir still not chosen
    refreshvifc(handles);
    disp('Saving data...');
    % save M
    save([get(handles.groupdir_choosebox,'String'),'LEAD_groupanalysis.mat'],'M','-v7.3');
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

for pt=1:length(M.patient.list)
    for side=1:2
        M.elstruct(pt).activecontacts{side}=M.S(pt).activecontacts{side};
    end
end
cuts=ea_writeplanes(options,M.elstruct(get(handles.patientlist,'Value')));

ea_busyaction('off',gcf,'group');



% --- Executes on button press in lc_SPM.
function lc_SPM_Callback(hObject, eventdata, handles)
% hObject    handle to lc_SPM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

spmdir=[M.ui.groupdir,'connectomics',filesep,get(handles.labelpopup,'String'),filesep,'graph',filesep,gecs,filesep,'SPM'];
rmdir(spmdir,'s');
mkdir([M.ui.groupdir,'connectomics']);
mkdir([M.ui.groupdir,'connectomics',filesep,get(handles.labelpopup,'String')]);
mkdir([M.ui.groupdir,'connectomics',filesep,get(handles.labelpopup,'String'),filesep,'graph']);
gecs=get(handles.lc_graphmetrics,'String');
[~,gecs]=gecs{M.ui.lc.graphmetrics};
mkdir([M.ui.groupdir,'connectomics',filesep,get(handles.labelpopup,'String'),filesep,'graph',filesep,gecs]);
mkdir(spmdir);


for sub=1:length(M.patient.list)
    if M.ui.lc.normalize==1;
        zzz='';
    elseif M.ui.lc.normalize==2;
        zzz='z';
        ea_histnormalize([M.patient.list{sub},'connectomics',filesep,get(handles.labelpopup,'String'),filesep,'graph',filesep,gecs,'.nii,1'],2);
    elseif M.ui.lc.normalize==3;
        zzz='k';
        ea_histnormalize([M.patient.list{sub},'connectomics',filesep,get(handles.labelpopup,'String'),filesep,'graph',filesep,gecs,'.nii,1'],3);
    end
    if M.ui.lc.smooth
        sss='s';
        ea_smooth([M.patient.list{sub},'connectomics',filesep,get(handles.labelpopup,'String'),filesep,'graph',filesep,zzz,gecs,'.nii,1']);
    else
        sss='';
    end
    fis{sub}=[M.patient.list{sub},'connectomics',filesep,get(handles.labelpopup,'String'),filesep,'graph',filesep,sss,zzz,gecs,'.nii,1'];
end

%% model specification:
matlabbatch{1}.spm.stats.factorial_design.dir = {spmdir};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = fis';
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
spm_jobman('run',{matlabbatch});
clear matlabbatch

%% model estimation:
matlabbatch{1}.spm.stats.fmri_est.spmmat = {spmdir};
matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
spm_jobman('run',{matlabbatch});
clear matlabbatch

%% contrast manager:
matlabbatch{1}.spm.stats.con.spmmat = {spmdir};
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'main effect';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;
spm_jobman('run',{matlabbatch});
clear matlabbatch


% --- Executes on selection change in lc_normalization.
function lc_normalization_Callback(hObject, eventdata, handles)
% hObject    handle to lc_normalization (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns lc_normalization contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lc_normalization
M=getappdata(gcf,'M');
M.ui.lc.normalization=get(handles.lc_normalization,'Value');
setappdata(gcf,'M',M);

% --- Executes during object creation, after setting all properties.
function lc_normalization_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lc_normalization (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in lc_graphmetric.
function lc_graphmetric_Callback(hObject, eventdata, handles)
% hObject    handle to lc_graphmetric (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns lc_graphmetric contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lc_graphmetric
M=getappdata(gcf,'M');
M.ui.lc.graphmetric=get(handles.lc_graphmetric,'Value');
setappdata(gcf,'M',M);

% --- Executes during object creation, after setting all properties.
function lc_graphmetric_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lc_graphmetric (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





% --- Executes on button press in lc_smooth.
function lc_smooth_Callback(hObject, eventdata, handles)
% hObject    handle to lc_smooth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of lc_smooth
M=getappdata(gcf,'M');
M.ui.lc.smooth=get(handles.lc_smooth,'Value');
setappdata(gcf,'M',M);



function ea_histnormalize(fname,zzz)
nii=ea_load_nii(fname);
[pth,fn,ext]=fileparts(fname);
vals=nii.img(~isnan(nii.img));
switch zzz
    case 2 % zscore
        nii.fname=[pth,'z',fn,ext];
        vals=zscore(vals(:));
    case 3 % albada
        nii.fname=[pth,'k',fn,ext];
        vals=ea_normal(vals(:));
end
nii.img(~isnan(nii.img))=vals;
spm_write_vol(nii,nii.img);

function ea_smooth(fname)
[pth,fn,ext]=fileparts(fname);

spm_smooth(fname,[pth,fn,'s',ext],[8 8 8]);





% --- Executes on button press in specify2doptions.
function specify2doptions_Callback(hObject, eventdata, handles)
% hObject    handle to specify2doptions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ea_spec2dwrite;


% --- Executes on button press in calcgroupconnectome.
function calcgroupconnectome_Callback(hObject, eventdata, handles)
% hObject    handle to calcgroupconnectome (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

options.prefs=ea_prefs('tmp');
M=getappdata(gcf,'M');

normalized_fibers_mm={}; % combined connectome
ea_dispercent(0,'Concatenating connectome');
for sub=1:length(M.patient.list)
    ea_dispercent(sub/length(M.patient.list));
    fs=load([M.patient.list{sub},filesep,options.prefs.FTR_normalized]);
    if isfield(fs,'normalized_fibers_mm')
        nfibs=fs.normalized_fibers_mm;
    else
        fn = fieldnames(fs);
        eval(sprintf('nfibs = fs.%s;',fn{1}));
    end
    if size(nfibs,1)>size(nfibs,2)
        nfibs=nfibs';
    end
    nfibs=nfibs(1:20000);
    normalized_fibers_mm=[normalized_fibers_mm,nfibs];
end
ea_dispercent(1,'end');
save([M.ui.groupdir,options.prefs.FTR_normalized],'normalized_fibers_mm','-v7.3');
load([M.ui.groupdir,options.prefs.FTR_normalized]);

% export to trackvis

% convert to vox format
options.root=[fileparts(fileparts(M.ui.groupdir)),filesep];
[~,options.patientname]=fileparts(fileparts(M.ui.groupdir));
options.earoot=[ea_getearoot];
specs.origin=[0,0,0];

nii=ea_load_nii([options.earoot,'templates',filesep,'mni_hires.nii']);
specs.dim=nii.dim;
specs.affine=nii.mat;
ea_dispercent(0,'Converting fibers to voxel format');
for fib=1:length(normalized_fibers_mm);
    ea_dispercent(fib/length(normalized_fibers_mm));
    normalized_fibers_mm{fib}=[normalized_fibers_mm{fib},ones(size(normalized_fibers_mm{fib},1),1)]';
    normalized_fibers_mm{fib}=nii.mat\normalized_fibers_mm{fib};
    normalized_fibers_mm{fib}=normalized_fibers_mm{fib}(1:3,:)';
end
ea_dispercent(1,'end');
normalized_fibers_vox=normalized_fibers_mm;
clear normalized_fibers_mm


ea_ftr2trk({normalized_fibers_vox,options.prefs.FTR_normalized},M.ui.groupdir,specs,options); % export normalized ftr to .trk

options.prefs=ea_prefs(M.ui.groupdir);



function lc_contrast_Callback(hObject, eventdata, handles)
% hObject    handle to lc_contrast (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lc_contrast as text
%        str2double(get(hObject,'String')) returns contents of lc_contrast as a double


% --- Executes during object creation, after setting all properties.
function lc_contrast_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lc_contrast (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in lc_stattest.
function lc_stattest_Callback(hObject, eventdata, handles)
% hObject    handle to lc_stattest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns lc_stattest contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lc_stattest


% --- Executes during object creation, after setting all properties.
function lc_stattest_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lc_stattest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in lc_metric.
function lc_metric_Callback(hObject, eventdata, handles)
% hObject    handle to lc_metric (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns lc_metric contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lc_metric


% --- Executes during object creation, after setting all properties.
function lc_metric_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lc_metric (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lc_threshold_Callback(hObject, eventdata, handles)
% hObject    handle to lc_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lc_threshold as text
%        str2double(get(hObject,'String')) returns contents of lc_threshold as a double


% --- Executes during object creation, after setting all properties.
function lc_threshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lc_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton37.
function pushbutton37_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton37 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton38.
function pushbutton38_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton38 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in lc_nbs.
function lc_nbs_Callback(hObject, eventdata, handles)
% hObject    handle to lc_nbs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

clearvars -global nbs
global nbs

earoot=getappdata(handles.lg_figure,'earoot');
UI.method.ui = 'Run NBS';
UI.test.ui = get(handles.lc_stattest,'String');
UI.test.ui=UI.test.ui{get(handles.lc_stattest,'Value')};
UI.thresh.ui = get(handles.lc_threshold,'String');
UI.contrast.ui = get(handles.lc_contrast,'String');
ea_preparenbs(handles)
root=get(handles.groupdir_choosebox,'String');
UI.design.ui = [root,'NBSdesignMatrix.mat'];
UI.matrices.ui = [root,'NBSdataMatrix.mat'];
UI.node_coor.ui = '';
UI.node_label.ui = '';

% get advanced options:
try
    lc=load([earoot,'connectomics',filesep,'lc_options.mat']);
catch
    lc=ea_initlcopts([]);
end
if ~isfield(lc,'nbs') % compatibility with older stored userdata (<v1.4.9)
    lc=ea_initlcopts([],lc); % will merely add the nbs stuff
end
save([earoot,'connectomics',filesep,'lc_options.mat'],'-struct','lc');
switch lc.nbs.adv.compsize
    case 1
        UI.size.ui = 'Extent';
    case 2
        UI.size.ui = 'Intensity';
end

UI.perms.ui = num2str(lc.nbs.adv.perm);
UI.alpha.ui = num2str(lc.nbs.adv.alpha);
UI.exchange.ui = lc.nbs.adv.exch;

ea_NBSrun(UI,[]);

load([root,'NBSdataMatrix.mat']);
load([root,'NBSdesignMatrix.mat']);

switch UI.test.ui
    case 't-test'
        c=eval(UI.contrast.ui);
        if length(c)>2
            ea_error('Only two-sample t-tests are fully supported at present.');
        end
        [~,~,~,tstat]=ttest2(permute(allX(:,:,logical(mX(:,find(c==1)))),[3,1,2]),...
            permute(allX(:,:,logical(mX(:,find(c==-1)))),[3,1,2]));

        T=squeeze(tstat.tstat);
        clear tstat
        clear pmask
        pmask=zeros(size(T));
        for network=1:nbs.NBS.n;
            X=full(nbs.NBS.con_mat{network});
            X=X+X';
            pmask(network,:,:)=X;
            save([root,'sig_',num2str(network)],'X');
        end
        pmask=squeeze(sum(pmask,1));
        T(~pmask)=nan;
    otherwise

        for network=1:nbs.NBS.n;
            X=full(nbs.NBS.con_mat{network});
            X=X+X';
            save([root,'sig_',num2str(network)],'X');
        end
        ea_error('Only t-test fully supported at present. Please process results manually for other tests.');
end

thisparc=get(handles.labelpopup,'String');
thisparc=thisparc{get(handles.labelpopup,'Value')};
thismetr=get(handles.lc_metric,'String');
thismetr=thismetr{get(handles.lc_metric,'Value')};
eval([thismetr,'=T;']);
expfn=[root,thisparc,'_',thismetr,'p<',UI.alpha.ui,'.mat'];
save(expfn,thismetr);

disp('** NBS done.');
if network
    disp('NBS found at least one significant network.');
    disp(['It has been stored in: ',expfn,'.']);

    disp(['To display the network(s), please load this file in the 3D viewer''s "Connectivity Visualization" under the "Matrix Level" panel.']);
else
    disp('NBS found no significant networks.');
end


function ea_preparenbs(handles)

% prepare designmatrix:

gstr=(get(handles.grouplist,'String'));

for pt=1:length(gstr);
   gv(pt)=str2double(gstr(pt));
end
mX=zeros(length(gv),max(gv));
for g=1:max(gv)
   mX(:,g)=gv==g;
end
save([get(handles.groupdir_choosebox,'String'),'NBSdesignMatrix'],'mX');

% prepare data matrix:

M=getappdata(handles.lg_figure,'M');
thisparc=get(handles.labelpopup,'String');
thisparc=thisparc{get(handles.labelpopup,'Value')};
thismetr=get(handles.lc_metric,'String');
thismetr=thismetr{get(handles.lc_metric,'Value')};

for pt=1:length(M.patient.list)
    X=load([M.patient.list{pt},filesep,'connectomics',filesep,thisparc,filesep,thismetr,'.mat']);
    fn=fieldnames(X);
    if ~exist('allX','var')
       allX=nan([size(X.(fn{1})),length(M.patient.list)]);
    end
    allX(:,:,pt)=X.(fn{1});
end
save([get(handles.groupdir_choosebox,'String'),'NBSdataMatrix'],'allX','-v7.3');



% --- Executes on button press in lc_nbsadvanced.
function lc_nbsadvanced_Callback(hObject, eventdata, handles)
% hObject    handle to lc_nbsadvanced (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ea_nbs_advanced;
