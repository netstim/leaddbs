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

% Last Modified by GUIDE v2.5 15-Jul-2015 14:23:43

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
options.earoot=[fileparts(which('lead')),filesep];
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



% setup modelselect popup

cnt=1;
earoot=[fileparts(which('lead')),filesep];
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

set(handles.modelselect,'String',ndc);



% get electrode model specs and place in popup
set(handles.elmodelselect,'String',[{'Patient specified'},ea_resolve_elspec]);

% set background image
set(gcf,'color','w');
im=imread('ea_logo.png');
image(im);
axis off;
axis equal;


% Fibers:
% Fibers:


fibd=dir([options.earoot,'fibers',filesep,'*.mat']);
fiberscell{1}='Patient-specific DTI-Data';

for fd=2:length(fibd)+1
    [~,fn]=fileparts(fibd(fd-1).name);
    fiberscell{fd}=fn;
end

set(handles.fiberspopup,'String',fiberscell);

try
    priorselection=find(ismember(fiberscell,stimparams.usefiberset)); % retrieve prior selection of fiberset.
    set(handles.fiberspopup,'Value',priorselection);
    
catch    % reinitialize using third entry.
    set(handles.fiberspopup,'Value',4);
    
    
end

% Labels:


ll=dir([options.earoot,'templates',filesep,'labeling',filesep,'*.nii']);
for lab=1:length(ll)
    [~,n]=fileparts(ll(lab).name);
    labelcell{lab}=n;
end

set(handles.labelpopup,'String',labelcell);
set(handles.lc_parcellation,'String',labelcell);

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
vnum=ea_getvsn('local');
set(handles.versiontxt,'String',['v ',num2str(vnum(1))]);


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
save([get(handles.groupdir_choosebox,'String'),'LEAD_groupanalysis.mat'],'M');





% --- Executes on button press in removeptbutton.
function removeptbutton_Callback(hObject, eventdata, handles)
% hObject    handle to removeptbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
M=getappdata(gcf,'M');

deleteentry=get(handles.patientlist,'Value');

M.patient.list(deleteentry)=[];

M.patient.group(deleteentry)=[];

try M.elstruct(deleteentry)=[];
    M.stimparams(deleteentry)=[]; end

for cvar=1:length(M.clinical.vars)
    M.clinical.vars{cvar}(deleteentry)=[];
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
    save([get(handles.groupdir_choosebox,'String'),'LEAD_groupanalysis.mat'],'M');
    disp('Done.');
end


resultfig=ea_render_view(options,M.elstruct(get(handles.patientlist,'Value')));



% --- Executes on button press in corrbutton.
function corrbutton_Callback(hObject, eventdata, handles)
% hObject    handle to corrbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

stats=preparedataanalysis(handles);


assignin('base','stats',stats);


% perform correlations:


if ~isempty(stats.vicorr.both)
    ea_corrplot([stats.corrcl,stats.vicorr.both],'Volume Intersections, both hemispheres',stats.vc_labels);
    ea_corrplot([stats.corrcl,stats.vicorr.nboth],'Volume Intersections, normalized, both hemispheres',stats.vc_labels);
end
if ~isempty(stats.vicorr.right)
    ea_corrplot([stats.corrcl,stats.vicorr.right],'Volume Intersections, right hemisphere',stats.vc_labels);
    ea_corrplot([stats.corrcl,stats.vicorr.nright],'Volume Intersections, normalized, right hemisphere',stats.vc_labels);
end
if ~isempty(stats.vicorr.left)
    ea_corrplot([stats.corrcl,stats.vicorr.left],'Volume Intersections, left hemisphere',stats.vc_labels);
    ea_corrplot([stats.corrcl,stats.vicorr.nleft],'Volume Intersections, normalized, left hemisphere',stats.vc_labels);
end


if ~isempty(stats.fccorr.both)
    ea_corrplot([stats.corrcl,stats.fccorr.both],'Fibercounts, both hemispheres',stats.fc_labels);
    
    ea_corrplot([stats.corrcl,stats.fccorr.nboth],'Fibercounts, normalized, both hemispheres',stats.fc_labels);
end


if ~isempty(stats.fccorr.left)
    ea_corrplot([stats.corrcl,stats.fccorr.right],'Fibercounts, right hemisphere',stats.fc_labels);
    
    ea_corrplot([stats.corrcl,stats.fccorr.nright],'Fibercounts, normalized, right hemisphere',stats.fc_labels);
end



if ~isempty(stats.fccorr.left)
    ea_corrplot([stats.corrcl,stats.fccorr.left],'Fibercounts, left hemisphere',stats.fc_labels);
    
    ea_corrplot([stats.corrcl,stats.fccorr.nleft],'Fibercounts, normalized, left hemisphere',stats.fc_labels);
end




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
set(gcf,'name','LEAD-DBS Groupanalysis (updating...)');
drawnow

% get model data

M=getappdata(gcf,'M');

if strcmp(get(handles.groupdir_choosebox,'String'),'Choose Group Directory') % not set yet.
    set(gcf,'name','LEAD-DBS Groupanalysis');
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


thisparc=get(handles.lc_parcellation,'String');
thisparc=thisparc{get(handles.lc_parcellation,'Value')};
gmdir=dir([M.patient.list{1},'connectomics',filesep,thisparc,filesep,'graph',filesep,'*.nii']);

gms{1}='';
for gm=1:length(gmdir)
    gms{gm}=gmdir(gm).name;
end
set(handles.lc_graphmetric,'String',gms);

% update UI

% 2D-Viz options:

if isfield(M.ui,'bbsize');    set(handles.bbsize,'String',M.ui.bbsize); end
if isfield(M.ui,'tdcolorscheck');    set(handles.tdcolorscheck,'Value',M.ui.tdcolorscheck); end
if isfield(M.ui,'tdcontourcheck');    set(handles.tdcontourcheck,'Value',M.ui.tdcontourcheck); end
if isfield(M.ui,'tdlabelcheck');    set(handles.tdlabelcheck,'Value',M.ui.tdlabelcheck); end
if isfield(M.ui,'tdlegendcheck');    set(handles.tdlegendcheck,'Value',M.ui.tdlegendcheck); end
if isfield(M.ui,'tdcontourcolor');    setappdata(handles.tdcontourcolor,'color',M.ui.tdcontourcolor); end



% make setstimparams button green if set.
if isfield(M,'stimparams')
    if isfield(M.stimparams,'U')
        set(handles.setstimparamsbutton,'BackgroundColor',[0.1;0.8;0.1]);
    else
        set(handles.setstimparamsbutton,'BackgroundColor',[0.93,0.93,0.93]);
    end
else
    set(handles.setstimparamsbutton,'BackgroundColor',[0.93,0.93,0.93]);
end




% make chooscolors button green if chosen.
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
try set(handles.modelselect,'Value',M.ui.modelselect); end
try set(handles.elmodelselect,'Value',M.ui.elmodelselect); end
try set(handles.normregpopup,'Value',M.ui.normregpopup); end
try set(handles.lc_parcellation,'Value',M.ui.lc.parcellation); end
try set(handles.lc_normalization,'Value',M.ui.lc.normalization); end
try set(handles.lc_graphmetric,'Value',M.ui.lc.graphmetric); end

% update enable-disable-dependencies:
if M.ui.elrendering==3
    try set(handles.colorpointcloudcheck,'Enable','on'); end
else
    try set(handles.colorpointcloudcheck,'Enable','off'); end
end
% hide detachbutton if already detached:
try
    if M.ui.detached
        set(handles.detachbutton,'Visible','off');
    end
end


if ~isempty(M.patient.list)
    for pt=1:length(M.patient.list)
        % set stimparams based on values provided by user
        for side=1:2
            
            M.stimparams(pt,side).usefiberset=get(handles.fiberspopup,'String');
            
            M.stimparams(pt,side).usefiberset=M.stimparams(pt,side).usefiberset{M.ui.fiberspopup};
            
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
        
        try
            load([M.patient.list{pt},filesep,'ea_reconstruction']);
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
        catch
            %warning('No reconstruction present in folder. Using information stored in group-file.');
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
            setappdata(gcf,'M',M);
            set(gcf,'name','LEAD-DBS Groupanalysis');
            return
        end
        
        priorvilist=M.vilist;
        try
            M.vilist=M.stats(pt).ea_stats.atlases.names;
        end
        % check and compare with prior atlas intersection list.
        
        if ~isempty(priorvilist) && ~isequal(priorvilist,M.vilist)
            
            warning('Patient stats are inhomogeneous. Please re-run group analysis (Section Prepare DBS stats).');
        end
        
        
        
        
        try
            priorfclist=M.fclist;
            M.fclist=ea_stats.stimulation(1).ft(1).labels{1};
            fcdone=1;
        catch
            fcdone=0;
            
        end
        
        % check and compare with prior fibertracking list.
        
        if fcdone
            
            if ~isempty(priorfclist) && ~isequal(priorfclist,M.fclist)
                warning('Trying to analyse inhomogeneous patient group. Please re-run single subject lead analysis with patients using always the same labeling atlas.');
            end
            
        end
    end
    try
        setappdata(gcf,'elstruct',elstruct);
    end
else
    M.vilist={};
    M.fclist={};
end


% store everything in Model
setappdata(gcf,'M',M);

% refresh UI

set(handles.vilist,'String',M.vilist);
set(handles.fclist,'String',M.fclist);




set(gcf,'name','LEAD-DBS Groupanalysis');


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
        for side=1:size(M.stats(pt).ea_stats.stimulation(usewhichstim).vat,1)
            for vat=1:size(M.stats(pt).ea_stats.stimulation(usewhichstim).vat,2);
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
M=getappdata(gcf,'M');

% set options
options=ea_setopts_local(handles);

for pt=1:length(M.patient.list)
    
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
    options.numcontacts=size(M.elstruct(pt).coords_mm{1},1);
    options.elmodel=M.elstruct(pt).elmodel;
    options=ea_resolve_elspec(options);
    options.prefs=ea_prefs(options.patientname);
    options.d3.verbose='off';
    
    
    
    
    options.d3.elrendering=M.ui.elrendering;
    options.d3.hlactivecontacts=get(handles.highlightactivecontcheck,'Value');
    options.d3.showactivecontacts=get(handles.showactivecontcheck,'Value');
    options.d3.showpassivecontacts=get(handles.showpassivecontcheck,'Value');
    try options.d3.isomatrix=M.isomatrix; end
    
    options.d3.isovscloud=M.ui.isovscloudpopup;
    options.d3.showisovolume=M.ui.showisovolumecheck;
    
    options.expstatvat.do=M.ui.statvat;
    try
        options.expstatvat.vars=M.clinical.vars(M.ui.clinicallist);
        options.expstatvat.labels=M.clinical.labels(M.ui.clinicallist);
        options.expstatvat.pt=pt;
    end
    options.expstatvat.dir=M.ui.groupdir;
    processlocal=0;
    try
        if M.ui.detached
            processlocal=1;
            mkdir([M.ui.groupdir,'tmp']);
            options.root=M.ui.groupdir;
            options.patientname='tmp';
            
            ea_stats=M.stats(pt).ea_stats;
            coords_mm=M.elstruct(pt).coords_mm;
            trajectory=M.elstruct(pt).trajectory;
            save([M.ui.groupdir,'tmp',filesep,'ea_stats'],'ea_stats');
            save([M.ui.groupdir,'tmp',filesep,'ea_reconstruction'],'coords_mm','trajectory');
        end
    catch
        if ~exist(options.root,'file') % data is not there. Act as if detached. Process in tmp-dir.
            processlocal=1;
            warning('on');
            warning('Data has been detached from group-directory. Will process locally. Please be aware that you might loose this newly-processed data once you re-attach the single-patient data to the analysis!');
            warning('off');
            mkdir([M.ui.groupdir,'tmp']);
            options.root=M.ui.groupdir;
            options.patientname='tmp';
            
            ea_stats=M.stats(pt).ea_stats;
            coords_mm=M.elstruct.coords_mm;
            trajectory=M.elstruct.trajectory;
            save([M.ui.groupdir,'tmp',filesep,'ea_stats'],'ea_stats');
            save([M.ui.groupdir,'tmp',filesep,'ea_reconstruction'],'coords_mm','trajectory');
        end
    end
    
        delete([options.root,options.patientname,filesep,'ea_stats.mat']);
    
    % Step 1: Re-calculate closeness to subcortical atlases.
    resultfig=ea_render_view(options);
    % save scene as matlab figure
    
    
    
    
    % Step 2: Re-calculate Fibertracts / VAT
    setappdata(resultfig,'stimparams',M.stimparams(pt,:));
    ea_showfibres_volume(resultfig,options);
    
    close(resultfig);
    
    if processlocal % gather stats and recos to M
        load([M.ui.groupdir,'tmp',filesep,'ea_stats']);
        load([M.ui.groupdir,'tmp',filesep,'ea_reconstruction']);
        
        M.stats(pt).ea_stats=ea_stats;
        M.elstruct(pt).coords_mm=coords_mm;
        M.elstruct(pt).trajectory=trajectory;
        setappdata(gcf,'M',M);
        
        save([M.ui.groupdir,'LEAD_groupanalysis.mat'],'M');
        try      movefile([options.root,options.patientname,filesep,'LEAD_scene.fig'],[M.ui.groupdir,'LEAD_scene_',num2str(pt),'.fig']); end
        rmdir([M.ui.groupdir,'tmp'],'s');
        
        
    end
end
%% processing done here.

% calculate group-results for expstatvat if required:
if options.expstatvat.do
    disp('Averaging VAT-Stat files to a group result.');
    for side=1:2
        switch side
            case 1
                si='rh';
            case 2
                si='lh';
        end
        fis=dir([M.ui.groupdir,'statvat_results',filesep,'*_',si,'.nii']);
        ea_dispercent(0,['Reading in files, ',si]);
        for fi=1:length(fis);
            ea_dispercent(fi/length(fis));
            vV=spm_vol([M.ui.groupdir,'statvat_results',filesep,fis(fi).name]);
            if ~exist('thisvat','var')
                tmp=spm_read_vols(vV);
                thisvat=zeros([size(tmp),length(fis)]);
                thisvat(:,:,:,1)=tmp;
                clear tmp
            else
                thisvat(:,:,:,fi)=spm_read_vols(vV);
            end
        end
        ea_dispercent(100,'end');
        disp('Averaging...');
        thisvat(thisvat==0)=nan;
        
        thisvat=nanmean(thisvat,4);
        disp('Done.');
        disp('Saving file...');
        %thisvat=thisvat/fi; % simple mean here.
        vV.fname=[M.ui.groupdir,'statvat_results',filesep,si,'_mean.nii'];
        vV.dt=[64,0];
        spm_write_vol(vV,thisvat);
        disp('Done.');
    end
end

refreshvifc(handles);


% --- Executes on selection change in modelselect.
function modelselect_Callback(hObject, eventdata, handles)
% hObject    handle to modelselect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns modelselect contents as cell array
%        contents{get(hObject,'Value')} returns selected item from modelselect
M=getappdata(gcf,'M');
M.ui.modelselect=get(handles.modelselect,'Value');
setappdata(gcf,'M',M);
refreshvifc(handles);


% --- Executes during object creation, after setting all properties.
function modelselect_CreateFcn(hObject, eventdata, handles)
% hObject    handle to modelselect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


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

options.earoot=[fileparts(which('lead')),filesep];
options.verbose=3;
options.sides=1:2; % re-check this later..
options.atlasset=get(handles.atlassetpopup,'String');
options.atlasset=options.atlasset{get(handles.atlassetpopup,'Value')};
options.fiberthresh=1;
options.writeoutstats=1;
options.labelatlas=get(handles.labelpopup,'String');
options.labelatlas=options.labelatlas{get(handles.labelpopup,'Value')};
options.writeoutpm=1;
options.colormap=jet;
options.d3.write=1;
options.d3.prolong_electrode=2;
options.d3.writeatlases=1;


% --- Executes on button press in groupdir_choosebox.
function groupdir_choosebox_Callback(hObject, eventdata, handles)
% hObject    handle to groupdir_choosebox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

nudir=[uigetdir];

if ~nudir % user pressed cancel
    return
end
nudir=[nudir,filesep];
M=initializeM;


set(handles.groupdir_choosebox,'String',nudir);

try % if file already exists, load it (and overwrite M).
    load([nudir,'LEAD_groupanalysis.mat']);
catch % if not, store it saving M.
    save([nudir,'LEAD_groupanalysis.mat'],'M');
end

M.ui.groupdir=nudir;
setappdata(gcf,'M',M);

refreshvifc(handles);



% --- Executes on button press in opensubgui.
function opensubgui_Callback(hObject, eventdata, handles)
% hObject    handle to opensubgui (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
M=getappdata(gcf,'M');
lead('loadsubs',M.patient.list);


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
M.ui.modelselect=1;
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

[pU,pIm]=getstimparams(M); % local function, gets out U and Im in cell format and returns zeros/thousands if not set.


for pt=1:length(M.patient.list)
    [~,ptname]=fileparts(M.patient.list{pt});
    Urproperties(pt) =    PropertyGridField(['double',num2str(pt)], pU{1}(pt,:), ...
        'Category', 'Voltage Right', ...
        'DisplayName', ptname, ...
        'Description', 'Double Matrix.');
    Ulproperties(pt) =    PropertyGridField(['double',num2str(pt)], pU{2}(pt,:), ...
        'Category', 'Voltage Left', ...
        'DisplayName', ptname, ...
        'Description', 'Double Matrix.');
    Irproperties(pt) =    PropertyGridField(['double',num2str(pt)], pIm{1}(pt,:), ...
        'Category', 'Impedance Right', ...
        'DisplayName', ptname, ...
        'Description', 'Double Matrix.');
    Ilproperties(pt) =    PropertyGridField(['double',num2str(pt)], pIm{2}(pt,:), ...
        'Category', 'Impedance Left', ...
        'DisplayName', ptname, ...
        'Description', 'Double Matrix.');
end




% arrange flat list into a hierarchy based on qualified names
Urproperties = Urproperties.GetHierarchy();
Irproperties = Irproperties.GetHierarchy();
Ulproperties = Ulproperties.GetHierarchy();
Ilproperties = Ilproperties.GetHierarchy();

% create figure
f = figure( ...
    'MenuBar', 'none', ...
    'Name', 'Please enter stimulation parameters for the group.', ...
    'NumberTitle', 'off', ...
    'Toolbar', 'none');

% procedural usage
gU{1} = PropertyGrid(f, ...            % add property pane to figure
    'Properties', Urproperties,'Position', [0 0 0.25 1]);     % set properties explicitly;
gU{2} = PropertyGrid(f, ...            % add property pane to figure
    'Properties', Ulproperties,'Position', [0.25 0 0.25 1]);     % set properties explicitly;
gI{1} = PropertyGrid(f, ...            % add property pane to figure
    'Properties', Irproperties,'Position', [0.5 0 0.25 1]);     % set properties explicitly;
gI{2} = PropertyGrid(f, ...            % add property pane to figure
    'Properties', Ilproperties,'Position', [0.75 0 0.25 1]);     % set properties explicitly;


% declarative usage, bind object to grid
obj = SampleObject;  % a value object

% update the type of a property assigned with type autodiscovery
%userproperties = PropertyGridField.GenerateFrom(obj);
%userproperties.FindByName('IntegerMatrix').Type = PropertyType('denserealdouble', 'matrix');
%disp(userproperties.FindByName('IntegerMatrix').Type);

% wait for figure to close
uiwait(f);

%disp(gUr.GetPropertyValues());





options=ea_setopts_local(handles);

for pt=1:length(M.patient.list)
    
    % set pt specific options
    options.root=[fileparts(M.patient.list{pt}),filesep];
    [~,options.patientname]=fileparts(M.patient.list{pt});
    if ~isempty(M.elstruct(pt).elmodel)
        options.elmodel=M.elstruct(pt).elmodel;
    else
        options.elmodel='Medtronic 3389';
    end
    options=ea_resolve_elspec(options);
    options.prefs=ea_prefs(options.patientname);
    options.d3.verbose='off';
    
    % assign correct .m-file to function.
    genvatfunctions=getappdata(gcf,'genvatfunctions');
    ea_genvat=eval(['@',genvatfunctions{get(handles.modelselect,'Value')}]);
    
    % set stimparams based on values provided by user
    
    
    for side=1:2
        M.stimparams(pt,side).U=gU{side}.Properties(pt).Value;
        M.stimparams(pt,side).Im=gI{side}.Properties(pt).Value;
        
        M.stimparams(pt,side).usefiberset=get(handles.fiberspopup,'String');
        
        M.stimparams(pt,side).usefiberset=M.stimparams(pt,side).usefiberset{M.ui.fiberspopup};
        
        M.stimparams(pt,side).labelatlas={options.labelatlas};
        M.stimparams(pt,side).showfibers=1;
        M.stimparams(pt,side).fiberthresh=1;
        [M.stimparams(pt,side).VAT.VAT,radius,volume]=feval(ea_genvat,M.elstruct(pt).coords_mm,M.stimparams(pt,:),side,options);
        
        M.stimparams(pt,side).radius=radius;
        M.stimparams(pt,side).volume=volume;
        M.stimparams(pt,side).showconnectivities=1;
        M.elstruct(pt).activecontacts{side}=find(M.stimparams(pt,side).U);
    end
    
    
    
end
try M.stimparams(pt+1:end,:)=[]; end % remove any additional entries if present.
setappdata(gcf,'M',M);
refreshvifc(handles);




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





% --- Executes on button press in stimoldbutn.
function stimoldbutn_Callback(hObject, eventdata, handles)
% hObject    handle to stimoldbutn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
M=getappdata(gcf,'M');

try
    uicell=inputdlg('Enter Variable name for Voltage-Parameters','Enter Stimulation Settings...',1);
    uidata.U=evalin('base',uicell{1});
catch
    warning('Stim-Params could not be evaluated. Please Try again.');
end
try
    uicell=inputdlg('Enter Variable name for Impedance-Parameters','Enter Stimulation Settings...',1);
    uidata.Im=evalin('base',uicell{1});
catch
    warning('Stim-Params could not be evaluated. Please Try again.');
end

options=ea_setopts_local(handles);
for pt=1:length(M.patient.list)
    
    % set pt specific options
    options.root=[fileparts(M.patient.list{pt}),filesep];
    [~,options.patientname]=fileparts(M.patient.list{pt});
    options.elmodel=M.elstruct(pt).elmodel;
    options=ea_resolve_elspec(options);
    options.prefs=ea_prefs(options.patientname);
    options.d3.verbose='off';
    % set stimparams based on values provided by user
    
    % assign correct .m-file to function.
    genvatfunctions=getappdata(gcf,'genvatfunctions');
    ea_genvat=eval(['@',genvatfunctions{get(handles.modelselect,'Value')}]);
    
    
    for side=1:2
        M.stimparams(pt,side).U=uidata.U{side}(pt,1:options.elspec.numel);
        M.stimparams(pt,side).Im=uidata.Im{side}(pt,1:options.elspec.numel);
        
        M.stimparams(pt,side).usefiberset=get(handles.fiberspopup,'String');
        
        M.stimparams(pt,side).usefiberset=M.stimparams(pt,side).usefiberset{get(handles.fiberspopup,'Value')};
        M.stimparams(pt,side).labelatlas={options.labelatlas};
        M.stimparams(pt,side).showfibers=1;
        M.stimparams(pt,side).fiberthresh=1;
        
        M.stimparams(pt,side).VAT.VAT=feval(ea_genvat,M.elstruct(pt).coords_mm,M.stimparams(pt,side),side,options);
        M.stimparams(pt,side).showconnectivities=1;
        M.elstruct(pt).activecontacts{side}=find(M.stimparams(pt,side).U);
    end
    
    
    
end
setappdata(gcf,'M',M);
refreshvifc(handles);


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


if ~strcmp(get(handles.groupdir_choosebox,'String'),'Choose Group Directory') % group dir still not chosen
    disp('Saving data...');
    % save M
    M=getappdata(hObject,'M');
    try
        save([get(handles.groupdir_choosebox,'String'),'LEAD_groupanalysis.mat'],'M');
    catch
        warning('Data could not be saved.');
    end
    disp('Bye for now.');
end

delete(hObject);


% --- Executes on button press in targetreport.
function targetreport_Callback(hObject, eventdata, handles)
% hObject    handle to targetreport (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
M=getappdata(gcf,'M');
ea_gentargetreport(M);


% --- Executes on button press in viz2dbutton.
function viz2dbutton_Callback(hObject, eventdata, handles)
% hObject    handle to viz2dbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc;
M=getappdata(gcf,'M');
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


options.d2.showlegend=M.ui.tdlegendcheck;


options.d3.isovscloud=M.ui.isovscloudpopup;
options.d3.showisovolume=M.ui.showisovolumecheck;
options.d3.colorpointcloud=M.ui.colorpointcloudcheck;
options.normregressor=M.ui.normregpopup;

options.d2.write=1;

options.d2.atlasopacity=0.15;
options.d2.col_overlay=get(handles.tdcolorscheck,'Value');
options.d2.con_overlay=get(handles.tdcontourcheck,'Value');
options.d2.con_color=getappdata(handles.tdcontourcolor,'color');
if isempty(options.d2.con_color)
    options.d2.con_color=[1,1,1]; % white
end

options.d2.lab_overlay=get(handles.tdlabelcheck,'Value');


options.d2.bbsize=str2double(get(handles.bbsize,'String'));

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
    disp('Saving data...');
    % save M
    save([get(handles.groupdir_choosebox,'String'),'LEAD_groupanalysis.mat'],'M');
    disp('Done.');
end


% export isovolume
if options.d3.showisovolume % export to nifti volume
    ea_exportisovolume(M.elstruct(get(handles.patientlist,'Value')),options);
end



cuts=ea_writeplanes(options,M.elstruct(get(handles.patientlist,'Value')));



function bbsize_Callback(hObject, eventdata, handles)
% hObject    handle to bbsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bbsize as text
%        str2double(get(hObject,'String')) returns contents of bbsize as a double
M=getappdata(gcf,'M');
M.ui.bbsize=get(handles.bbsize,'String');
setappdata(gcf,'M',M);
refreshvifc(handles);

% --- Executes during object creation, after setting all properties.
function bbsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bbsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in tdcolorscheck.
function tdcolorscheck_Callback(hObject, eventdata, handles)
% hObject    handle to tdcolorscheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tdcolorscheck
M=getappdata(gcf,'M');
M.ui.tdcolorscheck=get(handles.tdcolorscheck,'Value');
setappdata(gcf,'M',M);
refreshvifc(handles);

% --- Executes on button press in tdcontourcheck.
function tdcontourcheck_Callback(hObject, eventdata, handles)
% hObject    handle to tdcontourcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tdcontourcheck
M=getappdata(gcf,'M');
M.ui.tdcontourcheck=get(handles.tdcontourcheck,'Value');
setappdata(gcf,'M',M);
refreshvifc(handles);

% --- Executes on button press in tdlabelcheck.
function tdlabelcheck_Callback(hObject, eventdata, handles)
% hObject    handle to tdlabelcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tdlabelcheck
M=getappdata(gcf,'M');
M.ui.tdlabelcheck=get(handles.tdlabelcheck,'Value');
setappdata(gcf,'M',M);
refreshvifc(handles);

% --- Executes on button press in tdcontourcolor.
function tdcontourcolor_Callback(hObject, eventdata, handles)
% hObject    handle to tdcontourcolor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tdcol=uisetcolor;
setappdata(hObject,'color',tdcol);
M=getappdata(gcf,'M');
M.ui.tdcontourcolor=tdcol;
setappdata(gcf,'M',M);
refreshvifc(handles);


% --- Executes on button press in lc_SPM.
function lc_SPM_Callback(hObject, eventdata, handles)
% hObject    handle to lc_SPM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

spmdir=[M.ui.groupdir,'connectomics',filesep,get(handles.lc_parcellation,'String'),filesep,'graph',filesep,gcs,filesep,'SPM'];
rmdir(spmdir,'s');
mkdir([M.ui.groupdir,'connectomics']);
mkdir([M.ui.groupdir,'connectomics',filesep,get(handles.lc_parcellation,'String')]);
mkdir([M.ui.groupdir,'connectomics',filesep,get(handles.lc_parcellation,'String'),filesep,'graph']);
gcs=get(handles.lc_graphmetrics,'String');
[~,gcs]=gcs{M.ui.lc.graphmetrics};
mkdir([M.ui.groupdir,'connectomics',filesep,get(handles.lc_parcellation,'String'),filesep,'graph',filesep,gcs]);
mkdir(spmdir);


for sub=1:length(M.patient.list)
    if M.ui.lc.normalize==1;
        zzz='';
    elseif M.ui.lc.normalize==2;
        zzz='z';
        ea_histnormalize([M.patient.list{sub},'connectomics',filesep,get(handles.lc_parcellation,'String'),filesep,'graph',filesep,gcs,'.nii,1'],2);
    elseif M.ui.lc.normalize==3;
        zzz='k';
        ea_histnormalize([M.patient.list{sub},'connectomics',filesep,get(handles.lc_parcellation,'String'),filesep,'graph',filesep,gcs,'.nii,1'],3);
    end
    if M.ui.lc.smooth
        sss='s';
        ea_smooth([M.patient.list{sub},'connectomics',filesep,get(handles.lc_parcellation,'String'),filesep,'graph',filesep,zzz,gcs,'.nii,1']);
    else
        sss='';
    end
    fis{sub}=[M.patient.list{sub},'connectomics',filesep,get(handles.lc_parcellation,'String'),filesep,'graph',filesep,sss,zzz,gcs,'.nii,1'];
end

%% model specification:
matlabbatch{1}.spm.stats.factorial_design.dir = {spmdir};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = fis;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
cfg_util('run',{matlabbatch});
clear matlabbatch

%% model estimation:
matlabbatch{1}.spm.stats.fmri_est.spmmat = {spmdir};
matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
cfg_util('run',{matlabbatch});
clear matlabbatch

%% contrast manager:
matlabbatch{1}.spm.stats.con.spmmat = {spmdir};
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'main effect';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;
cfg_util('run',{matlabbatch});
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


% --- Executes on selection change in lc_parcellation.
function lc_parcellation_Callback(hObject, eventdata, handles)
% hObject    handle to lc_parcellation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns lc_parcellation contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lc_parcellation
M=getappdata(gcf,'M');

M.ui.lc.parcellation=get(handles.lc_parcellation,'Value');
setappdata(gcf,'M',M);

% --- Executes during object creation, after setting all properties.
function lc_parcellation_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lc_parcellation (see GCBO)
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


% --- Executes on button press in tdlegendcheck.
function tdlegendcheck_Callback(hObject, eventdata, handles)
% hObject    handle to tdlegendcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tdlegendcheck
M=getappdata(gcf,'M');
M.ui.tdlegendcheck=get(handles.tdlegendcheck,'Value');
setappdata(gcf,'M',M);
refreshvifc(handles);