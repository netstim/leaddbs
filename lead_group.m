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

% Last Modified by GUIDE v2.5 21-Dec-2016 16:19:08

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
% uiwait(handles.leadfigure);


% Build popup tables:

% atlassets:
options.earoot=ea_getearoot;
options.prefs=ea_prefs('');
setappdata(handles.leadfigure,'earoot',options.earoot);
as=dir([ea_space(options,'atlases')]);
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
set(handles.atlassetpopup,'String',asc);
[~,defix]=ismember(options.prefs.atlases.default,asc);
set(handles.atlassetpopup,'Value',defix);



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

setappdata(handles.leadfigure,'genvatfunctions',genvatfunctions);
setappdata(handles.leadfigure,'vatfunctionnames',ndc);




% get electrode model specs and place in popup
set(handles.elmodelselect,'String',[{'Patient specified'},ea_resolve_elspec]);

% set background image
set(gcf,'color','w');
im=imread([earoot,'icons',filesep,'logo_lead_group.png']);
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


ll=dir([ea_space(options,'labeling'),'*.nii']);
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
    M=ea_initializeM;
end
setappdata(gcf,'M',M);
ea_refresh_lg(handles);

handles.prod='group';
ea_firstrun(handles,options);

ea_menu_initmenu(handles,{'prefs','transfer'});

ea_processguiargs(handles,varargin)


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
M=getappdata(handles.leadfigure,'M');

folders=ea_uigetdir(ea_startpath,'Select Patient folders..');
M.patient.list=[M.patient.list;folders'];
M.patient.group=[M.patient.group;ones(length(folders),1)];

setappdata(handles.leadfigure,'M',M);
ea_refresh_lg(handles);
% save M
M=getappdata(handles.leadfigure,'M');
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
try M.stimparams(deleteentry)=[]; end

for cvar=1:length(M.clinical.vars)
    M.clinical.vars{cvar}(deleteentry,:)=[];
end

try
    M.stats(deleteentry)=[];
end


setappdata(gcf,'M',M);
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
options.d3.mirrorsides=get(handles.mirrorsides,'Value');
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

for reg=1:length(options.d3.isomatrix)
try options.d3.isomatrix{reg}=ea_reformat_isomatrix(options.d3.isomatrix{reg},M,options); end
end
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
options.groupmode=1;
options.patient_list=M.patient.list;
resultfig=ea_elvis(options,M.elstruct(get(handles.patientlist,'Value')));

ea_busyaction('off',handles.leadfigure,'group');


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
ea_refresh_lg(handles);






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
ea_refresh_lg(handles);



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
ea_refresh_lg(handles);

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
        S.label='gs';
        [ea_stats,usewhichstim]=ea_assignstimcnt(M.stats(pt).ea_stats,S);
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

ea_refresh_lg(handles);


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

ea_refresh_lg(handles);

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
%stimname=ea_detstimname();

% determine if fMRI or dMRI
mods=get(handles.fiberspopup,'String');
mod=mods{get(handles.fiberspopup,'Value')};
switch mod
    case {'Patient-specific fiber tracts','rest_tc'}
        fibersfile=mod;
    case 'Do not calculate connectivity stats'
    otherwise % load fibertracts once and for all subs here.
        [fibersfile.fibers,fibersfile.fibersidx]=ea_loadfibertracts([ea_getconnectomebase('dmri'),mod,filesep,'data.mat']);
end

[selection]=ea_groupselectorwholelist(M.ui.listselect,M.patient.list);

for pt=selection

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
    options.d3.elrendering=1; % hard code to viz electrodes in this setting.
    options.d3.colorpointcloud=0;
    options.native=0;

    options.d3.hlactivecontacts=get(handles.highlightactivecontcheck,'Value');
    options.d3.showactivecontacts=get(handles.showactivecontcheck,'Value');
    options.d3.showpassivecontacts=get(handles.showpassivecontcheck,'Value');
    try
        options.d3.isomatrix=M.isomatrix;
    catch
        options.d3.isomatrix={};
    end

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
        vfnames=getappdata(handles.leadfigure,'vatfunctionnames');

        [~,ix]=ismember(M.vatmodel,vfnames);
        vfs=getappdata(handles.leadfigure,'genvatfunctions');
        try
            ea_genvat=eval(['@',vfs{ix}]);
        catch
            keyboard
        end
        setappdata(handles.leadfigure,'resultfig',resultfig);
        for side=1:2
            setappdata(resultfig,'elstruct',M.elstruct(pt));
            setappdata(resultfig,'elspec',options.elspec);
 %           try
                [stimparams(1,side).VAT(1).VAT,volume]=feval(ea_genvat,M.elstruct(pt).coords_mm,M.S(pt),side,options,'gs',0.2,handles.leadfigure);
%            catch
%                ea_error(['Error while creating VTA of ',M.patient.list{pt},'.']);
%            end
            stimparams(1,side).volume=volume;
        end

        setappdata(resultfig,'stimparams',stimparams(1,:));
    end
    % this will add the volume stats (atlasIntersections) to stats file:
    ea_showfibres_volume(resultfig,options);


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
                    pV=spm_vol([ea_space(options,'labeling'),selectedparc,'.nii']);
                    pX=spm_read_vols(pV);
                    ea_cvshowvatfmri(resultfig,pX,directory,filesare,handles,pV,selectedparc,options);
                otherwise
                    ea_cvshowvatdmri(resultfig,directory,{mod,'gs'},selectedparc,options);
            end
        else
            ea_cvshowvatdmri(resultfig,directory,{fibersfile,'gs'},selectedparc,options);
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
%         ea_dispercent(1,'end');
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

ea_refresh_lg(handles);





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
ea_refresh_lg(handles);

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
ea_refresh_lg(handles);

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
ea_refresh_lg(handles);

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
ea_busyaction('on',handles.leadfigure,'group');

nudir=[nudir,filesep];
M=ea_initializeM;


set(handles.groupdir_choosebox,'String',nudir);

try % if file already exists, load it (and overwrite M).
    load([nudir,'LEAD_groupanalysis.mat']);
catch % if not, store it saving M.
    save([nudir,'LEAD_groupanalysis.mat'],'M','-v7.3');
end

M.ui.groupdir=nudir;
setappdata(handles.leadfigure,'M',M);
try
    setappdata(handles.leadfigure,'S',M.S);
    setappdata(handles.leadfigure,'vatmodel',M.S(1).model);
end

ea_busyaction('off',handles.leadfigure,'group');

ea_refresh_lg(handles);



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
        uisetcolor(M.groups.color(ismember(M.groups.group,g),:),['Group ',num2str(g),':']);
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
ea_refresh_lg(handles);

ea_stimparams(M.elstruct,handles.leadfigure,options);



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
ea_refresh_lg(handles);


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
ea_refresh_lg(handles);

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
ea_refresh_lg(handles);

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
ea_refresh_lg(handles);


% --- Executes on button press in colorpointcloudcheck.
function colorpointcloudcheck_Callback(hObject, eventdata, handles)
% hObject    handle to colorpointcloudcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of colorpointcloudcheck
M=getappdata(gcf,'M');
M.ui.colorpointcloudcheck=get(handles.colorpointcloudcheck,'Value');
setappdata(gcf,'M',M);
ea_refresh_lg(handles);

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
if ~strcmp(get(handles.groupdir_choosebox,'String'),'Choose Group Directory') % group dir still not chosen
    disp('Saving data...');
    % save M
    ea_refresh_lg(handles);
    M=getappdata(hObject,'M');
    disp('Saving data to disk...');
    try
        save([get(handles.groupdir_choosebox,'String'),'LEAD_groupanalysis.mat'],'M','-v7.3');
    catch
        warning('Data could not be saved.');
        keyboard
    end
    disp('Done.');
    disp('Bye for now.');
end
ea_busyaction('off',gcf,'group');
delete(hObject);


% --- Executes on button press in targetreport.
function targetreport_Callback(hObject, eventdata, handles)
% hObject    handle to targetreport (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ea_refresh_lg(handles);
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
try options.d3.isomatrix=M.isomatrix;
catch
    options.d3.isomatrix={};
end
try options.d3.isomatrix_name=M.isomatrix_name;
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
options.modality=3; % use template image
options=ea_amendtoolboxoptions(options);

if strcmp(options.atlasset,'Use none');
    options.d2.writeatlases=1;
else
    options.d2.writeatlases=1;
end

% Prior Results are loaded here inside the function (this way, function
% can be called just by giving the patient directory.

% Prepare isomatrix (includes a normalization step if M.ui.normregpopup
% says so:


if options.d3.showisovolume || options.expstatvat.do % regressors be used ? iterate through all
    allisomatrices=options.d3.isomatrix;
    allisonames=options.d3.isomatrix_name;
    for reg=1:length(allisomatrices)
        options.d3.isomatrix=allisomatrices{reg};
        options.d3.isomatrix_name=allisonames{reg};
        M.isomatrix=allisomatrices{reg};
        M.isomatrix_name=allisonames{reg};
        options.shifthalfup=0;
        try options.d3.isomatrix=ea_reformat_isomatrix(options.d3.isomatrix,M,options);
        if size(options.d3.isomatrix{1},2)==3 % pairs
        options.shifthalfup=1;
        end
        end

        if ~strcmp(get(handles.groupdir_choosebox,'String'),'Choose Group Directory') % group dir still not chosen
            ea_refresh_lg(handles);
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


% --- Executes on button press in lc_SPM.
function lc_SPM_Callback(hObject, eventdata, handles)
% hObject    handle to lc_SPM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
M=getappdata(handles.leadfigure,'M');

gecs=get(handles.lc_graphmetric,'String');
[~,gecs]=fileparts(gecs{M.ui.lc.graphmetric});
parc=get(handles.labelpopup,'String');
parc=parc{get(handles.labelpopup,'Value')};

spmdir=[M.ui.groupdir,'connectomics',filesep,parc,filesep,'graph',filesep,gecs,filesep,'SPM'];
try
rmdir(spmdir,'s');
end
mkdir([M.ui.groupdir,'connectomics']);
mkdir([M.ui.groupdir,'connectomics',filesep,parc]);
mkdir([M.ui.groupdir,'connectomics',filesep,parc,filesep,'graph']);

mkdir([M.ui.groupdir,'connectomics',filesep,parc,filesep,'graph',filesep,gecs]);
mkdir(spmdir);


for sub=1:length(M.patient.list)
    if M.ui.lc.normalization==1;
        zzz='';
    elseif M.ui.lc.normalization==2;
        zzz='z';
        ea_histnormalize([M.patient.list{sub},filesep,'connectomics',filesep,parc,filesep,'graph',filesep,gecs,'.nii,1'],2);
    elseif M.ui.lc.normalization==3;
        zzz='k';
        ea_histnormalize([M.patient.list{sub},filesep,'connectomics',filesep,parc,filesep,'graph',filesep,gecs,'.nii,1'],3);
    end
    if M.ui.lc.smooth
        sss='s';
        ea_smooth([M.patient.list{sub},filesep,'connectomics',filesep,parc,filesep,'graph',filesep,zzz,gecs,'.nii,1']);
    else
        sss='';
    end
    fis{sub}=[M.patient.list{sub},filesep,'connectomics',filesep,parc,filesep,'graph',filesep,sss,zzz,gecs,'.nii,1'];
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

msk=isnan(nii.img);
msk=msk+(nii.img==0);
vals=nii.img(~msk);
ext=ext(1:end-2);
switch zzz
    case 2 % zscore
        nii.fname=[pth,filesep,'z',fn,ext];
        vals=zscore(vals(:));
    case 3 % albada
        nii.fname=[pth,filesep,'k',fn,ext];
        vals=ea_normal(vals(:));
end
nii.img(~msk)=vals;
spm_write_vol(nii,nii.img);

function ea_smooth(fname)
[pth,fn,ext]=fileparts(fname);

spm_smooth(fname,[pth,fn,'s',ext],[8 8 8]);





% --- Executes on button press in specify2doptions.
function specify2doptions_Callback(hObject, eventdata, handles)
% hObject    handle to specify2doptions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
options.prefs=ea_prefs('');
options.groupmode=1;
options.native=0;
ea_spec2dwrite(options);


% --- Executes on button press in calcgroupconnectome.
function calcgroupconnectome_Callback(hObject, eventdata, handles)
% hObject    handle to calcgroupconnectome (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

options.prefs=ea_prefs('tmp');
M=getappdata(gcf,'M');

normalized_fibers_mm=[]; % combined connectome
allidx=[];
ea_dispercent(0,'Concatenating connectome');
maxfibno=0;
for sub=1:length(M.patient.list)
    ea_dispercent(sub/length(M.patient.list));

    [nfibs,idx]=ea_loadfibertracts([M.patient.list{sub},filesep,'connectomes',filesep,'dMRI',filesep,options.prefs.FTR_normalized]);
    idx=idx(1:20000); % only use first 20k fibers of each subject.
    sumidx=sum(idx);
    nfibs=nfibs(1:sumidx,:);
    nfibs(:,4)=nfibs(:,4)+maxfibno; % add offset
    normalized_fibers_mm=[normalized_fibers_mm;nfibs];
    allidx=[allidx;idx];
    maxfibno=max(normalized_fibers_mm(:,4));
end
ea_dispercent(1,'end');
ea_savefibertracts([M.ui.groupdir,'connectomes',filesep,'dMRI',filesep,options.prefs.FTR_normalized],normalized_fibers_mm,allidx,'mm');




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




% --- Executes on button press in lc_nbs.
function lc_nbs_Callback(hObject, eventdata, handles)
% hObject    handle to lc_nbs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

clearvars -global nbs
global nbs

earoot=getappdata(handles.leadfigure,'earoot');
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
        uT=T;
        clear tstat
        clear pmask
        pmask=zeros([nbs.NBS.n,size(T)]);
        for network=1:nbs.NBS.n;
            X=full(nbs.NBS.con_mat{network});
            X=X+X';
            pmask(network,:,:)=X;
            save([root,'sig_',num2str(network)],'X');
        end

        if network
        pmask=squeeze(sum(pmask,1));
        end
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

if ismember('&',thismetr) % clean from & for variable
    [~,ix]=ismember('&',thismetr);
    thismetr(ix)='';
end
eval([thismetr,'=T;']);
expfn=[root,thisparc,'_',thismetr,'T_p<',UI.alpha.ui,'.mat'];
save(expfn,thismetr);
eval([thismetr,'=uT;']);
expfn=[root,thisparc,'_',thismetr,'T_unthresholded.mat'];
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

M=getappdata(handles.leadfigure,'M');
thisparc=get(handles.labelpopup,'String');
thisparc=thisparc{get(handles.labelpopup,'Value')};

thismetr=get(handles.lc_metric,'String');
thismetr=thismetr{get(handles.lc_metric,'Value')};

if ismember('&',thismetr)
    % e.g. ON vs OFF metric - need to load two matrices per subject!
    [~,plusix]=ismember('&',thismetr);
    [scix]=strfind(thismetr,'_');
    suffx1=thismetr(scix(1)+1:plusix-1);
    suffx2=thismetr(plusix+1:scix(2)-1);
    base=thismetr(1:scix(1)-1);
    cap=thismetr(scix(2)+1:end);
    clear thismetr
    metrix{1}=[base,'_',suffx1,'_',cap];
    metrix{2}=[base,'_',suffx2,'_',cap];

    % inflate design matrix
    if ~all(mX==1)
        ea_error('Multi subject measurements only supported for a single cohort. All group values must be 1.');
    end

    mX=repmat(eye(2),length(mX),1);
    save([get(handles.groupdir_choosebox,'String'),'NBSdesignMatrix'],'mX');
else
 metrix{1}=thismetr;
end

Xcnt=1;
for pt=1:length(M.patient.list)
    for metr=1:length(metrix)
        X=load([M.patient.list{pt},filesep,'connectomics',filesep,thisparc,filesep,metrix{metr},'.mat']);
        fn=fieldnames(X);
        if ~exist('allX','var')
            allX=nan([size(X.(fn{1})),length(M.patient.list)]);
        end

        X=X.(fn{1});
        switch get(handles.normregpopup,'Value');
            case 2
                X(:)=ea_nanzscore(X(:));
            case 3
                X(:)=ea_normal(X(:));
        end

        allX(:,:,Xcnt)=X;
        Xcnt=Xcnt+1;
    end
end
save([get(handles.groupdir_choosebox,'String'),'NBSdataMatrix'],'allX','-v7.3');



% --- Executes on button press in lc_nbsadvanced.
function lc_nbsadvanced_Callback(hObject, eventdata, handles)
% hObject    handle to lc_nbsadvanced (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ea_nbs_advanced;


% --- Executes on button press in mirrorsides.
function mirrorsides_Callback(hObject, eventdata, handles)
% hObject    handle to mirrorsides (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of mirrorsides

M=getappdata(gcf,'M');
M.ui.mirrorsides=get(handles.mirrorsides,'Value');


setappdata(gcf,'M',M);
ea_refresh_lg(handles);
