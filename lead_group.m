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

% Last Modified by GUIDE v2.5 24-Jun-2014 15:46:36

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
% uiwait(handles.figure1);


% Build popup tables:

% atlassets:
options.earoot=[fileparts(which('lead')),filesep];
as=dir([options.earoot,'atlases',filesep]);
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


% Fibers:

fiberscell{1}='Gibbsconnectome';
fiberscell{2}='Gibbsconnectome5';
fiberscell{3}='Gibbsconnectome10';
fiberscell{4}='Gibbsconnectome20';
fiberscell{5}='Gibbsconnectome50';
fiberscell{6}='Gibbsconnectome100';
fiberscell{7}='Gibbsconnectome500';
fiberscell{8}='Patient-specific DTI-Data';

set(handles.fiberspopup,'String',fiberscell);

try
    priorselection=find(ismember(fiberscell,stimparams.usefiberset)); % retrieve prior selection of fiberset.
    set(handles.fiberspopup,'Value',priorselection);

catch    % reinitialize using third entry.
    set(handles.fiberspopup,'Value',3);
    
    
end

% Labels:


ll=dir([options.earoot,'templates',filesep,'labeling',filesep,'*.nii']);
for lab=1:length(ll)
    [~,n]=fileparts(ll(lab).name);
    labelcell{lab}=n;
end
labelcell{lab+1}='Use all';

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
set(handles.versiontxt,'String',ea_getvsn);


% make listboxes multiselectable:

set(handles.patientlist,'Max',100,'Min',0);
set(handles.grouplist,'Max',100,'Min',0);
set(handles.analysislist,'Max',100,'Min',0);
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
set(handles.analysislist,'Value',M.ui.listselect);

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
M.patient.analysis=[M.patient.analysis;zeros(length(folders),1)];

setappdata(gcf,'M',M);
refreshvifc(handles);





% --- Executes on button press in removeptbutton.
function removeptbutton_Callback(hObject, eventdata, handles)
% hObject    handle to removeptbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
M=getappdata(gcf,'M');

deleteentry=get(handles.patientlist,'Value');

M.patient.list(deleteentry)=[];

M.patient.group(deleteentry)=[];
M.patient.analysis(deleteentry)=[];

M.elstruct(deleteentry)=[];


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

options.d3.isovscloud=M.ui.isovscloudpopup;
options.d3.showisovolume=M.ui.showisovolumecheck;
resultfig=ea_render_view(options,M.elstruct(get(handles.patientlist,'Value')));



% --- Executes on button press in corrbutton.
function corrbutton_Callback(hObject, eventdata, handles)
% hObject    handle to corrbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[corrcl,vicorr,fccorr,vc_labels,fc_labels]=preparedataanalysis(handles);


% perform correlations:


if ~isempty(vicorr.both)
ea_corrplot([corrcl,vicorr.both],'Volume Intersections, both hemispheres',vc_labels);
end
if ~isempty(vicorr.right)
ea_corrplot([corrcl,vicorr.right],'Volume Intersections, right hemisphere',vc_labels);
end
if ~isempty(vicorr.left)
ea_corrplot([corrcl,vicorr.left],'Volume Intersections, left hemisphere',vc_labels);
end
if ~isempty(fccorr)
ea_corrplot([corrcl,fccorr],'Fibercounts',fc_labels);
end

try
    X=[corrcl,vicorr.both,vicorr.right,vicorr.left,fccorr];
    assignin('base','X',X);
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

% get user input
prompt = {'Enter title for comparison:','Enter clinical vector:'};
cvarst = inputdlg(prompt,'Enter Clinical Variable',[1;20]);

if isempty(cvarst)
    return
end
cvar = str2num(cvarst{2});

% store in model as variables
M.clinical.vars{end+1}=cvar;

M.clinical.labels{end+1}=cvarst{1};

% store model and refresh UI
setappdata(gcf,'M',M);

refreshvifc(handles);



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

% refresh analysis list
set(handles.analysislist,'String',M.patient.analysis);
try set(handles.analysislist,'Value',M.ui.listselect); end

% refresh patient list
set(handles.patientlist,'String',M.patient.list);
try set(handles.patientlist,'Value',M.ui.listselect); end


% refresh clinical list
set(handles.clinicallist,'String',M.clinical.labels);
try set(handles.clinicallist,'Value',M.ui.clinicallist); end


if get(handles.clinicallist,'Value')>length(get(handles.clinicallist,'String'))
   set(handles.clinicallist,'Value',length(get(handles.clinicallist,'String'))); 
end


% refresh selections on VI and FC Lists:
try
if M.ui.volumeintersections>length(get(handles.vilist,'String'))
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


% update UI

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

% make isomatrix button green if set.
if isfield(M,'isomatrix')
    set(handles.setisomatrixbutton,'BackgroundColor',[0.1;0.8;0.1]);
else
    set(handles.setisomatrixbutton,'BackgroundColor',[0.93,0.93,0.93]);
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
try set(handles.savefigscheck,'Value',M.ui.savefigscheck); end
try set(handles.removepriorstatscheck,'Value',M.ui.removepriorstatscheck); end


% update selectboxes:
try set(handles.elrenderingpopup,'Value',M.ui.elrendering); end
try set(handles.atlassetpopup,'Value',M.ui.atlassetpopup); end
try set(handles.fiberspopup,'Value',M.ui.fiberspopup); end
try set(handles.labelpopup,'Value',M.ui.labelpopup); end
try set(handles.modelselect,'Value',M.ui.modelselect); end





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
            M.elstruct(pt).elmodel='Medtronic 3389';
            
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
            % save M
            save([get(handles.groupdir_choosebox,'String'),'LEAD_groupanalysis.mat'],'M');
                set(gcf,'name','LEAD-DBS Groupanalysis');
            return
        end
        
        priorvilist=M.vilist;
        try
        M.vilist=M.stats(pt).ea_stats.atlases.names;
        end
        % check and compare with prior atlas intersection list.
        
        if ~isempty(priorvilist) && ~isequal(priorvilist,M.vilist)
            
            warning('Trying to analyse inhomogeneous patient group. Please re-run single subject lead analysis with patients using always the same atlas set.');
        end
        
        
        
        
        try
            priorfclist=M.fclist;
            M.fclist=ea_stats.ft(1).labels{1};
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

% save M
if ~strcmp(get(handles.groupdir_choosebox,'String'),'Choose Group Directory') % group dir still not chosen                
            save([get(handles.groupdir_choosebox,'String'),'LEAD_groupanalysis.mat'],'M');
end

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
set(handles.analysislist,'Value',M.ui.listselect);

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


[corrcl,vicorr,fccorr,vc_labels,fc_labels]=preparedataanalysis(handles);


% perform t-tests:
keyboard
if ~isempty(vicorr.both)
ea_ttest(vicorr.both(repmat(logical(corrcl),1,size(vicorr.both,2))),vicorr.both(~repmat(logical(corrcl),1,size(vicorr.both,2))),'Volume Intersections, both hemispheres',vc_labels);
end
if ~isempty(vicorr.right)
ea_ttest(vicorr.right(repmat(logical(corrcl),1,size(vicorr.right,2))),vicorr.right(~repmat(logical(corrcl),1,size(vicorr.right,2))),'Volume Intersections, right hemisphere',vc_labels);
end
if ~isempty(vicorr.left)
ea_ttest(vicorr.left(repmat(logical(corrcl),1,size(vicorr.left,2))),vicorr.left(~repmat(logical(corrcl),1,size(vicorr.left,2))),'Volume Intersections, left hemisphere',vc_labels);
end

if ~isempty(fccorr)
ea_ttest(fccorr(repmat(logical(corrcl),1,size(fccorr,2))),fccorr(~repmat(logical(corrcl),1,size(fccorr,2))),'Fibercounts',vc_labels);
end





function [corrcl,vicorr,fccorr,vc_labels,fc_labels]=preparedataanalysis(handles)

M=getappdata(gcf,'M');


%M.stats(get(handles.vilist,'Value'))

% Get volume intersections:
vicnt=1; ptcnt=1;
vicorr_right=[]; vicorr_left=[]; vicorr_both=[];
vc_labels={};
for vi=get(handles.vilist,'Value') % get volume interactions for each patient from stats
    for pt=get(handles.patientlist,'Value')
        for vat=M.stats(pt).ea_stats.vatanalyses(end+M.patient.analysis(pt)).vatsused;
            if M.stats(pt).ea_stats.vat(vat).Side==1 % right hemisphere
                try
                    pval=vicorr_right(ptcnt,vicnt); % prior val ? since there might be more than one VAT on one side
                catch
                    pval=0;
                end
                vicorr_right(ptcnt,vicnt)=pval+M.stats(pt).ea_stats.vat(vat).AtlasIntersection(vi);
            elseif M.stats(pt).ea_stats.vat(vat).Side==2 % left hemisphere
                try
                    pval=vicorr_left(ptcnt,vicnt); % prior val ? since there might be more than one VAT on one side
                catch
                    pval=0;
                end
                
                vicorr_left(ptcnt,vicnt)=pval+M.stats(pt).ea_stats.vat(vat).AtlasIntersection(vi);
          
            end
            try
                pval=vicorr_both(ptcnt,vicnt); % prior val ? since there might be more than one VAT on one side
            catch
                pval=0;
            end
            vicorr_both(ptcnt,vicnt)=pval+M.stats(pt).ea_stats.vat(vat).AtlasIntersection(vi);

        end
        
        
        % check if all three values have been served. if not, set to zero
        % (e.g. if there was no stimulation at all on one hemisphere, this
        % could happen.
        if size(vicorr_right,1)<ptcnt; vicorr_right(ptcnt,vicnt)=0; end
        if size(vicorr_left,1)<ptcnt; vicorr_left(ptcnt,vicnt)=0; end
        if size(vicorr_both,1)<ptcnt; vicorr_both(ptcnt,vicnt)=0; end
        
        ptcnt=ptcnt+1;
        
    end
    vc_labels{end+1}=M.stats(pt).ea_stats.atlases.names{vi};

    ptcnt=1;
    vicnt=vicnt+1;
end


% Get fibercounts (here no difference between sides. fts are always calculated for the whole brain):

fccnt=1; ptcnt=1;
fccorr=[];
fc_labels={};
for fc=get(handles.fclist,'Value') % get volume interactions for each patient from stats
    for pt=get(handles.patientlist,'Value')
        fibersused=M.stats(pt).ea_stats.vatanalyses(end+M.patient.analysis(pt)).fibersused;
            
                try
                    pval=fccorr(ptcnt,fccnt); % prior val ? since there might be more than one VAT on one side
                catch
                    pval=0;
                end
                fccorr(ptcnt,fccnt)=pval+M.stats(pt).ea_stats.ft(fibersused).fibercounts{1}(fc);
            
        
        ptcnt=ptcnt+1;
        
    end
    ptcnt=1;
    fccnt=fccnt+1;
    fc_labels{end+1}=M.stats(pt).ea_stats.ft(fibersused).labels{1}{fc};

end


% prepare outputs:

vicorr.both=vicorr_both;
vicorr.left=vicorr_left;
vicorr.right=vicorr_right;


% clinical vector:
corrcl=M.clinical.vars{get(handles.clinicallist,'Value')};

clinstrs=get(handles.clinicallist,'String');
vc_labels=[clinstrs(get(handles.clinicallist,'Value')),vc_labels]; % add name of clinical vector to labels
fc_labels=[clinstrs(get(handles.clinicallist,'Value')),fc_labels]; % add name of clinical vector to labels


% --- Executes on selection change in analysislist.
function analysislist_Callback(hObject, eventdata, handles)
% hObject    handle to analysislist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns analysislist contents as cell array
%        contents{get(hObject,'Value')} returns selected item from analysislist
M=getappdata(gcf,'M');

M.ui.listselect=get(handles.analysislist,'Value');

set(handles.patientlist,'Value',M.ui.listselect);
set(handles.grouplist,'Value',M.ui.listselect);

setappdata(gcf,'M',M);
refreshvifc(handles);

% --- Executes during object creation, after setting all properties.
function analysislist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to analysislist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in plusanalysisbutton.
function plusanalysisbutton_Callback(hObject, eventdata, handles)
% hObject    handle to plusanalysisbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
M=getappdata(gcf,'M');
if M.patient.analysis(get(handles.patientlist,'Value'))<0
M.patient.analysis(get(handles.patientlist,'Value'))=M.patient.analysis(get(handles.patientlist,'Value'))+1;
end
setappdata(gcf,'M',M);
refreshvifc(handles);

% --- Executes on button press in minusanalysisbutton.
function minusanalysisbutton_Callback(hObject, eventdata, handles)
% hObject    handle to minusanalysisbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
M=getappdata(gcf,'M');
M.patient.analysis(get(handles.patientlist,'Value'))=M.patient.analysis(get(handles.patientlist,'Value'))-1;
setappdata(gcf,'M',M);
refreshvifc(handles);


% --- Executes on button press in reviewvarbutton.
function reviewvarbutton_Callback(hObject, eventdata, handles)
% hObject    handle to reviewvarbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
M=getappdata(gcf,'M');

% get user input
prompt = {'Enter title for comparison:','Enter clinical vector:'};
defans{1}=M.clinical.labels{get(handles.clinicallist,'Value')};
defans{2}=num2str(M.clinical.vars{get(handles.clinicallist,'Value')});
cvarst = inputdlg(prompt,'Enter Clinical Variable',[1;20],defans);
if ~isempty(cvarst) % cancel.
cvar = str2num(cvarst{2});

% store in model as variables
M.clinical.vars{get(handles.clinicallist,'Value')}=cvar;

M.clinical.labels{get(handles.clinicallist,'Value')}=cvarst{1};

% store model and refresh UI
setappdata(gcf,'M',M);

refreshvifc(handles);
end

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
M.patient.analysis=M.patient.analysis(ix);
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
M.patient.analysis=M.patient.analysis(ix);
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
    options.patientname=M.patient.list{pt}(slashes(end)+1:end);
    
    options.root=M.patient.list{pt}(1:slashes(end));

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
    
    
    if ~exist(options.root,'file') % data is not there. Process in tmp-dir.
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
        
    else
        processlocal=0;
    end
    
    if get(handles.removepriorstatscheck,'Value')
        delete([options.root,options.patientname,filesep,'ea_stats.mat']);
    end
    
    % Step 1: Re-calculate closeness to subcortical atlases.
    resultfig=ea_render_view(options);
    % save scene as matlab figure
    
    
 
    
    % Step 2: Re-calculate Fibertracts / VAT
    setappdata(resultfig,'stimparams',M.stimparams(pt,:));
    ea_showfibres_volume(resultfig,options);
    if get(handles.savefigscheck,'Value')
    saveas(resultfig,[options.root,options.patientname,filesep,'LEAD_scene.fig']);
    end
    close(resultfig);
    
     if processlocal % gather stats and recos to M
         load([M.ui.groupdir,'tmp',filesep,'ea_stats']);
         load([M.ui.groupdir,'tmp',filesep,'ea_reconstruction']);
         M.stats(pt).ea_stats=ea_stats;
        M.elstruct(pt).coords_mm=coords_mm;
        M.elstruct(pt).trajectory=trajectory;
        setappdata(gcf,'M',M);
        
        save([M.ui.groupdir,'LEAD_groupanalysis.mat'],'M');
        
        rmdir([M.ui.groupdir,'tmp'],'s');
        
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


% --- Executes on button press in savefigscheck.
function savefigscheck_Callback(hObject, eventdata, handles)
% hObject    handle to savefigscheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of savefigscheck
M=getappdata(gcf,'M');
M.ui.savefigscheck=get(handles.savefigscheck,'Value');


setappdata(gcf,'M',M);
refreshvifc(handles);

% --- Executes on button press in removepriorstatscheck.
function removepriorstatscheck_Callback(hObject, eventdata, handles)
% hObject    handle to removepriorstatscheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of removepriorstatscheck
M=getappdata(gcf,'M');
M.ui.removepriorstatscheck=get(handles.removepriorstatscheck,'Value');


setappdata(gcf,'M',M);
refreshvifc(handles);

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
lead;


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
M.patient.analysis=[];

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
M.ui.savefigscheck=0;
M.ui.showisovolumecheck=0;
M.ui.isovscloudpopup=1;
M.ui.removepriorstatscheck=0;
M.ui.atlassetpopup=1;
M.ui.fiberspopup=3;
M.ui.labelpopup=1;
M.ui.modelselect=1;
M.ui.volumeintersections=1;
M.ui.fibercounts=1;
M.ui.elrendering=1;




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
    Urproperties(pt) =    PropertyGridField('double', pU{1}(pt,:), ...
        'Category', 'Voltage Right', ...
        'DisplayName', ptname, ...
        'Description', 'Double Matrix.');
    Ulproperties(pt) =    PropertyGridField('double', pU{2}(pt,:), ...
        'Category', 'Voltage Left', ...
        'DisplayName', ptname, ...
        'Description', 'Double Matrix.');
    Irproperties(pt) =    PropertyGridField('double', pIm{1}(pt,:), ...
        'Category', 'Impedance Right', ...
        'DisplayName', ptname, ...
        'Description', 'Double Matrix.');
    Ilproperties(pt) =    PropertyGridField('double', pIm{2}(pt,:), ...
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
    options.elmodel=M.elstruct(pt).elmodel;
    options=ea_resolve_elspec(options);
    options.prefs=ea_prefs(options.patientname);
    options.d3.verbose='off';
   % set stimparams based on values provided by user

   
   for side=1:2
       M.stimparams(pt,side).U=gU{side}.Properties(pt).Value;
       M.stimparams(pt,side).Im=gI{side}.Properties(pt).Value;
       
       M.stimparams(pt,side).usefiberset=get(handles.fiberspopup,'String');
       
       M.stimparams(pt,side).usefiberset=M.stimparams(pt,side).usefiberset{M.ui.fiberspopup};
       
       M.stimparams(pt,side).labelatlas={options.labelatlas};
       M.stimparams(pt,side).showfibers=1;
       M.stimparams(pt,side).fiberthresh=1;
       
       M.stimparams(pt,side).VAT.VAT=ea_genvat(M.elstruct(pt).coords_mm,M.stimparams(pt,side),side,options);
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

% --- Executes on button press in setisomatrixbutton.
function setisomatrixbutton_Callback(hObject, eventdata, handles)
% hObject    handle to setisomatrixbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
M=getappdata(gcf,'M');

try
    uicell=inputdlg('Enter Variable name for Isomatrix','Enter Isomatrix...',1);
    if isempty(uicell); return; end
    M.isomatrix=evalin('base',uicell{1});
catch
    try
    M=rmfield(M,'isomatrix');
    end
    warning('Isomatrix could not be evaluated. Please Try again.');
end

% check if isomatrix needs to be expanded:



iso=M.isomatrix;
if ~iscell(iso)
    if min(size(iso))==1 && length(size(iso))==2 % single vector
    
        for side=1:2
           stimmat{side}=cat(1,M.stimparams(:,1).U); 
            stimmat{side}=bsxfun(@times,stimmat{side}>0,iso);
        end
        
    end
end
M.isomatrix=stimmat;


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
    for side=1:size(M.stimparams,2)
   pU{side}(pt,:)=M.stimparams(pt,side).U;
   pIm{side}(pt,:)=M.stimparams(pt,side).Im;
    end
end
catch
    pU{1}=zeros(length(M.patient.list),size(M.elstruct(1).coords_mm{1},1));
    pIm{1}=repmat(1000,length(M.patient.list),size(M.elstruct(1).coords_mm{1},1));
end


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

   
   for side=1:2
       M.stimparams(pt,side).U=uidata.U{side}(pt,1:options.elspec.numel);
       M.stimparams(pt,side).Im=uidata.Im{side}(pt,1:options.elspec.numel);
       
       M.stimparams(pt,side).usefiberset=get(handles.fiberspopup,'String');
       
       M.stimparams(pt,side).usefiberset=M.stimparams(pt,side).usefiberset{get(handles.fiberspopup,'Value')};
       M.stimparams(pt,side).labelatlas={options.labelatlas};
       M.stimparams(pt,side).showfibers=1;
       M.stimparams(pt,side).fiberthresh=1;
       
       M.stimparams(pt,side).VAT.VAT=ea_genvat(M.elstruct(pt).coords_mm,M.stimparams(pt,side),side,options);
       M.stimparams(pt,side).showconnectivities=1;
       M.elstruct(pt).activecontacts{side}=find(M.stimparams(pt,side).U);
   end
    
   
    
end
setappdata(gcf,'M',M);
refreshvifc(handles);
