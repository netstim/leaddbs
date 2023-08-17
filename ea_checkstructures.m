function varargout = ea_checkstructures(varargin)
% EA_CHECKSTRUCTURES MATLAB code for ea_checkstructures.fig
%      EA_CHECKSTRUCTURES, by itself, creates a new EA_CHECKSTRUCTURES or raises the existing
%      singleton*.
%
%      H = EA_CHECKSTRUCTURES returns the handle to a new EA_CHECKSTRUCTURES or the handle to
%      the existing singleton*.
%
%      EA_CHECKSTRUCTURES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EA_CHECKSTRUCTURES.M with the given input arguments.
%
%      EA_CHECKSTRUCTURES('Property','Value',...) creates a new EA_CHECKSTRUCTURES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ea_checkstructures_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ea_checkstructures_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ea_checkstructures

% Last Modified by GUIDE v2.5 21-Sep-2021 10:48:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @ea_checkstructures_OpeningFcn, ...
    'gui_OutputFcn',  @ea_checkstructures_OutputFcn, ...
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


% --- Executes just before ea_checkstructures is made visible.
function ea_checkstructures_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ea_checkstructures (see VARARGIN)

% Choose default command line output for ea_checkstructures
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% add atlases contextmenu
options=varargin{1};
set(handles.checkstructures,'Name',['Check registration of specific structures: ',options.patientname]);
options.atlasset = 'DISTAL Minimal (Ewert 2017)'; % Force atlas to Distal for checkstructures (so it does not matter which atlas is set in the lead DBS window)

setappdata(handles.checkstructures,'options',options);
setappdata(handles.checkstructures,'hemisphere',2);
setappdata(handles.checkstructures,'offset',[0.5,0.5,0.5]);
presentfiles = fieldnames(options.subj.preopAnat);

c = uicontextmenu(handles.checkstructures);
handles.otherstructures.UIContextMenu = c;
atlases = dir(ea_space(options,'atlases'));
atlases = {atlases(cell2mat({atlases.isdir})).name};    % only keep folders
atlases = atlases(cellfun(@(x) ~strcmp(x(1),'.'), atlases));  % also remove '.', '..' and '.*' folders from dir results
atlmenu = cell(length(presentfiles),length(atlases));
warning('off');
set(handles.refinestatus,'String','Click on panels to add detail-corrections relative to the atlas of choice.');

for atl=1:length(atlases)
    if ~exist([ea_space(options,'atlases'),atlases{atl},filesep,'atlas_index.mat'],'file')
        continue
    end
    atlmenu{atl}=uimenu('Parent',c,'Label',atlases{atl});
    clear a

    a=load([ea_space(options,'atlases'),atlases{atl},filesep,'atlas_index.mat'],'structures');
    if isempty(fieldnames(a)) % old format
        disp(['Re-indexing ',atlases{atl},'...']);
        a=load([ea_space(options,'atlases'),atlases{atl},filesep,'atlas_index.mat']);
        a.structures=a.atlases.names;
        save([ea_space(options,'atlases'),atlases{atl},filesep,'atlas_index.mat'],'-struct','a','-v7.3');
    end
    a.structures=ea_rmext(a.structures);
    structmenu{atl,1}=uimenu('Parent',atlmenu{atl},'Label','ALL','Callback',{@ea_setnewatlas,options,handles});
    try
        for strct=1:length(a.structures)
            structmenu{atl,strct+1}=uimenu('Parent',atlmenu{atl},'Label',a.structures{strct},'Callback',{@ea_setnewatlas,options,handles});
        end
    catch
        keyboard
    end
end

warning('on');
axes(handles.tra);
imshow(zeros(10,10,3));
axis equal;
axis off;
axes(handles.cor);
imshow(zeros(10,10,3));
axis equal;
axis off;
axes(handles.sag);
imshow(zeros(10,10,3));
axis equal;
axis off;
drawnow
% add preop acquisitions to popup
cellentr = presentfiles;
set(handles.anat_select, 'String', cellentr);
modality = get(handles.anat_select, 'String');
modality = modality{get(handles.anat_select, 'Value')};
setappdata(handles.checkstructures, 'modality', modality);
options.prefs = ea_prefs(options.patientname);
setappdata(handles.checkstructures, 'options', options);
try
    switch(options.prefs.machine.checkreg.default)
        case 'DISTAL Minimal (Ewert 2017)@STN'
            ea_preset_stn(handles)
        case 'DISTAL Minimal (Ewert 2017)@GPi'
            ea_preset_gpi(handles)
        otherwise
            parts=ea_strsplit(options.prefs.machine.checkreg.default,'@');
            h.Parent.Label=parts{1};
            h.Label=parts{2};
            ea_setnewatlas(h,[],options,handles);
    end
catch % default (e.g. when changing to a different space
    sd=load([ea_space,'spacedef.mat']);
    defaultnucleus=sd.spacedef.defaultnucleus;
    parts=ea_strsplit(defaultnucleus,'@');
    h.Parent.Label=parts{1};
    h.Label=parts{2};
    ea_setnewatlas(h,[],options,handles);
end

% UIWAIT makes ea_checkstructures wait for user response (see UIRESUME)
% uiwait(handles.checkstructures);


function ea_preset_stn(handles)
set(handles.stn,'Value',1); set(handles.gpi,'Value',0);
stnmods={'T2w','T2starw','FGATIR','QSM'};
mods=get(handles.anat_select,'String');
[is,idx]=ismember(mods,stnmods);
if any(is) % only change modality if theres a suitable one available.
    [~,sorted]=sort(idx,'ascend');
    for mod=sorted'
        if is(mod)
            bestmod=mod;
            break
        end
    end
    set(handles.anat_select,'Value',bestmod);
    ea_setnewbackdrop(handles,1);
end
options=getappdata(handles.checkstructures,'options');
h.Parent.Label='DISTAL Minimal (Ewert 2017)';
h.Label='STN';
ea_setnewatlas(h,[],options,handles)


function ea_preset_gpi(handles)
set(handles.stn,'Value',0); set(handles.gpi,'Value',1);
gpimods={'FGATIR','IR','T1','PD','QSM'};
mods=get(handles.anat_select,'String');
[is,idx]=ismember(mods,gpimods);
if any(is) % only change modality if theres a suitable one available.
    [~,sorted]=sort(idx,'ascend');
    for mod=sorted'
        if is(mod)
            bestmod=mod;
            break
        end
    end
    set(handles.anat_select,'Value',bestmod);
    ea_setnewbackdrop(handles,1);
end
options=getappdata(handles.checkstructures,'options');

h.Parent.Label='DISTAL Minimal (Ewert 2017)';
h.Label='GPi';
ea_setnewatlas(h,[],options,handles)


function ea_setnewbackdrop(handles,dontupdate)
modality=get(handles.anat_select,'String');
modality=modality{get(handles.anat_select,'Value')};
setappdata(handles.checkstructures,'modality',modality);
if ~exist('dontupdate','var')
    options=getappdata(handles.checkstructures,'options');
    ea_updateviews(options,handles,1:3)
end


function ea_setnewhemisphere(handles,dontupdate)
hemisphere=get(handles.rh,'Value');
if ~hemisphere % LH
    hemisphere=2;
end
setappdata(handles.checkstructures,'hemisphere',hemisphere);
if ~exist('dontupdate','var')
    options=getappdata(handles.checkstructures,'options');
    ea_setnewatlas([],[],options,handles,1);
    ea_updateviews(options,handles,3)
end


function ea_setnewatlas(h,gf,options,handles,dontupdate)
if isempty(h)
    h=getappdata(handles.checkstructures,'h');
end

ea_setprefs('checkreg.default',[h.Parent.Label,'@',h.Label]);

options.atlasset=h.Parent.Label;
load([ea_space(options,'atlases'),options.atlasset,filesep,'atlas_index.mat']);
if strcmp(h.Label,'ALL')
    roi=atlases.roi;
    pixdim=atlases.pixdim;
    xyz=[];
    for i=1:numel(roi)
        try % some are empty in case of midline/mixed structures
            xyz=[xyz;roi{i}.fv.vertices];
        end
    end
    pixdim=pixdim{1};
else
    [~,idx]=ismember(h.Label,ea_rmext(atlases.names));
    roi=atlases.roi(idx,:);
    pixdim=atlases.pixdim(idx,:);
    if length(roi)>1
        xyz=[];
        for i=1:numel(roi)
            try % some are empty in case of midline/mixed structures
                xyz=[xyz;roi{i}.fv.vertices];
            end
        end
        pixdim=mean([pixdim{1};pixdim{2}]);
    else
        xyz=roi{1}.fv.vertices;
        pixdim=pixdim{1};
    end
end

mz=mean(xyz);
vmz=abs(max(xyz)-min(xyz));
hemisphere=getappdata(handles.checkstructures,'hemisphere');
if hemisphere==1
    xyz=xyz(xyz(:,1)>0,:);
elseif hemisphere==2
    xyz=xyz(xyz(:,1)<0,:);
end
mzsag=mean(xyz);
vmzsag=abs(max(xyz)-min(xyz));

setappdata(handles.checkstructures,'h',h);
setappdata(handles.checkstructures,'roi',roi);
setappdata(handles.checkstructures,'atlases',atlases);
setappdata(handles.checkstructures,'pixdim',pixdim);
setappdata(handles.checkstructures,'mz',mz);
setappdata(handles.checkstructures,'mzsag',mzsag);
setappdata(handles.checkstructures,'vmz',vmz);
setappdata(handles.checkstructures,'vmzsag',vmzsag);
setappdata(handles.checkstructures,'options',options);
if ~exist('dontupdate','var')
    ea_updateviews(options,handles,1:3)
end


function ea_updateviews(options,handles,cortrasag)
fv=getappdata(handles.checkstructures,'roi');
atlases=getappdata(handles.checkstructures,'atlases');
pixdim=getappdata(handles.checkstructures,'pixdim');
mz=getappdata(handles.checkstructures,'mz');
mzsag=getappdata(handles.checkstructures,'mzsag');
vmz=getappdata(handles.checkstructures,'vmz');
vmzsag=getappdata(handles.checkstructures,'vmzsag');

h=getappdata(handles.checkstructures,'h');
views={'tra','cor','sag'};
for cts=cortrasag
    options.d2.writeatlases=1;
    options.d2.col_overlay=0;
    options.d2.con_color=[0.8,0.1,0];
    options.d2.atlasopacity=0.2;
    options.d2.tracor=cts;
    options.d2.bbsize=(max(vmz)/1.7);
    offset=getappdata(handles.checkstructures,'offset')-0.5;
    if cts<3
        mz(ea_view2coord(cts))=mz(ea_view2coord(cts))+offset(ea_view2coord(cts)).*vmz(ea_view2coord(cts));
        options.d2.depth=mz;
    else
        mzsag(ea_view2coord(cts))=mzsag(ea_view2coord(cts))+offset(ea_view2coord(cts)).*vmzsag(ea_view2coord(cts));
        options.d2.depth=mzsag;
    end
    options.d2.showlegend=0;
    if strcmp(h.Label,'ALL')
        options.d2.showstructures=ea_rmext(atlases.names);
    else
        options.d2.showstructures={h.Label};
    end
    modality=getappdata(handles.checkstructures,'modality');

    [Vtra,Vcor,Vsag]=ea_assignbackdrop(['Patient Pre-OP (',modality,')'],options,'Patient');
    Vs={Vtra,Vcor,Vsag};
    options.sides=1;
    evalin('base','custom_cont=2;');
    contour=getappdata(handles.checkstructures,'contour');
    if isempty(contour)
        clear contour
    end
    options.d2.con_overlay=1;
    options.d2.lab_overlay=1;
    options.d2.col_overlay=0;
    [hf,img,bb,contour{cts}]=ea_writeplanes(options,options.d2.depth,options.d2.tracor,Vs{options.d2.tracor},'off',2);
    ea_delete([options.subj.subjDir, filesep, 'export', filesep, '2D', '*viewplane.txt']);
    bbs=getappdata(handles.checkstructures,'bbs');
    if isempty(bbs)
        clear bbs
    end
    bbs(cts).mm=bb;
    bbs(cts).imgdim=size(img);
    setappdata(handles.checkstructures,'bbs',bbs);
    setappdata(handles.checkstructures,'contour',contour);
    axes(handles.(views{cts}));
    voxdepth=Vs{1}.mat\[options.d2.depth,1]';
    voxz=voxdepth(ea_view2coord(cts));
    him=image(img);

    hold on
%     set(handles.(views{cts}),'ButtonDownFcn', @(h,e) ea_getmouse(handles.(views{cts}),handles,[voxz,ea_view2coord(cts),cts],'Color',[1,1,0.6],'linewidth',2));
    set(handles.(views{cts}),'xtick',[],'ytick',[],'xlabel',[],'ylabel',[]);
    set(him, 'HitTest', 'off');
end


function ea_getmouse(varargin)

hit=evalin('caller','e');
switch hit.Button
    case 1 % left click
        ea_refinestructs_freehanddraw(varargin{:});
    case 3 % right click
        ea_refinestructs_twopoints(varargin{:});
end


function coord=ea_view2coord(view)
switch view
    case 1
        coord=3;
    case 2
        coord=2;
    case 3
        coord=1;
end


% --- Outputs from this function are returned to the command line.
function varargout = ea_checkstructures_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in stn.
function stn_Callback(hObject, eventdata, handles)
% hObject    handle to stn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(hObject,'Value')
    ea_preset_stn(handles);
end


% --- Executes on button press in gpi.
function gpi_Callback(hObject, eventdata, handles)
% hObject    handle to gpi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(hObject,'Value')
    ea_preset_gpi(handles);
end


% --- Executes on button press in otherstructures.
function otherstructures_Callback(hObject, eventdata, handles)
% hObject    handle to otherstructures (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.stn,'Value',0); set(handles.gpi,'Value',0);

CurPos = get(0, 'PointerLocation');
figPos = get(gcf,'Position');
handles.otherstructures.UIContextMenu.Position = CurPos - figPos(1:2);
handles.otherstructures.UIContextMenu.Visible='on';


% --- Executes on selection change in anat_select.
function anat_select_Callback(hObject, eventdata, handles)
% hObject    handle to anat_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns anat_select contents as cell array
%        contents{get(hObject,'Value')} returns selected item from anat_select
ea_setnewbackdrop(handles);


% --- Executes during object creation, after setting all properties.
function anat_select_CreateFcn(hObject, eventdata, handles)
% hObject    handle to anat_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function traslide_Callback(hObject, eventdata, handles)
% hObject    handle to traslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
setappdata(handles.checkstructures,'offset',[get(handles.sagslide,'Value'),get(handles.corslide,'Value'),get(handles.traslide,'Value')]);
options=getappdata(handles.checkstructures,'options');
ea_updateviews(options,handles,1);


% --- Executes during object creation, after setting all properties.
function traslide_CreateFcn(hObject, eventdata, handles)
% hObject    handle to traslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function corslide_Callback(hObject, eventdata, handles)
% hObject    handle to corslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
setappdata(handles.checkstructures,'offset',[get(handles.sagslide,'Value'),get(handles.corslide,'Value'),get(handles.traslide,'Value')]);
options=getappdata(handles.checkstructures,'options');
ea_updateviews(options,handles,2);


% --- Executes during object creation, after setting all properties.
function corslide_CreateFcn(hObject, eventdata, handles)
% hObject    handle to corslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function sagslide_Callback(hObject, eventdata, handles)
% hObject    handle to sagslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
setappdata(handles.checkstructures,'offset',[get(handles.sagslide,'Value'),get(handles.corslide,'Value'),get(handles.traslide,'Value')]);
options=getappdata(handles.checkstructures,'options');
ea_updateviews(options,handles,3);


% --- Executes during object creation, after setting all properties.
function sagslide_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sagslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in lh.
function lh_Callback(hObject, eventdata, handles)
% hObject    handle to lh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of lh
set(handles.rh,'Value',0);
ea_setnewhemisphere(handles);


% --- Executes on button press in rh.
function rh_Callback(hObject, eventdata, handles)
% hObject    handle to rh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rh
set(handles.lh,'Value',0);
ea_setnewhemisphere(handles);


% --- Executes on button press in approvefiducial.
function approvefiducial_Callback(hObject, eventdata, handles)
% hObject    handle to approvefiducial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.approvefiducial,'visible','off');
set(handles.discardfiducial,'visible','off');
set(handles.refinestatus,'String','Great, correction added! Remember to re-run an ANTs-based normalization for changes to take effect.');
ea_busyaction('on',handles.checkstructures,'normcheckstructures');
linefiducial=getappdata(handles.checkstructures,'linefiducial');
clinefiducial=getappdata(handles.checkstructures,'clinefiducial');
fiducialview=getappdata(handles.checkstructures,'fiducialview');
bbs=getappdata(handles.checkstructures,'bbs'); % bounding boxes of views
thisbb=bbs(fiducialview);

if size(linefiducial,1)==3 && all(isinf(linefiducial(3,:))) % this is a two-point correction fiducial
    twopoint=1;
    clinefiducial=linefiducial(2,:);
    linefiducial=linefiducial(1,:);
else
    twopoint=0;
end



%linefiducial(:,1)=thisbb.imgdim(1)-linefiducial(:,1);
[planedim,onedim,secdim]=ea_getdims(fiducialview,1);

linefiducial(:,secdim)=thisbb.imgdim(1)-linefiducial(:,secdim);
expmm=zeros(size(linefiducial));
expmm(:,planedim)=thisbb.mm{planedim}(1); % all entries should be equal.

expmm(:,onedim)=linefiducial(:,onedim)./thisbb.imgdim(1)... % scale from 0 to 1
    *abs(diff([thisbb.mm{onedim}(1),thisbb.mm{onedim}(end)]))... % scale from 0 to fov size in mm
    +smallestentry([thisbb.mm{onedim}(1),thisbb.mm{onedim}(end)]); % shift to correct bb

expmm(:,secdim)=linefiducial(:,secdim)./thisbb.imgdim(2)... % scale from 0 to 1
    *diff([thisbb.mm{secdim}(1),thisbb.mm{secdim}(end)])... % scale from 0 to fov size in mm
    +smallestentry([thisbb.mm{secdim}(1),thisbb.mm{secdim}(end)]); % shift to correct bb

clinefiducial(:,secdim)=thisbb.imgdim(1)-clinefiducial(:,secdim);
cexpmm=zeros(size(clinefiducial));
cexpmm(:,planedim)=thisbb.mm{planedim}(1); % all entries should be equal.

cexpmm(:,onedim)=clinefiducial(:,onedim)./thisbb.imgdim(1)... % scale from 0 to 1
    *abs(diff([thisbb.mm{onedim}(1),thisbb.mm{onedim}(end)]))... % scale from 0 to fov size in mm
    +smallestentry([thisbb.mm{onedim}(1),thisbb.mm{onedim}(end)]); % shift to correct bb

cexpmm(:,secdim)=clinefiducial(:,secdim)./thisbb.imgdim(2)... % scale from 0 to 1
    *diff([thisbb.mm{secdim}(1),thisbb.mm{secdim}(end)])... % scale from 0 to fov size in mm
    +smallestentry([thisbb.mm{secdim}(1),thisbb.mm{secdim}(end)]); % shift to correct bb

uuid=getappdata(handles.checkstructures,'fidguiid');
if isempty(uuid)
    uuid=ea_generate_uuid;
    setappdata(handles.checkstructures,'fidguiid',uuid);
    tp_uuid=ea_generate_uuid;
    setappdata(handles.checkstructures,'tp_fidguiid',tp_uuid);
end

if 0 % export fiducial in template space - could be done for debugging.
    ea_mkdir([ea_space,'fiducials']);
    nii=ea_load_nii([ea_space,'t1.nii']);
    expvx=nii.mat\[expmm,ones(size(expmm,1),1)]';
    expvx=round(expvx(1:3,:))';
    nii.img(:)=0;
    nii.img(sub2ind(size(nii.img),(expvx(:,1)),(expvx(:,2)),(expvx(:,3))))=1;
    nii.fname=[uuid,'.nii'];
    ea_write_nii(nii);
end
map3d=0;
if map3d % this could be used to map in 3D instead of 2D - then could be incongruent to visualization which is in 2D
    % map points to closest point on atlas:
    atlfv=getappdata(handles.checkstructures,'roi'); % get current atlas
    allatlcoords=[];
    for entry=1:length(atlfv)
        try allatlcoords=[allatlcoords;atlfv{entry}.vertices]; end
    end
    [idx,d]=knnsearch(allatlcoords,expmm);
    cexpmm=allatlcoords(idx(d<3.5),:); % ignore fiducials further away than 3.5 mm
end

if twopoint % fiducials based on right-click two point info
    % append corrections to list.
    tp_allcexpmm=getappdata(handles.checkstructures,'tp_allcexpmm');
    tp_allexpmm=getappdata(handles.checkstructures,'tp_allexpmm');
    tp_allcexpmm=[tp_allcexpmm;cexpmm];
    tp_allexpmm=[tp_allexpmm;expmm];
    setappdata(handles.checkstructures,'tp_allcexpmm',tp_allcexpmm);
    setappdata(handles.checkstructures,'tp_allexpmm',tp_allexpmm);
else % fiducials based on left-click free hand drawings
    % append corrections to list.
    allcexpmm=getappdata(handles.checkstructures,'allcexpmm');
    allexpmm=getappdata(handles.checkstructures,'allexpmm');
    allcexpmm=[allcexpmm;cexpmm];
    allexpmm=[allexpmm;expmm];
    setappdata(handles.checkstructures,'allcexpmm',allcexpmm);
    setappdata(handles.checkstructures,'allexpmm',allexpmm);
end

ea_csremovedrawings(handles);
%ea_updateviews(options,handles,1:3)
ea_busyaction('off',handles.checkstructures,'normcheckstructures');

funs=getappdata(handles.checkstructures,'bdfuns');
set(handles.tra,'ButtonDownFcn', funs{1});
set(handles.cor,'ButtonDownFcn', funs{2});
set(handles.sag,'ButtonDownFcn', funs{3});


function v=smallestentry(ay)
ay=sort(ay,'ascend');
v=ay(1);


% --- Executes on button press in discardfiducial.
function discardfiducial_Callback(hObject, eventdata, handles)
% hObject    handle to discardfiducial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.approvefiducial,'visible','off');
set(handles.discardfiducial,'visible','off');
ea_csremovedrawings(handles);
funs=getappdata(handles.checkstructures,'bdfuns');
set(handles.tra,'ButtonDownFcn', funs{1});
set(handles.cor,'ButtonDownFcn', funs{2});
set(handles.sag,'ButtonDownFcn', funs{3});


% --- Executes when user attempts to close checkstructures.
function checkstructures_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to checkstructures (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
uuid=getappdata(handles.checkstructures,'fidguiid'); % uuid only used for freehandfiducials (not for twopoint fiducials)
if ~isempty(uuid)
    disp('Adding corrections to fiducial markers...');
    disp('Remember that changes only take effect if you rerun ANTs-based normalization!');

    options=getappdata(handles.checkstructures,'options');
    cexpmm=getappdata(handles.checkstructures,'allcexpmm');
    expmm=getappdata(handles.checkstructures,'allexpmm');

    tp_cexpmm=getappdata(handles.checkstructures,'tp_allcexpmm');
    tp_expmm=getappdata(handles.checkstructures,'tp_allexpmm');
    delete(hObject);

    directory=[options.root,options.patientname,filesep];
    ea_mkdir([directory,'fiducials']);
    ea_mkdir([directory,'fiducials',filesep,'native']);
    ea_mkdir([directory,'fiducials',filesep,ea_getspace]);

    options=ea_assignpretra(options);
    if ~isempty(cexpmm)
        % export this mapping in template space:


        if ~exist([directory,'fiducials',filesep,ea_getspace,filesep,uuid,'.nii'],'file')
            nii=ea_load_nii([ea_space,'t1.nii']);
            nii.fname=[directory,'fiducials',filesep,ea_getspace,filesep,uuid,'.nii'];
            nii.dt(1) = 16;
            nii.img(:)=0;
        else
            nii=ea_load_nii([directory,'fiducials',filesep,ea_getspace,filesep,uuid,'.nii']);
        end
        atlvx=nii.mat\[cexpmm,ones(size(cexpmm,1),1)]';
        atlvx=round(atlvx(1:3,:))';

        nii.img(sub2ind(size(nii.img),(atlvx(:,1)),(atlvx(:,2)),(atlvx(:,3))))=1;

        ea_write_nii(nii);


        % now project fids back to native space and export mapping there:
        expvx=nii.mat\[expmm,ones(size(expmm,1),1)]';
        [~,subcvx]=ea_map_coords(expvx,[ea_space,'t1.nii'],[directory,'forwardTransform'],[directory,options.prefs.prenii_unnormalized]);


        if ~exist([directory,'fiducials',filesep,'native',filesep,uuid,'.nii'],'file')
            nii=ea_load_nii([directory,options.prefs.prenii_unnormalized]);
            nii.fname=[directory,'fiducials',filesep,'native',filesep,uuid,'.nii'];
            nii.dt(1) = 16;
            nii.img(:)=0;
        else
            nii=ea_load_nii([directory,'fiducials',filesep,'native',filesep,uuid,'.nii']);
        end
        subcvx=round(subcvx(1:3,:))';
        nii.img(sub2ind(size(nii.img),(subcvx(:,1)),(subcvx(:,2)),(subcvx(:,3))))=1;
        ea_write_nii(nii);

        for pttemp=1:2
            switch pttemp
                case 1 % native
                    subdir='native';
                case 2 % template
                    subdir=ea_getspace;
            end
            [pathn,filen]=fileparts([directory,'fiducials',filesep,subdir,filesep,uuid,'.nii']);
            filen=[filen,'.nii'];
                        clear matlabbatch

            matlabbatch{1}.spm.spatial.smooth.data = {fullfile(pathn,filen)};
            matlabbatch{1}.spm.spatial.smooth.fwhm = [0.5 0.5 0.5];
            matlabbatch{1}.spm.spatial.smooth.dtype = 512;
            matlabbatch{1}.spm.spatial.smooth.im = 0;
            matlabbatch{1}.spm.spatial.smooth.prefix = 's';
            spm_jobman('run',{matlabbatch});
            movefile(fullfile(pathn,['s',filen]),fullfile(pathn,filen));
            gzip(fullfile(pathn,filen));
            delete(fullfile(pathn,filen));
        end
    end
    if ~isempty(tp_cexpmm)
        spacedef=ea_getspacedef;
        for pointpair=1:size(tp_cexpmm,1)
            % define in template:
            tps_uuid{pointpair}=ea_generate_uuid;

            ea_spherical_roi([directory,'fiducials',filesep,ea_getspace,filesep,tps_uuid{pointpair},'.nii'],tp_cexpmm(pointpair,:),ea_species_adjustsize(2),0,[ea_space,spacedef.templates{1},'.nii']);
            tfis{pointpair}=[directory,'fiducials',filesep,ea_getspace,filesep,tps_uuid{pointpair},'.nii'];

            % define in pt:
            ea_spherical_roi([directory,'fiducials',filesep,'native',filesep,tps_uuid{pointpair},'.nii'],tp_expmm(pointpair,:),ea_species_adjustsize(2),0,[directory,options.prefs.prenii_unnormalized]);
            pfis{pointpair}=[directory,'fiducials',filesep,'native',filesep,tps_uuid{pointpair},'.nii'];
        end


        if length(tfis)>1
            fguid=ea_generate_uuid;
            clear matlabbatch
            matlabbatch{1}.spm.util.imcalc.input = tfis';
            matlabbatch{1}.spm.util.imcalc.output = [fguid,'.nii'];
            matlabbatch{1}.spm.util.imcalc.outdir = {[directory,'fiducials',filesep,ea_getspace]};
            matlabbatch{1}.spm.util.imcalc.expression = 'mean(X)';
            matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
            matlabbatch{1}.spm.util.imcalc.options.dmtx = 1;
            matlabbatch{1}.spm.util.imcalc.options.mask = 0;
            matlabbatch{1}.spm.util.imcalc.options.interp = 1;
            matlabbatch{1}.spm.util.imcalc.options.dtype = 512;
            spm_jobman('run',{matlabbatch});
            smoothgzip([directory,'fiducials',filesep,ea_getspace],[fguid,'.nii']);
        else
            [pathn,filenn]=fileparts(tfis{1});
            smoothgzip(pathn,[filenn,'.nii']);
        end
        ea_delete(tfis);

        if length(pfis)>1
                        clear matlabbatch
            matlabbatch{1}.spm.util.imcalc.input = pfis';
            matlabbatch{1}.spm.util.imcalc.output = [fguid,'.nii'];
            matlabbatch{1}.spm.util.imcalc.outdir = {[directory,'fiducials',filesep,'native']};
            matlabbatch{1}.spm.util.imcalc.expression = 'mean(X)';
            matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
            matlabbatch{1}.spm.util.imcalc.options.dmtx = 1;
            matlabbatch{1}.spm.util.imcalc.options.mask = 0;
            matlabbatch{1}.spm.util.imcalc.options.interp = 1;
            matlabbatch{1}.spm.util.imcalc.options.dtype = 512;
            spm_jobman('run',{matlabbatch});
            smoothgzip([directory,'fiducials',filesep,'native'],[fguid,'.nii']);
        else
            [pathn,filenn]=fileparts(pfis{1});
            smoothgzip(pathn,[filenn,'.nii']);
        end
        ea_delete(pfis);


    end
    disp('Done.');

    % unapprove normalization since should be redone:
    directory=[options.root,options.patientname,filesep];
    approved=load([directory,'ea_coreg_approved.mat']);
    approved.(ea_stripext(options.prefs.gprenii))=0;
    save([directory,'ea_coreg_approved.mat'],'-struct','approved');

    if ismac
        system(['xattr -wx com.apple.FinderInfo "0000000000000000000C00000000000000000000000000000000000000000000" ',ea_path_helper([directory,ea_stripext(options.prefs.gprenii),'.nii'])]);
    end

    %% add methods dump:
    cits={
        'Horn, A., & Kuehn, A. A. (2015). Lead-DBS: a toolbox for deep brain stimulation electrode localizations and visualizations. NeuroImage, 107, 127?135. http://doi.org/10.1016/j.neuroimage.2014.12.002'
        'Horn, A., Li, N., Dembek, T. A., Kappel, A., Boulay, C., Ewert, S., et al. (2019). Lead-DBS v2: Towards a comprehensive pipeline for deep brain stimulation imaging. NeuroImage, 184, 293?316. http://doi.org/10.1016/j.neuroimage.2018.08.068'
        };
    ea_methods(options,['Fit to atlas structures was manually corrected using a custom-built tool implemented in Lead-DBS (www.lead-dbs.org; Horn & Kuehn 2015; Horn & Li et al. 2018).'],...
        cits);
else
    delete(hObject);
end


function smoothgzip(pathn,filen)
kernel=ea_species_adjustsize(2);
matlabbatch{1}.spm.spatial.smooth.data = {fullfile(pathn,filen)};
matlabbatch{1}.spm.spatial.smooth.fwhm = [kernel kernel kernel];
matlabbatch{1}.spm.spatial.smooth.dtype = 512;
matlabbatch{1}.spm.spatial.smooth.im = 0;
matlabbatch{1}.spm.spatial.smooth.prefix = 's';
spm_jobman('run',{matlabbatch});
movefile(fullfile(pathn,['s',filen]),fullfile(pathn,filen));
gzip(fullfile(pathn,filen));
delete(fullfile(pathn,filen));


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
