function varargout = ea_convis(varargin)
% EA_CONVIS MATLAB code for ea_convis.fig
%      EA_CONVIS, by itself, creates a new EA_CONVIS or raises the existing
%      singleton*.
%
%      H = EA_CONVIS returns the handle to a new EA_CONVIS or the handle to
%      the existing singleton*.
%
%      EA_CONVIS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EA_CONVIS.M with the given input arguments.
%
%      EA_CONVIS('Property','Value',...) creates a new EA_CONVIS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ea_convis_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ea_convis_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ea_convis

% Last Modified by GUIDE v2.5 09-Sep-2015 14:32:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @ea_convis_OpeningFcn, ...
    'gui_OutputFcn',  @ea_convis_OutputFcn, ...
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


% --- Executes just before ea_convis is made visible.
function ea_convis_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ea_convis (see VARARGIN)

set(gcf,'Name','Connectome Results');


% Choose default command line output for ea_anatomycontrol
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ea_anatomycontrol wait for user response (see UIRESUME)
% uiwait(handles.figure1);
resultfig=varargin{1};
options=varargin{2};
setappdata(gcf,'resultfig',resultfig);
setappdata(gcf,'options',options);
setappdata(resultfig,'convis',gcf);


refreshcv(handles);


function refreshcv(varargin)
handles=varargin{1};
hold=0; % do refresh seed coordinates
if nargin>1
    hold=varargin{2};
end


options=getappdata(gcf,'options');

if isempty(options)
   convis=getappdata(gcf,'convis'); 
    options=getappdata(convis,'options');
    set(0,'CurrentFigure',convis);
end

%% init/modify UI controls:

% parcellation popup:

directory=[options.root,options.patientname,filesep];

pdirs=dir([directory,'connectomics',filesep]);
cnt=1;

for pdir=1:length(pdirs)
    if pdirs(pdir).isdir && ~strcmp(pdirs(pdir).name,'.') && ~strcmp(pdirs(pdir).name,'..')
        parcs{cnt}=pdirs(pdir).name;
        cnt=cnt+1;
    end
end
if ~exist('parcs','var')
    cv_disableall(handles);
    return
end
set(handles.labelpopup,'String',parcs);
if get(handles.labelpopup,'Value')>length(get(handles.labelpopup,'String'));
    set(handles.labelpopup,'Value',length(get(handles.labelpopup,'String')));
end

selectedparc=parcs{get(handles.labelpopup,'Value')};

pdirectory=[options.root,options.patientname,filesep,'connectomics',filesep,selectedparc,filesep];

%% init matrix level controls:

pmdirs=dir([pdirectory,'*_CM.mat']);

for pmdir=1:length(pmdirs)
    [~,pmc{pmdir}]=fileparts(pmdirs(pmdir).name);
end

tcdirs=dir([pdirectory,'*_tc.mat']);

for tcdir=1:length(tcdirs)
    [~,pmc{end+1}]=fileparts(tcdirs(tcdir).name);
end

if ~exist('pmc','var')
    cv_disablemats(handles);
else
    cv_enablemats(handles);
end



set(handles.matmodality,'String',pmc);

if get(handles.matmodality,'Value')>length(get(handles.matmodality,'String'));
    set(handles.matmodality,'Value',length(get(handles.matmodality,'String')));
end

if ~isempty(strfind(pmc{get(handles.matmodality,'Value')},'_tc')) % no timeseries selected.
    cv_enabletime(handles);
else
    cv_disabletime(handles);
end

% parcellation scheme
aID = fopen([options.earoot,'templates',filesep,'labeling',filesep,selectedparc,'.txt']);
atlas_lgnd=textscan(aID,'%d %s');

% store selected parcellation in figure:
% store pV and pX in figure
pV=spm_vol([options.earoot,'templates',filesep,'labeling',filesep,selectedparc,'.nii']);
pX=spm_read_vols(pV);
setappdata(gcf,'pV',pV);
setappdata(gcf,'pX',pX);

%d=length(atlas_lgnd{1}); % how many ROI.
set(handles.matseed,'String',atlas_lgnd{2});
if get(handles.matseed,'Value')>length(get(handles.matseed,'String'));
    set(handles.matseed,'Value',length(get(handles.matseed,'String')));
end

%% init voxel level controls:
% Metric:
if exist([pdirectory,'graph'],'file')
    testits={'deg_','eig_','eff_','sfs_'};
    labelits={'degree centrality','eigenvector centrality','nodal efficiency','structure function similarity'};
    cnt=1;
    for ti=1:length(testits)
        tdir=dir([pdirectory,'graph',filesep,testits{ti},'*.nii']);
        if ~isempty(tdir)
            filesare{cnt}=testits{ti}; % present filetypes
            labelsare{cnt}=labelits{ti}; % present labelnames
            cnt=cnt+1;
        end
    end
    
    if ~isempty(labelsare)
        set(handles.voxmetric,'String',labelsare);
        if get(handles.voxmetric,'Value')>length(get(handles.voxmetric,'String'));
            set(handles.voxmetric,'Value',length(get(handles.voxmetric,'String')));
        end
        cv_enablevoxs(handles);
    else
        cv_disablevoxs(handles);
    end
else
    cv_disablevoxs(handles);
end

% Modality:
if exist('filesare','var')
    selectedmetric=get(handles.voxmetric,'Value');
    selectedprefix=filesare{selectedmetric}; % deg_, eig_, eff_ or sfs_
    
    fis=dir([pdirectory,'graph',filesep,selectedprefix,'*.nii']);
    cnt=1;
    for fi=1:length(fis)
        mods{cnt}=fis(fi).name(5:end);
        [~,mods{cnt}]=fileparts(mods{cnt}); % remove .nii extension
        cnt=cnt+1;
    end
    set(handles.voxmodality,'String',mods);
    if get(handles.voxmodality,'Value')>length(get(handles.voxmodality,'String'));
        set(handles.voxmodality,'Value',length(get(handles.voxmodality,'String')));
    end
end

%% recruit handles from prior results from figure
resultfig=getappdata(gcf,'resultfig');
matsurf=getappdata(resultfig,'matsurf');
seedsurf=getappdata(resultfig,'seedsurf');
graphsurf=getappdata(resultfig,'graphsurf');

%% delete any prior results
try delete(matsurf); end
try delete(seedsurf); end
try delete(graphsurf); end

if ~hold
    pV=getappdata(gcf,'pV');
    pX=getappdata(gcf,'pX');
    [xmm,ymm,zmm]=getcoordinates(pV,pX,get(handles.matseed,'Value'));
    set(handles.xmm,'String',num2str(xmm)); set(handles.ymm,'String',num2str(ymm)); set(handles.zmm,'String',num2str(zmm));
    set(handles.matseed,'ForegroundColor',[0,0,0]);
end

%% now show results
if get(handles.vizgraph,'Value'); % show voxel-level results
    mo_ds=get(handles.voxmodality,'String');
    mo_d=mo_ds{get(handles.voxmodality,'Value')};
    gV=spm_vol([directory,'connectomics',filesep,selectedparc,filesep,'graph',filesep,filesare{get(handles.voxmetric,'Value')},mo_d,'.nii']);
    gX=spm_read_vols(gV);
    thresh=get(handles.voxthresh,'String');
    if strcmp(thresh,'auto');
        thresh=nanmean(gX(:))+1*nanstd(gX(:));
    else
        thresh=str2double(thresh);
    end
    gX=gX>thresh;
    
    bb=[0,0,0;size(gX)];
    
    bb=map_coords_proxy(bb,pV);
    gv=cell(3,1);
    for dim=1:3
        gv{dim}=linspace(bb(1,dim),bb(2,dim),size(gX,dim));
    end
    [X,Y,Z]=meshgrid(gv{1},gv{2},gv{3});
    fv=isosurface(X,Y,Z,permute(gX,[2,1,3]),0.5); % graph_metric
    set(0,'CurrentFigure',resultfig)
    graphsurf=patch(fv,'FaceColor',options.prefs.lc.graphsurfc,'facealpha',0.7,'EdgeColor','none','facelighting','phong');
    
    setappdata(resultfig,'graphsurf',graphsurf);
end

if get(handles.vizmat,'Value'); % show matrix-level results
    %mV=pV; % duplicate labeling handle
    %mX=pX; % duplicate labeling data
    pX=round(pX);
    mms=get(handles.matmodality,'String');
    parcs=get(handles.labelpopup,'String');
    CM=load([directory,'connectomics',filesep,parcs{get(handles.labelpopup,'Value')},filesep,mms{get(handles.matmodality,'Value')}]);
    fn=fieldnames(CM);
    CM=eval(['CM.',fn{1},';']);
    
    if ~isempty(strfind(mms{get(handles.matmodality,'Value')},'_tc'))
        % timecourses selected: need to create a CM first. In this case, the variable CM is
        % not a connectivity matrix but time-courses!
        timedim=size(CM,1);
        tiwindow=get(handles.timewindow,'String');
        tiframe=get(handles.timeframe,'String');
        
        if strcmp(tiwindow,'all') || strcmp(tiframe,'all')
            % use whole CM
            CM=corrcoef(CM);
        else
            tiframe=str2double(tiframe);         tiwindow=str2double(tiwindow);
            % check if selected time window is possible:
            if (tiframe+tiwindow)>timedim || tiframe<1 % end is reached
                set(handles.timeframe,'String','1'); tiframe=1; % reset timeframe to 1
                if tiwindow>size(CM,1)
                    set(handles.timewindow,'String','1'); tiwindow=1;
                end
            end
            CM=corrcoef(CM(tiframe:tiframe+tiwindow,:)); % actual correlation
            
            if get(handles.timecircle,'Value')
                % make a step to next timeframe (prepare next iteration).
                if (tiframe+tiwindow+1)>timedim
                    set(handles.timeframe,'String','1')
                else
                    set(handles.timeframe,'String',num2str(tiframe+1))
                end
            end
        end
    end
    currentseed=get(handles.matseed,'Value');
    seedcon=CM(currentseed,:);
    thresh=get(handles.matthresh,'String');
    if strcmp(thresh,'auto');
        thresh=nanmean(seedcon)+1*nanstd(seedcon);
    else
        thresh=str2double(thresh);
    end
    tseedcon=seedcon>thresh;
    tseedcon(currentseed)=0;
    mX=ismember(round(pX),find(tseedcon));
    sX=ismember(round(pX),currentseed);
    bb=[0,0,0;size(mX)];
    
    bb=map_coords_proxy(bb,pV);
    gv=cell(3,1);
    for dim=1:3
        gv{dim}=linspace(bb(1,dim),bb(2,dim),size(mX,dim));
    end
    [X,Y,Z]=meshgrid(gv{1},gv{2},gv{3});
    
    
    fv=isosurface(X,Y,Z,permute(mX,[2,1,3]),0.5); % connected regions
    fvs=isosurface(X,Y,Z,permute(sX,[2,1,3]),0.5); % seed
    set(0,'CurrentFigure',resultfig)
    matsurf=patch(fv,'FaceColor',options.prefs.lc.matsurfc,'facealpha',0.7,'EdgeColor','none','facelighting','phong');
    
    
    %threshold seedcon
    seedcon=((seedcon-thresh)/(max(seedcon)-thresh))*255;
    
    cX=pX;
    for i=1:length(seedcon)
       cX(pX==i)=seedcon(i);
    end
    %nc=isonormals(X,Y,Z,permute(mX,[2,1,3]),matsurf);
    %nc=isocolors(X,Y,Z,permute(cX,[2,1,3]),matsurf);
    %matsurf.FaceColor='interp';

    seedsurf=patch(fvs,'FaceColor',options.prefs.lc.seedsurfc,'facealpha',0.7,'EdgeColor','none','facelighting','phong');
    
    setappdata(resultfig,'matsurf',matsurf);
    setappdata(resultfig,'seedsurf',seedsurf);
    
end

if get(handles.vizmat,'Value') && get(handles.timecircle,'Value') && strcmp(get(handles.timecircle,'Enable'),'on') % cycle over time..
    pause(0.1);
    refreshcv(handles);
end

%% save result handles to figure


function coords=map_coords_proxy(XYZ,V)

XYZ=[XYZ';ones(1,size(XYZ,1))];

coords=V.mat*XYZ;
coords=coords(1:3,:)';

%% helperfunctions to enable/disable GUI parts.

function cv_disableall(handles)
set(handles.labelpopup,'Enable','off');
set(handles.vizgraph,'Enable','off');
set(handles.voxmodality,'Enable','off');

set(handles.voxmetric,'Enable','off');
set(handles.voxthresh,'Enable','off');
set(handles.vizmat,'Enable','off');
set(handles.matmodality,'Enable','off');
set(handles.matseed,'Enable','off');
set(handles.xmm,'Enable','off');
set(handles.ymm,'Enable','off');
set(handles.zmm,'Enable','off');
set(handles.matthresh,'Enable','off');
set(handles.timewindow,'Enable','off');
set(handles.timeframe,'Enable','off');
set(handles.timecircle,'Enable','off');

function cv_disablevoxs(handles)
set(handles.vizgraph,'Enable','off');
set(handles.voxmodality,'Enable','off');
set(handles.voxmetric,'Enable','off');
set(handles.voxthresh,'Enable','off');

function cv_enablevoxs(handles)
set(handles.vizgraph,'Enable','on');
set(handles.voxmodality,'Enable','on');
set(handles.voxmetric,'Enable','on');
set(handles.voxthresh,'Enable','on');

function cv_disablemats(handles)
set(handles.matmodality,'Enable','off');
set(handles.matseed,'Enable','off');
set(handles.xmm,'Enable','off');
set(handles.ymm,'Enable','off');
set(handles.zmm,'Enable','off');
set(handles.matthresh,'Enable','off');
set(handles.timewindow,'Enable','off');
set(handles.timeframe,'Enable','off');
set(handles.timecircle,'Enable','off');

function cv_enablemats(handles)
set(handles.matmodality,'Enable','on');
set(handles.matseed,'Enable','on');
set(handles.xmm,'Enable','on');
set(handles.ymm,'Enable','on');
set(handles.zmm,'Enable','on');
set(handles.matthresh,'Enable','on');
set(handles.timewindow,'Enable','on');
set(handles.timeframe,'Enable','on');
set(handles.timecircle,'Enable','on');

function cv_disabletime(handles)
set(handles.matthresh,'Enable','off');
set(handles.timewindow,'Enable','off');
set(handles.timeframe,'Enable','off');
set(handles.timecircle,'Enable','off');

function cv_enabletime(handles)
set(handles.matthresh,'Enable','on');
set(handles.timewindow,'Enable','on');
set(handles.timeframe,'Enable','on');
set(handles.timecircle,'Enable','on');

% --- Outputs from this function are returned to the command line.
function varargout = ea_convis_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in voxmetric.
function voxmetric_Callback(hObject, eventdata, handles)
% hObject    handle to voxmetric (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns voxmetric contents as cell array
%        contents{get(hObject,'Value')} returns selected item from voxmetric
refreshcv(handles);


% --- Executes during object creation, after setting all properties.
function voxmetric_CreateFcn(hObject, eventdata, handles)
% hObject    handle to voxmetric (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in vizgraph.
function vizgraph_Callback(hObject, eventdata, handles)
% hObject    handle to vizgraph (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of vizgraph
refreshcv(handles);



function voxthresh_Callback(hObject, eventdata, handles)
% hObject    handle to voxthresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of voxthresh as text
%        str2double(get(hObject,'String')) returns contents of voxthresh as a double
refreshcv(handles);


% --- Executes during object creation, after setting all properties.
function voxthresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to voxthresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in matseed.
function matseed_Callback(hObject, eventdata, handles)
% hObject    handle to matseed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns matseed contents as cell array
%        contents{get(hObject,'Value')} returns selected item from matseed

pV=getappdata(gcf,'pV');
pX=getappdata(gcf,'pX');
[xmm,ymm,zmm]=getcoordinates(pV,pX,get(handles.matseed,'Value'));
set(handles.xmm,'String',num2str(xmm)); set(handles.ymm,'String',num2str(ymm)); set(handles.zmm,'String',num2str(zmm));
set(handles.matseed,'ForegroundColor',[0,0,0]);
refreshcv(handles);

function [xmm,ymm,zmm]=getcoordinates(pV,pX,ix)

[xx,yy,zz]=ind2sub(size(pX),find(round(pX)==ix));
XYZ=[xx,yy,zz];
centrvx=[mean(XYZ,1),1];
centrmm=pV.mat*centrvx';
xmm=centrmm(1); ymm=centrmm(2); zmm=centrmm(3);

% --- Executes during object creation, after setting all properties.
function matseed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to matseed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in matmodality.
function matmodality_Callback(hObject, eventdata, handles)
% hObject    handle to matmodality (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns matmodality contents as cell array
%        contents{get(hObject,'Value')} returns selected item from matmodality
refreshcv(handles);

% --- Executes during object creation, after setting all properties.
function matmodality_CreateFcn(hObject, eventdata, handles)
% hObject    handle to matmodality (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function matthresh_Callback(hObject, eventdata, handles)
% hObject    handle to matthresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of matthresh as text
%        str2double(get(hObject,'String')) returns contents of matthresh as a double
refreshcv(handles);

% --- Executes during object creation, after setting all properties.
function matthresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to matthresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in vizmat.
function vizmat_Callback(hObject, eventdata, handles)
% hObject    handle to vizmat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of vizmat
refreshcv(handles);


function timewindow_Callback(hObject, eventdata, handles)
% hObject    handle to timewindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of timewindow as text
%        str2double(get(hObject,'String')) returns contents of timewindow as a double
refreshcv(handles);


% --- Executes during object creation, after setting all properties.
function timewindow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to timewindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in timecircle.
function timecircle_Callback(hObject, eventdata, handles)
% hObject    handle to timecircle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of timecircle
refreshcv(handles);



function timeframe_Callback(hObject, eventdata, handles)
% hObject    handle to timeframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of timeframe as text
%        str2double(get(hObject,'String')) returns contents of timeframe as a double
refreshcv(handles);


% --- Executes during object creation, after setting all properties.
function timeframe_CreateFcn(hObject, eventdata, handles)
% hObject    handle to timeframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in voxmodality.
function voxmodality_Callback(hObject, eventdata, handles)
% hObject    handle to voxmodality (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns voxmodality contents as cell array
%        contents{get(hObject,'Value')} returns selected item from voxmodality
refreshcv(handles);


% --- Executes during object creation, after setting all properties.
function voxmodality_CreateFcn(hObject, eventdata, handles)
% hObject    handle to voxmodality (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function xmm_Callback(hObject, eventdata, handles)
% hObject    handle to xmm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xmm as text
%        str2double(get(hObject,'String')) returns contents of xmm as a double
setcoordinates(handles);
refreshcv(handles,1);

function [ix,err]=setcoordinates(handles)
% set seed selection based on manual coordinate entry.
pV=getappdata(gcf,'pV');
pX=getappdata(gcf,'pX');
xmm=str2double(get(handles.xmm,'String'));
ymm=str2double(get(handles.ymm,'String'));
zmm=str2double(get(handles.zmm,'String'));


err=0;
XYZmm=[xmm,ymm,zmm,1]';
XYZvox=pV.mat\XYZmm;
ix=0;
try
    ix=pX(round(XYZvox(1)),round(XYZvox(2)),round(XYZvox(3)));
end
if ~ix
    ix=nan;
    err=1;
end

if ~err
    set(handles.matseed,'Value',ix);
    set(handles.matseed,'ForegroundColor',[0,0,0]);
else
    set(handles.matseed,'ForegroundColor',[1,0,0]);
end

% --- Executes during object creation, after setting all properties.
function xmm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xmm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ymm_Callback(hObject, eventdata, handles)
% hObject    handle to ymm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ymm as text
%        str2double(get(hObject,'String')) returns contents of ymm as a double
setcoordinates(handles);
refreshcv(handles,1);

% --- Executes during object creation, after setting all properties.
function ymm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ymm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function zmm_Callback(hObject, eventdata, handles)
% hObject    handle to zmm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zmm as text
%        str2double(get(hObject,'String')) returns contents of zmm as a double
setcoordinates(handles);
refreshcv(handles,1);

% --- Executes during object creation, after setting all properties.
function zmm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zmm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
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
refreshcv(handles);

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
