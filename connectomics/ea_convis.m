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

% Last Modified by GUIDE v2.5 01-Mar-2018 11:37:26

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

set(hObject,'Name','Connectome Results');


% Choose default command line output for ea_anatomycontrol
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ea_anatomycontrol wait for user response (see UIRESUME)
% uiwait(handles.convis);
resultfig=varargin{1};
options=varargin{2};
setappdata(handles.convis,'resultfig',resultfig);
setappdata(handles.convis,'options',options);
setappdata(resultfig,'convis',handles.convis);

refreshcv(handles);

% Initialize the coordinate boxes.
pV=getappdata(handles.convis,'pV');
pX=getappdata(handles.convis,'pX');
[xmm,ymm,zmm]=getcoordinates(pV,pX,get(handles.matseed,'Value'));
set(handles.xmm,'String',num2str(xmm,'%.2f'));
set(handles.ymm,'String',num2str(ymm,'%.2f'));
set(handles.zmm,'String',num2str(zmm,'%.2f'));
set(handles.matseed,'ForegroundColor',[0,0,0]);


function refreshcv(varargin)
handles=varargin{1};
hold=0; % do refresh seed coordinates
apply=0;
if nargin>1
    hold=varargin{2};
end
if nargin>2
    apply=varargin{3};
end

options=getappdata(handles.convis,'options');

if isempty(options)
    convis=getappdata(handles.convis,'convis');
    options=getappdata(convis,'options');
    set(0,'CurrentFigure',convis);
else
    convis=handles.convis;
end
set(convis,'name','Processing...');
drawnow
%% init figure
[directory,pdirectory,selectedparc]=ea_cvinitgui(handles,options);

%% initialize controls
ea_initvatlevel(handles,directory,selectedparc,options);
ea_initseedlevel(handles,directory,pdirectory,selectedparc,options);
filesare=ea_initvoxellevel(handles,pdirectory);
ea_initmatlevel(handles,directory,pdirectory,selectedparc,options);

%% check if sliding window view is possible (timecourses are set as a modality somewhere)
mmc=get(handles.matmodality,'String');
vmc=get(handles.vatmodality,'String');
vomc=get(handles.voxmodality,'String');

if isempty(mmc)
   cv_disabletime(handles);
else
    mmc = mmc{get(handles.matmodality,'Value')};
    vmc = vmc{get(handles.vatmodality,'Value')};
    vomc = vomc{get(handles.voxmodality,'Value')};
    if isempty(regexp(mmc, '^Patient''s fMRI - ', 'once')) && ...
       isempty(regexp(vmc, '^Patient''s fMRI - ', 'once')) && ...
       isempty(regexp(vomc, '_fMRI$', 'once')) % no timeseries selected.

       cv_disabletime(handles);
    else
       cv_enabletime(handles);
    end
end
if apply % update elvis
    %% retrieve and delete prior results
    resultfig=ea_cvcleanup(handles);
    pV=getappdata(handles.convis,'pV');
    pX=getappdata(handles.convis,'pX');
    ea_deleteML(resultfig,handles,options);
    if ~hold
        [xmm,ymm,zmm]=getcoordinates(pV,pX,get(handles.matseed,'Value'));
        set(handles.xmm,'String',num2str(xmm,'%.2f'));
        set(handles.ymm,'String',num2str(ymm,'%.2f'));
        set(handles.zmm,'String',num2str(zmm,'%.2f'));
        set(handles.matseed,'ForegroundColor',[0,0,0]);
    end

    %% now show results
    if get(handles.vizvat,'Value') % show voxel-level results
        ea_cvshowvatresults(resultfig,pX,directory,filesare,handles,pV,selectedparc,options);
    else
        ea_deletePL(resultfig,'PL','vat');
    end

    if get(handles.vizgraph,'Value') % show voxel-level results
        ea_cvshowvoxresults(resultfig,directory,filesare,handles,pV,selectedparc,options);
    end

    if get(handles.vizmat,'Value') % show seed-based-connectivity results
        ea_cvshowseedbasedresults(resultfig,directory,pV,pX,selectedparc,handles,options);
    else
        ea_deletePL(resultfig,'PL','mat');
    end


    if get(handles.wmatedges,'Value') || get(handles.wmatnodes,'Value'); % show matrix-level results
        ea_cvshowmatresults(resultfig,directory,selectedparc,handles,options);
    end

    if (get(handles.vizmat,'Value') || get(handles.vizvat,'Value')) && get(handles.timecircle,'Value') && strcmp(get(handles.timecircle,'Enable'),'on') % cycle over time..
        pause(0.1);
        refreshcv(handles,0,1);
    end
end
set(convis,'name','Connectome Results');


function ea_deleteML(resultfig,handles,options)

ML=getappdata(resultfig,'ML');
if isempty(ML)
    return
end
try delete(ML.pedges); end
try delete(ML.pnodes); end
try delete(ML.hn); end
try delete(ML.nnodes); end
% delete matrix level stuff here.


function ea_cvshowmatresults(resultfig,directory,selectedparc,handles,options)

% show matrix level stuff here.
ML=cv_build_graph_fig(resultfig,directory,selectedparc,handles,options);
setappdata(resultfig,'ML',ML);


function ea_initmatlevel(handles,directory,pdirectory,selectedparc,options)

pmdirs=dir([pdirectory,'*_CM.mat']);
cnt=1;
modlist=cell(0);
for pmdir=1:length(pmdirs)
    [~,modlist{cnt}]=fileparts(pmdirs(pmdir).name);
    cnt=cnt+1;
end

if ~exist('modlist','var')
    %cv_disablewmats(handles);
else
    %cv_enablewmats(handles);
    modlist{end+1}='Choose...';
    set(handles.wmatmodality,'String',modlist);
end

if get(handles.wmatmodality,'Value')>length(get(handles.wmatmodality,'String'))
    set(handles.wmatmodality,'Value',length(get(handles.wmatmodality,'String')));
end

mods=get(handles.wmatmodality,'String');
if ~strcmp('Choose...',mods{get(handles.wmatmodality,'Value')})
    set(handles.chosenwmatstr,'String','');
    set(handles.chosenwmatstr,'TooltipString','');
end


function ea_cvshowvoxresults(resultfig,directory,filesare,handles,pV,selectedparc,options)

mo_ds=get(handles.voxmodality,'String');
mo_d=mo_ds{get(handles.voxmodality,'Value')};
gV=spm_vol([directory,'connectomics',filesep,selectedparc,filesep,'graph',filesep,filesare{get(handles.voxmetric,'Value')},mo_d,'.nii']);
gX=spm_read_vols(gV);
thresh=get(handles.voxthresh,'String');
if strcmp(thresh,'auto')
    thresh=nanmean(gX(:))+1*nanstd(gX(:));
else
    thresh=str2double(thresh);
end
set(handles.voxthreshis,'String',num2str(thresh));
graphsurf=ea_showconnectivitypatch(resultfig,gV,gX,thresh,[],[],[],1,0);

setappdata(resultfig,'graphsurf',graphsurf);


function ea_cvshowseedbasedresults(resultfig,directory,pV,pX,selectedparc,handles,options)
%mV=pV; % duplicate labeling handle
%mX=pX; % duplicate labeling data

% determine if CM/TC or fiberset is selected
matmodality=get(handles.matmodality,'String');
matmodality=matmodality{get(handles.matmodality,'Value')};
if regexp(matmodality, '^Patient''s fMRI - ')
    ea_deletePL(resultfig,'PL','mat');
    ea_cvshowmatresultsCMTC(resultfig,directory,pV,pX,handles,options);
else % use fiberset

    % fibers filename
    switch matmodality
        case 'Patient''s fiber tracts'
            fibersfile=[directory,'connectomes',filesep,'dMRI',filesep,options.prefs.FTR_normalized];
        otherwise
            fibersfile=[ea_getconnectomebase('dmri'),matmodality,filesep,'data.mat'];
    end

    % seed filename
    seed=ea_load_nii([ea_space(options,'labeling'),selectedparc,'.nii']);
    % delete everything but set selected parcellation to 1.
    oseed=seed.img;
    seed.img(:)=0;
    seed.img(round(oseed)==get(handles.matseed,'Value'))=1;

    targetsfile=ea_load_nii([ea_space(options,'labeling'),selectedparc,'.nii']);
    targetsfile.img(round(targetsfile.img)==get(handles.matseed,'Value'))=0;
    thresh=get(handles.matthresh,'String');
    options.writeoutstats=0;
    options.writeoutpm = 0;


    [changedstates,ret]=ea_checkfschanges(resultfig,fibersfile,seed,targetsfile,thresh,'mat');

    if ~ret % something has changed since last time.
        ea_deletePL(resultfig,'PL','mat');

        [~,thresh]=ea_cvshowfiberconnectivities(resultfig,fibersfile,seed,targetsfile,thresh,1,options,'',changedstates,'mat',get(handles.vizmat_regs,'Value'),get(handles.vizmat_labs,'Value'));
        set(handles.matthreshis,'String',num2str(thresh));
    end
end


function ea_cvshowmatresultsCMTC(resultfig,directory,pV,pX,handles,options)

pX = round(pX);
mods = get(handles.matmodality,'String');
mod = mods{get(handles.matmodality,'Value')};
rest_CM = [strrep(mod, 'Patient''s fMRI - ', ''), '_fMRI_CM'];
parcs = get(handles.labelpopup,'String');
parc = parcs{get(handles.labelpopup,'Value')};
CM=load([directory,'connectomics',filesep,parc,filesep,rest_CM]);
fn=fieldnames(CM);
CM=eval(['CM.',fn{1},';']);

if get(handles.timecircle, 'Value')
    % timecourses selected: need to create a CM first. In this case, the variable CM is
    % not a connectivity matrix but time-courses!

    timedim=size(CM,1);
    tiwindow=get(handles.timewindow,'String');
    tiframe=get(handles.timeframe,'String');

    if strcmp(tiwindow,'all') || strcmp(tiframe,'all')
        % use whole CM
        CM=corrcoef(CM);
    else
        tiframe=str2double(tiframe);
        tiwindow=str2double(tiwindow);
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
if strcmp(thresh,'auto')
    thresh=nanmean(seedcon)+1*nanstd(seedcon);
else
    thresh=str2double(thresh);
end
tseedcon=seedcon;
tseedcon(tseedcon<thresh)=0;
tseedcon(currentseed)=0;

mX=pX;
for cs=1:length(tseedcon) % assign each voxel of the corresponding cluster with the entries in tseedcon. Fixme, this should be doable wo forloop..
    mX(ismember(round(pX),cs))=tseedcon(cs);
end
sX=ismember(round(pX),currentseed);
set(handles.matthreshis,'String',num2str(thresh));
matsurf=ea_showconnectivitypatch(resultfig,pV,mX,thresh,[],[],[],1,0);
seedsurf=ea_showseedpatch(resultfig,pV,sX,options);

setappdata(resultfig,'matsurf',matsurf);
setappdata(resultfig,'seedsurf',seedsurf);


function [directory,pdirectory,selectedparc]=ea_cvinitgui(handles,options)
%% init/modify UI controls:

% parcellation popup:
pdirs=dir([ea_space(options,'labeling'),'*.nii']);
cnt=1;

for pdir=1:length(pdirs)
    [~,parcs{cnt}]=fileparts(pdirs(pdir).name);
    cnt=cnt+1;
end

set(handles.labelpopup,'String',parcs);
if get(handles.labelpopup,'Value')>length(get(handles.labelpopup,'String'))
    set(handles.labelpopup,'Value',length(get(handles.labelpopup,'String')));
end

selectedparc=parcs{get(handles.labelpopup,'Value')};
pdirectory=[options.root,options.patientname,filesep,'connectomics',filesep,selectedparc,filesep];

% check if selected parcellation exists for patient data as well:
directory=[options.root,options.patientname,filesep];

pdirs=dir([directory,'connectomics',filesep]);

cnt=1;
parcs=cell(1);
parcs{1}='';
for pdir=1:length(pdirs)
    if pdirs(pdir).isdir && ~strcmp(pdirs(pdir).name,'.') && ~strcmp(pdirs(pdir).name,'..')
        parcs{cnt}=pdirs(pdir).name;
        cnt=cnt+1;
    end
end

if ~ismember({selectedparc},parcs)
    cv_enablevats(handles);
end


function ea_initseedlevel(handles,directory,pdirectory,selectedparc,options)

%% init matrix level controls:
modlist=ea_genmodlist(directory,selectedparc,options);
% % dMRI:
% cnt=1;
% % check if pat-specific fibertracts are present:
% if exist([directory,'connectomes',filesep,'dMRI',filesep,options.prefs.FTR_normalized],'file');
%     modlist{cnt}='Patient''s fiber tracts';
%     cnt=cnt+1;
% end
%
% % check for canonical fiber sets
% fdfibs=dir([ea_getconnectomebase('dmri'),'*.mat']);
% for fdf=1:length(fdfibs)
%     [~,fn]=fileparts(fdfibs(fdf).name);
%     modlist{cnt}=fn;
%     cnt=cnt+1;
% end
%
%
% pmdirs=dir([pdirectory,'*_CM.mat']);
%
% for pmdir=1:length(pmdirs)
%     [~,modlist{end+1}]=fileparts(pmdirs(pmdir).name);
% end
%
% tcdirs=dir([pdirectory,'*_tc.mat']);
%
% for tcdir=1:length(tcdirs)
%     [~,modlist{end+1}]=fileparts(tcdirs(tcdir).name);
% end

if ~exist('modlist','var')
    cv_disablevats(handles);
    cv_disablemats(handles);
else
    cv_enablevats(handles);
    cv_enablemats(handles);
end

set(handles.matmodality,'String',modlist);

if get(handles.matmodality,'Value')>length(get(handles.matmodality,'String'))
    set(handles.matmodality,'Value',length(get(handles.matmodality,'String')));
end

% parcellation scheme
aID = fopen([ea_space(options,'labeling'),selectedparc,'.txt']);
atlas_lgnd=textscan(aID,'%d %s');

% store selected parcellation in figure:
% store pV and pX in figure
pV=spm_vol([ea_space(options,'labeling'),selectedparc,'.nii']);
pX=spm_read_vols(pV);
setappdata(handles.convis,'pV',pV);
setappdata(handles.convis,'pX',pX);

%d=length(atlas_lgnd{1}); % how many ROI.

mm=get(handles.matmodality,'String');

set(handles.matseed,'String',[atlas_lgnd{2}]);

if get(handles.matseed,'Value')>length(get(handles.matseed,'String'));
    set(handles.matseed,'Value',length(get(handles.matseed,'String')));
end


function [filesare,labelsare,mods]=ea_initvoxellevel(handles,pdirectory)

%% init voxel level controls:
% Metric:
filesare{1}='';
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
        if get(handles.voxmetric,'Value')>length(get(handles.voxmetric,'String'))
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
    mods{1}='';
    for fi=1:length(fis)
        mods{cnt}=fis(fi).name(5:end);
        [~,mods{cnt}]=fileparts(mods{cnt}); % remove .nii extension
        cnt=cnt+1;
    end
    set(handles.voxmodality,'String',mods);
    if get(handles.voxmodality,'Value')>length(get(handles.voxmodality,'String'))
        set(handles.voxmodality,'Value',length(get(handles.voxmodality,'String')));
    end
end


function ea_initvatlevel(handles,directory,selectedparc,options)

modlist=ea_genmodlist(directory,selectedparc,options);

%% VATs:

vdirs=dir([directory,'stimulations',filesep,ea_nt(options)]);
cnt=1;
vdicell=cell(0); %empty cell
for vdir=1:length(vdirs)
    if vdirs(vdir).isdir && ~strcmp(vdirs(vdir).name(1),'.')
        vdicell{cnt}=vdirs(vdir).name;
        cnt=cnt+1;
    end
end

%% set popup strings
if isempty(modlist)
    set(handles.vatmodality,'String','No connectivity data found...');
else
    set(handles.vatmodality,'String',modlist);
end

%% correct for wrong selections in popup menus.
if get(handles.vatmodality,'Value')>length(modlist) % probably user has changed parcellation..
    set(handles.vatmodality,'Value',1);
end
if get(handles.vatseed,'Value')>length(get(handles.vatseed,'String'))
    set(handles.vatseed,'Value',1);
end

%% handle empty popup cases:
if isempty(vdicell)
    set(handles.vatseed,'String','No stimulation found...');
else
    set(handles.vatseed,'String',vdicell);
end

if isempty(vdicell) || isempty(modlist)
    cv_disablevats(handles)
else
    cv_enablevats(handles)

    %% check if left/right VATs are present
    stimfolder=vdicell{get(handles.vatseed,'Value')};
    vatdir=dir([directory,'stimulations',filesep,ea_nt(options),stimfolder,filesep,'*.nii']);
    for vt=1:length(vatdir)
        vatcell{vt}=vatdir(vt).name;
    end

    set(handles.rvatcheck,'Enable', ea_getonofftruefalse(ismember('vat_right.nii',vatcell)));
    set(handles.lvatcheck,'Enable', ea_getonofftruefalse(ismember('vat_left.nii',vatcell)));
    if strcmp(get(handles.rvatcheck,'Enable'),'off')
        set(handles.rvatcheck,'Value', 0);
    end
    if strcmp(get(handles.lvatcheck,'Enable'),'off')
        set(handles.lvatcheck,'Value', 0);
    end
end


function oo=ea_getonofftruefalse(tf)
oo='off';
if tf
    oo='on';
end

function resultfig=ea_cvcleanup(handles)
%% recruit handles from prior results from figure
resultfig=getappdata(handles.convis,'resultfig');
matsurf=getappdata(resultfig,'matsurf');
vatsurf=getappdata(resultfig,'vatsurf');
seedsurf=getappdata(resultfig,'seedsurf');
vatseedsurf=getappdata(resultfig,'vatseedsurf');
graphsurf=getappdata(resultfig,'graphsurf');
%% delete any prior results

try delete(matsurf); catch
    for l=1:length(matsurf)
        try delete(matsurf{l}); end
    end
end
try delete(vatseedsurf); catch
    for l=1:length(vatseedsurf)
        try delete(vatseedsurf{l}); end
    end
end
try delete(vatsurf); catch
    for l=1:length(vatsurf)
        try delete(vatsurf{l}); end
    end
end
try delete(seedsurf); end
try delete(graphsurf); end

%% fiber results are only cleaned if really triggered (since they take quite long).

%% helperfunctions to enable/disable GUI parts.

function cv_enablevats(handles)
% set(handles.labelpopup,'Enable','on');
% set(handles.vizgraph,'Enable','off');
% set(handles.voxmodality,'Enable','off');
%
% set(handles.voxmetric,'Enable','off');
% set(handles.voxthresh,'Enable','off');
% set(handles.vizmat,'Enable','off');
% set(handles.matmodality,'Enable','off');
% set(handles.matseed,'Enable','off');
% set(handles.xmm,'Enable','off');
% set(handles.ymm,'Enable','off');
% set(handles.zmm,'Enable','off');
% set(handles.matthresh,'Enable','off');
% set(handles.timewindow,'Enable','off');
% set(handles.timeframe,'Enable','off');
% set(handles.timecircle,'Enable','off');

set(handles.vatmodality,'Enable','on');
set(handles.vatseed,'Enable','on');
set(handles.vatthresh,'Enable','on');
set(handles.vatthreshis,'Enable','on');
set(handles.lvatcheck,'Enable','on');
set(handles.rvatcheck,'Enable','on');
set(handles.vizvat,'Enable','on');


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


function cv_disablewmats(handles)
set(handles.vizwmat,'Enable','off');
set(handles.wmatmodality,'Enable','off');
set(handles.wmatnodes,'Enable','off');
set(handles.wmatedges,'Enable','off');
set(handles.wmatthresh,'Enable','off');


function cv_enablewmats(handles)
set(handles.vizwmat,'Enable','on');
set(handles.wmatmodality,'Enable','on');
set(handles.wmatnodes,'Enable','on');
set(handles.wmatedges,'Enable','on');
set(handles.wmatthresh,'Enable','on');


function cv_enablemats(handles)
set(handles.labelpopup,'Enable','on');
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
%set(handles.matthresh,'Enable','off');
set(handles.timewindow,'Enable','off');
set(handles.timeframe,'Enable','off');
set(handles.timecircle,'Enable','off');


function cv_enabletime(handles)
%set(handles.matthresh,'Enable','on');
set(handles.timewindow,'Enable','on');
set(handles.timeframe,'Enable','on');
set(handles.timecircle,'Enable','on');


function cv_disablevats(handles)
set(handles.vizvat,'Enable','off');
set(handles.vatmodality,'Enable','off');
set(handles.vatseed,'Enable','off');
set(handles.lvatcheck,'Enable','off');
set(handles.rvatcheck,'Enable','off');
set(handles.vatthresh,'Enable','off');


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

pV=getappdata(handles.convis,'pV');
pX=getappdata(handles.convis,'pX');
[xmm,ymm,zmm]=getcoordinates(pV,pX,get(handles.matseed,'Value'));
set(handles.xmm,'String',num2str(xmm,'%.2f'));
set(handles.ymm,'String',num2str(ymm,'%.2f'));
set(handles.zmm,'String',num2str(zmm,'%.2f'));
set(handles.matseed,'ForegroundColor',[0,0,0]);
refreshcv(handles);

function [xmm,ymm,zmm]=getcoordinates(pV,pX,ix)
if ix<1
    xmm=nan; ymm=nan; zmm=nan;
    return
end
[xx,yy,zz]=ind2sub(size(pX),find(round(pX)==ix));
XYZ=[xx,yy,zz];
centrvx=[mean(XYZ,1),1];
centrmm=pV.mat*centrvx';
xmm=centrmm(1);
ymm=centrmm(2);
zmm=centrmm(3);


function [ix,err]=setcoordinates(handles)
% set seed selection based on manual coordinate entry.
pV=getappdata(handles.convis,'pV');
pX=getappdata(handles.convis,'pX');
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


% --- Executes on selection change in vatseed.
function vatseed_Callback(hObject, eventdata, handles)
% hObject    handle to vatseed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns vatseed contents as cell array
%        contents{get(hObject,'Value')} returns selected item from vatseed

refreshcv(handles);


% --- Executes during object creation, after setting all properties.
function vatseed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vatseed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function vatthresh_Callback(hObject, eventdata, handles)
% hObject    handle to vatthresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vatthresh as text
%        str2double(get(hObject,'String')) returns contents of vatthresh as a double

refreshcv(handles);


% --- Executes during object creation, after setting all properties.
function vatthresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vatthresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in vizvat.
function vizvat_Callback(hObject, eventdata, handles)
% hObject    handle to vizvat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of vizvat
refreshcv(handles);


% --- Executes on selection change in vatmodality.
function vatmodality_Callback(hObject, eventdata, handles)
% hObject    handle to vatmodality (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns vatmodality contents as cell array
%        contents{get(hObject,'Value')} returns selected item from vatmodality
refreshcv(handles);


% --- Executes during object creation, after setting all properties.
function vatmodality_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vatmodality (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in lvatcheck.
function lvatcheck_Callback(hObject, eventdata, handles)
% hObject    handle to lvatcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of lvatcheck
refreshcv(handles);

% --- Executes on button press in rvatcheck.
function rvatcheck_Callback(hObject, eventdata, handles)
% hObject    handle to rvatcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rvatcheck
refreshcv(handles);


function wmatthresh_Callback(hObject, eventdata, handles)
% hObject    handle to wmatthresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of wmatthresh as text
%        str2double(get(hObject,'String')) returns contents of wmatthresh as a double
refreshcv(handles);


% --- Executes during object creation, after setting all properties.
function wmatthresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wmatthresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in wmatmodality.
function wmatmodality_Callback(hObject, eventdata, handles)
% hObject    handle to wmatmodality (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns wmatmodality contents as cell array
%        contents{get(hObject,'Value')} returns selected item from wmatmodality

chosenmod=get(hObject,'String');
if strcmp('Choose...',chosenmod{get(hObject,'Value')})
    [fn,pth]=uigetfile('*.mat','Choose connectivity matrix...');
    if fn
        set(handles.chosenwmatstr,'String',[pth,filesep,fn]);
        set(handles.chosenwmatstr,'TooltipString',[pth,filesep,fn]);
    end
else
    set(handles.chosenwmatstr,'String','');
    set(handles.chosenwmatstr,'TooltipString','');
end


% --- Executes during object creation, after setting all properties.
function wmatmodality_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wmatmodality (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in wmatnodes.
function wmatnodes_Callback(hObject, eventdata, handles)
% hObject    handle to wmatnodes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of wmatnodes
refreshcv(handles);


% --- Executes on button press in wmatedges.
function wmatedges_Callback(hObject, eventdata, handles)
% hObject    handle to wmatedges (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of wmatedges
refreshcv(handles);


% --- Executes on button press in apply.
function apply_Callback(hObject, eventdata, handles)
% hObject    handle to apply (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
refreshcv(handles,0,1);


% --- Executes on button press in vizvat_regs.
function vizvat_regs_Callback(hObject, eventdata, handles)
% hObject    handle to vizvat_regs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of vizvat_regs
refreshcv(handles);


% --- Executes on button press in vizvat_labs.
function vizvat_labs_Callback(hObject, eventdata, handles)
% hObject    handle to vizvat_labs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of vizvat_labs
refreshcv(handles);


% --- Executes on button press in vizmat_regs.
function vizmat_regs_Callback(hObject, eventdata, handles)
% hObject    handle to vizmat_regs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of vizmat_regs
refreshcv(handles);


% --- Executes on button press in vizmat_labs.
function vizmat_labs_Callback(hObject, eventdata, handles)
% hObject    handle to vizmat_labs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of vizmat_labs
refreshcv(handles);


% --- Executes on button press in wmatlabs.
function wmatlabs_Callback(hObject, eventdata, handles)
% hObject    handle to wmatlabs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of wmatlabs
refreshcv(handles);


% --- Executes on button press in savefibers.
function savefibers_Callback(hObject, eventdata, handles)
% hObject    handle to savefibers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of savefibers
refreshcv(handles);
