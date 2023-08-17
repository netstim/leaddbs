function varargout = ea_atlasselect(varargin)
% EA_ATLASSELECT MATLAB code for ea_atlasselect.fig
%      EA_ATLASSELECT, by itself, creates a new EA_ATLASSELECT or raises the existing
%      singleton*.
%
%      H = EA_ATLASSELECT returns the handle to a new EA_ATLASSELECT or the handle to
%      the existing singleton*.
%
%      EA_ATLASSELECT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EA_ATLASSELECT.M with the given input arguments.
%
%      EA_ATLASSELECT('Property','Value',...) creates a new EA_ATLASSELECT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ea_atlasselect_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ea_atlasselect_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ea_atlasselect

% Last Modified by GUIDE v2.5 04-May-2017 19:22:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @ea_atlasselect_OpeningFcn, ...
    'gui_OutputFcn',  @ea_atlasselect_OutputFcn, ...
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


% --- Executes just before ea_atlasselect is made visible.
function ea_atlasselect_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ea_atlasselect (see VARARGIN)

% Choose default command line output for ea_atlasselect
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ea_atlasselect wait for user response (see UIRESUME)
% uiwait(handles.atlasselect);

% ea_atlasselect(colorbuttons, atlassurfs, atlases, options, resultfig, labelbutton, atlaslabels)

atlases=varargin{3};
options=varargin{4};
resultfig=varargin{5};

if ~isfield(options,'native')
    options.native=0;
end

movegui(hObject,'northeast'); drawnow;

set(handles.atlasselect,'Visible',options.d3.verbose); % set invisible if called from lead group

ea_listatlassets(options,handles,options.native);

[~,handles.atlassetpopup.Value]=ismember(options.atlasset,handles.atlassetpopup.String);
if handles.atlassetpopup.Value==0
handles.atlassetpopup.String=[handles.atlassetpopup.String;{options.atlasset}];
[~,handles.atlassetpopup.Value]=ismember(options.atlasset,handles.atlassetpopup.String);
end
setappdata(handles.atlasselect,'handles',handles);
setappdata(handles.atlasselect,'atlases',atlases);
setappdata(handles.atlasselect,'options',options);
setappdata(handles.atlasselect,'resultfig',resultfig);

axis off
ea_createpcmenu(handles);
setappdata(handles.atlasselect,'treeinit',1);
setuptree([{handles}, varargin])

ea_ListBoxRenderer(handles.atlassetpopup);
ea_ListBoxRenderer(handles.presets);


function setuptree(varargin)

% IO handling
handles=varargin{1}{1};

ea_busyaction('on',handles.atlasselect,'atlcontrol');

togglebuttons=varargin{1}{2};
atlassurfs=varargin{1}{3};
atlases=varargin{1}{4};
labelbutton=varargin{1}{end-1};
atlaslabels=varargin{1}{end};

setappdata(handles.atlasselect,'colorbuttons',togglebuttons);
setappdata(handles.atlasselect,'atlassurfs',atlassurfs);
setappdata(handles.atlasselect,'atlases',atlases);
setappdata(handles.atlasselect,'labelbutton',labelbutton);
setappdata(handles.atlasselect,'atlaslabels',atlaslabels);

% check if other labelsets are available
if isfield(atlases,'labelnames')
    handles.namingscheme.String=[atlases.labelnames,{'NIfTI filenames'}];
else
    handles.namingscheme.String={'NIfTI filenames'};
end

try
    if ~isfield(atlases,'subgroups')
        atlases.subgroups(1).label='Structures';
        isdiscfibers = cellfun(@(x) ischar(x) && strcmp(x, 'discfibers'), atlases.pixdim);
        atlases.subgroups(1).entries=find(~sum(isdiscfibers,2)');
    end
catch
    keyboard
end

for tb=1:length(togglebuttons)
    tbcell{tb}=togglebuttons(tb).Tag;
end

atlchecks=cell(length(atlases.names),1);
import com.mathworks.mwswing.checkboxtree.*
if handles.namingscheme.Value>length(handles.namingscheme.String)
    handles.namingscheme.Value=1;
end

if handles.presets.Value>length(handles.presets.String)
    handles.presets.Value=1;
end

switch handles.namingscheme.String{handles.namingscheme.Value}
    case 'NIfTI filenames'  % use atlases.names
        uselabelname = 0;
    otherwise  % use atlases.labels
        [~,uselabelname] = ismember(handles.namingscheme.String{handles.namingscheme.Value},atlases.labelnames);
        if uselabelname == 0
            uselabelname = 1;
        end
end

for subgroup=1:length(atlases.subgroups)
    h.sg{subgroup} = DefaultCheckBoxNode(atlases.subgroups(subgroup).label);
    for node=1:length(atlases.subgroups(subgroup).entries)
        [~,thisatlfname]=fileparts(atlases.names{atlases.subgroups(subgroup).entries(node)});
        [~,thisatlfname]=fileparts(thisatlfname);

        if uselabelname == 0
            thisatlname = thisatlfname;
        else
            thisatlname=atlases.labels{uselabelname}{node};
        end

        % nii(.gz) support
        if regexp(thisatlname, '\.nii(\.gz)?$')
            [~,thisatlname]=ea_niifileparts(thisatlname);
        end

        color = round(squeeze(atlases.colormap(ceil(atlases.colors(atlases.subgroups(subgroup).entries(node))),:))*256);
        color = sprintf('rgb(%d,%d,%d)', color(1),color(2),color(3));

        atlaslabel = ['<HTML><BODY>' ...
            '<FONT color=',color,' bgcolor=',color,'>ico</FONT>' ...
            '<FONT color="black">&nbsp;&nbsp;',thisatlname,'</FONT>' ...
            '</BODY></HTML>'];
        h.sgsub{subgroup}{node}=DefaultCheckBoxNode(atlaslabel,true);
        h.sg{subgroup}.add(h.sgsub{subgroup}{node});
        switch atlases.types(atlases.subgroups(subgroup).entries(node))
            case {3,4,6,10} % need lh and rh entries
                [~,thistb]=ismember([thisatlfname,'_right'],tbcell);
                try
                    checked=onoff2bool(togglebuttons(thistb).State);
                catch
                    keyboard
                end
                lrlabel = ['<HTML><BODY>' ...
                    '<FONT color=',color,' bgcolor=',color,'>ico</FONT>' ...
                    '<FONT color="black">&nbsp;&nbsp;RH</FONT>' ...
                    '</BODY></HTML>'];
                h.sgsubside{subgroup}{node}{1}=DefaultCheckBoxNode(lrlabel,checked);
                h.sgsub{subgroup}{node}.add(h.sgsubside{subgroup}{node}{1});
                h.sgsubfi{subgroup}{node}=thisatlfname;

                [~,thistb]=ismember([thisatlfname,'_left'],tbcell);
                checked=onoff2bool(togglebuttons(thistb).State);
                lrlabel = ['<HTML><BODY>' ...
                    '<FONT color=',color,' bgcolor=',color,'>ico</FONT>' ...
                    '<FONT color="black">&nbsp;&nbsp;LH</FONT>' ...
                    '</BODY></HTML>'];
                h.sgsubside{subgroup}{node}{2}=DefaultCheckBoxNode(lrlabel,checked);
                h.sgsub{subgroup}{node}.add(h.sgsubside{subgroup}{node}{2});
            case 1 % RH only
                [~,thistb]=ismember([thisatlfname,'_right'],tbcell);
                checked=onoff2bool(togglebuttons(thistb).State);
                lrlabel = ['<HTML><BODY>' ...
                    '<FONT color=',color,' bgcolor=',color,'>ico</FONT>' ...
                    '<FONT color="black">&nbsp;&nbsp;RH</FONT>' ...
                    '</BODY></HTML>'];
                h.sgsubside{subgroup}{node}{1}=DefaultCheckBoxNode(lrlabel,checked);
                h.sgsub{subgroup}{node}.add(h.sgsubside{subgroup}{node}{1});
                h.sgsubfi{subgroup}{node}=thisatlfname;
            case 2 % LH only
                [~,thistb]=ismember([thisatlfname,'_left'],tbcell);
                checked=onoff2bool(togglebuttons(thistb).State);
                lrlabel = ['<HTML><BODY>' ...
                    '<FONT color=',color,' bgcolor=',color,'>ico</FONT>' ...
                    '<FONT color="black">&nbsp;&nbsp;LH</FONT>' ...
                    '</BODY></HTML>'];
                h.sgsubside{subgroup}{node}{1}=DefaultCheckBoxNode(lrlabel,checked);
                h.sgsub{subgroup}{node}.add(h.sgsubside{subgroup}{node}{1});
                h.sgsubfi{subgroup}{node}=thisatlfname;
            case 5 % Midline
                [~,thistb]=ismember([thisatlfname,'_midline'],tbcell);
                checked=onoff2bool(togglebuttons(thistb).State);
                lrlabel = ['<HTML><BODY>' ...
                    '<FONT color=',color,' bgcolor=',color,'>ico</FONT>' ...
                    '<FONT color="black">&nbsp;&nbsp;Midline</FONT>' ...
                    '</BODY></HTML>'];
                h.sgsubside{subgroup}{node}{1}=DefaultCheckBoxNode(lrlabel,checked);
                h.sgsub{subgroup}{node}.add(h.sgsubside{subgroup}{node}{1});
                h.sgsubfi{subgroup}{node}=thisatlfname;
        end
        atlchecks{atlases.subgroups(subgroup).entries(node)}=...
            [atlchecks{atlases.subgroups(subgroup).entries(node)},h.sgsub{subgroup}{node}];

    end
end

ea_cleanpriortree(handles);

% Create a standard MJTree:
jTree = com.mathworks.mwswing.MJTree(h.sg{1});

% Now present the CheckBoxTree:
jCheckBoxTree = CheckBoxTree(jTree.getModel);

jScrollPane = com.mathworks.mwswing.MJScrollPane(jCheckBoxTree);
treeinit=getappdata(handles.atlasselect,'treeinit');
setappdata(handles.atlasselect,'treeinit',0);

atlN=length(atlases.names);
height=(atlN+1.5)*18;
norm=360; % max height if full size figure shown.
if height>360
    height=360;
end
if height<100
    height=100;
end

[jComp,hc] = ea_javacomponent(jScrollPane,[10,5,285,height],handles.atlasselect);
setappdata(handles.atlasselect,'uitree',jComp);

ea_busyaction('del',handles.atlasselect,'atlcontrol');

h.togglebuttons=togglebuttons;
h.atlassurfs=atlassurfs;
h.atlases=atlases;
h.atlchecks=atlchecks;
h.uselabelname = uselabelname;
h.labelbutton = labelbutton;
h.atlaslabels = reshape(atlaslabels',1,[]);
set(jCheckBoxTree, 'MouseReleasedCallback', {@mouseReleasedCallback, h})
setappdata(handles.atlasselect,'h',h);
setappdata(handles.atlasselect,'jtree',jCheckBoxTree);
sels=ea_storeupdatemodel(jCheckBoxTree,h);

if treeinit
    if handles.namingscheme.Value>length(handles.namingscheme.String)
        handles.namingscheme.Value=1;
    end
    if handles.presets.Value>length(handles.presets.String)
        handles.presets.Value=1;
    end
    %handles.atlasselect.Position(2)=handles.atlasselect.Position(2)-(450);
    handles.atlasselect.Position(4)=(534-(360-height));

    handles.atlstructxt.Position(2)=handles.atlasselect.Position(4)-25;

    handles.atlassetpopup.Position(2)=handles.atlasselect.Position(4)-75;
    handles.atlassetstatic.Position(2)=handles.atlassetpopup.Position(2)+28;

    handles.presets.Position(2)=handles.atlasselect.Position(4)-120;
    handles.presetsstatic.Position(2)=handles.presets.Position(2)+28;

    handles.namingscheme.Position(2)=handles.atlasselect.Position(4)-165;
    handles.namingstatic.Position(2)=handles.namingscheme.Position(2)+28;

    set(0,'CurrentFigure',handles.atlasselect);
    axis off
    movegui(handles.atlasselect,'northeast');
end


function ea_cleanpriortree(handles)

jComp=getappdata(handles.atlasselect,'uitree');
if ~isempty(jComp)
    delete(jComp);
end


function tf=bin2bool(t)
if t
    tf=true;
else
    tf=false;
end

function tf=onoff2bool(t)
switch t
    case 'on'
        tf=true;
    case 'off'
        tf=false;
end


function ea_showhideatlases(jtree,h)

sels=ea_storeupdatemodel(jtree,h);
for branch=1:length(sels.branches)
    for leaf=1:length(sels.leaves{branch})
        if ~isempty(sels.sides{branch}{leaf}) % has side children
            for side=1:length(sels.sides{branch}{leaf})

                sidec=getsidec(length(sels.sides{branch}{leaf}),side,h.atlases.types(leaf));

                [ixs,ixt]=ea_getsubindex(h.sgsub{branch}{leaf}.toString,sidec,h.atlassurfs,h.togglebuttons,h.uselabelname,h.atlases);
                if strcmp(sels.sides{branch}{leaf}{side},'selected')
                    if ~strcmp(h.atlassurfs{ixs}.Visible,'on')
                        h.atlassurfs{ixs}.Visible='on';
                        h.togglebuttons(ixt).State='on';
                    end

                    if strcmp(h.labelbutton.State, 'on')
                        h.atlaslabels(ixs).Visible='on';
                    end
                elseif strcmp(sels.sides{branch}{leaf}{side},'not selected')
                    if ~strcmp(h.atlassurfs{ixs}.Visible,'off')
                        h.atlassurfs{ixs}.Visible='off';
                        h.togglebuttons(ixt).State='off';
                    end

                    if strcmp(h.labelbutton.State, 'on')
                        h.atlaslabels(ixs).Visible='off';
                    end
                end
            end

        else
            keyboard
        end
    end

end


function sidec=getsidec(sel,side,type)

if sel==2
    switch side
        case 1
            sidec='_right';
        case 2
            sidec='_left';
    end
elseif sel==1
    switch type
        case 1
            sidec='_right';
        case 2
            sidec='_left';
        case 5
            sidec='_midline';
    end
end


% Set the mouse-press callback
function mouseReleasedCallback(jtree, eventData, h)

clickX = eventData.getX;
clickY = eventData.getY;
treePath = jtree.getPathForLocation(clickX, clickY);

oldselstate=getappdata(jtree,'selectionstate');
newselstate=ea_storeupdatemodel(jtree,h);
if ~isequal(oldselstate,newselstate)
    ea_showhideatlases(jtree,h);
end


% --- Outputs from this function are returned to the command line.
function varargout = ea_atlasselect_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in presets.
function presets_Callback(hObject, eventdata, handles)
% hObject    handle to presets (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns presets contents as cell array
%        contents{get(hObject,'Value')} returns selected item from presets

%pcmenu=getappdata(handles.presets,'uimenu');

%showcontextmenu(hObject,pcmenu);

presetactions=getappdata(handles.presets,'presetactions');
ea_makeselection([],[],handles,presetactions{handles.presets.Value});


function ea_createpcmenu(handles)
atlases=getappdata(handles.atlasselect,'atlases');

def.presets(1).label='Show all structures';
def.presets(1).show=1:length(atlases.names);
def.presets(1).hide=[];
def.presets(1).default='absolute';
def.presets(2).label='Hide all structures';
def.presets(2).show=[];
def.presets(2).hide=1:length(atlases.names);
def.presets(2).default='absolute';
prescell=cell(0);
presetactions=cell(0);
% add defaults:
for ps=1:length(def.presets)
    prescell{end+1}=def.presets(ps).label;
    presetactions{end+1}=def.presets(ps);
    %    uimenu(pcmenu, 'Label',,'Callback',{@ea_makeselection,handles,def.presets(ps)});
end

% add from atlas index:
if isfield(atlases,'presets')
    for ps=1:length(atlases.presets)
        prescell{end+1}=atlases.presets(ps).label;
        presetactions{end+1}=atlases.presets(ps);

        %        uimenu(pcmenu, 'Label',atlases.presets(ps).label,'Callback',{@ea_makeselection,handles,atlases.presets(ps)});
    end
end

handles.presets.String=prescell;
handles.presets.Value=1;

if isfield(atlases,'defaultset')
    if length(handles.presets.String)>2 % custom sets available
        handles.presets.Value=atlases.defaultset+2;
    end
end
setappdata(handles.presets,'presetactions',presetactions);


function ea_saveselection(~,~,handles,options)
ea_busyaction('on',handles.atlasselect,'atlcontrol');

jtree=getappdata(handles.atlasselect,'jtree');
h=getappdata(handles.atlasselect,'h');
sels=ea_storeupdatemodel(jtree,h);
atlases=getappdata(handles.atlasselect,'atlases');
pres.default='absolute';
pres.show=[];
pres.hide=[];
[~,atlases.names]=cellfun(@fileparts,atlases.names,'Uniformoutput',0);
[~,atlases.names]=cellfun(@fileparts,atlases.names,'Uniformoutput',0);

for branch=1:length(sels.branches)
    for leaf=1:length(sels.leaves{branch})
        for side=1:length(sels.sides{branch}{leaf})

            sidec=getsidec(length(sels.sides{branch}{leaf}),side,atlases.types(leaf));

            %[ixs,ixt]=getsubindex(h.sgsub{branch}{leaf},sidec,h.atlassurfs,h.togglebuttons);

            [~,ix]=ismember(char(h.sgsub{branch}{leaf}),atlases.names);
            if strcmp(sels.sides{branch}{leaf}{side},'selected')
                pres.show=[pres.show,ix];


            elseif strcmp(sels.sides{branch}{leaf}{side},'not selected')
                pres.hide=[pres.hide,ix];
            end
        end

    end
end

try WinOnTop(handles.atlasselect,false); end

tag=inputdlg('Please enter a name for the preset:','Preset name');
pres.label=tag{1};
try WinOnTop(handles.atlasselect,true); end

% refresh content menu.
ea_createpcmenu(handles)

ea_busyaction('off',handles.atlasselect,'atlcontrol');


function str=getridofspaces(str)
str=strrep(str,'(','');
str=strrep(str,')','');
str=strrep(str,' ','');
str=strrep(str,'-','');


function ea_makeselection(~,~,handles,preset)

ea_busyaction('on',handles.atlasselect,'atlcontrol');

h=getappdata(handles.atlasselect,'h');
jtree=getappdata(handles.atlasselect,'jtree');
atlases=getappdata(handles.atlasselect,'atlases');
onatlasnames=atlases.names(preset.show);
offatlasnames=atlases.names(preset.hide);
% get rid of file extensions:
[~,onatlasnames]=cellfun(@fileparts,onatlasnames,'Uniformoutput',0);
[~,onatlasnames]=cellfun(@fileparts,onatlasnames,'Uniformoutput',0);
[~,offatlasnames]=cellfun(@fileparts,offatlasnames,'Uniformoutput',0);
[~,offatlasnames]=cellfun(@fileparts,offatlasnames,'Uniformoutput',0);

% iterate through jTree to set selection according to preset:
sels=ea_storeupdatemodel(jtree,h);
for branch=1:length(sels.branches)
    for leaf=1:length(sels.leaves{branch})
        for side=1:length(sels.sides{branch}{leaf})

            sidec=getsidec(length(sels.sides{branch}{leaf}),side,atlases.types(leaf));
            [ixs,ixt]=ea_getsubindex(h.sgsub{branch}{leaf}.toString,sidec,h.atlassurfs,h.togglebuttons,h.uselabelname,h.atlases);

            if ismember(char(h.sgsubfi{branch}{leaf}),onatlasnames)
                h.atlassurfs{ixs}.Visible='on';
                if strcmp(h.labelbutton.State, 'on')
                    h.atlaslabels(ixs).Visible='on';
                end
                h.togglebuttons(ixt).State='on';
            elseif ismember(char(h.sgsubfi{branch}{leaf}),offatlasnames)
                h.atlassurfs{ixs}.Visible='off';
                h.atlaslabels(ixs).Visible='off';
                h.togglebuttons(ixt).State='off';
            else % not explicitly mentioned
                switch preset.default
                    case 'absolute'
                        h.atlassurfs{ixs}.Visible='off';
                        h.atlaslabels(ixs).Visible='off';
                        h.togglebuttons(ixt).State='off';
                    case 'relative'
                        % leave state as is.
                end
            end
        end

    end
end
ea_busyaction('off',handles.atlasselect,'atlcontrol');

ea_synctree(handles)


% --- Executes during object creation, after setting all properties.
function presets_CreateFcn(hObject, eventdata, handles)
% hObject    handle to presets (see GCBO)
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


ea_busyaction('on',handles.atlasselect,'atlcontrol');

% retrieve necessary info from atlasselect figure:
resultfig=getappdata(handles.atlasselect,'resultfig');
options=getappdata(handles.atlasselect,'options');

% Clear atlas toolbar
atlastoolbar = getappdata(resultfig,'addht');
arrayfun(@delete, atlastoolbar.Children);

% surfaces
atlassurfs=getappdata(resultfig,'atlassurfs');
for atl=1:numel(atlassurfs)
    delete(atlassurfs{atl});
end

% labels
atlaslabels=getappdata(resultfig,'atlaslabels');
for atl=1:numel(atlaslabels)
    delete(atlaslabels(atl))
end

elstruct=getappdata(resultfig,'elstruct');
options.atlasset=get(handles.atlassetpopup,'String'); %{get(handles.atlassetpopup,'Value')}
options.atlasset=options.atlasset{get(handles.atlassetpopup,'Value')};
options.atlassetn=get(handles.atlassetpopup,'Value');
setappdata(resultfig,'options',options); % update options in resultfig for VAT model
[atlases,colorbuttons,atlassurfs,atlaslabels]=ea_showatlas(resultfig,elstruct,options);
setappdata(handles.atlasselect,'atlases',atlases);
setappdata(handles.atlasselect,'treeinit',1);
labelbutton = getappdata(resultfig,'labelbutton');

setuptree({handles,colorbuttons,atlassurfs,atlases,labelbutton,atlaslabels});
ea_createpcmenu(handles);
ea_busyaction('off',handles.atlasselect,'atlcontrol');


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


% --- Executes on selection change in namingscheme.
function namingscheme_Callback(hObject, eventdata, handles)
% hObject    handle to namingscheme (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns namingscheme contents as cell array
%        contents{get(hObject,'Value')} returns selected item from namingscheme


colorbuttons=getappdata(handles.atlasselect,'colorbuttons');
atlassurfs=getappdata(handles.atlasselect,'atlassurfs');
atlases=getappdata(handles.atlasselect,'atlases');
labelbutton=getappdata(handles.atlasselect,'labelbutton');
atlaslabels=getappdata(handles.atlasselect,'atlaslabels');

setuptree({handles,colorbuttons,atlassurfs,atlases,labelbutton,atlaslabels});


% --- Executes during object creation, after setting all properties.
function namingscheme_CreateFcn(hObject, eventdata, handles)
% hObject    handle to namingscheme (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
