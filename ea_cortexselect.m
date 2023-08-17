function varargout = ea_cortexselect(varargin)
% EA_CORTEXSELECT MATLAB code for ea_cortexselect.fig
%      EA_CORTEXSELECT, by itself, creates a new EA_CORTEXSELECT or raises the existing
%      singleton*.
%
%      H = EA_CORTEXSELECT returns the handle to a new EA_CORTEXSELECT or the handle to
%      the existing singleton*.
%
%      EA_CORTEXSELECT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EA_CORTEXSELECT.M with the given input arguments.
%
%      EA_CORTEXSELECT('Property','Value',...) creates a new EA_CORTEXSELECT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ea_cortexselect_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ea_cortexselect_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ea_cortexselect

% Last Modified by GUIDE v2.5 25-Nov-2017 08:08:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @ea_cortexselect_OpeningFcn, ...
    'gui_OutputFcn',  @ea_cortexselect_OutputFcn, ...
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


% --- Executes just before ea_cortexselect is made visible.
function ea_cortexselect_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ea_cortexselect (see VARARGIN)

% Choose default command line output for ea_cortexselect
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ea_cortexselect wait for user response (see UIRESUME)
% uiwait(handles.cortexselect);

% cswin = ea_cortexselect(cortex, annot, atlases, colorindex, struct_names, options, resultfig);

cortex=varargin{1};
annot=varargin{2};
atlases=varargin{3};
colorindex=varargin{4};
struct_names=varargin{5};
options=varargin{6};
resultfig=varargin{7};

if ~isfield(options,'native')
    options.native=0;
end

set(handles.cortexselect,'Visible',options.d3.verbose); % set invisible if called from lead group

movegui(hObject,'southeast');

% available atlases
 set(handles.atlaspopup,'String',atlases);

 % opacity
set(handles.alphaslider,'Value',options.prefs.d3.cortexalpha);
set(handles.alphaedit,'String',num2str(options.prefs.d3.cortexalpha));

setappdata(handles.cortexselect,'handles',handles);
setappdata(handles.cortexselect,'atlases',atlases);
setappdata(handles.cortexselect,'annot',annot);
setappdata(handles.cortexselect,'struct_names',struct_names);
setappdata(handles.cortexselect,'options',options);
setappdata(handles.cortexselect,'resultfig',resultfig);

setappdata(handles.alphaslider,'annot',annot);
setappdata(handles.alphaslider,'struct_names',struct_names);
setappdata(handles.alphaslider,'options',options);
setappdata(handles.alphaslider,'resultfig',resultfig);


axis off
ea_createpcmenu(handles);
setappdata(handles.cortexselect,'treeinit',1);
setuptree([{handles}, varargin])


function setuptree(varargin)

% IO handling
handles=varargin{1}{1};

ea_busyaction('on',handles.cortexselect,'cortex_control');

cortex=varargin{1}{2};
annot=varargin{1}{3};
atlases=varargin{1}{4};
colorindex=varargin{1}{5};
struct_names=varargin{1}{6};
options=varargin{1}{7};
resultfig=varargin{1}{8};

setappdata(handles.cortexselect,'annot',annot);
setappdata(handles.cortexselect,'atlases',atlases);
setappdata(handles.cortexselect,'colorindex',colorindex);
setappdata(handles.cortexselect,'struct_names',struct_names);

cortchecks=cell(length(struct_names),1);
import com.mathworks.mwswing.checkboxtree.*
uselabelname = 1;

% for s = 1:2  % side 1=Right, 2=Left
%     if s==1
%         h.sg{s} = DefaultCheckBoxNode('Right');
%     elseif s==2
%         h.sg{s} = DefaultCheckBoxNode('Left');
%     end
%     for node=1:length(struct_names)
%         thisstruct=struct_names{node};
%
%         color = annot(s).colortable.table(node,1:3);
%         color = sprintf('rgb(%d,%d,%d)', color(1),color(2),color(3));
%
%         structlabel = ['<HTML><BODY>' ...
%             '<FONT color=',color,' bgcolor=',color,'>ico</FONT>' ...
%             '<FONT color="black">&nbsp;&nbsp;',thisstruct,'</FONT>' ...
%             '</BODY></HTML>'];
%         h.sgsub{s}{node}=DefaultCheckBoxNode(structlabel,true);
%         h.sg{s}.add(h.sgsub{s}{node});
%
%          cortchecks{node}=...
%             [cortchecks{node},h.sgsub{s}{node}];
%
%     end
% end
h.sg{1} = DefaultCheckBoxNode('Hemisphere',true);
for s = 1:2  % side 1=Right, 2=Left
    if s==1
        h.sgsub{s} = DefaultCheckBoxNode('Right',true);
        h.sg{1}.add(h.sgsub{s});
    elseif s==2
        h.sgsub{s} = DefaultCheckBoxNode('Left',true);
        h.sg{1}.add(h.sgsub{s});
    end
end

for s=1:2
    for node=1:length(struct_names)
        thisstruct=struct_names{node};

        color = annot(s).colortable.table(node,1:3);
        color = sprintf('rgb(%d,%d,%d)', color(1),color(2),color(3));

        structlabel = ['<HTML><BODY>' ...
            '<FONT color=',color,' bgcolor=',color,'>ico</FONT>' ...
            '<FONT color="black">&nbsp;&nbsp;',thisstruct,'</FONT>' ...
            '</BODY></HTML>'];
        h.sgsubfi{s}{node}=DefaultCheckBoxNode(structlabel,true);
        h.sgsub{s}.add(h.sgsubfi{s}{node});

    end
end

    ea_cleanpriortree(handles);

    % Create a standard MJTree:
    jTree = com.mathworks.mwswing.MJTree(h.sg{1});

    % Now present the CheckBoxTree:
    jCheckBoxTree = CheckBoxTree(jTree.getModel);

    jScrollPane = com.mathworks.mwswing.MJScrollPane(jCheckBoxTree);
    treeinit=getappdata(handles.cortexselect,'treeinit');
    setappdata(handles.cortexselect,'treeinit',0);

    atlN=length(struct_names);
    height=(atlN+1.5)*18;
    norm=360; % max height if full size figure shown.
    if height>360
        height=360;
    end
    if height<100
        height=100;
    end

    [jComp,hc] = ea_javacomponent(jScrollPane,[10,5,285,height],handles.cortexselect);
    setappdata(handles.cortexselect,'uitree',jComp);

    ea_busyaction('del',handles.cortexselect,'atlcontrol');

    h.uselabelname = uselabelname;
    h.cortex = cortex;
    h.annot = annot;
    h.struct_names = struct_names;
    h.colorindex = colorindex;
    h.options = options;
    h.resultfig = resultfig;
    set(jCheckBoxTree, 'MouseReleasedCallback', {@mouseReleasedCallback, h})
    setappdata(handles.cortexselect,'h',h);
    setappdata(handles.cortexselect,'jtree',jCheckBoxTree);
    %sels=ea_storeupdatecortex(jCheckBoxTree,h);

    if treeinit
        if handles.alphaslider.Value>length(handles.alphaslider.String)
            handles.alphaslider.Value=1;
        end
        if handles.togglepopup.Value>length(handles.togglepopup.String)
            handles.togglepopup.Value=1;
        end
        %handles.cortexselect.Position(2)=handles.cortexselect.Position(2)-(450);
        handles.cortexselect.Position(4)=(534-(360-height));

        handles.cortstructxt.Position(2)=handles.cortexselect.Position(4)-25;

        handles.atlaspopup.Position(2)=handles.cortexselect.Position(4)-75;
        handles.atlasstatic.Position(2)=handles.atlaspopup.Position(2)+28;

        handles.togglepopup.Position(2)=handles.cortexselect.Position(4)-120;
        handles.togglestatic.Position(2)=handles.togglepopup.Position(2)+28;

        handles.alphaslider.Position(2)=handles.cortexselect.Position(4)-165;
        handles.alphastatic.Position(2)=handles.alphaslider.Position(2)+28;
        handles.alphaedit.Position(2)=handles.alphastatic.Position(2);

        set(0,'CurrentFigure',handles.cortexselect);
        axis off
        movegui(handles.cortexselect,'southeast');
    end


function ea_cleanpriortree(handles)

jComp=getappdata(handles.cortexselect,'uitree');
if ~isempty(jComp)
    delete(jComp);
end


function [structures,labelidx] = ea_showhideatlases(jtree,h)

updatematrix = ea_updatematrix(jtree,h);
for s = 1:size(updatematrix,2)
    X{s} = mat2cell([1:length(updatematrix{s})]',ones(length(updatematrix{s}),1));
    labelidx{s} = X{s}(logical(updatematrix{s}));
    structures{s} = h.struct_names(logical(updatematrix{s}));
end
ea_updatecortex(h.options,h.resultfig,1:2,structures,labelidx);

selstate.structures=structures;
selstate.labelidx=labelidx;
setappdata(h.resultfig,'cortsurfs',selstate)


function updatematrix = ea_updatematrix(jtree,h)

sels=ea_storeupdatecortex(jtree,h);
for branch=1:length(sels.branches) % Hemisphere
    for leaf=1:length(sels.leaves{branch}) % Label

        % Turn Visibility On
        if strcmp(sels.branches{branch},'mixed') && strcmp(sels.leaves{branch}{leaf},'selected')
            updatematrix{branch}(leaf) = 1;

        elseif strcmp(sels.branches{branch},'selected') && strcmp(sels.leaves{branch}{leaf},'selected')
            updatematrix{branch}(leaf) = 1;

        % Turn Visibility Off
        elseif strcmp(sels.branches{branch},'mixed') && strcmp(sels.leaves{branch}{leaf},'not selected')
            updatematrix{branch}(leaf) = 0;

        elseif strcmp(sels.branches{branch},'not selected')
            updatematrix{branch}(leaf) = 0;
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
newselstate=ea_storeupdatecortex(jtree,h);
if ~isequal(oldselstate,newselstate)
    ea_showhideatlases(jtree,h);
end


% --- Outputs from this function are returned to the command line.
function varargout = ea_cortexselect_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in togglepopup.
function togglepopup_Callback(hObject, eventdata, handles)
% hObject    handle to togglepopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns togglepopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from togglepopup

%pcmenu=getappdata(handles.togglepopup,'uimenu');

%showcontextmenu(hObject,pcmenu);

presetactions=getappdata(handles.togglepopup,'presetactions');
ea_makeselection([],[],handles,presetactions{handles.togglepopup.Value});


function ea_createpcmenu(handles)
annot=getappdata(handles.cortexselect,'annot');
struct_names=getappdata(handles.cortexselect,'struct_names');

def.presets(1).label='Show all structures';
def.presets(1).show=1:length(struct_names);
def.presets(1).hide=[];
def.presets(1).default='absolute';
def.presets(2).label='Hide all structures';
def.presets(2).show=[];
def.presets(2).hide=1:length(struct_names);
def.presets(2).default='absolute';
prescell=cell(0);
presetactions=cell(0);
% add defaults:
for ps=1:length(def.presets)
    prescell{end+1}=def.presets(ps).label;
    presetactions{end+1}=def.presets(ps);
    %    uimenu(pcmenu, 'Label',,'Callback',{@ea_makeselection,handles,def.togglepopup(ps)});
end

% add from atlas index:
if isfield(annot,'presets')
    for ps=1:length(annot.presets)
        prescell{end+1}=annot.presets(ps).label;
        presetactions{end+1}=annot.presets(ps);

        %        uimenu(pcmenu, 'Label',atlases.togglepopup(ps).label,'Callback',{@ea_makeselection,handles,atlases.togglepopup(ps)});
    end
end

handles.togglepopup.String=prescell;
handles.togglepopup.Value=1;

if isfield(annot,'defaultset')
    if length(handles.togglepopup.String)>2 % custom sets available
        handles.togglepopup.Value=annot.defaultset+2;
    end
end
setappdata(handles.togglepopup,'presetactions',presetactions);


function ea_saveselection(~,~,handles,options)
ea_busyaction('on',handles.cortexselect,'atlcontrol');

jtree=getappdata(handles.cortexselect,'jtree');
h=getappdata(handles.cortexselect,'h');
sels=ea_storeupdatecortex(jtree,h);
atlases=getappdata(handles.cortexselect,'atlases');
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

try WinOnTop(handles.cortexselect,false); end

tag=inputdlg('Please enter a name for the preset:','Preset name');
pres.label=tag{1};
try WinOnTop(handles.cortexselect,true); end

% refresh content menu.
ea_createpcmenu(handles)

ea_busyaction('off',handles.cortexselect,'atlcontrol');


function str=getridofspaces(str)
str=strrep(str,'(','');
str=strrep(str,')','');
str=strrep(str,' ','');
str=strrep(str,'-','');


function ea_makeselection(~,~,handles,preset)

ea_busyaction('on',handles.cortexselect,'atlcontrol');

h=getappdata(handles.cortexselect,'h');
jtree=getappdata(handles.cortexselect,'jtree');
atlases=getappdata(handles.cortexselect,'atlases');
onatlasnames=atlases.names(preset.show);
offatlasnames=atlases.names(preset.hide);
% get rid of file extensions:
[~,onatlasnames]=cellfun(@fileparts,onatlasnames,'Uniformoutput',0);
[~,onatlasnames]=cellfun(@fileparts,onatlasnames,'Uniformoutput',0);
[~,offatlasnames]=cellfun(@fileparts,offatlasnames,'Uniformoutput',0);
[~,offatlasnames]=cellfun(@fileparts,offatlasnames,'Uniformoutput',0);

% iterate through jTree to set selection according to preset:
sels=ea_storeupdatecortex(jtree,h);
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
ea_busyaction('off',handles.cortexselect,'atlcontrol');

ea_synctree(handles)


% --- Executes during object creation, after setting all properties.
function togglepopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to togglepopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in atlaspopup.
function atlaspopup_Callback(hObject, eventdata, handles)
% hObject    handle to atlaspopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns atlaspopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from atlaspopup


ea_busyaction('on',handles.cortexselect,'atlcontrol');

% retrieve necessary info from cortexselect figure:
resultfig=getappdata(handles.cortexselect,'resultfig');
options=getappdata(handles.cortexselect,'options');

% surfaces
atlassurfs=getappdata(resultfig,'atlassurfs');
for atl=1:numel(atlassurfs)
    delete(atlassurfs{atl})
end

% labels
atlaslabels=getappdata(resultfig,'atlaslabels');
for atl=1:numel(atlaslabels)
    delete(atlaslabels(atl))
end

elstruct=getappdata(resultfig,'elstruct');
options.atlasset=get(handles.atlaspopup,'String'); %{get(handles.atlaspopup,'Value')}
options.atlasset=options.atlasset{get(handles.atlaspopup,'Value')};
options.atlassetn=get(handles.atlaspopup,'Value');
setappdata(resultfig,'options',options); % update options in resultfig for VAT model
[atlases,colorbuttons,atlassurfs,atlaslabels]=ea_showatlas(resultfig,elstruct,options);
setappdata(handles.cortexselect,'atlases',atlases);
setappdata(handles.cortexselect,'treeinit',1);
labelbutton = getappdata(resultfig,'labelbutton');

setuptree({handles,colorbuttons,atlassurfs,atlases,labelbutton,atlaslabels});
ea_createpcmenu(handles);
ea_busyaction('off',handles.cortexselect,'atlcontrol');


% --- Executes during object creation, after setting all properties.
function atlaspopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to atlaspopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in alphaslider.
function alphaslider_Callback(hObject, eventdata, handles)
% hObject    handle to alphaslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns alphaslider contents as cell array
%        contents{get(hObject,'Value')} returns selected item from alphaslider

% ea_updatecortex(options,resultfig,sides,structures,labelidx,alpha)
alpha = handles.alphaslider.Value;
resultfig = getappdata(hObject,'resultfig');
options = getappdata(hObject,'options');
annot = getappdata(hObject,'annot');
selstate = getappdata(resultfig,'cortsurfs');

ea_updatecortex(options,resultfig,1:2,selstate.structures,selstate.labelidx,alpha);
set(handles.alphaedit,'String',num2str(alpha));
% setuptree({handles,colorbuttons,atlassurfs,atlases,labelbutton,atlaslabels});


% --- Executes during object creation, after setting all properties.
function alphaslider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alphaslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function alphaedit_Callback(hObject, eventdata, handles)
% hObject    handle to alphaedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of alphaedit as text
%        str2double(get(hObject,'String')) returns contents of alphaedit as a double


% --- Executes during object creation, after setting all properties.
function alphaedit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alphaedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
