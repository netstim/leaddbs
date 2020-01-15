function varargout = ea_lg_3dsetting(varargin)
% EA_LG_3DSETTING MATLAB code for ea_lg_3dsetting.fig
%      EA_LG_3DSETTING, by itself, creates a new EA_LG_3DSETTING or raises the existing
%      singleton*.
%
%      H = EA_LG_3DSETTING returns the handle to a new EA_LG_3DSETTING or the handle to
%      the existing singleton*.
%
%      EA_LG_3DSETTING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EA_LG_3DSETTING.M with the given input arguments.
%
%      EA_LG_3DSETTING('Property','Value',...) creates a new EA_LG_3DSETTING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ea_lg_3dsetting_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ea_lg_3dsetting_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ea_lg_3dsetting

% Last Modified by GUIDE v2.5 14-Jan-2020 12:34:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ea_lg_3dsetting_OpeningFcn, ...
                   'gui_OutputFcn',  @ea_lg_3dsetting_OutputFcn, ...
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


% --- Executes just before ea_lg_3dsetting is made visible.
function ea_lg_3dsetting_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ea_lg_3dsetting (see VARARGIN)

% Choose default command line output for ea_lg_3dsetting
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ea_lg_3dsetting wait for user response (see UIRESUME)
% uiwait(handles.lg_3dsetting);

leadfigure = varargin{1};

setappdata(handles.lg_3dsetting, 'leadfigure', leadfigure);

M = getappdata(leadfigure, 'M');
set(handles.elmodelselect,'String',[{'Patient specified'},ea_resolve_elspec]);
try
	set(handles.elmodelselect, 'Value', M.ui.elmodelselect);
catch
    set(handles.elmodelselect, 'Value', 1);
end

try
	set(handles.elrenderingpopup, 'Value', M.ui.elrendering);
catch
    set(handles.elrenderingpopup, 'Value', 1);
end

try
	set(handles.isovscloudpopup, 'Value', M.ui.isovscloudpopup);
catch
    set(handles.isovscloudpopup, 'Value', 1);
end

try
	set(handles.colorpointcloudcheck, 'Value', M.ui.colorpointcloudcheck);
catch
    set(handles.colorpointcloudcheck, 'Value', 1);
end

if get(handles.elrenderingpopup, 'Value') == 3
    set(handles.colorpointcloudcheck, 'Enable', 'on');
else
    set(handles.colorpointcloudcheck, 'Enable', 'off');
end

try
	set(handles.mirrorsides, 'Value', M.ui.mirrorsides);
catch
    set(handles.mirrorsides, 'Value', 0);
end

try
	set(handles.showdiscfibers, 'Value', M.ui.showdiscfibers);
catch
    set(handles.showdiscfibers, 'Value', 0);
end

prefs = ea_prefs('');
discfibers = prefs.machine.lg.discfibers;
switch discfibers.showfibersset
    case 'positive'
        set(handles.showposonly, 'Value', 1);
        set(handles.pospredthreshold, 'Enable', 'on');
        set(handles.negpredthreshold, 'Enable', 'off');
    case 'negative'
        set(handles.shownegonly, 'Value', 1);
        set(handles.pospredthreshold, 'Enable', 'off');
        set(handles.negpredthreshold, 'Enable', 'on');
    case 'both'
        set(handles.showboth, 'Value', 1);
        set(handles.pospredthreshold, 'Enable', 'on');
        set(handles.negpredthreshold, 'Enable', 'on');
end
set(handles.pospredthreshold, 'String', num2str(discfibers.pospredthreshold));
set(handles.negpredthreshold, 'String', num2str(discfibers.negpredthreshold));
set(handles.statmetric,'Value',discfibers.statmetric);


% --- Outputs from this function are returned to the command line.
function varargout = ea_lg_3dsetting_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in elmodelselect.
function elmodelselect_Callback(hObject, eventdata, handles)
% hObject    handle to elmodelselect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns elmodelselect contents as cell array
%        contents{get(hObject,'Value')} returns selected item from elmodelselect
leadfigure = getappdata(handles.lg_3dsetting, 'leadfigure');

M = getappdata(leadfigure, 'M');

if eventdata.Source.Value ~= M.ui.elmodelselect
    if ~isempty(M.patient.list)
        disp('Updating electrode model...');
        ea_busyaction('on', leadfigure, 'group');
        M.ui.elmodelselect = eventdata.Source.Value;
        for pt=1:length(M.patient.list)
            % load localization
            [~, patientname] = fileparts(M.patient.list{pt});

            M.elstruct(pt).group = M.patient.group(pt);
            M.elstruct(pt).groupcolors = M.groups.color;
            M.elstruct(pt).groups = M.groups.group;

            options.sides = 1:2;
            options.native = 0;
            try
                [options.root,options.patientname] = fileparts(M.patient.list{pt});
                options.root = [options.root, filesep];
                options = ea_resolve_elspec(options);
                if exist([options.root,options.patientname,filesep,'ea_reconstruction.mat'],'file')
                    [~,~,markers,elmodel,~,coords_acpc] = ea_load_reconstruction(options);

                    if M.ui.elmodelselect == 1 % use patient specific elmodel
                        if exist('elmodel','var')
                            M.elstruct(pt).elmodel = elmodel;
                        else % use default for older reconstructions that did not store elmodel.
                            M.elstruct(pt).elmodel = 'Medtronic 3389';
                        end
                    else
                        elmodels = [{'Patient specified'},ea_resolve_elspec];
                        M.elstruct(pt).elmodel = elmodels{M.ui.elmodelselect};
                    end

                    % make sure coords_mm is congruent to coded electrode model
                    poptions = options;
                    poptions.native = 0;
                    poptions.elmodel = M.elstruct(pt).elmodel;
                    poptions = ea_resolve_elspec(poptions);
                    [coords_mm,trajectory,markers] = ea_resolvecoords(markers,poptions,0);

                    M.elstruct(pt).coords_mm = coords_mm;
                    M.elstruct(pt).coords_acpc = coords_acpc;
                    M.elstruct(pt).trajectory = trajectory;

                    M.elstruct(pt).name = patientname;
                    if ~exist('markers','var') % backward compatibility to old recon format
                        for side=1:2
                            markers(side).head = coords_mm{side}(1,:);
                            markers(side).tail = coords_mm{side}(4,:);
                            normtrajvector = (markers(side).tail-markers(side).head)./norm(markers(side).tail-markers(side).head);
                            orth = null(normtrajvector)*(options.elspec.lead_diameter/2);
                            markers(side).x = coords_mm{side}(1,:)+orth(:,1)';
                            markers(side).y = coords_mm{side}(1,:)+orth(:,2)'; % corresponding points in reality
                        end
                    end
                    M.elstruct(pt).markers = markers;
                end
            catch
                if pt>1 % first patient has worked but some other patient seems not to have worked.
                    try
                        if ~M.ui.detached
                            M.elstruct(1).coords_mm; % probe if error happens in pt. 1 ? if not show warning
                            warning(['No reconstruction present for ',patientname,'. Please check.']);
                        end
                    end
                end
            end
        end

        disp('Storing everything in model...');
        timestamp = datetime('now');
        timestamp.Format = 'uuuMMddHHmmss';
        M.ui.lastupdated = str2double(char(timestamp));
        setappdata(leadfigure, 'M', M);
        
        ea_busyaction('off', leadfigure, 'group');
    end
end


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


% --- Executes on selection change in elrenderingpopup.
function elrenderingpopup_Callback(hObject, eventdata, handles)
% hObject    handle to elrenderingpopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns elrenderingpopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from elrenderingpopup
if get(handles.elrenderingpopup, 'Value') == 3
    set(handles.colorpointcloudcheck, 'Enable', 'on');
else
    set(handles.colorpointcloudcheck, 'Enable', 'off');
end


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


% --- Executes on selection change in isovscloudpopup.
function isovscloudpopup_Callback(hObject, eventdata, handles)
% hObject    handle to isovscloudpopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns isovscloudpopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from isovscloudpopup


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


% --- Executes on button press in colorpointcloudcheck.
function colorpointcloudcheck_Callback(hObject, eventdata, handles)
% hObject    handle to colorpointcloudcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of colorpointcloudcheck


% --- Executes on button press in mirrorsides.
function mirrorsides_Callback(hObject, eventdata, handles)
% hObject    handle to mirrorsides (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of mirrorsides


% --- Executes on button press in showdiscfibers.
function showdiscfibers_Callback(hObject, eventdata, handles)
% hObject    handle to showdiscfibers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showdiscfibers


% --- Executes on button press in save3dsetting.
function save3dsetting_Callback(hObject, eventdata, handles)
% hObject    handle to save3dsetting (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
leadfigure = getappdata(handles.lg_3dsetting, 'leadfigure');

M = getappdata(leadfigure, 'M');
M.ui.elmodelselect = get(handles.elmodelselect, 'Value');
M.ui.elrendering = get(handles.elrenderingpopup, 'Value');
M.ui.isovscloudpopup = get(handles.isovscloudpopup, 'Value');
M.ui.colorpointcloudcheck = get(handles.colorpointcloudcheck, 'Value');
M.ui.mirrorsides = get(handles.mirrorsides, 'Value');
M.ui.showdiscfibers = get(handles.showdiscfibers, 'Value');

setappdata(leadfigure, 'M', M);

prefs=ea_prefs('');
discfibers = prefs.machine.lg.discfibers;
switch get(get(handles.showfiberssetpanel, 'SelectedObject'), 'Tag')
    case 'showposonly'
        discfibers.showfibersset = 'positive';
    case 'shownegonly'
        discfibers.showfibersset = 'negative';
    case 'showboth'
        discfibers.showfibersset = 'both';
end
discfibers.pospredthreshold = str2double(get(handles.pospredthreshold,'String'));
discfibers.negpredthreshold = str2double(get(handles.negpredthreshold,'String'));
discfibers.statmetric = get(handles.statmetric,'Value');
ea_setprefs('lg.discfibers', discfibers);

delete(handles.lg_3dsetting);


function pospredthreshold_Callback(hObject, eventdata, handles)
% hObject    handle to pospredthreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pospredthreshold as text
%        str2double(get(hObject,'String')) returns contents of pospredthreshold as a double


% --- Executes during object creation, after setting all properties.
function pospredthreshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pospredthreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function negpredthreshold_Callback(hObject, eventdata, handles)
% hObject    handle to negpredthreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of negpredthreshold as text
%        str2double(get(hObject,'String')) returns contents of negpredthreshold as a double


% --- Executes during object creation, after setting all properties.
function negpredthreshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to negpredthreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in statmetric.
function statmetric_Callback(hObject, eventdata, handles)
% hObject    handle to statmetric (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns statmetric contents as cell array
%        contents{get(hObject,'Value')} returns selected item from statmetric


% --- Executes during object creation, after setting all properties.
function statmetric_CreateFcn(hObject, eventdata, handles)
% hObject    handle to statmetric (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in showposonly.
function showposonly_Callback(hObject, eventdata, handles)
% hObject    handle to showposonly (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showposonly
set(handles.pospredthreshold, 'Enable', 'on');
set(handles.negpredthreshold, 'Enable', 'off');


% --- Executes on button press in shownegonly.
function shownegonly_Callback(hObject, eventdata, handles)
% hObject    handle to shownegonly (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of shownegonly
set(handles.pospredthreshold, 'Enable', 'off');
set(handles.negpredthreshold, 'Enable', 'on');


% --- Executes on button press in showboth.
function showboth_Callback(hObject, eventdata, handles)
% hObject    handle to showboth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showboth
set(handles.pospredthreshold, 'Enable', 'on');
set(handles.negpredthreshold, 'Enable', 'on');
