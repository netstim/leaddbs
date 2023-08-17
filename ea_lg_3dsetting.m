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

% Last Modified by GUIDE v2.5 24-Mar-2020 09:01:32

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
set(handles.elmodelselect,'String',[{'Patient specified'};ea_resolve_elspec]);
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

ea_ListBoxRenderer(handles.elmodelselect);


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
                if isfile([options.root,options.patientname,filesep,'reconstruction',filesep,options.patientname,'_desc-reconstruction.mat'])
                    [~,~,markers,elmodel,~,coords_acpc] = ea_load_reconstruction(options);

                    if M.ui.elmodelselect == 1 % use patient specific elmodel
                        if exist('elmodel','var')
                            M.elstruct(pt).elmodel = elmodel;
                        else % use default for older reconstructions that did not store elmodel.
                            M.elstruct(pt).elmodel = 'Medtronic 3389';
                        end
                    else
                        elmodels = [{'Patient specified'};ea_resolve_elspec];
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
                            [xunitv, yunitv] = ea_calcxy(markers(side).head, markers(side).tail);
                            markers(side).x = coords_mm{side}(1,:) + xunitv*(options.elspec.lead_diameter/2);
                            markers(side).y = coords_mm{side}(1,:) + yunitv*(options.elspec.lead_diameter/2);
                        end
                    end
                    M.elstruct(pt).markers = markers;
                end
            catch
                if pt>1 % first patient has worked but some other patient seems not to have worked.
                    try
                        M.elstruct(1).coords_mm; % probe if error happens in pt. 1 ? if not show warning
                        warning(['No reconstruction present for ',patientname,'. Please check.']);
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

setappdata(leadfigure, 'M', M);

delete(handles.lg_3dsetting);


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


% --- Executes on button press in regressorcolormap.
function regressorcolormap_Callback(hObject, eventdata, handles)
% hObject    handle to regressorcolormap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
leadfigure = getappdata(handles.lg_3dsetting, 'leadfigure');

M = getappdata(leadfigure, 'M');
if isfield(M.ui, 'regressorcolormap')
    M.ui.regressorcolormap = ea_select_colormap(length(gray),'custom2',M.ui.regressorcolormap);
else
    M.ui.regressorcolormap = ea_select_colormap(length(gray),'custom2');
end
setappdata(leadfigure, 'M', M);
