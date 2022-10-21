function varargout = ea_lg_stats(varargin)
% EA_LG_STATS MATLAB code for ea_lg_stats.fig
%      EA_LG_STATS, by itself, creates a new EA_LG_STATS or raises the existing
%      singleton*.
%
%      H = EA_LG_STATS returns the handle to a new EA_LG_STATS or the handle to
%      the existing singleton*.
%
%      EA_LG_STATS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EA_LG_STATS.M with the given input arguments.
%
%      EA_LG_STATS('Property','Value',...) creates a new EA_LG_STATS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ea_lg_stats_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ea_lg_stats_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ea_lg_stats

% Last Modified by GUIDE v2.5 13-Nov-2020 13:12:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ea_lg_stats_OpeningFcn, ...
                   'gui_OutputFcn',  @ea_lg_stats_OutputFcn, ...
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


% --- Executes just before ea_lg_stats is made visible.
function ea_lg_stats_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ea_lg_stats (see VARARGIN)

% Choose default command line output for ea_lg_stats
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ea_lg_stats wait for user response (see UIRESUME)
% uiwait(handles.lg_stats);

set(handles.vilist,'Max',100,'Min',0);
set(handles.fclist,'Max',100,'Min',0);

leadfigure = varargin{1};
setappdata(handles.lg_stats, 'leadfigure', leadfigure);
M = getappdata(leadfigure, 'M');

if ~isempty(M.clinical.labels)
    set(handles.clinicalvars, 'String', M.clinical.labels);
end

% Disable group comparison button for variables per hemisphere or
% non-logical variable, or when no clinical variable exists
if isempty(M.clinical.vars)
    set(handles.clinicalvars, 'Enable', 'off');
    set(handles.corrbutton_vta, 'Enable', 'off');
    set(handles.ttestbutton_vta, 'Enable', 'off');
    set(handles.fclist, 'Enable', 'off');
    set(handles.corrbutton_ft, 'Enable', 'off');
    set(handles.ttestbutton_ft, 'Enable', 'off');
else
    clinicvar = M.clinical.vars{get(handles.clinicalvars, 'Value')};
    if size(clinicvar,2)==2 || numel(unique(clinicvar(~isnan(clinicvar))))~=2
        set(handles.ttestbutton_vta, 'Enable', 'off');
        set(handles.ttestbutton_ft, 'Enable', 'off');
    end
end

% refresh UI
if ~isempty(M.vilist)
    set(handles.vilist, 'String', M.vilist);
end

if ~isempty(M.fclist)
    set(handles.fclist, 'String', M.fclist);
end

% refresh selections on VI and FC Lists:
try
    vlist = get(handles.vilist, 'String');
    if max(M.ui.volumeintersections) > length(vlist) || ...
        (ischar(vlist) && strcmp(vlist, 'subcortical_regions'))
        set(handles.vilist,'Value',1);
    else
        set(handles.vilist,'Value',M.ui.volumeintersections);
    end
end

try
    fclist = get(handles.fclist, 'String');
    if M.ui.fibercounts > length(fclist) || ...
        (ischar(fclist) && strcmp(fclist, 'whole-brain_regions'))
        set(handles.fclist,'Value',1);
    else
        set(handles.fclist,'Value',M.ui.fibercounts);
    end
end


% --- Outputs from this function are returned to the command line.
function varargout = ea_lg_stats_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in vilist.
function vilist_Callback(hObject, eventdata, handles)
% hObject    handle to vilist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns vilist contents as cell array
%        contents{get(hObject,'Value')} returns selected item from vilist
leadfigure = getappdata(handles.lg_stats, 'leadfigure');
M = getappdata(leadfigure, 'M');

M.ui.volumeintersections = get(handles.vilist,'Value');
setappdata(leadfigure, 'M', M);


% --- Executes during object creation, after setting all properties.
function vilist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vilist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
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
leadfigure = getappdata(handles.lg_stats, 'leadfigure');
M = getappdata(leadfigure, 'M');

M.ui.fibercounts = get(handles.fclist,'Value');
setappdata(leadfigure, 'M', M);


% --- Executes during object creation, after setting all properties.
function fclist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fclist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in corrbutton_vta.
function corrbutton_vta_Callback(hObject, eventdata, handles)
% hObject    handle to corrbutton_vta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
leadfigure = getappdata(handles.lg_stats, 'leadfigure');
ea_busyaction('on', leadfigure, 'group');

stats=preparedataanalysis_vta(handles);
assignin('base','stats',stats);

M=getappdata(leadfigure,'M');

colors = M.groups.color(unique(M.patient.group(M.ui.listselect)),:);
% Colors in Lead Group are the same, fallback to gramm default colormap
if size(colors,1)>1 && size(unique(colors, 'rows'),1) == 1
    colors = [1.0000  0.3686  0.4118];
end

% perform correlations:
if size(stats.corrcl,2)==1 % one value per patient
    if ~isempty(stats.vicorr.both)
        ea_corrplot(stats.corrcl,stats.vicorr.both,'permutation',...
            [{'VTA intersection, both hemispheres'},stats.vc_labels],...
            M.patient.group(M.ui.listselect),'',colors);
        ea_corrplot(stats.corrcl,stats.vicorr.nboth,'permutation',...
            [{'Normalized VTA intersection, both hemispheres'},stats.vc_labels],...
            M.patient.group(M.ui.listselect),'',colors);
    end

    if ~isempty(stats.vicorr.right)
        ea_corrplot(stats.corrcl,stats.vicorr.right,'permutation',...
            [{'VTA intersection, right hemisphere'},stats.vc_labels],...
            M.patient.group(M.ui.listselect),'',colors);
        ea_corrplot(stats.corrcl,stats.vicorr.nright,'permutation',...
            [{'Normalized VTA intersection, right hemisphere'},stats.vc_labels],...
            M.patient.group(M.ui.listselect),'',colors);
    end

    if ~isempty(stats.vicorr.left)
        ea_corrplot(stats.corrcl,stats.vicorr.left,'permutation',...
            [{'VTA intersection, left hemisphere'},stats.vc_labels],...
            M.patient.group(M.ui.listselect),'',colors);
        ea_corrplot(stats.corrcl,stats.vicorr.nleft,'permutation',...
            [{'Normalized VTA intersection, left hemisphere'},stats.vc_labels],...
            M.patient.group(M.ui.listselect),'',colors);
    end
elseif size(stats.corrcl,2)==2 % one value per hemisphere
    if ~isempty(stats.vicorr.both)
        ea_corrplot(stats.corrcl(:),[stats.vicorr.right;stats.vicorr.left],'permutation',...
            [{'VTA intersection, both hemispheres'},stats.vc_labels],...
            repmat(M.patient.group(M.ui.listselect),2,1),'',colors);
        ea_corrplot(stats.corrcl(:),[stats.vicorr.nright;stats.vicorr.nleft],'permutation',...
            [{'Normalized VTA intersection, both hemispheres'},stats.vc_labels],...
            repmat(M.patient.group(M.ui.listselect),2,1),'',colors);
    end

    if ~isempty(stats.vicorr.right)
        ea_corrplot(stats.corrcl(:,1),stats.vicorr.right,'permutation',...
            [{'VTA intersection, right hemisphere'},stats.vc_labels],...
            M.patient.group(M.ui.listselect),'',colors);
        ea_corrplot(stats.corrcl(:,1),stats.vicorr.nright,'permutation',...
            [{'Normalized VTA intersection, right hemispheres'},stats.vc_labels],...
            M.patient.group(M.ui.listselect),'',colors);
    end

    if ~isempty(stats.vicorr.left)
        ea_corrplot(stats.corrcl(:,2),stats.vicorr.left,'permutation',...
            [{'VTA intersection, left hemisphere'},stats.vc_labels],...
            M.patient.group(M.ui.listselect),'',colors);
        ea_corrplot(stats.corrcl(:,2),stats.vicorr.nleft,'permutation',...
            [{'Normalized VTA intersection, left hemispheres'},stats.vc_labels],...
            M.patient.group(M.ui.listselect),'',colors);
    end
else
    ea_error('Please select a regressor with one value per patient or per hemisphere to perform this correlation.');
end
ea_busyaction('off', leadfigure, 'group');


% --- Executes on button press in ttestbutton_vta.
function ttestbutton_vta_Callback(hObject, eventdata, handles)
% hObject    handle to ttestbutton_vta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
leadfigure = getappdata(handles.lg_stats, 'leadfigure');
ea_busyaction('on', leadfigure, 'group');

stats = preparedataanalysis_vta(handles);
assignin('base','stats',stats);

% perform t-tests:
if ~isempty(stats.vicorr.both)
    ea_ttest(stats.vicorr.both(~repmat(logical(stats.corrcl),1,size(stats.vicorr.both,2))),...
        stats.vicorr.both(repmat(logical(stats.corrcl),1,size(stats.vicorr.both,2))),...
        'VTA intersection, both hemispheres',stats.vc_labels{end});
end

if ~isempty(stats.vicorr.right)
    ea_ttest(stats.vicorr.right(~repmat(logical(stats.corrcl),1,size(stats.vicorr.right,2))),...
        stats.vicorr.right(repmat(logical(stats.corrcl),1,size(stats.vicorr.right,2))),...
        'VTA intersection, right hemisphere',stats.vc_labels{end});
end

if ~isempty(stats.vicorr.left)
    ea_ttest(stats.vicorr.left(~repmat(logical(stats.corrcl),1,size(stats.vicorr.left,2))),...
        stats.vicorr.left(repmat(logical(stats.corrcl),1,size(stats.vicorr.left,2))),...
        'VTA intersection, left hemisphere',stats.vc_labels{end});
end

if ~isempty(stats.vicorr.nboth)
    ea_ttest(stats.vicorr.nboth(~repmat(logical(stats.corrcl),1,size(stats.vicorr.both,2))),...
        stats.vicorr.nboth(repmat(logical(stats.corrcl),1,size(stats.vicorr.both,2))),...
        'Normalized VTA intersection, both hemispheres',stats.vc_labels{end});
end

if ~isempty(stats.vicorr.nright)
    ea_ttest(stats.vicorr.nright(~repmat(logical(stats.corrcl),1,size(stats.vicorr.right,2))),...
        stats.vicorr.nright(repmat(logical(stats.corrcl),1,size(stats.vicorr.right,2))),...
        'Normalized VTA intersection, right hemisphere',stats.vc_labels{end});
end

if ~isempty(stats.vicorr.nleft)
    ea_ttest(stats.vicorr.nleft(~repmat(logical(stats.corrcl),1,size(stats.vicorr.left,2))),...
        stats.vicorr.nleft(repmat(logical(stats.corrcl),1,size(stats.vicorr.left,2))),...
        'Normalized VTA intersection, left hemisphere',stats.vc_labels{end});
end

ea_busyaction('off', leadfigure, 'group');


% --- Executes on button press in corrbutton_ft.
function corrbutton_ft_Callback(hObject, eventdata, handles)
% hObject    handle to corrbutton_ft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
leadfigure = getappdata(handles.lg_stats, 'leadfigure');
ea_busyaction('on', leadfigure, 'group');

stats = preparedataanalysis_ft(handles);
assignin('base','stats',stats);

M=getappdata(leadfigure,'M');

colors = M.groups.color(unique(M.patient.group(M.ui.listselect)),:);
% Colors in Lead Group are the same, fallback to gramm default colormap
if size(colors,1)>1 && size(unique(colors, 'rows'),1) == 1
    colors = [];
end

% perform correlations:
if size(stats.corrcl,2)==1 % one value per patient
    if ~isempty(stats.fccorr.both)
        ea_corrplot(stats.corrcl,stats.fccorr.both,'permutation',...
            [{'Fiber connection, both hemispheres'},stats.fc_labels],...
            M.patient.group(M.ui.listselect),'',colors);
        ea_corrplot(stats.corrcl,stats.fccorr.nboth,'permutation',...
            [{'Normalized Fiber connection, both hemispheres'},stats.fc_labels],...
            M.patient.group(M.ui.listselect),'',colors);
    end

    if ~isempty(stats.fccorr.right)
        ea_corrplot(stats.corrcl,stats.fccorr.right,'permutation',...
            [{'Fiber connection, right hemispheres'},stats.fc_labels],...
            M.patient.group(M.ui.listselect),'',colors);
        ea_corrplot(stats.corrcl,stats.fccorr.nright,'permutation',...
            [{'Normalized Fiber connection, right hemispheres'},stats.fc_labels],...
            M.patient.group(M.ui.listselect),'',colors);
    end

    if ~isempty(stats.fccorr.left)
        ea_corrplot(stats.corrcl,stats.fccorr.left,'permutation',...
            [{'Fiber connection, left hemispheres'},stats.fc_labels],...
            M.patient.group(M.ui.listselect),'',colors);
        ea_corrplot(stats.corrcl,stats.fccorr.nleft,'permutation',...
            [{'Normalized Fiber connection, left hemispheres'},stats.fc_labels],...
            M.patient.group(M.ui.listselect),'',colors);
    end
elseif size(stats.corrcl,2)==2 % one value per hemisphere
    if ~isempty(stats.fccorr.both)
        ea_corrplot(stats.corrcl(:),[stats.fccorr.right;stats.fccorr.left],'permutation',...
            [{'Fiber connection, both hemispheres'},stats.fc_labels],...
            repmat(M.patient.group(M.ui.listselect),2,1),'',colors);
        ea_corrplot(stats.corrcl(:),[stats.fccorr.nright;stats.fccorr.nleft],'permutation',...
            [{'Normalized Fiber connection, both hemispheres'},stats.fc_labels],...
            repmat(M.patient.group(M.ui.listselect),2,1),'',colors);
    end

    if ~isempty(stats.fccorr.right)
        ea_corrplot(stats.corrcl(:,1),stats.fccorr.right,'permutation',...
            [{'Fiber connection, right hemispheres'},stats.fc_labels],...
            M.patient.group(M.ui.listselect),'',colors);
        ea_corrplot(stats.corrcl(:,1),stats.fccorr.nright,'permutation',...
            [{'Normalized Fiber connection, right hemispheres'},stats.fc_labels],...
            M.patient.group(M.ui.listselect),'',colors);
    end

    if ~isempty(stats.fccorr.left)
        ea_corrplot(stats.corrcl(:,2),stats.fccorr.left,'permutation',...
            [{'Fiber connection, left hemispheres'},stats.fc_labels],...
            M.patient.group(M.ui.listselect),'',colors);
        ea_corrplot(stats.corrcl(:,2),stats.fccorr.nleft,'permutation',...
            [{'Normalized Fiber connection, left hemispheres'},stats.fc_labels],...
            M.patient.group(M.ui.listselect),'',colors);
    end
else
    ea_error('Please select a regressor with one value per patient or per hemisphere to perform this correlation.');
end

ea_busyaction('off', leadfigure, 'group');


% --- Executes on button press in ttestbutton_ft.
function ttestbutton_ft_Callback(hObject, eventdata, handles)
% hObject    handle to ttestbutton_ft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
leadfigure = getappdata(handles.lg_stats, 'leadfigure');
ea_busyaction('on', leadfigure, 'group');

stats = preparedataanalysis_ft(handles);
assignin('base','stats',stats);

% perform t-tests:
if ~isempty(stats.fccorr.both)
    ea_ttest(stats.fccorr.both(repmat(logical(stats.corrcl),1,size(stats.fccorr.both,2))),...
        stats.fccorr.both(~repmat(logical(stats.corrcl),1,size(stats.fccorr.both,2))),...
        'Fiber connection, both hemispheres',stats.fc_labels{end});
end

if ~isempty(stats.fccorr.nboth)
    ea_ttest(stats.fccorr.nboth(repmat(logical(stats.corrcl),1,size(stats.fccorr.nboth,2))),...
        stats.fccorr.nboth(~repmat(logical(stats.corrcl),1,size(stats.fccorr.nboth,2))),...
        'Normalized Fiber connection, both hemispheres',stats.fc_labels{end});
end

if ~isempty(stats.fccorr.right)
    ea_ttest(stats.fccorr.right(repmat(logical(stats.corrcl),1,size(stats.fccorr.right,2))),...
        stats.fccorr.right(~repmat(logical(stats.corrcl),1,size(stats.fccorr.right,2))),...
        'Fiber connection, right hemispheres',stats.fc_labels{end});
end

if ~isempty(stats.fccorr.nright)
    ea_ttest(stats.fccorr.nright(repmat(logical(stats.corrcl),1,size(stats.fccorr.nright,2))),...
        stats.fccorr.nright(~repmat(logical(stats.corrcl),1,size(stats.fccorr.nright,2))),...
        'Normalized Fiber connection, right hemispheres',stats.fc_labels{end});
end

if ~isempty(stats.fccorr.left)
    ea_ttest(stats.fccorr.left(repmat(logical(stats.corrcl),1,size(stats.fccorr.left,2))),...
        stats.fccorr.left(~repmat(logical(stats.corrcl),1,size(stats.fccorr.left,2))),...
        'Fiber connection, left hemispheres',stats.fc_labels{end});
end

if ~isempty(stats.fccorr.nleft)
    ea_ttest(stats.fccorr.nleft(repmat(logical(stats.corrcl),1,size(stats.fccorr.nleft,2))),...
        stats.fccorr.nleft(~repmat(logical(stats.corrcl),1,size(stats.fccorr.nleft,2))),...
        'Normalized Fiber connection, left hemispheres',stats.fc_labels{end});
end

ea_busyaction('off', leadfigure, 'group');

% --- Executes on selection change in VTAvsEfield.
function VTAvsEfield_Callback(hObject, eventdata, handles)
% hObject    handle to VTAvsEfield (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns VTAvsEfield contents as cell array
%        contents{get(hObject,'Value')} returns selected item from VTAvsEfield


% --- Executes during object creation, after setting all properties.
function VTAvsEfield_CreateFcn(hObject, eventdata, handles)
% hObject    handle to VTAvsEfield (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function [stats] = preparedataanalysis_vta(handles)
leadfigure = getappdata(handles.lg_stats, 'leadfigure');
M=getappdata(leadfigure, 'M');

% Get VTA intersection:
vicnt=1; ptcnt=1;

howmanyvis=length(get(handles.vilist,'Value'));

patientlist = leadfigure.findobj('Tag','patientlist');
howmanypts=length(get(patientlist,'Value'));

vicorr_right=zeros(howmanypts,howmanyvis);
vicorr_left=zeros(howmanypts,howmanyvis);
vicorr_both=zeros(howmanypts,howmanyvis);
nvicorr_right=zeros(howmanypts,howmanyvis);
nvicorr_left=zeros(howmanypts,howmanyvis);
nvicorr_both=zeros(howmanypts,howmanyvis);
vc_labels={};

switch get(handles.VTAvsEfield,'value')
    case 1 % VTA
        vtavsefield='vat';
        vtavsefieldLabel='VTA';
    case 2 % E-field
        vtavsefield='efield';
        vtavsefieldLabel = 'E-field';
end

for vi=get(handles.vilist,'Value') % get volume interactions for each patient from stats
    for pt=get(patientlist,'Value')
        S.label=['gs_',M.guid];
        try
            [ea_stats,usewhichstim]=ea_assignstimcnt(M.stats(pt).ea_stats,S);

            for side=1:size(M.stats(pt).ea_stats.stimulation(usewhichstim).vat,1)
                for vat=1
                    if side==1 % right hemisphere
                        vicorr_right(ptcnt,vicnt)=vicorr_right(ptcnt,vicnt)+M.stats(pt).ea_stats.stimulation(usewhichstim).(vtavsefield)(side,vat).AtlasIntersection(vi);
                        nvicorr_right(ptcnt,vicnt)=nvicorr_right(ptcnt,vicnt)+M.stats(pt).ea_stats.stimulation(usewhichstim).(vtavsefield)(side,vat).nAtlasIntersection(vi);

                    elseif side==2 % left hemisphere
                        vicorr_left(ptcnt,vicnt)=vicorr_left(ptcnt,vicnt)+M.stats(pt).ea_stats.stimulation(usewhichstim).(vtavsefield)(side,vat).AtlasIntersection(vi);
                        nvicorr_left(ptcnt,vicnt)=nvicorr_left(ptcnt,vicnt)+M.stats(pt).ea_stats.stimulation(usewhichstim).(vtavsefield)(side,vat).nAtlasIntersection(vi);
                    end
                    vicorr_both(ptcnt,vicnt)=vicorr_both(ptcnt,vicnt)+M.stats(pt).ea_stats.stimulation(usewhichstim).(vtavsefield)(side,vat).AtlasIntersection(vi);
                    nvicorr_both(ptcnt,vicnt)=nvicorr_both(ptcnt,vicnt)+M.stats(pt).ea_stats.stimulation(usewhichstim).(vtavsefield)(side,vat).nAtlasIntersection(vi);
                end
            end
        catch
            ea_error(['DBS stats for patient ',M.patient.list{pt},' need to be calculated.']);
        end

        % check if all three values have been served. if not, set to zero
        % (e.g. if there was no stimulation at all on one hemisphere, this
        % could happen.
        ptcnt=ptcnt+1;
    end
    vc_labels{end+1}=[ea_stripext(M.stats(pt).ea_stats.atlases.names{vi}),': ',vtavsefieldLabel,' intersection'];

    ptcnt=1;
    vicnt=vicnt+1;
end

% prepare outputs:
vicorr.both=vicorr_both;
vicorr.left=vicorr_left;
vicorr.right=vicorr_right;
vicorr.nboth=nvicorr_both;
vicorr.nleft=nvicorr_left;
vicorr.nright=nvicorr_right;

% clinical vector:
corrcl=M.clinical.vars{get(handles.clinicalvars,'Value')};
corrcl=corrcl(get(patientlist,'Value'),:);

clinstrs = get(handles.clinicalvars,'String');
vc_labels = [clinstrs(get(handles.clinicalvars,'Value')),vc_labels]; % add name of clinical vector to labels

stats.corrcl=corrcl;
stats.vicorr=vicorr;
stats.vc_labels=vc_labels;


function [stats]=preparedataanalysis_ft(handles)
leadfigure = getappdata(handles.lg_stats, 'leadfigure');
M=getappdata(leadfigure, 'M');

% Get VTA intersection:
vicnt=1; ptcnt=1;

patientlist = leadfigure.findobj('Tag','patientlist');
howmanypts=length(get(patientlist,'Value'));

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
    for pt=get(patientlist,'Value')
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
    fc_labels{end+1} = [M.stats(pt).ea_stats.stimulation(usewhichstim).ft(1).labels{1}{fc},':  fiber connection'];
end

% prepare outputs:
fccorr.both=fccorr_both;
fccorr.nboth=nfccorr_both;
fccorr.right=fccorr_right;
fccorr.nright=nfccorr_right;
fccorr.left=fccorr_left;
fccorr.nleft=nfccorr_left;

% clinical vector:
corrcl=M.clinical.vars{get(handles.clinicalvars,'Value')};
corrcl=corrcl(get(patientlist,'Value'),:);

clinstrs = get(handles.clinicalvars,'String');
fc_labels=[clinstrs(get(handles.clinicalvars,'Value')),fc_labels]; % add name of clinical vector to labels

stats.corrcl=corrcl;
stats.fccorr=fccorr;
stats.fc_labels=fc_labels;


% --- Executes on selection change in clinicalvars.
function clinicalvars_Callback(hObject, eventdata, handles)
% hObject    handle to clinicalvars (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns clinicalvars contents as cell array
%        contents{get(hObject,'Value')} returns selected item from clinicalvars

leadfigure = getappdata(handles.lg_stats, 'leadfigure');
M = getappdata(leadfigure, 'M');

% Disable group comparison button for variables per hemisphere or
% non-logical variable
clinicvar = M.clinical.vars{get(handles.clinicalvars, 'Value')};
if size(clinicvar,2)==2 || numel(unique(clinicvar(~isnan(clinicvar))))~=2
    set(handles.ttestbutton_vta, 'Enable', 'off');
    set(handles.ttestbutton_ft, 'Enable', 'off');
else
    set(handles.ttestbutton_vta, 'Enable', 'on');
    set(handles.ttestbutton_ft, 'Enable', 'on');
end


% --- Executes during object creation, after setting all properties.
function clinicalvars_CreateFcn(hObject, eventdata, handles)
% hObject    handle to clinicalvars (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in targetreport.
function targetreport_Callback(hObject, eventdata, handles)
% hObject    handle to targetreport (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
leadfigure = getappdata(handles.lg_stats, 'leadfigure');
M = getappdata(leadfigure, 'M');
ea_gentargetreport(M);
