function varargout = ea_edit_regressor(varargin)
% EA_EDIT_REGRESSOR MATLAB code for ea_edit_regressor.fig
%      EA_EDIT_REGRESSOR, by itself, creates a new EA_EDIT_REGRESSOR or raises the existing
%      singleton*.
%
%      H = EA_EDIT_REGRESSOR returns the handle to a new EA_EDIT_REGRESSOR or the handle to
%      the existing singleton*.
%
%      EA_EDIT_REGRESSOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EA_EDIT_REGRESSOR.M with the given input arguments.
%
%      EA_EDIT_REGRESSOR('Property','Value',...) creates a new EA_EDIT_REGRESSOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ea_edit_regressor_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ea_edit_regressor_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ea_edit_regressor

% Last Modified by GUIDE v2.5 11-Oct-2015 15:47:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @ea_edit_regressor_OpeningFcn, ...
    'gui_OutputFcn',  @ea_edit_regressor_OutputFcn, ...
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


% --- Executes just before ea_edit_regressor is made visible.
function ea_edit_regressor_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ea_edit_regressor (see VARARGIN)

set(gcf,'Name','Edit Regressor');

% Choose default command line output for ea_edit_regressor
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

M=varargin{1};

setappdata(handles.editregressor, 'maxnumel', max([M.S.numContacts]));

[~,ptnames]=cellfun(@fileparts,M.patient.list,'UniformOutput',0);
set(handles.datatable,'RowName',ptnames);

try
    regressor=M.clinical.vars{M.ui.clinicallist};
catch % new variable
    regressor=nan(length(ptnames),1);
end

if iscell(regressor)
   regressor=cell2mat(regressor);
end

set(handles.datatable,'Data',regressor);
set(handles.datatable,'ColumnEditable',true(1,size(regressor,2)));
if isperpatient(regressor)
    switchperpatient(handles);
elseif isperhemisphere(regressor)
    switchperhemisphere(handles);
elseif ispercontactpair(regressor)
    switchpercontactpair(handles);
elseif ispercontact(regressor)
    switchpercontact(handles);
end

try
    set(handles.varname,'String',M.clinical.labels{M.ui.clinicallist});
end

% UIWAIT makes ea_edit_regressor wait for user response (see UIRESUME)
uiwait(hObject);


function switchperpatient(handles)
reg=get(handles.datatable,'Data');
if ~isempty(reg)
   if ~isperpatient(reg)
       answ=questdlg('Warning: switching variable type will delete/modify variable! Are you sure you want this?','Warning','Yes','No','No');
   switch answ
       case 'No'
           return
   end
   end
end
set(handles.variabletype, 'Value', 1);
set(handles.datatable,'ColumnName',{'Value'});

reg=reg(1:size(reg,1),1);
set(handles.datatable,'Data',reg);
set(handles.datatable,'ColumnEditable',true(1,1));

function switchperhemisphere(handles)
reg=get(handles.datatable,'Data');
if ~isempty(reg)
   if ~isperhemisphere(reg)
       answ=questdlg('Warning: switching variable type will delete/modify variable! Are you sure you want this?','Warning','Yes','No','No');
   switch answ
       case 'No'
           return
   end
   end
end
set(handles.variabletype, 'Value', 2);
set(handles.datatable,'ColumnName',{'Right Hem.','Left Hem.'});
try
    reg=reg(1:size(reg,1),1:2);
catch
   reg=[reg(1:size(reg,1),1),nan(size(reg,1),1)];
end
set(handles.datatable,'Data',reg);
set(handles.datatable,'ColumnEditable',true(1,1));


function switchpercontact(handles)
reg=get(handles.datatable,'Data');
if ~isempty(reg)
   if ~ispercontact(reg)
       answ=questdlg('Warning: switching variable type will delete/modify variable! Are you sure you want this?','Warning','Yes','No','No');
   switch answ
       case 'No'
           return
   end
   end
end

set(handles.variabletype, 'Value', 3);

maxnumel = getappdata(handles.editregressor, 'maxnumel');
conInd = arrayfun(@num2str, 1:maxnumel, 'Uni', 0);
set(handles.datatable,'ColumnName',[strcat('k', conInd, 'R'), strcat('k', conInd, 'L')]);

nreg=nan(size(reg,1),8);
nreg(1:size(reg,1),1:size(reg,2))=reg;

reg=nreg; clear('nreg');
set(handles.datatable,'Data',reg);
set(handles.datatable,'ColumnEditable',true(1,8));


function switchpercontactpair(handles)
reg=get(handles.datatable,'Data');
if ~isempty(reg)
   if ~ispercontactpair(reg)
       answ=questdlg('Warning: switching variable type will delete/modify variable! Are you sure you want this?','Warning','Yes','No','No');
   switch answ
       case 'No'
           return
   end
   end
end

set(handles.variabletype, 'Value', 4);

maxnumel = getappdata(handles.editregressor, 'maxnumel');
comb = nchoosek(1:maxnumel, 2);
conComb = cellstr(strcat(string(comb(:,1)), '-', string(comb(:,2))))';
set(handles.datatable,'ColumnName',[strcat('K',conComb,'R'),strcat('K',conComb,'L')]);

nreg=nan(size(reg,1),6);
try
    nreg(1:size(reg,1),1:size(nreg,2))=reg(1:size(reg,1),1:size(nreg,2));
catch
    nreg(1:size(reg,1),1:size(reg,2))=reg(1:size(reg,1),1:size(reg,2));
end
reg=nreg; clear('nreg');
set(handles.datatable,'Data',reg);
set(handles.datatable,'ColumnEditable',true(1,6));


function yn=ispercontactpair(regressor)
yn=(iscell(regressor) && size(regressor{1},2)==3) || (~iscell(regressor) && size(regressor,2)==6);
function yn=ispercontact(regressor)
yn=(iscell(regressor) && size(regressor{1},2)==4) || (~iscell(regressor) && size(regressor,2)==8);
function yn=isperpatient(regressor)
yn=~iscell(regressor) && size(regressor,2)==1;
function yn=isperhemisphere(regressor)
yn=~iscell(regressor) && size(regressor,2)==2;

% --- Outputs from this function are returned to the command line.
function varargout = ea_edit_regressor_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

if getappdata(gcf,'save')
    varargout{1} = get(handles.datatable,'Data');
    varargout{2}=get(handles.varname,'String');
else
    varargout{1}=[];
    varargout{2}=[];
end

delete(hObject);

% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setappdata(gcf,'save',1);
close(gcf);

% --- Executes on button press in cancel.
function cancel_Callback(hObject, eventdata, handles)
% hObject    handle to cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setappdata(gcf,'save',0);
close(gcf);


% --- Executes when user attempts to close editregressor.
function editregressor_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to editregressor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
%delete(hObject);
uiresume(hObject);



function varname_Callback(hObject, eventdata, handles)
% hObject    handle to varname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of varname as text
%        str2double(get(hObject,'String')) returns contents of varname as a double


% --- Executes during object creation, after setting all properties.
function varname_CreateFcn(hObject, eventdata, handles)
% hObject    handle to varname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in variabletype.
function variabletype_Callback(hObject, eventdata, handles)
% hObject    handle to variabletype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns variabletype contents as cell array
%        contents{get(hObject,'Value')} returns selected item from variabletype
popvals=get(hObject,'String');
variabletype=popvals{get(hObject,'Value')};
switch variabletype
    case '1 variable per patient'
        switchperpatient(handles);
    case '1 variable per hemisphere'
        switchperhemisphere(handles);
    case '1 variable per contact'
        switchpercontact(handles);
    case '1 variable per contact pair'
        switchpercontactpair(handles);
end




% --- Executes during object creation, after setting all properties.
function variabletype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to variabletype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
