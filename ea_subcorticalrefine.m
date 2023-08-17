function varargout = ea_subcorticalrefine(varargin)
% EA_SUBCORTICALREFINE MATLAB code for ea_subcorticalrefine.fig
%      EA_SUBCORTICALREFINE, by itself, creates a new EA_SUBCORTICALREFINE or raises the existing
%      singleton*.
%
%      H = EA_SUBCORTICALREFINE returns the handle to a new EA_SUBCORTICALREFINE or the handle to
%      the existing singleton*.
%
%      EA_SUBCORTICALREFINE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EA_SUBCORTICALREFINE.M with the given input arguments.
%
%      EA_SUBCORTICALREFINE('Property','Value',...) creates a new EA_SUBCORTICALREFINE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ea_subcorticalrefine_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ea_subcorticalrefine_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ea_subcorticalrefine

% Last Modified by GUIDE v2.5 30-Jan-2019 09:00:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ea_subcorticalrefine_OpeningFcn, ...
                   'gui_OutputFcn',  @ea_subcorticalrefine_OutputFcn, ...
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


% --- Executes just before ea_subcorticalrefine is made visible.
function ea_subcorticalrefine_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ea_subcorticalrefine (see VARARGIN)

options = varargin{1};

setappdata(handles.scrf, 'options', options);
handles.patientname.String = options.subj.subjId;

set(handles.scrf, 'name', ['Brainshift Correction: ', options.subj.subjId]);
options.init = 1;
ispresent = ea_refreshscrf(options,handles);

switch options.scrf.mask
    case 'No mask'
        handles.mask0.Value = 1;
        handles.mask1.Value = 0;
        handles.mask2.Value = 0;
        if ~ispresent || isfield(options, 'autobrainshift')
            ea_compute_scrf(handles)
        end
    case 'Coarse + Fine mask (Schönecker 2008)'
        handles.mask0.Value = 0;
        handles.mask1.Value = 0;
        handles.mask2.Value = 1;
        if ~ispresent || isfield(options, 'autobrainshift')
            ea_compute_scrf(handles)
        end
    otherwise % 'Coarse mask (Schönecker 2008)'
        handles.mask0.Value = 0;
        handles.mask1.Value = 1;
        handles.mask2.Value = 0;
        if ~ispresent || isfield(options, 'autobrainshift')
            ea_compute_scrf(handles)
        end
end

% Choose default command line output for ea_subcorticalrefine
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ea_subcorticalrefine wait for user response (see UIRESUME)
if options.d2.write || options.d3.write
    uiwait(handles.scrf);
end

if isfield(options,'autobrainshift') && options.autobrainshift
    saveandclose([], [], handles); % close figure again.
end


% --- Outputs from this function are returned to the command line.
function varargout = ea_subcorticalrefine_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
%varargout{1} = handles.output;


% --- Executes on selection change in methodbutton.
function methodbutton_Callback(hObject, eventdata, handles)
% hObject    handle to methodbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns methodbutton contents as cell array
%        contents{get(hObject,'Value')} returns selected item from methodbutton


% --- Executes during object creation, after setting all properties.
function methodbutton_CreateFcn(hObject, eventdata, handles)
% hObject    handle to methodbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in computebutn.
function computebutn_Callback(hObject, eventdata, handles)
% hObject    handle to computebutn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ea_compute_scrf(handles)


% --- Executes on button press in approvebutn.
function approvebutn_Callback(hObject, eventdata, handles)
% hObject    handle to approvebutn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

options = getappdata(handles.scrf,'options');
if ~isfile(options.subj.brainshift.transform.instore)
	msgbox('Please generate a transform first (Click on "Compute subcortical refine transform"). If you don''t want to compute a transform, simply click on "Continue without subcortical transform".');
else
    copyfile(options.subj.brainshift.transform.converted, options.subj.brainshift.transform.scrf);
    if isfile(options.subj.recon.recon)
        ea_recalc_reco([],[],options.subj.subjDir);
    end

    % Set approval status
    if isfile(options.subj.brainshift.log.method)
        json = loadjson(options.subj.brainshift.log.method);
        json.approval = 1;
        savejson('', json, options.subj.brainshift.log.method);
    end

    ea_methods(options,...
        ['DBS electrode localizations were corrected for brainshift in postoperative acquisitions by applying a refined affine transform calculated between ',...
        'pre- and postoperative acquisitions that were restricted to a subcortical area of interest as implemented in the brainshift-correction module of Lead-DBS software',...
        ' (Horn & Kuehn 2005; SCR_002915; https://www.lead-dbs.org).'],...
        {'Horn, A., & Kuehn, A. A. (2015). Lead-DBS: a toolbox for deep brain stimulation electrode localizations and visualizations. NeuroImage, 107, 127?135. http://doi.org/10.1016/j.neuroimage.2014.12.002'});

    closescrf(handles);
end


function saveandclose(hObject, eventdata, handles)
% hObject    handle to approvebutn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

options=getappdata(handles.scrf,'options');
if ~isfile(options.subj.brainshift.transform.instore)
	msgbox('Please generate a transform first (Click on "Compute subcortical refine transform"). If you don''t want to compute a transform, simply click on "Continue without subcortical transform".');
else
    copyfile(options.subj.brainshift.transform.converted, options.subj.brainshift.transform.scrf);
    if isfile(options.subj.recon.recon)
        ea_recalc_reco([],[],options.subj.subjDir);
    end

    ea_methods(options,...
        ['DBS electrode localizations were corrected for brainshift in postoperative acquisitions by applying a refined affine transform calculated between ',...
        'pre- and postoperative acquisitions that were restricted to a subcortical area of interest as implemented in the brainshift-correction module of Lead-DBS software',...
        ' (Horn & Kuehn 2005; SCR_002915; https://www.lead-dbs.org).'],...
        {'Horn, A., & Kuehn, A. A. (2015). Lead-DBS: a toolbox for deep brain stimulation electrode localizations and visualizations. NeuroImage, 107, 127?135. http://doi.org/10.1016/j.neuroimage.2014.12.002'});

    closescrf(handles);
end

% --- Executes on button press in disapprovebutn.
function disapprovebutn_Callback(hObject, eventdata, handles)
% hObject    handle to disapprovebutn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

options=getappdata(handles.scrf,'options');
ea_delete(options.subj.brainshift.transform.converted);
ea_delete(options.subj.brainshift.transform.scrf);

% Set approval status
if isfile(options.subj.brainshift.log.method)
    json = loadjson(options.subj.brainshift.log.method);
    json.approval = 0;
    savejson('', json, options.subj.brainshift.log.method);
end

closescrf(handles);


function closescrf(handles)

options=getappdata(handles.scrf,'options');
if options.d2.write || options.d3.write
    uiresume(handles.scrf);
end
delete(handles.scrf);


% --- Executes on button press in mask0.
function mask0_Callback(hObject, eventdata, handles)
% hObject    handle to mask0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of mask0
handles.mask1.Value=0;
handles.mask2.Value=0;
% --- Executes on button press in mask1.


% --- Executes on button press in mask1.
function mask1_Callback(hObject, eventdata, handles)
% hObject    handle to mask1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.mask0.Value=0;
handles.mask2.Value=0;
% Hint: get(hObject,'Value') returns toggle state of mask1


% --- Executes on button press in mask2.
function mask2_Callback(hObject, eventdata, handles)
% hObject    handle to mask2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.mask0.Value=0;
handles.mask1.Value=0;
% Hint: get(hObject,'Value') returns toggle state of mask2


% --- Executes on button press in back.
function back_Callback(hObject, eventdata, handles)
% hObject    handle to back (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
options=getappdata(handles.scrf,'options');
ea_checkreg(options);
closescrf(handles);


% --- Executes on button press in openpatientdir.
function openpatientdir_Callback(hObject, eventdata, handles)
% hObject    handle to openpatientdir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
options = getappdata(handles.scrf,'options');
ea_opendir(options.subj.subjDir);
