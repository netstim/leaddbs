function varargout = ea_stimparams(varargin)
% EA_STIMPARAMS MATLAB code for ea_stimparams.fig
%      EA_STIMPARAMS, by itself, creates a new EA_STIMPARAMS or raises the existing
%      singleton*.
%
%      H = EA_STIMPARAMS returns the handle to a new EA_STIMPARAMS or the handle to
%      the existing singleton*.
%
%      EA_STIMPARAMS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EA_STIMPARAMS.M with the given input arguments.
%
%      EA_STIMPARAMS('Property','Value',...) creates a new EA_STIMPARAMS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ea_stimparams_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ea_stimparams_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ea_stimparams

% Last Modified by GUIDE v2.5 03-Mar-2021 16:26:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @ea_stimparams_OpeningFcn, ...
    'gui_OutputFcn',  @ea_stimparams_OutputFcn, ...
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


% --- Executes just before ea_stimparams is made visible.
function ea_stimparams_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ea_stimparams (see VARARGIN)

% Change name
set(handles.stimfig,'Name','Stimulation Parameters');

% store input arguments in figure to make it available to subroutines.
elstruct = varargin{1};
resultfig = varargin{2};
options = varargin{3};

setappdata(handles.stimfig,'elstruct',elstruct);

if options.native
    set(handles.estimateInTemplate,'Visible','off');
end

set(handles.estimateInTemplate,'Value',ea_getprefs('vatsettings.estimateInTemplate'));

if strcmp(options.leadprod, 'group')
    groupmode=1;
    M=getappdata(resultfig,'M');
    try % priorly loaded stim params.
        gS=M.S;
        setappdata(handles.stimfig,'gS',gS);
        gSv.vatmodel=M.vatmodel;
        setappdata(handles.stimfig,'gSv',gSv);
    end
    set(handles.prevpt,'visible','on');
    set(handles.nextpt,'visible','on');
    set(handles.saveparams,'visible','on');

    set(handles.stimlabel,'visible','off');
    set(handles.stimlab,'visible','off');
    set(handles.stimulate,'visible','off');
    if length(elstruct) == 1
        set(handles.prevpt,'enable','off');
        set(handles.nextpt,'enable','off');
    end

    % UI adaption
    Yshift = 24;
    set(handles.stimfig, 'Position', handles.stimfig.Position - [0 0 0 Yshift]);
    set(handles.headertxt, 'Position', handles.headertxt.Position - [0 Yshift 0 0]);
    set(handles.settings, 'Position', handles.settings.Position - [0 Yshift+2 0 0]);
else
    groupmode=0;

    set(handles.headertxt, 'String', ['Patient: ', elstruct(1).name]);

    set(handles.prevpt,'visible','off');
    set(handles.nextpt,'visible','off');
    set(handles.saveparams,'visible','off');

    set(handles.stimlabel,'visible','on');
    set(handles.stimlab,'visible','on');
    set(handles.stimulate,'visible','on');
end

setappdata(handles.stimfig,'groupmode',groupmode);
setappdata(handles.stimfig,'actpt',1);

setappdata(handles.stimfig,'elstruct',elstruct);
setappdata(handles.stimfig,'resultfig',resultfig);
setappdata(handles.stimfig,'options',options);

stimparams=getappdata(resultfig,'stimparams'); % get info from resultfig.
if isempty(stimparams)
    stimparams = struct();
end
setappdata(handles.stimfig,'stimparams',stimparams); % store stimulation settings from resultfig to stim (this) fig for subroutines.

% setup modelselect popup
if strcmp(options.leadprod, 'group')
    isdirected=0; % for now allow everything in lead group
else
    e=load(fullfile(ea_getearoot,'templates','electrode_models',options.elspec.matfname));
    if isfield(e.electrode,'isdirected')
        isdirected=e.electrode.isdirected;
    else
        isdirected=0;
    end
end

% Hide OSS-DBS option in case non-dev env or elmodel not available
if strcmp(options.leadprod, 'dbs')
    funcs = ea_regexpdir(ea_getearoot, 'ea_genvat_.*\.m$', 0);
    funcs = regexp(funcs, '(ea_genvat_.*)(?=\.m)', 'match', 'once');
    [names, supportDirected] = cellfun(@(x) eval([x, '(''prompt'');']), funcs, 'Uni', 0);
    if isdirected
        funcs = funcs(cell2mat(supportDirected));
        names = names(cell2mat(supportDirected));
    end
    if ~options.prefs.env.dev || ~ismember(options.elmodel,ea_ossdbs_elmodel)
        ossdbsInd = find(contains(names,'OSS-DBS'));
        funcs(ossdbsInd) = [];
        names(ossdbsInd) = [];
    end
else % Call in lead 'group'
    funcs = getappdata(resultfig, 'genvatfunctions');
    names = getappdata(resultfig, 'vatfunctionnames');
end

setappdata(gcf, 'genvatfunctions', funcs);
value = find(contains(names, handles.modelselect.String(handles.modelselect.Value)));
set(handles.modelselect, 'String', names);
set(handles.modelselect, 'Value', value);

% if ~isempty(stimparams) % stimfigure has been used before..
%     for side=1:2
%         for el=1:4
%             %keyboard
%             set(eval(['handles.k',num2str(((side-1)*4)+el-1),'u']),'String', num2str(stimparams(side).U(el)));
%             set(eval(['handles.k',num2str(((side-1)*4)+el-1),'im']),'String',num2str(stimparams(side).Im(el)));
%         end
%     end
%
%     set(handles.fiberthresh,'String',num2str(stimparams(1).fiberthresh))
%     set(handles.showfibs,'Value',stimparams(1).showfibers);
%     set(handles.showconns,'Value',stimparams(1).showconnectivities);
% end

pos=get(handles.stimfig,'position');
set(handles.stimfig,'position',[51,51,pos(3),pos(4)]);

ea_refreshguisp(handles,options);

if ~strcmp(options.leadprod, 'group')
    label =handles.stimlabel.String{handles.stimlabel.Value};
    label(strfind(label, ' ')) = '';
    stimDir = fullfile(options.subj.stimDir, ea_nt(options), label);
    filePrefix = ['sub-', options.subj.subjId, '_sim-'];
    visualizeVAT = 1;

    stimParams = ea_regexpdir(stimDir, 'stimparameters\.mat$', 0);
    load(stimParams{1}, 'S');
    modelLabel = ea_simModel2Label(S.model);

    if visualizeVAT
        if isfile([stimDir, filesep, filePrefix, 'binary_model-', modelLabel, '_hemi-R.mat']) && isfile([stimDir, filesep, filePrefix, 'binary_model-', modelLabel, '_hemi-L.mat'])
            load([stimDir, filesep, filePrefix, 'binary_model-', modelLabel, '_hemi-R.mat']);
            stimparams(1,1).VAT.VAT = vatfv;
            stimparams(1,1).volume = vatvolume;
            if exist('vatgrad','var')
                vatgradtemp(1) = vatgrad;
            end
            load([stimDir, filesep, filePrefix, 'binary_model-', modelLabel, '_hemi-L.mat']);
            stimparams(1,2).VAT.VAT = vatfv;
            stimparams(1,2).volume = vatvolume;
            if exist('vatgrad','var')
                vatgradtemp(2) = vatgrad;
                vatgrad = vatgradtemp;
            end
        elseif isfile([stimDir, filesep, filePrefix, 'binary_model-', modelLabel, '_hemi-R.mat'])
            load([stimDir, filesep, filePrefix, 'binary_model-', modelLabel, '_hemi-R.mat']);
            stimparams(1,1).VAT.VAT = vatfv;
            stimparams(1,1).volume = vatvolume;
        elseif isfile([stimDir, filesep, filePrefix, 'binary_model-', modelLabel, '_hemi-L.mat'])
            load([stimDir, filesep, filePrefix, 'binary_model-', modelLabel, '_hemi-L.mat']);
            %For consistency, left is always on 2nd element of stimparams
            stimparams(1,2).VAT.VAT = vatfv;
            stimparams(1,2).volume = vatvolume;
        else
            if isfile([stimDir, filesep, filePrefix, 'binary_model-', modelLabel, '_hemi-R.nii']) && isfile([stimDir, filesep, filePrefix, 'binary_model-', modelLabel, '_hemi-L.nii'])
                nii = ea_load_nii([stimDir, filesep, filePrefix, 'binary_model-', modelLabel, '_hemi-R.nii']);
                vatfv = ea_niiVAT2fvVAT(nii);
                stimparams(1,1).VAT.VAT = vatfv;
                stimparams(1,1).volume = sum(nii.img(:))*nii.voxsize(1)*nii.voxsize(2)*nii.voxsize(3);
                nii = ea_load_nii([stimDir, filesep, filePrefix, 'binary_model-', modelLabel, '_hemi-L.nii']);
                vatfv = ea_niiVAT2fvVAT(nii);
                stimparams(1,2).VAT.VAT = vatfv;
                stimparams(1,2).volume = sum(nii.img(:))*nii.voxsize(1)*nii.voxsize(2)*nii.voxsize(3);
            elseif isfile([stimDir, filesep, filePrefix, 'binary_model-', modelLabel, '_hemi-R.nii'])
                nii = ea_load_nii([stimDir, filesep, filePrefix, 'binary_model-', modelLabel, '_hemi-R.nii']);
                vatfv = ea_niiVAT2fvVAT(nii);
                stimparams(1,1).VAT.VAT = vatfv;
                stimparams(1,1).volume = sum(nii.img(:))*nii.voxsize(1)*nii.voxsize(2)*nii.voxsize(3);
            elseif isfile([stimDir, filesep, filePrefix, 'binary_model-', modelLabel, '_hemi-L.nii'])
                nii = ea_load_nii([stimDir, filesep, filePrefix, 'binary_model-', modelLabel, '_hemi-L.nii']);
                vatfv = ea_niiVAT2fvVAT(nii);
                %For consistency, left is always on 2nd element of stimparams
                stimparams(1,2).VAT.VAT = vatfv;
                stimparams(1,2).volume = sum(nii.img(:))*nii.voxsize(1)*nii.voxsize(2)*nii.voxsize(3);
            else
                visualizeVAT = 0;
            end
        end

        if isfile([stimDir, filesep, filePrefix, 'fiberactivation_model-', modelLabel, '_hemi-R.mat']) ...
                && isfile([stimDir, filesep, filePrefix, 'fiberactivation_model-', modelLabel, '_hemi-L.mat'])
            resultfig = getappdata(handles.stimfig,'resultfig');
            PL=getappdata(resultfig,'PL');
            for group=1:length(PL)
                ea_deletePL(PL(group));
            end
            clear PL
            ea_fiberactivation_viz([stimDir, filesep, filePrefix, 'fiberactivation_model-', modelLabel, '_hemi-R.mat'], resultfig);
            ea_fiberactivation_viz([stimDir, filesep, filePrefix, 'fiberactivation_model-', modelLabel, '_hemi-L.mat'], resultfig);
        elseif isfile([stimDir, filesep, filePrefix, 'fiberactivation_model-', modelLabel, '_hemi-R.mat'])
            resultfig = getappdata(handles.stimfig,'resultfig');
            PL=getappdata(resultfig,'PL');
            for group=1:length(PL)
                ea_deletePL(PL(group));
            end
            clear PL
            ea_fiberactivation_viz([stimDir, filesep, filePrefix, 'fiberactivation_model-', modelLabel, '_hemi-R.mat'], resultfig);
        elseif isfile([stimDir, filesep, filePrefix, 'fiberactivation_model-', modelLabel, '_hemi-L.mat'])
            resultfig = getappdata(handles.stimfig,'resultfig');
            PL=getappdata(resultfig,'PL');
            for group=1:length(PL)
                ea_deletePL(PL(group));
            end
            clear PL
            ea_fiberactivation_viz([stimDir, filesep, filePrefix, 'fiberactivation_model-', modelLabel, '_hemi-L.mat'], resultfig);
        end

        if visualizeVAT
            setappdata(handles.stimfig,'stimparams',stimparams);
            resultfig = getappdata(handles.stimfig,'resultfig');
            PL=getappdata(resultfig,'PL');
            for group=1:length(PL)
                ea_deletePL(PL(group));
            end
            clear PL
            if exist('vatgrad')
                setappdata(resultfig,'vatgrad',vatgrad);
            end
            setappdata(resultfig,'stimparams',stimparams(1,:));
            S=ea_loadstimulation(label,options);
            setappdata(resultfig,'curS',S(1))
            options.writeoutstats = 1;
            ea_calc_vatstats(resultfig,options);
        end
    end
end

% Choose default command line output for ea_stimparams
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ea_stimparams wait for user response (see UIRESUME)
% uiwait(handles.stimfig);


% --- Outputs from this function are returned to the command line.
function varargout = ea_stimparams_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function k0u_Callback(hObject, eventdata, handles)
% hObject    handle to k0u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k0u as text
%        str2double(get(hObject,'String')) returns contents of k0u as a double
S=getappdata(handles.stimfig,'S');
options=getappdata(handles.stimfig,'options');

eval(['S.Rs',num2str(S.active(1)),'.k0.perc=',num2str(get(hObject,'String')),';']);

setappdata(handles.stimfig,'S',S);
ea_refreshguisp(handles,options,hObject);


% --- Executes during object creation, after setting all properties.
function k0u_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k0u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function k1u_Callback(hObject, eventdata, handles)
% hObject    handle to k1u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k1u as text
%        str2double(get(hObject,'String')) returns contents of k1u as a double
S=getappdata(handles.stimfig,'S');
options=getappdata(handles.stimfig,'options');

eval(['S.Rs',num2str(S.active(1)),'.k1.perc=',num2str(get(hObject,'String')),';']);

setappdata(handles.stimfig,'S',S);
ea_refreshguisp(handles,options,hObject);


% --- Executes during object creation, after setting all properties.
function k1u_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k1u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function k2u_Callback(hObject, eventdata, handles)
% hObject    handle to k2u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k2u as text
%        str2double(get(hObject,'String')) returns contents of k2u as a double
S=getappdata(handles.stimfig,'S');
options=getappdata(handles.stimfig,'options');

eval(['S.Rs',num2str(S.active(1)),'.k2.perc=',num2str(get(hObject,'String')),';']);

setappdata(handles.stimfig,'S',S);
ea_refreshguisp(handles,options,hObject);


% --- Executes during object creation, after setting all properties.
function k2u_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k2u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function k3u_Callback(hObject, eventdata, handles)
% hObject    handle to k3u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k3u as text
%        str2double(get(hObject,'String')) returns contents of k3u as a double
S=getappdata(handles.stimfig,'S');
options=getappdata(handles.stimfig,'options');

eval(['S.Rs',num2str(S.active(1)),'.k3.perc=',num2str(get(hObject,'String')),';']);

setappdata(handles.stimfig,'S',S);
ea_refreshguisp(handles,options,hObject);


% --- Executes during object creation, after setting all properties.
function k3u_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k3u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function k4u_Callback(hObject, eventdata, handles)
% hObject    handle to k3u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k3u as text
%        str2double(get(hObject,'String')) returns contents of k3u as a double
S=getappdata(handles.stimfig,'S');
options=getappdata(handles.stimfig,'options');

eval(['S.Rs',num2str(S.active(1)),'.k4.perc=',num2str(get(hObject,'String')),';']);

setappdata(handles.stimfig,'S',S);
ea_refreshguisp(handles,options,hObject);


% --- Executes during object creation, after setting all properties.
function k4u_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k3u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function k5u_Callback(hObject, eventdata, handles)
% hObject    handle to k3u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k3u as text
%        str2double(get(hObject,'String')) returns contents of k3u as a double
S=getappdata(handles.stimfig,'S');
options=getappdata(handles.stimfig,'options');

eval(['S.Rs',num2str(S.active(1)),'.k5.perc=',num2str(get(hObject,'String')),';']);

setappdata(handles.stimfig,'S',S);
ea_refreshguisp(handles,options,hObject);


% --- Executes during object creation, after setting all properties.

function k5u_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k3u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function k6u_Callback(hObject, eventdata, handles)
% hObject    handle to k3u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k3u as text
%        str2double(get(hObject,'String')) returns contents of k3u as a double
S=getappdata(handles.stimfig,'S');
options=getappdata(handles.stimfig,'options');

eval(['S.Rs',num2str(S.active(1)),'.k6.perc=',num2str(get(hObject,'String')),';']);

setappdata(handles.stimfig,'S',S);
ea_refreshguisp(handles,options,hObject);


% --- Executes during object creation, after setting all properties.
function k6u_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k3u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function k7u_Callback(hObject, eventdata, handles)
% hObject    handle to k3u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k3u as text
%        str2double(get(hObject,'String')) returns contents of k3u as a double
S=getappdata(handles.stimfig,'S');
options=getappdata(handles.stimfig,'options');

eval(['S.Rs',num2str(S.active(1)),'.k7.perc=',num2str(get(hObject,'String')),';']);

setappdata(handles.stimfig,'S',S);
ea_refreshguisp(handles,options,hObject);


% --- Executes during object creation, after setting all properties.
function k7u_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k3u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function k4im_Callback(hObject, eventdata, handles)
% hObject    handle to k3u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k3u as text
%        str2double(get(hObject,'String')) returns contents of k3u as a double
S=getappdata(handles.stimfig,'S');
options=getappdata(handles.stimfig,'options');

eval(['S.Rs',num2str(S.active(1)),'.k4.imp=',num2str(get(hObject,'String')),';']);

setappdata(handles.stimfig,'S',S);
ea_refreshguisp(handles,options);


% --- Executes during object creation, after setting all properties.
function k4im_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k3u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function k5im_Callback(hObject, eventdata, handles)
% hObject    handle to k3u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k3u as text
%        str2double(get(hObject,'String')) returns contents of k3u as a double
S=getappdata(handles.stimfig,'S');
options=getappdata(handles.stimfig,'options');

eval(['S.Rs',num2str(S.active(1)),'.k5.imp=',num2str(get(hObject,'String')),';']);

setappdata(handles.stimfig,'S',S);
ea_refreshguisp(handles,options);

% --- Executes during object creation, after setting all properties.

function k5im_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k3u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function k6im_Callback(hObject, eventdata, handles)
% hObject    handle to k3u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k3u as text
%        str2double(get(hObject,'String')) returns contents of k3u as a double
S=getappdata(handles.stimfig,'S');
options=getappdata(handles.stimfig,'options');

eval(['S.Rs',num2str(S.active(1)),'.k6.imp=',num2str(get(hObject,'String')),';']);

setappdata(handles.stimfig,'S',S);
ea_refreshguisp(handles,options);


% --- Executes during object creation, after setting all properties.
function k6im_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k3u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function k7im_Callback(hObject, eventdata, handles)
% hObject    handle to k3u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k3u as text
%        str2double(get(hObject,'String')) returns contents of k3u as a double
S=getappdata(handles.stimfig,'S');
options=getappdata(handles.stimfig,'options');

eval(['S.Rs',num2str(S.active(1)),'.k7.imp=',num2str(get(hObject,'String')),';']);

setappdata(handles.stimfig,'S',S);
ea_refreshguisp(handles,options);


% --- Executes during object creation, after setting all properties.
function k7im_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k3u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function k0im_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k0im (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function k0im_Callback(hObject, eventdata, handles)
% hObject    handle to k1im (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k1im as text
%        str2double(get(hObject,'String')) returns contents of k1im as a double
S=getappdata(handles.stimfig,'S');
options=getappdata(handles.stimfig,'options');

eval(['S.Rs',num2str(S.active(1)),'.k0.imp=',num2str(get(hObject,'String')),';']);

setappdata(handles.stimfig,'S',S);
ea_refreshguisp(handles,options);


function k1im_Callback(hObject, eventdata, handles)
% hObject    handle to k1im (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k1im as text
%        str2double(get(hObject,'String')) returns contents of k1im as a double
S=getappdata(handles.stimfig,'S');
options=getappdata(handles.stimfig,'options');

eval(['S.Rs',num2str(S.active(1)),'.k1.imp=',num2str(get(hObject,'String')),';']);

setappdata(handles.stimfig,'S',S);
ea_refreshguisp(handles,options);


% --- Executes during object creation, after setting all properties.
function k1im_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k1im (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function k2im_Callback(hObject, eventdata, handles)
% hObject    handle to k2im (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k2im as text
%        str2double(get(hObject,'String')) returns contents of k2im as a double
S=getappdata(handles.stimfig,'S');
options=getappdata(handles.stimfig,'options');

eval(['S.Rs',num2str(S.active(1)),'.k2.imp=',num2str(get(hObject,'String')),';']);

setappdata(handles.stimfig,'S',S);
ea_refreshguisp(handles,options);


% --- Executes during object creation, after setting all properties.
function k2im_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k2im (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function k3im_Callback(hObject, eventdata, handles)
% hObject    handle to k3im (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k3im as text
%        str2double(get(hObject,'String')) returns contents of k3im as a double
S=getappdata(handles.stimfig,'S');
options=getappdata(handles.stimfig,'options');

eval(['S.Rs',num2str(S.active(1)),'.k3.imp=',num2str(get(hObject,'String')),';']);

setappdata(handles.stimfig,'S',S);
ea_refreshguisp(handles,options);


% --- Executes during object creation, after setting all properties.
function k3im_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k3im (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function k8u_Callback(hObject, eventdata, handles)
% hObject    handle to k8u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k8u as text
%        str2double(get(hObject,'String')) returns contents of k8u as a double
S=getappdata(handles.stimfig,'S');
options=getappdata(handles.stimfig,'options');

eval(['S.Ls',num2str(S.active(2)),'.k8.perc=',num2str(get(hObject,'String')),';']);

setappdata(handles.stimfig,'S',S);
ea_refreshguisp(handles,options,hObject);


% --- Executes during object creation, after setting all properties.
function k8u_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k8u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function k9u_Callback(hObject, eventdata, handles)
% hObject    handle to k9u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k9u as text
%        str2double(get(hObject,'String')) returns contents of k9u as a double
S=getappdata(handles.stimfig,'S');
options=getappdata(handles.stimfig,'options');

eval(['S.Ls',num2str(S.active(2)),'.k9.perc=',num2str(get(hObject,'String')),';']);

setappdata(handles.stimfig,'S',S);
ea_refreshguisp(handles,options,hObject);


% --- Executes during object creation, after setting all properties.
function k9u_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k9u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function k10u_Callback(hObject, eventdata, handles)
% hObject    handle to k10u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k10u as text
%        str2double(get(hObject,'String')) returns contents of k10u as a double
S=getappdata(handles.stimfig,'S');
options=getappdata(handles.stimfig,'options');

eval(['S.Ls',num2str(S.active(2)),'.k10.perc=',num2str(get(hObject,'String')),';']);

setappdata(handles.stimfig,'S',S);
ea_refreshguisp(handles,options,hObject);


% --- Executes during object creation, after setting all properties.
function k10u_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k10u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function k11u_Callback(hObject, eventdata, handles)
% hObject    handle to k11u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k11u as text
%        str2double(get(hObject,'String')) returns contents of k11u as a double
S=getappdata(handles.stimfig,'S');
options=getappdata(handles.stimfig,'options');

eval(['S.Ls',num2str(S.active(2)),'.k11.perc=',num2str(get(hObject,'String')),';']);

setappdata(handles.stimfig,'S',S);
ea_refreshguisp(handles,options,hObject);


% --- Executes during object creation, after setting all properties.
function k11u_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k11u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function k8im_Callback(hObject, eventdata, handles)
% hObject    handle to k8im (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k8im as text
%        str2double(get(hObject,'String')) returns contents of k8im as a double
S=getappdata(handles.stimfig,'S');
options=getappdata(handles.stimfig,'options');

eval(['S.Ls',num2str(S.active(2)),'.k8.imp=',num2str(get(hObject,'String')),';']);

setappdata(handles.stimfig,'S',S);
ea_refreshguisp(handles,options);


% --- Executes during object creation, after setting all properties.
function k8im_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k8im (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function k9im_Callback(hObject, eventdata, handles)
% hObject    handle to k9im (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k9im as text
%        str2double(get(hObject,'String')) returns contents of k9im as a double
S=getappdata(handles.stimfig,'S');
options=getappdata(handles.stimfig,'options');

eval(['S.Ls',num2str(S.active(2)),'.k9.imp=',num2str(get(hObject,'String')),';']);

setappdata(handles.stimfig,'S',S);
ea_refreshguisp(handles,options);


% --- Executes during object creation, after setting all properties.
function k9im_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k9im (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function k10im_Callback(hObject, eventdata, handles)
% hObject    handle to k10im (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k10im as text
%        str2double(get(hObject,'String')) returns contents of k10im as a double
S=getappdata(handles.stimfig,'S');
options=getappdata(handles.stimfig,'options');

eval(['S.Ls',num2str(S.active(2)),'.k10.imp=',num2str(get(hObject,'String')),';']);

setappdata(handles.stimfig,'S',S);
ea_refreshguisp(handles,options);


% --- Executes during object creation, after setting all properties.
function k10im_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k10im (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function k11im_Callback(hObject, eventdata, handles)
% hObject    handle to k11im (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k11im as text
%        str2double(get(hObject,'String')) returns contents of k11im as a double
S=getappdata(handles.stimfig,'S');
options=getappdata(handles.stimfig,'options');

eval(['S.Ls',num2str(S.active(2)),'.k11.imp=',num2str(get(hObject,'String')),';']);

setappdata(handles.stimfig,'S',S);
ea_refreshguisp(handles,options);


% --- Executes during object creation, after setting all properties.
function k11im_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k11im (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in modelselect.
function modelselect_Callback(hObject, eventdata, handles)
% hObject    handle to modelselect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns modelselect contents as cell array
%        contents{get(hObject,'Value')} returns selected item from modelselect

% Set model
S = getappdata(handles.stimfig,'S');
models = get(handles.modelselect,'String');
S.model = models{get(handles.modelselect,'Value')};
setappdata(handles.stimfig,'S',S);

% Handle special case for call from LeadGroup
groupmode=getappdata(handles.stimfig,'groupmode');
if groupmode
    choice = questdlg('Changing VAT model will delete stimulation parameters of all patients! Continue?', ...
        'Warning', ...
        'Yes, sure','No','No');

    gSv=getappdata(handles.stimfig,'gSv');

    switch choice
        case 'No'
            % Keep current model
            currentModelInd = find(ismember(get(hObject,'String'),gSv.vatmodel),1);
            set(hObject,'Value',currentModelInd);
            return
        case 'Yes, sure'
            setappdata(handles.stimfig,'gS',[]);

            % Set new model
            models = get(hObject,'String');
            gSv.vatmodel = models{get(hObject,'Value')};
            setappdata(handles.stimfig,'gSv',gSv);

            % Clear stimulation parameters
            setappdata(handles.stimfig,'S',[]);
    end
end

% Refresh GUI
options = getappdata(handles.stimfig,'options');
ea_refreshguisp(handles,options);

% Save stimulation parameters
if ~groupmode
    S = getappdata(handles.stimfig,'S');
    ea_savestimulation(S,options);
end


% --- Executes during object creation, after setting all properties.
function modelselect_CreateFcn(hObject, eventdata, handles)
% hObject    handle to modelselect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in stimulate.
function stimulate_Callback(hObject, eventdata, handles)
% hObject    handle to stimulate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ea_busyaction('on',handles.stimfig,'stim');
elstruct = getappdata(handles.stimfig,'elstruct');
resultfig = getappdata(handles.stimfig,'resultfig');
options = getappdata(handles.stimfig,'options');
% refresh prefs:
options.prefs = ea_prefs;
setappdata(resultfig,'options',options);
setappdata(handles.stimfig,'options',options);
S = getappdata(handles.stimfig,'S');
S = ea_activecontacts(S);

options = getappdata(resultfig,'options'); % selected atlas could have refreshed.
options.orignative = options.native;

if strcmp('on',get(handles.estimateInTemplate,'Visible')) % only allowed for specific VTA functions
    switch get(handles.estimateInTemplate,'Value')
        case 0
            S.template = 'warp';
            options.native = 1;
        case 1
            S.template = 'direct';
    end
end

if ~isfield(options.subj, 'norm') && options.native
    ea_cprintf('CmdWinWarnings', 'Calculating VTA in template space since patient folder %s is incomplete.\n', options.subj.subjId);
    options.native = 0;
end

ea_savestimulation(S,options);
setappdata(handles.stimfig,'S',S);

if isfield(elstruct,'group') % group analysis, more than one electrode set
    % this should not happen, in this case the stim button is
    % hidden.
    keyboard
end

% assign correct .m-file to function.
genvatfunctions = getappdata(handles.stimfig,'genvatfunctions');
ea_genvat = eval(['@',genvatfunctions{get(handles.modelselect,'Value')}]);
stimname = S.label;

for el = 1:length(elstruct)
    % Load stim coordinates
    if options.native % Reload native space coordinates
        coords = ea_load_reconstruction(options);
    else
        coords = elstruct(el).coords_mm;
    end

    if strcmp(S.model, 'OSS-DBS (Butenko 2020)') % For OSS-DBS, side iteration is within the genvat function
        % Set stimSetMode flag to options
        % (avoid additional parameter or setting appdata, to make it scriptable)
        if handles.addStimSet.Value
            options.stimSetMode = 1;
        else
            options.stimSetMode = 0;
        end
        if options.prefs.machine.vatsettings.butenko_calcAxonActivation
            feval(ea_genvat,getappdata(handles.stimfig,'S'),options,handles.stimfig);
            ea_busyaction('off',handles.stimfig,'stim');
            return;
        else
            [~, stimparams] = feval(ea_genvat,getappdata(handles.stimfig,'S'),options,handles.stimfig);
            flix = 1;
        end
    else
        stimparams = struct();
        for iside = 1:length(options.sides)
            side = options.sides(iside);
            [vatfv, vatvolume] = feval(ea_genvat,coords,getappdata(handles.stimfig,'S'),side,options,stimname,handles.stimfig);
            stimparams(1,side).VAT(el).VAT = vatfv;
            stimparams(1,side).volume = vatvolume;
            flix = 1;
        end
    end
end

options.native = options.orignative;
PL = getappdata(resultfig, 'PL');
for group = 1:length(PL)
    ea_deletePL(PL(group));
end
clear PL

for group = flix
    setappdata(resultfig,'stimparams',stimparams(group,:));
    setappdata(resultfig,'curS',S(group));

    if ~exist('hmchanged','var')
        hmchanged = 1;
    end
    ea_calc_vatstats(resultfig,options,hmchanged);

    try % TODO: fix dir
        copyfile([options.root,options.patientname,filesep,'ea_pm.nii'],[options.root,options.patientname,filesep,'ea_pm_group_',num2str(group),'.nii']);
    end

    try
        PL(group) = getappdata(resultfig,'PL');
    catch
        keyboard
    end
end
setappdata(resultfig,'PL',PL);

ea_busyaction('off',handles.stimfig,'stim');





function k12u_Callback(hObject, eventdata, handles)
% hObject    handle to k12u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k12u as text
%        str2double(get(hObject,'String')) returns contents of k12u as a double
S=getappdata(handles.stimfig,'S');
options=getappdata(handles.stimfig,'options');

eval(['S.Ls',num2str(S.active(2)),'.k12.perc=',num2str(get(hObject,'String')),';']);

setappdata(handles.stimfig,'S',S);
ea_refreshguisp(handles,options,hObject);


% --- Executes during object creation, after setting all properties.
function k12u_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k12u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function k13u_Callback(hObject, eventdata, handles)
% hObject    handle to k13u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k13u as text
%        str2double(get(hObject,'String')) returns contents of k13u as a double
S=getappdata(handles.stimfig,'S');
options=getappdata(handles.stimfig,'options');

eval(['S.Ls',num2str(S.active(2)),'.k13.perc=',num2str(get(hObject,'String')),';']);

setappdata(handles.stimfig,'S',S);
ea_refreshguisp(handles,options,hObject);


% --- Executes during object creation, after setting all properties.
function k13u_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k13u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function k14u_Callback(hObject, eventdata, handles)
% hObject    handle to k14u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k14u as text
%        str2double(get(hObject,'String')) returns contents of k14u as a double
S=getappdata(handles.stimfig,'S');
options=getappdata(handles.stimfig,'options');

eval(['S.Ls',num2str(S.active(2)),'.k14.perc=',num2str(get(hObject,'String')),';']);

setappdata(handles.stimfig,'S',S);
ea_refreshguisp(handles,options,hObject);


% --- Executes during object creation, after setting all properties.
function k14u_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k14u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function k15u_Callback(hObject, eventdata, handles)
% hObject    handle to k15u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k15u as text
%        str2double(get(hObject,'String')) returns contents of k15u as a double
S=getappdata(handles.stimfig,'S');
options=getappdata(handles.stimfig,'options');

eval(['S.Ls',num2str(S.active(2)),'.k15.perc=',num2str(get(hObject,'String')),';']);

setappdata(handles.stimfig,'S',S);
ea_refreshguisp(handles,options,hObject);


% --- Executes during object creation, after setting all properties.
function k15u_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k15u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function k12im_Callback(hObject, eventdata, handles)
% hObject    handle to k12im (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k12im as text
%        str2double(get(hObject,'String')) returns contents of k12im as a double
S=getappdata(handles.stimfig,'S');
options=getappdata(handles.stimfig,'options');

eval(['S.Ls',num2str(S.active(2)),'.k12.imp=',num2str(get(hObject,'String')),';']);

setappdata(handles.stimfig,'S',S);
ea_refreshguisp(handles,options);


% --- Executes during object creation, after setting all properties.
function k12im_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k12im (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function k13im_Callback(hObject, eventdata, handles)
% hObject    handle to k13im (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k13im as text
%        str2double(get(hObject,'String')) returns contents of k13im as a double
S=getappdata(handles.stimfig,'S');
options=getappdata(handles.stimfig,'options');

eval(['S.Ls',num2str(S.active(2)),'.k13.imp=',num2str(get(hObject,'String')),';']);

setappdata(handles.stimfig,'S',S);
ea_refreshguisp(handles,options);


% --- Executes during object creation, after setting all properties.
function k13im_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k13im (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function k14im_Callback(hObject, eventdata, handles)
% hObject    handle to k14im (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k14im as text
%        str2double(get(hObject,'String')) returns contents of k14im as a double
S=getappdata(handles.stimfig,'S');
options=getappdata(handles.stimfig,'options');

eval(['S.Ls',num2str(S.active(2)),'.k14.imp=',num2str(get(hObject,'String')),';']);

setappdata(handles.stimfig,'S',S);
ea_refreshguisp(handles,options);


% --- Executes during object creation, after setting all properties.
function k14im_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k14im (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function k15im_Callback(hObject, eventdata, handles)
% hObject    handle to k15im (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k15im as text
%        str2double(get(hObject,'String')) returns contents of k15im as a double
S=getappdata(handles.stimfig,'S');
options=getappdata(handles.stimfig,'options');

eval(['S.Ls',num2str(S.active(2)),'.k15.imp=',num2str(get(hObject,'String')),';']);

setappdata(handles.stimfig,'S',S);
ea_refreshguisp(handles,options);


% --- Executes during object creation, after setting all properties.
function k15im_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k15im (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function RCu_Callback(hObject, eventdata, handles)
% hObject    handle to RCu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RCu as text
%        str2double(get(hObject,'String')) returns contents of RCu as a double
S=getappdata(handles.stimfig,'S');
options=getappdata(handles.stimfig,'options');

eval(['S.Rs',num2str(S.active(1)),'.case.perc=',num2str(get(hObject,'String')),';']);

setappdata(handles.stimfig,'S',S);
ea_refreshguisp(handles,options,hObject);


% --- Executes during object creation, after setting all properties.
function RCu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RCu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Rs1am_Callback(hObject, eventdata, handles)
% hObject    handle to Rs1am (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Rs1am as text
%        str2double(get(hObject,'String')) returns contents of Rs1am as a double
S=getappdata(handles.stimfig,'S');
options=getappdata(handles.stimfig,'options');
S.active(1)=1;
S.Rs1.amp=str2double(get(hObject,'String'));

setappdata(handles.stimfig,'S',S);
ea_refreshguisp(handles,options);


% --- Executes during object creation, after setting all properties.
function Rs1am_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Rs1am (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Rs1va.
function Rs1va_Callback(hObject, eventdata, handles)
% hObject    handle to Rs1va (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Rs1va contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Rs1va
S=getappdata(handles.stimfig,'S');
options=getappdata(handles.stimfig,'options');
S.active(1)=1;
S.Rs1.va=get(hObject,'Value');

setappdata(handles.stimfig,'S',S);
ea_refreshguisp(handles,options);


% --- Executes during object creation, after setting all properties.
function Rs1va_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Rs1va (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Rs2am_Callback(hObject, eventdata, handles)
% hObject    handle to Rs2am (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Rs2am as text
%        str2double(get(hObject,'String')) returns contents of Rs2am as a double
S=getappdata(handles.stimfig,'S');
options=getappdata(handles.stimfig,'options');
S.active(1)=2;
S.Rs2.amp=str2double(get(hObject,'String'));

setappdata(handles.stimfig,'S',S);
ea_refreshguisp(handles,options);


% --- Executes during object creation, after setting all properties.
function Rs2am_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Rs2am (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Rs2va.
function Rs2va_Callback(hObject, eventdata, handles)
% hObject    handle to Rs2va (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Rs2va contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Rs2va
S=getappdata(handles.stimfig,'S');
options=getappdata(handles.stimfig,'options');
S.active(1)=2;
S.Rs2.va=get(hObject,'Value');
setappdata(handles.stimfig,'S',S);
ea_refreshguisp(handles,options);


% --- Executes during object creation, after setting all properties.
function Rs2va_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Rs2va (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function LCu_Callback(hObject, eventdata, handles)
% hObject    handle to LCu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LCu as text
%        str2double(get(hObject,'String')) returns contents of LCu as a double
S=getappdata(handles.stimfig,'S');
options=getappdata(handles.stimfig,'options');

eval(['S.Ls',num2str(S.active(2)),'.case.perc=',num2str(get(hObject,'String')),';']);

setappdata(handles.stimfig,'S',S);
ea_refreshguisp(handles,options,hObject);


% --- Executes during object creation, after setting all properties.
function LCu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LCu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Ls1am_Callback(hObject, eventdata, handles)
% hObject    handle to Ls1am (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ls1am as text
%        str2double(get(hObject,'String')) returns contents of Ls1am as a double
S=getappdata(handles.stimfig,'S');
options=getappdata(handles.stimfig,'options');
S.active(2)=1;
S.Ls1.amp=str2double(get(hObject,'String'));

setappdata(handles.stimfig,'S',S);
ea_refreshguisp(handles,options);


% --- Executes during object creation, after setting all properties.
function Ls1am_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ls1am (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Ls1va.
function Ls1va_Callback(hObject, eventdata, handles)
% hObject    handle to Ls1va (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Ls1va contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Ls1va
S=getappdata(handles.stimfig,'S');
options=getappdata(handles.stimfig,'options');
S.active(2)=1;
S.Ls1.va=get(hObject,'Value');

setappdata(handles.stimfig,'S',S);
ea_refreshguisp(handles,options);


% --- Executes during object creation, after setting all properties.
function Ls1va_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ls1va (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Ls2am_Callback(hObject, eventdata, handles)
% hObject    handle to Ls2am (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ls2am as text
%        str2double(get(hObject,'String')) returns contents of Ls2am as a double
S=getappdata(handles.stimfig,'S');
options=getappdata(handles.stimfig,'options');
S.active(2)=2;
S.Ls2.amp=str2double(get(hObject,'String'));

setappdata(handles.stimfig,'S',S);
ea_refreshguisp(handles,options);


% --- Executes during object creation, after setting all properties.
function Ls2am_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ls2am (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Ls2va.
function Ls2va_Callback(hObject, eventdata, handles)
% hObject    handle to Ls2va (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Ls2va contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Ls2va
S=getappdata(handles.stimfig,'S');
options=getappdata(handles.stimfig,'options');
S.active(2)=2;
S.Ls2.va=get(hObject,'Value');

setappdata(handles.stimfig,'S',S);
ea_refreshguisp(handles,options);


% --- Executes during object creation, after setting all properties.
function Ls2va_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ls2va (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Linterleaved.
function Linterleaved_Callback(hObject, eventdata, handles)
% hObject    handle to Linterleaved (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Linterleaved


function stimlabel_Callback(hObject, eventdata, handles)
% hObject    handle to stimlabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stimlabel as text
%        str2double(get(hObject,'String')) returns contents of stimlabel as a double
S = getappdata(handles.stimfig,'S');
options = getappdata(handles.stimfig,'options');
sel = get(handles.stimlabel,'String');
sel = sel{get(handles.stimlabel,'Value')};
if startsWith(sel,' => ') % command, not entry
    switch sel(5:end)
        case 'New stimulation'
            resultfig = getappdata(handles.stimfig,'resultfig');

            PL = getappdata(resultfig, 'PL');
            for group = 1:length(PL)
                ea_deletePL(PL(group));
            end
            clear PL

            ea_savestimulation(S,options);
            S = []; % this will create the prompt to generate a new S.
            options.gen_newstim = 1;
            setappdata(handles.stimfig,'options',options);
            setappdata(handles.stimfig,'S',S);
            ea_refreshguisp(handles,options);
            S = getappdata(handles.stimfig,'S');
            ea_savestimulation(S,options);
            options.gen_newstim = 0; % reset new stim flag
            setappdata(handles.stimfig,'options',options);
        case 'Rename stimulation'
            stimlabel = getappdata(handles.stimfig,'stimlabel');

            [~,ix] = ismember(stimlabel,get(handles.stimlabel,'String'));
            set(handles.stimlabel,'Value',ix);
            stimc = inputdlg('Please enter a label for this stimulation','Stimulation Label',1,{stimlabel});
            if isfolder([options.subj.stimDir,filesep,ea_nt(0),stimlabel])
                movefile([options.subj.stimDir,filesep,ea_nt(0),stimlabel],[options.subj.stimDir,filesep,ea_nt(0),stimc{1}]);
            end
            if isfolder([options.subj.stimDir,filesep,ea_nt(1),stimlabel])
                movefile([options.subj.stimDir,filesep,ea_nt(1),stimlabel],[options.subj.stimDir,filesep,ea_nt(1),stimc{1}]);
            end
            slabelc = get(handles.stimlabel,'String');
            slabelc{ix} = stimc{1};
            set(handles.stimlabel,'String',slabelc);
            S.label = stimc{1};
            setappdata(handles.stimfig,'S',S);
            setappdata(handles.stimfig,'stimlabel',S.label);
            ea_refreshguisp(handles,options);
            ea_savestimulation(S,options);
        case 'Delete stimulation'
            answ = questdlg(['Are you sure you wish to delete the stimulation parameters for ',...
                S.label,'?'],'Delete stimulation parameters','Sure','No','No');
            if strcmp(answ,'No')
                set(handles.stimlabel,'Value',1);
            else % truly delete Stimulation parameters
                resultfig = getappdata(handles.stimfig,'resultfig');
                PL = getappdata(resultfig, 'PL');
                for group = 1:length(PL)
                    ea_deletePL(PL(group));
                end
                clear PL
                ea_delete([options.subj.stimDir,filesep,ea_nt(0),S.label]);
                ea_delete([options.subj.stimDir,filesep,ea_nt(1),S.label]);
                S = []; % this will create the prompt to generate a new S.
                setappdata(handles.stimfig,'S',S);
                set(handles.stimlabel,'Value',1);
                setappdata(handles.stimfig,'stimlabel','');
                options.gen_newstim = 0;
                setappdata(handles.stimfig,'S',S);
                ea_refreshguisp(handles,options);
                S = getappdata(handles.stimfig,'S');
                ea_savestimulation(S,options);
            end
    end
else
    label = handles.stimlabel.String{handles.stimlabel.Value};
    label(strfind(label,' ')) = '';
    S = ea_loadstimulation(label,options);
    S.label = label;
    setappdata(handles.stimfig,'S',S);
    setappdata(handles.stimfig,'stimlabel',S.label);
    setappdata(handles.stimfig,'S',S);
    ea_refreshguisp(handles,options);

    %% stuff by Till for visualizing VATs by selecting them from the stimlabel list
    % tries to load .mat-files which are now created by ea_genvat_horn.m
    % and contain the VAT as well as the quiver. In case no .mat-files are
    % available the vat_xxx.nii is loaded and visualized

    stimDir = fullfile(options.subj.stimDir, ea_nt(options), label);
    filePrefix = ['sub-', options.subj.subjId, '_sim-'];

    stimParams = ea_regexpdir(stimDir, 'stimparameters\.mat$', 0);
    load(stimParams{1}, 'S');
    modelLabel = ea_simModel2Label(S.model);

    visualizeVAT = 1;
    stimparams = struct();
    if isfile([stimDir, filesep, filePrefix, 'binary_model-', modelLabel, '_hemi-R.mat']) && isfile([stimDir, filesep, filePrefix, 'binary_model-', modelLabel, '_hemi-L.mat'])
        load([stimDir, filesep, filePrefix, 'binary_model-', modelLabel, '_hemi-R.mat']);
        stimparams(1,1).VAT.VAT = vatfv;
        stimparams(1,1).volume = vatvolume;
        if exist('vatgrad','var')
            vatgradtemp(1) = vatgrad;
        end
        load([stimDir, filesep, filePrefix, 'binary_model-', modelLabel, '_hemi-L.mat']);
        stimparams(1,2).VAT.VAT = vatfv;
        stimparams(1,2).volume = vatvolume;
        if exist('vatgrad','var')
            vatgradtemp(2) = vatgrad;
            vatgrad = vatgradtemp;
        end
    elseif isfile([stimDir, filesep, filePrefix, 'binary_model-', modelLabel, '_hemi-R.mat'])
        load([stimDir, filesep, filePrefix, 'binary_model-', modelLabel, '_hemi-R.mat']);
        stimparams(1,1).VAT.VAT = vatfv;
        stimparams(1,1).volume = vatvolume;
    elseif isfile([stimDir, filesep, filePrefix, 'binary_model-', modelLabel, '_hemi-L.mat'])
        load([stimDir, filesep, filePrefix, 'binary_model-', modelLabel, '_hemi-L.mat']);
        stimparams(1,2).VAT.VAT = vatfv;
        stimparams(1,2).volume = vatvolume;
    else
        if isfile([stimDir, filesep, filePrefix, 'binary_model-', modelLabel, '_hemi-R.nii']) && isfile([stimDir, filesep, filePrefix, 'binary_model-', modelLabel, '_hemi-L.nii'])
            nii = ea_load_nii([stimDir, filesep, filePrefix, 'binary_model-', modelLabel, '_hemi-R.nii']);
            vatfv = ea_niiVAT2fvVAT(nii);
            stimparams(1,1).VAT.VAT = vatfv;
            stimparams(1,1).volume = sum(nii.img(:))*nii.voxsize(1)*nii.voxsize(2)*nii.voxsize(3);
            nii = ea_load_nii([stimDir, filesep, filePrefix, 'binary_model-', modelLabel, '_hemi-L.nii']);
            vatfv = ea_niiVAT2fvVAT(nii);
            stimparams(1,2).VAT.VAT = vatfv;
            stimparams(1,2).volume = sum(nii.img(:))*nii.voxsize(1)*nii.voxsize(2)*nii.voxsize(3);
        elseif isfile([stimDir, filesep, filePrefix, 'binary_model-', modelLabel, '_hemi-R.nii'])
            nii = ea_load_nii([stimDir, filesep, filePrefix, 'binary_model-', modelLabel, '_hemi-R.nii']);
            vatfv = ea_niiVAT2fvVAT(nii);
            stimparams(1,1).VAT.VAT = vatfv;
            stimparams(1,1).volume = sum(nii.img(:))*nii.voxsize(1)*nii.voxsize(2)*nii.voxsize(3);
        elseif isfile([stimDir, filesep, filePrefix, 'binary_model-', modelLabel, '_hemi-L.nii'])
            nii = ea_load_nii([stimDir, filesep, filePrefix, 'binary_model-', modelLabel, '_hemi-L.nii']);
            vatfv = ea_niiVAT2fvVAT(nii);
            stimparams(1,2).VAT.VAT = vatfv;
            stimparams(1,2).volume = sum(nii.img(:))*nii.voxsize(1)*nii.voxsize(2)*nii.voxsize(3);
        else
            visualizeVAT = 0;
        end
    end

    visualizeFiberActivation = 1;
    if isfile([stimDir, filesep, filePrefix, 'fiberactivation_model-', modelLabel, '_hemi-R.mat']) ...
            && isfile([stimDir, filesep, filePrefix, 'fiberactivation_model-', modelLabel, '_hemi-L.mat'])
        resultfig = getappdata(handles.stimfig,'resultfig');
        PL=getappdata(resultfig,'PL');
        for group=1:length(PL)
            ea_deletePL(PL(group));
        end
        clear PL
        ea_fiberactivation_viz([stimDir, filesep, filePrefix, 'fiberactivation_model-', modelLabel, '_hemi-R.mat'], resultfig);
        ea_fiberactivation_viz([stimDir, filesep, filePrefix, 'fiberactivation_model-', modelLabel, '_hemi-L.mat'], resultfig);
    elseif isfile([stimDir, filesep, filePrefix, 'fiberactivation_model-', modelLabel, '_hemi-R.mat'])
        resultfig = getappdata(handles.stimfig,'resultfig');
        PL=getappdata(resultfig,'PL');
        for group=1:length(PL)
            ea_deletePL(PL(group));
        end
        clear PL
        ea_fiberactivation_viz([stimDir, filesep, filePrefix, 'fiberactivation_model-', modelLabel, '_hemi-R.mat'], resultfig);
    elseif isfile([stimDir, filesep, filePrefix, 'fiberactivation_model-', modelLabel, '_hemi-L.mat'])
        resultfig = getappdata(handles.stimfig,'resultfig');
        PL=getappdata(resultfig,'PL');
        for group=1:length(PL)
            ea_deletePL(PL(group));
        end
        clear PL
        ea_fiberactivation_viz([stimDir, filesep, filePrefix, 'fiberactivation_model-', modelLabel, '_hemi-L.mat'], resultfig);
    else
        visualizeFiberActivation = 0;
    end

    if visualizeVAT
        setappdata(handles.stimfig,'stimparams',stimparams);
        resultfig = getappdata(handles.stimfig,'resultfig');
        PL=getappdata(resultfig,'PL');
        for group=1:length(PL)
            ea_deletePL(PL(group));
        end
        clear PL
        if exist('vatgrad')
            setappdata(resultfig,'vatgrad',vatgrad);
        end
        setappdata(resultfig,'stimparams',stimparams(1,:));
        setappdata(resultfig,'curS',S(1))
        options.writeoutstats = 1;
        ea_calc_vatstats(resultfig,options);
    end

    if ~visualizeVAT && ~visualizeFiberActivation
        fprintf('\n');
        warning('off', 'backtrace');
        warning('Nothing to be visualized, please rerun stimulation!!');
        warning('on', 'backtrace');
        fprintf('\n');
    end
end


% --- Executes during object creation, after setting all properties.
function stimlabel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stimlabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function ea_refreshguisp(varargin)

handles=varargin{1};
options=varargin{2};

elstruct=getappdata(handles.stimfig,'elstruct');
S=getappdata(handles.stimfig,'S');
groupmode=getappdata(handles.stimfig,'groupmode');
actpt=getappdata(handles.stimfig,'actpt');

if isempty(actpt) || length(actpt)>1
    actpt=1;
end

if groupmode
    grouploaded=getappdata(handles.stimfig,'grouploaded');
    if isempty(grouploaded) % this is done only once and gets the selection info from lead_group initially (which patient shown).
        lgfig=getappdata(handles.stimfig,'resultfig');
        M=getappdata(lgfig,'M');
        actpt=M.ui.listselect;

        if length(actpt)>1 % more than one entry selected
            actpt=1;
        end

        % Ensure active patient is non empty
        % This can happen if you delete a patient, then add a new one, without clicking on the patient window
        if isempty(actpt)
            actpt=1;
        end
        setappdata(handles.stimfig,'actpt',actpt);
        % set grouploaded true is being done below.
    end

    elstruct=getappdata(handles.stimfig,'elstruct');

    set(handles.headertxt,'String',['Patient (',num2str(actpt),'/', num2str(length(elstruct)),'): ',elstruct(actpt).name]);

    gSv=getappdata(handles.stimfig,'gSv');
    if isfield(gSv,'vatmodel')
        if isempty(gSv.vatmodel)
            nms=get(handles.modelselect,'String');
            try
                gSv.vatmodel=nms{get(handles.modelselect,'Value')};
            catch
                keyboard
            end
            setappdata(handles.stimfig,'gSv',gSv);
        else
            [~,ind]=ismember(gSv.vatmodel,get(handles.modelselect,'String'));
            set(handles.modelselect,'Value',ind);
        end
    else
        nms=get(handles.modelselect,'String');
        try
            gSv.vatmodel=nms{get(handles.modelselect,'Value')};
        catch
            keyboard
        end
        setappdata(handles.stimfig,'gSv',gSv);
    end

    % load gS - updated with each refresh:
    gS=getappdata(handles.stimfig,'gS');
    if isempty(grouploaded)
        if ~isempty(gS)
            % determine stimlabel from priorly set gS:
            for sub=1:length(gS)
                stimlabel=['gs_',M.guid];
                if ~isempty(stimlabel)
                    break
                end
            end
            setappdata(handles.stimfig,'stimlabel',stimlabel);

            % if gS is defined but group has just now been loaded
            try
                if ~isempty(gS(actpt).Rs1) % current patient is defined -> set S to gS of this patient.
                    S = gS(actpt);
                end
            end
        end

        % now tell everyone that the figure has been opened for a while already:
        setappdata(handles.stimfig,'grouploaded',1);
    end
end

stimlabel=getappdata(handles.stimfig,'stimlabel');

if isempty(S)
    S=ea_initializeS(stimlabel,options,handles);
    try S.model=gSv.vatmodel; end
    setappdata(handles.stimfig,'stimlabel',S.label);
else
    if isempty(S.Rs1)
        S=ea_initializeS(stimlabel,options,handles);
        try S.model=gSv.vatmodel; end
        setappdata(handles.stimfig,'stimlabel',S.label);
    end
end

if isfield(S, 'model')
    [~,ix]=ismember(S.model,get(handles.modelselect,'String'));
    if ~ix
        ea_error('The model of the selected stimulation is not available.');
    else
        set(handles.modelselect,'Value',ix);
    end
else
    set(handles.modelselect,'Value',1);
end

Ractive=S.active(1);
Lactive=S.active(2);

if nargin==3
    if ischar(varargin{3})
        switch varargin{3}
            case {'Rcase'}
                ks={'k0','k1','k2','k3','k4','k5','k6','k7'};
                sidec='R'; side=1;
                S=ea_redistribute_voltage(S,varargin{3});
                if S.(['Rs',num2str(Ractive)]).case.pol==1

                    S.(['Rs',num2str(Ractive)]).case.pol=2;
                    S=ea_redistribute_voltage(S,varargin{3});
                    for k=0:7
                        if S.([sidec,'s',num2str(S.active(side))]).(['k',num2str(k)]).pol==2
                            S.([sidec,'s',num2str(S.active(side))]).(['k',num2str(k)]).pol=0;
                            S=ea_redistribute_voltage(S,['k',num2str(k)]);
                        end
                    end
                end
            case {'Lcase'}
                ks={'k8','k9','k10','k11','k12','k13','k14','k15'};
                sidec='L'; side=2;
                S=ea_redistribute_voltage(S,varargin{3});
                if S.(['Ls',num2str(Lactive)]).case.pol==1
                    S.(['Ls',num2str(Lactive)]).case.pol=2;
                    S=ea_redistribute_voltage(S,varargin{3});
                    for k=8:15
                        if S.([sidec,'s',num2str(S.active(side))]).(['k',num2str(k)]).pol==2
                            S.([sidec,'s',num2str(S.active(side))]).(['k',num2str(k)]).pol=0;
                            S=ea_redistribute_voltage(S,['k',num2str(k)]);
                        end
                    end
                end
            otherwise
                S=ea_redistribute_voltage(S,varargin{3});

                switch varargin{3}
                    case {'k0','k1','k2','k3','k4','k5','k6','k7'}
                        sidec='R'; side=1;
                    case {'k8','k9','k10','k11','k12','k13','k14','k15'}
                        sidec='L'; side=2;
                end
                if S.([sidec,'s',num2str(S.active(side))]).(varargin{3}).pol==2 && S.([sidec,'s',num2str(S.active(side))]).case.pol==2
                    S.([sidec,'s',num2str(S.active(side))]).case.pol=0;
                    S=ea_redistribute_voltage(S,[sidec,'case']);
                end
        end
    else
        S=ea_redistribute_voltage(S,varargin{3});
    end
end

setappdata(handles.stimfig,'S',S);

% set stim amplitudes
for source=1:4
    S.amplitude{1}(source)=S.(['Rs',num2str(source)]).amp;

    set(eval(['handles.Rs',num2str(source),'am']),'String',num2str(S.amplitude{1}(source)));
    set(eval(['handles.Rs',num2str(source),'va']),'Value',eval(['S.Rs',num2str(source),'.va']));

    %if eval(['S.Rs',num2str(source),'.amp']) % check if a valid +/- combination is active, if not set defaults.
    anycontactpositive=0; anycontactnegative=0;
    for k=0:7
        if eval(['S.Rs',num2str(source),'.k',num2str(k),'.pol==1'])
            anycontactnegative=1;
        elseif eval(['S.Rs',num2str(source),'.k',num2str(k),'.pol==2'])
            anycontactpositive=2;
        end
    end

    if ~anycontactnegative
        eval(['S.Rs',num2str(source),'.k1.pol=1;']);
        eval(['S.Rs',num2str(source),'.k1.perc=100;']);
    end

    if ~anycontactpositive
        eval(['S.Rs',num2str(source),'.case.pol=2;']);
        eval(['S.Rs',num2str(source),'.case.perc=100;']);
    end
    %end
end

for source=1:4
    S.amplitude{2}(source)=S.(['Ls',num2str(source)]).amp;

    set(eval(['handles.Ls',num2str(source),'am']),'String',num2str(S.amplitude{2}(source)));
    set(eval(['handles.Ls',num2str(source),'va']),'Value',eval(['S.Ls',num2str(source),'.va']));

    % if eval(['S.Ls',num2str(source),'.amp']) % check if a valid +/- combination is active, if not set defaults.
    anycontactpositive=0; anycontactnegative=0;
    for k=8:15
        if eval(['S.Ls',num2str(source),'.k',num2str(k),'.pol==1'])
            anycontactnegative=1;
        elseif eval(['S.Ls',num2str(source),'.k',num2str(k),'.pol==2'])
            anycontactpositive=1;
        end
    end

    if ~anycontactnegative
        eval(['S.Ls',num2str(source),'.k9.pol=1;']);
        eval(['S.Ls',num2str(source),'.k9.perc=100;']);
    end
    if ~anycontactpositive
        eval(['S.Ls',num2str(source),'.case.pol=2;']);
        eval(['S.Ls',num2str(source),'.case.perc=100;']);
    end
    % end
end

%% model to handles: all GUI elements.

source=Ractive;
for k=0:7
    val=eval(['S.Rs',num2str(source),'.k',num2str(k),'.perc']);
    set(eval(['handles.k',num2str(k),'u']),'String',num2str(val));

    val=eval(['S.Rs',num2str(source),'.k',num2str(k),'.imp']);
    set(eval(['handles.k',num2str(k),'im']),'String',num2str(val));
end
% set case

set(handles.RCu,'String',num2str(eval(['S.Rs',num2str(source),'.case.perc'])));

source=Lactive;
for k=8:15
    val=eval(['S.Ls',num2str(source),'.k',num2str(k),'.perc']);
    set(eval(['handles.k',num2str(k),'u']),'String',num2str(val));

    val=eval(['S.Ls',num2str(source),'.k',num2str(k),'.imp']);
    set(eval(['handles.k',num2str(k),'im']),'String',num2str(val));
end

% set case

set(handles.LCu,'String',num2str(eval(['S.Ls',num2str(source),'.case.perc'])));

%% model to handles: Axes objects:
for k=0:7
    if eval(['S.Rs',num2str(Ractive),'.k',num2str(k),'.pol==0']) % off
        im=ea_get_icn(['empty',num2str(Ractive)]);
    elseif eval(['S.Rs',num2str(Ractive),'.k',num2str(k),'.pol==1']) % negative S1
        im=ea_get_icn(['minus',num2str(Ractive)]);
    elseif eval(['S.Rs',num2str(Ractive),'.k',num2str(k),'.pol==2']) % positive S1
        im=ea_get_icn(['plus',num2str(Ractive)]);
    end
    set(0,'CurrentFigure',handles.stimfig);
    set(handles.stimfig,'CurrentAxes',eval(['handles.k',num2str(k),'ax']));
    h=image(im);
    set(h,'ButtonDownFcn',{@ea_inc_polarity,handles,options,['k',num2str(k)]});
    axis off;
    axis equal;
end

for k=8:15
    if eval(['S.Ls',num2str(Lactive),'.k',num2str(k),'.pol==0']) % off
        im=ea_get_icn(['empty',num2str(Lactive)]);
    elseif eval(['S.Ls',num2str(Lactive),'.k',num2str(k),'.pol==1']) % negative S1
        im=ea_get_icn(['minus',num2str(Lactive)]);
    elseif eval(['S.Ls',num2str(Lactive),'.k',num2str(k),'.pol==2']) % positive S1
        im=ea_get_icn(['plus',num2str(Lactive)]);
    end
    set(0,'CurrentFigure',handles.stimfig);
    set(handles.stimfig,'CurrentAxes',eval(['handles.k',num2str(k),'ax']));
    h=image(im);
    set(h,'ButtonDownFcn',{@ea_inc_polarity,handles,options,['k',num2str(k)]});
    axis off;
    axis equal;
end

% right case:
if eval(['S.Rs',num2str(Ractive),'.case.pol==0']) % off
    im=ea_get_icn(['empty',num2str(Ractive)]);
elseif eval(['S.Rs',num2str(Ractive),'.case.pol==1']) % negative
    im=ea_get_icn(['minus',num2str(Ractive)]);
elseif eval(['S.Rs',num2str(Ractive),'.case.pol==2']) % positive
    im=ea_get_icn(['plus',num2str(Ractive)]);
end
set(0,'CurrentFigure',handles.stimfig);
set(handles.stimfig,'CurrentAxes',handles.RCax);

h=image(im);
set(h,'ButtonDownFcn',{@ea_inc_polarity,handles,options,'Rcase'});
axis off;
axis equal;

% left case:
if eval(['S.Ls',num2str(Lactive),'.case.pol==0']) % off
    im=ea_get_icn(['empty',num2str(Lactive)]);
elseif eval(['S.Ls',num2str(Lactive),'.case.pol==1']) % negative
    im=ea_get_icn(['minus',num2str(Lactive)]);
elseif eval(['S.Ls',num2str(Lactive),'.case.pol==2']) % positive
    im=ea_get_icn(['plus',num2str(Lactive)]);
end
set(0,'CurrentFigure',handles.stimfig);
set(handles.stimfig,'CurrentAxes',handles.LCax);

h=image(im);
set(h,'ButtonDownFcn',{@ea_inc_polarity,handles,options,'Lcase'});
axis off;
axis equal;

%% add label

%set(handles.stimlabel,'String',S.label);

%% check consistency with chosen VAT model.
%% check consistency with chosen electrode model.
if ~isfield(options,'elspec')
    if isempty(elstruct(actpt).elmodel)
        error('Model is empty. Was the electrode segmentation fully executed in single patient mode?')
    end
    toptions=ea_resolve_elspec(elstruct(actpt));
    try
        options.elspec=toptions.elspec;
    catch
        keyboard
    end
end

models=get(handles.modelselect,'String');
try
    model=models{get(handles.modelselect,'Value')};
catch
    set(handles.modelselect,'Value',1);
    model=models{1};
end

if options.elspec.numel > 8
    warning('Only electrode with less than 8 contacts are fully supported.');
else
    ea_toggle_contacts(handles, options.elspec.numel);
end

%if strcmp(options.elspec.matfname,'boston_vercise_directed')
%    ea_error('VTA modeling for directed leads is not yet supported.');
%end

if get(handles.(['Rs',num2str(Ractive),'va']),'Value')==1 % Volt
    ea_show_percent(handles,options,1,'off'); % right hemisphere
else % Ampere
    ea_show_percent(handles,options,1,'on'); % right hemisphere
end
if get(handles.(['Ls',num2str(Ractive),'va']),'Value')==1 % Volt
    ea_show_percent(handles,options,2,'off'); % left hemisphere
else % Ampere
    ea_show_percent(handles,options,2,'on'); % left hemisphere
end

%% enable/disable panel based on sides that are present
is_side_present=arrayfun(@(xside) ~ea_arenopoints4side(elstruct(actpt).trajectory, xside), [1,2]);%First element is R, second is L
if is_side_present(1)>0%check if R side is present
    set(findall(handles.uipanel2, '-property', 'enable'), 'enable', 'on')
    %fix color (ensure they are reloaded correctly)
    try
        handles.Rs1am.BackgroundColor=[1 1 1];%force redraw color
        handles.Rs1am.BackgroundColor=[0.953 0.871 0.733];%orange
        handles.Rs2am.BackgroundColor=[1 1 1];%force redraw color
        handles.Rs2am.BackgroundColor=[0.729 0.831 0.957];%blue
        handles.Rs3am.BackgroundColor=[1 1 1];%force redraw color
        handles.Rs3am.BackgroundColor=[0.925 0.839 0.839];%red
        handles.Rs4am.BackgroundColor=[1 1 1];%force redraw color
        handles.Rs4am.BackgroundColor=[0.757 0.867 0.776];%green
    catch
        %for older versions of matlab (remove if dropping support of older Matlabs)
        set(findall(handles.Rs1am, '-property', 'BackgroundColor'), 'BackgroundColor', [1 1 1]);%force redraw color
        set(findall(handles.Rs1am, '-property', 'BackgroundColor'), 'BackgroundColor', [0.953 0.871 0.733]);%orange
        set(findall(handles.Rs2am, '-property', 'BackgroundColor'), 'BackgroundColor', [1 1 1]);%force redraw color
        set(findall(handles.Rs2am, '-property', 'BackgroundColor'), 'BackgroundColor', [0.729 0.831 0.957]);%blue
        set(findall(handles.Rs3am, '-property', 'BackgroundColor'), 'BackgroundColor', [1 1 1]);%force redraw color
        set(findall(handles.Rs3am, '-property', 'BackgroundColor'), 'BackgroundColor', [0.925 0.839 0.839]);%red
        set(findall(handles.Rs4am, '-property', 'BackgroundColor'), 'BackgroundColor', [1 1 1]);%force redraw color
        set(findall(handles.Rs4am, '-property', 'BackgroundColor'), 'BackgroundColor', [0.757 0.867 0.776]);%green
    end
else
    set(findall(handles.uipanel2, '-property', 'enable'), 'enable', 'off')
end
if is_side_present(2)>0%check if L side is present
    set(findall(handles.uipanel3, '-property', 'enable'), 'enable', 'on')

    %fix color (ensure they are reloaded correctly)
    try
        handles.Ls1am.BackgroundColor=[1 1 1];%force redraw color
        handles.Ls1am.BackgroundColor=[0.953 0.871 0.733];%orange
        handles.Ls2am.BackgroundColor=[1 1 1];%force redraw color
        handles.Ls2am.BackgroundColor=[0.729 0.831 0.957];%blue
        handles.Ls3am.BackgroundColor=[1 1 1];%force redraw color
        handles.Ls3am.BackgroundColor=[0.925 0.839 0.839];%red
        handles.Ls4am.BackgroundColor=[1 1 1];%force redraw color
        handles.Ls4am.BackgroundColor=[0.757 0.867 0.776];%green
    catch
        %for older versions of matlab (remove if dropping support of older Matlabs)
        set(findall(handles.Ls1am, '-property', 'BackgroundColor'), 'BackgroundColor', [1 1 1]);%force redraw color
        set(findall(handles.Ls1am, '-property', 'BackgroundColor'), 'BackgroundColor', [0.953 0.871 0.733]);%orange
        set(findall(handles.Ls2am, '-property', 'BackgroundColor'), 'BackgroundColor', [1 1 1]);%force redraw color
        set(findall(handles.Ls2am, '-property', 'BackgroundColor'), 'BackgroundColor', [0.729 0.831 0.957]);%blue
        set(findall(handles.Ls3am, '-property', 'BackgroundColor'), 'BackgroundColor', [1 1 1]);%force redraw color
        set(findall(handles.Ls3am, '-property', 'BackgroundColor'), 'BackgroundColor', [0.925 0.839 0.839]);%red
        set(findall(handles.Ls4am, '-property', 'BackgroundColor'), 'BackgroundColor', [1 1 1]);%force redraw color
        set(findall(handles.Ls4am, '-property', 'BackgroundColor'), 'BackgroundColor', [0.757 0.867 0.776]);%green
    end
else
    set(findall(handles.uipanel3, '-property', 'enable'), 'enable', 'off')
end

switch model
    case 'SimBio/FieldTrip (see Horn 2017)'
        ea_hide_impedance(handles);
        set(handles.estimateInTemplate,'Visible','on');
        S.monopolarmodel=0;
        ea_enable_vas(handles,options);
        set(handles.betawarning,'visible','on');
        set(handles.settings,'visible','on');
        set(handles.addStimSet,'visible','off');
    case 'Maedler 2012'
        ea_show_impedance(handles);
        set(handles.estimateInTemplate,'Visible','off');
        S.monopolarmodel=1;
        ea_disable_vas(handles,options);
        set(handles.betawarning,'visible','off');
        set(handles.settings,'visible','off');
        set(handles.addStimSet,'visible','off');
    case 'Kuncel 2008'
        ea_hide_impedance(handles);
        set(handles.estimateInTemplate,'Visible','off');
        S.monopolarmodel=1;
        ea_disable_vas(handles,options);
        set(handles.betawarning,'visible','off');
        set(handles.settings,'visible','off');
        set(handles.addStimSet,'visible','off');
    case 'Dembek 2017'
        ea_show_impedance(handles);
        set(handles.estimateInTemplate,'Visible','off');
        S.monopolarmodel=1;
        ea_enable_vas(handles,options);
        set(handles.betawarning,'visible','off');
        set(handles.settings,'visible','on');
        set(handles.addStimSet,'visible','off');
    case 'Fastfield (Baniasadi 2020)'
        ea_show_impedance(handles);
        set(handles.estimateInTemplate,'Visible','off');
        S.monopolarmodel=0;
        ea_enable_vas(handles,options);
        set(handles.betawarning,'visible','off');
        set(handles.settings,'visible','on');
        set(handles.addStimSet,'visible','off');
    case 'OSS-DBS (Butenko 2020)'
        ea_hide_impedance(handles);
        set(handles.estimateInTemplate,'Visible','on');
        S.monopolarmodel=0;
        ea_enable_vas(handles,options);
        set(handles.betawarning,'visible','on');
        set(handles.settings,'visible','on');
        set(handles.addStimSet,'visible','on');

end
S.model=model;

ea_savestimulation(S,options);
setappdata(handles.stimfig,'S',S);


function ea_show_percent(handles,options,side,onoff)

switch side
    case 1
        sel=0:7;
        sidestr='R';
        ptval=1;
    case 2
        sel=8:15;
        sidestr='L';
        ptval=3;
end

% Only support up to 8 contacts for now
if options.elspec.numel<=8
    sel=sel(1:options.elspec.numel);
else
    sel=sel(1:8);
end

for k=sel
    set(handles.(['k',num2str(k),'u']),'visible',onoff);
end

set(handles.([sidestr,'Cu']),'visible',onoff);

set(handles.(['perctext',num2str(ptval)]),'visible',onoff);
if options.elspec.numel>4
    set(handles.(['perctext',num2str(ptval+1)]),'visible',onoff);
end


function ea_toggle_contacts(handles, numel)

if numel == 8
    status = 'on';
else
    status = 'off';
end


for k=[numel:7,numel+8:15]
    set(handles.(['k',num2str(k),'u']),'visible',status);
    set(handles.(['k',num2str(k),'im']),'visible',status);
    set(handles.(['k',num2str(k),'txt']),'visible',status);
    handles2hide = [get(handles.(['k',num2str(k),'ax']),'Children')];
    set(handles2hide,'visible',status)
end

set(handles.perctext2,'visible',status);
set(handles.kohmtext2,'visible',status);
set(handles.perctext4,'visible',status);
set(handles.kohmtext4,'visible',status);


function ea_disable_vas(handles,options)

RL={'R','L'};
for iside=1:length(options.sides)
    side=options.sides(iside);

    for Rva=1:4
        set(handles.([RL{side},'s',num2str(Rva),'va']),'enable','off');
        set(handles.([RL{side},'s',num2str(Rva),'va']),'value',1);
    end
end


function ea_enable_vas(handles,options)

RL={'R','L'};
for iside=1:length(options.sides)
    side=options.sides(iside);

    for Rva=1:4
        set(handles.([RL{side},'s',num2str(Rva),'va']),'enable','on');
    end
end


function ea_hide_impedance(handles)

for k=0:15
    eval(['set(handles.k',num2str(k),'im,''visible'',''off'');']);
end

for ohm=1:4
    eval(['set(handles.kohmtext',num2str(ohm),',''visible'',''off'');']);
end


function ea_show_impedance(handles)

for k=0:15
    eval(['set(handles.k',num2str(k),'im,''visible'',''on'');']);
end

for ohm=1:4
    eval(['set(handles.kohmtext',num2str(ohm),',''visible'',''on'');']);
end


function S=ea_redistribute_voltage(S,changedobj)
Rconts={'k0','k1','k2','k3','k4','k5','k6','k7'};
Lconts={'k8','k9','k10','k11','k12','k13','k14','k15'};
LcontsCase=[Lconts,{'case'}];
RcontsCase=[Rconts,{'case'}];
if ischar(changedobj) % different polarity on the block
    switch changedobj
        case Rconts
            conts=Rconts;
            contsCase=RcontsCase;
            sidec='R';
            side=1;
        case Lconts
            conts=Lconts;
            contsCase=LcontsCase;
            sidec='L';
            side=2;
        case 'Rcase'
            conts=Rconts;
            changedobj='case';
            contsCase=RcontsCase;

            side=1;
            sidec='R';
        case 'Lcase'
            conts=Lconts;
            contsCase=LcontsCase;

            changedobj='case';
            side=2;
            sidec='L';
    end

    % check polarity of changed object:
    polchanged=eval(['S.',sidec,'s',num2str(S.active(side)),'.',changedobj,'.pol']);

    % check for monopolar models:
    if S.monopolarmodel % these allow only 1 active anode contact per model.
        for c=1:length(conts)
            eval(['S.',sidec,'s',num2str(S.active(side)),'.',conts{c},'.pol=0;']);
            eval(['S.',sidec,'s',num2str(S.active(side)),'.',conts{c},'.perc=0;']);
        end
        eval(['S.',sidec,'s',num2str(S.active(side)),'.',changedobj,'.pol=1;']);
        eval(['S.',sidec,'s',num2str(S.active(side)),'.',changedobj,'.perc=100;']);

        return

    else
        %         if S.([sidec,'s',num2str(S.active(side))]).va==2 % ampere only allows one anode and one cathode
        %             for c=1:length(contsCase)
        %
        %                 if S.([sidec,'s',num2str(S.active(side))]).(contsCase{c}).pol==polchanged % same polarity as changed object
        %                     S.([sidec,'s',num2str(S.active(side))]).(contsCase{c}).pol=ea_swappol(polchanged);
        %                     S.([sidec,'s',num2str(S.active(side))]).(contsCase{c}).perc=100;
        %                 else
        %                     S.([sidec,'s',num2str(S.active(side))]).(contsCase{c}).pol=0;
        %                     S.([sidec,'s',num2str(S.active(side))]).(contsCase{c}).perc=0;
        %                 end
        %             end
        %             S.([sidec,'s',num2str(S.active(side))]).(changedobj).pol=1;
        %             S.([sidec,'s',num2str(S.active(side))]).(changedobj).perc=100;
        %         end
    end

    if polchanged==0
        % set changed contacts percentage to zero:
        eval(['S.',sidec,'s',num2str(S.active(side)),'.',changedobj,'.perc=0;']);
    else
        % determine how many other nodes with this polarity exist:
        divby=1;
        contacts={};
        for con=1:length(conts)
            if eval(['S.',sidec,'s',num2str(S.active(side)),'.',conts{con},'.pol==polchanged'])
                if ~strcmp(conts{con},changedobj)
                    %voltages{divby}=eval(['S.Rs',num2str(S.active(side)),'.',Rconts{con},'.perc']);
                    contacts{divby}=conts{con};
                    divby=divby+1;
                end
            end
        end

        if eval(['S.',sidec,'s',num2str(S.active(side)),'.case.pol==polchanged'])
            if ~strcmp(changedobj,'case')
                contacts{divby}='case';
                divby=divby+1;
            end
        end
        % add case to calculation.

        % set changed contacts percentage:
        eval(['S.',sidec,'s',num2str(S.active(side)),'.',changedobj,'.perc=100/divby;']);

        % reduce all other contacts percentages:

        try divby=divby/length(contacts); end
        for c=1:length(contacts)
            eval(['S.',sidec,'s',num2str(S.active(side)),'.',contacts{c},'.perc=',...
                'S.',sidec,'s',num2str(S.active(side)),'.',contacts{c},'.perc/divby;']);
        end
    end

    % now clean up mess from polarity that the contact used to have..

    polchanged=ea_polminus(polchanged);
    sumpercs=0;

    if polchanged % polarization has changed from negative to positive. clean up negatives. or changed from positive to off. clean up positives.
        contacts={};
        cnt=0;
        for con=1:length(conts)
            if eval(['S.',sidec,'s',num2str(S.active(side)),'.',conts{con},'.pol==polchanged'])
                if ~strcmp(conts{con},changedobj)
                    %voltages{divby}=eval(['S.Rs',num2str(S.active(side)),'.',Rconts{con},'.perc']);

                    cnt=cnt+1;
                    contacts{cnt}=conts{con};
                    sumpercs=sumpercs+eval(['S.',sidec,'s',num2str(S.active(side)),'.',conts{con},'.perc']);
                end
            end
        end
        % add case to calculation:
        if eval(['S.',sidec,'s',num2str(S.active(side)),'.case.pol==polchanged'])
            if ~strcmp(changedobj,'case')
                cnt=cnt+1;
                contacts{cnt}='case';
                sumpercs=sumpercs+eval(['S.',sidec,'s',num2str(S.active(side)),'.case.perc']);
            end
        end

        multby=(100/sumpercs);
        if cnt
            for c=1:length(contacts)
                eval(['S.',sidec,'s',num2str(S.active(side)),'.',contacts{c},'.perc=',...
                    'S.',sidec,'s',num2str(S.active(side)),'.',contacts{c},'.perc*multby;']);
            end
        end
    end

else % voltage percentage changed
    changedobj=get(changedobj,'Tag');
    changedobj=changedobj(1:end-1);

    switch changedobj
        case Rconts
            conts=Rconts;
            sidec='R';
            side=1;
        case Lconts
            conts=Lconts;
            sidec='L';
            side=2;
        case 'RC'
            conts=Rconts;
            changedobj='case';
            side=1;
            sidec='R';
        case 'LC'
            conts=Lconts;
            changedobj='case';
            side=2;
            sidec='L';
    end

    % check for monopolar models:
    if S.monopolarmodel % these allow only 1 active anode contact per model.
        for c=1:length(conts)
            eval(['S.',sidec,'s',num2str(S.active(side)),'.',conts{c},'.pol=0;']);
            eval(['S.',sidec,'s',num2str(S.active(side)),'.',conts{c},'.perc=0;']);
        end
        eval(['S.',sidec,'s',num2str(S.active(side)),'.',changedobj,'.pol=1;']);
        eval(['S.',sidec,'s',num2str(S.active(side)),'.',changedobj,'.perc=100;']);

        return
    end

    % check polarity of changed object:
    try
        polchanged=eval(['S.',sidec,'s',num2str(S.active(side)),'.',changedobj,'.pol']);
    catch
        keyboard
    end

    if polchanged==0 % set changed contacts polarity to negative
        eval(['S.',sidec,'s',num2str(S.active(side)),'.',changedobj,'.pol=1;']);
        polchanged=1;
    end

    % determine how many other nodes with this polarity exist:
    divby=1;
    contacts={};
    sumpercent=0;
    for con=1:length(conts)
        if eval(['S.',sidec,'s',num2str(S.active(side)),'.',conts{con},'.pol==polchanged'])
            if ~strcmp(conts{con},changedobj)
                sumpercent=sumpercent+eval(['S.',sidec,'s',num2str(S.active(side)),'.',conts{con},'.perc']);
                contacts{divby}=conts{con};
                divby=divby+1;
            end
        end
    end

    % add case to calculation.
    if eval(['S.',sidec,'s',num2str(S.active(side)),'.case.pol==polchanged'])
        if ~strcmp(changedobj,'case')
            contacts{divby}='case';
            divby=divby+1;
        end
    end

    if divby==1 % only one contact -> set to 100 percent.
        eval(['S.',sidec,'s',num2str(S.active(side)),'.',changedobj,'.perc=100;']);
    end

    % reduce all other contacts percentages:
    divby=sumpercent/(100-eval(['S.',sidec,'s',num2str(S.active(side)),'.',changedobj,'.perc']));

    for c=1:length(contacts)
        eval(['S.',sidec,'s',num2str(S.active(side)),'.',contacts{c},'.perc=',...
            'S.',sidec,'s',num2str(S.active(side)),'.',contacts{c},'.perc/divby;']);
    end
end


function opol=ea_polminus(pol)
if pol==1
    opol=0;
elseif pol==2
    opol=1;
elseif pol==0
    opol=2;
end


function opol=ea_swappol(pol)
if pol==1
    opol=2;
elseif pol==2
    opol=1;
else
    opol=0;
end


function ea_inc_polarity(h,h2,handles,options,ID)

S=getappdata(handles.stimfig,'S');

switch ID
    case {'k0','k1','k2','k3','k4','k5','k6','k7'}
        side=1;
        sidec='R';
        gID=ID;
        ks=0:7;
    case {'k8','k9','k10','k11','k12','k13','k14','k15'}
        side=2;
        sidec='L';
        gID=ID;
        ks=8:15;
    case 'Rcase'
        gID='case';
        side=1;
        sidec='R';
        ks=0:7;
    case 'Lcase'
        gID='case';
        side=2;
        sidec='L';
        ks=8:15;
end

cycles=[0,1,2];

try
    oldval=eval(['S.',sidec,'s',num2str(S.active(side)),'.',gID,'.pol']);
catch
    keyboard
end
[~,newval]=ismember(oldval,cycles);
newval=newval+1;
if newval>length(cycles)
    newval=1;
end

newval=cycles(newval);

eval(['S.',sidec,'s',num2str(S.active(side)),'.',gID,'.pol=',num2str(newval),';']);

% now check if any other contact is left with the old polarity

anycontactpositive=0; anycontactnegative=0;
for k=ks
    if eval(['S.',sidec,'s',num2str(num2str(S.active(side))),'.k',num2str(k),'.pol==1'])
        anycontactnegative=1;
    elseif eval(['S.',sidec,'s',num2str(num2str(S.active(side))),'.k',num2str(k),'.pol==2'])
        anycontactpositive=1;
    end
end

% also check case
if eval(['S.',sidec,'s',num2str(num2str(S.active(side))),'.case.pol==1'])
    anycontactnegative=1;
elseif eval(['S.',sidec,'s',num2str(num2str(S.active(side))),'.case.pol==2'])
    anycontactpositive=1;
end

if anycontactnegative && anycontactpositive % only then save results..
    setappdata(handles.stimfig,'S',S);
end
ea_refreshguisp(handles,options,ID);


function Ls3am_Callback(hObject, eventdata, handles)
% hObject    handle to Ls3am (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ls3am as text
%        str2double(get(hObject,'String')) returns contents of Ls3am as a double

S=getappdata(handles.stimfig,'S');
options=getappdata(handles.stimfig,'options');
S.active(2)=3;
S.Ls3.amp=str2double(get(hObject,'String'));

setappdata(handles.stimfig,'S',S);
ea_refreshguisp(handles,options);


% --- Executes during object creation, after setting all properties.
function Ls3am_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ls3am (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Ls3va.
function Ls3va_Callback(hObject, eventdata, handles)
% hObject    handle to Ls3va (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Ls3va contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Ls3va
S=getappdata(handles.stimfig,'S');
options=getappdata(handles.stimfig,'options');
S.active(2)=3;
S.Ls3.va=get(hObject,'Value');

setappdata(handles.stimfig,'S',S);
ea_refreshguisp(handles,options);


% --- Executes during object creation, after setting all properties.
function Ls3va_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ls3va (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Ls4am_Callback(hObject, eventdata, handles)
% hObject    handle to Ls4am (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ls4am as text
%        str2double(get(hObject,'String')) returns contents of Ls4am as a double
S=getappdata(handles.stimfig,'S');
options=getappdata(handles.stimfig,'options');
S.active(2)=4;
S.Ls4.amp=str2double(get(hObject,'String'));

setappdata(handles.stimfig,'S',S);
ea_refreshguisp(handles,options);


% --- Executes during object creation, after setting all properties.
function Ls4am_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ls4am (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Ls4va.
function Ls4va_Callback(hObject, eventdata, handles)
% hObject    handle to Ls4va (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Ls4va contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Ls4va
S=getappdata(handles.stimfig,'S');
options=getappdata(handles.stimfig,'options');
S.active(2)=4;
S.Ls4.va=get(hObject,'Value');

setappdata(handles.stimfig,'S',S);
ea_refreshguisp(handles,options);


% --- Executes during object creation, after setting all properties.
function Ls4va_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ls4va (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Rs3am_Callback(hObject, eventdata, handles)
% hObject    handle to Rs3am (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Rs3am as text
%        str2double(get(hObject,'String')) returns contents of Rs3am as a double
S=getappdata(handles.stimfig,'S');
options=getappdata(handles.stimfig,'options');
S.active(1)=3;
S.Rs3.amp=str2double(get(hObject,'String'));

setappdata(handles.stimfig,'S',S);
ea_refreshguisp(handles,options);


% --- Executes during object creation, after setting all properties.
function Rs3am_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Rs3am (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Rs3va.
function Rs3va_Callback(hObject, eventdata, handles)
% hObject    handle to Rs3va (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Rs3va contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Rs3va
S=getappdata(handles.stimfig,'S');
options=getappdata(handles.stimfig,'options');
S.active(1)=3;
S.Rs3.va=get(hObject,'Value');
setappdata(handles.stimfig,'S',S);
ea_refreshguisp(handles,options);


% --- Executes during object creation, after setting all properties.
function Rs3va_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Rs3va (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Rs4am_Callback(hObject, eventdata, handles)
% hObject    handle to Rs4am (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Rs4am as text
%        str2double(get(hObject,'String')) returns contents of Rs4am as a double
S=getappdata(handles.stimfig,'S');
options=getappdata(handles.stimfig,'options');
S.active(1)=4;
S.Rs4.amp=str2double(get(hObject,'String'));

setappdata(handles.stimfig,'S',S);
ea_refreshguisp(handles,options);


% --- Executes during object creation, after setting all properties.
function Rs4am_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Rs4am (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Rs4va.
function Rs4va_Callback(hObject, eventdata, handles)
% hObject    handle to Rs4va (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Rs4va contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Rs4va
S=getappdata(handles.stimfig,'S');
options=getappdata(handles.stimfig,'options');
S.active(1)=4;
S.Rs4.va=get(hObject,'Value');

setappdata(handles.stimfig,'S',S);
ea_refreshguisp(handles,options);


% --- Executes during object creation, after setting all properties.
function Rs4va_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Rs4va (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in nextpt.
function nextpt_Callback(hObject, eventdata, handles)
% hObject    handle to nextpt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

S=getappdata(handles.stimfig,'S');
gS=getappdata(handles.stimfig,'gS');

gSv=getappdata(handles.stimfig,'gSv');

actpt=getappdata(handles.stimfig,'actpt');
elstruct=getappdata(handles.stimfig,'elstruct');
options=getappdata(handles.stimfig,'options');

%ensure active patient is non empty
if isempty(actpt)
    actpt=1;
end

if isempty(gS)
    clear gS
end

S=ea_activecontacts(S);
try
    gS(actpt)=S;
catch
    S.sources=1:4;
    S.volume=[0,0];
    gS(actpt)=S;
end
setappdata(handles.stimfig,'gS',gS);

if (actpt+1)>length(elstruct)
    setto=1;
else
    setto=actpt+1;
end
try
    setappdata(handles.stimfig,'S',gS(setto));
catch % S for this next pt not defined yet.
    setappdata(handles.stimfig,'S',[]);
end
setappdata(handles.stimfig,'actpt',setto);

ea_refreshguisp(handles,options);


% --- Executes on button press in prevpt.
function prevpt_Callback(hObject, eventdata, handles)
% hObject    handle to prevpt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

S=getappdata(handles.stimfig,'S');
gS=getappdata(handles.stimfig,'gS');
actpt=getappdata(handles.stimfig,'actpt');
elstruct=getappdata(handles.stimfig,'elstruct');
options=getappdata(handles.stimfig,'options');

%ensure active patient is non empty
if isempty(actpt)
    actpt=1;
end

if isempty(gS)
    clear gS
end

S=ea_activecontacts(S);
gS(actpt)=S;
setappdata(handles.stimfig,'gS',gS);

if (actpt-1)<1
    setto=length(elstruct);
else
    setto=actpt-1;
end

try
    setappdata(handles.stimfig,'S',gS(setto));
catch
    setappdata(handles.stimfig,'S',[]);
end
setappdata(handles.stimfig,'actpt',setto);

ea_refreshguisp(handles,options);


% --- Executes on button press in saveparams.
function saveparams_Callback(hObject, eventdata, handles)
% hObject    handle to saveparams (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
S = getappdata(handles.stimfig, 'S');
gS = getappdata(handles.stimfig, 'gS');
actpt = getappdata(handles.stimfig, 'actpt');
elstruct = getappdata(handles.stimfig, 'elstruct');
options = getappdata(handles.stimfig, 'options');

if isempty(gS)
    clear gS
end

S = ea_activecontacts(S);

try
    gS(actpt)=S;
catch
    S.sources=1:4;
    S.volume=[0,0];
    gS(actpt)=S;
end
setappdata(handles.stimfig, 'gS', gS);

gSv = getappdata(handles.stimfig, 'gSv');
lgfig = getappdata(handles.stimfig, 'resultfig');

setappdata(lgfig, 'S', gS);
setappdata(lgfig, 'vatmodel', gSv.vatmodel);
setstimbutton = lgfig.findobj('Tag','setstimparamsbutton');
if ~isempty(gS)
    set(setstimbutton,'BackgroundColor',[0.1;0.8;0.1]);
else
    set(setstimbutton,'BackgroundColor',[0.93,0.93,0.93]);
end

ea_setprefs('vatsettings.estimateInTemplate',get(handles.estimateInTemplate,'Value'));

close(handles.stimfig);


% --- Executes on key press with focus on Rs2am and none of its controls.
function Rs2am_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to Rs2am (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
S=getappdata(handles.stimfig,'S');
S.active(1)=2;
S.Rs2.amp=str2double(get(hObject,'String'));
setappdata(handles.stimfig,'S',S);
options=getappdata(handles.stimfig,'options');
ea_refreshguisp(handles,options);


% --- Executes on key press with focus on Rs1am and none of its controls.
function Rs1am_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to Rs1am (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
S=getappdata(handles.stimfig,'S');
S.active(1)=1;
S.Rs1.amp=str2double(get(hObject,'String'));
setappdata(handles.stimfig,'S',S);
options=getappdata(handles.stimfig,'options');
ea_refreshguisp(handles,options);


% --- Executes on key press with focus on Rs3am and none of its controls.
function Rs3am_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to Rs3am (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
S=getappdata(handles.stimfig,'S');
S.active(1)=3;
S.Rs3.amp=str2double(get(hObject,'String'));
setappdata(handles.stimfig,'S',S);
options=getappdata(handles.stimfig,'options');
ea_refreshguisp(handles,options);


% --- Executes on key press with focus on Rs4am and none of its controls.
function Rs4am_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to Rs4am (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
S=getappdata(handles.stimfig,'S');
S.active(1)=4;
S.Rs4.amp=str2double(get(hObject,'String'));
setappdata(handles.stimfig,'S',S);
options=getappdata(handles.stimfig,'options');
ea_refreshguisp(handles,options);


% --- Executes on key press with focus on Ls1am and none of its controls.
function Ls1am_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to Ls1am (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
S=getappdata(handles.stimfig,'S');
S.active(2)=1;
S.Ls1.amp=str2double(get(hObject,'String'));
setappdata(handles.stimfig,'S',S);
options=getappdata(handles.stimfig,'options');
ea_refreshguisp(handles,options);


% --- Executes on key press with focus on Ls2am and none of its controls.
function Ls2am_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to Ls2am (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
S=getappdata(handles.stimfig,'S');
S.active(2)=2;
S.Ls2.amp=str2double(get(hObject,'String'));
setappdata(handles.stimfig,'S',S);
options=getappdata(handles.stimfig,'options');
ea_refreshguisp(handles,options);


% --- Executes on key press with focus on Ls3am and none of its controls.
function Ls3am_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to Ls3am (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
S=getappdata(handles.stimfig,'S');
S.active(2)=3;
S.Ls3.amp=str2double(get(hObject,'String'));
setappdata(handles.stimfig,'S',S);
options=getappdata(handles.stimfig,'options');
ea_refreshguisp(handles,options);


% --- Executes on key press with focus on Ls4am and none of its controls.
function Ls4am_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to Ls4am (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
S=getappdata(handles.stimfig,'S');
S.active(2)=4;
S.Ls4.amp=str2double(get(hObject,'String'));
setappdata(handles.stimfig,'S',S);
options=getappdata(handles.stimfig,'options');
ea_refreshguisp(handles,options);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over Rs1am.
function Rs1am_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to Rs1am (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
S=getappdata(handles.stimfig,'S');
S.active(1)=1;
S.Rs1.amp=str2double(get(hObject,'String'));
setappdata(handles.stimfig,'S',S);
options=getappdata(handles.stimfig,'options');
ea_refreshguisp(handles,options);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over Rs2am.
function Rs2am_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to Rs2am (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
S=getappdata(handles.stimfig,'S');
S.active(1)=2;
S.Rs1.amp=str2double(get(hObject,'String'));
setappdata(handles.stimfig,'S',S);
options=getappdata(handles.stimfig,'options');
ea_refreshguisp(handles,options);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over Rs3am.
function Rs3am_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to Rs3am (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
S=getappdata(handles.stimfig,'S');
S.active(1)=3;
S.Rs3.amp=str2double(get(hObject,'String'));
setappdata(handles.stimfig,'S',S);
options=getappdata(handles.stimfig,'options');
ea_refreshguisp(handles,options);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over Rs4am.
function Rs4am_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to Rs4am (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
S=getappdata(handles.stimfig,'S');
S.active(1)=4;
S.Rs4.amp=str2double(get(hObject,'String'));
setappdata(handles.stimfig,'S',S);
options=getappdata(handles.stimfig,'options');
ea_refreshguisp(handles,options);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over Ls1am.
function Ls1am_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to Ls1am (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
S=getappdata(handles.stimfig,'S');
S.active(2)=1;
S.Ls1.amp=str2double(get(hObject,'String'));
setappdata(handles.stimfig,'S',S);
options=getappdata(handles.stimfig,'options');
ea_refreshguisp(handles,options);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over Ls2am.
function Ls2am_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to Ls2am (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
S=getappdata(handles.stimfig,'S');
S.active(2)=2;
S.Ls2.amp=str2double(get(hObject,'String'));
setappdata(handles.stimfig,'S',S);
options=getappdata(handles.stimfig,'options');
ea_refreshguisp(handles,options);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over Ls3am.
function Ls3am_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to Ls3am (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
S=getappdata(handles.stimfig,'S');
S.active(2)=3;
S.Ls3.amp=str2double(get(hObject,'String'));
setappdata(handles.stimfig,'S',S);
options=getappdata(handles.stimfig,'options');
ea_refreshguisp(handles,options);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over Ls4am.
function Ls4am_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to Ls4am (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
S=getappdata(handles.stimfig,'S');
S.active(2)=4;
S.Ls4.amp=str2double(get(hObject,'String'));
setappdata(handles.stimfig,'S',S);
options=getappdata(handles.stimfig,'options');
ea_refreshguisp(handles,options);


% --- Executes on button press in predictstim.
function predictstim_Callback(hObject, eventdata, handles)
% hObject    handle to predictstim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in settings.
function settings_Callback(hObject, eventdata, handles)
% hObject    handle to settings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
models=get(handles.modelselect,'String');
try
    model=models{get(handles.modelselect,'Value')};
catch
    set(handles.modelselect,'Value',1);
    model=models{1};
end
switch model
    case 'SimBio/FieldTrip (see Horn 2017)'
        ea_vatsettings_horn;
    case 'Dembek 2017'
        ea_vatsettings_dembek;
    case 'Fastfield (Baniasadi 2020)'
        ea_vatsettings_fastfield;
    case 'OSS-DBS (Butenko 2020)'
        ea_vatsettings_butenko;
end


% --- Executes on button press in estimateInTemplate.
function estimateInTemplate_Callback(hObject, eventdata, handles)
% hObject    handle to estimateInTemplate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of estimateInTemplate


% --- Executes when user attempts to close stimfig.
function stimfig_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to stimfig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
ea_setprefs('vatsettings.estimateInTemplate',get(handles.estimateInTemplate,'Value'));
delete(hObject);


% --- Executes on button press in addStimSet.
function addStimSet_Callback(hObject, eventdata, handles)
% hObject    handle to addStimSet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of addStimSet
if hObject.Value
    options = getappdata(handles.stimfig, 'options');
    stimLabel = getappdata(handles.stimfig, 'stimlabel');

    if strcmp(options.leadprod, 'dbs')
        patdir = [options.root, options.patientname];
        numContacts = options.elspec.numel;
    else % Special case when call from Lead Group
        actpt = getappdata(handles.stimfig, 'actpt');
        resultfig = getappdata(handles.stimfig, 'resultfig');
        M = getappdata(resultfig, 'M');
        patdir = M.patient.list{actpt};
        numContacts = size(M.elstruct(actpt).coords_mm{1}, 1);
    end

    stimFolder = fullfile(patdir ,'stimulations', ea_nt(~handles.estimateInTemplate.Value), stimLabel);
    ea_mkdir(stimFolder);
    ea_addStimSet(numContacts, stimFolder, hObject);

    % Set stimSetMode flag for Lead Group
    if hObject.Value && strcmp(options.leadprod, 'group')
        M.ui.stimSetMode = 1;
        setappdata(resultfig, 'M');
    end
end
