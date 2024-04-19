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

% Last Modified by GUIDE v2.5 16-Feb-2024 15:14:00

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
    if ~ismember(options.elmodel,ea_ossdbs_elmodel)
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

if handles.estimateInTemplate.Visible % only allowed for specific VTA functions
    switch handles.estimateInTemplate.Value
        case 0
            options.native = 1;
        case 1
            options.native = 0;
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
        if options.prefs.machine.vatsettings.butenko_calcPAM
            feval(ea_genvat,getappdata(handles.stimfig,'S'),options,handles.stimfig);
            ea_busyaction('off',handles.stimfig,'stim');
            return;
        else
            getappdata(handles.stimfig,'S')
            S_old = S;
            file_path = '/Users/savirmadan/Downloads/exportedData.json';
            new_data = fileread(file_path);
            importedS = jsondecode(new_data);
            S = importedS.S;
            S.label = S_old.label;
%             S.model = S_old.model;
            S.amplitude = {S.amplitude.rightAmplitude.', S.amplitude.leftAmplitude.'};
%             S.amplitude{2}(2:)
%             S.activecontacts = S.activecontacts{1};
            S.monopolarmodel = 0;
%             S.activecontacts = {cell2mat(S.activecontacts(6)).', cell2mat(S.activecontacts(2)).'};
            % Encode new data to JSON format
            [~, stimparams] = feval(ea_genvat,S,options,handles.stimfig);
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
        ea_vatsettings_butenko(handles);
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


% --- Executes during object creation, after setting all properties.
function pulseWidthTextbox_R_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pulseWidthTextbox_R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function pulseWidthTextbox_R_Callback(hObject, eventdata, handles)
% hObject    handle to pulseWidthTextbox_R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pulseWidthTextbox_R as text
%        str2double(get(hObject,'String')) returns contents of pulseWidthTextbox_R as a double
S=getappdata(handles.stimfig,'S');
options=getappdata(handles.stimfig,'options');

eval(['S.Rs',num2str(S.active(1)),'.pulseWidth=',hObject.String,';']);

setappdata(handles.stimfig,'S',S);
ea_refreshguisp(handles,options,hObject);


% --- Executes during object creation, after setting all properties.
function pulseWidthTextbox_L_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pulseWidthTextbox_L (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function pulseWidthTextbox_L_Callback(hObject, eventdata, handles)
% hObject    handle to pulseWidthTextbox_L (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pulseWidthTextbox_L as text
%        str2double(get(hObject,'String')) returns contents of pulseWidthTextbox_L as a double
S=getappdata(handles.stimfig,'S');
options=getappdata(handles.stimfig,'options');

eval(['S.Ls',num2str(S.active(2)),'.pulseWidth=',hObject.String,';']);

setappdata(handles.stimfig,'S',S);
ea_refreshguisp(handles,options,hObject);
