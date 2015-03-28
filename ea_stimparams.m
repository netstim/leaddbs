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

% Last Modified by GUIDE v2.5 11-Mar-2014 17:55:53

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
set(gcf,'Name','Stimulation Parameters');

% store input arguments in figure to make it available to subroutines.
elstruct=varargin{1};
resultfig=varargin{2};
options=varargin{3};
setappdata(gcf,'elstruct',elstruct);
setappdata(gcf,'resultfig',resultfig);
setappdata(gcf,'options',options);

stimparams=getappdata(resultfig,'stimparams'); % get info from resultfig.
setappdata(gcf,'stimparams',stimparams); % store stimulation settings from resultfig to stim (this) fig for subroutines.



% setup modelselect popup

cnt=1;
earoot=[fileparts(which('lead')),filesep];
ndir=dir([earoot,'ea_genvat_*.m']);
for nd=length(ndir):-1:1
    [~,methodf]=fileparts(ndir(nd).name);
    try
        [thisndc]=eval([methodf,'(','''prompt''',')']);
        ndc{cnt}=thisndc;
        genvatfunctions{cnt}=methodf;
        cnt=cnt+1;
    end
end
setappdata(gcf,'genvatfunctions',genvatfunctions);

set(handles.modelselect,'String',ndc);

%



if ~isempty(stimparams) % stimfigure has been used before..
    for side=1:2
        for el=1:4
            %keyboard
            
            set(eval(['handles.k',num2str(((side-1)*4)+el-1),'u']),'String', num2str(stimparams(side).U(el)));
            set(eval(['handles.k',num2str(((side-1)*4)+el-1),'im']),'String',num2str(stimparams(side).Im(el)));
            
            
        end
    end
    
    
    set(handles.fiberthresh,'String',num2str(stimparams(1).fiberthresh))
    
    set(handles.showfibs,'Value',stimparams(1).showfibers);
    set(handles.showconns,'Value',stimparams(1).showconnectivities);
end


% Build popup tables:
% Fibers:


fibd=dir([options.earoot,'fibers',filesep,'*.mat']);
fiberscell{1}='Patient-specific DTI-Data';

for fd=2:length(fibd)+1
[~,fn]=fileparts(fibd(fd-1).name);
fiberscell{fd}=fn;
end

set(handles.fiberspopup,'String',fiberscell);

try
    priorselection=find(ismember(fiberscell,stimparams.usefiberset)); % retrieve prior selection of fiberset.
    set(handles.fiberspopup,'Value',priorselection);

catch    % reinitialize using third entry.
    set(handles.fiberspopup,'Value',4);
    
    
end

% Labels:


ll=dir([options.earoot,'templates',filesep,'labeling',filesep,'*.nii']);
for lab=1:length(ll)
    [~,n]=fileparts(ll(lab).name);
    labelcell{lab}=n;
end
% historical part that supported more than one labelatlas:
%labelcell{lab+1}='Use all';

set(handles.labelpopup,'String',labelcell);


try
    priorselection=find(ismember(labelcell,stimparams.labelatlas)); % retrieve prior selection of fiberset.
    if length(priorselection)==1
    set(handles.labelpopup,'Value',priorselection); % set to prior selection
    else % if priorselection was a cell array with more than one entry, set to use all
            set(handles.labelpopup,'Value',lab+1); % set to use all

    end
catch    % reinitialize using third entry.
    set(handles.labelpopup,'Value',1);
    
    
end




% Choose default command line output for ea_stimparams
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);








% UIWAIT makes ea_stimparams wait for user response (see UIRESUME)
% uiwait(handles.figure1);


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



function k0im_Callback(hObject, eventdata, handles)
% hObject    handle to k0im (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k0im as text
%        str2double(get(hObject,'String')) returns contents of k0im as a double


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



function k1im_Callback(hObject, eventdata, handles)
% hObject    handle to k1im (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k1im as text
%        str2double(get(hObject,'String')) returns contents of k1im as a double


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



function k4u_Callback(hObject, eventdata, handles)
% hObject    handle to k4u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k4u as text
%        str2double(get(hObject,'String')) returns contents of k4u as a double


% --- Executes during object creation, after setting all properties.
function k4u_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k4u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function k5u_Callback(hObject, eventdata, handles)
% hObject    handle to k5u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k5u as text
%        str2double(get(hObject,'String')) returns contents of k5u as a double


% --- Executes during object creation, after setting all properties.
function k5u_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k5u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function k6u_Callback(hObject, eventdata, handles)
% hObject    handle to k6u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k6u as text
%        str2double(get(hObject,'String')) returns contents of k6u as a double


% --- Executes during object creation, after setting all properties.
function k6u_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k6u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function k7u_Callback(hObject, eventdata, handles)
% hObject    handle to k7u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k7u as text
%        str2double(get(hObject,'String')) returns contents of k7u as a double


% --- Executes during object creation, after setting all properties.
function k7u_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k7u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function k4im_Callback(hObject, eventdata, handles)
% hObject    handle to k4im (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k4im as text
%        str2double(get(hObject,'String')) returns contents of k4im as a double


% --- Executes during object creation, after setting all properties.
function k4im_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k4im (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function k5im_Callback(hObject, eventdata, handles)
% hObject    handle to k5im (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k5im as text
%        str2double(get(hObject,'String')) returns contents of k5im as a double


% --- Executes during object creation, after setting all properties.
function k5im_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k5im (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function k6im_Callback(hObject, eventdata, handles)
% hObject    handle to k6im (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k6im as text
%        str2double(get(hObject,'String')) returns contents of k6im as a double


% --- Executes during object creation, after setting all properties.
function k6im_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k6im (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function k7im_Callback(hObject, eventdata, handles)
% hObject    handle to k7im (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k7im as text
%        str2double(get(hObject,'String')) returns contents of k7im as a double


% --- Executes during object creation, after setting all properties.
function k7im_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k7im (see GCBO)
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







elstruct=getappdata(gcf,'elstruct');
resultfig=getappdata(gcf,'resultfig');
options=getappdata(gcf,'options');
if isfield(elstruct,'group')
    
    gcnt=ones(length(elstruct(1).groups),1);
    
end

% assign correct .m-file to function.
genvatfunctions=getappdata(gcf,'genvatfunctions');
ea_genvat=eval(['@',genvatfunctions{get(handles.modelselect,'Value')}]);

for el=1:length(elstruct)
    for side=1:length(elstruct.coords_mm)
    if isfield(elstruct,'group') % group analysis, more than one electrode set
        for elin=1:4
            %keyboard
            stimparams(elstruct(el).group,side).U(elin)=str2double(get(eval(['handles.k',num2str(elin-1+((side-1)*options.elspec.numel)),'u']),'String'));
            stimparams(elstruct(el).group,side).Im(elin)=str2double(get(eval(['handles.k',num2str(elin-1+((side-1)*options.elspec.numel)),'im']),'String'));
            
        end
        stimparams(elstruct(el).group,side).group=elstruct(el).group;
        stimparams(elstruct(el).group,side).groupcolors=elstruct(el).groupcolors;
        stimparams(elstruct(el).group,side).groups=elstruct(el).groups;
        flix=elstruct(el).groups;
        
        
        

        
        [stimparams(elstruct(el).group,side).VAT(gcnt(elstruct(el).group)).VAT,radius,volume]=feval(ea_genvat,elstruct(el).coords_mm,stimparams,options);
        stimparams(elstruct(el).group,side).radius=radius;
        stimparams(elstruct(el).group,side).volume=volume;
        
        gcnt(elstruct(el).group)=gcnt(elstruct(el).group)+1;
    else % single patient
        
        for elin=1:4
            %keyboard
            stimparams(1,side).U(elin)=str2double(get(eval(['handles.k',num2str(elin-1+((side-1)*options.elspec.numel)),'u']),'String'));
            stimparams(1,side).Im(elin)=str2double(get(eval(['handles.k',num2str(elin-1+((side-1)*options.elspec.numel)),'im']),'String'));
            
        end
        
        [stimparams(1,side).VAT(el).VAT,radius,volume]=feval(ea_genvat,elstruct(el).coords_mm,stimparams,side,options);
           stimparams(1,side).radius=radius;
        stimparams(1,side).volume=volume;
        flix=1;
    end
    end
end

for group=1:length(stimparams)
    fiberscell=get(handles.fiberspopup,'String');
    stimparams(group).usefiberset=fiberscell{get(handles.fiberspopup,'Value')};
    labelcell=get(handles.labelpopup,'String');
    stimparams(group).labelatlas=labelcell(get(handles.labelpopup,'Value'));
    stimparams(group).showfibers=(get(handles.showfibs,'Value') == get(handles.showfibs,'Max'));
    stimparams(group).showconnectivities=(get(handles.showconns,'Value') == get(handles.showconns,'Max'));
    stimparams(group).fiberthresh=str2double(get(handles.fiberthresh,'String'));

    % historical part that supported more than one labelatlas
%     if strcmp(stimparams(group).labelatlas{1},'Use all')
%         ll=dir([options.earoot,'templates',filesep,'labeling',filesep,'*.nii']);
%         for lab=1:length(ll)
%             [~,n]=fileparts(ll(lab).name);
%             stimparams(group).labelatlas{lab}=n;
%         end
%     end
end









PL=getappdata(resultfig,'PL');
for group=1:length(PL)
        deletePL(PL(group));
    
end
clear PL

figtitle=get(resultfig,'Name');
set(resultfig,'Name',[figtitle,'...building...']);

for group=flix
    setappdata(resultfig,'stimparams',stimparams(group,:));
    
    ea_showfibres_volume(resultfig,options);
    copyfile([options.root,options.patientname,filesep,'ea_stats.mat'],[options.root,options.patientname,filesep,'ea_stats_group_',num2str(group),'.mat']);
    try
    copyfile([options.root,options.patientname,filesep,'ea_pm.nii'],[options.root,options.patientname,filesep,'ea_pm_group_',num2str(group),'.nii']);
    end
    try
    PL(group)=getappdata(resultfig,'PL');
    catch
        keyboard
    end
end
setappdata(resultfig,'PL',PL);
set(resultfig,'Name',figtitle);


function deletePL(PL)
if verLessThan('matlab','8.5') % ML <2014a support
    
    
    for p=1:length(PL)
        
        
        if isfield(PL(p),'vatsurfs')
            delete(PL(p).vatsurfs(logical(PL(p).vatsurfs)));
        end
        if isfield(PL(p),'fib_plots')
            if isfield(PL(p).fib_plots,'fibs')
                delete(PL(p).fib_plots.fibs(logical(PL(p).fib_plots.fibs)));
            end
            
            if isfield(PL(p).fib_plots,'dcfibs')
                todelete=PL(p).fib_plots.dcfibs((PL(p).fib_plots.dcfibs(:)>0));
                delete(todelete(:));
                
            end
        end
        if isfield(PL(p),'regionsurfs')
            todelete=PL(p).regionsurfs(logical(PL(p).regionsurfs));
            delete(todelete(:));
        end
        if isfield(PL(p),'conlabels')
            todelete=PL(p).conlabels(logical(PL(p).conlabels));
            delete(todelete(:));
        end
        if isfield(PL(p),'ht')
            delete(PL(p).ht);
        end
    end
    
    
else
    for p=1:length(PL) 
        if isfield(PL(p),'vatsurfs')
            delete(PL(p).vatsurfs);
        end
        if isfield(PL(p),'fib_plots')
            if isfield(PL(p).fib_plots,'fibs')
                delete(PL(p).fib_plots.fibs);
            end
            
            if isfield(PL(p).fib_plots,'dcfibs')
                delete(PL(p).fib_plots.dcfibs);
            end
        end
        if isfield(PL(p),'regionsurfs')
            delete(PL(p).regionsurfs);
        end
        if isfield(PL(p),'conlabels')
            delete(PL(p).conlabels);
        end
        if isfield(PL(p),'ht')
            delete(PL(p).ht);
        end
    end
    
end



function fiberthresh_Callback(hObject, eventdata, handles)
% hObject    handle to fiberthresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fiberthresh as text
%        str2double(get(hObject,'String')) returns contents of fiberthresh as a double


% --- Executes during object creation, after setting all properties.
function fiberthresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fiberthresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in showfibs.
function showfibs_Callback(hObject, eventdata, handles)
% hObject    handle to showfibs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showfibs


% --- Executes on button press in showconns.
function showconns_Callback(hObject, eventdata, handles)
% hObject    handle to showconns (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showconns


% --- Executes on selection change in fiberspopup.
function fiberspopup_Callback(hObject, eventdata, handles)
% hObject    handle to fiberspopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns fiberspopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from fiberspopup


% --- Executes during object creation, after setting all properties.
function fiberspopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fiberspopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
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
