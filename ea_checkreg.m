function varargout = ea_checkreg(varargin)
% EA_CHECKREG MATLAB code for ea_checkreg.fig
%      EA_CHECKREG, by itself, creates a new EA_CHECKREG or raises the existing
%      singleton*.
%
%      H = EA_CHECKREG returns the handle to a new EA_CHECKREG or the handle to
%      the existing singleton*.
%
%      EA_CHECKREG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EA_CHECKREG.M with the given input arguments.
%
%      EA_CHECKREG('Property','Value',...) creates a new EA_CHECKREG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ea_checkreg_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ea_checkreg_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ea_checkreg

% Last Modified by GUIDE v2.5 03-Sep-2021 19:47:16

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @ea_checkreg_OpeningFcn, ...
    'gui_OutputFcn',  @ea_checkreg_OutputFcn, ...
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


% --- Executes just before ea_checkreg is made visible.
function ea_checkreg_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ea_checkreg (see VARARGIN)

options = varargin{1};

handles.leadfigure.Name = [options.subj.subjId, ': Check Coregistration'];
handles.patientname.String = options.subj.subjId;

set(handles.normsettings, 'Visible', 'off');

% Get coregistered pre-op images (except for the anchor image)
preopCoregImages = struct2cell(options.subj.coreg.anat.preop);
preopCoregImages = preopCoregImages(2:end);

% Get coregistered post-op images
if strcmp(options.subj.postopModality, 'CT')
    postopCoregImages = options.subj.coreg.anat.postop.tonemapCT;
elseif strcmp(options.subj.postopModality, 'MRI')
    postopCoregImages = struct2cell(options.subj.coreg.anat.postop);
else
    postopCoregImages = [];
end

% Get normalized pre-op anchor image
preopNormImage = options.subj.norm.anat.preop.(options.subj.AnchorModality);

% Get brain shift corrected image
if isfield(options.subj, 'brainshift')
    brainshiftImage = options.subj.brainshift.anat.scrf;
else
    brainshiftImage = [];
end

% List pf images for checkreg
checkregImages = [preopCoregImages; postopCoregImages; preopNormImage; brainshiftImage];
checkregImages = checkregImages(cellfun(@(f) ea_reglocked(options, f)~=1 & isfile(f), checkregImages));

% fMRI
restfiles = dir([options.root,options.patientname,filesep,options.prefs.rest_searchstring]);
options.prefs.n_rest = numel(restfiles);
b0restanchor = cell(length(checkregImages),1);
for irest = 1:options.prefs.n_rest
    % Set filenames for this iteration
    if exist(['r',ea_stripext(restfiles(irest).name),'_t1'],'file')
        if ~ea_reglocked(options,['r',ea_stripext(restfiles(irest).name),'_t1'])
            checkregImages = [checkregImages;{['r',ea_stripext(restfiles(irest).name),'_t1']}];
            b0restanchor{length(checkregImages)} = ['mean',restfiles(irest).name];
        end
    end
end

% b0 image
if exist([ea_stripext(options.prefs.b0),'_t1'],'file')
    if ~ea_reglocked(options,[ea_stripext(options.prefs.b0),'_t1'])
        checkregImages = [checkregImages;{[ea_stripext(options.prefs.b0),'_t1']}];
        b0restanchor{length(checkregImages)} = [options.prefs.b0];
    end
end

if isempty(checkregImages)
    evalin('base','checkregempty=1;');
    close(handles.leadfigure)
    return
else
    evalin('base','checkregempty=0;');
end

options.overwriteapproved = 0;

%set(handles.previous,'visible','off'); set(handles.next,'visible','off');
setappdata(handles.leadfigure, 'checkregImages', checkregImages)
setappdata(handles.leadfigure, 'b0restanchor', b0restanchor)
setappdata(handles.leadfigure, 'activevolume', 1);
setappdata(handles.leadfigure, 'options', options);

set(handles.checkatl, 'Visible', 'off');

ea_mrcview(handles);

if isvalid(hObject)
    % Choose default command line output for ea_checkreg
    handles.output = hObject;

    % Update handles structure
    guidata(hObject, handles);
end
% UIWAIT makes ea_checkreg wait for user response (see UIRESUME)


function ea_mrcview(handles)

options = getappdata(handles.leadfigure, 'options');
checkregImages = getappdata(handles.leadfigure, 'checkregImages');
activevolume = getappdata(handles.leadfigure, 'activevolume');
b0restanchor = getappdata(handles.leadfigure, 'b0restanchor');

if activevolume==length(checkregImages)
    set(handles.disapprovebutn, 'String', 'Disapprove & Close');
    set(handles.approvebutn, 'String', 'Approve & Close');
else
    set(handles.approvebutn, 'String', 'Approve & Next >>');
    set(handles.disapprovebutn, 'String', 'Disapprove & Next >>');
end

currvol = checkregImages{activevolume};

% Brain shift corrected image
if isfield(options.subj, 'brainshift') && strcmp(currvol, options.subj.brainshift.anat.scrf)
    ea_subcorticalrefine(options);
    close(handles.leadfigure);
    return
end

if strcmp(currvol, options.subj.norm.anat.preop.(options.subj.AnchorModality))
    handles.checkatl.Visible='on';
    anchor = options.primarytemplate;

    ea_init_normpopup(handles, options.prefs.normalize.default, 'coregmrmethod');

    if isfile(options.subj.norm.log.method)
        json = loadjson(options.subj.norm.log.method);
    end

    json.approval = 0;
    ea_mkdir(fileparts(options.subj.coreg.log.method));
    savejson('', json, options.subj.norm.log.method);

    checkregFig = options.subj.norm.checkreg.preop.(options.subj.AnchorModality);

    set(handles.anchortxt, 'String', 'Template (red wires):');
    set(handles.coregresultstxt, 'String', 'Normalization results');
    set(handles.normsettings, 'Visible', 'on');
    set(handles.recomputebutn, 'String', '(Re-) compute normalization using...');
    set(handles.coregmrmethod, 'TooltipString', 'Choose a normalization method');
    set(handles.leadfigure, 'Name', [options.subj.subjId, ': Check Normalization']);
else
    handles.checkatl.Visible = 'off';
    anchor = options.subj.AnchorModality;
    set(handles.anchortxt, 'String', 'Anchor modality (red wires):');
    set(handles.coregresultstxt, 'String', 'Coregistration results');
    set(handles.leadfigure, 'Name', [options.subj.subjId, ': Check Coregistration']);

    ea_mkdir(fileparts(options.subj.coreg.log.method));

    if strcmp(options.subj.postopModality, 'CT') && strcmp(currvol, options.subj.coreg.anat.postop.tonemapCT)
        ea_init_coregctpopup(handles, options.prefs.ctcoreg.default, 'coregmrmethod');

        if isfile(options.subj.coreg.log.method)
            json = loadjson(options.subj.coreg.log.method);
        end

        json.approval.CT = 0;
        savejson('', json, options.subj.coreg.log.method);

        checkregFig = options.subj.coreg.checkreg.postop.tonemapCT;
    else % MR
        ea_init_coregmrpopup(handles, options.prefs.mrcoreg.default);

        if isfile(options.subj.coreg.log.method)
            json = loadjson(options.subj.coreg.log.method);
        end

        % Extract image modality
        modality = ea_getmodality(currvol);

        json.approval.(modality) = 0;
        savejson('', json, options.subj.coreg.log.method);

        if contains(currvol, 'ses-preop')
            checkregFig = options.subj.coreg.checkreg.preop.(modality);
        else
            checkregFig = options.subj.coreg.checkreg.postop.(modality);
        end
    end

    set(handles.normsettings, 'Visible', 'off');
    set(handles.recomputebutn, 'String', '(Re-) compute coregistration using...');
    set(handles.coregmrmethod, 'TooltipString', 'Choose a coregistration method');
end

% show result:
if ~isempty(b0restanchor) && ~isempty(b0restanchor{activevolume}) % rest or b0 registration
    set(handles.substitute, 'Visible', 'on');
    set(handles.substitute, 'String', ea_getsubstitutes(options));
    checkregFig = [directory, 'checkreg', filesep,ea_stripext(currvol),'2',strrep(ea_stripext(b0restanchor{activevolume}),'mean','r'),'_',method,'.png'];
    set(handles.anchormod,'String',ea_stripext(b0restanchor{activevolume}));
else % normal anatomical 2 anatomical registration
    set(handles.substitute, 'Visible', 'off');
    set(handles.anchormod, 'String', anchor);
end

set(handles.imgfn, 'Visible', 'on');
set(handles.imgfn, 'String', strrep(checkregFig, [fileparts(options.subj.subjDir), filesep], ''));
set(handles.imgfn, 'TooltipString', strrep(checkregFig, [fileparts(options.subj.subjDir), filesep], ''));

switch currvol
    case options.subj.norm.anat.preop.(options.subj.AnchorModality)
        anchorPath = [ea_space, options.primarytemplate];
    otherwise
        if ~isempty(b0restanchor) && ~isempty(b0restanchor{activevolume}) % rest or b0 registration
            anchorPath = ea_stripext(b0restanchor{activevolume});
        else
            anchorPath = options.subj.coreg.anat.preop.(options.subj.AnchorModality);
        end
end

if ~isfile(checkregFig)
    ea_gencheckregpair(currvol, anchorPath, checkregFig);

    if ~isfile(checkregFig)
        checkregFig = fullfile(ea_getearoot,'helpers','gui','coreg_msg.png');
        set(handles.imgfn, 'String', '');
        set(handles.imgfn, 'Visible', 'off');
    end
end

im = imread(checkregFig);
set(0,'CurrentFigure',handles.leadfigure);
set(handles.leadfigure,'CurrentAxes',handles.standardax);

imagesc(im);
axis off
axis equal

% textfields:
set(handles.depvolume, 'String', [ea_stripext(currvol),'.nii']);
set(handles.depvolume, 'Tooltip', [ea_stripext(currvol),'.nii']);


function [pretras]=ea_getsubstitutes(options)

[~, checkregImages] = ea_assignpretra(options);
for fi=1:length(checkregImages)
    if fi==1
        pretras{fi} = ['Use ',checkregImages{fi}, ' (default)'];
    else
        pretras{fi} = ['Substitute moving file with ',checkregImages{fi}];
    end
end


% --- Outputs from this function are returned to the command line.
function varargout = ea_checkreg_OutputFcn(hObject, eventdata, handles)
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


% --- Executes on button press in recomputebutn.
function recomputebutn_Callback(hObject, eventdata, handles)
% hObject    handle to recomputebutn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ea_busyaction('on', handles.leadfigure, 'coreg');

options = getappdata(handles.leadfigure, 'options');
options.overwriteapproved = 1;

checkregImages = getappdata(handles.leadfigure, 'checkregImages');
activevolume = getappdata(handles.leadfigure, 'activevolume');
currvol = checkregImages{activevolume};

anchorImage = options.subj.coreg.anat.preop.(options.subj.AnchorModality);

if strcmp(currvol, options.subj.norm.anat.preop.(options.subj.AnchorModality))
    ea_delete(struct2cell(options.subj.norm.anat.preop));
    if isfield(options.subj.norm.anat, 'postop')
        ea_delete(struct2cell(options.subj.norm.anat.postop));
    end

    % Get normalization method
    options.normalize.method = handles.coregmrmethod.String{handles.coregmrmethod.Value};

    % Re-run normlization
    ea_normalize(options);

    % Create checkreg fig
    fprintf('\nRegenerating checkreg figure...\n\n');
    ea_gencheckregfigs(options, 'norm');

    % Disapprove after recompute
    json = loadjson(options.subj.norm.log.method);
    json.approval = 0;
    savejson('', json, options.subj.norm.log.method);

    if isfolder(options.subj.brainshiftDir)
        ea_cprintf('CmdWinWarnings', 'Normalization has been rerun. Please also rerun brain shift correction!\n');
    end

elseif strcmp(options.subj.postopModality, 'CT') && strcmp(currvol, options.subj.coreg.anat.postop.tonemapCT)
    % Get CT coregistration method
    options.coregct.method = handles.coregmrmethod.String{handles.coregmrmethod.Value};

    % Run CT coregistration
    ea_coregpostopct(options);

    % Create checkreg fig
    fprintf('\nRegenerating checkreg figure...\n\n');
    ea_gencheckregpair(currvol, anchorImage, options.subj.coreg.checkreg.postop.tonemapCT);

    % Disapprove after recompute
    json = loadjson(options.subj.coreg.log.method);
    json.approval.CT = 0;
    savejson('', json, options.subj.coreg.log.method);

    if isfolder(options.subj.brainshiftDir)
        ea_cprintf('CmdWinWarnings', 'CT coregistration has been rerun. Please also rerun brain shift correction!\n');
    end

elseif strcmp(ea_stripext(options.prefs.fa2anat), 'FA') % FA
    options.coregmr.method=get(handles.coregmrmethod,'String');
    options.coregmr.method=options.coregmr.method{get(handles.coregmrmethod,'Value')};
    ea_backuprestore([directory,options.prefs.fa]);
    ea_coregimages(options,[directory,options.prefs.fa],[directory,anchor],[directory,checkregImages{activevolume}],{},0);

else % MR
    options.coregmr.method = handles.coregmrmethod.String{handles.coregmrmethod.Value};

    b0restanchor = getappdata(handles.leadfigure, 'b0restanchor');
    if ~isempty(b0restanchor{activevolume})  % b0 or rest files
        thisrest=strrep(ea_stripext(b0restanchor{activevolume}),'mean','r');

        ea_delete([directory,thisrest,'2',ea_stripext(anchor),'_',ea_matext(options.coregmr.method)]);
        ea_delete([directory,ea_stripext(anchor),'2',thisrest,'_',ea_matext(options.coregmr.method)]);

        substitute=get(handles.substitute,'Value');
        [~,pf]=ea_assignpretra(options);
        useasanchor=pf{substitute};

        % in following line correct that useasanchor is the *moving*
        % image (since we're going from anchor to rest/b0.
        ea_coregimages(options,[directory,useasanchor],[directory,b0restanchor{activevolume}],[directory,checkregImages{activevolume}],{},1);
        if ~isequal([directory,ea_stripext(b0restanchor{activevolume}),'2',ea_stripext(useasanchor),'_',ea_matext(options.coregmr.method)],...
                [directory,thisrest,'2',ea_stripext(anchor),'_',ea_matext(options.coregmr.method)])
            movefile([directory,ea_stripext(b0restanchor{activevolume}),'2',ea_stripext(useasanchor),'_',ea_matext(options.coregmr.method)],...
                [directory,thisrest,'2',ea_stripext(anchor),'_',ea_matext(options.coregmr.method)]);
        end
        if ~isequal([directory,ea_stripext(useasanchor),'2',ea_stripext(b0restanchor{activevolume}),'_',ea_matext(options.coregmr.method)],...
                [directory,ea_stripext(anchor),'2',thisrest,'_',ea_matext(options.coregmr.method)])
            movefile([directory,ea_stripext(useasanchor),'2',ea_stripext(b0restanchor{activevolume}),'_',ea_matext(options.coregmr.method)],...
                [directory,ea_stripext(anchor),'2',thisrest,'_',ea_matext(options.coregmr.method)]);
        end
        ea_cleandownstream(directory,thisrest);
    else  % other images
        session = regexp(currvol, '(?<=_ses-)(preop|postop)', 'match', 'once');
        modality = ea_getmodality(currvol);
        preprocImage = options.subj.preproc.anat.(session).(modality);
        ea_coregimages(options, preprocImage, anchorImage, currvol, {}, 0);

        % Create checkreg fig
        fprintf('\nRegenerating checkreg figure...\n\n');
        ea_gencheckregpair(currvol, anchorImage, options.subj.coreg.checkreg.(session).(modality));

        % Disapprove after recompute
        json = loadjson(options.subj.coreg.log.method);
        json.approval.(modality) = 0;
        savejson('', json, options.subj.coreg.log.method);

        if strcmp(session, 'postop') && isfolder(options.subj.brainshiftDir)
            ea_cprintf('CmdWinWarnings', 'Postop MR coregistration has been rerun. Please also rerun brain shift correction!\n');
        elseif strcmp(session, 'preop')
            ea_cprintf('CmdWinWarnings', 'Preop MR coregistration has been rerun. Please also rerun normalization!\n');
        end
    end
end

ea_mrcview(handles);
title = get(handles.leadfigure, 'Name');
ea_chirp(options);
ea_busyaction('off', handles.leadfigure, 'coreg');
set(handles.leadfigure, 'Name', title);


function ea_cleandownstream(directory, thisrest)
% cleanup fibertracking mask
ea_delete([directory,'trackingmask.nii']);
ea_delete([directory,'ttrackingmask.nii']);

% cleanup /templates/labelings (these need to be recalculated)
ea_delete([directory,'templates',filesep,'labeling',filesep,thisrest,'*.nii']);
parcdirs=dir([directory,'connectomics',filesep]);

% cleanup /connectomics results (these need to be recalculated):
for pd=1:length(parcdirs)
    if ~strcmp(parcdirs(pd).name(1),'.')
        ea_delete([directory,'connectomics',filesep,parcdirs(pd).name,filesep,thisrest(2:end),'*.*']);
    end
end

for tn=0:1
    options.native=tn;
    stimdirs=dir([directory,'stimulations',filesep,ea_nt(options)]);
    % cleanup /connectomics results (these need to be recalculated):
    for pd=1:length(stimdirs)
        if ~strcmp(stimdirs(pd).name(1),'.')
            connfolders=dir([directory,'stimulations',ea_nt(options),filesep,stimdirs(pd).name]);
            for connfolder=1:length(connfolders)
                if ~strcmp(connfolders(connfolder).name(1),'.')

                    if connfolders(connfolder).isdir && contains(connfolders(connfolder).name,thisrest(2:end))
                        rmdir([directory,'stimulations',filesep,ea_nt(options),stimdirs(pd).name,filesep,connfolders(connfolder).name],'s');
                    end

                    ea_delete([directory,'stimulations',filesep,ea_nt(options),stimdirs(pd).name,filesep,'*',thisrest(2:end),'*.nii']);
                end
            end
        end
    end
end


function ext=ea_matext(method)

switch upper(method)
    case 'SPM'
        ext='spm.mat';
    case 'FSL FLIRT'
        ext='flirt1.mat';
    case 'FSL BBR'
        ext='flirtbbr.mat';
    case 'ANTS'
        ext='ants1.mat';
    case 'BRAINSFIT'
        ext='brainsfit.h5';
end


% --- Executes on button press in approvebutn.
function approvebutn_Callback(hObject, eventdata, handles)
% hObject    handle to approvebutn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ea_busyaction('on', handles.leadfigure, 'coreg');

options = getappdata(handles.leadfigure, 'options');

checkregImages = getappdata(handles.leadfigure, 'checkregImages');
activevolume = getappdata(handles.leadfigure, 'activevolume');
currvol = checkregImages{activevolume};

% Get coregistered pre-op images (except for the anchor image)
preopCoregImages = struct2cell(options.subj.coreg.anat.preop);
preopCoregImages = preopCoregImages(2:end);

if strcmp(currvol, options.subj.norm.anat.preop.(options.subj.AnchorModality))
    json = loadjson(options.subj.norm.log.method);
    if all(cellfun(@(f) ea_reglocked(options, f), preopCoregImages))
        json.approval = 1;
    else
        warning('off', 'backtrace');
        warning('You approved normalization before all preoperative co-registrations were approved. Lead-DBS will still override / redo normalization if applying a multispectral method.');
        warning('on', 'backtrace');
        json.approval = 0.5;
    end

    savejson('', json, options.subj.norm.log.method);

elseif strcmp(options.subj.postopModality, 'CT') && strcmp(currvol, options.subj.coreg.anat.postop.tonemapCT)
    json = loadjson(options.subj.coreg.log.method);
    json.approval.CT = 1;
    savejson('', json, options.subj.coreg.log.method);

else
    json = loadjson(options.subj.coreg.log.method);
    modality = ea_getmodality(currvol);

    if ismember(currvol, preopCoregImages)
        try % Field might not exist.
            coregWasApproved = json.approval.(modality);
        catch
            coregWasApproved = 0;
        end

        json.approval.(modality) = 1;
        savejson('', json, options.subj.coreg.log.method);

        if ~coregWasApproved && isfile(options.subj.norm.log.method)
            json = loadjson(options.subj.norm.log.method);
            if isfield(json, 'approval') && json.approval==1
                % In this situation we had the normalization approved before
                % all coregistrations were approved. This could lead to suboptimal
                % normalization when a multispectral normalization method is used.
                % Thus we set the normalization approval flag to 0.5. This way,
                % LeadDBS will redo the normalization when using a multispectral
                % method.
                warning('off', 'backtrace');
                warning('Normalization had been approved before all pre-op coregistrations were approved. Lead-DBS will redo normalization when using a multispectral method.');
                warning('on', 'backtrace');
                json.approval = 0.5;
                savejson('', json, options.subj.norm.log.method);
            end
        end
    else
        json.approval.(modality) = 1;
        savejson('', json, options.subj.coreg.log.method);
    end
end

if ismac
    system(['xattr -wx com.apple.FinderInfo "0000000000000000000400000000000000000000000000000000000000000000" ',currvol]);
end

if activevolume == length(checkregImages)
    close(handles.leadfigure);
    return
elseif activevolume==length(checkregImages)-1 ...
        && isfield(options.subj, 'brainshift') ...
        && strcmp(checkregImages{end}, options.subj.brainshift.anat.scrf)
    close(handles.leadfigure);
    ea_subcorticalrefine(options);
    return
else
    activevolume = activevolume+1;
end

setappdata(handles.leadfigure, 'activevolume', activevolume);

ea_mrcview(handles);
title = get(handles.leadfigure, 'Name');
ea_busyaction('off', handles.leadfigure, 'coreg');
set(handles.leadfigure, 'Name', title);


% --- Executes on selection change in coregmrmethod.
function coregmrmethod_Callback(hObject, eventdata, handles)
% hObject    handle to coregmrmethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns coregmrmethod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from coregmrmethod

options=getappdata(handles.leadfigure,'options');

checkregImages = getappdata(handles.leadfigure,'checkregImages');
activevolume = getappdata(handles.leadfigure,'activevolume');
currvol = checkregImages{activevolume};

% Set callback and Enable status of normalization setting button
if strcmp(currvol, options.subj.norm.anat.preop.(options.subj.AnchorModality))
    ea_normsettings(handles, 'coregmrmethod');
end


% --- Executes during object creation, after setting all properties.
function coregmrmethod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to coregmrmethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in openviewer.
function openviewer_Callback(hObject, eventdata, handles)
% hObject    handle to openviewer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
options = getappdata(handles.leadfigure,'options');
checkregImages = getappdata(handles.leadfigure,'checkregImages');
activevolume = getappdata(handles.leadfigure,'activevolume');
b0restanchor = getappdata(handles.leadfigure,'b0restanchor');

currvol = checkregImages{activevolume};
if strcmp(currvol, options.subj.norm.anat.preop.(options.subj.AnchorModality))
	ea_show_normalization(options);
else
    options.moving = checkregImages{activevolume};
    if ~isempty(b0restanchor{activevolume})
        options.fixed = b0restanchor{activevolume};
        options.tag = [checkregImages{activevolume},' & ',b0restanchor{activevolume}];
    else
        options.fixed = options.subj.coreg.anat.preop.(options.subj.AnchorModality);
        options.tag = [ea_getmodality(options.moving), ' & ', ea_getmodality(options.fixed)];
    end

    ea_show_coregistration(options);
end


% --- Executes on button press in normsettings.
function normsettings_Callback(hObject, eventdata, handles)
% hObject    handle to normsettings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
normsettingsfunc = getappdata(handles.normsettings,'normsettingsfunc');
feval(normsettingsfunc, handles);


% --- Executes on button press in disapprovebutn.
function disapprovebutn_Callback(hObject, eventdata, handles)
% hObject    handle to disapprovebutn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% hObject    handle to approvebutn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ea_busyaction('on',handles.leadfigure,'coreg');

options = getappdata(handles.leadfigure,'options');

checkregImages = getappdata(handles.leadfigure,'checkregImages');
activevolume = getappdata(handles.leadfigure,'activevolume');
currvol = checkregImages{activevolume};

% Get coregistered pre-op images (except for the anchor image)
preopCoregImages = struct2cell(options.subj.coreg.anat.preop);
preopCoregImages = preopCoregImages(2:end);

if strcmp(currvol, options.subj.norm.anat.preop.(options.subj.AnchorModality))
    json = loadjson(options.subj.norm.log.method);
    if all(cellfun(@(f) ea_reglocked(options, f), preopCoregImages))
        json.approval = 0.5;
    else
        json.approval = 0;
    end

    savejson('', json, options.subj.norm.log.method);

elseif strcmp(options.subj.postopModality, 'CT') && strcmp(currvol, options.subj.coreg.anat.postop.tonemapCT)
    json = loadjson(options.subj.coreg.log.method);
    json.approval.CT = 0;
    savejson('', json, options.subj.coreg.log.method);

else
    json = loadjson(options.subj.coreg.log.method);
    modality = ea_getmodality(currvol);

    if ismember(currvol, preopCoregImages)
        try % Field might not exist.
            coregWasApproved = json.approval.(modality);
        catch
            coregWasApproved = 0;
        end

        json.approval.(modality) = 0;
        savejson('', json, options.subj.coreg.log.method);

        if coregWasApproved && isfile(options.subj.norm.log.method)
            json = loadjson(options.subj.norm.log.method);
            if isfield(json, 'approval') && json.approval==1
                % In this situation we had the normalization approved but the
                % coregistration is disapproved. This could lead to suboptimal
                % normalization when a multispectral normalization method is
                % used. Thus we set the normalization approval flag to 0.5.
                % This way, LeadDBS will redo the normalization when using
                % a multispectral method.
                warning('off', 'backtrace');
                warning('Normalization had been approved but pre-op coregistration is disapproved now. Lead-DBS will redo normalization when using a multispectral method.');
                warning('on', 'backtrace');
                json.approval = 0.5;
                savejson('', json, options.subj.norm.log.method);
            end
        end
    else
        json.approval.(modality) = 0;
        savejson('', json, options.subj.coreg.log.method);
    end
end

if ismac
    system(['xattr -wx com.apple.FinderInfo "0000000000000000000C00000000000000000000000000000000000000000000" ', currvol]);
end

b0restanchor = getappdata(handles.leadfigure,'b0restanchor');
if ~isempty(b0restanchor{activevolume})
    thisrest = strrep(ea_stripext(b0restanchor{activevolume}),'mean','r');
    ea_cleandownstream(directory,thisrest)
end

checkregImages = getappdata(handles.leadfigure,'checkregImages');
activevolume = getappdata(handles.leadfigure,'activevolume');

if activevolume == length(checkregImages)
    close(handles.leadfigure);
    return
else
    activevolume = activevolume+1;
end

setappdata(handles.leadfigure, 'activevolume', activevolume);
ea_mrcview(handles);
try
    title = get(handles.leadfigure, 'Name');
    ea_busyaction('off',handles.leadfigure,'coreg');
    set(handles.leadfigure, 'Name', title);
end


% --- Executes on button press in back.
function back_Callback(hObject, eventdata, handles)
% hObject    handle to back (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ea_busyaction('on', handles.leadfigure, 'coreg');

activevolume = getappdata(handles.leadfigure, 'activevolume');

if activevolume==1
    ea_busyaction('off',handles.leadfigure, 'coreg');
    return
else
    activevolume = activevolume-1;
    setappdata(handles.leadfigure, 'activevolume', activevolume);
end

ea_mrcview(handles);
title = get(handles.leadfigure, 'Name');
ea_busyaction('off', handles.leadfigure, 'coreg');
set(handles.leadfigure, 'Name', title);


% --- Executes on button press in refreshview.
function refreshview_Callback(hObject, eventdata, handles)
% hObject    handle to refreshview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ea_busyaction('on',handles.leadfigure,'coreg');

options = getappdata(handles.leadfigure, 'options');

checkregImages = getappdata(handles.leadfigure, 'checkregImages');
activevolume = getappdata(handles.leadfigure, 'activevolume');
currvol = checkregImages{activevolume};

anchorImage = options.subj.coreg.anat.preop.(options.subj.AnchorModality);

if strcmp(currvol, options.subj.norm.anat.preop.(options.subj.AnchorModality))
     % Create checkreg fig
    fprintf('\nRegenerating checkreg figure...\n\n');
    ea_gencheckregfigs(options, 'norm');

elseif strcmp(options.subj.postopModality, 'CT') && strcmp(currvol, options.subj.coreg.anat.postop.tonemapCT)
    % Create checkreg fig
    fprintf('\nRegenerating checkreg figure...\n\n');
    ea_gencheckregpair(currvol, anchorImage, options.subj.coreg.checkreg.postop.tonemapCT);
else
    b0restanchor = getappdata(handles.leadfigure,'b0restanchor');
    if ~isempty(b0restanchor{activevolume})
       anchorImage = b0restanchor{activevolume};
    end

    session = regexp(currvol, '(?<=_ses-)(preop|postop)', 'match', 'once');
    modality = ea_getmodality(currvol);

    % Create checkreg fig
    fprintf('\nRegenerating checkreg figure...\n\n');
    ea_gencheckregpair(currvol, anchorImage, options.subj.coreg.checkreg.(session).(modality));
end

ea_mrcview(handles);
title = get(handles.leadfigure, 'Name');
ea_busyaction('off', handles.leadfigure, 'coreg');
set(handles.leadfigure, 'Name', title);


% --- Executes on button press in checkatl.
function checkatl_Callback(hObject, eventdata, handles)
% hObject    handle to checkatl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
options=getappdata(handles.leadfigure,'options');

ea_checkstructures(options);


% --- Executes on button press in openpatientdir.
function openpatientdir_Callback(hObject, eventdata, handles)
% hObject    handle to openpatientdir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
options = getappdata(handles.leadfigure, 'options');
ea_opendir(options.subj.subjDir);


% --- Executes on selection change in substitute.
function substitute_Callback(hObject, eventdata, handles)
% hObject    handle to substitute (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns substitute contents as cell array
%        contents{get(hObject,'Value')} returns selected item from substitute


% --- Executes during object creation, after setting all properties.
function substitute_CreateFcn(hObject, eventdata, handles)
% hObject    handle to substitute (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
