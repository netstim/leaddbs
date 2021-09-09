function varargout = ea_checkcoreg(varargin)
% EA_CHECKCOREG MATLAB code for ea_checkcoreg.fig
%      EA_CHECKCOREG, by itself, creates a new EA_CHECKCOREG or raises the existing
%      singleton*.
%
%      H = EA_CHECKCOREG returns the handle to a new EA_CHECKCOREG or the handle to
%      the existing singleton*.
%
%      EA_CHECKCOREG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EA_CHECKCOREG.M with the given input arguments.
%
%      EA_CHECKCOREG('Property','Value',...) creates a new EA_CHECKCOREG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ea_checkcoreg_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ea_checkcoreg_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ea_checkcoreg

% Last Modified by GUIDE v2.5 03-Sep-2021 19:47:16

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @ea_checkcoreg_OpeningFcn, ...
    'gui_OutputFcn',  @ea_checkcoreg_OutputFcn, ...
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


% --- Executes just before ea_checkcoreg is made visible.
function ea_checkcoreg_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ea_checkcoreg (see VARARGIN)

options=varargin{1};
%ea_init_coregmrpopup(handles,1);
set(handles.leadfigure,'Name',[options.patientname, ': Check Coregistration']);

directory=[options.root,options.patientname,filesep];

setappdata(handles.leadfigure,'options',options);
setappdata(handles.leadfigure,'directory',directory);

[~, patientname]=fileparts(fileparts(directory));
handles.patientname.String=patientname;

set(handles.leadfigure,'Name',[patientname, ': MR-Coregistration']);

presentfiles=ea_getall_coregcheck(options);
anchor=presentfiles{1};
presentfiles(1)=[];

set(handles.normsettings,'Visible','off');
if exist([directory,options.prefs.gprenii],'file') && ~ea_reglocked(options,options.prefs.gprenii)
    presentfiles=[presentfiles;{[directory,options.prefs.gprenii]}];
end

% if exist([directory,options.prefs.gtranii],'file') && ~ea_reglocked(options,options.prefs.gtranii)
%     presentfiles=[presentfiles;{[directory,options.prefs.gtranii]}];
% end
%
% if exist([directory,options.prefs.tp_gctnii],'file') && ~ea_reglocked(options,options.prefs.tp_gctnii)
%     presentfiles=[presentfiles;{[directory,options.prefs.tp_gctnii]}];
% end

% add coregchecks for b0 and rest:
% get files with rs-fMRI data
restfiles = dir([options.root,options.patientname,filesep,options.prefs.rest_searchstring]);

% get number of files with rs-fMRI data
options.prefs.n_rest = numel(restfiles);

b0restanchor=cell(length(presentfiles),1);
for irest = 1:options.prefs.n_rest
    % set filenames for this iteration
    if exist([directory,'r',ea_stripext(restfiles(irest).name),'_',anchor],'file')
        if ~ea_reglocked(options,['r',ea_stripext(restfiles(irest).name),'_',anchor])
            presentfiles=[presentfiles;{['r',ea_stripext(restfiles(irest).name),'_',anchor]}];
            b0restanchor{length(presentfiles)} = ['mean',restfiles(irest).name];
        end
    end
end

% add b0:
if exist([directory,ea_stripext(options.prefs.b0),'_',anchor],'file')
    if ~ea_reglocked(options,[directory,ea_stripext(options.prefs.b0),'_',anchor])
        presentfiles=[presentfiles;{[ea_stripext(options.prefs.b0),'_',anchor]}];
        b0restanchor{length(presentfiles)} = [options.prefs.b0];
    end
end

if exist([directory,'scrf',filesep,'scrf_instore_converted.mat'],'file')
    if ~ea_reglocked(options,'brainshift')
        presentfiles=[presentfiles; {'brainshift'}];
    end
end

if isempty(presentfiles)
    evalin('base','checkregempty=1;');
    close(handles.leadfigure)
    return
else
    evalin('base','checkregempty=0;');
end

%set(handles.previous,'visible','off'); set(handles.next,'visible','off');
setappdata(handles.leadfigure,'presentfiles',presentfiles)
setappdata(handles.leadfigure,'anchor',anchor)
setappdata(handles.leadfigure,'b0restanchor',b0restanchor)
setappdata(handles.leadfigure,'activevolume',1);
setappdata(handles.leadfigure,'options',options);

set(handles.checkatl,'Visible','off');


ea_mrcview(handles);

if isvalid(hObject)
    % Choose default command line output for ea_checkcoreg
    handles.output = hObject;

    % Update handles structure
    guidata(hObject, handles);
end

% UIWAIT makes ea_checkcoreg wait for user response (see UIRESUME)


function ea_mrcview(handles)

options=getappdata(handles.leadfigure,'options');

presentfiles=getappdata(handles.leadfigure,'presentfiles');
activevolume=getappdata(handles.leadfigure,'activevolume');
directory=getappdata(handles.leadfigure,'directory');
b0restanchor=getappdata(handles.leadfigure,'b0restanchor');

if activevolume==length(presentfiles)
    set(handles.disapprovebutn,'String','Disapprove & Close');
    set(handles.approvebutn,'String','Approve & Close');
else
    set(handles.approvebutn,'String','Approve & Next >>');
    set(handles.disapprovebutn,'String','Disapprove & Next >>');
end

currvol=presentfiles{activevolume};
if strcmp(currvol,'brainshift')
    ea_subcorticalrefine(options);
    close(handles.leadfigure);
    return
end
switch ea_stripext(currvol)
    case ea_stripext(options.prefs.gprenii)
    % case ea_stripext({options.prefs.gprenii, options.prefs.gtranii, options.prefs.tp_gctnii})
        handles.checkatl.Visible='on';
        [options] = ea_assignpretra(options);
        anchor=[ea_space,options.primarytemplate,'.nii'];
        set(handles.leadfigure,'Name',[options.patientname, ': Check Normalization']);

        ea_addnormmethods(handles,options,'coregmrmethod');

        if ~exist([directory,'ea_normmethod_applied.mat'],'file')
            method='';
        else
            method=load([directory,'ea_normmethod_applied.mat']);
            method=method.norm_method_applied{end};
        end

        set(handles.anchortxt,'String','Template (red wires):');
        set(handles.coregresultstxt,'String','Normalization results');
        set(handles.normsettings,'Visible','on');
        set(handles.recomputebutn,'String','(Re-) compute normalization using...');
        set(handles.coregmrmethod,'TooltipString','Choose a normalization method');
        set(handles.leadfigure,'Name',[options.patientname, ': Check Normalization']);
        set(gcf,'Name',[options.patientname, ': Check Normalization']);
    otherwise
        handles.checkatl.Visible='off';
        anchor=getappdata(handles.leadfigure,'anchor');
        set(handles.anchortxt,'String','Anchor modality (red wires):');
        set(handles.coregresultstxt,'String','Coregistration results');
        set(handles.leadfigure,'Name',[options.patientname, ': Check Coregistration']);

            switch currvol
                case ['tp_',options.prefs.ctnii_coregistered] % CT
                    ea_init_coregctpopup(handles, options, 'coregmrmethod');
                    if ~exist([directory,'ea_coregctmethod_applied.mat'],'file')
                        method='';
                    else
                        method=load([directory,'ea_coregctmethod_applied.mat']);
                        if iscell(method.coregct_method_applied)
                            method=method.coregct_method_applied{end};
                        else
                            method=method.coregct_method_applied;
                        end
                    end
                otherwise % MR
                    ea_init_coregmrpopup(handles,1);
                    if ~exist([directory,'ea_coregmrmethod_applied.mat'],'file')
                        method='';
                    else
                        method=load([directory,'ea_coregmrmethod_applied.mat']);
                        if isfield(method,ea_stripext(currvol)) % specific method used for this modality
                            method=method.(ea_stripext(currvol));
                        else
                            if isfield(method,'coregmr_method_applied')
                                if iscell(method.coregmr_method_applied)
                                    method=method.coregmr_method_applied{end};
                                else
                                    method=method.coregmr_method_applied;
                                end
                            else
                                method='';
                            end
                        end
                    end
            end

        set(handles.normsettings,'Visible','off');
        set(handles.recomputebutn,'String','(Re-) compute coregistration using...');
        set(handles.coregmrmethod,'TooltipString','Choose a coregistration method');
end

if ~exist([directory,'ea_coreg_approved.mat'],'file') % init
    for vol=1:length(presentfiles)
        approved.(ea_stripext(presentfiles{vol}))=0;
    end
    save([directory,'ea_coreg_approved.mat'],'-struct','approved');
else
    approved=load([directory,'ea_coreg_approved.mat']);
end
setappdata(handles.leadfigure,'method',method);

% show result:
if ~isempty(b0restanchor) && ~isempty(b0restanchor{activevolume}) % rest or b0 registration
    set(handles.substitute,'Visible','on');
    set(handles.substitute,'String',ea_getsubstitutes(options));
    checkfig=[directory,'checkreg',filesep,ea_stripext(currvol),'2',strrep(ea_stripext(b0restanchor{activevolume}),'mean','r'),'_',method,'.png'];
    set(handles.anchormod,'String',ea_stripext(b0restanchor{activevolume}));

else % normal anatomical 2 anatomical registration
    set(handles.substitute,'Visible','off');
    checkfig=[directory,'checkreg',filesep,ea_stripext(currvol),'2',ea_stripext(anchor),'_',method,'.png'];
    set(handles.anchormod,'String',ea_stripext(anchor));

end

set(handles.imgfn,'Visible','on');
set(handles.imgfn,'String',checkfig);
set(handles.imgfn,'TooltipString',checkfig);
switch ea_stripext(currvol)
    case ea_stripext(options.prefs.gprenii)
    % case ea_stripext({options.prefs.gprenii, options.prefs.gtranii, options.prefs.tp_gctnii})
        options=ea_assignpretra(options);
        anchorpath=[ea_space,options.primarytemplate];
    otherwise
        if ~isempty(b0restanchor) && ~isempty(b0restanchor{activevolume}) % rest or b0 registration
            anchorpath=[directory,ea_stripext(b0restanchor{activevolume})];
        else
            anchorpath=[directory,ea_stripext(anchor)];
        end
end

if ~exist(checkfig,'file')
    ea_gencheckregpair([directory,ea_stripext(currvol)],anchorpath,checkfig);

    if ~exist(checkfig,'file')
        checkfig=fullfile(ea_getearoot,'helpers','gui','coreg_msg.png');
        set(handles.imgfn,'String','');
        set(handles.imgfn,'Visible','off');
    end
end

setappdata(handles.leadfigure,'anchorpath',anchorpath);
im=imread(checkfig);
set(0,'CurrentFigure',handles.leadfigure);
set(handles.leadfigure,'CurrentAxes',handles.standardax);

imagesc(im);
axis off
axis equal

% textfields:
set(handles.depvolume,'String',[ea_stripext(currvol),'.nii']);


function [pretras]=ea_getsubstitutes(options)

[~,presentfiles]=ea_assignpretra(options);
for fi=1:length(presentfiles)
    if fi==1
        pretras{fi}=['Use ',presentfiles{fi}, ' (default)'];
    else
        pretras{fi}=['Substitute moving file with ',presentfiles{fi}];
    end
end


function presentfiles=ea_getall_coregcheck(options)
directory=[options.root,options.patientname,filesep];
[options,presentfiles]=ea_assignpretra(options);
% add postoperative volumes:
switch options.modality
    case 1 % MR
        if exist([directory,options.prefs.tranii_unnormalized],'file')
            presentfiles=[presentfiles;options.prefs.tranii_unnormalized];
        end
        if exist([directory,options.prefs.cornii_unnormalized],'file')
            presentfiles=[presentfiles;options.prefs.cornii_unnormalized];
        end
        if exist([directory,options.prefs.sagnii_unnormalized],'file')
            presentfiles=[presentfiles;options.prefs.sagnii_unnormalized];
        end
    case 2 % CT
        if exist([directory,'tp_',options.prefs.ctnii_coregistered],'file')
            presentfiles=[presentfiles;['tp_',options.prefs.ctnii_coregistered]];
        end
end

if exist([directory,options.prefs.fa2anat],'file')
    presentfiles=[presentfiles;options.prefs.fa2anat];
end

% now check if those are already approved (then don't show again):
todel=[];
for pf=1:length(presentfiles)
    if ea_reglocked(options,presentfiles{pf})
        todel=[todel,pf];
    end
end
presentfiles(todel)=[];


% --- Outputs from this function are returned to the command line.
function varargout = ea_checkcoreg_OutputFcn(hObject, eventdata, handles)
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
ea_busyaction('on',handles.leadfigure,'coreg');

options=getappdata(handles.leadfigure,'options');
options.overwriteapproved=1;
presentfiles=getappdata(handles.leadfigure,'presentfiles');
anchor=getappdata(handles.leadfigure,'anchor');
b0restanchor=getappdata(handles.leadfigure,'b0restanchor');

activevolume=getappdata(handles.leadfigure,'activevolume');
directory=getappdata(handles.leadfigure,'directory');
currvol=presentfiles{activevolume};

switch ea_stripext(currvol)
    case ea_stripext(options.prefs.gprenii)
        ea_delete([options.root,options.patientname,filesep,options.prefs.gprenii]);
        [options,anatspresent]=ea_assignpretra(options);
        for fi=2:length(anatspresent)
            ea_delete([options.root,options.patientname,filesep,'gl',anatspresent{fi}]);
        end

        options.normalize.method=getappdata(handles.leadfigure,'normmethod');
        options.normalize.method=options.normalize.method{get(handles.coregmrmethod,'Value')};
        options.normalize.methodn=get(handles.coregmrmethod,'Value');
        ea_dumpnormmethod(options,options.normalize.method,'normmethod'); % has to come first due to applynormalization.
        eval([options.normalize.method,'(options)']); % triggers the normalization function and passes the options struct to it.

        if options.modality == 2 % (Re-) compute tonemapped (normalized) CT
            ea_tonemapct_file(options,'mni');
        end

    case ea_stripext(['tp_',options.prefs.ctnii_coregistered]) % CT
        % Get CT coregistration method
        options.coregct.method = handles.coregmrmethod.String{handles.coregmrmethod.Value};

        % Run CT coregistration
        ea_coregpostopct(options);

        % Dump method
        ea_dumpmethod(options, 'coreg');

        ea_tonemapct_file(options,'native'); % (Re-) compute tonemapped (native space) CT
        ea_gencheckregfigs(options, 'coreg'); % generate checkreg figures

    case ea_stripext(options.prefs.fa2anat) % FA
        options.coregmr.method=get(handles.coregmrmethod,'String');
        options.coregmr.method=options.coregmr.method{get(handles.coregmrmethod,'Value')};
        ea_backuprestore([directory,options.prefs.fa]);
        ea_coregimages(options,[directory,options.prefs.fa],[directory,anchor],[directory,presentfiles{activevolume}],{},0);
        ea_dumpspecificmethod(handles,options.coregmr.method)

    otherwise % MR
        options.coregmr.method=get(handles.coregmrmethod,'String');
        options.coregmr.method=options.coregmr.method{get(handles.coregmrmethod,'Value')};
        if ~isempty(b0restanchor{activevolume})  % b0 or rest files

            thisrest=strrep(ea_stripext(b0restanchor{activevolume}),'mean','r');

            ea_delete([directory,thisrest,'2',ea_stripext(anchor),'_',ea_matext(options.coregmr.method)]);
            ea_delete([directory,ea_stripext(anchor),'2',thisrest,'_',ea_matext(options.coregmr.method)]);

            substitute=get(handles.substitute,'Value');
            [~,pf]=ea_assignpretra(options);
            useasanchor=pf{substitute};

            % in following line correct that useasanchor is the *moving*
            % image (since we're going from anchor to rest/b0.
            ea_coregimages(options,[directory,useasanchor],[directory,b0restanchor{activevolume}],[directory,presentfiles{activevolume}],{},1);
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
            ea_backuprestore([directory,presentfiles{activevolume}]);
            ea_coregimages(options,[directory,presentfiles{activevolume}],[directory,anchor],[directory,presentfiles{activevolume}],{},0);
        end
        ea_dumpspecificmethod(handles,options.coregmr.method)
end

% regenerate checkfig.
anchorpath=getappdata(handles.leadfigure,'anchorpath');

method=getappdata(handles.leadfigure,'method');

if ~isempty(b0restanchor{activevolume})
   anchor=b0restanchor{activevolume};
end

checkfig=[directory,'checkreg',filesep,ea_stripext(currvol),'2',strrep(ea_stripext(anchor),'mean','r'),'_',method,'.png'];

ea_gencheckregpair([directory,ea_stripext(currvol)],anchorpath,checkfig);
% now disapprove again since this new computation hasn't been approved yet.
approved=load([directory,'ea_coreg_approved.mat']);
approved.(ea_stripext(currvol))=0;
save([directory,'ea_coreg_approved.mat'],'-struct','approved');

ea_mrcview(handles)
title = get(handles.leadfigure, 'Name');    % Fix title
ea_chirp(options);
ea_busyaction('off',handles.leadfigure,'coreg');
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


function ea_dumpspecificmethod(handles,method)
options=getappdata(handles.leadfigure,'options');

presentfiles=getappdata(handles.leadfigure,'presentfiles');
anchor=getappdata(handles.leadfigure,'anchor');
activevolume=getappdata(handles.leadfigure,'activevolume');
directory=getappdata(handles.leadfigure,'directory');
try
    m=load([directory,'ea_coregmrmethod_applied.mat']);
catch
    m=struct;
end
m.(ea_stripext(presentfiles{activevolume}))=method;
save([directory,'ea_coregmrmethod_applied.mat'],'-struct','m');


% --- Executes on button press in approvebutn.
function approvebutn_Callback(hObject, eventdata, handles)
% hObject    handle to approvebutn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ea_busyaction('on',handles.leadfigure,'coreg');

options=getappdata(handles.leadfigure,'options');
presentfiles=getappdata(handles.leadfigure,'presentfiles');
activevolume=getappdata(handles.leadfigure,'activevolume');
directory=getappdata(handles.leadfigure,'directory');
currvol=presentfiles{activevolume};

switch ea_stripext(currvol)
    case ea_stripext(options.prefs.gprenii)
    case ea_stripext(['tp_',options.prefs.ctnii_coregistered])
    otherwise % make sure method gets logged for specific volume.
        method=getappdata(handles.leadfigure,'method');
        if exist([directory,'ea_coregmrmethod_applied.mat'],'file')
            m=load([directory,'ea_coregmrmethod_applied.mat']);
        end
        m.(ea_stripext(currvol))=method;
        save([directory,'ea_coregmrmethod_applied.mat'],'-struct','m');
end

approved=load([directory,'ea_coreg_approved.mat']);
try
    wasapprovedalready=approved.(ea_stripext(currvol));
catch
    wasapprovedalready=0;
end
approved.(ea_stripext(currvol))=1;
if strcmp(ea_stripext(currvol),ea_stripext(options.prefs.gprenii))
    [options,preniis]=ea_assignpretra(options); % get all preop versions
    allcoreg=1; % check if all preniis are already approved
    for pn=2:length(preniis)
        if ~approved.(ea_stripext(preniis{pn}))
            allcoreg=0;
        end
    end
    if allcoreg
        approved.(ea_stripext(currvol))=2; % set to permanent approved =2 normalization. This will not be overriden no matter what (as long is override flag is not set).
    else
        ea_warning('You approved normalization before all preoperative co-registrations were approved. Lead-DBS will still override / redo normalization if applying a multispectral method.');
    end
else
    if isfield(approved,ea_stripext(options.prefs.gprenii))
        [~,preopfiles]=ea_assignpretra(options);
        if ismember([ea_stripext(currvol),'.nii'],preopfiles)
            if approved.(ea_stripext(options.prefs.gprenii))==2
                % now in this situation we had the normalization approved before
                % all coregistrations were approved. This could lead to suboptimal
                % normalizations *only* if a multispectral protocol is used. Thus
                % we set the normalization approval rate to 1. This way, it will
                % still be overriden in case of running a multispectral
                % normalization.
                if ~wasapprovedalready
                    ea_warning('Normalization had been approved before all preoperative co-registrations were approved. Lead-DBS will still override / redo normalization if applying a multispectral method.');
                    approved.stripex(options.prefs.gprenii)=1; % this will be overriden when using a multispectral normalization.
                end
            end
        end
    end
end

save([directory,'ea_coreg_approved.mat'],'-struct','approved');
if strcmp(computer('arch'),'maci64')
    system(['xattr -wx com.apple.FinderInfo "0000000000000000000400000000000000000000000000000000000000000000" ',ea_path_helper([directory,ea_stripext(currvol),'.nii'])]);
end

presentfiles=getappdata(handles.leadfigure,'presentfiles');
anchor=getappdata(handles.leadfigure,'anchor');
activevolume=getappdata(handles.leadfigure,'activevolume');

if activevolume==length(presentfiles)
    close(handles.leadfigure); % make an exit
    return
elseif (activevolume==length(presentfiles)-1 && strcmp(presentfiles{end},'brainshift'))
    close(handles.leadfigure); % make an exit
    ea_subcorticalrefine(options);
    return
else
    activevolume=activevolume+1;
end
setappdata(handles.leadfigure,'activevolume',activevolume);

ea_mrcview(handles);
title = get(handles.leadfigure, 'Name');    % Fix title
ea_busyaction('off',handles.leadfigure,'coreg');
set(handles.leadfigure, 'Name', title);


% --- Executes on selection change in coregmrmethod.
function coregmrmethod_Callback(hObject, eventdata, handles)
% hObject    handle to coregmrmethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns coregmrmethod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from coregmrmethod

options=getappdata(handles.leadfigure,'options');

presentfiles=getappdata(handles.leadfigure,'presentfiles');
activevolume=getappdata(handles.leadfigure,'activevolume');
currvol=presentfiles{activevolume};
% init retry popup:
if strcmp(currvol,'glanat.nii')
    ea_switchnormmethod(handles,'coregmrmethod');
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
options=getappdata(handles.leadfigure,'options');
presentfiles=getappdata(handles.leadfigure,'presentfiles');
activevolume=getappdata(handles.leadfigure,'activevolume');
b0restanchor=getappdata(handles.leadfigure,'b0restanchor');

currvol=presentfiles{activevolume};
switch ea_stripext(currvol)
    case ea_stripext(options.prefs.gprenii)
        ea_show_normalization(options);
    otherwise
        presentfiles=getappdata(handles.leadfigure,'presentfiles');
        anchor=getappdata(handles.leadfigure,'anchor');
        activevolume=getappdata(handles.leadfigure,'activevolume');
        directory=getappdata(handles.leadfigure,'directory');

        options.moving=[directory,presentfiles{activevolume}];
        if ~isempty(b0restanchor{activevolume})
            options.fixed=[directory,b0restanchor{activevolume}];
            options.tag=[presentfiles{activevolume},' & ',b0restanchor{activevolume}];
        else
            options.fixed=[directory,anchor];
            options.tag=[presentfiles{activevolume},' & ',anchor];
        end

        ea_show_coregistration(options);
end


% --- Executes on button press in normsettings.
function normsettings_Callback(hObject, eventdata, handles)
% hObject    handle to normsettings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
currentNormMethod=getappdata(handles.normsettings,'currentNormMethod');
ea_shownormsettings(currentNormMethod,handles)


% --- Executes on button press in disapprovebutn.
function disapprovebutn_Callback(hObject, eventdata, handles)
% hObject    handle to disapprovebutn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% hObject    handle to approvebutn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ea_busyaction('on',handles.leadfigure,'coreg');

options=getappdata(handles.leadfigure,'options');
presentfiles=getappdata(handles.leadfigure,'presentfiles');
activevolume=getappdata(handles.leadfigure,'activevolume');
directory=getappdata(handles.leadfigure,'directory');
b0restanchor=getappdata(handles.leadfigure,'b0restanchor');
currvol=presentfiles{activevolume};

approved=load([directory,'ea_coreg_approved.mat']);

approved.(ea_stripext(currvol))=0;
save([directory,'ea_coreg_approved.mat'],'-struct','approved');
if strcmp(computer('arch'),'maci64')
    system(['xattr -wx com.apple.FinderInfo "0000000000000000000C00000000000000000000000000000000000000000000" ',ea_path_helper([directory,ea_stripext(currvol),'.nii'])]);
end

switch ea_stripext(currvol)
    case ea_stripext(options.prefs.gprenii)

    case ea_stripext(['tp_',options.prefs.ctnii_coregistered])

    otherwise % make sure method gets unlogged for specific volume.
        method=getappdata(handles.leadfigure,'method');
        if exist([directory,'ea_coregmrmethod_applied.mat'],'file')
            m=load([directory,'ea_coregmrmethod_applied.mat']);
        else
            m=struct;
        end
        if isfield(m,ea_stripext(currvol))
            m=rmfield(m,ea_stripext(currvol));
        end
        save([directory,'ea_coregmrmethod_applied.mat'],'-struct','m');
end

if ~isempty(b0restanchor{activevolume})
    thisrest=strrep(ea_stripext(b0restanchor{activevolume}),'mean','r');
    ea_cleandownstream(directory,thisrest)
end

presentfiles=getappdata(handles.leadfigure,'presentfiles');
anchor=getappdata(handles.leadfigure,'anchor');
activevolume=getappdata(handles.leadfigure,'activevolume');

if activevolume==length(presentfiles)
    close(handles.leadfigure); % make an exit
    return
else
    activevolume=activevolume+1;
end

setappdata(handles.leadfigure,'activevolume',activevolume);
ea_mrcview(handles);
try
    title = get(handles.leadfigure, 'Name');    % Fix title
    ea_busyaction('off',handles.leadfigure,'coreg');
    set(handles.leadfigure, 'Name', title);
end


% --- Executes on button press in back.
function back_Callback(hObject, eventdata, handles)
% hObject    handle to back (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ea_busyaction('on',handles.leadfigure,'coreg');

activevolume=getappdata(handles.leadfigure,'activevolume');

if activevolume==1
    ea_busyaction('off',handles.leadfigure,'coreg');
    return
else
    activevolume=activevolume-1;
    setappdata(handles.leadfigure,'activevolume',activevolume);
end

ea_mrcview(handles);
title = get(handles.leadfigure, 'Name');    % Fix title
ea_busyaction('off',handles.leadfigure,'coreg');
set(handles.leadfigure, 'Name', title);


% --- Executes on button press in refreshview.
function refreshview_Callback(hObject, eventdata, handles)
% hObject    handle to refreshview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ea_busyaction('on',handles.leadfigure,'coreg');

options=getappdata(handles.leadfigure,'options');
presentfiles=getappdata(handles.leadfigure,'presentfiles');
activevolume=getappdata(handles.leadfigure,'activevolume');
directory=getappdata(handles.leadfigure,'directory');
currvol=presentfiles{activevolume};
anchorpath=getappdata(handles.leadfigure,'anchorpath');
method=getappdata(handles.leadfigure,'method');

anchor=getappdata(handles.leadfigure,'anchor');
switch ea_stripext(currvol)
    case ea_stripext(options.prefs.gprenii)
        options=ea_assignpretra(options);
        anchor=[ea_space,options.primarytemplate,'.nii'];
    case ea_stripext(['tp_',options.prefs.ctnii_coregistered]) % make sure tp matches rpostop_ct.
        ea_tonemapct_file(options,'native');
end
b0restanchor=getappdata(handles.leadfigure,'b0restanchor');
if ~isempty(b0restanchor{activevolume})
   anchor=b0restanchor{activevolume};
end

checkfig=[directory,'checkreg',filesep,ea_stripext(currvol),'2',strrep(ea_stripext(anchor),'mean','r'),'_',method,'.png'];
ea_delete(checkfig);
ea_gencheckregpair([directory,ea_stripext(currvol)],anchorpath,checkfig);
ea_mrcview(handles); % refresh
title = get(handles.leadfigure, 'Name');    % Fix title
ea_busyaction('off',handles.leadfigure,'coreg');
set(handles.leadfigure, 'Name', title);


% --- Executes on button press in checkatl.
function checkatl_Callback(hObject, eventdata, handles)
% hObject    handle to checkatl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
options=getappdata(handles.leadfigure,'options');
presentfiles=getappdata(handles.leadfigure,'presentfiles');
directory=getappdata(handles.leadfigure,'directory');

ea_checkstructures(options);


% --- Executes on button press in openpatientdir.
function openpatientdir_Callback(hObject, eventdata, handles)
% hObject    handle to openpatientdir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ea_opendir(getappdata(handles.leadfigure,'directory'));


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
