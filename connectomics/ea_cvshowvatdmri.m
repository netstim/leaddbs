function ea_cvshowvatdmri(resultfig,directory,handles,selectedparc,options)

% get everything we need to know from handles:
if isstruct(handles) % called from GUI
    vatmodality=get(handles.vatmodality,'String');
    vatmodality=vatmodality{get(handles.vatmodality,'Value')};
    vs=get(handles.vatseed,'String'); % dont need this below
    vss=get(handles.vatseed,'Value'); % dont need this below
    vsname=vs{vss};
    [usevat,dimensionality,~,sides]=ea_checkvatselection(handles);
    thresh=get(handles.vatthresh,'String');
elseif iscell(handles) % called from lead_group
    vatmodality=handles{1};
    vsname=[handles{2},'_',options.groupid];
    thresh='auto';
    usevat={'right','left'};
    dimensionality=2; % how many ROI.
    sides=[1,2];
end

% fibers filename
if isstruct(vatmodality)
    fibersfile=vatmodality;
else
    switch vatmodality
        case 'Patient''s fiber tracts'
            fibersfile=[directory,'connectomes',filesep,'dMRI',filesep,options.prefs.FTR_normalized];
        otherwise
            fibersfile=[ea_getconnectomebase('dmri'),vatmodality,filesep,'data.mat'];
    end
end

subPrefix = ['sub-', options.subj.subjId];
hemiTag = regexprep(usevat, {'right', 'left'}, {'R', 'L'});

try
    load([directory,'stimulations',filesep,ea_nt(options),vsname,filesep,subPrefix,'_desc-stimparameters.mat'], 'S');
catch
    ea_error(['Could not find stimulation parameters for ',directory,ea_nt(options),vsname,'.']);
end

modelLabel = ea_simModel2Label(S.model);
vatPrefix = [subPrefix, '_sim-binary_model-', modelLabel, '_hemi-'];

% seed filename
seedfile=cell(length(hemiTag),1);
for v=1:length(hemiTag)
    seedfile{v}=[directory,'stimulations',filesep,ea_nt(options),vsname,filesep,vatPrefix,hemiTag{options.sides(v)},'.nii'];
end

targetsfile=[ea_space(options,'labeling'),selectedparc,'.nii'];

options.writeoutstats=1;
options.writeoutpm = 0;

try
    load([directory,'connvisfibers/','fiberstate.mat'])
end

[changedstates,ret]=ea_checkfschanges(resultfig,fibersfile,seedfile,targetsfile,thresh,'vat');

% Option to save fibers and check for workspace to load, pass via options
if isstruct(handles) && get(handles.savefibers,'Value')
    options.savefibers.save = 1;
    savedir = ea_getfsstatesavedir(directory,vsname,fibersfile,seedfile,targetsfile,thresh,'vat');

    if ~exist(savedir,'dir')
    	mkdir(savedir);
    end

    if ~exist([savedir,'workspace.mat'],'file')
        mode = 'vat'; save([savedir,'fiberstate.mat'],'fibersfile','seedfile','targetsfile','thresh','mode')
        options.savefibers.load=0;
    elseif exist([savedir,'workspace.mat'],'file') && exist([savedir,'fiberstate.mat'],'file')
        options.savefibers.load=1;
    end

    options.savefibers.dir=savedir;
end


if ~ret % something has changed since last time.
    ea_deletePL(resultfig,'PL','vat');
    if dimensionality % one of the vat checkboxes is active
        if isstruct(handles) % call from GUI
            [~,thresh]=ea_cvshowfiberconnectivities(resultfig,fibersfile,seedfile,targetsfile,thresh,sides,options,S,changedstates,'vat',get(handles.vizvat_regs,'Value'),get(handles.vizvat_labs,'Value')); % 'vat' only used for storage of changes.
        elseif iscell(handles) % call from Lead Group
            [~,thresh]=ea_cvshowfiberconnectivities(resultfig,fibersfile,seedfile,targetsfile,thresh,sides,options,S,changedstates,'vat',0,0); % 'vat' only used for storage of changes.
        end
        if isstruct(handles)
            set(handles.vatthreshis,'String',num2str(thresh));
        end
    end
end
