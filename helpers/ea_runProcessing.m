function ea_runProcessing(dataset, patientID, opts)
arguments
    dataset         {mustBeFolder}  % dataset folder
    patientID       {mustBeText}    % patient ID
    opts.process    {mustBeMember(opts.process, {'pre', 'post', 'both'})} = 'both'
end

%% get bids structure
options = ea_getptopts(fullfile(dataset, 'derivatives', 'leaddbs', ['sub-', patientID]));
options.overwriteapproved = 0;

%% Set electrode model
uiprefs = load(options.subj.uiprefs);
options.elmodel = uiprefs.elmodel;
options = ea_resolve_elspec(options);

%% Preprocessing
% Pre-op images
if ismember(opts.process, {'pre', 'both'})
    fields = fieldnames(options.subj.preopAnat);
    for i=1:length(fields)
        if ~isfile(options.subj.preopAnat.(fields{i}).preproc)
            % Copy files
            ea_mkdir(fileparts(options.subj.preopAnat.(fields{i}).preproc));
    
            % Copy file to preproc, take care of .nii.gz raw image
            if strcmp(options.prefs.niiFileExt, '.nii')
                copyfile(options.subj.preopAnat.(fields{i}).raw, [options.subj.preopAnat.(fields{i}).preproc, '.gz']);
                gunzip([options.subj.preopAnat.(fields{i}).preproc, '.gz']);
                delete([options.subj.preopAnat.(fields{i}).preproc, '.gz']);
            else
                copyfile(options.subj.preopAnat.(fields{i}).raw, options.subj.preopAnat.(fields{i}).preproc);
            end
    
            % Run reorientation, cropping and bias field correction
            ea_anatpreprocess(options.subj.preopAnat.(fields{i}).preproc);
    
            % Preprocessing steps only for pre-op anchor image
            % Skip reslicing anchor image to 0.7^3 resolution when only pre-op images exist
            if isfield(options.subj, 'postopAnat')
                ea_resliceanat(options.subj.preopAnat.(fields{i}).preproc);
            end
        end
    end
end

% Post-op images
if ismember(opts.process, {'post', 'both'})
    fields = fieldnames(options.subj.postopAnat);
    for i=1:length(fields)
        if ~isfile(options.subj.postopAnat.(fields{i}).preproc)
            ea_mkdir(fileparts(options.subj.postopAnat.(fields{i}).preproc));
            % Copy file to preproc, take care of .nii.gz raw image
            if strcmp(options.prefs.niiFileExt, '.nii')
                copyfile(options.subj.postopAnat.(fields{i}).raw, [options.subj.postopAnat.(fields{i}).preproc, '.gz']);
                gunzip([options.subj.postopAnat.(fields{i}).preproc, '.gz']);
                delete([options.subj.postopAnat.(fields{i}).preproc, '.gz']);
            else
                copyfile(options.subj.postopAnat.(fields{i}).raw, options.subj.postopAnat.(fields{i}).preproc);
            end
        end
    end
end

ea_cprintf('*Comments', '%s: preprocessing finished.\n', patientID);

%% Co-registration
% Set primary template
subjAnchor = regexprep(options.subj.AnchorModality, '[^\W_]+_', '');
if ismember(subjAnchor, fieldnames(options.bids.spacedef.norm_mapping))
    options.primarytemplate = options.bids.spacedef.norm_mapping.(subjAnchor);
else
    options.primarytemplate = options.bids.spacedef.misfit_template;
end

options.coregmr.method = 'ANTs (Avants 2008)'; % 'SPM (Friston 2007)'
options.coregct.method = 'ANTs (Avants 2008)';

% Pre-coregister pre-op anchor image
if ismember(opts.process, {'pre', 'both'})
    if isfield(options.subj, 'preopAnat')
        if ~isfile(options.subj.preopAnat.(options.subj.AnchorModality).coreg)
            ea_precoreg(options.subj.preopAnat.(options.subj.AnchorModality).preproc, ... % Input anchor image
                options.primarytemplate, ... % Template to use
                options.subj.preopAnat.(options.subj.AnchorModality).coreg, ... % Output pre-coregistered image
                options.subj.coreg.transform.(options.subj.AnchorModality)); % % Pre-coregistration transform
    
            % Check if anchor image has been properly pre-coregistered. If
            % not, fallback to preproc image.
            try
                load_nii(options.subj.preopAnat.(options.subj.AnchorModality).coreg);
            catch
                ea_cprintf('CmdWinWarnings', 'Anchor image was not properly pre-coregistered. Fallback to preproc image instead.\n');
                ea_delete(options.subj.coreg.transform.(options.subj.AnchorModality));
                ea_mkdir(fileparts(options.subj.preopAnat.(options.subj.AnchorModality).coreg));
                copyfile(options.subj.preopAnat.(options.subj.AnchorModality).preproc, options.subj.preopAnat.(options.subj.AnchorModality).coreg);
            end
        end
    end
    
    ea_coregpreopmr(options);
end

if ismember(opts.process, {'post', 'both'})
    if strcmp(options.subj.postopModality, 'MRI')
        % Coregister post-op MRI to pre-op MRI
        ea_coregpostopmr(options);
    elseif strcmp(options.subj.postopModality, 'CT')
        % Coregister post-op CT to pre-op MRI
        ea_coregpostopct(options);
    end
end

% Generate checkreg figures for coregistration
ea_gencheckregfigs(options, 'coreg');

ea_cprintf('*Comments', '%s: coregistration finished.\n', patientID)

%% Normalization
options.normalize.method = 'ANTs (Avants 2008)';

if ismember(opts.process, {'pre', 'both'})
    ea_normalize(options);
end

if ismember(opts.process, {'post'})
    ea_apply_normalization(options);
end

ea_gencheckregfigs(options, 'norm');

ea_cprintf('*Comments', '%s: normalization finished.\n', patientID)

%% Brainshift correction
if ismember(opts.process, {'post', 'both'})
    options.autobrainshift = 1;
    options.d2.write = 0;
    options.d3.write = 0;
    options.scrf.mask = 'Coarse';
    ea_subcorticalrefine(options);
    
    ea_cprintf('*Comments', '%s: brain shift correction finished.\n', patientID);
end

%% Pre-reconstruction
if ismember(opts.process, {'post', 'both'})
    options.hybridsave=1;
    options.reconmethod = 'PaCER (Husch 2017)';
    ea_mkdir(options.subj.reconDir);
    try
        [coords_mm,trajectory,markers]=ea_runpacer(options);
        options.native = 1;
        ea_save_reconstruction(coords_mm, trajectory, markers, options.elmodel, 0, options);
    
        ea_cprintf('*Comments', '%s: pre-reconstruction finished.\n', patientID);
    catch ME
        ea_cprintf('CmdWinErrors', '%s: pre-reconstruction failed.\n%s\n', patientID, ME.message);
    end
end
