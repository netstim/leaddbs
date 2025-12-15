function ea_preprocess_fmri(options)

directory=[options.root,options.patientname,filesep];

V=spm_vol([directory,options.prefs.rest]);
signallength=length(V);

%% run sequence of proxyfunctions (below):
ea_realign_fmri(signallength,options); % realign fMRI

% BIDS FIX: Check if segmentation already exists before running ea_newseg
[anatDir, anatName] = fileparts(fullfile(directory, options.prefs.prenii_unnormalized));
c1File = fullfile(anatDir, ['c1', anatName, '.nii']);
c2File = fullfile(anatDir, ['c2', anatName, '.nii']);
c3File = fullfile(anatDir, ['c3', anatName, '.nii']);

if ~exist(c1File, 'file') || ~exist(c2File, 'file') || ~exist(c3File, 'file')
    % Segmentation not done yet - but use standard SPM12 segmentation instead of ea_newseg
    % to avoid dependency on Lorio-Draganski templates
    disp('Running SPM12 segmentation for anat...');
    
    anatFile = fullfile(directory, options.prefs.prenii_unnormalized);
    matlabbatch{1}.spm.spatial.preproc.channel.vols = {[anatFile, ',1']};
    matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
    matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
    matlabbatch{1}.spm.spatial.preproc.channel.write = [0 0];
    
    % Tissue classes: GM, WM, CSF
    for ti = 1:3
        matlabbatch{1}.spm.spatial.preproc.tissue(ti).tpm = {[spm('Dir'), '/tpm/TPM.nii,', num2str(ti)]};
        matlabbatch{1}.spm.spatial.preproc.tissue(ti).ngaus = ti;
        matlabbatch{1}.spm.spatial.preproc.tissue(ti).native = [1 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(ti).warped = [0 0];
    end
    
    % Other tissues
    for ti = 4:6
        matlabbatch{1}.spm.spatial.preproc.tissue(ti).tpm = {[spm('Dir'), '/tpm/TPM.nii,', num2str(ti)]};
        matlabbatch{1}.spm.spatial.preproc.tissue(ti).ngaus = 2;
        matlabbatch{1}.spm.spatial.preproc.tissue(ti).native = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(ti).warped = [0 0];
    end
    
    matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
    matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
    matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
    matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
    matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
    matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
    matlabbatch{1}.spm.spatial.preproc.warp.write = [0 0];
    
    spm_jobman('run', {matlabbatch});
    clear matlabbatch
    disp('Done.');
else
    disp('Segmentation already exists, skipping...');
end

ea_coreg_pre2fmri(options); % register pre 2 fmri (for timecourse-extraction).

ea_smooth_fmri(signallength,options); % slightly smooth fMRI data


function ea_realign_fmri(signallength,options)
%% realign fmri.
directory=[options.root,options.patientname,filesep];

% BIDS FIX: Add prefix to filename only
rest_r = ea_prependFilename(options.prefs.rest, 'r');

if ~exist([directory,rest_r],'file')

    restbackup = ea_niifileparts([directory,options.prefs.rest]);
    copyfile([directory,options.prefs.rest], restbackup);

    disp('Realignment of rs-fMRI data...');
    filetimepts = ea_appendVolNum([directory,options.prefs.rest], 1:signallength);
    matlabbatch{1}.spm.spatial.realign.estwrite.data = {filetimepts};
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 1;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 0;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = {''};
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
    spm_jobman('run',{matlabbatch});
    clear matlabbatch;
    disp('Done.');
%    ea_reslice_nii([directory, options.prefs.rest], [directory, 'r',options.prefs.rest], [1,1,1], 0, 0, 1, [], [],0)
    movefile(restbackup, [directory, options.prefs.rest]);
end


function ea_coreg_pre2fmri(options)
directory=[options.root,options.patientname,filesep];
% Disable Hybrid coregistration
coregmethod = options.coregmr.method;
options.coregmr.method = strrep(coregmethod, 'Hybrid SPM & ', '');

% BIDS FIX: Helper to add prefix to filename (not path)
% For BIDS: 'preprocessing/func/sub-XXX_bold.nii' -> 'preprocessing/func/rsub-XXX_bold.nii'
% For classic: 'rest.nii' -> 'rrest.nii'
rest_r = ea_prependFilename(options.prefs.rest, 'r');
rest_mean = ea_prependFilename(options.prefs.rest, 'mean');

% Re-calculate mean re-aligned image if not found
if ~exist([directory, rest_mean], 'file')
    ea_meanimage([directory, rest_r], rest_mean);
end

if isfield(options, 'overwriteapproved') && options.overwriteapproved
    overwrite = 1;
else
    overwrite = 0;
end

anatfname=ea_stripext(options.prefs.prenii_unnormalized);
rest_r_noext = ea_stripext(rest_r);
refname = ea_stripext(rest_r);  % Get basename without path for matching
if contains(refname, filesep)
    [~, refname] = fileparts(refname);
end
reference=rest_mean; % okay here to not use the hd version of the image since this is about the csf/wm masks.

% BIDS FIX: Check for coregistration method log (in both locations)
coregLogFile = fullfile(directory, 'coregistration', 'log', 'ea_coregmrmethod_applied.mat');
if ~exist(coregLogFile, 'file')
    coregLogFile = fullfile(directory, 'ea_coregmrmethod_applied.mat');
end

if exist(coregLogFile, 'file')
    % For this pair of approved coregistations, find out which method to use
    coregmethodsused = load(coregLogFile);
    fn = fieldnames(coregmethodsused);
    for field = 1:length(fn)
        if contains(fn{field}, ea_stripext(options.prefs.rest))
            if ~isempty(coregmethodsused.(fn{field}))
                disp(['For this pair of coregistrations, the user specifically approved the ',coregmethodsused.(fn{field}),' method, so we will overwrite the current global options and use this transform.']);
                options.coregmr.method = coregmethodsused.(fn{field});
            end
            break
        end
    end
end

% Disable Hybrid coregistration
coregmethod = strrep(options.coregmr.method, 'Hybrid SPM & ', '');
options.coregmr.method = coregmethod;

% BIDS FIX: Use dir() instead of ea_regexpdir to find transformation
% Check if the corresponding transform already exists
[~, anatBaseName] = ea_niifileparts(options.prefs.prenii_unnormalized);
pattern = [anatBaseName, '2', refname, '_', lower(coregmethod), '*.mat'];

% Search in preprocessing/anat first, then root
searchDirs = {fullfile(directory, 'preprocessing', 'anat'), directory};
transform = {};

for iDir = 1:length(searchDirs)
    if exist(searchDirs{iDir}, 'dir')
        files = dir(fullfile(searchDirs{iDir}, pattern));
        if ~isempty(files)
            transform = {fullfile(searchDirs{iDir}, files(end).name)};
            break;
        end
    end
end

% BIDS FIX: Build output name correctly (needed for both coregister and apply)
% Get basename of anatomical (without path)
[~, anatBaseName] = fileparts(ea_stripext(options.prefs.prenii_unnormalized));
rest_r_anat = ea_prependFilename(rest_r_noext, '', ['_', anatBaseName, '.nii']);

if numel(transform) == 0 || overwrite
    if numel(transform) == 0
        warning('Transformation not found! Running coregistration now!');
    end
    
    transform = ea_coregimages(options,[directory,options.prefs.prenii_unnormalized],...
        [directory,reference],...
        [directory,rest_r_anat],...
        [],1,[],1);
    
    transform = transform{1}; % Forward transformation
else
    if numel(transform) > 1
        warning(['Multiple transformations of the same type found! ' ...
            'Will use the last one:\n%s'], transform{end});
    end
    transform = transform{end};
end

ea_apply_coregistration([directory,reference], ...
    [directory,options.prefs.prenii_unnormalized], ...
    [directory,rest_r_anat], ...
    transform, 'linear');

% segmented anat images registered to mean rest image
for i=1:3
    % BIDS FIX: Add c prefix to filename only
    anat_c = ea_prependFilename(options.prefs.prenii_unnormalized, ['c', num2str(i)]);
    rest_r_c_anat = ea_prependFilename(rest_r_noext, '', ['_c', num2str(i), anatBaseName, '.nii']);
    ea_apply_coregistration([directory,reference], ...
        [directory,anat_c], ...
        [directory,rest_r_c_anat], ...
        transform, 'linear');
end


function ea_smooth_fmri(signallength,options)
directory=[options.root,options.patientname,filesep];

% BIDS FIX: Add prefix to filename only
rest_r = ea_prependFilename(options.prefs.rest, 'r');
rest_sr = ea_prependFilename(options.prefs.rest, 'sr');

filetimepts = ea_appendVolNum([directory,rest_r], 1:signallength);
if ~exist([directory,rest_sr],'file')
    matlabbatch{1}.spm.spatial.smooth.data = filetimepts;
    matlabbatch{1}.spm.spatial.smooth.fwhm = [6 6 6];
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    matlabbatch{1}.spm.spatial.smooth.im = 0;
    matlabbatch{1}.spm.spatial.smooth.prefix = 's';
    jobs{1}=matlabbatch;
    spm_jobman('run',jobs);
    clear jobs matlabbatch
end
