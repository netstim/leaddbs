function ea_preprocess_fmri(options)

directory=[options.root,options.patientname,filesep];

V=spm_vol([directory,options.prefs.rest]);
signallength=length(V);

%% run sequence of proxyfunctions (below):
ea_realign_fmri(signallength,options); % realign fMRI

ea_newseg(fullfile(directory,options.prefs.prenii_unnormalized),0,1); % Segment anat

ea_coreg_pre2fmri(options); % register pre 2 fmri (for timecourse-extraction).

ea_smooth_fmri(signallength,options); % slightly smooth fMRI data


function ea_realign_fmri(signallength,options)
%% realign fmri.
directory=[options.root,options.patientname,filesep];
if ~exist([directory,'r',options.prefs.rest],'file')

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

% Re-calculate mean re-aligned image if not found
if ~exist([directory, 'mean', options.prefs.rest], 'file')
    ea_meanimage([directory, 'r', options.prefs.rest], ['mean', options.prefs.rest]);
end

if isfield(options, 'overwriteapproved') && options.overwriteapproved
    overwrite = 1;
else
    overwrite = 0;
end

anatfname=ea_stripext(options.prefs.prenii_unnormalized);
refname=['r',ea_stripext([options.prefs.rest])];
reference=['mean',options.prefs.rest]; % okay here to not use the hd version of the image since this is about the csf/wm masks.

% For this pair of approved coregistations, find out which method to use -
% irrespective of the current selection in coregmethod.
coregmethodsused=load([directory,'ea_coregmrmethod_applied.mat']);
fn=fieldnames(coregmethodsused);
for field=1:length(fn)
    if contains(fn{field},ea_stripext(options.prefs.rest))
        if ~isempty(coregmethodsused.(fn{field}))
            disp(['For this pair of coregistrations, the user specifically approved the ',coregmethodsused.(fn{field}),' method, so we will overwrite the current global options and use this transform.']);
            options.coregmr.method=coregmethodsused.(fn{field});
        end
        break
    end
end

% Disable Hybrid coregistration
coregmethod = strrep(options.coregmr.method, 'Hybrid SPM & ', '');
options.coregmr.method = coregmethod;

% Check if the corresponding transform already exists
xfm = [anatfname, '2', refname, '_', lower(coregmethod), '\d*\.(mat|h5)$'];
transform = ea_regexpdir(directory, xfm, 0);

if numel(transform) == 0 || overwrite
    if numel(transform) == 0
        warning('Transformation not found! Running coregistration now!');
    end

    transform = ea_coregimages(options,[directory,options.prefs.prenii_unnormalized],...
        [directory,reference],...
        [directory,'r',ea_stripext(options.prefs.rest),'_',options.prefs.prenii_unnormalized],...
        [],1,[],1);
    % Fix transformation names, replace 'mean' by 'r' for fMRI
    if strcmp(reference, ['mean', options.prefs.rest])
        cellfun(@(f) movefile(f, strrep(f, 'mean', 'r')), transform);
        transform = strrep(transform, 'mean', 'r');
    end
    transform = transform{1}; % Forward transformation
else
    if numel(transform) > 1
        warning(['Multiple transformations of the same type found! ' ...
            'Will use the last one:\n%s'], transform{end});
    end
    transform = transform{end};
end

ea_apply_coregistration([directory,'mean',options.prefs.rest], ...
    [directory,options.prefs.prenii_unnormalized], ...
    [directory,'r',ea_stripext(options.prefs.rest),'_',options.prefs.prenii_unnormalized], ...
    transform, 'linear');

% segmented anat images registered to mean rest image
for i=1:3
    ea_apply_coregistration([directory,'mean',options.prefs.rest], ...
        [directory,'c',num2str(i),options.prefs.prenii_unnormalized], ...
        [directory,'r',ea_stripext(options.prefs.rest),'_c',num2str(i),options.prefs.prenii_unnormalized], ...
        transform, 'linear');
end


function ea_smooth_fmri(signallength,options)
directory=[options.root,options.patientname,filesep];

filetimepts = ea_appendVolNum([directory,'r',options.prefs.rest], 1:signallength);
if ~exist([directory,'sr',options.prefs.rest],'file')
    matlabbatch{1}.spm.spatial.smooth.data = filetimepts;
    matlabbatch{1}.spm.spatial.smooth.fwhm = [6 6 6];
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    matlabbatch{1}.spm.spatial.smooth.im = 0;
    matlabbatch{1}.spm.spatial.smooth.prefix = 's';
    jobs{1}=matlabbatch;
    spm_jobman('run',jobs);
    clear jobs matlabbatch
end
