function ea_warp_parcellation(reference,options)
directory=[options.root,options.patientname,filesep];

% Regenerate the template or not
if isfield(options, 'overwriteapproved') && options.overwriteapproved
    overwrite = 1;
else
    overwrite = 0;
end

if ~exist([directory,'templates',filesep,'labeling',filesep, ...
        'w',options.lc.general.parcellation,'.nii'],'file') ...
        || overwrite

    %% warp atlas into pre_tra-space:
    if ~exist([directory,'templates'],'dir')
        mkdir([directory,'templates']);
    end
    if ~exist([directory,'templates',filesep,'labeling'],'dir')
        mkdir([directory,'templates',filesep,'labeling']);
    end

    whichnormmethod=ea_whichnormmethod([directory]);
    switch whichnormmethod
        case ea_getantsnormfuns
            useinterp='GenericLabel';
            parc=ea_load_nii([ea_space(options,'labeling'),options.lc.general.parcellation,'.nii']);
            if length(unique(parc.img(:)))>500
                useinterp='NearestNeighbor'; % both GenericLabel and MultiLabel take ages on high dimensional parcellations.
            end
            ea_ants_apply_transforms(options, ...
                {[ea_space(options,'labeling'),options.lc.general.parcellation,'.nii']}, ...
                {[directory,'templates',filesep,'labeling',filesep,'w',options.lc.general.parcellation,'.nii']},...
                1,'','',useinterp);

        case ea_getfslnormfuns

            ea_fsl_apply_normalization(options, ...
                {[ea_space(options,'labeling'),options.lc.general.parcellation,'.nii']}, ...
                {[directory,'templates',filesep,'labeling',filesep,'w',options.lc.general.parcellation,'.nii']},...
                1,'','','nn');

        otherwise
            switch spm('ver')
                case 'SPM8'
                    matlabbatch{1}.spm.util.defs.comp{1}.def = {[directory,'y_ea_inv_normparams.nii']};
                    matlabbatch{1}.spm.util.defs.ofname = '';
                    matlabbatch{1}.spm.util.defs.fnames = {[ea_space(options,'labeling'),options.lc.general.parcellation,'.nii,1']};
                    matlabbatch{1}.spm.util.defs.savedir.saveusr = {[directory,'templates',filesep,'labeling',filesep]};
                    matlabbatch{1}.spm.util.defs.interp = 0;
                    spm_jobman('run',{matlabbatch});
                    clear matlabbatch
                case 'SPM12'
                    matlabbatch{1}.spm.util.defs.comp{1}.def = {[directory,'y_ea_inv_normparams.nii']};
                    matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = {[ea_space(options,'labeling'),options.lc.general.parcellation,'.nii']};
                    matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {[directory,'templates',filesep,'labeling',filesep]};
                    matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 0;
                    matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
                    matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
                    spm_jobman('run',{matlabbatch});
                    clear matlabbatch
            end
    end
end

[~,refname]=fileparts(reference);
[~,anatfname]=fileparts(options.prefs.prenii_unnormalized);

% For fMRI, the real reference image is 'meanrest_*.nii' rather than 'rrest_*.nii'
if strcmp(reference, ['r', options.prefs.rest])
    reference = ['hdmean', options.prefs.rest];
    % Re-calculate mean re-aligned image if not found
    if ~exist([directory, reference], 'file')
        if ~exist(['mean', options.prefs.rest],'file')
            ea_meanimage([directory, 'r', options.prefs.rest], ['mean', options.prefs.rest]);
        end
        ea_reslice_nii([directory,'mean', options.prefs.rest],[directory,reference],[0.7,0.7,0.7],0,0,1,[],[],0);
    end

end

if ~exist([directory,'templates',filesep,'labeling',filesep,refname,'w', ...
        options.lc.general.parcellation,'.nii'],'file') ...
        || overwrite

    % for this pair of approved coregistations, find out which method to use -
    % irrespective of the current selection in coregmethod.

    coregmethodsused=load([directory,'ea_coregmrmethod_applied.mat']);
    fn=fieldnames(coregmethodsused);
    for field=1:length(fn)
        if contains(fn{field},refname)
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

        transform = ea_coreg2images(options,[directory,options.prefs.prenii_unnormalized],...
            [directory,reference],...
            [directory,refname,'_',options.prefs.prenii_unnormalized],...
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

    ea_apply_coregistration([directory,reference], ...
        [directory,'templates',filesep,'labeling',filesep,'w',options.lc.general.parcellation,'.nii'], ...
        [directory,'templates',filesep,'labeling',filesep,refname,'w',options.lc.general.parcellation,'.nii'], ...
        transform, 'label');

    ea_gencheckregpair([directory,'templates',filesep,'labeling',filesep,refname,'w',options.lc.general.parcellation],...
        [directory,refname],...
        [directory,'checkreg',filesep,options.lc.general.parcellation,'2',refname,'.png']);
end
