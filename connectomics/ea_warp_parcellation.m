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
        ea_reslice_nii([directory,'mean', options.prefs.rest],[directory,reference],[0.7,0.7,0.7],0,0,1,[],[],3);
    end
end

if ~exist([directory,'templates',filesep,'labeling',filesep,refname,'w', ...
        options.lc.general.parcellation,'.nii'],'file') ...
        || overwrite

    if exist([directory,'ea_coregmrmethod_applied.mat'], 'file')
        % Check coreg method used
        coregmrmethod = load([directory,'ea_coregmrmethod_applied.mat']);
        images = fieldnames(coregmrmethod);
        if any(contains(images, refname))
            options.coregmr.method = coregmrmethod.(images{contains(images, refname)});
            fprintf(['Coregistration already done using ', options.coregmr.method, '.\n'...
                'Will try to use the existing transform.\n']);
        end
    else
        fprintf(['Unable to determine the coregistration method used.\n', ...
                 'Will redo the coregistration using the method chosen from GUI.\n']);
        coregmrmethod.([refname, '_', anatfname]) = checkCoregMethod(options);
        save([directory,'ea_coregmrmethod_applied.mat'],'-struct','coregmrmethod');
    end

    % Disable Hybrid coregistration
    options.coregmr.method = strrep(options.coregmr.method, 'Hybrid SPM & ', '');

    % Check if the corresponding transform already exists
    if strcmp(options.coregmr.method, 'ANTsSyN')
        transform = [refname, '2', anatfname, 'InverseComposite\.nii\.gz$'];
    else
        transform = [anatfname, '2', refname, '_', lower(options.coregmr.method), '\d*\.(mat|h5)$'];
    end

    transform = ea_regexpdir(directory, transform, 0);

    if numel(transform) == 0 || overwrite
        if numel(transform) == 0
            warning('Transformation not found! Running coregistration now!');
        end

        if strcmp(options.coregmr.method, 'ANTsSyN') % ANTs nonlinear case, only for b0 coreg
            ea_ants_nonlinear_coreg([directory, options.prefs.prenii_unnormalized],...
                [directory, reference],...
                [directory, refname, '2', options.prefs.prenii_unnormalized]);
            ea_delete([directory, refname, '2', options.prefs.prenii_unnormalized]);
            transform = [directory, refname, '2', anatfname, 'InverseComposite.nii.gz'];
        else
            transform = ea_coregimages(options,[directory,options.prefs.prenii_unnormalized],...
                [directory,reference],...
                [directory,refname,'_',options.prefs.prenii_unnormalized],...
                [],1,[],1);

            % Fix transformation names, replace 'mean' by 'r' for fMRI
            if strcmp(reference, ['mean', options.prefs.rest])
                cellfun(@(f) movefile(f, strrep(f, 'mean', 'r')), transform);
                transform = strrep(transform, 'mean', 'r');
            end
            transform = transform{1}; % Forward transformation
        end
    else
        if numel(transform) > 1
            warning(['Multiple transformations of the same type found! ' ...
                'Will use the last one:\n%s'], transform{end});
        end
        transform = transform{end};
    end

    if strcmp(options.coregmr.method, 'ANTsSyN') % ANTs nonlinear case, only for b0 coreg
        ea_ants_apply_transforms(struct,...
            [directory,'templates',filesep,'labeling',filesep,'w',options.lc.general.parcellation,'.nii'],...
            [directory,'templates',filesep,'labeling',filesep,refname,'w',options.lc.general.parcellation,'.nii'],...
            0,[directory,reference],...
            transform,'NearestNeighbor');
    else
        ea_apply_coregistration([directory,reference], ...
            [directory,'templates',filesep,'labeling',filesep,'w',options.lc.general.parcellation,'.nii'], ...
            [directory,'templates',filesep,'labeling',filesep,refname,'w',options.lc.general.parcellation,'.nii'], ...
            transform, 'label');
    end

    ea_gencheckregpair([directory,'templates',filesep,'labeling',filesep,refname,'w',options.lc.general.parcellation],...
        [directory,refname],...
        [directory,'checkreg',filesep,options.lc.general.parcellation,'2',refname,'.png']);
end


function coregmethod = checkCoregMethod(options)
if strcmp(options.coregmr.method, 'ANTs') && options.coregb0.addSyN
    coregmethod = 'ANTsSyN';
else
    coregmethod = options.coregmr.method;
end
