function ea_warp_parcellation(reference,options)
directory=[options.root,options.patientname,filesep];

% Regenerate the template or not
if isfield(options, 'overwriteapproved') && options.overwriteapproved
    overwrite = 1;
else
    overwrite = 0;
end

% Require the subject's active normalization; do not estimate a second one.
ea_gettransformfiles(options);

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

    ea_apply_normalization_tofile(options, ...
        {[ea_space(options,'labeling'),options.lc.general.parcellation,'.nii']}, ...
        {[directory,'templates',filesep,'labeling',filesep,'w',options.lc.general.parcellation,'.nii']}, ...
        1, 0);
end

[~,refname]=fileparts(reference);
[~,anatfname]=fileparts(options.prefs.prenii_unnormalized);

% BIDS FIX: Determine modality and use short names for template files
% This avoids too-long BIDS filenames in template directory
rest_r = ea_prependFilename(options.prefs.rest, 'r');
isFMRI = false;
isDWI = false;

% Check if this is fMRI reference
if strcmp(reference, rest_r)
    rest_hdmean = ea_prependFilename(options.prefs.rest, 'hdmean');
    reference = rest_hdmean;
    isFMRI = true;
    % Re-calculate mean re-aligned image if not found
    if ~exist([directory, reference], 'file')
        rest_mean = ea_prependFilename(options.prefs.rest, 'mean');
        if ~exist([directory, rest_mean],'file')
            ea_meanimage([directory, rest_r], rest_mean);
        end
        ea_reslice_nii([directory, rest_mean],[directory,reference],[0.7,0.7,0.7],0,0,1,[],[],3);
    end
    refname = 'rest';
% Check if this is b0 reference (DWI)
elseif contains(reference, 'b0')
    isDWI = true;
    refname = 'b0';
end

if ~exist([directory,'templates',filesep,'labeling',filesep,refname,'w', ...
        options.lc.general.parcellation,'.nii'],'file') ...
        || overwrite

    coregLogFile = fullfile(directory, 'coregistration', 'log', 'ea_coregmrmethod_applied.mat');
    if ~isfile(coregLogFile)
        coregLogFile = fullfile(directory, 'ea_coregmrmethod_applied.mat');
    end
    if exist(coregLogFile, 'file')
        % Check coreg method used
        coregmrmethod = load(coregLogFile);
        images = fieldnames(coregmrmethod);
        if any(contains(images, refname))
            options.coregmr.method = coregmrmethod.(images{contains(images, refname)});
            fprintf(['Coregistration already done using ', options.coregmr.method, '.\n'...
                'Will try to use the existing transform.\n']);
        end
    else
        fprintf(['Unable to determine the coregistration method used.\n', ...
                 'Will redo the coregistration using the method chosen from GUI.\n']);
        % BIDS FIX: Sanitize field name for MATLAB (max 63 chars)
        fieldName = matlab.lang.makeValidName([refname, '_', anatfname], 'ReplacementStyle', 'delete');
        if length(fieldName) > 63
            fieldName = fieldName(1:63);
        end
        coregmrmethod.(fieldName) = checkCoregMethod(options);
        coregLogDir = fullfile(directory, 'coregistration', 'log');
        if ~isfolder(coregLogDir), mkdir(coregLogDir); end
        save(fullfile(coregLogDir, 'ea_coregmrmethod_applied.mat'), '-struct', 'coregmrmethod');
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

        % BIDS FIX: Build output name correctly (get basename of anat)
        [~, anatBaseName] = fileparts(ea_stripext(options.prefs.prenii_unnormalized));
        [~, refBaseName] = fileparts(refname);
        outputName = [refBaseName, '_', anatBaseName, '.nii'];
        
        % For BIDS paths, construct full output path
        if contains(reference, filesep)
            [refDir, ~] = fileparts(reference);
            outputPath = fullfile(refDir, outputName);
        else
            outputPath = outputName;
        end

        if strcmp(options.coregmr.method, 'ANTsSyN') % ANTs nonlinear case, only for b0 coreg
            ea_ants_nonlinear_coreg([directory, options.prefs.prenii_unnormalized],...
                [directory, reference],...
                [directory, outputPath]);
            ea_delete([directory, outputPath]);
            transform = [directory, refBaseName, '2', anatBaseName, 'InverseComposite.nii.gz'];
        else
            transform = ea_coregimages(options,[directory,options.prefs.prenii_unnormalized],...
                [directory,reference],...
                [directory,outputPath],...
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

    ea_gencheckregpair([directory,'templates',filesep,'labeling',filesep,refname,'w',options.lc.general.parcellation,'.nii'],...
        [directory,reference],...
        [directory,'checkreg',filesep,options.lc.general.parcellation,'2',refname,'.png']);
end


function coregmethod = checkCoregMethod(options)
if strcmp(options.coregmr.method, 'ANTs') && options.coregb0.addSyN
    coregmethod = 'ANTsSyN';
else
    coregmethod = options.coregmr.method;
end
