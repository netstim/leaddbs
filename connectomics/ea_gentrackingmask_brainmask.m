function ea_gentrackingmask_brainmask(options,threshold)
directory=[options.root,options.patientname,filesep];

% BIDS: Place masks in preprocessing/dwi/ subdirectory
maskDir = fullfile(directory, 'preprocessing', 'dwi');
if ~isfolder(maskDir)
    mkdir(maskDir);
end

% BIDS FIX: Generate brainmask directly from native anat image
% This avoids dependency on potentially outdated normalization transformations
brainmaskFile = fullfile(maskDir, 'brainmask.nii');

if ~exist(brainmaskFile, 'file') || (isfield(options, 'overwriteapproved') && options.overwriteapproved)
    fprintf('Generating brain mask from native anatomical image...\n');
    
    % Use SPM segmentation to get GM, WM, CSF
    anatFile = [directory, options.prefs.prenii_unnormalized];
    
    % Run SPM12 segmentation if tissue classes don't exist
    [anatDir, anatName] = fileparts(anatFile);
    c1File = fullfile(anatDir, ['c1', anatName, '.nii']);
    c2File = fullfile(anatDir, ['c2', anatName, '.nii']);
    c3File = fullfile(anatDir, ['c3', anatName, '.nii']);
    
    if ~exist(c1File, 'file') || ~exist(c2File, 'file') || ~exist(c3File, 'file')
        fprintf('Running SPM segmentation...\n');
        matlabbatch{1}.spm.spatial.preproc.channel.vols = {[anatFile, ',1']};
        matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
        matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
        matlabbatch{1}.spm.spatial.preproc.channel.write = [0 0];
        
        % Tissue classes: GM, WM, CSF
        for ti = 1:3
            matlabbatch{1}.spm.spatial.preproc.tissue(ti).tpm = {[spm('Dir'), '/tpm/TPM.nii,', num2str(ti)]};
            matlabbatch{1}.spm.spatial.preproc.tissue(ti).ngaus = ti;
            matlabbatch{1}.spm.spatial.preproc.tissue(ti).native = [1 0];  % Native space
            matlabbatch{1}.spm.spatial.preproc.tissue(ti).warped = [0 0];
        end
        
        % Other tissues (skull, soft tissue, air)
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
    end
    
    % Combine GM + WM + CSF to create brain mask
    c1 = ea_load_nii(c1File);
    c2 = ea_load_nii(c2File);
    c3 = ea_load_nii(c3File);
    
    % Create binary mask (threshold at 0.5)
    brainmask = c1;
    brainmask.img = (c1.img + c2.img + c3.img) > 0.5;
    brainmask.fname = brainmaskFile;
    brainmask.dt = [2 0];  % uint8
    ea_write_nii(brainmask);
    
    fprintf('Brain mask generated successfully.\n');
else
    fprintf('Brain mask already exists, skipping generation.\n');
end

% Create valid MATLAB field name from b0 and anat filenames (BIDS-safe)
b0_anat = matlab.lang.makeValidName([ea_stripext(options.prefs.b0),'_',ea_stripext(options.prefs.prenii_unnormalized)]);

if exist([directory,'ea_coreg_approved.mat'], 'file')
    isapproved = load([directory,'ea_coreg_approved.mat']);
    if isfield(isapproved, b0_anat)
        isapproved = isapproved.(b0_anat);
    else
        isapproved = 0;
    end
else
    isapproved = 0;
end

if ~isapproved || isfield(options, 'overwriteapproved') && options.overwriteapproved
    fprintf(['\nCoregistration between b0 and anat not done/not approved yet or overwrite mode enabled!\n', ...
             'Will do the coregistration using the method chosen from GUI.\n\n'])

    coregmrmethod.(b0_anat) = checkCoregMethod(options);

    % BIDS: Save coregistration method log in coregistration/log/
    coregLogDir = fullfile(directory, 'coregistration', 'log');
    if ~exist(coregLogDir, 'dir')
        mkdir(coregLogDir);
    end
    save(fullfile(coregLogDir, 'ea_coregmrmethod_applied.mat'),'-struct','coregmrmethod');
    docoreg = 1;
else
    coregLogFile = fullfile(directory, 'coregistration', 'log', 'ea_coregmrmethod_applied.mat');
    % Fallback: try root for backwards compatibility
    if ~exist(coregLogFile, 'file')
        coregLogFile = fullfile(directory, 'ea_coregmrmethod_applied.mat');
    end
    if exist(coregLogFile, 'file')
        coregmrmethod = load(coregLogFile);
        if isfield(coregmrmethod, b0_anat) && ~isempty(coregmrmethod.(b0_anat))
            if strcmp(coregmrmethod.(b0_anat), 'ANTsSyN')
                transform = [ea_stripext(options.prefs.b0), '2', ea_stripext(options.prefs.prenii_unnormalized), ...
                             '(Inverse)?Composite\.nii\.gz$'];
            else
                transform = [ea_stripext(options.prefs.prenii_unnormalized), '2', ea_stripext(options.prefs.b0), ...
                             '_', lower(strrep(coregmrmethod.(b0_anat), 'Hybrid SPM & ', '')), '\d*\.(mat|h5)$'];
            end

            transform = ea_regexpdir(directory, transform, 0);

            if isempty(transform) % Approved but the transformation file not found
                fprintf(['\nThe coregistration between b0 and anat using ''%s'' has been approved.\n', ...
                         'But the tranformation file doesn''t exist!\n', ...
                         'Will redo the coregistration using ''%s''.\n\n'], ...
                         coregmrmethod.(b0_anat), coregmrmethod.(b0_anat));
                docoreg = 1;
            else % Approved and the transformation file found
                transform = transform{end};
                fprintf(['\nThe coregistration between b0 and anat using ''%s'' has been approved.\n', ...
                         'Will use this transform to generate the fiber tracking mask:\n', ...
                         '%s\n\n'], coregmrmethod.(b0_anat), transform);
                docoreg = 0;
            end
            options.coregmr.method = checkCoregMethod(options);
        else % 'ea_coregmrmethod_applied.mat' found but method not set, redo coreg use method from GUI
            fprintf(['\nThe coregistration between b0 and anat has been approved.\n', ...
                     'But it is not clear which method has been used.\n', ...
                     'Will redo the coregistration using the method chosen from GUI.\n\n']);
            coregmrmethod.(b0_anat) = checkCoregMethod(options);
            coregLogDir = fullfile(directory, 'coregistration', 'log');
            if ~exist(coregLogDir, 'dir')
                mkdir(coregLogDir);
            end
            save(fullfile(coregLogDir, 'ea_coregmrmethod_applied.mat'),'-struct','coregmrmethod');
            docoreg = 1;
        end
    else % 'ea_coregmrmethod_applied.mat' not found, redo coreg use method from GUI
        fprintf(['\nThe coregistration between b0 and anat has been approved.\n', ...
                 'But it is not clear which method has been used.\n', ...
                 'Will redo the coregistration using the method chosen from GUI.\n\n']);
        coregmrmethod.(b0_anat) = checkCoregMethod(options);
        coregLogDir = fullfile(directory, 'coregistration', 'log');
        if ~exist(coregLogDir, 'dir')
            mkdir(coregLogDir);
        end
        save(fullfile(coregLogDir, 'ea_coregmrmethod_applied.mat'),'-struct','coregmrmethod');
        docoreg = 1;
    end
end

% Initialize affinefile
affinefile = [];

if docoreg
    if strcmp(coregmrmethod.(b0_anat), 'ANTsSyN')
        % BIDS: Place ANTs transformation files in coregistration/transformations/
        xfmDir = fullfile(directory, 'coregistration', 'transformations');
        if ~isfolder(xfmDir)
            mkdir(xfmDir);
        end
        
        [~, b0Base] = fileparts(options.prefs.b0);
        [~, anatBase] = fileparts(options.prefs.prenii_unnormalized);
        xfmBase = fullfile(xfmDir, [b0Base, '2', anatBase]);
        
        ea_ants_nonlinear_coreg([directory,options.prefs.prenii_unnormalized],...
            [directory,options.prefs.b0], xfmBase);
        ea_delete(xfmBase);
        ea_ants_apply_transforms(struct, fullfile(maskDir, 'brainmask.nii'),...
            fullfile(maskDir, 'trackingmask.nii'), 0, [directory,options.prefs.b0],...
            [xfmBase, 'InverseComposite.nii.gz'], 'Linear');
        affinefile = []; % ANTs doesn't use affinefile
    else
        copyfile(fullfile(maskDir, 'brainmask.nii'), fullfile(maskDir, 'cbrainmask.nii'));
        
        % Construct output path with 'c' prefix on basename (BIDS-safe)
        [anatdir, anatname, anatext] = fileparts([directory,options.prefs.prenii_unnormalized]);
        outAnat = fullfile(anatdir, ['c', anatname, anatext]);
        
        affinefile = ea_coregimages(options, ...
            [directory,options.prefs.prenii_unnormalized], ... % moving
            [directory,options.prefs.b0], ... % fixed
            outAnat, ... % out
            {fullfile(maskDir, 'cbrainmask.nii')}, ... % other
            1); % writeout transform
        
        % Rename saved transform for further use: canat2b0*.mat to anat2b0*.mat,
        % and b02canat*.mat to b02anat*.mat
        [~, anat] = ea_niifileparts(options.prefs.prenii_unnormalized);
        fprintf('DEBUG: affinefile from ea_coregimages: ');
        disp(affinefile);
        fprintf('DEBUG: anat basename: %s\n', anat);
        
        if ~isempty(affinefile) && iscell(affinefile) && ~isempty(affinefile{1})
            for i = 1:numel(affinefile)
                oldPath = affinefile{i};
                newPath = strrep(oldPath, ['c',anat], anat);
                fprintf('DEBUG: Checking transform %d: %s\n', i, oldPath);
                if ~strcmp(oldPath, newPath)
                    if exist(oldPath, 'file')
                        fprintf('DEBUG: Moving %s -> %s\n', oldPath, newPath);
                        movefile(oldPath, newPath);
                        affinefile{i} = newPath;
                    else
                        fprintf('DEBUG: WARNING - file does not exist: %s\n', oldPath);
                    end
                else
                    fprintf('DEBUG: No rename needed (no c-prefix found)\n');
                end
            end
        else
            fprintf('DEBUG: ERROR - affinefile is empty or invalid!\n');
        end
        fprintf('DEBUG: Final affinefile: ');
        disp(affinefile);

        movefile(fullfile(maskDir, 'cbrainmask.nii'), fullfile(maskDir, 'trackingmask.nii'));
        ea_delete(outAnat);
    end
else
    % Reuse approved coregistration
    if strcmp(coregmrmethod.(b0_anat), 'ANTsSyN')
        % BIDS: ANTs transformation files in coregistration/transformations/
        xfmDir = fullfile(directory, 'coregistration', 'transformations');
        [~, b0Base] = fileparts(options.prefs.b0);
        [~, anatBase] = fileparts(options.prefs.prenii_unnormalized);
        xfmBase = fullfile(xfmDir, [b0Base, '2', anatBase]);
        
        ea_ants_apply_transforms(struct, fullfile(maskDir, 'brainmask.nii'),...
            fullfile(maskDir, 'trackingmask.nii'), 0, [directory,options.prefs.b0],...
            [xfmBase, 'InverseComposite.nii.gz'],'Linear');
        affinefile = []; % ANTs doesn't use affinefile
    else
        copyfile(fullfile(maskDir, 'brainmask.nii'), fullfile(maskDir, 'cbrainmask.nii'));
        ea_apply_coregistration([directory,options.prefs.b0], fullfile(maskDir, 'cbrainmask.nii'), ...
            fullfile(maskDir, 'trackingmask.nii'), transform);
        ea_delete(fullfile(maskDir, 'cbrainmask.nii'));
        % Store transform for checkreg generation
        affinefile = {transform};
    end
end

tr=ea_load_nii(fullfile(maskDir, 'trackingmask.nii'));
if threshold
    tr.img=tr.img>0.1;
    tr.fname=fullfile(maskDir, 'ttrackingmask.nii');
    ea_write_nii(tr);
    % Legacy: also save in root for DSI Studio when not BIDS
    if ~contains(directory, 'derivatives') && ~contains(directory, 'leaddbs')
        tr.fname = [directory, 'ttrackingmask.nii'];
        ea_write_nii(tr);
    end
end

% Transform anat to b0 to generate checkreg image (BIDS: place in coregistration/dwi/)
coregDwiDir = fullfile(directory, 'coregistration', 'dwi');
if ~isfolder(coregDwiDir)
    mkdir(coregDwiDir);
end

if strcmp(coregmrmethod.(b0_anat), 'ANTsSyN')
    % BIDS: ANTs transformation files in coregistration/transformations/
    xfmDir = fullfile(directory, 'coregistration', 'transformations');
    [~, b0Base] = fileparts(options.prefs.b0);
    [~, anatBase] = fileparts(options.prefs.prenii_unnormalized);
    xfmBase = fullfile(xfmDir, [b0Base, '2', anatBase]);
    
    ea_ants_apply_transforms(struct,[directory,options.prefs.prenii_unnormalized],...
            fullfile(coregDwiDir, [b0_anat, '.nii']), 0, [directory,options.prefs.b0],...
            [xfmBase, 'InverseComposite.nii.gz'],'BSpline');
else
    % Use the transformation file (SPM/FSL)
    % Pass the transformation file directly - ea_apply_coregistration now handles _spm.mat format
    if ~isempty(affinefile) && iscell(affinefile) && ~isempty(affinefile{1})
        ea_apply_coregistration([directory,options.prefs.b0], [directory,options.prefs.prenii_unnormalized], ...
                                fullfile(coregDwiDir, [b0_anat, '.nii']), affinefile{1});
    else
        warning('No transformation file found, skipping checkreg image generation');
    end
end


function coregmethod = checkCoregMethod(options)
if strcmp(options.coregmr.method, 'ANTs') && options.coregb0.addSyN
    coregmethod = 'ANTsSyN';
else
    coregmethod = options.coregmr.method;
end
