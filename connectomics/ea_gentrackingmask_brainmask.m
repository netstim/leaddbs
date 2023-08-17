function ea_gentrackingmask_brainmask(options,threshold)
directory=[options.root,options.patientname,filesep];

copyfile([ea_space,'brainmask.nii.gz'],[directory,'brainmask.nii.gz']);
gunzip([directory,'brainmask.nii.gz']);
ea_delete([directory,'brainmask.nii.gz']);
ea_apply_normalization_tofile(options, {[directory,'brainmask.nii']}, {[directory,'brainmask.nii']}, 1, 0);

b0_anat = [ea_stripext(options.prefs.b0),'_',ea_stripext(options.prefs.prenii_unnormalized)];

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

    save([directory,'ea_coregmrmethod_applied.mat'],'-struct','coregmrmethod');
    docoreg = 1;
else
    if exist([directory,'ea_coregmrmethod_applied.mat'], 'file')
        coregmrmethod = load([directory,'ea_coregmrmethod_applied.mat']);
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
            save([directory,'ea_coregmrmethod_applied.mat'],'-struct','coregmrmethod');
            docoreg = 1;
        end
    else % 'ea_coregmrmethod_applied.mat' not found, redo coreg use method from GUI
        fprintf(['\nThe coregistration between b0 and anat has been approved.\n', ...
                 'But it is not clear which method has been used.\n', ...
                 'Will redo the coregistration using the method chosen from GUI.\n\n']);
        coregmrmethod.(b0_anat) = checkCoregMethod(options);
        save([directory,'ea_coregmrmethod_applied.mat'],'-struct','coregmrmethod');
        docoreg = 1;
    end
end

if docoreg
    if strcmp(coregmrmethod.(b0_anat), 'ANTsSyN')
        ea_ants_nonlinear_coreg([directory,options.prefs.prenii_unnormalized],...
            [directory,options.prefs.b0],...
            [directory,ea_stripext(options.prefs.b0), '2', options.prefs.prenii_unnormalized]);
        ea_delete([directory,ea_stripext(options.prefs.b0), '2', options.prefs.prenii_unnormalized]);
        ea_ants_apply_transforms(struct,[directory,'brainmask.nii'],...
            [directory,'trackingmask.nii'],0,[directory,options.prefs.b0],...
            [directory,ea_stripext(options.prefs.b0), '2', ea_stripext(options.prefs.prenii_unnormalized), 'InverseComposite.nii.gz'], 'Linear');
    else
        copyfile([directory,'brainmask.nii'],[directory,'cbrainmask.nii']);
        affinefile = ea_coregimages(options, ...
            [directory,options.prefs.prenii_unnormalized], ... % moving
            [directory,options.prefs.b0], ... % fix
            [directory,'c',options.prefs.prenii_unnormalized], ... % out
            {[directory,'cbrainmask.nii']}, ... % other
            1); % writeout transform

        % Rename saved transform for further use: canat2b0*.mat to anat2b0*.mat,
        % and b02canat*.mat to b02anat*.mat
        [~, anat] = ea_niifileparts(options.prefs.prenii_unnormalized);
        if ~isempty(affinefile)
            for i = 1:numel(affinefile)
                if ~strcmp(affinefile{i},strrep(affinefile{i}, ['c',anat], anat))
                    movefile(affinefile{i}, strrep(affinefile{i}, ['c',anat], anat));
                end
            end
        end

        movefile([directory,'cbrainmask.nii'],[directory,'trackingmask.nii']);
        delete([directory,'c',options.prefs.prenii_unnormalized]);
    end
else
    % Reuse approved coregistration
    if strcmp(coregmrmethod.(b0_anat), 'ANTsSyN')
        ea_ants_apply_transforms(struct,[directory,'brainmask.nii'],...
            [directory,'trackingmask.nii'],0,[directory,options.prefs.b0],...
            [directory,ea_stripext(options.prefs.b0), '2', ea_stripext(options.prefs.prenii_unnormalized), 'InverseComposite.nii.gz'],'Linear');
    else
        copyfile([directory,'brainmask.nii'],[directory,'cbrainmask.nii']);
        ea_apply_coregistration([directory,options.prefs.b0], [directory,'cbrainmask.nii'], ...
            [directory,'trackingmask.nii'], transform);
        ea_delete([directory,'cbrainmask.nii']);
    end
end

tr=ea_load_nii([options.root,options.patientname,filesep,'trackingmask.nii']);
if threshold
    tr.img=tr.img>0.1;
    tr.fname=[options.root,options.patientname,filesep,'ttrackingmask.nii'];
    ea_write_nii(tr);
end

% Transform anat to b0 to generate checkreg image
if strcmp(coregmrmethod.(b0_anat), 'ANTsSyN')
    ea_ants_apply_transforms(struct,[directory,options.prefs.prenii_unnormalized],...
            [directory,ea_stripext(options.prefs.b0), '_', options.prefs.prenii_unnormalized],0,[directory,options.prefs.b0],...
            [directory,ea_stripext(options.prefs.b0), '2', ea_stripext(options.prefs.prenii_unnormalized), 'InverseComposite.nii.gz'],'BSpline');
else
    ea_apply_coregistration([directory,options.prefs.b0], [directory,options.prefs.prenii_unnormalized], ...
                            [directory,b0_anat,'.nii'], strrep(options.coregmr.method, 'Hybrid SPM & ', ''));
end


function coregmethod = checkCoregMethod(options)
if strcmp(options.coregmr.method, 'ANTs') && options.coregb0.addSyN
    coregmethod = 'ANTsSyN';
else
    coregmethod = options.coregmr.method;
end
