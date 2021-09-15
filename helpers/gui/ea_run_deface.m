function ea_run_deface(hobj, evt, handles)

answ=questdlg(['Warning: this method only works for the MNI ICBM 2009b NLIN ASYM Space!' sprintf('\n\n') 'Do you want to run brain masking or defacing to anonymize patient data?'], ...
    'Anonymize base files', ...
    'Mask','Deface','Cancel', 'Cancel');

if ~strcmp(answ,'Cancel')
    if strcmp( ea_getspace, 'MNI_ICBM_2009b_NLIN_ASYM')      % check if space is actually MNI ICBM 2009b
        dirs=getappdata(handles.leadfigure, 'uipatdir');                 % get current list of patients loaded in the gui
        
        for pt=1:length(dirs)
            
            [whichnormmethod, template] = ea_whichnormmethod(dirs{pt});
            options=ea_getptopts(dirs{pt});
            switch whichnormmethod                   % before running, check whether patient has been normalized or not
                case ''
                    fprintf('\nAnonymization not possible for patient %s, because normalization was not performed!\nPlease perform normalization first.\n', options.patientname);
                otherwise
                    fprintf('\nRunning anonymization with the option ''%s''...\n', answ);
                    ea_run_anonymize_pt(dirs{pt}, answ, options);               % run deface function with desired methode (deface or brain extract)
                    fprintf('** Process Done.\n')
            end
            
        end
    else
        fprintf('Space is %s, please switch to MNI_ICBM_2009b_NLIN_ASYM space to enable anonymization of patient data!\n', ea_getspace)
    end
end

    function  ea_run_anonymize_pt(directory, method, options)
        
        [~,presentfiles] = ea_assignpretra(options);     % get present pre-op files
        presentfiles=stripstn(presentfiles);                    % strip stn files if present
        
        presentfiles=ea_add_postop_files(presentfiles, options);   % add postop files
        
        % create cell without the endings
        for fi = 1:length(presentfiles)
            [p, f, e] = fileparts(presentfiles{fi});
            presentfiles_fnames{fi, 1} = f;
        end
        
        if ~strcmp(directory(end), filesep)
            directory=[directory, filesep];
        end
        
        % make export directory if it does not exist
        exportfolder = fullfile(directory, 'export', 'anonymized');
        
        if ~exist(exportfolder, 'dir')
            mkdir(exportfolder);
        end
        
        switch method
            % simple brain masking using the MNI ICBM 2009b brain mask
            case 'Mask'
                mask_file = fullfile(options.earoot, 'templates', 'space', ea_getspace, 'brainmask.nii.gz');
                anonymized_file_suffix = '_brain_mask';
                mask_file_pt_space = fullfile(exportfolder, 'mask_brain.nii.gz');
                % defacing by defaced MNI ICBM 2009b brain mask
            case 'Deface'
                mask_file = fullfile(options.earoot, 'templates', 'space', ea_getspace, 'brainmask_defaced.nii.gz');
                anonymized_file_suffix = '_defaced';
                mask_file_pt_space = fullfile(exportfolder, 'mask_deface.nii.gz');
        end
        
        % normalize MNI brainmask to patient space
        ea_apply_normalization_tofile(options, mask_file, mask_file_pt_space, 1, 0, fullfile(directory, presentfiles{1}));
        
        % now go through all the files and multiply them with the mask to
        % create brain extracted images
        for fi = 1:length(presentfiles)
            if exist(fullfile(directory, presentfiles{fi}), 'file')
                applyMask(fullfile(directory, presentfiles{fi}), fullfile(exportfolder, [presentfiles_fnames{fi}, anonymized_file_suffix, '.nii']), mask_file_pt_space);
                gzip(fullfile(exportfolder, [presentfiles_fnames{fi}, anonymized_file_suffix, '.nii']));
                delete(fullfile(exportfolder, [presentfiles_fnames{fi}, anonymized_file_suffix, '.nii']));
            else
                fprintf('Did not find file %s, skipping.\n', fullfile(directory, presentfiles{fi}));
            end
        end
    end
end

function applyMask(infilename, outfilename, mask_filename)
mask_pt_space = ea_load_nii(mask_filename);         % load mask .nii
mask_pt_space.img=logical(mask_pt_space.img);    % convert to logical

anatfile = ea_load_nii(infilename);
outfile = anatfile;
outfile.img(~logical(mask_pt_space.img)) =0;
outfile.fname=outfilename;
ea_write_nii(outfile);
fprintf('Written anonymized image: %s\n', outfilename);
end

function presentfiles=stripstn(presentfiles)

[mem, ~]=ismember(presentfiles,{'anat_STN.nii','anat_RN.nii','anat_GPi.nii','anat_GPe.nii'});
presentfiles=presentfiles(~mem);
end

function presentfiles=ea_add_postop_files(presentfiles,options)
switch options.modality
    case 1 % MR
        % postop MR
        presentfiles{end+1}=options.prefs.tranii_unnormalized;
        presentfiles{end+1}=options.prefs.cornii_unnormalized;
        presentfiles{end+1}=options.prefs.sagnii_unnormalized;
    case 2 % CT
        % postop CT
        presentfiles{end+1}=options.prefs.tp_ctnii_coregistered;
end
end







