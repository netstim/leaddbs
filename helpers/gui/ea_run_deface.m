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
                    ea_run_anonymize_pt(dirs{pt}, answ, options);               % run deface function with desired methode (deface or brain extract)
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

if ~strcmp(directory(end), filesep)
    directory=[directory, filesep];
end

% make export directory if it does not exist
exportfolder = fullfile(directory, 'export');

if ~exist(exportfolder, 'dir')
    mkdir(fullfile(directory, 'export'))
end


switch method
    
    % simple brain masking using the MNI ICBM 2009b brain mask
    case 'Mask'
        mask_file = fullfile(options.earoot, 'templates', 'space', ea_getspace, 'brainmask.nii.gz');
        
        % normalize MNI brainmask to patient space
        ea_apply_normalization(options, mask_file, fullfile(ea_getleadtempdir, 'brainmask_'))
        for fi = 1:length(presentfiles)
            ea_apply_normalization_tofile(options, )
        end
end

tempdir=ea_getleadtempdir;
copyfile([ea_space,'mask.nii.gz'],[tempdir,'mask.nii.gz']);
gunzip([tempdir,'mask.nii.gz']);
delete([tempdir,'mask.nii.gz']);
ea_apply_normalization_tofile(options,{[tempdir,'mask.nii']},{[directory,'brainmask.nii']},directory,1,0,[directory,presentfiles{1}]);

mask=ea_load_nii([directory,'brainmask.nii']);
for fi=1:length(presentfiles)
    
    
    % also deface raw pre-/postop files
    if exist([directory,'raw_',presentfiles{fi}],'file');
        defaceraw([directory,presentfiles{fi}],[directory,'raw_',presentfiles{fi}],directory,options);
    end
    
    % also deface raw postop-ct file
    if exist([directory,presentfiles{fi}(2:end)],'file') && strcmp(presentfiles{fi}(2:end),options.prefs.rawctnii_unnormalized);
        defaceraw([directory,presentfiles{fi}],[directory,presentfiles{fi}(2:end)],directory,options);
    end

    
    if exist([directory,presentfiles{fi}],'file');
        thisanat=ea_load_nii([directory,presentfiles{fi}]);
        thisanat.img(~logical(mask.img))=0;
        ea_write_nii(thisanat);
    end
end

function defaceraw(registered,rawfile,directory,options)
 % going back by backregistering here. We should change this to export
    % .mats by default when registering raw to anchor modality in the first
    % place, in the future.

    copyfile([directory,'brainmask.nii'],[directory,'rbrainmask.nii']); % create brainmask for raw file copy.
    ea_coreg2images(options,registered,rawfile,[directory,'tmp.nii'],{[directory,'rbrainmask.nii']},0,[],0);
    
    rmask=ea_load_nii([directory,'rbrainmask.nii']);
    rmask.img=logical(rmask.img);
    
    thisanat=ea_load_nii(rawfile);
    thisanat.img(~logical(rmask.img))=0;
    ea_write_nii(thisanat);



function presentfiles=stripstn(presentfiles)

[mem,ix]=ismember(presentfiles,{'anat_STN.nii','anat_RN.nii','anat_GPi.nii','anat_GPe.nii'});
presentfiles=presentfiles(~mem);

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
                        
                        
                        
                        
                        
                        
                        
                        
                        
