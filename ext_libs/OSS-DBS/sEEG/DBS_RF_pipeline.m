cd([ea_getearoot, '/ext_libs/OSS-DBS/sEEG'])

%% Inputs
subjs = {};

anchor = 'T1';
removeElectrodeFromRF = false;

space = 'native';
BIDS_format = false;
cleanup = true; % otherwise some large .nii files might remain

Segment_SVD = false;
regenerate_segmask = false;

warp2MNI = true;   % warp from the anchor space to MNI, we assume ANTs
if strcmp(space,'mni')
    % not needed in this case
    warp2MNI = false;
end

flip_lh2rh = true;
reslice2ROI = true; % if true, all VTRs are stored in the same voxel space defined by their bounding box thresholded at ROI_RF_threshold
ROI_RF_threshold = 0.1; % because we solve for -1V
vox_size = [0.75,0.75,0.75];

% hardcoded for now
VTR = true;             % if true, computes recording fields (aka volume of tissue recorded)
StimMode = 'CC';       % CC - current-controlled, VC - voltage-controlled

% auto-definitions
if strcmp(space, 'native')
    image_folder = 'coregistration';
    stimFolderSfx = fullfile('stimulations','native','BSS');
else
    image_folder = 'normalization';
    stimFolderSfx = fullfile('stimulations','MNI152NLin2009bAsym','BSS');
end

%% compute RFs for specified patients

if reslice2ROI && (warp2MNI || strcmp(space,'mni'))
    % in this case, put RFs in the common voxel space across subjects!
    images2ROI = {};
end  % otherwise this will be done within the subject

for subj_inx = 1:size(subjs,2)

    % get the imaging
    if BIDS_format
        all_images = dir_without_dots(fullfile(subjs{1,subj_inx},image_folder,'anat'));
        allNames = {all_images.name};
        logicalIndex = contains(allNames, 'T2starw.nii') | contains(allNames, 'T1w.nii') | contains(allNames, 'T2w.nii') | contains(allNames, 'FLAIR.nii');
        images = all_images(logicalIndex);
        % re-order to have T2w and T1w before FLAIR
        images = flipud(images);
    
        if strcmp(anchor,'T2') || strcmp(anchor,'t2')
            idx = find(contains({images.name}, 'T2w.nii'));
        elseif strcmp(anchor,'T1') || strcmp(anchor,'t1')
            idx = find(contains({images.name}, 'T1w.nii'));
        elseif strcmp(anchor,'FLAIR') || strcmp(anchor,'flair')
            idx = find(contains({images.name}, 'FLAIR.nii'));  
        else
            warn_msg = sprintf("The anchor modality %s is not supported, please use T1, T2 or FLAIR",anchor);
            ea_warndlg(warn_msg)
        end
    else
        all_images = dir_without_dots(fullfile(subjs{1,subj_inx}));
        allNames = {all_images.name};
        if strcmp(space, 'native')
            logicalIndex = strcmp(allNames, 'anat_t1.nii') | strcmp(allNames, 'anat_t2.nii') | strcmp(allNames, 'anat_t2star.nii') | strcmp(allNames, 'anat_flair.nii');
        else
            logicalIndex = strcmp(allNames, 'glanat.nii') | strcmp(allNames, 'glanat_t2.nii') | strcmp(allNames, 'glanat_t2star.nii') | strcmp(allNames, 'glanat_flair.nii');
        end

        images = all_images(logicalIndex);
        % re-order to have T2w and T1w before FLAIR
        images = flipud(images);      

        if strcmp(space,'mni')
            if strcmp(anchor,'T2') || strcmp(anchor,'t2')
                idx = find(contains({images.name}, 'glanat_t2.nii'));
            elseif strcmp(anchor,'T1') || strcmp(anchor,'t1')
                idx = find(contains({images.name}, 'glanat.nii'));
            elseif strcmp(anchor,'FLAIR') || strcmp(anchor,'flair')
                idx = find(contains({images.name}, 'glanat_flair.nii'));  
            else
                warn_msg = sprintf("The anchor modality %s is not supported, please use T1, T2 or FLAIR",anchor);
                ea_warndlg(warn_msg)
            end    
        else
            if strcmp(anchor,'T2') || strcmp(anchor,'t2')
                idx = find(contains({images.name}, 'anat_t2.nii'));
            elseif strcmp(anchor,'T1') || strcmp(anchor,'t1')
                idx = find(contains({images.name}, 'anat_t1.nii'));
            elseif strcmp(anchor,'FLAIR') || strcmp(anchor,'flair')
                idx = find(contains({images.name}, 'anat_flair.nii'));  
            else
                warn_msg = sprintf("The anchor modality %s is not supported, please use T1, T2 or FLAIR",anchor);
                ea_warndlg(warn_msg)
            end

            if isempty(idx)
                if strcmp(anchor,'T1') || strcmp(anchor,'t1')
                    idx = find(contains({images.name}, 'anat_t2.nii'));
                else
                    idx = find(contains({images.name}, 'anat_t1.nii'));
                end
            end

        end
    end


    anchorImage = fullfile(images(idx).folder,images(idx).name);

    % get the transform (if necessary)
    if warp2MNI == true
        if BIDS_format
            all_transforms = dir_without_dots(fullfile(subjs{1,subj_inx},'normalization','transformations'));
            allNames = {all_transforms.name};
            % from native to MNI, you need the opposite warp, i.e. from-MNI152NLin2009bAsym_to-anchorNative
             % IMPORTANT: make sure ants is mentioned in the name and use single quote
            logicalIndex = contains(allNames, 'from-MNI152NLin2009bAsym_to-anchorNative');
            transform_file = all_images(logicalIndex);
            if size(transform_file,1) ~= 1
                ea_warndlg("One (only) transform from-MNI152NLin2009bAsym_to-anchorNative should be present")
            end
            transform = fullfile(transform_file(1).folder,transform_file(1).name);
        else
            transform = [subjs{1,subj_inx},filesep,'glanatInverseComposite.nii.gz'];
            if ~isfile(transform)
                transform = [path2patient,filesep,'glanatInverseComposite.mat'];
            end
        end
    end

    % create a segmask (if not available)
    ea_mkdir(fullfile(subjs{1,subj_inx},stimFolderSfx))
    segmaskFile = fullfile(subjs{1,subj_inx},stimFolderSfx,'segmask.nii');
    
    if ~isfile(segmaskFile) || regenerate_segmask
        if Segment_SVD
            ea_svd_segmentation(images, anchorImage)
        else
            % segment the anchor modality
            % TBD: multimodal segmentation with SynthSeg
            [workingDir,~,~] = fileparts(anchorImage); 
            multisegment = [workingDir,filesep,'segmask_synthSeg.nii'];
            if ~isfile(segmaskFile)
                if ~isfile(multisegment)
                    ea_synthseg(anchorImage, multisegment)
                end
                ea_convert_synthSeg2segmask(multisegment, segmaskFile);
            end
        end
    end

    % Set OSS-DBS python path
    env = ea_conda_env('OSS-DBSv2');
    ea_checkOSSDBSInstallv2(env);

    % call a python script to compute stim. volumes
    OSS4DBS_RFs_script = [ea_getearoot, 'ext_libs/OSS-DBS/sEEG/run_OSS4DBS_RFs.py'];
    
    % prepare stim_json
    stim_json.path2subject = subjs{1,subj_inx};
    stim_json.path2segmask = segmaskFile;
    stim_json.space = space;
    stim_json.stimFolder = fullfile(subjs{1,subj_inx},stimFolderSfx);
    stim_json.StimMode = StimMode;
    stim_json.remove_el = removeElectrodeFromRF;
    DictPath = fullfile(stim_json.stimFolder,'RF_Input.json');

    % get the electrode reconstructions
    [~,subj_label,~] = fileparts(subjs{1,subj_inx});

    if BIDS_format
        [~,subj_label,~] = fileparts(subjs{1,subj_inx});
        recon = [subjs{1,subj_inx},filesep,'reconstruction',filesep,subj_label,'_desc-reconstruction.mat'];
    else
        recon = [subjs{1,subj_inx},filesep,'ea_reconstruction.mat'];
    end

    load(recon, 'reco');

    if isfield(reco,'electrode')
        stim_json.el_type = reco.electrode(1).dbs.elmodel; % assuming same electrode type
        eleNum = size(reco.electrode,2); % Number of electrodes
    else
        stim_json.el_type = reco.props(1).elmodel;
        eleNum = size(reco.props,2); % Number of electrodes
    end

    settings.yMarker = nan(eleNum, 3);
    settings.head = nan(eleNum, 3);     % not gonna be used

    for i=1:eleNum
        if strcmp(space, 'native')
            if isfield(reco,'scrf')
                if ~isempty(reco.scrf.markers) && ~isempty(reco.scrf.markers(i).y)
    	            settings.yMarker(i,:) = reco.scrf.markers(i).y;
                end
                if ~isempty(reco.scrf.markers) && ~isempty(reco.scrf.markers(i).head)
    	            settings.head(i,:) = reco.scrf.markers(i).head;
                end
            else
                if ~isempty(reco.native.markers) && ~isempty(reco.native.markers(i).y)
    	            settings.yMarker(i,:) = reco.native.markers(i).y;
                end
                if ~isempty(reco.native.markers) && ~isempty(reco.native.markers(i).head)
    	            settings.head(i,:) = reco.native.markers(i).head;
                end                
            end
        else
            if ~isempty(reco.mni.markers) && ~isempty(reco.mni.markers(i).y)
    	        settings.yMarker(i,:) = reco.mni.markers(i).y;
            end
            if ~isempty(reco.mni.markers) && ~isempty(reco.mni.markers(i).head)
    	        settings.head(i,:) = reco.mni.markers(i).head;
            end
        end

        if i == 1
            stim_json.marker_y_coords_rh = settings.yMarker(i,:);
            stim_json.marker_head_coords_rh = settings.head(i,:);
        else
            stim_json.marker_y_coords_lh = settings.yMarker(i,:);
            stim_json.marker_head_coords_lh = settings.head(i,:);   
        end

    end

    if strcmp(space, 'native')
        if isfield(reco,'scrf')
            coords_mm = reco.scrf.coords_mm;
        else
            coords_mm = reco.native.coords_mm;
        end
    else
        coords_mm = reco.mni.coords_mm;
    end

    % ToDo: unilateral cases
    stim_json.contact_coords_rh_x = coords_mm{1,1}(:,1);
    stim_json.contact_coords_rh_y = coords_mm{1,1}(:,2);
    stim_json.contact_coords_rh_z = coords_mm{1,1}(:,3);

    stim_json.contact_coords_lh_x = coords_mm{1,2}(:,1);
    stim_json.contact_coords_lh_y = coords_mm{1,2}(:,2);
    stim_json.contact_coords_lh_z = coords_mm{1,2}(:,3);

    % save the input to a json
    jsonString = jsonencode(stim_json, 'PrettyPrint', true);
    fid = fopen(DictPath, 'w');
    if fid == -1
        error('ea_merge:FileWriteError', 'Could not open file %s for writing.', DictPath);
    end
    
    try
        fprintf(fid, '%s', jsonString);
        disp(['Successfully saved formatted JSON to: ', DictPath]);
    catch ME
        fclose(fid);
        rethrow(ME);
    end
    fclose(fid);

    % run OSS-DBS
    env.system(['python ', ea_path_helper(OSS4DBS_RFs_script), ' ', ea_path_helper(DictPath)])

    % create niftis from Lattice (.csv)
    if VTR
        Result_folders = dir_without_dots([stim_json.stimFolder,filesep,'Results_protocol_*']);
        ea_mkdir([stim_json.stimFolder,filesep,'RFs']); 
    else
        Result_folders = dir_without_dots([stim_json.stimFolder,filesep,'Results_VTA*']);
    end
    

    %% Postprocessing
    subj_RFs = {};
    counter_rh = 0;
    counter_lh = 0;
    for vtr_i = 1:size(Result_folders,1)
    
        % convert each field from .csv to the required .nii
        if VTR
            if strcmp(Result_folders(vtr_i).name(end-2:end),'_rh')
                hemi_sfx = '_right';
                counter_rh = counter_rh + 1; % we can do it this way, because Results_folders are ordered
                cnt_idx = counter_rh;
            else
                hemi_sfx = '_left';
                counter_lh = counter_lh + 1;
                cnt_idx = counter_lh;
            end


            file2save = [stim_json.stimFolder,filesep,'RFs',filesep,'1mA_c',num2str(cnt_idx),'_efield', hemi_sfx ,'.nii'];  % can be adjusted as you wish, this is just an output
            %file2save_MNI = [Result_folders(vtr_i).folder,filesep,'RFs',filesep,'1mA_c',num2str(cnt_idx),'_efield_wb_MNI.nii'];  % can be adjusted as you wish, this is just an output
            field_in_csv = [Result_folders(vtr_i).folder,filesep,Result_folders(vtr_i).name,filesep,'E_field_Lattice.csv'];

            %field_in_csv = [Result_folders(vtr_i).folder,filesep,Result_folders(vtr_i).name,filesep,'oss_potentials_Lattice.csv'];
            %file2save = [Result_folders(vtr_i).folder,filesep,'1V_c',num2str(cnt_idx),'_phi', hemi_sfx ,'.nii'];  % can be adjusted as you wish, this is just an output
            subj_RFs{vtr_i,1} = file2save;
            if strcmp(space,'native') && warp2MNI 
                stimFolderMNI = fullfile(subjs{1,subj_inx},'stimulations','MNI152NLin2009bAsym','BSS');
                if ~isfolder(stimFolderMNI)
                    ea_mkdir(stimFolderMNI)
                end
                if ~isfolder([stimFolderMNI,filesep,'RFs'])
                    ea_mkdir([stimFolderMNI,filesep,'RFs']);
                end
                %file2save_MNI = [stimFolderMNI,filesep,'1V_c',num2str(cnt_idx),'_phi', hemi_sfx ,'.nii'];  
                file2save_MNI = [stimFolderMNI,filesep,'RFs',filesep,'1mA_c',num2str(cnt_idx),'_efield', hemi_sfx ,'.nii'];  % can be adjusted as you wish, this is just an output
            end
           
        else
            file2save = [Result_folders(vtr_i).folder,filesep,Result_folders(vtr_i).name(13:end),'_efield.nii'];
            file2save_MNI = [Result_folders(vtr_i).folder,filesep,Result_folders(vtr_i).name(13:end),'_efield_MNI.nii'];  % can be adjusted as you wish, this is just an output
            field_in_csv = [Result_folders(vtr_i).folder,filesep,Result_folders(vtr_i).name,filesep,'E_field_Lattice.csv'];
        end
    
        %get_sEEG_field_from_csv(field_in_csv, file2save, VTR, 0, 0, false)
        get_sEEG_field_from_csv(field_in_csv, file2save, 0, 0, 0, true)
    
        if warp2MNI
            % ToDo: Test!
            if ~isfile(file2save_MNI)
                get_sEEG_field_in_MNI_from_csv(field_in_csv, file2save_MNI, VTR, anchorImage, transform, true)
                get_sEEG_field_in_MNI_from_csv(field_in_csv, file2save_MNI, 0, anchorImage, transform, true)
            end
            subj_RFs{vtr_i,1} = file2save_MNI;
            % if cleanup
            %     % maybe we actually wanna preserve native RFs
            %     ea_delete(file2save);  
            % end
        end
        
        % flip left RFs if required
        if flip_lh2rh && strcmp(hemi_sfx,'_left')
            if warp2MNI
                file2save_MNI_fl = [file2save_MNI(1:end-4),'_fl.nii.gz'];
                ea_flip_lr_nonlinear([file2save_MNI,'.gz'],file2save_MNI_fl)
                if cleanup
                    ea_delete(file2save_MNI);
                end
                subj_RFs{vtr_i,1} = file2save_MNI_fl;
            else
                if strcmp(space,'native')
                    ea_error("LH2RH flip is only allowed in MNI space!")
                else
                    file2save_fl = [file2save(1:end-4),'_fl.nii.gz'];
                    ea_flip_lr_nonlinear([file2save,'.gz'],file2save_fl)
                    if cleanup
                        ea_delete(file2save);
                    end
                    subj_RFs{vtr_i,1} = file2save_fl;
                end
            end
        end
    end

    if reslice2ROI && (warp2MNI || strcmp(space,'mni'))
        images2ROI = [images2ROI;subj_RFs];
        % will be resliced outside
        if cleanup
            ea_delete([stim_json.stimFolder,filesep,'Results_protocol_*'])
        end

        continue
    elseif reslice2ROI
        % reslice to the common voxel space within subject
        resliceNii2ROI(subj_RFs,ROI_RF_threshold,vox_size)     

        if cleanup
            ea_delete([stim_json.stimFolder,filesep,'Results_protocol_*'])
        end
    end

end

if reslice2ROI && (warp2MNI || strcmp(space,'mni'))
    % reslice to the common voxel space across subjects
    resliceNii2ROI(images2ROI,ROI_RF_threshold,vox_size)    
end

