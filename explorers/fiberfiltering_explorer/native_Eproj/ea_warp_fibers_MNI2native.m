function ftr = ea_warp_fibers_MNI2native(pt_folder, MNI_connectome, transform, anchor_img)
    % Warp normative connectome to the native space
    % By Butenko, konstantinmgtu@gmail.com
    
    arguments
        pt_folder            % full path to the patient in Lead-DBS dataset
        MNI_connectome       % full path to the (merged) connectome in MNI space
        transform            % full path to the ANTs transform file from native to MNI
        anchor_img           % full path tio the Anchore image in native space
    end

    % BIDS
    %anchor_img = options.subj.preopAnat.(options.subj.AnchorModality).coreg;
    %transform = [options.subj.subjDir, filesep, 'forwardTransform'];

    t1 = [ea_space,'t1.nii'];
    %anchor_img = [pt_folder, '/anat_t1.nii'];
    %transform = [pt_folder, '/glanatComposite.nii.gz'];

    %fprintf('Loading connectome: %s ...\n', MNI_connectome);
    conn = load(MNI_connectome);
    % Convert connectome fibers from MNI space to anchor space
    fprintf('Convert connectome into native space...\n\n');
    fibersMNIVox = ea_mm2vox(conn.fibers(:,1:3), t1)';
    conn.fibers(:,1:3)  = ea_map_coords(fibersMNIVox, ...
        t1, ...
        transform, ...
        anchor_img, 'ANTS')';

    ftr.fibers = conn.fibers;
    ftr.idx = conn.idx;
    ftr.ea_fibformat = conn.ea_fibformat;
    ftr.fourindex = conn.fourindex;
    % store metadata
    ftr.orig_connectome = 'placeholder';

    [connectomePath,connectomeFileName,connectomeExtension] = fileparts(MNI_connectome);

    % check the folder name if stored in data or merged_pathways.mat
    if strcmp('data',connectomeFileName) || strcmp('merged_pathways',connectomeFileName)
        [~,connectomeName,~] = fileparts(connectomePath);
    else
        connectomeName = connectomeFileName;
    end
     
    % create miscellaneous if in patient folder
        mkdir([pt_folder, filesep,'connectomes',filesep,'dMRI', filesep, connectomeName])
    save(strcat(pt_folder, filesep,'connectomes',filesep,'dMRI', filesep, connectomeName, filesep, connectomeFileName, connectomeExtension),'-struct','ftr');
end
