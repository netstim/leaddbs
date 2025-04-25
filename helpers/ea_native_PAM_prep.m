function ea_native_PAM_prep(pt_folder,tracts_folder)
% Warp patient's tracts to MNI. This is just needed for a simplicity of
% processing. The computations will be done on the actual native space tracts!
% By Butenko, konstantinmgtu@gmail.com

arguments
    pt_folder                   % full path to the patient folder within Lead-DBS
    tracts_folder               % full path to the folder containing patient tracts in the ftr format
end

% auto-definitions
anchor_img = ea_regexpdir([pt_folder,'/coregistration/anat'], '.*ses-preop_space-anchorNative.*',0,'f',0);
% yes, use inverse warp!
transform = ea_regexpdir([pt_folder, '/normalization/transformations'],'.*from-MNI152NLin2009bAsym_to-anchorNative_desc-ants.nii.gz',0,'f',0);
t1 = [ea_space,'t1.nii'];

[~,connectome_name,~] = fileparts(tracts_folder);
native_tracts = dir_without_dots(tracts_folder);

% create the connectome folder
if size(native_tracts,1) > 1
    connectome_MNI = [ea_getconnectomebase,'dMRI_MultiTract',filesep,connectome_name];
else
    connectome_MNI = [ea_getconnectomebase,'dMRI',filesep,connectome_name];
    if ~strcmp(native_tracts(1).name,'data.mat')
        ea_warndlg("Rename 'the one file connectome' as data.mat")
        return
    end
end
   
% always re-create
if isfolder(connectome_MNI)
    ea_delete(connectome_MNI)
end
mkdir(connectome_MNI)

for tr_i = 1:size(native_tracts,1)

    native_tract = [native_tracts(tr_i).folder,filesep,native_tracts(tr_i).name];

    %fprintf('Loading connectome: %s ...\n', MNI_connectome);
    conn = load(native_tract);
    % Convert connectome fibers from MNI space to anchor space
    fprintf('Convert connectome into MNI space...\n\n');
    fibersNativeVox = ea_mm2vox(conn.fibers(:,1:3), anchor_img{1})';
    conn.fibers(:,1:3)  = ea_map_coords(fibersNativeVox, ...
        anchor_img{1}, ...
        transform{1}, ...
        t1, 'ANTS')';

    ftr.fibers = conn.fibers;
    ftr.idx = conn.idx;
    ftr.ea_fibformat = conn.ea_fibformat;
    ftr.fourindex = conn.fourindex;
    % store metadata
    ftr.orig_connectome = connectome_name;
    ftr.patient_specific = 1;

    
    if size(native_tracts,1) == 1
        save(strcat(connectome_MNI, filesep, 'data.mat'),'-struct','ftr');
    else
        save(strcat(connectome_MNI, filesep, native_tracts(tr_i).name),'-struct','ftr');
    end
end