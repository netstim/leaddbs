
% Inputs
sEEG_table = '/home/forel/Documents/data/SEEG/tst/reconstruction_sheet_midsize.tsv';
sEEG_stim_protocols = {'/home/forel/Documents/data/SEEG/tst/stim_protocols_EL1.csv'
    '/home/forel/Documents/data/SEEG/tst/stim_protocols_EL2.csv'};
anchor = '/home/forel/Documents/data/SEEG/tst/anat_t1.nii.gz';
multisegment = '/home/forel/Documents/data/SEEG/tst/anat_t1_synthSeg.nii.gz';  % output of the segmentation alg.
segmaskFile = '/home/forel/Documents/data/SEEG/tst/segmask.nii';
reslice2segmask = false; % if true, all VTRs are stored in the same voxel space defined by segmask
                        % this is handy when voxelwise operations are used
                        % but each VTR nii is > 100 MBs
VTR = false;             % if true, computes recording fields (aka volume of tissue recorded)
Stim_Mode = 'VC';       % CC - current-controlled, VC - voltage-controlled
warp2MNI = false;   % warp from the anchor space to MNI
transform = 'JohnDoe_from-MNI152NLin2009bAsym_to-anchorNative-ANTS.nii.gz';      % from native to MNI, you need the opposite warp, i.e. from-MNI152NLin2009bAsym_to-anchorNative
                                                                                 % IMPORTANT: make sure ants is mentioned in the name and use single quotes
% auto-definitions
[basepath,~,~] = fileparts(sEEG_table);
if VTR
    OSS_sEEG_script = [ea_getearoot, '/ext_libs/OSS-DBS/sEEG/run_OSS4SEEG.py'];
else
    OSS_sEEG_script = [ea_getearoot, '/ext_libs/OSS-DBS/sEEG/run_OSS4SEEG_Stim_no_shift.py'];
end

% segment the anchor modality
% TBD: multimodal segmentation with SynthSeg
if ~isfile(segmaskFile)
    if ~isfile(multisegment)
        ea_synthseg(anchor, multisegment)
    end
    ea_convert_synthSeg2segmask(multisegment, segmaskFile);
end

% Set OSS-DBS python path
env = ea_conda_env('OSS-DBSv2');
ea_checkOSSDBSInstallv2(env);

% call a python script to compute stim. volumes
if VTR
    env.system(['python ', ea_path_helper(OSS_sEEG_script), ' ', ea_path_helper(sEEG_table),  ' ',Stim_Mode])
else
    % each electrode has a separate sitm protocol sheet
    for elec = 1:size(sEEG_stim_protocols,1)
        env.system(['python ', ea_path_helper(OSS_sEEG_script), ' ', ea_path_helper(sEEG_table), ' ',ea_path_helper(sEEG_stim_protocols{elec,1}),  ' ',Stim_Mode, ' ', num2str(elec)])
    end
end

% create niftis from Lattice (.csv)
if VTR
    Result_folders = dir_without_dots([basepath,filesep,'Results_VTR*']);
else
    Result_folders = dir_without_dots([basepath,filesep,'Results_VTA*']);
end

mono_field = {};
for vtr_i = 1:size(Result_folders,1)

    if VTR
        file2save = [Result_folders(vtr_i).folder,filesep,Result_folders(vtr_i).name(13:end),'_potential.nii'];  % can be adjusted as you wish, this is just an output
        file2save_MNI = [Result_folders(vtr_i).folder,filesep,Result_folders(vtr_i).name(13:end),'_potential_MNI.nii'];  % can be adjusted as you wish, this is just an output
        field_in_csv = [Result_folders(vtr_i).folder,filesep,Result_folders(vtr_i).name,filesep,'oss_potentials_Lattice.csv'];
    else
        file2save = [Result_folders(vtr_i).folder,filesep,Result_folders(vtr_i).name(13:end),'_efield.nii'];
        file2save_MNI = [Result_folders(vtr_i).folder,filesep,Result_folders(vtr_i).name(13:end),'_efield_MNI.nii'];  % can be adjusted as you wish, this is just an output
        field_in_csv = [Result_folders(vtr_i).folder,filesep,Result_folders(vtr_i).name,filesep,'E_field_Lattice.csv'];
    end

    if ~isfile(file2save)
        get_sEEG_field_from_csv(field_in_csv, file2save, VTR, reslice2segmask, segmaskFile)
    end

    if warp2MNI
        % ToDo: Test!
        if ~isfile(file2save_MNI)
            get_sEEG_field_in_MNI_from_csv(field_in_csv, file2save_MNI, VTR, anchor, transform)
        end
    end
end
