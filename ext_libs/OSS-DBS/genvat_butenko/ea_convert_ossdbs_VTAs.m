function [vatfv,vatvolume,source_efield_side_source,source_vta_side_source] = ea_convert_ossdbs_VTAs(options,settings,side,multiSourceMode, source_use_index,outputDir,outputBasePath,templateOutputBasePath)

anchorImage = options.subj.preopAnat.(options.subj.AnchorModality).coreg;
source_efield_side_source = [];
source_vta_side_source = [];
vatfv = [];
vatvolume = [];


switch side
    case 0
        sideLabel = 'R';
        sideCode = 'rh';
        sideStr = 'right';
    case 1
        sideLabel = 'L';
        sideCode = 'lh';
        sideStr = 'left';
end

if settings.removeElectrode
    % create nii for distorted grid
    if options.native
        ea_get_field_from_csv(anchorImage, [outputDir, filesep, 'Results_', sideCode, filesep,'E_field_Lattice.csv'], settings.Activation_threshold_VTA(side+1), sideLabel, outputBasePath, source_use_index)
    else
        ea_get_field_from_csv([ea_space, options.primarytemplate, '.nii'], [outputDir, filesep, 'Results_', sideCode, filesep,'E_field_Lattice.csv'], settings.Activation_threshold_VTA(side+1), sideLabel, outputBasePath, source_use_index)
    end
else
    % convert original OSS-DBS VTAs to BIDS in the corresponding space
    if ~multiSourceMode(side+1)
        copyfile(fullfile([outputDir, filesep, 'Results_', sideCode, filesep,'E_field_solution_Lattice.nii']), fullfile([outputBasePath, 'efield_model-ossdbs_hemi-', sideLabel, '.nii']));
        copyfile(fullfile([outputDir, filesep, 'Results_', sideCode, filesep,'VTA_solution_Lattice.nii']), fullfile([outputBasePath, 'binary_model-ossdbs_hemi-', sideLabel, '.nii']));
    else
        copyfile(fullfile([outputDir, filesep, 'Results_', sideCode, filesep,'E_field_solution_Lattice.nii']), fullfile([outputBasePath, 'efield_model-ossdbs_hemi-', sideLabel,'_S',num2str(source_use_index), '.nii']));
        copyfile(fullfile([outputDir, filesep, 'Results_', sideCode, filesep,'VTA_solution_Lattice.nii']), fullfile([outputBasePath, 'binary_model-ossdbs_hemi-', sideLabel,'_S',num2str(source_use_index), '.nii']));
    end
    %ea_autocrop([outputBasePath, 'binary_model-ossdbs_hemi-', sideLabel, '.nii'], margin=10);
    %ea_autocrop([outputBasePath, 'efield_model-ossdbs_hemi-', sideLabel, '.nii'], margin=10);
end

% always transform to MNI space
if options.native
    ea_get_MNI_field_from_csv(options, [outputDir, filesep, 'Results_', sideCode, filesep,'E_field_Lattice.csv'], settings.Activation_threshold_VTA(side+1), sideLabel, templateOutputBasePath, source_use_index)
end

if options.native && ~options.orignative &&  ~multiSourceMode(side+1)
    % Visualize MNI space VTA computed in native
    vatToViz = [templateOutputBasePath, 'binary_model-ossdbs_hemi-', sideLabel, '.nii'];
else
    vatToViz = [outputBasePath, 'binary_model-ossdbs_hemi-', sideLabel, '.nii'];
end

if ~multiSourceMode(side+1)
    % Calc vat fv and volume
    vat = ea_load_nii(vatToViz);
    vatfv = ea_niiVAT2fvVAT(vat,1,3);
    vatvolume = sum(vat.img(:))*vat.voxsize(1)*vat.voxsize(2)*vat.voxsize(3);
    save(strrep(vatToViz, '.nii', '.mat'), 'vatfv', 'vatvolume');
else
    source_efield_side_source = fullfile([outputBasePath, 'efield_model-ossdbs_hemi-', sideLabel,'_S',num2str(source_use_index), '.nii']);
    source_vta_side_source = fullfile([outputBasePath, 'binary_model-ossdbs_hemi-', sideLabel,'_S',num2str(source_use_index), '.nii']);
end