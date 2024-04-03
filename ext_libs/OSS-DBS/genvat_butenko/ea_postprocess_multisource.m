function stimparams = ea_postprocess_multisource(options,settings,side,source_efields,source_vtas,templateOutputBasePath,outputBasePath)

switch side
    case 1
        sideLabel = 'R';
    case 2
        sideLabel = 'L';
end

ea_merge_multisource_fields(outputBasePath,source_efields,side,settings.Activation_threshold_VTA(side),sideLabel)

% clean-up to avoid any misimport downstream
for i = 1:size(source_efields,2)
    if ~isempty(source_efields{side,i})
        ea_delete(source_efields{side,i});
        ea_delete(source_vtas{side,i});
    end
end

% always transform to MNI space
if options.native
    % also merge in MNI space

    for i = 1:size(source_efields,2)
        if ~isempty(source_efields{side,i})
            source_efields{side,i} = fullfile([templateOutputBasePath, 'efield_model-ossdbs_hemi-', sideLabel,'_S',num2str(i), '.nii']);
            source_vtas{side,i} = fullfile([templateOutputBasePath, 'binary_model-ossdbs_hemi-', sideLabel,'_S',num2str(i), '.nii']);
        end
    end

    ea_merge_multisource_fields(templateOutputBasePath,source_efields,side,settings.Activation_threshold_VTA(side),sideLabel)
    % clean-up to avoid any misimport downstream
    for i = 1:size(source_efields,2)
        if ~isempty(source_efields{side,i})
            ea_delete(fullfile([templateOutputBasePath, 'efield_model-ossdbs_hemi-', sideLabel,'_S',num2str(i), '.nii']));
            ea_delete(fullfile([templateOutputBasePath, 'binary_model-ossdbs_hemi-', sideLabel,'_S',num2str(i), '.nii']));
        end
    end
end

if options.native && ~options.orignative
    % Visualize MNI space VTA computed in native
    vatToViz = [templateOutputBasePath, 'binary_model-ossdbs_hemi-', sideLabel, '.nii'];
else
    vatToViz = [outputBasePath, 'binary_model-ossdbs_hemi-', sideLabel, '.nii'];
end

% Calc vat fv and volume
vat = ea_load_nii(vatToViz);
vatfv = ea_niiVAT2fvVAT(vat,1,3);
vatvolume = sum(vat.img(:))*vat.voxsize(1)*vat.voxsize(2)*vat.voxsize(3);
save(strrep(vatToViz, '.nii', '.mat'), 'vatfv', 'vatvolume');
stimparams(side).VAT.VAT = vatfv;
stimparams(side).volume = vatvolume;