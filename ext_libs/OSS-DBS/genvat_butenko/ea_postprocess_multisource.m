function stimparams = ea_postprocess_multisource(options,settings,side,source_efields,source_vtas,outputPaths)
% Merge multisource VATs.
% By Butenko and Li, konstantinmgtu@gmail.com

arguments
    options             % Lead-DBS options for electrode reconstruction and stimulation
    settings            % parameters for OSS-DBS simulation
    side                {mustBeNumeric} % hemisphere index (0 - rh, 1 - lh)
    source_efields      % cell array, full paths to the e-field computed for source_use_index 
    source_vtas         % cell array, full paths to the VATs computed for source_use_index 
    outputPaths         % various paths to conform with lead-dbs BIDS structure 
end

switch side
    case 1
        sideLabel = 'R';
    case 2
        sideLabel = 'L';
end

% use SPM imcalc to merge fields
ea_merge_multisource_fields(outputPaths.outputBasePath,source_efields,side,settings.Activation_threshold_VTA(side))

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
            source_efields{side,i} = fullfile([outputPaths.templateOutputBasePath, 'efield_model-ossdbs_hemi-', sideLabel,'_S',num2str(i), '.nii']);
            source_vtas{side,i} = fullfile([outputPaths.templateOutputBasePath, 'binary_model-ossdbs_hemi-', sideLabel,'_S',num2str(i), '.nii']);
        end
    end

    ea_merge_multisource_fields(outputPaths.templateOutputBasePath,source_efields,side,settings.Activation_threshold_VTA(side))
    % clean-up to avoid any misimport downstream
    for i = 1:size(source_efields,2)
        if ~isempty(source_efields{side,i})
            ea_delete(fullfile([outputPaths.templateOutputBasePath, 'efield_model-ossdbs_hemi-', sideLabel,'_S',num2str(i), '.nii']));
            ea_delete(fullfile([outputPaths.templateOutputBasePath, 'binary_model-ossdbs_hemi-', sideLabel,'_S',num2str(i), '.nii']));
        end
    end
end

if options.native && ~options.orignative
    % Visualize MNI space VTA computed in native
    vatToViz = [outputPaths.templateOutputBasePath, 'binary_model-ossdbs_hemi-', sideLabel, '.nii'];
else
    vatToViz = [outputPaths.outputBasePath, 'binary_model-ossdbs_hemi-', sideLabel, '.nii'];
end

% Calc vat fv and volume
vat = ea_load_nii(vatToViz);
vatfv = ea_niiVAT2fvVAT(vat,1,3);
vatvolume = sum(vat.img(:))*vat.voxsize(1)*vat.voxsize(2)*vat.voxsize(3);
save(strrep(vatToViz, '.nii', '.mat'), 'vatfv', 'vatvolume');
stimparams(side).VAT.VAT = vatfv;
stimparams(side).volume = vatvolume;