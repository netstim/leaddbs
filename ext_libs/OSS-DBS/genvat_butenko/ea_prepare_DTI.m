function DTI_data_name = ea_prepare_DTI(options,outputPaths)
% Warp and recompute DTI for native space modeling, otherwise just copy for
% templates/.
% By Butenko and Li, konstantinmgtu@gmail.com

arguments
    options     % Lead-DBS options for electrode reconstruction and stimulation
    outputPaths % various paths to conform with lead-dbs BIDS structure 
end

tensorName = options.prefs.machine.vatsettings.butenko_tensorFileName;
scalingMethod = options.prefs.machine.vatsettings.butenko_tensorScalingMethod;
scaledTensorName = strrep(tensorName, '.nii', ['_', scalingMethod, '.nii']);

ea_mkdir([options.subj.coregDir, filesep, 'dwi']);
nativeTensor = [options.subj.coregDir, filesep, 'dwi', filesep, outputPaths.subDescPrefix, tensorName];
nativeTensorScaled = [options.subj.coregDir, filesep, 'dwi', filesep, outputPaths.subDescPrefix, scaledTensorName];
templateTensor = [ea_space, tensorName];
templateTensorScaled = [ea_space, scaledTensorName];
tensorData = [outputPaths.outputDir, filesep, scaledTensorName]; % Final tensor data input for OSS-DBS

if options.prefs.machine.vatsettings.butenko_useTensorData
    if isfile(tensorData)
        % Scaled tensor data found in stimulation folder
        DTI_data_name = scaledTensorName;

    elseif ~options.native && isfile(templateTensorScaled)
        % MNI mode, scaled tensor data found in MNI space folder
        copyfile(templateTensorScaled, outputPaths.outputDir);
        DTI_data_name = scaledTensorName;

    elseif options.native && isfile(nativeTensorScaled)
        % native mode, scaled tensor data found in patient folder
        copyfile(nativeTensorScaled, tensorData);
        DTI_data_name = scaledTensorName;

    else
        if ~options.native
            % MNI mode, tensor data found
            if isfile(templateTensor)
                tensorDir = ea_space;
                tensorPrefix = '';
            end
        else
            % native mode, tensor data not found, warp template tensor data
            if ~isfile(nativeTensor) && isfile(templateTensor)
                % Warp tensor data only when ANTs was used for normalization
                json = loadjson(options.subj.norm.log.method);
                if contains(json.method, {'ANTs','EasyReg'})
                    fprintf('Warping tensor data into patient space...\n\n')
                    ea_ants_apply_transforms(options,...
                        [ea_space, tensorName],... % From
                        nativeTensor,... % To
                        1, ... % Useinverse is 1
                        '', ... % Reference, auto-detected
                        '', ... % Transformation, auto-detected
                        0, ... % NearestNeighbor interpolation
                        3, ... % Dimension
                        'tensor');
                else
                    warning('off', 'backtrace');
                    warning('Warping tensor data is only supported when ANTs was used for normalization! Skipping...');
                    warning('on', 'backtrace');
                end
            end

            if isfile(nativeTensor) % Scale tensor data
                tensorDir = fileparts(nativeTensor);
                tensorPrefix = outputPaths.subDescPrefix;
            end
        end

        % Scale tensor data
        if exist('tensorDir', 'var')
            fprintf('Scaling tensor data...\n\n')

            system(['python ', ea_getearoot, 'ext_libs/OSS-DBS/MRI_DTI_processing/Tensor_scaling.py ', tensorDir,filesep, tensorPrefix, tensorName, ' ', scalingMethod]);
            
            if ~isfile([tensorDir, filesep, tensorPrefix, scaledTensorName])
                disp('Parallel tensor scaling failed, trying a single thread...')
                system(['python ', ea_getearoot, 'ext_libs/OSS-DBS/MRI_DTI_processing/Tensor_scaling_one_thread.py ', tensorDir,filesep, tensorPrefix, tensorName, ' ', scalingMethod]);
            end

            % Copy scaled tensor data to stimulation directory, update setting
            copyfile([tensorDir, filesep, tensorPrefix, scaledTensorName], tensorData);
            DTI_data_name = scaledTensorName;
        end
    end

    fprintf('Scaled tensor data added: %s\n\n', DTI_data_name)
    % get the full path
    DTI_data_name = [outputPaths.outputDir, filesep, DTI_data_name];

else
    DTI_data_name = 'no dti';
end
