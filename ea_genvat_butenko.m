function varargout = ea_genvat_butenko(varargin)
% Wrapper for OSS-DBS for VTA calculation

time = datetime('now', 'TimeZone', 'local');
timezone = time.TimeZone;
setenv('TZ', timezone);

% set to 1 if you only want to prep files for cluster comp.
prepFiles_cluster = 0; % for now hardcoded

if nargin==2
    S=varargin{1};
    options=varargin{2};
elseif nargin==3
    S=varargin{1};
    options=varargin{2};
    hFigure=varargin{3};
elseif nargin==1 && ischar(varargin{1}) % return name of method.
    varargout{1} = 'OSS-DBS (Butenko 2020)';
    varargout{2} = true; % Support directed lead
    return
end

% Check OSS-DBS installation, set env
env = ea_conda_env('OSS-DBSv2');
ea_checkOSSDBSInstallv2(env);

binPath = getenv('PATH'); % Current PATH
if isunix
    pythonPath = [env.path, filesep, 'bin'];
    setenv('PATH', [pythonPath, ':', binPath]);
else
    pythonPath = [env.path,';',env.path,filesep,'Scripts'];
    setenv('PATH', [pythonPath, ';', binPath]);
end

% Double check if lead is supported by OSS-DBS.
if ~ismember(options.elmodel, ea_ossdbs_elmodel)
    ea_error([options.elmodel, 'is not supported by OSS-DBS yet!'], simpleStack = 1);
end

% new parameters
settings.butenko_segmAlg = options.prefs.machine.vatsettings.butenko_segmAlg;
settings.butenko_intersectStatus = options.prefs.machine.vatsettings.butenko_intersectStatus;
settings.removeElectrode = options.prefs.machine.vatsettings.butenko_removeElectrode;
settings.neuronModel = options.prefs.machine.vatsettings.butenko_axonModel;
settings.signalType = options.prefs.machine.vatsettings.butenko_signalType;
%settings.pulseWidth = options.prefs.machine.vatsettings.butenko_pulseWidth;
settings.biphasic = options.prefs.machine.vatsettings.butenko_biphasic;
settings.butenko_tensorData = options.prefs.machine.vatsettings.butenko_tensorData;
settings.AdaptiveRef = options.prefs.machine.vatsettings.butenko_AdaptiveRef;
settings.encapsulationType = options.prefs.machine.vatsettings.butenko_encapsulation;
settings.adaptive_threshold = options.prefs.machine.vatsettings.butenko_adaptive_ethresh;
settings.outOfCore = 0;

% Set output path
subDescPrefix = ['sub-', options.subj.subjId, '_desc-'];
subSimPrefix = ['sub-', options.subj.subjId, '_sim-'];
outputDir = [options.subj.stimDir, filesep, ea_nt(options.native), S.label];
outputBasePath = [outputDir, filesep, subSimPrefix];
ea_mkdir(outputDir);
if options.native
    templateOutputDir = [options.subj.stimDir, filesep, ea_nt(0), S.label];
    ea_mkdir(templateOutputDir);
    templateOutputBasePath = [templateOutputDir, filesep, subSimPrefix];
end

%% Set MRI_data_name
% Segment MRI
segmaskName = 'segmask.nii';
anchorImage = options.subj.preopAnat.(options.subj.AnchorModality).coreg;
switch settings.butenko_segmAlg
    case 'SPM'
        if options.native
            [anchorImageDir, anchorImageName] = fileparts(anchorImage);
            anchorImageDir = [anchorImageDir, filesep];
            anchorImageName = [anchorImageName, '.nii'];

            mod = replace(options.subj.AnchorModality, textBoundary('start') + alphanumericsPattern + "_", "");
            c1File = setBIDSEntity(anchorImage, 'mod', mod, 'label', 'GM', 'suffix', 'mask');
            c2File = setBIDSEntity(anchorImage, 'mod', mod, 'label', 'WM', 'suffix', 'mask');
            c3File = setBIDSEntity(anchorImage, 'mod', mod, 'label', 'CSF', 'suffix', 'mask');
            if ~isfile(c1File) || ~isfile(c2File) || ~isfile(c3File)
                ea_newseg(anchorImage, 0, 1);
                movefile([anchorImageDir, 'c1', anchorImageName], c1File);
                movefile([anchorImageDir, 'c2', anchorImageName], c2File);
                movefile([anchorImageDir, 'c3', anchorImageName], c3File);
            end

            segMaskPath = setBIDSEntity(anchorImage, 'label', 'C123', 'mod', options.subj.AnchorModality, 'suffix', 'mask');
        else
            c1File = [ea_space, 'c1mask.nii'];
            c2File = [ea_space, 'c2mask.nii'];
            c3File = [ea_space, 'c3mask.nii'];
            if ~isfile(c1File) || ~isfile(c2File) || ~isfile(c3File)
                ea_newseg(fullfile(ea_space, [options.primarytemplate, '.nii']), 0, 1);
                movefile([ea_space, 'c1', options.primarytemplate, '.nii'], c1File);
                movefile([ea_space, 'c2', options.primarytemplate, '.nii'], c2File);
                movefile([ea_space, 'c3', options.primarytemplate, '.nii'], c3File);
            end

            segMaskPath = [ea_space, segmaskName];
        end

        if ~isfile(segMaskPath)
            % Binarize segmentations
            c1 = ea_load_nii(c1File);
            c2 = ea_load_nii(c2File);
            c3 = ea_load_nii(c3File);
            c1.img = c1.img>0.5;
            c2.img = c2.img>0.5;
            c3.img = c3.img>0.5;

            % Fuse segmentations by voting in the order  CSF -> WM -> GM
            c2.img(c3.img) = 0;
            c1.img(c2.img | c3.img) = 0;
            c1.img = c1.img + c2.img*2 + c3.img*3;
            c1.fname = segMaskPath;
            c1.dt = [2 0]; % unit8 according to spm_type
            c1.descrip = 'Tissue 1 + 2 + 3';
            c1.pinfo(1:2) = [1,0]; % uint8 is enough for output values, no need for scaling
            ea_write_nii(c1);
        end
        copyfile(segMaskPath, [outputDir, filesep, segmaskName]);
    case 'Atlas Based'

        % always overwrite in this case
        if isfile([outputDir, filesep, segmaskName])
            ea_delete([outputDir, filesep, segmaskName])
        end

        if options.native
            segMaskPath = [options.subj.atlasDir,filesep,options.atlasset,filesep,'segmask_atlas.nii'];
            atlas_gm_mask_path = [options.subj.atlasDir,filesep,options.atlasset,filesep,'gm_mask.nii.gz'];
            ea_convert_atlas2segmask(atlas_gm_mask_path, segMaskPath, 0.5)
            copyfile(segMaskPath, [outputDir, filesep, segmaskName]);
        else
            % save directly to stim folder
            segMaskPath = [outputDir,filesep,'segmask.nii'];
            atlas_gm_mask_path = [ea_space,filesep,'atlases',filesep,options.atlasset,filesep,'gm_mask.nii.gz'];
            ea_convert_atlas2segmask(atlas_gm_mask_path, segMaskPath, 0.5)
        end
    case'SynthSeg'
        ea_warndlg("SynthSeg segmentations are not supported yet")
        return
        %SynthSeg_segmask_image = ea_smart_BIDS_function_to_find_SynthSeg;
        %ea_convert_synthSeg2segmask(SynthSeg_segmask_image, segmask_output);
end

%% Set patient folder
settings.Patient_folder = options.subj.subjDir;

%% Set native/MNI flag
settings.Estimate_In_Template = ~options.native;

%% Set MRI path
settings.MRI_data_name = [outputDir, filesep, segmaskName];

%% Check tensor data
tensorName = options.prefs.machine.vatsettings.butenko_tensorFileName;
scalingMethod = options.prefs.machine.vatsettings.butenko_tensorScalingMethod;
scaledTensorName = strrep(tensorName, '.nii', ['_', scalingMethod, '.nii']);

ea_mkdir([options.subj.coregDir, filesep, 'dwi']);
nativeTensor = [options.subj.coregDir, filesep, 'dwi', filesep, subDescPrefix, tensorName];
nativeTensorScaled = [options.subj.coregDir, filesep, 'dwi', filesep, subDescPrefix, scaledTensorName];
templateTensor = [ea_space, tensorName];
templateTensorScaled = [ea_space, scaledTensorName];
tensorData = [outputDir, filesep, scaledTensorName]; % Final tensor data input for OSS-DBS

% initialize
settings.DTI_data_name = 'no dti';

if options.prefs.machine.vatsettings.butenko_useTensorData
    if isfile(tensorData)
        % Scaled tensor data found in stimulation folder
        settings.DTI_data_name = scaledTensorName;

    elseif ~options.native && isfile(templateTensorScaled)
        % MNI mode, scaled tensor data found in MNI space folder
        copyfile(templateTensorScaled, outputDir);
        settings.DTI_data_name = scaledTensorName;

    elseif options.native && isfile(nativeTensorScaled)
        % native mode, scaled tensor data found in patient folder
        copyfile(nativeTensorScaled, tensorData);
        settings.DTI_data_name = scaledTensorName;

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
                tensorPrefix = subDescPrefix;
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
            settings.DTI_data_name = scaledTensorName;
        end
    end
end

if ~isempty(settings.DTI_data_name) && ~strcmp('no dti', settings.DTI_data_name)
    fprintf('Scaled tensor data added: %s\n\n', settings.DTI_data_name)
end

if ~strcmp(settings.DTI_data_name, 'no dti')
    % get the full path
    settings.DTI_data_name = [outputDir, filesep, settings.DTI_data_name];
end

%% Index of the tissue in the segmented MRI data
settings.GM_index = 1;
settings.WM_index = 2;
settings.CSF_index = 3;

settings.default_material = 'GM'; % GM, WM or CSF

%% Electrodes information
settings.Electrode_type = options.elmodel;

% Reload reco since we need to decide whether to use native or MNI coordinates.
coords_mm = ea_load_reconstruction(options);
settings.contactLocation = coords_mm;
eleNum = length(coords_mm); % Number of electrodes
conNum = options.elspec.numel; % Number of contacts per electrode

% Save both native and MNI space y and head markers for OSS-DBS
settings.yMarkerNative = nan(eleNum, 3);
settings.yMarkerMNI = nan(eleNum, 3);
settings.headNative = nan(eleNum, 3);
settings.headMNI = nan(eleNum, 3);
[markersNative, markersMNI] = ea_get_markers(options);
for i=1:eleNum
    if ~isempty(markersNative) && ~isempty(markersNative(i).y)
    	settings.yMarkerNative(i,:) = markersNative(i).y;
    end
    if ~isempty(markersMNI) && ~isempty(markersMNI(i).y)
    	settings.yMarkerMNI(i,:) = markersMNI(i).y;
    end
    if ~isempty(markersNative) && ~isempty(markersNative(i).head)
    	settings.headNative(i,:) = markersNative(i).head;
    end
    if ~isempty(markersMNI) && ~isempty(markersMNI(i).head)
    	settings.headMNI(i,:) = markersMNI(i).head;
    end
end

% Head
settings.Implantation_coordinate = nan(eleNum, 3);
for i=1:eleNum
    if ~isempty(coords_mm{i})
        settings.Implantation_coordinate(i,:) = coords_mm{i}(1,:);
    end
end

% Tail
settings.Second_coordinate = nan(eleNum, 3);
for i=1:eleNum
    if conNum == 1 % Exception for electrode with only one contact
        if options.native
            settings.Second_coordinate(i,:) = markersNative(i).tail;
        else
            settings.Second_coordinate(i,:) = markersMNI(i).tail;
        end
    elseif ~isempty(coords_mm{i})
        if contains(options.elmodel, 'DIXI D08')
            settings.Second_coordinate(i,:) = coords_mm{i}(4,:);
        else
            settings.Second_coordinate(i,:) = coords_mm{i}(end,:);
        end
    end
end

%% Stimulation Information
% Set stimSetMode flag
settings.stimSetMode = options.stimSetMode;
if settings.stimSetMode
    ea_warndlg("Not yet supported in V2")
    return
end

% if right electrode only, insert null Stim Protocol for the left
% this work around is not needed for the left only, handled by Lead-DBS
if eleNum == 1
    S = ea_add_StimVector_to_S(S, zeros(1, conNum),1);
    eleNum = 2;
elseif isempty(coords_mm)
    S = ea_add_StimVector_to_S(S, zeros(1, conNum),0);
end

% Initialize current control flag
settings.current_control = nan(eleNum, 1);

% Initalize Phi vector
settings.Phi_vector = nan(eleNum, conNum);

% Initialize grounding
settings.Case_grounding = zeros(eleNum, 1);

% Get the stimulation parameters from S in case stimSetMode is 0, otherwise
% they will be loaded directly from the Current_protocols_[0|1].csv files
% also get the center of the grid
settings.stim_center = nan(2, 3);

nActiveSources = zeros(2,1); % count sources for each side
activeSources = nan(2,4);
multiSourceMode = [0;0];
runStatusMultiSource = zeros(4,2);  % check status for each source
for side = 1:2
    for source_index = 1:4
        if S.amplitude{side}(source_index) ~= 0 && ~isnan(S.amplitude{side}(source_index))
            nActiveSources(side) = nActiveSources(side) + 1;
            activeSources(side,source_index) = source_index;
        end
    end

    if nActiveSources(side) > 1
        multiSourceMode(side) = 1;
        if settings.adaptive_threshold
            ea_warndlg('Adaptive Thresholding is not supported for MultiSource. Switching to machine.vatsettings.butenko_ethresh')
            settings.adaptive_threshold = 0;
        end
    end
end

if any(nActiveSources > 1)
    source_efields = cell(2,4);
    source_vtas = cell(2,4);
    if options.prefs.machine.vatsettings.butenko_calcPAM
        ea_warndlg('MultiSource Mode is not supported for PAM!')
        % Restore working directory and environment variables
        runStatus = [0 0];
        varargout{1} = runStatus;
        varargout{2} = struct(); % empty for stimparameters
    
        % Restore working directory and environment variables
        libpath = getenv('LD_LIBRARY_PATH');
        setenv('LD_LIBRARY_PATH', ''); % Clear LD_LIBRARY_PATH to resolve conflicts
        setenv('LD_LIBRARY_PATH', libpath);
        setenv('PATH', binPath);
        return
    else
        ea_warndlg('MultiSource Mode is used! Stimulation volumes will be computed separately and merged using max(|E|)')
    end
end

stimparams = struct();

for source_index = 1:4

    settings.Activation_threshold_VTA = [];

    if ~settings.stimSetMode
        [settings.Phi_vector, settings.current_control, settings.Case_grounding] = ea_get_OneSourceStimVector(S, eleNum, conNum,activeSources(:,source_index));
        settings.pulseWidth = [double(S.(['Rs', num2str(source_index)]).pulseWidth);double(S.(['Ls', num2str(source_index)]).pulseWidth)];
        for side = 1:2

            % estimate center of VAT grid
            if ~isnan(settings.current_control(side))

                % Threshold for Astrom VTA (V/m)
                if settings.adaptive_threshold
                    settings.Activation_threshold_VTA = [settings.Activation_threshold_VTA;ea_get_adaptiveEthreshold(settings.pulseWidth(side))];
                else
                    settings.Activation_threshold_VTA = [settings.Activation_threshold_VTA;options.prefs.machine.vatsettings.butenko_ethresh];
                end

                stimamp = sum(abs(settings.Phi_vector(side,:)),"all",'omitnan');
                settings.stim_center(side,:) = sum(settings.contactLocation{side}.*abs(settings.Phi_vector(side,:)')./stimamp,1,'omitnan');
            
                 % estimate extent of the stimulation along the lead
                phi_temp = settings.Phi_vector(side,:);
                phi_temp(isnan(phi_temp)) = 0;
                first_active = find(phi_temp,1,'first');
                last_active = find(phi_temp,1,'last');
                length_active_span = norm(settings.contactLocation{side}(last_active,:)-settings.contactLocation{side}(first_active,:));
                if strcmp('Boston Scientific Vercise', options.elmodel)
                    if length_active_span > 24.0
                        ea_warndlg("Large span of active contacts is detected. Consider extending VAT grid, see PointModel.Lattice.Shape in lead_settings.py")
                    end
                elseif length_active_span > 18.0
                    ea_warndlg("Large span of active contacts is detected. Consider extending VAT grid, see PointModel.Lattice.Shape in lead_settings.py")
                end
            else
                % add nonsense value as a placeholder
                settings.Activation_threshold_VTA = [settings.Activation_threshold_VTA;-1.0];
            end

        end
    else
        settings.stim_center = [NaN;NaN];
    end

    
    % Set stimulation protocol
    if settings.stimSetMode
        stimProtocol = ea_regexpdir(outputDir, '^Current_protocols_\d\.csv$', 0);
    else
        stimProtocol = S;
    end
    
    %% check what we simulate
    settings.calcAxonActivation = options.prefs.machine.vatsettings.butenko_calcPAM;
    settings.exportVAT = options.prefs.machine.vatsettings.butenko_calcVAT;
    
    
    % Axon activation setting
    if settings.calcAxonActivation
        %settings.AxonModel = options.prefs.machine.vatsettings.butenko_AxonModel;

        settings.pulseWidth = [S.Rs1.pulseWidth;S.Ls1.pulseWidth];
        settings.connectome = options.prefs.machine.vatsettings.butenko_connectome;
        settings.axonLength = options.prefs.machine.vatsettings.butenko_axonLength;
        settings.fiberDiameter = options.prefs.machine.vatsettings.butenko_fiberDiameter;
    
    
        preopAnchor = options.subj.preopAnat.(options.subj.AnchorModality).coreg;
        if ~startsWith(settings.connectome, 'Multi-Tract: ') % Normal connectome
            fprintf('Loading connectome: %s ...\n', settings.connectome);
            conn = load([ea_getconnectomebase, 'dMRI', filesep, settings.connectome, filesep, 'data.mat']);
            if options.native
                originalFib = conn;
                % Convert connectome fibers from MNI space to anchor space
                fprintf('Convert connectome into native space...\n\n');
                fibersMNIVox = ea_mm2vox(conn.fibers(:,1:3), [ea_space, options.primarytemplate, '.nii'])';
                conn.fibers(:,1:3)  = ea_map_coords(fibersMNIVox, ...
                    [ea_space, options.primarytemplate, '.nii'], ...
                    [options.subj.subjDir, filesep, 'forwardTransform'], ...
                    preopAnchor)';
            end
    
            % Filter fibers based on the spherical ROI
            if options.native
        	    fiberFiltered = ea_filterfiber_stim(conn, coords_mm, stimProtocol, 'kuncel', 2, preopAnchor);
            else
                fiberFiltered = ea_filterfiber_stim(conn, coords_mm, stimProtocol, 'kuncel', 2, [ea_space, options.primarytemplate, '.nii']);
            end
    
            % Filter fibers based on the minimal length
            fiberFiltered = ea_filterfiber_len(fiberFiltered, settings.axonLength);
    
            % Move original fiber id to the 5th column, the 4th column will be 1:N
            fibersFound = zeros(size(fiberFiltered));
            for i=1:length(fiberFiltered)
                if ~isempty(fiberFiltered{i}.fibers)
                    fibers = zeros(size(fiberFiltered{i}.fibers,1),5);
                    fibers(:,[1,2,3,5]) = fiberFiltered{i}.fibers;
                    fibers(:,4) = repelem(1:length(fiberFiltered{i}.idx), fiberFiltered{i}.idx)';
                    fiberFiltered{i}.fibers = fibers;
                    fibersFound(i) = 1;
                end
            end
    
            settings.connectomePath = [outputDir, filesep, settings.connectome];
            ea_mkdir(settings.connectomePath);
            for i=1:length(fiberFiltered)
                % store the original number of fibers
                % to compute percent activation
                fiberFiltered{i}.origNum = size(conn.idx,1);
                buffer = fiberFiltered{i};
                save([settings.connectomePath, filesep, 'data', num2str(i), '.mat'], '-struct', 'buffer', '-v7.3');
            end
        else % Multi-Tract connectome
            % Extract connectome name
            connName = strrep(settings.connectome, 'Multi-Tract: ', '');
    
            % Create output folder
            settings.connectomePath = [outputDir, filesep, connName];
            ea_mkdir(settings.connectomePath);
    
            % Get paths of tracts
            connFolder = [ea_getconnectomebase, 'dMRI_MultiTract', filesep, connName];
            tracts = ea_regexpdir(connFolder, '\.mat$', 0);
    
            settings.connectomeTractNames = cell(size(tracts));
            data1 = struct;
            data2 = struct;
            fibersFound = zeros(numel(tracts),2);
            for t=1:numel(tracts)
                tract = tracts{t};
                [~, tractName] = fileparts(tract);
                settings.connectomeTractNames{t} = tractName;
                fprintf('Loading connectome: %s, Tract: %s ...\n', connName, tractName);
                conn = load(tract);
                if options.native
                    originalFib = conn;
                    % Convert connectome fibers from MNI space to anchor space
                    fprintf('Convert connectome into native space...\n\n');
                    fibersMNIVox = ea_mm2vox(conn.fibers(:,1:3), [ea_space, options.primarytemplate, '.nii'])';
                    conn.fibers(:,1:3)  = ea_map_coords(fibersMNIVox, ...
                        [ea_space, options.primarytemplate, '.nii'], ...
                        [options.subj.subjDir, filesep, 'forwardTransform'], ...
                        preopAnchor)';
                end
    
                % Filter fibers based on the spherical ROI
                if options.native
        	        fiberFiltered = ea_filterfiber_stim(conn, coords_mm, stimProtocol, 'kuncel', 2, preopAnchor);
                else
                    fiberFiltered = ea_filterfiber_stim(conn, coords_mm, stimProtocol, 'kuncel', 2, [ea_space, options.primarytemplate, '.nii']);
                end
    
                % Filter fibers based on the minimal length
                fiberFiltered = ea_filterfiber_len(fiberFiltered, settings.axonLength(t));
    
                % Move original fiber id to the 5th column, the 4th column will be 1:N
                for i=1:length(fiberFiltered)
                    if ~isempty(fiberFiltered{i}.fibers)
                        fibers = zeros(size(fiberFiltered{i}.fibers,1),5);
                        fibers(:,[1,2,3,5]) = fiberFiltered{i}.fibers;
                        fibers(:,4) = repelem(1:length(fiberFiltered{i}.idx), fiberFiltered{i}.idx)';
                        fiberFiltered{i}.fibers = fibers;
    
                        % store the original number of fibers
                        % to compute percent activation
                        fiberFiltered{i}.origNum = size(conn.idx,1);
    
                        fibersFound(t,i) = 1;
                    end
                end
                data1.tractName = eval('fiberFiltered{1}');
                data2.tractName = eval('fiberFiltered{2}');
            end
    
            % Save filtered fibers
            save([settings.connectomePath, filesep, 'data1.mat'], '-struct', 'data1', '-v7.3');
            save([settings.connectomePath, filesep, 'data2.mat'], '-struct', 'data2', '-v7.3');
        end
    else
        % for now, hardcode 60 us pw for all VATs
        settings.pulseWidth = [60.0;60.0];
    end
    
    %% Save settings for OSS-DBS
    if any(~isnan(settings.current_control))
        parameterFile = fullfile(outputDir, 'oss-dbs_parameters.mat');
        save(parameterFile, 'settings', '-v7.3');
        ea_savestimulation(S, options);
        if options.native
            poptions = options;
            poptions.native = 0;
            ea_savestimulation(S, poptions);
        end
        
        % Delete previous results from stimSetMode
        ea_delete([outputDir, filesep, 'Result_StimProt_*']);
        if options.native
            ea_delete([templateOutputDir, filesep, 'Result_StimProt_*']);
        end
        
        % full clean-up for V2
        ea_delete([outputDir, filesep, 'Results_*']);
        ea_delete([outputDir, filesep, 'oss-dbs_parameters.json']);
        ea_delete([outputDir, filesep, 'Allocated_axons.h5'])
        
        %% Run OSS-DBS
        libpath = getenv('LD_LIBRARY_PATH');
        setenv('LD_LIBRARY_PATH', ''); % Clear LD_LIBRARY_PATH to resolve conflicts
        
        % Delete flag files before running
        ea_delete([outputDir, filesep, 'success_rh.txt']);
        ea_delete([outputDir, filesep, 'fail_rh.txt']);
        ea_delete([outputDir, filesep, 'skip_rh.txt']);
        ea_delete([outputDir, filesep, 'success_lh.txt']);
        ea_delete([outputDir, filesep, 'fail_lh.txt']);
        ea_delete([outputDir, filesep, 'skip_lh.txt']);
    end
    
    if prepFiles_cluster == 1
        % Restore working directory and environment variables
        %runStatusMultiSource(source_index,:) = [0 0];
        varargout{1} = [0 0];
        varargout{2} = struct(); % empty for stimparameters
    
        % Restore working directory and environment variables
        setenv('LD_LIBRARY_PATH', libpath);
        setenv('PATH', binPath);
        return
    end
    
    % Iterate sides, index side: 0 - rh , 1 - lh
    for side = 0:1
    
        if ~multiSourceMode(side+1)
            % not relevant in this case, terminate after one iterationrce_vtas = {};
            source_use_index = 5;  
        else
            source_use_index = source_index; 
        end

        if ~multiSourceMode(side+1) && all(isnan(settings.current_control))
            % skip without message
            runStatusMultiSource(source_index,side+1) = 1;
            continue
        end

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
    
        if ~settings.stimSetMode && isnan(settings.current_control(side+1))
            warning('off', 'backtrace');
            warning('No stimulation exists for %s side! Skipping...\n', sideStr);
            warning('on', 'backtrace');
            fclose(fopen([outputDir, filesep, 'skip_', sideCode, '.txt'], 'w'));
            runStatusMultiSource(source_index,side+1) = 1;
            continue;
        end
    
        if settings.stimSetMode && ~isfile(strcat(outputDir, filesep, 'Current_protocols_',string(side),'.csv'))
            warning('off', 'backtrace');
            warning('No stimulation set for %s side! Skipping...\n', sideStr);
            warning('on', 'backtrace');
            fclose(fopen([outputDir, filesep, 'skip_', sideCode, '.txt'], 'w'));
            runStatusMultiSource(source_index,side+1) = 1;
            continue;
        end
    
        if settings.calcAxonActivation && ~any(fibersFound(:,side+1))
            warning('off', 'backtrace');
            warning('No fibers found for %s side! Skipping...\n', sideStr);
            warning('on', 'backtrace');
            fclose(fopen([outputDir, filesep, 'skip_', sideCode, '.txt'], 'w'));
            runStatusMultiSource(source_index,side+1) = 1;
            continue;
        end
    
        fprintf('\nRunning OSS-DBS for %s side stimulation...\n\n', sideStr);
    
        if settings.calcAxonActivation
            ea_delete([outputDir, filesep, 'Allocated_axons.h5']);
            system(['python ', ea_getearoot, 'ext_libs/OSS-DBS/Axon_Processing/axon_allocation.py ', outputDir,' ', num2str(side), ' ', parameterFile]);
            % call axon_allocation script
        end
    
        % use OSS-DBS v2 environment
        system(['leaddbs2ossdbs --hemi_side ', num2str(side), ' ', parameterFile, ...
            ' --output_path ', outputDir]);
        parameterFile_json = [parameterFile(1:end-3), 'json'];
        system(['ossdbs ' , parameterFile_json]);
    
        if settings.calcAxonActivation
            % copy NEURON folder to the stimulation folder 
            leaddbs_neuron = [ea_getearoot, 'ext_libs/OSS-DBS/Axon_Processing/Axon_files'];
            neuron_folder = fullfile(outputDir,'Axon_files');
            copyfile(leaddbs_neuron, neuron_folder)
    
            % call the NEURON module
            folder2save = [outputDir,filesep,'Results_', sideCode];
            timeDomainSolution = [outputDir,filesep,'Results_', sideCode, filesep, 'oss_time_result_PAM.h5'];
            pathwayParameterFile = [outputDir,filesep, 'Allocated_axons_parameters.json'];

            % check if the time domain results is available
            if ~isfile(timeDomainSolution)
                ea_warndlg('OSS-DBS failed to prepare a time domain solution. If RAM consumption exceeded the hardware limit, set settings.outOfCore to 1')
                return
            end
    
            system(['python ', ea_getearoot, 'ext_libs/OSS-DBS/Axon_Processing/PAM_caller.py ', neuron_folder, ' ', folder2save,' ', timeDomainSolution, ' ', pathwayParameterFile]);
        end
    
        % clean-up for outOfCore
        if settings.outOfCore == 1
            ea_delete([outputDir, filesep, 'Results_',sideCode,filesep,'oss_freq_domain_tmp_PAM.hdf5']);
        end
    
        % Check if OSS-DBS calculation is finished
        while ~isfile([outputDir, filesep, 'success_', sideCode, '.txt']) ...
                && ~isfile([outputDir, filesep, 'fail_', sideCode, '.txt'])
            continue;
        end
    
        % Get resultfig handle
        if exist('hFigure', 'var')
            resultfig = getappdata(hFigure,'resultfig');
        end
    
        if isfile([outputDir, filesep, 'success_', sideCode, '.txt'])
            runStatusMultiSource(source_index,side+1) = 1;
            fprintf('\nOSS-DBS calculation succeeded!\n\n')
    
            if settings.exportVAT
    
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
                    stimparams(side+1).VAT.VAT = vatfv;
                    stimparams(side+1).volume = vatvolume;
                else
                    source_efields{side+1,source_use_index} = fullfile([outputBasePath, 'efield_model-ossdbs_hemi-', sideLabel,'_S',num2str(source_use_index), '.nii']);
                    source_vtas{side+1,source_use_index} = fullfile([outputBasePath, 'binary_model-ossdbs_hemi-', sideLabel,'_S',num2str(source_use_index), '.nii']);
                end
            end
    
    
            axonState = ea_regexpdir([outputDir, filesep, 'Results_', sideCode], 'Axon_state.*\.mat', 0);
            if ~isempty(axonState)
                for f=1:length(axonState)
    
                    % Determine tract name
                    if startsWith(settings.connectome, 'Multi-Tract: ')
                        tractName = regexp(axonState{f}, '(?<=Axon_state_).+(?=\.mat$)', 'match', 'once');
                    end
    
                    % If stimSetMode, extract the index from tractName (but axonState is still checked on the indexed file)
                    if settings.stimSetMode
                        if startsWith(settings.connectome, 'Multi-Tract: ')
                            stimProt_index = regexp(tractName, '(?<=_)\d+$', 'match', 'once');
                            tractName = regexp(tractName, '.+(?=_\d+$)', 'match', 'once');
                        else
                            stimProt_index = regexp(axonState{f}, '(?<=Axon_state_)\d+(?=\.mat$)', 'match', 'once');
                        end
                    end
    
                    % Get fiber id and state from OSS-DBS result
                    ftr = load(axonState{f});
                    [fibId, ind] = unique(ftr.fibers(:,4));
    
                    fibState = ftr.fibers(ind,5);
    
                    % Restore full length fiber (as in original filtered fiber)
                    if startsWith(settings.connectome, 'Multi-Tract: ')
                        ftr = load([settings.connectomePath, filesep, 'data', num2str(side+1), '.mat'], tractName);
                        ftr = ftr.(tractName);
                        ftr.connectome_name = connName;
                    else
                        ftr = load([settings.connectomePath, filesep, 'data', num2str(side+1), '.mat']);
                        ftr.connectome_name = settings.connectome;
                    end
                    ftr.fibers = ftr.fibers(ismember(ftr.fibers(:,4), fibId), :);
                    originalFibID = ftr.fibers(:,5);
    
                    % Extract original conn fiber id and idx, needed in case
                    % calculation is done in native space
                    [connFibID, idx] = unique(ftr.fibers(:,5));
    
                    % Set fiber state
                    for fib=1:length(fibId)
                        ftr.fibers(ftr.fibers(:,4)==fibId(fib),5) = fibState(fib);
                    end
    
                    % Extract state of original conn fiber, needed in case
                    % calculation is  done in native space
                    connFibState = ftr.fibers(idx, 5);
    
                    % Reset original fiber id as in the connectome
                    ftr.fibers(:,4) = originalFibID;
    
                    if strcmp(settings.butenko_intersectStatus,'activated')
                        ftr.fibers(ftr.fibers(:,5) == -1 | ftr.fibers(:,5) == -3, 5) = 1;
                    elseif strcmp(settings.butenko_intersectStatus,'activated_at_active_contacts')
                        ftr.fibers = OSS_DBS_Damaged2Activated(settings,ftr.fibers,ftr.idx,side+1);
                    end
    
                    % Save result for visualization
    
                    % If stimSets, save to a corresponding folder
                    if settings.stimSetMode
                        resultProtocol = [outputDir, filesep, 'Result_StimProt_', sideStr, '_', stimProt_index];
                        ea_mkdir(resultProtocol);
                        if startsWith(settings.connectome, 'Multi-Tract: ')
                            fiberActivation = [resultProtocol, filesep, subSimPrefix, 'fiberActivation_model-ossdbs_hemi-', sideLabel, '_tract-', tractName,'_prot-', stimProt_index, '.mat'];
                        else
                            fiberActivation = [resultProtocol, filesep, subSimPrefix, 'fiberActivation_model-ossdbs_hemi-', sideLabel,'_prot-', stimProt_index, '.mat'];
                        end
                    else
                        if startsWith(settings.connectome, 'Multi-Tract: ')
                            fiberActivation = [outputBasePath, 'fiberActivation_model-ossdbs_hemi-', sideLabel, '_tract-', tractName, '.mat'];
                        else
                            fiberActivation = [outputBasePath, 'fiberActivation_model-ossdbs_hemi-', sideLabel, '.mat'];
                        end
                    end
    
                    save(fiberActivation, '-struct', 'ftr');
    
                    if options.native % Generate fiber activation file in MNI space
                        fprintf('Restore connectome in MNI space: %s ...\n\n', settings.connectome);
    
                        if startsWith(settings.connectome, 'Multi-Tract: ')
                            % load the particular pathway in MNI
                            [atlas_folder,~] = fileparts(tract);
                            originalFib = load([atlas_folder,filesep,tractName]);
                        end
                        conn = originalFib;
    
                        fprintf('Convert fiber activation result into MNI space...\n\n');
                        conn.fibers = conn.fibers(ismember(conn.fibers(:,4), connFibID), :);
                        % Set fiber state
                        conn.fibers = [conn.fibers, zeros(size(conn.fibers,1),1)];
                        for fib=1:length(connFibID)
                            conn.fibers(conn.fibers(:,4)==connFibID(fib),5) = connFibState(fib);
                        end
    
                        % Recreate fiber idx
                        [~, ~, idx] = unique(conn.fibers(:,4), 'stable');
                        conn.idx = accumarray(idx,1);
    
                        % Reset original fiber id as in the connectome
                        ftr.fibers(:,4) = originalFibID;
    
                        % Save MNI space fiber activation result
    
                        % If stimSets, save to a corresponding folder
                        if settings.stimSetMode
                            resultProtocol = [templateOutputDir, filesep, 'Result_StimProt_', sideStr, '_', stimProt_index];
                            ea_mkdir(resultProtocol);
                            if startsWith(settings.connectome, 'Multi-Tract: ')
                                fiberActivationMNI = [resultProtocol, filesep, subSimPrefix, 'fiberActivation_model-ossdbs_hemi-', sideLabel, '_tract-', tractName,'_prot-', stimProt_index, '.mat'];
                            else
                                fiberActivationMNI = [resultProtocol, filesep, subSimPrefix, 'fiberActivation_model-ossdbs_hemi-', sideLabel, '_prot-', stimProt_index, '.mat'];
                            end
                        else
                            if startsWith(settings.connectome, 'Multi-Tract: ')
                                fiberActivationMNI = [templateOutputBasePath, 'fiberActivation_model-ossdbs_hemi-', sideLabel, '_tract-', tractName, '.mat'];
                            else
                                fiberActivationMNI = [templateOutputBasePath, 'fiberActivation_model-ossdbs_hemi-', sideLabel, '.mat'];
                            end
                        end
    
                        if startsWith(settings.connectome, 'Multi-Tract: ')
                            conn.connectome_name = connName;
                        else
                            conn.connectome_name = settings.connectome;
                        end
    
                        save(fiberActivationMNI, '-struct', 'conn');
    
                        if ~options.orignative % Visualize MNI space result
                            fiberActivation = fiberActivationMNI;
                        end
                    end
    
                    % Visualize fiber activation, but not for stimSetMode
                    if ~settings.stimSetMode
                        if exist('resultfig', 'var')
                            set(0, 'CurrentFigure', resultfig);
                            ea_fiberactivation_viz(fiberActivation, resultfig);
                        end
                    end
                end
            end
        elseif isfile([outputDir, filesep, 'fail_', sideCode, '.txt'])
            fprintf('\n')
            warning('off', 'backtrace');
            warning('OSS-DBS calculation failed for %s side!\n', sideStr);
            warning('on', 'backtrace');
            %runStatus(side+1) = 0;
        end
    
        % Clean up
        ea_delete([outputDir, filesep, 'Brain_substitute.brep']);
        %ea_delete(ea_regexpdir(outputDir, '^(?!Current_protocols_).*\.csv$', 0));
        % ea_delete([outputDir, filesep, '*.py']);
    
        % Delete this folder in MATLAB since shutil.rmtree may raise
        % I/O error
        % ea_delete([outputDir, filesep, 'Axons_in_time']);
    end

    % check only the first source for PAM
    if settings.calcAxonActivation
        runStatusMultiSource(2:end,:) = 1;
        break
    end

end

% here we merge and display multiSourceModes
for side = 1:2
    if multiSourceMode(side)

        switch side
            case 1
                sideLabel = 'R';
            case 2
                sideLabel = 'L';
        end

        if nActiveSources(side,:) > 0
            % don't call it if all zeros
            ea_merge_multisource_fields(outputBasePath,source_efields,side,settings.Activation_threshold_VTA(side),sideLabel)
        else
            continue
        end

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

    end
end

% unused sources are set to 1 above
runStatus = [all(runStatusMultiSource(:,1)==1), all(runStatusMultiSource(:,2)==1)];
varargout{1} = runStatus;

if ~settings.calcAxonActivation && exist('stimparams', 'var')
    varargout{2} = stimparams;
end

% Restore working directory and environment variables
setenv('LD_LIBRARY_PATH', libpath);
setenv('PATH', binPath);


%% Helper function to get markers in bothe native and MNI space
function [markersNative, markersMNI] = ea_get_markers(options)
options.native = 1;
try
    [~, ~, markersNative] = ea_load_reconstruction(options);
catch
    markersNative = [];
    fprintf('\n')
    warning('off', 'backtrace');
    warning('Failed to load native reconstruction!');
    warning('on', 'backtrace');
end

options.native = 0;
try
    [~, ~, markersMNI] = ea_load_reconstruction(options);
catch
    markersMNI = [];
    fprintf('\n')
    warning('off', 'backtrace');
    warning('Failed to load MNI reconstruction!');
    warning('on', 'backtrace');
end
