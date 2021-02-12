function varargout = ea_genvat_butenko(varargin)
% Wrapper for OSS-DBS for VTA calculation

if nargin==2
    S=varargin{1};
    options=varargin{2};
elseif nargin==3
    S=varargin{1};
    options=varargin{2};
    lgfigure=varargin{3};
elseif nargin==1 && ischar(varargin{1}) % return name of method.
    varargout{1} = 'OSS-DBS (Butenko 2020)';
    return
end

% Check OSS-DBS installation, set env
if ~options.prefs.machine.vatsettings.oss_dbs.installed
    ea_checkOSSDBSInstall;
else
    binPath = getenv('PATH'); % Current PATH
    pythonPath = options.prefs.env.pythonPath;
    if isunix
        setenv('PATH', [pythonPath, ':', binPath]);
        pythonBinName = 'python3';
    else
        setenv('PATH', [pythonPath, ';', binPath]);
        pythonBinName = 'python';
    end
end

% Double check if lead is supported by OSS-DBS.
if ~ismember(options.elmodel, ea_ossdbs_elmodel)
    ea_error([options.elmodel, 'is not supported by OSS-DBS yet!'], 'Error', dbstack)
end

directory = [options.root, options.patientname, filesep];

if ~exist([directory,'stimulations',filesep,ea_nt(options.native),S.label],'dir')
    mkdir([directory,'stimulations',filesep,ea_nt(options.native),S.label]);
end

options = ea_assignpretra(options);

%% Set MRI_data_name
% Segment MRI
if options.native
    if ~isfile([directory, 'c1', options.prefs.prenii_unnormalized]) ...
            || ~isfile([directory, 'c2', options.prefs.prenii_unnormalized]) ...
            || ~isfile([directory, 'c3', options.prefs.prenii_unnormalized])
        ea_newseg(directory, options.prefs.prenii_unnormalized, 0, options, 1);
    end

    segMaskDir = directory;
    segFileSuffix = options.prefs.prenii_unnormalized;
else
    if ~isfile([ea_space, 'c1mask.nii']) ...
            || ~isfile([ea_space, 'c2mask.nii']) ...
            || ~isfile([ea_space, 'c3mask.nii'])
        ea_newseg(ea_space, 't1.nii', 0, options, 1);
        movefile([ea_space, 'c1t1.nii'], [ea_space, 'c1mask.nii']);
        movefile([ea_space, 'c2t1.nii'], [ea_space, 'c2mask.nii']);
        movefile([ea_space, 'c3t1.nii'], [ea_space, 'c3mask.nii']);
    end

    segMaskDir = ea_space;
    segFileSuffix = 'mask.nii';
end

if ~isfile([segMaskDir, 'segmask.nii'])
    % Binarize segmentations
    c1 = ea_load_nii([segMaskDir, 'c1', segFileSuffix]);
    c2 = ea_load_nii([segMaskDir, 'c2', segFileSuffix]);
    c3 = ea_load_nii([segMaskDir, 'c3', segFileSuffix]);
    c1.img = c1.img>0.5;
    c2.img = c2.img>0.5;
    c3.img = c3.img>0.5;

    % Fuse segmentations by voting in the order  CSF -> WM -> GM
    c2.img(c3.img) = 0;
    c1.img(c2.img | c3.img) = 0;
    c1.fname = [segMaskDir, 'segmask.nii'];
    c1.dt = [4 0];
    c1.img = int16(c1.img) + int16(c2.img)*2 + int16(c3.img)*3;
    ea_write_nii(c1);
end

%% Set patient folder
settings.Patient_folder = directory;

%% Set native/MNI flag
settings.Estimate_In_Template = options.prefs.machine.vatsettings.estimateInTemplate;

%% Set MRI path
% Put the MRI file in stimulation folder
copyfile([segMaskDir, 'segmask.nii'], [directory,'stimulations',filesep,ea_nt(options.native),S.label]);
settings.MRI_data_name = [directory,'stimulations',filesep,ea_nt(options.native),S.label,filesep,'segmask.nii'];

%% Scaled tensor data
settings.DTI_data_name = ''; % 'dti_tensor.nii';

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
    if ~isempty(coords_mm{i})
        settings.Second_coordinate(i,:) = coords_mm{i}(end,:);
    end
end

% Rotation around the lead axis in degrees
settings.Rotation_Z = 0.0;

%% Stimulation Information
source = nan(eleNum,1);
for i=1:eleNum
    if any(S.amplitude{i})
        source(i) = find(S.amplitude{i},1);
    end
end

% 0 - Voltage Control; 1 - Current Control
settings.current_control = nan(eleNum,1);
for i=1:eleNum
   	switch i
        case 1
            sideCode = 'R';
        case 2
            sideCode = 'L';
    end

    if ~isnan(source(i))
        settings.current_control(i) = uint8(S.([sideCode, 's', num2str(source(i))]).va==2);
    end
end

% Get stimulation amplitude
amp = nan(eleNum,1);
for i=1:eleNum
    if ~isnan(source(i))
        amp(i) = S.amplitude{i}(source(i));
    end
end

% Initalize Phi vector
settings.Phi_vector = nan(eleNum,options.elspec.numel);

% Initialize grounding
settings.Case_grounding = zeros(eleNum,1);

for i = 1:eleNum
    switch i
        case 1
            sideCode = 'R';
            cntlabel = {'k0','k1','k2','k3','k4','k5','k6','k7'};
        case 2
            sideCode = 'L';
            cntlabel = {'k8','k9','k10','k11','k12','k13','k14','k15'};
    end

    if ~isnan(source(i))
        stimSource = S.([sideCode, 's', num2str(source(i))]);
        for cnt = 1:options.elspec.numel
            if S.activecontacts{i}(cnt)
                switch stimSource.(cntlabel{cnt}).pol
                    case 1 % Negative, cathode
                        settings.Phi_vector(i, cnt) = -amp(i)*stimSource.(cntlabel{cnt}).perc/100;
                    case 2 % Postive, anode
                        settings.Phi_vector(i, cnt) = amp(i)*stimSource.(cntlabel{cnt}).perc/100;
                end
            end
        end
        if stimSource.case.perc == 100
            settings.Case_grounding(i) = 1;
        end
    end
end

% Threshold for Astrom VTA (V/mm)
settings.Activation_threshold_VTA = options.prefs.machine.vatsettings.butenko_ethresh;

% Set output path
outputPath = [directory, 'stimulations', filesep, ea_nt(options.native), S.label];
if options.native
    MNIoutputPath = [directory, 'stimulations', filesep, ea_nt(0), S.label];
    ea_mkdir(MNIoutputPath);
end

% Axon activation setting
settings.calcAxonActivation = options.prefs.machine.vatsettings.butenko_calcAxonActivation;
if settings.calcAxonActivation
    settings.connectome = options.prefs.machine.vatsettings.butenko_connectome;
    settings.axonLength = options.prefs.machine.vatsettings.butenko_axonLength;
    settings.fiberDiameter = options.prefs.machine.vatsettings.butenko_fiberDiameter;
    fprintf('Loading connectome: %s ...\n', settings.connectome);
    conn = load([ea_getconnectomebase, 'dMRI', filesep, settings.connectome, filesep, 'data.mat']);
    if options.native
        originalFib = conn.fibers;
        % Convert connectome fibers from MNI space to anchor space
        fprintf('Convert connectome into native space...\n\n');
        fibersMNIVox = ea_mm2vox(conn.fibers(:,1:3), [ea_space, 't1.nii'])';
        conn.fibers(:,1:3)  = ea_map_coords(fibersMNIVox, ...
            [ea_space, 't1.nii'], ...
            [directory, 'y_ea_normparams.nii'], ...
            [directory, options.prefs.prenii_unnormalized])';
    end

    % Filter fibers based on the spherical ROI
    fiberFiltered = ea_filterfiber_stim(conn, coords_mm, S, 'kuncel', 2);

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

    settings.connectomePath = [outputPath, filesep, settings.connectome];
    ea_mkdir(settings.connectomePath);
    for i=1:length(fiberFiltered)
        buffer = fiberFiltered{i};
        save([settings.connectomePath, filesep, 'data', num2str(i), '.mat'], '-struct', 'buffer', '-v7.3');
    end
end

% Interactive mode setting
settings.interactiveMode = options.prefs.machine.vatsettings.butenko_interactive;

%% Save settings for OSS-DBS
parameterFile = [outputPath, filesep, 'oss-dbs_parameters.mat'];
save(parameterFile, 'settings', '-v7.3');

%% Run OSS-DBS
libpath = getenv('LD_LIBRARY_PATH');
setenv('LD_LIBRARY_PATH', ''); % Clear LD_LIBRARY_PATH to resolve conflicts

% Delete flag files before running
ea_delete([outputPath, filesep, 'success_rh.txt']);
ea_delete([outputPath, filesep, 'fail_rh.txt']);
ea_delete([outputPath, filesep, 'skip_rh.txt']);
ea_delete([outputPath, filesep, 'success_lh.txt']);
ea_delete([outputPath, filesep, 'fail_lh.txt']);
ea_delete([outputPath, filesep, 'skip_lh.txt']);

if ispc || ismac
    dockerImage = 'ningfei/oss-dbs';
else % Linux
    dockerImage = 'custom_oss-dbs';
end

% Iterate sides, index side: 0 - rh , 1 - lh
runStatus = [0 0]; % Succeed or not
stimparams = struct();
for side=0:1
    % Stop and Remove running docker container on start
    if isempty(getenv('SINGULARITY_NAME')) % Only do it when using docker
        [~, containerID] = system(['docker ps -qf ancestor=', dockerImage]);
        if ~isempty(containerID)
            containerID = strsplit(strip(containerID));
            fprintf('\nStop running container...\n')
            cellfun(@(id) system(['docker stop ', id, newline]), containerID);
            % fprintf('\nClean up running container...\n')
            % cellfun(@(id) system(['docker rm ', id, newline]), containerID);
        end
    end

    switch side
        case 0
            sideCode = 'rh';
            sideStr = 'right';
        case 1
            sideCode = 'lh';
            sideStr = 'left';
    end

    if isnan(settings.current_control(side+1))
        warning('off', 'backtrace');
        warning('No stimulation exists for %s side! Skipping...\n', sideStr);
        warning('on', 'backtrace');
        fclose(fopen([outputPath, filesep, 'skip_', sideCode, '.txt'], 'w'));
        continue;
    end

    if settings.calcAxonActivation && ~fibersFound(side+1)
        warning('off', 'backtrace');
        warning('No fibers found for %s side! Skipping...\n', sideStr);
        warning('on', 'backtrace');
        fclose(fopen([outputPath, filesep, 'skip_', sideCode, '.txt'], 'w'));
        continue;
    end

    fprintf('\nRunning OSS-DBS for %s side stimulation...\n\n', sideStr);

    % Calculate axon allocation when option enabled
    if settings.calcAxonActivation
            fprintf('Calculating axon allocation for %s side stimulation...\n\n', sideStr);

            % Make sure to clean up, useful in manually interruption
            ea_delete([outputPath, filesep, 'Brain_substitute.brep']);
            ea_delete([outputPath, filesep,'Allocated_axons.h5']);
            ea_delete([outputPath, filesep,'*.csv']);
            ea_delete([outputPath, filesep,'*.py']);

            % Delete this folder in MATLAB since shutil.rmtree may raise
            % I/O error
            ea_delete([outputPath, filesep,'Axons_in_time']);

            if isempty(getenv('SINGULARITY_NAME')) % Docker
                system(['docker run ', ...
                        '--volume ', ea_getearoot, 'ext_libs/OSS-DBS:/opt/OSS-DBS ', ...
                        '--volume ', outputPath, ':/opt/Patient ', ...
                        '--rm ', dockerImage, ' ', ...
                        'python3 /opt/OSS-DBS/OSS_platform/Axon_allocation.py /opt/Patient ', num2str(side)]);
            else % Singularity
                system(['python3 ', ea_getearoot, 'ext_libs/OSS-DBS/OSS_platform/Axon_allocation.py ', outputPath, ' ', num2str(side)]);
            end
    end

    % Call OSS-DBS GUI to start calculation
    system([pythonBinName, ' ', ea_getearoot, 'ext_libs/OSS-DBS/OSS_platform/OSS-DBS_LeadDBS_integrator.py ', ...
            parameterFile, ' ', num2str(side)]);	% 0 is right side, 1 is the left side here

    % Check if OSS-DBS calculation is finished
    while ~isfile([outputPath, filesep, 'success_', sideCode, '.txt']) ...
            && ~isfile([outputPath, filesep, 'fail_', sideCode, '.txt'])
        continue;
    end

    % Get resultfig handle
    if exist('lgfigure', 'var')
        resultfig = getappdata(lgfigure,'resultfig');
    end

    if isfile([outputPath, filesep, 'success_', sideCode, '.txt'])
        runStatus(side+1) = 1;
        fprintf('\nOSS-DBS calculation succeeded!\n\n')
        % Copy VAT files
        if isfile([outputPath, filesep, 'Results_', sideCode, filesep, 'E_field_solution.nii'])
            copyfile([outputPath, filesep, 'Results_', sideCode, filesep, 'E_field_solution.nii'], ...
                     [outputPath, filesep, 'vat_efield_', sideStr, '.nii'])
            if options.native % Transform to MNI space
                ea_apply_normalization_tofile(options,...
                    [outputPath, filesep, 'vat_efield_', sideStr, '.nii'],... % from
                    [MNIoutputPath, filesep, 'vat_efield_', sideStr, '.nii'],... % to
                    directory,... % patient directory
                    0, ... % useinverse is 0
                    1, ... % linear interpolation
                    [ea_space, 't1.nii']);
                ea_autocrop([MNIoutputPath, filesep, 'vat_efield_', sideStr, '.nii']);
            end
        end

        if isfile([outputPath, filesep, 'Results_', sideCode, filesep, 'VTA_solution.nii'])
            copyfile([outputPath, filesep, 'Results_', sideCode, filesep, 'VTA_solution.nii'], ...
                     [outputPath, filesep, 'vat_', sideStr, '.nii'])

            vatToViz = [outputPath, filesep, 'vat_', sideStr, '.nii'];
            if options.native % Transform to MNI space
                ea_apply_normalization_tofile(options,...
                    [outputPath, filesep, 'vat_', sideStr, '.nii'],... % from
                    [MNIoutputPath, filesep, 'vat_', sideStr, '.nii'],... % to
                    directory,... % patient directory
                    0, ... % useinverse is 0
                    0, ... % nn interpolation
                    [ea_space, 't1.nii']);
                ea_autocrop([MNIoutputPath, filesep, 'vat_', sideStr, '.nii']);

                if ~options.orignative % Visualize MNI space VTA
                    vatToViz = [MNIoutputPath, filesep, 'vat_', sideStr, '.nii'];
                end
            end

            % Calc vat fv and volume
            vat = ea_load_nii(vatToViz);
            vatfv = ea_niiVAT2fvVAT(vat);
            vatvolume = sum(vat.img(:))*vat.voxsize(1)*vat.voxsize(2)*vat.voxsize(3);
            save([outputPath, filesep, 'vat_', sideStr, '.mat'], 'vatfv', 'vatvolume');
            stimparams(side+1).VAT.VAT = vatfv;
            stimparams(side+1).volume = vatvolume;
        end

        if isfile([outputPath, filesep, 'Results_', sideCode, filesep, 'Axon_state_data',num2str(side+1),'.mat'])
            % Get fiber id and state from OSS-DBS result
            ftr = load([outputPath, filesep, 'Results_', sideCode, filesep, 'Axon_state_data',num2str(side+1),'.mat']);
            [fibId, ind] = unique(ftr.fibers(:,4));
            fibState = ftr.fibers(ind,5);

            % Restore full length fiber (as in original filtered fiber)
            ftr = load([settings.connectomePath, filesep, 'data', num2str(side+1), '.mat']);
            ftr.fibers = ftr.fibers(ismember(ftr.fibers(:,4), fibId), :);
            originalFibID = ftr.fibers(:,5);

            % Extract original conn fiber id and idx, needed in case
            % calculation is done in native space
            [connFibID, idx] = unique(ftr.fibers(:,5));

            % Set fiber state
            for f=1:length(fibId)
                ftr.fibers(ftr.fibers(:,4)==fibId(f),5) = fibState(f);
            end

            % Extract state of original conn fiber, needed in case
            % calculation is  done in native space
            connFibState = ftr.fibers(idx, 5);

            % Reset original fiber id as in the connectome
            ftr.fibers(:,4) = originalFibID;

            % Save result for visualization
            save([outputPath, filesep, 'fiberActivation_', sideStr, '.mat'], '-struct', 'ftr');
            fiberActivation = [outputPath, filesep, 'fiberActivation_', sideStr, '.mat'];

            if options.native % Generate fiber activation file in MNI space
                fprintf('Restore connectome in MNI space: %s ...\n\n', settings.connectome);
                conn.fibers = originalFib;

                fprintf('Convert fiber activation result into MNI space...\n\n');
                conn.fibers = conn.fibers(ismember(conn.fibers(:,4), connFibID), :);
                % Set fiber state
                conn.fibers = [conn.fibers, zeros(size(conn.fibers,1),1)];
                for f=1:length(connFibID)
                    conn.fibers(conn.fibers(:,4)==connFibID(f),5) = connFibState(f);
                end

                % Recreate fiber idx
                [~, ~, idx] = unique(conn.fibers(:,4), 'stable');
                conn.idx = accumarray(idx,1);

                % Reset original fiber id as in the connectome
                ftr.fibers(:,4) = originalFibID;

                % Save MNI space fiber activation result
                save([MNIoutputPath, filesep, 'fiberActivation_', sideStr, '.mat'], '-struct', 'conn');

                if ~options.orignative % Visualize MNI space result
                    fiberActivation = [MNIoutputPath, filesep, 'fiberActivation_', sideStr, '.mat'];
                end
            end

            % Visualize fiber activation
            if exist('resultfig', 'var')
                set(0, 'CurrentFigure', resultfig);
                ea_fiberactivation_viz(fiberActivation, resultfig);
            end
        end
    elseif isfile([outputPath, filesep, 'fail_', sideCode, '.txt'])
        fprintf('\n')
        warning('off', 'backtrace');
        warning('OSS-DBS calculation failed for %s side!\n', sideStr);
        warning('on', 'backtrace');
    end

    % Clean up
    ea_delete([outputPath, filesep, 'Brain_substitute.brep']);
    ea_delete([outputPath, filesep,'Allocated_axons.h5']);
    ea_delete([outputPath, filesep,'*.csv']);
    ea_delete([outputPath, filesep,'*.py']);

    % Delete this folder in MATLAB since shutil.rmtree may raise
    % I/O error
    ea_delete([outputPath, filesep,'Axons_in_time']);
end

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
