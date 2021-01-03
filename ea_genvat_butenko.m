function varargout = ea_genvat_butenko(varargin)
% Wrapper for OSS-DBS for VTA calculation

if nargin==5
    acoords=varargin{1};
    S=varargin{2};
    side=varargin{3};
    options=varargin{4};
    stimname=varargin{5};
elseif nargin==6
    acoords=varargin{1};
    S=varargin{2};
    side=varargin{3};
    options=varargin{4};
    stimname=varargin{5};
    lgfigure=varargin{6};
elseif nargin==1 && ischar(varargin{1}) % return name of method.
    varargout{1} = 'OSS-DBS (Butenko 2020)';
    return
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
[~, ~, markers] = ea_load_reconstruction(options);
coords_mm = ea_resolvecoords(markers, options);
settings.contactLocation = coords_mm;
eleNum = length(coords_mm); % Number of electrodes

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

% Set grounding
settings.Case_grounding = nan(eleNum,1);
for i=1:eleNum
    if ~isnan(source(i))
        settings.Case_grounding(i) = 0;
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

% Axon activation setting
settings.calcAxonActivation = options.prefs.machine.vatsettings.butenko_calcAxonActivation;
if settings.calcAxonActivation
    settings.connectome = options.prefs.machine.vatsettings.butenko_connectome;
    settings.axonLength = options.prefs.machine.vatsettings.butenko_axonLength;
    settings.fiberDiameter = options.prefs.machine.vatsettings.butenko_fiberDiameter;
    conn = load([ea_getconnectomebase, 'dMRI', filesep, settings.connectome, filesep, 'data.mat']);
    if options.native
        % Convert connectome fibers from MNI space to anchor space
        fibersMNIVox = ea_mm2vox(conn.fibers(:,1:3), [ea_space, 't1.nii'])';
        conn.fibers(:,1:3)  = ea_map_coords(fibersMNIVox, ...
            [ea_space, 't1.nii'], ...
            [directory, 'y_ea_normparams.nii'], ...
            [directory, options.prefs.prenii_unnormalized])';
    end

    % Filter fibers based on the spherical ROI
    fiberFiltered = ea_filterfiber_stim(conn, coords_mm, S, 'kuncel');

    % Filter fibers based on the minimal length
    fiberFiltered = ea_filterfiber_len(fiberFiltered, settings.axonLength);

    settings.connectomePath = [outputPath, filesep, settings.connectome];
    ea_mkdir(settings.connectomePath);
    for i=1:length(fiberFiltered)
        buffer = fiberFiltered{i};
        save([settings.connectomePath, filesep, 'data', num2str(i), '.mat'], '-struct', 'buffer');
    end
end

%% Save settings for OSS-DBS
parameterFile = [outputPath, filesep, 'oss-dbs_parameters.mat'];
save(parameterFile, 'settings', '-v7.3');

%% Run OSS-DBS
currentPath = pwd;
libpath = getenv('LD_LIBRARY_PATH');
setenv('LD_LIBRARY_PATH', ''); % Clear LD_LIBRARY_PATH to resolve conflicts

% Delete flag files before running
ea_delete([outputPath, filesep, 'success_rh.txt']);
ea_delete([outputPath, filesep, 'fail_rh.txt']);
ea_delete([outputPath, filesep, 'success_lh.txt']);
ea_delete([outputPath, filesep, 'fail_lh.txt']);

% Iterate sides, index side: 0 - rh , 1 - lh
for side=0:1
    switch side
        case 0
            disp('Running OSS-DBS for right side stimulation...');
            sideCode = 'rh';
            sideStr = 'right';
        case 1
            disp('Running OSS-DBS for left side stimulation...');
            sideCode = 'lh';
            sideStr = 'left';
    end

    % Calculate axon allocation when option enabled
    if settings.calcAxonActivation
            switch side
                case 0
                    fprintf('Calculating axon allocation for right side stimulation...\n');
                case 1
                    fprintf('Calculating axon allocation for left side stimulation...\n');
            end

            system(['docker run ', ...
                    '--volume ', ea_getearoot, 'ext_libs/OSS-DBS:/opt/OSS-DBS ', ...
                    '--volume ', outputPath, ':/opt/Patient ', ...
                    '-it --rm sfbelaine/oss_dbs:python_latest ', ...
                    'python3 /opt/OSS-DBS/OSS_platform/Axon_allocation.py ', num2str(side)]);
    end

    % Call OSS-DBS GUI to start calculation
    cd([ea_getearoot, 'ext_libs/OSS-DBS/OSS_platform']);
    system(['cd "', ea_getearoot, 'ext_libs/OSS-DBS/OSS_platform";', ...
            'python3 ', ea_getearoot, 'ext_libs/OSS-DBS/OSS_platform/OSS-DBS_LeadDBS_integrator.py ', ...
            parameterFile,' ', num2str(side)]);	% 0 is right side, 1 is the left side here

    % Check if OSS-DBS calculation is finished
    while ~isfile([outputPath, filesep, 'success_', sideCode, '.txt']) ...
            && ~isfile([outputPath, filesep, 'fail_', sideCode, '.txt'])
        continue;
    end

    if isfile([outputPath, filesep, 'success_', sideCode, '.txt'])
        disp('OSS-DBS calculation succeeded!')
        % Copy VAT files
        if isfile([outputPath, filesep, 'Results_', sideCode, filesep, 'E_field_solution.nii'])
            copyfile([outputPath, filesep, 'Results_', sideCode, filesep, 'E_field_solution.nii'], ...
                     [outputPath, filesep, 'vat_efield_', sideStr, '.nii'])
        end

        if isfile([outputPath, filesep, 'Results_', sideCode, filesep, 'VAT_solution.nii'])
            copyfile([outputPath, filesep, 'Results_', sideCode, filesep, 'VAT_solution.nii'], ...
                     [outputPath, filesep, 'vat_', sideStr, '.nii'])
        end
    elseif isfile([outputPath, filesep, 'fail_', sideCode, '.txt'])
        warning('off', 'backtrace');
        warning('OSS-DBS calculation failed for %s side!', sideStr);
        warning('on', 'backtrace');
    end
end

% Restore working directory and env
cd(currentPath);
setenv('LD_LIBRARY_PATH', libpath);
