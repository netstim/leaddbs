function varargout = ea_genvat_butenko(varargin)
% Wrapper for OSS-DBS for VTA calculation

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


time = datetime('now', 'TimeZone', 'local');
timezone = time.TimeZone;
setenv('TZ', timezone);
% set to 1 if you only want to prep files for cluster comp.
prepFiles_cluster = 0; % for now hardcoded

settings = ea_prepare_ossdbs(options);

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

% Segment MRI image
settings = ea_segment_MRI(options, settings, outputDir);

% Prepare tensor data
settings.DTI_data_name = ea_prepare_DTI(options);

% get electrode reconstruction parameters in OSS-DBS format
[settings,eleNum] = ea_get_oss_reco(options, settings);

%% Stimulation Information
[settings, runStatusMultiSource, activeSources] = ea_check_stimSources(S,settings,eleNum);
nActiveSources = [nnz(~isnan(activeSources(1,:))), nnz(~isnan(activeSources(2,:)))];
stimparams = struct();

% Set stimulation protocol
if settings.stimSetMode
    stimProtocol = ea_regexpdir(outputDir, '^Current_protocols_\d\.csv$', 0);
else
    stimProtocol = S;
end

if any(nActiveSources > 1)
    % files to store results for each source
    source_efields = cell(2,4);
    source_vtas = cell(2,4);
    if options.prefs.machine.vatsettings.butenko_calcPAM
        ea_warndlg('MultiSource Mode is not supported for PAM!')
        % Restore working directory and environment variables
        [varargout{1}, varargout{2}] = ea_exit_genvat_butenko();
        return
    else
        ea_warndlg('MultiSource Mode is used! Stimulation volumes will be computed separately and merged using max(|E|)')
    end
end

% if single source, we will run only one iteration
for source_index = 1:4

    % get stim settings for particular source    
    settings = ea_get_stimProtocol(S,settings,activeSources);

    % check what we simulate
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
                eval(['data1.', tractName, ' = fiberFiltered{1};']);
                eval(['data2.', tractName, ' = fiberFiltered{2};']);
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
        [varargout{1}, varargout{2}] = ea_exit_genvat_butenko();
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
setenv('LD_LIBRARY_PATH', getenv('LD_LIBRARY_PATH'));
setenv('PATH', getenv('PATH'));


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


%% Helper function to fail exit
function [runStatus, stimparameters] = ea_exit_genvat_butenko()

    runStatus = [0 0];  % runStatus
    stimparameters = struct(); % empty for stimparameters

    % Restore working directory and environment variables
    setenv('LD_LIBRARY_PATH', ''); % Clear LD_LIBRARY_PATH to resolve conflicts
    setenv('LD_LIBRARY_PATH', getenv('LD_LIBRARY_PATH'));
    setenv('PATH', getenv('PATH'));


