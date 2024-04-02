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

% some hardcoded parameters, can be added to GUI later
prepFiles_cluster = 0; % set to 1 if you only want to prep files for cluster comp.
true_VTA = 0; % set to 1 to compute classic VAT using axonal grids
prob_PAM = 0; % set to 1 to compute PAM over uncertain parameters (e.g. fiber diameter)

% import settings from Lead-DBS GUI
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
% Get resultfig handle
if exist('hFigure', 'var')
    resultfig = getappdata(hFigure,'resultfig');
end

% Segment MRI image
settings = ea_segment_MRI(options, settings, outputDir);

% Prepare tensor data
settings.DTI_data_name = ea_prepare_DTI(options,subDescPrefix,outputDir);

% get electrode reconstruction parameters in OSS-DBS format
[settings,eleNum,conNum] = ea_get_oss_reco(options, settings);

%% Stimulation Information
[settings, runStatusMultiSource, activeSources, multiSourceMode] = ea_check_stimSources(options,S,settings,eleNum,conNum);
nActiveSources = [nnz(~isnan(activeSources(1,:))), nnz(~isnan(activeSources(2,:)))];
stimparams = struct();

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
    settings = ea_get_stimProtocol(options,S,settings,activeSources,conNum,source_index);
    
    % Axon activation setting
    if settings.calcAxonActivation
        % will exit after the first source
        if true_VTA
            % custom case of PAM is original VTA
            % should be tested for unilateral!
            fibersFound = [[0],[0]];
            if ~isnan(activeSources(1,source_index))
                settings = ea_switch2VATgrid(options,settings,0);
                fibersFound(:,1) = 1; 
            end
            if ~isnan(activeSources(side+1,source_index))
                settings = ea_switch2VATgrid(options,settings,1);
                fibersFound(:,2) = 1;
            end
        else
            [fibersFound] = ea_prepare_fibers(options, S, settings);
        end
    end
    
    % Save settings for OSS-DBS
    parameterFile = ea_save_ossdbs_settings(options, S, settings, outputDir, templateOutputDir);
    
    if prepFiles_cluster == 1
        % now you can run OSS-DBS externally
        [varargout{1}, varargout{2}] = ea_exit_genvat_butenko();
        return
    end
    
    % Iterate sides, index side: 0 - rh , 1 - lh
    for side = 0:1
    
        if ~multiSourceMode(side+1)
            % not relevant in this case, terminate after one iteration;
            source_use_index = 5;  
        else
            source_use_index = source_index; 
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
    
        % skip non-active sources when using single source
        if ~multiSourceMode(side+1) && isnan(activeSources(side+1,source_index))
            runStatusMultiSource(source_index,side+1) = 1;
            continue
        end

        % skip if non-active source for this side
        if ~settings.stimSetMode && isnan(activeSources(side+1,source_index))
            if ~multiSourceMode(side+1)
                warning('off', 'backtrace');
                warning('No stimulation exists for %s side! Skipping...\n', sideStr);
                warning('on', 'backtrace');
            end
            fclose(fopen([outputDir, filesep, 'skip_', sideCode, '.txt'], 'w'));
            runStatusMultiSource(source_index,side+1) = 1;
            continue;
        end

        % skip stimSets if not provided for this side
        if settings.stimSetMode && ~isfile(strcat(outputDir, filesep, 'Current_protocols_',string(side),'.csv'))
            warning('off', 'backtrace');
            warning('No stimulation set for %s side! Skipping...\n', sideStr);
            warning('on', 'backtrace');
            fclose(fopen([outputDir, filesep, 'skip_', sideCode, '.txt'], 'w'));
            runStatusMultiSource(source_index,side+1) = 1;
            continue;
        end
    
        % skip PAM if no fibers were preserved for the stim protocol
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
            if prob_PAM 

                % parameters for "probabilistic" run
                % here we iterate over different fiber diameters 2.0 - 4.0 um
                N_samples = 5;
                parameter_limits = [2.0,4.0];
                parameter_step = (parameter_limits(2)-parameter_limits(1)) / (N_samples - 1);
                scaling = -1.0;
                % one can swap it to probabilistic 
            else
                N_samples = 1;
            end
        else
            N_samples = 1;
        end
    
        %% OSS-DBS part (using the corresponding conda environment)
        for i = 1:N_samples

            if settings.calcAxonActivation
                if prob_PAM 
                    settings.fiberDiameter = options.prefs.machine.vatsettings.butenko_fiberDiameter;
                    settings.fiberDiameter(:) = parameter_limits(1);
                    settings.fiberDiameter(:) = settings.fiberDiameter(:) + (i-1)*parameter_step;

                    parameterFile = fullfile(outputDir, 'oss-dbs_parameters.mat');
                    save(parameterFile, 'settings', '-v7.3');
                end
        
                % clean-up
                ea_delete([outputDir, filesep, 'Allocated_axons.h5']);
                ea_delete([folder2save,filesep,'oss_time_result.h5'])
    
                system(['python ', ea_getearoot, 'ext_libs/OSS-DBS/Axon_Processing/axon_allocation.py ', outputDir,' ', num2str(side), ' ', parameterFile]);
            end

            % prepare OSS-DBS input as oss-dbs_parameters.json
            system(['leaddbs2ossdbs --hemi_side ', num2str(side), ' ', parameterFile, ...
                ' --output_path ', outputDir]);
            parameterFile_json = [parameterFile(1:end-3), 'json'];
    
            % run OSS-DBS
            system(['ossdbs ', parameterFile_json]);
        
            % prepare NEURON simulation
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

                if prob_PAM
                    system(['python ', ea_getearoot, 'ext_libs/OSS-DBS/Axon_Processing/PAM_caller.py ', neuron_folder, ' ', folder2save,' ', timeDomainSolution, ' ', pathwayParameterFile, ' ', num2str(scaling), ' ', num2str(i)]);
                else
                    system(['python ', ea_getearoot, 'ext_libs/OSS-DBS/Axon_Processing/PAM_caller.py ', neuron_folder, ' ', folder2save,' ', timeDomainSolution, ' ', pathwayParameterFile]);
                end
            end
        end

        if prob_PAM
            % convert binary PAM status over uncertain parameter to "probabilistic activation"
            ea_get_probab_axon_state(folder2save,1,strcmp(settings.butenko_intersectStatus,'activated'));
        end

        %% Postprocessing in Lead-DBS

        % Check if OSS-DBS calculation is finished
        while ~isfile([outputDir, filesep, 'success_', sideCode, '.txt']) ...
                && ~isfile([outputDir, filesep, 'fail_', sideCode, '.txt'])
            continue;
        end

        % clean-up if outOfCore was used
        if settings.outOfCore == 1
            ea_delete([outputDir, filesep, 'Results_',sideCode,filesep,'oss_freq_domain_tmp_PAM.hdf5']);
        end
    
        if isfile([outputDir, filesep, 'success_', sideCode, '.txt'])
            runStatusMultiSource(source_index,side+1) = 1;
            fprintf('\nOSS-DBS calculation succeeded!\n\n')
    
            if settings.exportVAT
                [stimparams(side+1).VAT.VAT,stimparams(side+1).volume,source_efields{side+1,source_use_index},source_vtas{side+1,source_use_index}] = ea_convert_ossdbs_VTAs(options,settings,side,multiSourceMode,source_use_index,outputDir,outputBasePath,templateOutputBasePath);
            end

            if settings.calcAxonActivation
                ea_convert_ossdbs_axons(settings,side,resultfig,outputDir,outputBasePath,templateOutputBasePath)
            end

        elseif isfile([outputDir, filesep, 'fail_', sideCode, '.txt'])
            fprintf('\n')
            warning('off', 'backtrace');
            warning('OSS-DBS calculation failed for %s side!\n', sideStr);
            warning('on', 'backtrace');
            %runStatus(side+1) = 0;
        end
   
    end

    % check only the first source for PAM
    if settings.calcAxonActivation
        runStatusMultiSource(2:end,:) = 1;
        break
    end

end

% here we merge and display multiSourceModes
for side = 0:1
    if multiSourceMode(side+1) && nActiveSources(side+1,:) > 0
        stimparams = ea_postprocess_multisource(options,settings,side+1,source_efields,source_vtas);
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

%% Helper function to fail exit
function [runStatus, stimparameters] = ea_exit_genvat_butenko()

    runStatus = [0 0];  % runStatus
    stimparameters = struct(); % empty for stimparameters

    % Restore working directory and environment variables
    setenv('LD_LIBRARY_PATH', ''); % Clear LD_LIBRARY_PATH to resolve conflicts
    setenv('LD_LIBRARY_PATH', getenv('LD_LIBRARY_PATH'));
    setenv('PATH', getenv('PATH'));


