function varargout = ea_genvat_butenko(varargin)
% Wrapper for OSS-DBS simulations.
% By Butenko and Li, konstantinmgtu@gmail.com

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
% get resultfig handle
if exist('hFigure', 'var')
    resultfig = getappdata(hFigure,'resultfig');
else
    resultfig = [];
end

time = datetime('now', 'TimeZone', 'local');
timezone = time.TimeZone;
setenv('TZ', timezone);

% import settings from Lead-DBS GUI
settings = ea_prepare_ossdbs(options);

% some hardcoded parameters, can be added to GUI later
prepFiles_cluster = 0; % set to 1 if you only want to prep files for cluster comp.
true_VTA = 0; % set to 1 to compute classic VAT using axonal grids
prob_PAM = 0; % set to 1 to compute PAM over an uncertain parameter (e.g. fiber diameter)
settings.outOfCore = 0;  % set to 1 if RAM capacity is exceeded during PAM

if settings.stimSetMode
    ea_warndlg("Not yet supported in V2")
    return
end

% set outputs
outputPaths = ea_get_oss_outputPaths(options,S);

% segment MRI image
settings = ea_segment_MRI(options, settings, outputPaths);

% prepare tensor data
settings.DTI_data_name = ea_prepare_DTI(options,outputPaths);

% get electrode reconstruction parameters in OSS-DBS format
settings = ea_get_oss_reco(options, settings);

%% Stimulation Information
% Check which sources are active and insert blank stimulation for unilateral cases
[S, settings, activeSources] = ea_check_stimSources(options,S,settings);
nActiveSources = [nnz(~isnan(activeSources(1,:))), nnz(~isnan(activeSources(2,:)))];
multiSourceMode = [nActiveSources(1) > 1; nActiveSources(2) > 1];
runStatusMultiSource = zeros(4,2);  % check status for each source
source_efields = cell(2,4);  % temp files to store results for each source
source_vtas = cell(2,4);
stimparams = struct();

% multiple sources are not supported for PAM
if any(nActiveSources > 1)
    if options.prefs.machine.vatsettings.butenko_calcPAM
        ea_warndlg('MultiSource Mode is not supported for PAM!')
        % Restore working directory and environment variables
        [varargout{1}, varargout{2}] = ea_exit_genvat_butenko();
        return
    else
        ea_warndlg('MultiSource Mode is used! Stimulation volumes will be computed separately and merged using max(||E||)')
    end
end

% if single source, we will run only one iteration
for source_index = 1:4

    % get stim settings for particular source    
    settings = ea_get_stimProtocol(options,S,settings,activeSources,source_index);

    % to optimize monopolar
    % solve for 10 mA (anodic to avoid sign confusion later)
    settings.Phi_vector(:,1) = 10;
    settings.Phi_vector(:,2:end) = 0;
    
    if settings.calcAxonActivation
        % will exit after the first source
        if true_VTA
            % custom case of PAM is original VAT
            % should be tested for unilateral!
            fibersFound = [0,0];
            if ~isnan(activeSources(1,source_index))
                settings = ea_switch2VATgrid(options,S,settings,0,outputPaths);
                fibersFound(:,1) = 1; 
            end
            if ~isnan(activeSources(2,source_index))
                settings = ea_switch2VATgrid(options,S,settings,1,outputPaths);
                fibersFound(:,2) = 1;
            end
        else
            % warp fibers, remove too short and too far away ones
            [settings,fibersFound] = ea_prepare_fibers(options, S, settings, outputPaths);
        end
    end
    
    % Save settings for OSS-DBS
    if ~all(isnan(activeSources(:,source_index)))
        parameterFile = ea_save_ossdbs_settings(options, S, settings, outputPaths);
    end
    
    if prepFiles_cluster == 1
        % now you can run OSS-DBS externally
        [varargout{1}, varargout{2}] = ea_exit_genvat_butenko();
        return
    end
    
    % Iterate over hemispheres: 0 - rh , 1 - lh
    for side = 0:1

        switch side
            case 0
                sideCode = 'rh';
                sideStr = 'right';
            case 1
                sideCode = 'lh';
                sideStr = 'left';
        end

        if ~multiSourceMode(side+1)
            % not relevant in this case, terminate after one iteration;
            source_use_index = 5;  
        else
            source_use_index = source_index; 
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

            fclose(fopen([outputPaths.outputDir, filesep, 'skip_', sideCode, '.txt'], 'w'));
            runStatusMultiSource(source_index,side+1) = 1;
            continue;
        end

        % skip stimSets if not provided for this side
        if settings.stimSetMode && ~isfile(strcat(outputPaths.outputDir, filesep, 'Current_protocols_',string(side),'.csv'))
            warning('off', 'backtrace');
            warning('No stimulation set for %s side! Skipping...\n', sideStr);
            warning('on', 'backtrace');
            fclose(fopen([outputPaths.outputDir, filesep, 'skip_', sideCode, '.txt'], 'w'));
            runStatusMultiSource(source_index,side+1) = 1;
            continue;
        end
    
        % skip PAM if no fibers were preserved for the stim protocol
        if settings.calcAxonActivation && ~any(fibersFound(:,side+1))
            warning('off', 'backtrace');
            warning('No fibers found for %s side! Skipping...\n', sideStr);
            warning('on', 'backtrace');
            fclose(fopen([outputPaths.outputDir, filesep, 'skip_', sideCode, '.txt'], 'w'));
            runStatusMultiSource(source_index,side+1) = 1;
            continue;
        end
    
        fprintf('\nRunning OSS-DBS for %s side stimulation...\n\n', sideStr);
    
        N_samples = 1;  % by default, run one instance
        if settings.calcAxonActivation && prob_PAM
            % "probabilistic" run
            % here we iterate over different fiber diameters 2.0 - 4.0 um
            N_samples = 5;  % hardcoded for now
            % one can swap it to probabilistic values
        end
    
        %% OSS-DBS part (using the corresponding conda environment)
        for i = 1:N_samples

            if settings.calcAxonActivation
                if prob_PAM 
                    settings = ea_updatePAM_parameter(options,settings,N_samples,outputPaths,i);
                    scaling = 1.0; % same current scaling across the parameter sweep
                end
        
                % clean-up
                folder2save = [outputPaths.outputDir,filesep,'Results_', sideCode];
                ea_delete([outputPaths.outputDir, filesep, 'Allocated_axons.h5']);
                ea_delete([folder2save,filesep,'oss_time_result.h5'])
    
                % allocate computational axons on fibers
                system(['python ', ea_getearoot, 'ext_libs/OSS-DBS/Axon_Processing/axon_allocation.py ', outputPaths.outputDir,' ', num2str(side), ' ', parameterFile]);
            end

            % prepare OSS-DBS input as oss-dbs_parameters.json
            system(['leaddbs2ossdbs --hemi_side ', num2str(side), ' ', parameterFile, ...
                ' --output_path ', outputPaths.outputDir]);
            parameterFile_json = [parameterFile(1:end-3), 'json'];
    
            % run OSS-DBS
            system(['ossdbs ', parameterFile_json]);
        
            % prepare NEURON simulation
            if settings.calcAxonActivation
                % copy NEURON folder to the stimulation folder 
                leaddbs_neuron = [ea_getearoot, 'ext_libs/OSS-DBS/Axon_Processing/Axon_files'];
                neuron_folder = fullfile(outputPaths.outputDir,'Axon_files');
                copyfile(leaddbs_neuron, neuron_folder)
        
                % call the NEURON module
                timeDomainSolution = [outputPaths.outputDir,filesep,'Results_', sideCode, filesep, 'oss_time_result_PAM.h5'];
                pathwayParameterFile = [outputPaths.outputDir,filesep, 'Allocated_axons_parameters.json'];
    
                % check if the time domain results is available
                if ~isfile(timeDomainSolution)
                    ea_warndlg('OSS-DBS failed to prepare a time domain solution. If RAM consumption exceeded the hardware limit, set settings.outOfCore to 1')
                    return
                end

                if prob_PAM
                    system(['python ', ea_getearoot, 'ext_libs/OSS-DBS/Axon_Processing/PAM_caller.py ', neuron_folder, ' ', folder2save,' ', timeDomainSolution, ' ', pathwayParameterFile, ' ', num2str(scaling), ' ', num2str(i)]);
                else
                    % instead of PAM_caller, we call optimization algorithm here
                    PAM_caller_script = [ea_getearoot, 'ext_libs/OSS-DBS/Axon_Processing/PAM_caller.py'];
                    optim_settings_dict = [outputPaths.outputDir,filesep,'netblend_dict.json'];
                    system(['python ', ea_getearoot, 'cleartune/PathwayTune/pam_optimizer.py ', PAM_caller_script, ' ', neuron_folder, ' ', optim_settings_dict, ' ', outputPaths.outputDir, ' ', num2str(side)])

                    %system(['python ', ea_getearoot, 'ext_libs/OSS-DBS/Axon_Processing/PAM_caller.py ', neuron_folder, ' ', folder2save,' ', timeDomainSolution, ' ', pathwayParameterFile]);
                end
            end
        end

        %% Postprocessing in Lead-DBS

        % Check if OSS-DBS calculation is finished
        while ~isfile([outputPaths.outputDir, filesep, 'success_', sideCode, '.txt']) ...
                && ~isfile([outputPaths.outputDir, filesep, 'fail_', sideCode, '.txt'])
            continue;
        end

        if prob_PAM
            % convert binary PAM status over uncertain parameter to "probabilistic activations"
            ea_get_probab_axon_state(folder2save,1,strcmp(settings.butenko_intersectStatus,'activated'));
        end

        % clean-up time domain solution if outOfCore was used
        if settings.outOfCore == 1
            ea_delete([outputPaths.outputDir, filesep, 'Results_',sideCode,filesep,'oss_freq_domain_tmp_PAM.hdf5']);
        end
    
        if isfile([outputPaths.outputDir, filesep, 'success_', sideCode, '.txt'])
            runStatusMultiSource(source_index,side+1) = 1;
            fprintf('\nOSS-DBS calculation succeeded!\n\n')
    
            % prepare Lead-DBS BIDS format VATs
            if settings.exportVAT
                [stimparams(side+1).VAT.VAT,stimparams(side+1).volume,source_efields{side+1,source_use_index},source_vtas{side+1,source_use_index}] = ea_convert_ossdbs_VTAs(options,settings,side,multiSourceMode,source_use_index,outputPaths);
            end

            % prepare Lead-DBS BIDS format fiber activations
            if settings.calcAxonActivation
                ea_convert_ossdbs_axons(options,settings,side,prob_PAM,resultfig,outputPaths)
            end

        elseif isfile([outputPaths.outputDir, filesep, 'fail_', sideCode, '.txt'])
            fprintf('\n')
            warning('off', 'backtrace');
            warning('OSS-DBS calculation failed for %s side!\n', sideStr);
            warning('on', 'backtrace');
            %runStatus(side+1) = 0;
        end
   
    end

    % check only the first source for PAM
    if settings.calcAxonActivation || all(~multiSourceMode)
        runStatusMultiSource(2:end,:) = 1;
        break
    end

end

% merge multisource VATs
for side = 0:1
    if multiSourceMode(side+1) && nActiveSources(1,side+1) > 0
        stimparams = ea_postprocess_multisource(options,settings,side+1,source_efields,source_vtas,outputPaths);
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


