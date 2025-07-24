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
    varargout{3} = true; % Support estimation in native space
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
%options.stimSetMode = 0;
[settings,S] = ea_prepare_ossdbs(options,S);
env = ea_conda_env('OSS-DBSv2');

% some hardcoded parameters, can be added to the GUI later
settings.reuse_warped_connectome = 0;  % set to 1 if the connectome was already processed for the given stim. settings
prepFiles_cluster = 0; % set to 1 if you only want to prep files for cluster comp.
true_VTA = 0; % set to 1 to compute classic VAT using axonal grids
settings.outOfCore = 0; % set to 1 if RAM capacity is exceeded during PAM

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

% check if multisource mode
if any(nActiveSources > 1)
    if options.prefs.machine.vatsettings.butenko_calcPAM
        % ea_warndlg('MultiSource Mode is not supported for PAM!')
        % % Restore working directory and environment variables
        % [varargout{1}, varargout{2}] = ea_exit_genvat_butenko();
        % return
        ea_warndlg('MultiSource Mode is used! Disjunction of PAM results across sources will be stored!')
    else
        ea_warndlg('MultiSource Mode is used! Stimulation volumes will be computed separately and merged using max(||E||)')
    end
    first_active_source = 1;  % check all sources
else
    first_active_source_rh = find(~isnan(activeSources(1,:)),1,'first');
    first_active_source_lh = find(~isnan(activeSources(2,:)),1,'first');
    first_active_source = min([first_active_source_rh,first_active_source_lh]);
end

% if single source, we will run only one iteration
for source_index = first_active_source:4

    % get stim settings for particular source
    settings = ea_get_stimProtocol(options,S,settings,activeSources,source_index);

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
            % TBD: skip inactive sources!
            [settings,fibersFound] = ea_prepare_fibers(options, S, settings, outputPaths,  source_index);
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

        % hemisphere specific data is stored in this subfolder
        outputPaths.HemiSimFolder = [outputPaths.outputDir, filesep, 'OSS_sim_files_', sideCode];

        if ~multiSourceMode(side+1)
            % not relevant in this case, terminate after one iteration;
            source_use_index = 5;
        else
            source_use_index = source_index;
        end

        % copy Current_protocols if generated externally
        if settings.optimizer || settings.trainANN
            if ~exist(outputPaths.HemiSimFolder,'dir')
                mkdir(outputPaths.HemiSimFolder)
            end
            copyfile([outputPaths.outputDir,filesep,'NB_',sideCode,filesep,'Current_protocols_',num2str(side),'.csv'],[outputPaths.HemiSimFolder,filesep,'Current_protocols_',num2str(side),'.csv'])
        elseif settings.stimSetMode
            if ~exist(outputPaths.HemiSimFolder,'dir')
                mkdir(outputPaths.HemiSimFolder)
            end
            copyfile([outputPaths.outputDir,filesep,'Current_protocols_',num2str(side),'.csv'],[outputPaths.HemiSimFolder,filesep,'Current_protocols_',num2str(side),'.csv'])
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
        if settings.stimSetMode && ~isfile(strcat(outputPaths.HemiSimFolder,  filesep, 'Current_protocols_',string(side),'.csv'))
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

        if ~exist(outputPaths.HemiSimFolder,'dir')
            mkdir(outputPaths.HemiSimFolder)
        end

        %% OSS-DBS part (using the corresponding conda environment)
        for i = 1:settings.N_samples  % mutiple samples if probablistic PAM is used, otherwise 1

            if settings.calcAxonActivation
                if settings.prob_PAM && (source_use_index == 1 || source_use_index == 5)
                    settings = ea_updatePAM_parameter(options,settings,outputPaths,i);
                    if any(multiSourceMode)
                        vatsettings = options.prefs.machine.vatsettings;
                        pparam = vatsettings.butenko_probabilistic_parameter;
                        copyfile([outputPaths.HemiSimFolder, filesep, pparam,'_samples.mat'], [outputPaths.outputDir, filesep, pparam,'_samples_',sideCode,'.mat'])
                    end
                elseif settings.prob_PAM
                    % for other sources, load already sampled parameters
                    settings = ea_load_prob_parameter(options, settings, outputPaths, sideCode, i);
                end

                % clean-up
                ea_delete([outputPaths.HemiSimFolder, filesep, 'Allocated_axons.h5']);
                ea_delete([outputPaths.HemiSimFolder, filesep, 'Results', filesep,'oss_time_result*'])

                % allocate computational axons on fibers
                %system(['python ', ea_getearoot, 'ext_libs/OSS-DBS/Axon_Processing/axon_allocation.py ', outputPaths.outputDir,' ', num2str(side), ' ', parameterFile]);
                env.system(['prepareaxonmodel ',ea_path_helper(outputPaths.outputDir),' --hemi_side ',num2str(side),' --description_file ', ea_path_helper(parameterFile)]);
            end

            % prepare OSS-DBS input as oss-dbs_parameters.json
            env.system(['leaddbs2ossdbs --hemi_side ', num2str(side), ' ', ea_path_helper(parameterFile), ...
                ' --output_path ', ea_path_helper(outputPaths.HemiSimFolder)]);
            [~,input_name,~] = fileparts(parameterFile);
            parameterFile_json = [outputPaths.HemiSimFolder, filesep, input_name, '.json'];

            % run OSS-DBS
            [~, cmdout] = env.system(['ossdbs ', ea_path_helper(parameterFile_json)])

            % detec error related to Bnd_Box
            if contains(cmdout, 'Bnd_Box is void')
                disp ('Error "Bnd_Box is void" detected, increasing the dimensions ...');

                % increase the Bnd_Box dimensions
                env.system(cell2mat(['python ' ea_regexpdir(ea_getearoot, 'BndBoxDimensionsEdits.py') ' ', ea_path_helper(parameterFile_json)]));

                % run OSS-DBS
                env.system(['ossdbs ', ea_path_helper(parameterFile_json)])
            end

            % prepare NEURON simulation
            if settings.calcAxonActivation

                if ~strcmp(settings.butenko_intersectStatus,'damaged')
                    % we additionally correct for the tissue push and
                    % downscale the solution (equivalent of pulling VTAs into the electrode volume)
                    scaling = 0.80;  % estimate for our default comp. domain, see eq. for el. potential in co-axial cables
                else
                    scaling = 1.0;
                end

                % check if the time domain results is available
                timeDomainSolution = [outputPaths.HemiSimFolder,filesep,'Results', filesep, 'oss_time_result_PAM.h5'];
                if ~isfile(timeDomainSolution) && ~settings.stimSetMode
                    ea_warndlg('OSS-DBS failed to prepare a time domain solution. If RAM consumption exceeded the hardware limit, set settings.outOfCore to 1')
                    return
                end

                if settings.optimizer
                    env.system(['python ', ea_getearoot,'ext_libs',filesep,'PathwayTune',filesep,'pam_optimizer.py ', settings.netblend_settings_file, ' ', ea_path_helper(outputPaths.outputDir), ' ', num2str(side), ' ', ea_path_helper(parameterFile_json), ' ', num2str(scaling)])
                else
                    if settings.prob_PAM
                        %system(['python ', ea_getearoot, 'ext_libs/OSS-DBS/Axon_Processing/PAM_caller.py ', neuron_folder, ' ', folder2save,' ', timeDomainSolution, ' ', pathwayParameterFile, ' ', num2str(scaling), ' ', num2str(i)]);
                        env.system(['run_pathway_activation ', ea_path_helper(parameterFile_json), ' --scaling_index ', num2str(i), ' --scaling ', num2str(scaling)]);
                    else
                        %system(['python ', ea_getearoot, 'ext_libs/OSS-DBS/Axon_Processing/PAM_caller.py ', neuron_folder, ' ', folder2save,' ', timeDomainSolution, ' ', pathwayParameterFile]);
                        env.system(['run_pathway_activation ', ea_path_helper(parameterFile_json), ' --scaling ', num2str(scaling)]);
                    end
                end

                % remove the large file containing the time-domain solution (but not for StimSets!)
                ea_delete([outputPaths.HemiSimFolder, filesep, 'Results', filesep,'oss_time_result*'])
            end
        end

        %% Postprocessing in Lead-DBS

        % Check if OSS-DBS calculation is finished
        while ~isfile([outputPaths.HemiSimFolder, filesep, 'success_', sideCode, '.txt']) ...
                && ~isfile([outputPaths.HemiSimFolder, filesep, 'fail_', sideCode, '.txt'])
            continue;
        end

        if settings.prob_PAM && all(~multiSourceMode)
            % convert binary PAM status over uncertain parameter to "probabilistic activations"
            ea_get_probab_axon_state([outputPaths.HemiSimFolder,filesep,'Results'],1,settings,side);
        end

        % clean-up time domain solution if outOfCore was used
        if settings.outOfCore == 1
            ea_delete([outputPaths.HemiSimFolder, filesep, 'Results', filesep, 'oss_freq_domain_tmp_PAM.hdf5']);
        end

        if isfile([outputPaths.HemiSimFolder, filesep, 'success_', sideCode, '.txt'])
            runStatusMultiSource(source_index,side+1) = 1;
            fprintf('\nOSS-DBS calculation succeeded!\n\n')

            % prepare Lead-DBS BIDS format VATs
            if settings.exportVAT && settings.optimizer
                % get 4-D unit niftis for the optimizer and exit
                ea_convert_ossdbs_StimSets_VTAs(settings,side,outputPaths)
                stimparams(side+1).VAT.VAT = -42.0;
                stimparams(side+1).volume = -42.0;
                ea_exit_genvat_butenko;
            elseif settings.exportVAT
                [stimparams(side+1).VAT.VAT,stimparams(side+1).volume,source_efields{side+1,source_use_index},source_vtas{side+1,source_use_index}] = ea_convert_ossdbs_VTAs(options,settings,side,multiSourceMode,source_use_index,outputPaths);
            end

            if settings.prob_PAM && any(multiSourceMode)
                % for multisource, we will convert in the external loop
                % we just need to add the source index to the Axon States
                axonStateFolder = ea_sourceIndex4AxonStates(outputPaths, side, source_use_index);
                continue
            end

            % prepare Lead-DBS BIDS format fiber activations
            if settings.calcAxonActivation && ~settings.optimizer
                ea_convert_ossdbs_axons(options,settings,side,settings.prob_PAM,resultfig,outputPaths,source_use_index);
            end

        elseif isfile([outputPaths.HemiSimFolder, filesep, 'fail_', sideCode, '.txt'])
            fprintf('\n')
            warning('off', 'backtrace');
            warning('OSS-DBS calculation failed for %s side!\n', sideStr);
            warning('on', 'backtrace');
            %runStatus(side+1) = 0;
        end

        % clean-up StimSets FEM solutions
        ea_delete([outputPaths.HemiSimFolder, filesep, 'ResultsE*']);

    end

    % check only the first source for PAM
    if all(~multiSourceMode)
        runStatusMultiSource(2:end,:) = 1;
        break
    end

end

% merge multisource VATs
for side = 0:1
    if multiSourceMode(side+1) && nActiveSources(1,side+1) > 0 && ~settings.calcAxonActivation
        %stimparams = ea_postprocess_multisource(options,settings,side+1,source_efields,source_vtas,outputPaths);
        [vatfv,vatvolume] = ea_postprocess_multisource(options,settings,side+1,source_efields,source_vtas,outputPaths);
        stimparams(side+1).VAT.VAT = vatfv;
        stimparams(side+1).volume = vatvolume;
    elseif multiSourceMode(side+1) && nActiveSources(1,side+1) > 0 && settings.calcAxonActivation && ~settings.prob_PAM
        ea_postprocess_multisource_pam(options,settings,side+1)
    end
end

% special case for multisource probabilistic PAM
if any(multiSourceMode) && settings.prob_PAM
    ea_get_prob_fiber_states_for_multisource(options,settings,outputPaths,axonStateFolder,resultfig)
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


