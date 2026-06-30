function parameterFile = ea_save_ossdbs_settings(options,S,settings,outputPaths)
% Save simulation parameters in oss-dbs_parameters.mat and do a clean-up.
% By Butenko and Li, konstantinmgtu@gmail.com

arguments
    options     % Lead-DBS options for electrode reconstruction and stimulation
    S           % Lead-DBS stimulation settings
    settings    % parameters for OSS-DBS simulation
    outputPaths % various paths to conform with lead-dbs BIDS structure
end

% Skip destructive cleanup (and the parameter-file overwrite) when only
% postprocessing existing cluster results — those folders contain the
% Axon_state_*.mat / Pathway_status_*.json we need to keep.
postprocess_cluster = isfield(settings,'postprocess_cluster') && settings.postprocess_cluster;

parameterFile = fullfile(outputPaths.outputDir, 'oss-dbs_parameters.mat');
if ~postprocess_cluster
    save(parameterFile, 'settings', '-v7.3');
end
ea_savestimulation(S, options);
if options.native
    poptions = options;
    poptions.native = 0;
    ea_savestimulation(S, poptions);
end

if ~postprocess_cluster
    % Delete previous results from stimSetMode
    ea_delete([outputPaths.outputDir, filesep, 'Result_StimProt_*']);
    if options.native
        ea_delete([outputPaths.templateOutputDir, filesep, 'Result_StimProt_*']);
    end

    % full clean-up for V2 (Lead-DBS outputs are still kept!)
    ea_delete([outputPaths.outputDir, filesep, 'OSS_sim_files*']);
    ea_delete([outputPaths.outputDir, filesep, 'skip_rh.txt']);
    ea_delete([outputPaths.outputDir, filesep, 'skip_lh.txt']);
end

%% Run OSS-DBS
setenv('LD_LIBRARY_PATH', ''); % Clear LD_LIBRARY_PATH to resolve conflicts

%end