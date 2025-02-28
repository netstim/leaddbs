function parameterFile = ea_save_ossdbs_settings(options,S,settings,outputPaths)
% Save simulation parameters in oss-dbs_parameters.mat and do a clean-up.
% By Butenko and Li, konstantinmgtu@gmail.com

arguments
    options     % Lead-DBS options for electrode reconstruction and stimulation
    S           % Lead-DBS stimulation settings
    settings    % parameters for OSS-DBS simulation
    outputPaths % various paths to conform with lead-dbs BIDS structure
end

parameterFile = fullfile(outputPaths.outputDir, 'oss-dbs_parameters.mat');
save(parameterFile, 'settings', '-v7.3');
ea_savestimulation(S, options);
if options.native
    poptions = options;
    poptions.native = 0;
    ea_savestimulation(S, poptions);
end

% Delete previous results from stimSetMode
ea_delete([outputPaths.outputDir, filesep, 'Result_StimProt_*']);
if options.native
    ea_delete([outputPaths.templateOutputDir, filesep, 'Result_StimProt_*']);
end

% full clean-up for V2 (Lead-DBS outputs are still kept!)
ea_delete([outputPaths.outputDir, filesep, 'OSS_sim_files*']);
ea_delete([outputPaths.outputDir, filesep, 'skip_rh.txt']);
ea_delete([outputPaths.outputDir, filesep, 'skip_lh.txt']);

%% Run OSS-DBS
setenv('LD_LIBRARY_PATH', ''); % Clear LD_LIBRARY_PATH to resolve conflicts

%end