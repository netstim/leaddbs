function parameterFile = ea_save_ossdbs_settings(options,S,settings,outputDir,templateOutputDir)
%if any(~isnan(activeSources))
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
%end