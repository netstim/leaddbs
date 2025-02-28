function ea_get_prob_fiber_states_for_multisource(options,settings,outputPaths,axonStateFolder,resultfig)

% convert PAM samples from multiple sources to prob. PAM
% IMPORTANT:this is a special case for multisource prob. PAM
% By Butenko, konstantinmgtu@gmail.com

arguments
    options             % Lead-DBS options for electrode reconstruction and stimulation
    settings            % parameters for OSS-DBS simulation
    outputPaths         % various paths to conform with lead-dbs BIDS structure 
    axonStateFolder     % temp. folder where we store AxonStates from different samples and sources
    resultfig           % figure handle
end

ea_postprocess_multisource_axonstates(axonStateFolder,options,settings)
% then create probabilistic fiber states
ea_get_probab_axon_state_for_multisource(axonStateFolder,strcmp(settings.butenko_intersectStatus,'activated'));
for side = 0:1
    switch side
        case 0
            sideCode = 'rh';
        case 1
            sideCode = 'lh';
    end
    outputPaths.HemiSimFolder = [outputPaths.outputDir, filesep, 'OSS_sim_files_', sideCode];
    if ~isfolder(outputPaths.HemiSimFolder)
        % if one source on one side
        mkdir(outputPaths.HemiSimFolder)
        mkdir([outputPaths.HemiSimFolder,filesep,'Results'])
    end
    % move them back to Results (or change hardwired folder in ea_convert_ossdbs_axons)
    axon_state_files = dir([axonStateFolder,filesep,'Axon_state_*']); 
    for file_i = 1:length(axon_state_files)
        % poor approach 
        if (side == 1 && (contains(axon_state_files(file_i).name, '_right_') || contains(axon_state_files(file_i).name, '_rh_'))) || (side == 0 && (contains(axon_state_files(file_i).name, '_left_') || contains(axon_state_files(file_i).name, '_lh_')))
            continue
        else
            oldfile = [axon_state_files(file_i).folder,filesep,axon_state_files(file_i).name];
            newfile = [outputPaths.HemiSimFolder,filesep,'Results',filesep,axon_state_files(file_i).name];
            movefile(oldfile,newfile)
        end
    end
    ea_convert_ossdbs_axons(options,settings,side,settings.prob_PAM,resultfig,outputPaths,5);
end
