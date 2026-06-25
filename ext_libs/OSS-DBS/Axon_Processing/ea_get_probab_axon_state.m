function ea_get_probab_axon_state(results_folder,plot_rates,settings,side)
% Estimate probability of axon activation based on a sweep of parameters
% (e.g. fiber diameters).
% By Butenko and Li, konstantinmgtu@gmail.com

arguments
    results_folder      % path to save the output
    plot_rates          {mustBeNumericOrLogical}% if true, plot activation rates over parameter sweep
    settings            % parameters for OSS-DBS simulation
    side                {mustBeNumeric} % hemisphere index (0 - rh, 1 - lh)
end

% check which pathways were simulated and their percent activations
rate_files = dir([results_folder,filesep,'Pathway_status*']);
pathways = {};
activations_over_pathways = {};
pt_counter = 0;
% we enforce correct loading order
for file_i = 1:length(rate_files)
    parts = strsplit(rate_files(file_i).name,'_');
    % very stupid way to recover the pathway name
    pathway_name_orig = [];
    for i = 3:length(parts)-1
        if i == 3
            pathway_name_orig = [pathway_name_orig,parts{i}];
        else
            pathway_name_orig = [pathway_name_orig,'_',parts{i}];
        end
    end
    if ~any(strcmp(pathways,pathway_name_orig))
        scaling_counter = 1;
    end
    rate_file = fullfile(rate_files(file_i).folder, ['Pathway_status_',pathway_name_orig,'_',num2str(scaling_counter),'.json']);
    scaling_counter = scaling_counter + 1;
    jsonText = fileread(rate_file);
    disp(rate_file)
    % Convert JSON formatted text to MATLAB data types 
    jsonDict = jsondecode(jsonText); 
    %pathway_name = strrep(pathway_name_orig,'_',' ');
    index = rate_files(file_i).name(end-5);
    if ~any(strcmp(pathways,pathway_name_orig))
        pt_counter = pt_counter + 1;
        pathways{pt_counter} = pathway_name_orig;
        % scalings are always ordered
        activations_over_pathways{pt_counter} = [jsonDict.percent_activated];
    else
        activations_over_pathways{pt_counter} = [activations_over_pathways{pt_counter},jsonDict.percent_activated];
    end
end

N_samples = length(activations_over_pathways{pt_counter});

% plot activation curves over the parameter sweep
if plot_rates
    figure()
    for pt_counter = 1:length(pathways)
        pathway_name = strrep(pathways{pt_counter},'_',' ');
        disp(pathway_name)
        plot(1:N_samples,activations_over_pathways{pt_counter},'DisplayName',pathway_name)
        ylim([0,100]);
        xlabel('Scaling')
        ylabel('Percent Activation')
        hold on 
    end
    legend('Location','eastoutside')
end

% now let's compute activation probabilty for each fiber
% we compute probabilities for the filtered fibers and will map to global
% indices later in ea_genvat_butenko
for pt_counter = 1:length(pathways)
    skip_pathway = false;
    % iterate over scalings (fiber diameters)
    % important: number of compartments can differ based on the fiber
    % diameter / length!
    for sample_i = 1:N_samples

        % If stimSetMode, extract the index from tractName (but axonState is still checked on the indexed file)
        if settings.stimSetMode && ~settings.optimizer
            if startsWith(settings.connectome, 'Multi-Tract: ')
                stimProt_index = regexp(tractName, '(?<=_)\d+$', 'match', 'once');
                tractName = regexp(tractName, '.+(?=_\d+$)', 'match', 'once');
            else
                stimProt_index = regexp(axonState{f}, '(?<=Axon_state_)\d+(?=\.mat$)', 'match', 'once');
            end

            % in this case, each sample is in a separate folder
            Axon_state_gpe2stn_sm_left_0   % without prob

        else
            Axon_state_file = fullfile(results_folder, ['Axon_state_',pathways{pt_counter},'_',num2str(sample_i),'.mat']);
        end

        try
            load(Axon_state_file, 'fibers','ea_fibformat','connectome_name');
        catch ME
            % Corrupt / non-MAT cluster output for this pathway. Skip rather
            % than crash -- downstream merge code treats a missing *_prob.mat
            % as zero activation.
            fprintf('Skipping pathway %s: failed to load sample %d (%s).\n     file: %s\n', ...
                pathways{pt_counter}, sample_i, ME.message, Axon_state_file);
            skip_pathway = true;
            break;
        end
        if isempty(fibers)
            % No axons survived for this pathway (e.g. all pre-filtered by
            % Kuncel VTA). Skip the whole pathway -- downstream merge code
            % treats a missing *_prob.mat as zero activation.
            fprintf('Skipping pathway %s: empty fibers in sample %d (no *_prob.mat will be written).\n', ...
                pathways{pt_counter}, sample_i);
            skip_pathway = true;
            break;
        end
        n_comp = sum(bsxfun(@eq, fibers(:,4), fibers(:,4)'), 2)';
        % the number is the same since these are OSS-DBS axons
        idx = n_comp(1) * ones(length(unique(fibers(:,4))),1);

        if strcmp(settings.butenko_intersectStatus,'activated_at_active_contacts')
            fibers = OSS_DBS_Damaged2Activated(settings,fibers,idx,side+1);
        end

        % intialize new axon state file with probabilistic activations
        % morphologz defined by the first scaling!
        if sample_i == 1
            fibers_prob = fibers;
            % unnecessary, but for clarity
            fibers_prob(:,5) = 0;
        end

        fibers_state = zeros(max(fibers(:,4)),1);
        for fiber_i = 1:max(fibers(:,4))   % relative index!
            idx_comp = find(fibers(:,4)==fiber_i);
            % we need to check only status of one compartment
            % here we only count really activated ones
            fibers_state(fiber_i) = (fibers(idx_comp(1),5) == 1);
            % can also add damaged
            if fibers_state(fiber_i) == 0 && strcmp(settings.butenko_intersectStatus,'activated')
                fibers_state(fiber_i) = (fibers(idx_comp(1),5) == -1 || fibers(idx_comp(1),5) == -3);
            end

            % The key line. Probability is estimated as the number of
            % activations across scalings divided by the number of the scalings
            idx_comp_orig = find(fibers_prob(:,4)==fiber_i);
            fibers_prob(idx_comp_orig,5) = fibers_prob(idx_comp_orig,5) + fibers_state(fiber_i) / N_samples;
        end
    end
    if skip_pathway
        % no axons survived for any sample -> nothing to write
        continue;
    end
    ftr.fibers = fibers_prob;
    ftr.ea_fibformat = ea_fibformat;
    ftr.connectome_name = connectome_name;
    Axon_state_file_prob = fullfile(results_folder, ['Axon_state_',pathways{pt_counter},'_prob.mat']);
    save(Axon_state_file_prob,'-struct','ftr')
end
