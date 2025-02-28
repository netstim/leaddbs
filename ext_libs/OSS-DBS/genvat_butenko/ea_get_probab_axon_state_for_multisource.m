function ea_get_probab_axon_state_for_multisource(results_folder,damaged_as_activated)
% Estimate probability of axon activation based on a sweep of parameters
% special case for multisource probabilistic PAM
% (e.g. fiber diameters).
% By Butenko and Li, konstantinmgtu@gmail.com

arguments
    results_folder          % path to save the output
    damaged_as_activated    {mustBeNumericOrLogical} = 0  % if true, interpret axons intersected with the electrode (damaged) as activated
end

% check which pathways were simulated and their percent activations
axon_state_files = dir([results_folder,filesep,'Axon_state_*']);
pathways = {};
activations_over_pathways = {};
pt_counter = 0;
% we enforce correct loading order
for file_i = 1:length(axon_state_files)
    parts = strsplit(axon_state_files(file_i).name,'_');
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
    scaling_counter = scaling_counter + 1;
    if ~any(strcmp(pathways,pathway_name_orig))
        pt_counter = pt_counter + 1;
        pathways{pt_counter} = pathway_name_orig;
        activations_over_pathways{pt_counter} = -42.0;  % dummy
    else
        activations_over_pathways{pt_counter} = [activations_over_pathways{pt_counter},-42.0];
    end
end

N_scalings = length(activations_over_pathways{pt_counter});

% now let's compute activation probabilty for each fiber
% we compute probabilities for the filtered fibers and will map to global
% indices later in ea_genvat_butenko
for pt_counter = 1:length(pathways)
    % iterate over scalings (fiber diameters)
    % important: number of compartments can differ based on the fiber
    % diameter / length!
    for scaling_i = 1:N_scalings

        Axon_state_file = fullfile(results_folder, ['Axon_state_',pathways{pt_counter},'_',num2str(scaling_i),'.mat']);
        load(Axon_state_file, 'fibers','ea_fibformat','connectome_name');

        % intialize new axon state file with probabilistic activations
        % morphologz defined by the first scaling!
        if scaling_i == 1
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
            if fibers_state(fiber_i) == 0 && damaged_as_activated
                fibers_state(fiber_i) = (fibers(idx_comp(1),5) == -1 || fibers(idx_comp(1),5) == -3);
            end

            % key line. Probability is estimated as number of success
            % across scaling divided by the number of scalings
            idx_comp_orig = find(fibers_prob(:,4)==fiber_i);
            fibers_prob(idx_comp_orig,5) = fibers_prob(idx_comp_orig,5) + fibers_state(fiber_i) / N_scalings;
        end
        ea_delete(Axon_state_file);
    end
    ftr.fibers = fibers_prob;
    ftr.ea_fibformat = ea_fibformat;
    ftr.connectome_name = connectome_name;
    Axon_state_file_prob = fullfile(results_folder, ['Axon_state_',pathways{pt_counter},'_prob.mat']);
    save(Axon_state_file_prob,'-struct','ftr')
end
