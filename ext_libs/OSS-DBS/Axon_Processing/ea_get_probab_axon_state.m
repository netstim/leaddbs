function ea_get_probab_axon_state(results_folder,plot_rates,damaged_as_activated)
% Estimate probability of axon activation based on a sweep of parameters
% (e.g. fiber diameters).
% By Butenko and Li, konstantinmgtu@gmail.com

arguments
    results_folder          % path to save the output
    plot_rates              {mustBeNumericOrLogical} = 0  % if true, plot activation rates over parameter sweep
    damaged_as_activated    {mustBeNumericOrLogical} = 0  % if true, interpret axons intersected with the electrode (damaged) as activated
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

N_scalings = length(activations_over_pathways{pt_counter});

% plot activation curves over the parameter sweep
if plot_rates
    figure()
    for pt_counter = 1:length(pathways)
        pathway_name = strrep(pathways{pt_counter},'_',' ');
        disp(pathway_name)
        plot(1:N_scalings,activations_over_pathways{pt_counter},'DisplayName',pathway_name)
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
            if damaged_as_activated
                fibers_state(fiber_i) = (fibers(idx_comp(1),5) == -1);
                fibers_state(fiber_i) = (fibers(idx_comp(1),5) == -3);
            end

            % key line. Probability is estimated as number of success
            % across scaling divided by the number of scalings
            idx_comp_orig = find(fibers_prob(:,4)==fiber_i);
            fibers_prob(idx_comp_orig,5) = fibers_prob(idx_comp_orig,5) + fibers_state(fiber_i) / N_scalings;
        end
    end
    ftr.fibers = fibers_prob;
    ftr.ea_fibformat = ea_fibformat;
    ftr.connectome_name = connectome_name;
    Axon_state_file_prob = fullfile(results_folder, ['Axon_state_',pathways{pt_counter},'_prob.mat']);
    save(Axon_state_file_prob,'-struct','ftr')
end
