function contact_perc_currents = ea_get_optimal_currents(N_contacts,stimfolder,side,plot_progress)

% check pam_optimizer results and return currents for the best solution

if side == 0
    side_suffix = '_rh';
else
    side_suffix = '_lh';
end

if isfile([stimfolder,filesep,'NB',side_suffix,filesep,'optim_iterations_CSE.csv'])
    optim_result = readtable([stimfolder,filesep,'NB',side_suffix,filesep,'optim_iterations_CSE.csv']);
    % remove iterations with predicted critical side-effects
    CSE_field = ['CSE', side_suffix];
    toDelete = optim_result.(CSE_field) == 1;
    optim_result(toDelete,:) = [];

    if isempty(optim_result)
        ea_warndlg("All iterations predicted CSE!")
        contact_perc_currents = false;
        return
    end
elseif isfile([stimfolder,filesep,'NB',side_suffix,filesep,'optim_iterations.csv'])
    optim_result = readtable([stimfolder,filesep,'NB',side_suffix,filesep,'optim_iterations.csv']);
else
    ea_warndlg("No optimization result found!")
    contact_perc_currents = false;
    return
end


% there are now entries with critical side-effects stored in optim_iterations.csv
% maximum because we maximize improvement
[~,max_ind] = max(optim_result.weighted_total_score);
contact_currents = [];
for contact_i = 0:N_contacts-1
    % convert to mA
    contact_currents = [contact_currents,optim_result.(['Contact_',num2str(contact_i)])(max_ind)*1000.0];

    % % just for now because we do monopolar
    %contact_currents(2:end) = 0.0; 
end

if any(contact_currents > 0.0) && any(contact_currents < 0.0)
    % bipolar case
    cathode_currents = sum(contact_currents(contact_currents < 0.0)); 
    anode_currents = sum(contact_currents(contact_currents > 0.0)); 
    amplitude = max([abs(cathode_currents),anode_currents]);
else
    amplitude = sum(abs(contact_currents));
end

grounded_current = -1*sum(contact_currents);
% compute percentages in relation to the amplitude
contact_perc_currents = 100 * contact_currents / amplitude;
contact_perc_currents = [amplitude,contact_perc_currents];

if plot_progress
    % Construct the field names for target and SSE
    target_field = ['target', side_suffix];
    SSE_field = ['SSE', side_suffix];
    
    % Check if the necessary SSE field exists
    Exist_Column = strcmp(SSE_field,optim_result.Properties.VariableNames);
    show_target_and_SSE = Exist_Column(Exist_Column==1);
    
    % --- Data Extraction (Using Logical Indexing and Pre-allocation) ---
    
    % Assume optim_result is a structure array where each element 'i' is an iteration.
    % The logic seems to track iterations where the weighted_total_score improves (is greater than the last recorded score),
    % or is the very first iteration.
    
    % 1. Find the indices of iterations that show improvement or are the first one
    num_iterations = size(optim_result, 1);
    weighted_scores = [optim_result.weighted_total_score]; % Convert to a vector for easy comparison
    
    % Initialize a logical vector for the selected indices
    improvement_indices = false(1, num_iterations);
    if num_iterations >= 1
        improvement_indices(1) = true; % Always include the first iteration
        current_max_score = weighted_scores(1);
        
        % Find subsequent improvements
        for i = 2:num_iterations
            if weighted_scores(i) > current_max_score
                improvement_indices(i) = true;
                current_max_score = weighted_scores(i);
            end
        end
    end
    
    % 2. Extract the data for the selected iterations
    iter_num = find(improvement_indices);
    goal_function = weighted_scores(improvement_indices); % Optimization progress scores
    
    % Check if goal_function is empty (only happens if optim_result is empty)
    if isempty(goal_function)
        warning('optim_result appears to be empty or invalid.');
        contact_perc_currents = false;
        return; % Exit if no data to plot
    end
    
    if show_target_and_SSE
        % Extract target and SSE values directly using the logical indices and dynamic field names
        target_values = [optim_result.(target_field)(improvement_indices)];
        SSE_values = [optim_result.(SSE_field)(improvement_indices)];
    end
    
    % --- Plotting ---
    
    figure;
    hold on; % Use hold on to plot multiple lines on the same axes
    
    % Plot Optimization Progress (Goal Function)
    plot(iter_num, goal_function, 'b-o', 'DisplayName', 'Weighted Total Score');
    
    if show_target_and_SSE
        % Plot Target Values
        plot(iter_num, target_values, 'r--x', 'DisplayName', ['Target (', side_suffix(2:end), ')']);
        
        % Plot SSE Values
        plot(iter_num, SSE_values, 'g:s', 'DisplayName', ['SSE (', side_suffix(2:end), ')']);
    end
    
    xlabel('Iteration Number');
    ylabel('Value'); % Changed to 'Value' since it plots multiple types of values
    title('Optimization Progress');
    legend('show', 'Location', 'NorthWest'); % Show legend to distinguish the lines
    grid on;
    
    hold off;
    saveas(gcf,[stimfolder,filesep,'NB',side_suffix,filesep,'Optimization_convergence',side_suffix,'.png'])
end

    

