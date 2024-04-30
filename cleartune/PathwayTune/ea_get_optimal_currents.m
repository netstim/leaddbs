function contact_perc_currents = ea_get_optimal_currents(N_contacts,stimfolder,side,plot_progress)

% check pam_optimizer results and return currents for the best solution

if side == 0
    side_suffix = '_rh';
else
    side_suffix = '_lh';
end

if isfile([stimfolder,filesep,'NB',side_suffix,filesep,'optim_iterations.csv'])
    optim_result = readtable([stimfolder,filesep,'NB',side_suffix,filesep,'optim_iterations.csv']);
else
    ea_warndlg("No optimization result found!")
    return
end


% there are now entries with critical side-effects stored in optim_iterations.csv
% maximum because we maximize improvement
[~,max_ind] = max(optim_result.weighted_total_score);
contact_currents = [];
for contact_i = 0:N_contacts-1
    % convert to mA
    contact_currents = [contact_currents,optim_result.(['Contact_',num2str(contact_i)])(max_ind)*1000.0];

    % just for now because we do monopolar
    contact_currents(2:end) = 0.0; 
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
    goal_function = [];
    iter_num = [];

    for i = 1:size(optim_result,1)
        if i == 1
            goal_function = [goal_function, optim_result.weighted_total_score(i)];
            iter_num = [iter_num, i];
        else
            if optim_result.weighted_total_score(i) > goal_function(end)
               goal_function = [goal_function, optim_result.weighted_total_score(i)];
               iter_num = [iter_num, i]; 
            end
        end
    end

    figure
    plot(iter_num,goal_function,'b--o');
    xlabel('Iteration Number')
    ylabel('Weighted General Improvement')
    saveas(gcf,[stimfolder,filesep,'NB',side_suffix,filesep,'Optimization_convergence.png'])
end

    

