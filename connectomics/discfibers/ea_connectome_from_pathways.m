function ea_connectome_from_pathways(obj, connectomeName, N_voters, N_sides, negative_vals)

    % check how many are duplicates
    % store in one folder in connectomes/dMRI_MultiTract
    % give original names to connectomes, the function will overwrite
    % files!
    
    connectome_folder = [ea_getconnectomebase('dMRI_multitract'), connectomeName];
    % add pathways to the connectome if exists (Symptoms + side-effects)
    if ~isfolder(connectome_folder)
        mkdir(connectome_folder)
        %rmdir(connectome_folder, 's')
    end

    all_indices = []; % potentially slow
    disp("Copying pathways to dMRI_MultiTract/...")
    for voter = 1:N_voters
        for side = 1:N_sides
    
            if side == 1
                prefix_side = '_rh';
            else
                prefix_side = '_lh';
            end

            % think about some standard way of symptom naming
            if N_voters > 1
                % think about a function here for multitracts,
                symptomName = [obj.subscore.labels{voter,1},prefix_side];
            else
                symptomName = [obj.responsevarlabel,prefix_side];
            end
        
            % change - to _
            symptomName = strrep(symptomName,'-','_');
    
            % prepare output folders for pathways
            [filepath,~,~] = fileparts(obj.leadgroup);
            symptomNameFolder = [filepath,filesep,symptomName(1:end-3)];
    
            % store soft side-effect pathways separately
            if negative_vals == 1
                symptomNameFolder = [symptomNameFolder,'_SE'];
            end
    
            pathways_files = dir([symptomNameFolder,filesep,'*.mat']);
            if side == 1
                pathways_files = pathways_files(contains({pathways_files.name}, '_rh_'));
            else
                pathways_files = pathways_files(contains({pathways_files.name}, '_lh_'));
            end

            for k = 1:length(pathways_files)
                pathway_file = fullfile(pathways_files(k).folder, pathways_files(k).name);
                % we only need the vector of global indices (defined for each point here)
                load(pathway_file, 'glob_ind');
                
                gl_indices = unique(glob_ind)';
                all_indices = [all_indices, gl_indices];
                % fibers are unique within the symptom, so we compare only with
                % other symptoms at the end
    
    
                pathway_file_conn = fullfile(connectome_folder, pathways_files(k).name);
                copyfile(pathway_file, pathway_file_conn)
            end


        end
    end

    disp("Checking number of repetitions")
    N_rep = 0;
    indices = unique(all_indices);
    for i = 1:length(indices)
        N_rep = N_rep + sum(all_indices == indices(i)) - 1;
    end
    disp(N_rep)

end
