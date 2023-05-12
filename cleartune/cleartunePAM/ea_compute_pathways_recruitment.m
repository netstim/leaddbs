function ea_compute_pathways_recruitment(connectome_name,stim_prot_file,stim_folder,side,VAT_thresh)

%% Compute pathway recruitment based on VAT
% The results are used to train ANN
% Iterate over Current_protocols_' + str(side) + '.csv
% Save results in Activations_over_StimSets
% Currently, OSS-DBS file structure is used

% stim_prot_file - array (N_protocols X N_contacts)
% side           - 0-rh (OSS-DBS notation)

if side == 0
    NB_folder = 'NB_0';
    res_folder = 'Results_rh';
else
    NB_folder = 'NB_1';
    res_folder = 'Results_lh';
end


pathways_files = dir([connectome_name,filesep,'*.mat']);
if side == 0
    pathways_files = pathways_files(contains({pathways_files.name}, '_rh_'));
else
    pathways_files = pathways_files(contains({pathways_files.name}, '_lh_'));
end

% will be stored in Activations_over_StimSets
percent_recruit = zeros(size(stim_prot_file,1),length(pathways_files) + 1);

%stim_protocols = load(stim_prot_file);

for stim_i = 1:size(stim_prot_file,1)

    percent_recruit(stim_i,1) = stim_i;  % store iteration (for consistency with OSS-DBS)

        
    % stim_vector = stim_protocols(stim_i,:)
    % convert mA to perc
    % stay with cathodic stim only for now
    % simply sum up all currents, recompute percentages
    % Efield = ea_get_magic_Efield(stim_vector_perc, options, reconst)

    % Threshold the vat efield
    vatInd = find(abs(Efield.img(:)) > VAT_thresh);

    % Trim connectome fibers
    [xvox, yvox, zvox] = ind2sub(size(Efield.img), vatInd);
    vatmm = ea_vox2mm([xvox, yvox, zvox], Efield.mat);
    
    for path_i = 1:length(pathways_files)

        pathway_file = fullfile(pathways_files(k).folder, pathways_files(k).name);
        load(pathway_file, 'fibers', 'idx');

        filter = all(fibers(:,1:3)>=min(vatmm),2) & all(fibers(:,1:3)<=max(vatmm), 2);

        % Skip further calculation in case VAT is totally not connected
        if ~any(filter)
            continue;
        end

        trimmedFiber = fibers(filter,:);

        % Map mm connectome fibers into VAT voxel space
        [trimmedFiberInd, ~, trimmedFiberID] = unique(trimmedFiber(:,4), 'stable');
        fibVoxInd = splitapply(@(fib) {ea_mm2uniqueVoxInd(fib, Efield)}, trimmedFiber(:,1:3), trimmedFiberID);

        % Remove outliers
        fibVoxInd(cellfun(@(x) any(isnan(x)), fibVoxInd)) = [];
        trimmedFiberInd(cellfun(@(x) any(isnan(x)), fibVoxInd)) = [];

        % Find connected fibers
        connected = cellfun(@(fib) any(ismember(fib, vatInd)), fibVoxInd);
        
        % check!!!
        percent_recruit(stim_i,path_i+1) = sum(connected);


        % save pathway name and number of fibers
        % only needed once
        if stim_i == 1
            jsonDict.vat_paths_dict.(genvarname(pathways_files(k).name)) = length(idx);
        end
    end
end

% save simulated pathway info
jsonText = jsonencode(jsonDict);
fid = fopen([stim_folder,filesep,NB_folder,filesep,'VAT_pathways.json'], 'w');
fprintf(fid, '%s', jsonText)
fclose(fid);

% save to .csv
if side == 0
    csvwrite([stim_folder,res_folder,filesep,'Activations_over_StimSets_rh.csv'],percent_recruit)
else
    csvwrite([stim_folder,res_folder,filesep,'Activations_over_StimSets_lh.csv'],percent_recruit)
end

