function ea_compute_pathways_recruitment(connectome_name,stim_prot_file,pat_folder,stim_folder,side,VAT_thresh)

%% Compute pathway recruitment based on VAT
% The results are used to train ANN
% Iterate over Current_protocols_' + str(side) + '.csv
% Save results in Activations_over_StimSets
% Currently, OSS-DBS file structure is used

% stim_prot_file - array (N_protocols X N_contacts)
% side           - 0-rh (OSS-DBS notation)

% load csv as array
stim_protocols = readtable(stim_prot_file);
stim_protocols = table2array(stim_protocols);

% output folders (OSS-DBS structure)
if side == 0
    NB_folder = 'NB_0';
    res_folder = 'Results_rh';
else
    NB_folder = 'NB_1';
    res_folder = 'Results_lh';
end

% get all pathways from the connectome
pathways_files = dir([connectome_name,filesep,'*.mat']);
if side == 0
    pathways_files = pathways_files(contains({pathways_files.name}, '_rh_'));
else
    pathways_files = pathways_files(contains({pathways_files.name}, '_lh_'));
end

% will be stored in Results / Activations_over_StimSets
N_recruited = zeros(size(stim_protocols,1),length(pathways_files) + 1);
save_pathways_to_json = 1;

for stim_i = 1:size(stim_protocols,1)

    N_recruited(stim_i,1) = stim_i - 1;  % store iteration (for consistency with OSS-DBS)

    % StimVector = stim_protocols(stim_i,:)
    % convert mA to perc
    % stay with cathodic stim only for now (pol = 1 for all)
    % simply sum up all currents, recompute percentages
    % S = ea_add_StimVector_to_S(S, StimVector,side)

    % just because all cathodic
    if any(stim_protocols(:) > 0.0)
        disp("Only cathodic stimulation is supported!")
        return
    else
        ampselect = sum(abs(stim_protocols(stim_i,:)));
    end
    perc_val = zeros(1,8);
    % keep original sign for perc here
    perc_val(1:size(stim_protocols,2)) = 100.0 * stim_protocols(stim_i,:)./ampselect;
    constcurr = 0;  % 0 - CC, 1 - VC (Lead-DBS notation)
    writeVTA = 0;
    modelVTA = 'FieldTrip';

    % Nanditha's script
    try
        [Efield,~] = ea_generate_vat(pat_folder,ampselect,perc_val,constcurr,side+1,writeVTA,modelVTA);
    catch
        disp("Failed to compute Efield for this stim protocol, it will be marked")
        N_recruited(stim_i,:) = -1.0;
        continue
    end

    % Threshold the vat efield
    vatInd = find(abs(Efield.img(:)) > VAT_thresh);

    % Trim connectome fibers
    [xvox, yvox, zvox] = ind2sub(size(Efield.img), vatInd);
    vatmm = ea_vox2mm([xvox, yvox, zvox], Efield.mat);
    
    for path_i = 1:length(pathways_files)

        pathway_file = fullfile(pathways_files(path_i).folder, pathways_files(path_i).name);
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
        N_recruited(stim_i,path_i+1) = sum(connected);


        if stim_i == 1 || save_pathways_to_json == 1
            % save pathway name and number of fibers
            % only needed once
            % also save index to resort json later as Activations_over_StimSets
            jsonDict.vat_paths_dict.(genvarname(pathways_files(path_i).name(1:end-4))) = [length(idx),path_i];
            save_pathways_to_json = 0;
        end
    end
end

% save simulated pathway info
jsonText = jsonencode(jsonDict);
fid = fopen([stim_folder,filesep,NB_folder,filesep,'VAT_pathways.json'], 'w');
fprintf(fid, '%s', jsonText)
fclose(fid);

if ~exist([stim_folder,filesep,res_folder], 'dir')
   mkdir([stim_folder,filesep,res_folder])
end

% save to .csv
if side == 0
    csvwrite([stim_folder,filesep,res_folder,filesep,'Activations_over_StimSets_rh.csv'],N_recruited)
else
    csvwrite([stim_folder,filesep,res_folder,filesep,'Activations_over_StimSets_lh.csv'],N_recruited)
end

