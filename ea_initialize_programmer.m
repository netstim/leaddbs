function ea_initialize_programmer(handles, bids, action)
    %% Set the lead DBS patient folder path here
    groupFolderPath = bids.datasetDir; 
    %% Initialize Parameters
    % handles = struct('atlasset', 'DISTAL Minimal (Ewert 2017)');
    groupOptions = BIDSFetcher(groupFolderPath);
    numPatients = numel(groupOptions.subjId);
    
    %% Prepare PLY Reconstructions
    preparePLYReconstructions(groupFolderPath, groupOptions, handles);
    
    %% Prepare Reconstruction JSON Files
    prepareReconstructionJSON(groupFolderPath, groupOptions);
    
    %% Prepare Atlases for Export
    prepareAtlasesForExport();
    
    %% Prepare Reconstructions and Orientations for Export
    % masterData = struct();
    % for i = 1:numPatients
    %     ptID = groupOptions.subjId{i};
    %     BidsID = strcat('sub-', ptID);
    %     patientFolderPath = fullfile(groupFolderPath, 'derivatives', 'leaddbs', BidsID);
    % 
    %     % Ensure BidsID is a valid field name
    %     validBidsID = matlab.lang.makeValidName(BidsID);
    % 
    %     % Get patient options and clinical directories
    %     options = ea_getptopts(patientFolderPath);
    %     clinicalDirectory = options.subj.clinicalDir;
    %     exportDirectory = options.subj.exportDir;
    % 
    %     % Initialize master data structure for this patient
    %     masterData.(validBidsID) = struct();
    %     masterData.(validBidsID).clinicalData = gatherClinicalData(clinicalDirectory, BidsID);
    % 
    %     % Process and store PLY paths
    %     masterData.(validBidsID).exportData = gatherPLYData(exportDirectory);
    % end
    % 
    % % Save the final master dataset as a JSON file
    % datasetMasterPath = fullfile(groupFolderPath, 'dataset_master.json');
    % savejson('', masterData, datasetMasterPath);

    %% Conditionally launch application
    if (action == "standalone")
        writeParticipantsJson(bids);
        currentOS = ea_getarch;
        releaseDir = fullfile(ea_getearoot, 'programmer', 'app', 'release');
        zipFile = fullfile(releaseDir, ['LeadDBSProgrammer_', currentOS, '.zip']);
        if ismac
            appFile = fullfile(ea_prefsdir, 'Programmer', 'LeadDBSProgrammer.app', 'Contents', 'MacOS', 'LeadDBSProgrammer');
            if ~isfile(appFile)
                unzip(zipFile, fullfile(ea_prefsdir, 'Programmer'));
                system(['xattr -cr ', ea_path_helper(fullfile(ea_prefsdir, 'Programmer', 'LeadDBSProgrammer.app'))]);
                savejson('', struct('LeadDBS_Path', ea_getearoot), fullfile(ea_prefsdir, 'Programmer', 'Preferences.json'));
            end
        elseif isunix
            appFile = fullfile(ea_prefsdir, 'Programmer', 'LeadDBSProgrammer', 'LeadDBSProgrammer');
            if ~isfile(appFile)
                unzip(zipFile, fullfile(ea_prefsdir, 'Programmer', 'LeadDBSProgrammer'));
                savejson('', struct('LeadDBS_Path', ea_getearoot), fullfile(ea_prefsdir, 'Programmer', 'Preferences.json'));
            end
        else
            appFile = fullfile(ea_prefsdir, 'Programmer', 'LeadDBSProgrammer', 'LeadDBSProgrammer.exe');
            if ~isfile(appFile)
                unzip(zipFile, fullfile(ea_prefsdir, 'Programmer', 'LeadDBSProgrammer'));
                savejson('', struct('LeadDBS_Path', ea_getearoot), fullfile(ea_prefsdir, 'Programmer', 'Preferences.json'));
            end
        end
        system([appFile, ' ', groupFolderPath, ' > /dev/null 2>&1 &']);
    end
end


%% Helper Functions

function preparePLYReconstructions(groupFolderPath, groupOptions, handles)
    for i = 1:numel(groupOptions.subjId)
        ptID = groupOptions.subjId{i};
        bidsPt = strcat('sub-', ptID);
        ptFilePath = fullfile(groupFolderPath, 'derivatives', 'leaddbs', bidsPt);
        plyDir = fullfile(ptFilePath, 'export', 'ply');
        
        if ~isfolder(plyDir)
            % Folder does not exist: create PLY reconstructions
            ea_pat2ply(ptFilePath, handles);
            disp(['Processed PLY reconstructions for: ', bidsPt]);
        else
            % Folder exists: check whether it is empty
            dirContents = dir(plyDir);
            % Remove '.' and '..' entries
            dirContents = dirContents(~ismember({dirContents.name}, {'.', '..'}));
            
            if isempty(dirContents)
                % Folder is empty, so run the processing
                ea_pat2ply(ptFilePath, handles);
                disp(['Processed PLY reconstructions for: ', bidsPt]);
            % else
            %     disp(['PLY reconstructions already exist for: ', bidsPt]);
            end
        end
    end
end


function prepareReconstructionJSON(groupFolderPath, groupOptions)
    for i = 1:numel(groupOptions.subjId)
        ptID = groupOptions.subjId{i};
        bidsPt = strcat('sub-', ptID);
        clinicalDirectory = fullfile(groupFolderPath, 'derivatives', 'leaddbs', bidsPt, 'clinical');
        recoOutputFile = fullfile(clinicalDirectory, strcat(bidsPt, '_desc-reconstruction.json'));
        
        if ~isfile(recoOutputFile)
            recoFileName = strcat(bidsPt, '_desc-reconstruction.mat');
            ptFilePath = fullfile(groupFolderPath, 'derivatives', 'leaddbs', bidsPt, 'reconstruction', recoFileName);
            
            if isfile(ptFilePath)
                reco_data = writeRecoToJSON(ptFilePath);
                if ~exist(clinicalDirectory, 'dir')
                    mkdir(clinicalDirectory);
                end
                savejson('', reco_data, recoOutputFile);
                disp(['Saved reconstruction JSON for: ', bidsPt]);
            else
                warning(['Reconstruction file not found for: ', bidsPt]);
            end
        else
            disp(['Reconstruction JSON already exists for: ', bidsPt]);
        end
    end
end


function prepareAtlasesForExport()
    leadPath = ea_getearoot();
    atlasDir = fullfile(leadPath, 'templates/space/MNI_ICBM_2009b_NLIN_ASYM/atlases');
    atlasList = dir(atlasDir);

    for i = 1:numel(atlasList)
        atlasName = atlasList(i).name;
        if any(strcmp(atlasName, {'.', '..', '.DS_Store'})), continue, end
    
        atlasPath = fullfile(atlasDir, atlasName);
        outputDir = fullfile(atlasDir, atlasName, strcat(atlasName, '.ply'));
        
        if ~isfile(outputDir)
            try
                ea_atlas2ply({atlasName}, outputDir);
                disp(['Successfully processed atlas: ', atlasName]);
            catch ME
                disp(['Failed to process atlas: ', atlasName]);
                disp(['Error: ', ME.message]);
            end
        % else
        %     disp(['Atlas already exported: ', atlasName]);
        end
    end
end


% function clinicalData = gatherClinicalData(clinicalDirectory, BidsID)
%     clinicalData = struct();
%     sessionDirs = dir(clinicalDirectory);
%     sessionDirs = sessionDirs([sessionDirs.isdir]);
%     sessionDirs = sessionDirs(~ismember({sessionDirs.name}, {'.', '..'}));
%     
%     for j = 1:numel(sessionDirs)
%         sessionName = sessionDirs(j).name;
%         sessionFieldName = matlab.lang.makeValidName(sessionName);
%         sessionFolderPath = fullfile(clinicalDirectory, sessionName);
%     
%         % List JSON files in the session directory
%         fileList = dir(fullfile(sessionFolderPath, '*.json'));
%         for k = 1:numel(fileList)
%             filePath = fullfile(sessionFolderPath, fileList(k).name);
%             clinicalData.(sessionFieldName){k} = filePath;
%         end
%     
%         % Add reconstruction file if available
%         reconstructionFilePattern = fullfile(clinicalDirectory, sprintf('%s_desc-reconstruction.json', BidsID));
%         if isfile(reconstructionFilePattern)
%             clinicalData.reconstructionJson = reconstructionFilePattern;
%         end
%     end
% end


% function exportData = gatherPLYData(exportDirectory)
%     exportData = struct();
%     plyDir = fullfile(exportDirectory, 'ply');
%     
%     if isfolder(plyDir)
%         anatomyFilePath = fullfile(plyDir, 'anatomy.ply');
%         combinedElectrodesFilePath = fullfile(plyDir, 'combined_electrodes.ply');
%     
%         if isfile(anatomyFilePath)
%             exportData.anatomyPly = anatomyFilePath;
%         else
%             warning('anatomy.ply not found');
%         end
%     
%         if isfile(combinedElectrodesFilePath)
%             exportData.combinedElectrodesPly = combinedElectrodesFilePath;
%         else
%             warning('combined_electrodes.ply not found');
%         end
%     else
%         warning('PLY directory not found');
%     end
% end


function [reco_output] = writeRecoToJSON(matFileName)
    % Load the .mat file
    data = load(matFileName);
    
    % Extract the necessary fields from reco
    reco = data.reco;
    
    % Extract the electrode model and markers
    elmodel = reco.props(1).elmodel;
    
    try
        % Extract marker data
        head1 = reco.mni.markers(1).head;
        head2 = reco.mni.markers(2).head;
        tail1 = reco.mni.markers(1).tail;
        tail2 = reco.mni.markers(2).tail;
        
        % Calculate the directionality (unit vector and roll out) for the first marker
        unitvector1 = reco.native.markers(1).y - reco.native.markers(1).head;
        unitvector1(3) = 0;  % Set the z-component to zero
        unitvector1 = unitvector1 / norm(unitvector1);  % Normalize the vector
        roll_out_1 = rad2deg(atan2(norm(cross([0 1 0], unitvector1)), dot([0 1 0], unitvector1)));
    
        % Calculate the directionality (unit vector and roll out) for the second marker
        unitvector2 = reco.native.markers(2).y - reco.native.markers(2).head;
        unitvector2(3) = 0;  % Set the z-component to zero
        unitvector2 = unitvector2 / norm(unitvector2);  % Normalize the vector
        roll_out_2 = rad2deg(atan2(norm(cross([0 1 0], unitvector2)), dot([0 1 0], unitvector2)));
        
        % Prepare data for JSON encoding
        jsonData = struct( ...
            'elmodel', elmodel, ...
            'markers', struct( ...
                'head1', head1, ...
                'head2', head2, ...
                'tail1', tail1, ...
                'tail2', tail2 ...
            ), ...
            'directionality', struct( ...
                'unitvector1', unitvector1, ...
                'roll_out_right', roll_out_1, ...
                'unitvector2', unitvector2, ...
                'roll_out_left', roll_out_2 ...
            ) ...
        );
    catch ME
        % If an error occurs, set jsonData to an empty structure
        jsonData = struct();
        % Optionally, log or display the error message
        fprintf('An error occurred: %s\n', ME.message);
    end

    reco_output = jsonData;
end


function writeParticipantsJson(bids)
    % Path to the participants.json file
    groupFolderPath = bids.datasetDir; 
    groupOptions = BIDSFetcher(groupFolderPath);
    numPatients = numel(groupOptions.subjId);
    participantsPath = fullfile(groupFolderPath, 'participants.json');
    
    % Initialize participantsData
    if isfile(participantsPath)
        % If the file exists, read the existing data
        participantsData = loadjson(participantsPath);
    else
        % Initialize an empty array if the file doesn't exist
        participantsData = [];
    end
    
    % Extract existing IDs from the cell array
    if iscell(participantsData)
        existingIds = cellfun(@(x) x.id, participantsData, 'UniformOutput', false);
    else
        existingIds = {};
    end
    idIndexMap = containers.Map(existingIds, 1:length(existingIds));
    
    % Update or add patient entries
    for i = 1:numPatients
        ptID = groupOptions.subjId{i};
        BidsID = strcat('sub-', ptID);
        patientFolderPath = fullfile(groupFolderPath, 'derivatives', 'leaddbs', BidsID);
        
        % Get patient options
        options = ea_getptopts(patientFolderPath);
        
        % Check if the patient ID already exists in the JSON file
        if isKey(idIndexMap, BidsID)
            % Update the existing entry
            participantsData{idIndexMap(BidsID)}.elmodel = options.elmodel;
        else
            % Add a new entry for this patient
            newEntry = struct();
            newEntry.id = BidsID;
            newEntry.elmodel = options.elmodel;
            participantsData = [participantsData; newEntry]; %#ok<AGROW>
        end
    end
    
    % Remove entries that are not in the current groupOptions
    % validIds = strcat('sub-', groupOptions.subjId);
    % participantsData = participantsData(ismember({participantsData.id}, validIds));
    % participantsPath = '/Users/savirmadan/Downloads/participants.json';
    % Save the updated data back to the JSON file
    savejson('', participantsData, participantsPath);
end
