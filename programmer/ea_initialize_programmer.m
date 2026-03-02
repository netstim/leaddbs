function ea_initialize_programmer(handles, bids, action)
    %% Set the lead DBS patient folder path here
    groupFolderPath = bids.datasetDir; 
    %% Initialize Parameters
    groupOptions = BIDSFetcher(groupFolderPath);
    
    %% Prepare PLY Reconstructions
    preparePLYReconstructions(groupFolderPath, groupOptions, handles);
    
    %% Prepare Reconstruction JSON Files
    prepareReconstructionJSON(groupFolderPath, groupOptions);
    
    %% Prepare Atlases for Export
    % prepareAtlasesForExport();

    %% Prepare participants json file
     writeParticipantsJson(bids);

    %% Conditionally launch application
    if (action == "standalone")
        prepareStimulations(groupFolderPath, groupOptions);

        % prepare_ossdbs(groupOptions, handles);
        system([ea_getProgrammerPath, ' ', groupFolderPath, ' > /dev/null 2>&1 &']);
    end
end


function prepare_ossdbs(groupOptions, handles)
    ptID = groupOptions.subjId{1};
    BidsID = strcat('sub-', ptID);
    groupFolderPath = groupOptions.datasetDir;
    patientFolderPath = fullfile(groupFolderPath, 'derivatives', 'leaddbs', BidsID);
    options = ea_getptopts(patientFolderPath);
    S = ea_initializeS('initialize', options);
    options.stimSetMode = 0;
    S.Ls1.k1.perc = 100;
    S.Ls1.k1.pol = 1;
    S.Ls1.amp = 2;
    S.Rs1.k1.perc = 100;
    S.Rs1.k1.pol = 1;
    S.Rs1.amp = 1;
    S.activecontacts{1} = [1,0,0,0,0,0,0,0];
    S.activecontacts{2} = S.activecontacts{1};
    S.amplitude{1} = [2,0,0,0];
    S.amplitude{2} = S.amplitude{1};
    options.optimizer = 1;
    ea_genvat_butenko(S, options);
end


%% Helper Functions

function prepareStimulations(groupFolderPath, groupOptions)
    for i = 1:numel(groupOptions.subjId)
        ptID = groupOptions.subjId{i};
        bidsPt = strcat('sub-', ptID);
        ptFilePath = fullfile(groupFolderPath, 'derivatives', 'leaddbs', bidsPt);
        stimDir = fullfile(ptFilePath, 'stimulations');
        clinicalDir = fullfile(ptFilePath, 'clinical');
        stimFileName = [bidsPt, '_desc-stimparameters'];
        stimMatFile = ea_regexpdir(stimDir, ['^', stimFileName, '\.mat$'], 1, 'f');
        [~, stimFolder] = fileparts(fileparts(stimMatFile));
        if ~iscell(stimFolder)
            stimFolder = {stimFolder};
        end
        % stimMatFile = stimMatFile(~startsWith(stimFolder, 'gs_'));
        % stimFolder = stimFolder(~startsWith(stimFolder, 'gs_'));
        for k=1:numel(stimMatFile)
            S = ea_checkStimParams(stimMatFile{k});
            sessionId = strcat('ses-', stimFolder{k});
            filename = strcat(bidsPt, '_', sessionId, '_stimparameters.json');
            stimJsonFile = fullfile(clinicalDir, sessionId, filename);
            if ~exist(stimJsonFile)
                ea_mkdir(fullfile(clinicalDir, sessionId));
                savejson('', S, stimJsonFile);
                disp(['Created ', stimJsonFile])
            end
        end
    end
end


function preparePLYReconstructions(groupFolderPath, groupOptions, handles)
    for i = 1:numel(groupOptions.subjId)
        ptID = groupOptions.subjId{i};
        bidsPt = strcat('sub-', ptID);
        ptFilePath = fullfile(groupFolderPath, 'derivatives', 'leaddbs', bidsPt);
        plyDir = fullfile(ptFilePath, 'export', 'ply');
        
        if ~isfolder(plyDir)
            try 
                ea_pat2ply(ptFilePath, handles);
                disp(['Processed PLY reconstructions for: ', bidsPt]);
            catch
                disp('failed');
            end
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
    atlasDir = fullfile(ea_space, 'atlases');
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
        coords1 = reco.mni.coords_mm{1,1};
        coords2 = reco.mni.coords_mm{1,2};
        % Prepare data for JSON encoding
        jsonData = struct( ...
            'elmodel', elmodel, ...
            'coords_right', coords1, ...
            'coords_left', coords2, ...
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
        participantsData = loadjson(participantsPath);
        if isempty(participantsData)
            participantsData = {}; % Ensure it remains a cell array
        elseif ~iscell(participantsData)
            participantsData = {participantsData}; % Convert struct to cell array
        end
    else
        participantsData = {}; % No file exists, start fresh
    end
    % Extract existing IDs from the cell array, ensuring it's not empty
    if ~isempty(participantsData)
%         existingIds = cellfun(@(x) x.id, participantsData, 'UniformOutput', false);
%         idIndexMap = containers.Map(existingIds, 1:length(existingIds));
        validParticipants = participantsData(~cellfun(@(x) ~isfield(x, 'id') || isempty(x.id), participantsData));
        
        % Extract IDs from the valid participants
        existingIds = cellfun(@(x) x.id, validParticipants, 'UniformOutput', false);
        
        % Create the map using only valid IDs
        idIndexMap = containers.Map(existingIds, 1:length(existingIds));

    else
        existingIds = {}; % Prevent empty keys
        idIndexMap = containers.Map('KeyType', 'char', 'ValueType', 'double'); % Empty map
    end
    % Update or add patient entries
    for i = 1:numPatients
        % disp(i);
        ptID = groupOptions.subjId{i};
        BidsID = strcat('sub-', ptID);
        patientFolderPath = fullfile(groupFolderPath, 'derivatives', 'leaddbs', BidsID);
        % Get patient options
        options = ea_getptopts(patientFolderPath);
        % Ensure elmodel exists
        if ~isfield(options, 'elmodel') || isempty(options.elmodel)
            options.elmodel = 'Medtronic 3387';
        end
        % Check if the patient ID already exists in the JSON file
        if isKey(idIndexMap, BidsID)
            % Update the existing entry
            participantsData{idIndexMap(BidsID)}.elmodel = options.elmodel;
        else
            % Add a new entry for this patient
            newEntry = struct();
            newEntry.id = BidsID;
            newEntry.elmodel = options.elmodel;
            participantsData{end+1} = newEntry; % Append to cell array
        end
    end
    % Save the updated data back to the JSON file
    savejson('', participantsData, participantsPath);
end
