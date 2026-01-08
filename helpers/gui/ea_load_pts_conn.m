function ea_load_pts_conn(handles, uipatdir)
% Simplified version for Lead-Connectome: load patient and check BIDS status

if ~iscell(uipatdir)
    uipatdir = {uipatdir};
end

%% BIDS Detection & Mapper Integration
for i = 1:length(uipatdir)
    patdir = uipatdir{i};
    
    % Detect if this is BIDS or non-BIDS data
    isBIDS = ea_detect_bids(patdir);
    
    if ~isBIDS
        fprintf('\n=== Non-BIDS data detected ===\n');
        fprintf('Patient directory: %s\n', patdir);
        fprintf('Opening BIDS mapper to convert to BIDS structure...\n\n');
        
        % Open BIDS mapper
        try
            % Determine dataset root and subject ID
            BIDSRoot = ea_getdataset;
            [~, subjDirName] = fileparts(patdir);
            subjId = validateSubjId(regexprep(subjDirName, '^sub-', ''));

            % For non-BIDS data, gather all NIfTI files from the patient directory
            niiFiles = dir(fullfile(patdir, '*.nii'));
            niiGzFiles = dir(fullfile(patdir, '*.nii.gz'));
            allNiiFiles = [niiFiles; niiGzFiles];
            
            if ~isempty(allNiiFiles)
                % Build full paths
                niiFilePaths = cell(length(allNiiFiles), 1);
                for j = 1:length(allNiiFiles)
                    niiFilePaths{j} = fullfile(allNiiFiles(j).folder, allNiiFiles(j).name);
                end
                
                % Call ea_nifti_to_bids directly with the files
                fprintf('Opening BIDS mapping GUI with %d NIfTI files...\n', length(niiFilePaths));
                
                % Wrap in try-catch to handle GUI close gracefully
                try
                    ea_nifti_to_bids(niiFilePaths, BIDSRoot, ['sub-', subjId]);
                catch ME
                    % User might have closed the GUI - that's OK
                    fprintf('BIDS mapping GUI closed.\n');
                end
                
                % After mapping, files should be in rawdata/sub-XXX/
                % Lead-Connectome will copy them to derivatives/preprocessing/ as needed
            else
                warning('No NIfTI files found in: %s', patdir);
            end
            
            % Create derivatives folder structure
            derivDir = fullfile(BIDSRoot, 'derivatives', 'leaddbs', ['sub-', subjId]);
            ea_mkdir(fullfile(derivDir, 'preprocessing', 'dwi'));
            ea_mkdir(fullfile(derivDir, 'preprocessing', 'func'));
            ea_mkdir(fullfile(derivDir, 'connectomics'));
            
            % Create participants.tsv if it doesn't exist (required for BIDSFetcher)
            participantsFile = fullfile(BIDSRoot, 'participants.tsv');
            if ~exist(participantsFile, 'file')
                fprintf('Creating participants.tsv...\n');
                fid = fopen(participantsFile, 'w');
                fprintf(fid, 'participant_id\n');
                fprintf(fid, 'sub-%s\n', subjId);
                fclose(fid);
            else
                % Add subject to existing participants.tsv if not already there
                participants = readtable(participantsFile, 'FileType', 'text', 'Delimiter', '\t');
                subjIdFull = ['sub-', subjId];
                if ~any(strcmp(participants.participant_id, subjIdFull))
                    fprintf('Adding subject to participants.tsv...\n');
                    fid = fopen(participantsFile, 'a');
                    fprintf(fid, '%s\n', subjIdFull);
                    fclose(fid);
                end
            end
            
            if isfolder(derivDir)
                fprintf('BIDS structure ready at: %s\n', derivDir);
                uipatdir{i} = derivDir;
                
                % IMPORTANT: Update the GUI to point to the new BIDS location
                fprintf('\n*** Patient directory updated to BIDS location ***\n');
                fprintf('New location: %s\n\n', derivDir);
            else
                warning('Failed to create BIDS structure at: %s', newPatdir);
            end
        catch ME
            warning(ME.identifier, 'BIDS mapper failed: %s', ME.message);
            fprintf('Attempting to proceed with classic file structure...\n');
        end
    else
        % Determine dataset root and subject ID
        BIDSRoot = regexp(patdir, ['.*(?=\',filesep,'(derivatives|derivatives\',filesep,'leaddbs|rawdata|sourcedata))'], 'match', 'once');
        [~, subjDirName] = fileparts(patdir);
        subjId = regexprep(subjDirName, '^sub-', '');
        % Create derivatives folder structure
        derivDir = fullfile(BIDSRoot, 'derivatives', 'leaddbs', ['sub-', subjId]);
        ea_mkdir(fullfile(derivDir, 'preprocessing', 'dwi'));
        ea_mkdir(fullfile(derivDir, 'preprocessing', 'func'));
        ea_mkdir(fullfile(derivDir, 'connectomics'));
        fprintf('BIDS-compliant data detected: %s\n', derivDir);
        uipatdir{i} = derivDir;
    end
end

if exist('BIDSRoot', 'var')
    bids = BIDSFetcher(BIDSRoot);
    setappdata(handles.leadfigure, 'bids', bids);
    setappdata(handles.leadfigure, 'subjId', subjId);
end

% Store in appdata for use by other functions
setappdata(handles.leadfigure, 'uipatdir', uipatdir);

% Update the patdir choosebox
if isfield(handles, 'patdir_choosebox')
    if length(uipatdir) > 1
        handles.patdir_choosebox.String = ['Multiple (', num2str(length(uipatdir)), ')'];
        handles.patdir_choosebox.TooltipString = strjoin(uipatdir, '\n');
    else
        handles.patdir_choosebox.String = uipatdir{1};
        handles.patdir_choosebox.TooltipString = uipatdir{1};
    end
end

% Add to recent patients
if isfield(handles, 'leadfigure')
    ea_addrecent(handles, uipatdir, 'patients');
end

% Update guidata to ensure the new path is used
guidata(handles.leadfigure, handles);

fprintf('Subject loaded.\n\n');


function subjId = validateSubjId(subjId)
if ~isempty(regexp(subjId, '[\W_]', 'once'))
    subjId = regexprep(subjId, '[\W_]', '');
    ea_cprintf('CmdWinWarnings', 'It looks like you have special chars in your subj folder name.\nWe will use a cleaned name ''%s'' for the BIDS dataset. Please check manually.\n', subjId);
end 
