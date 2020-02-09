function options = ea_read_bids(options, bids_subject_folder)

%% get prefs
if ~exist('options', 'var') || isempty(options)
    options.prefs = ea_prefs('');
end

%% Choose the BIDS subject folder to be imported and read folder structure
if ~exist('bids_subject_folder', 'var')
    bids_subject_folder = uigetdir;
end

% remove filesep if required
if strcmp(bids_subject_folder(end), filesep)
    bids_subject_folder(end)=[];
end

% separate bids root folder from id folder and get id:
[bids_rawdata, id]=fileparts(bids_subject_folder);

% at this point lead checks that the subject folder begins with 'sub', not sure this is necessary:
if ~strcmp(id(1:3), 'sub')
    error('Not a valid BIDS subject folder')
end

%% Get BIDS session names - subfolders where preop and postop images are stored
% Note that users can specify the same folder for both, or even
% '' for one or both if no session-specific subfolder is used, by modifying
% their .ea_prefs file.
if ~isfield(options.prefs,'bids_session_preop') 
    options.prefs.bids_session_preop = 'ses-preDBS';
end
if ~isfield(options.prefs,'bids_session_postop') 
    options.prefs.bids_session_postop = 'ses-postDBS';
end
ses_preop = options.prefs.bids_session_preop;
ses_postop = options.prefs.bids_session_postop;


%% Create a lead derivatives folder in the BIDS root
% If subject's parent folder (bids_rawdata) is named 'rawdata',
% then the parent of _that_ is the bids dataset root, and we want our
% derivates folder to be a sister of rawdata. Otherwise we assume
% the subject's parent folder is in fact the dataset root.
[bids_dataset, rawdata] = fileparts(bids_rawdata);
if ~strcmpi(rawdata, 'rawdata')
    bids_dataset = bids_rawdata;
end
derivatives_folder = fullfile(bids_dataset, 'derivatives', 'leaddbs', id);
if ~exist(derivatives_folder, 'dir')
    mkdir(derivatives_folder)
end

%% Lookup table to map BIDS names to lead-dbs names.
% Note that the .nii.gz extension is assumed if not provided.
lead2bids_lookup = {
    % lead_dbs_sequence_name        BIDS_ses    BIDS_subfolder file_search_expression
    'prenii_unnormalized'           ses_preop    'anat'         '.*T2w'
    'prenii_unnormalized_t1'        ses_preop    'anat'         '.*T1w'
    'prenii_unnormalized_pd'        ses_preop    'anat'         '.*PD'
    'tranii_unnormalized'           ses_postop   'anat'         '.*tra.*T2w' % needs refinement
    'sagnii_unnormalized'           ses_postop   'anat'         '.*sag.*T2w' % needs refinement
    'cornii_unnormalized'           ses_postop   'anat'         '.*cor.*T2w' % needs refinement
    'rawctnii_unnormalized'         ses_postop   'ct'           '.*ct'
    'rest'                          ses_preop    'func'         '.*bold'
    'b0'                            ses_preop    'dwi'          '.*b0'
    'fa'                            ses_preop    'dwi'          '.*fa'
    'fa2anat'                       ses_preop    'dwi'          '.*fa2anat'
    'dti'                           ses_preop    'dwi'          '.*diff'
    };

% Override/append with options.prefs.lead2bids_lookup.
if isfield(options.prefs, 'lead2bids_lookup')
    % Override
    [~, ia, ib] = intersect(lead2bids_lookup(:, 1), options.prefs.lead2bids_lookup(:, 1));
    for ov_ix = 1:length(ia)
        for col_ix = 2:size(lead2bids_lookup, 2)
            lead2bids_lookup{ia(ov_ix), col_ix} = options.prefs.lead2bids_lookup{ib(ov_ix), col_ix};
        end
    end
    % Append
    [~, ia] = setdiff(options.prefs.lead2bids_lookup(:, 1), lead2bids_lookup(:, 1));
    lead2bids_lookup = cat(1, lead2bids_lookup, options.prefs.lead2bids_lookup(ia, :));
end

% If the user appended a file extension then we should assume they know
% what they want. Otherwise, and for the entries from the table above, we
% need to add the file extension manually.
for a = 1:size(lead2bids_lookup, 1)
    if ~contains(lead2bids_lookup{a, 4}, '.nii')
        % Must have .nii and optionally .gz
        lead2bids_lookup{a, 4} = [lead2bids_lookup{a, 4} '\.nii(\.gz)?'];
    end
end


%% Define auxiliary files to copy
% The first column is the sequence type from the above table,
% the second column is a cell array of file extensions (without .)
% of files to be copied along with the nifti.
copy_aux_files = {
    'dti'   {'bvec', 'bval'}
};

%% Find images in BIDS structure, copy and gunzip to derivatives_folder
% search for each type of file defined in lead2bids_lookup
for a = 1:size(lead2bids_lookup, 1)
    % Identify all session folders. Usually only 1.
    ses_folder = dir([bids_subject_folder filesep lead2bids_lookup{a, 2} '*']);
    ses_folder = ses_folder(~ismember({ses_folder.name}, {'.', '..'}));
    if ~isempty(ses_folder)
        ses_folder = strcat([bids_subject_folder filesep], {ses_folder([ses_folder(:).isdir]).name}');
    else 
        ses_folder = {bids_subject_folder};
    end
    
    % Search all session folders, storing matching files
    files = [];
    for b = 1:length(ses_folder)
        search_dir = fullfile(ses_folder{b}, lead2bids_lookup{a, 3});
        found_nifti = dir(fullfile(search_dir, '*.nii*'));
        match_nifti = regexpi({found_nifti.name}, lead2bids_lookup{a, 4}, 'match', 'once');
        match_nifti = found_nifti(~cellfun(@isempty, match_nifti));
        if ~isempty(match_nifti)
            files = [files match_nifti];
        end
    end
        
    % Process (copy and gunzip) matching files
    if length(files) > 1
        warning(['More than 1 file for ' lead2bids_lookup{a, 1} ' found!'])
    end
    for c = 1:length(files)
        [~, fname, ext] = fileparts(files(c).name);
        is_gzipped = strcmpi(ext, '.gz');
        if is_gzipped
            [~, fname, ~] = fileparts(fname);
        end
        % rename the file for lead dbs
        suffix = '';
        if c > 1
            suffix = ['_' num2str(c-1)];  % Multiple-matches get appended with _number.
        end
        outname_noext = options.prefs.(lead2bids_lookup{a,1})(1:end-4);
        outname = fullfile(derivatives_folder, [outname_noext suffix '.nii']);
        if is_gzipped
            outname = [outname '.gz'];
        end
        
        % Copy, gunzip, delete gzip
        copyfile(fullfile(files(c).folder, files(c).name), outname);
        if is_gzipped
            gunzip(outname, derivatives_folder)
        end
        delete(outname)

        % Copy auxiliary files. See copy_aux_files table above.
        for aux_ix = 1:size(copy_aux_files, 1)
            if strcmpi(lead2bids_lookup{a, 1}, copy_aux_files{aux_ix, 1})
                for aux_f_ix = 1:length(copy_aux_files{aux_ix, 2})
                    aux_f_ext = ['.' copy_aux_files{aux_ix, 2}{aux_f_ix}];
                    aux_f_src = fullfile(files(c).folder, [fname aux_f_ext]);
                    if exist(aux_f_src, 'file')
                        aux_f_dest = fullfile(options.prefs.patientdir, [outname_noext aux_f_ext]);
                        copyfile(aux_f_src, aux_f_dest);
                    else
                        warning(['No ' aux_f_ext ' file found for ' lead2bids_lookup{a, 1} '.']);
                    end
                end
            end
        end
        
    end
    
end

%% set patientdir to the newly created derivatives folder
options.uipatdirs = {derivatives_folder};
% Change to the derivatives folder. Still not working.
% Andy wanted to correct this.
cd(options.uipatdirs{1})
