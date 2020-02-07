function options = ea_read_bids(options, bids_subject_folder, run)

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

%% Get BIDS session names for preop and postop
% Note that user's can specify the same folder for both, or even
% '' for one or both if no session subfolder is used, by modifying
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

%% Define default lookup table to map BIDS names to lead-dbs names.
% Each row is a lead-dbs<-->BIDS mapping
% The columns are:
% lead_dbs_sequence_name BIDS_ses BIDS_subfolder file_search_expression
% Note that the .nii.gz extension is assumed if not provided.
lead2bids_lookup = ...
    {'prenii_unnormalized'          ses_preop    'anat'  '.*T2w'
    'prenii_unnormalized_t1'        ses_preop    'anat'  '.*T1w'
    'prenii_unnormalized_pd'        ses_preop    'anat'  '.*PD'
    'prenii_unnormalized_t2star'    ses_preop    'anat'  '.*T2star'
    'prenii_unnormalized_swi'       ses_preop    'anat'  '.*SWI'
    'prenii_unnormalized_fgatir'    ses_preop    'anat'  '.*FGATIR'
    'tranii_unnormalized'           ses_postop   'anat'  '.*tra.*T2w' % needs refinement
    'sagnii_unnormalized'           ses_postop   'anat'  '.*sag.*T2w' % needs refinement
    'cornii_unnormalized'           ses_postop   'anat'  '.*cor.*T2w' % needs refinement
    'rawctnii_unnormalized'         ses_postop   'ct'	'.*ct'
    'rest'                          ses_preop    'func'	'.*bold'
    'b0'                            ses_preop    'dwi'	'.*b0'
    'fa'                            ses_preop    'dwi'	'.*fa'
    'fa2anat'                       ses_preop    'dwi'   '.*fa2anat'
    'dti'                           ses_preop    'dwi'   '.*diff'
    };

% Any sequences that are not already defined by lead-dbs should be added
% to the prefs.
options.prefs.prenii_unnormalized_t2star = 'anat_t2star.nii';
options.prefs.prenii_unnormalized_swi = 'anat_swi.nii';
options.prefs.prenii_unnormalized_fgatir = 'anat_fgatir.nii';

% Append/override with options.prefs.lead2bids_lookup.
% Triple-nested-loop feels bad. This could be a perf problem.
if isfield(options.prefs, 'lead2bids_lookup')
    for a = 1:size(lead2bids_lookup, 1)
        for b = 1:size(options.prefs.lead2bids_lookup)
            if strcmpi(lead2bids_lookup{a, 1}, options.prefs.lead2bids_lookup{b, 1})
                for col_ix = 1:length({options.prefs.lead2bids_lookup{b, :}})
                    lead2bids_lookup{a, col_ix} = options.prefs.lead2bids_lookup{b, col_ix};
                end
            end
        end
    end
end

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
    
    % Fixup the extension on the filename regexp, if necessary.
    pattern = lead2bids_lookup{a, 4};
    if (length(pattern) > 7)  && (strcmpi(pattern(end-6:end), '.nii.gz'))
        pattern = pattern(1:end-7);
    end
    if (length(pattern) < 9) || (~strcmpi(pattern(end-8:end), '\.nii\.gz'))
        pattern = [pattern '\.nii\.gz'];
    end
    
    % Search all session folders, storing matching files
    files = [];
    for b = 1:length(ses_folder)
        search_dir = fullfile(ses_folder{b}, lead2bids_lookup{a, 3});
        found_nifti = dir(fullfile(search_dir, '*.nii.gz'));
        match_nifti = regexpi({found_nifti.name}, pattern, 'match', 'once');
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
        inname_noext = files(c).name(1:end-7);
        % rename the file for lead dbs
        suffix = '';
        if c > 1
            suffix = ['_' num2str(c-1)];  % Multiple-matches get appended with _number.
        end
        outname_noext = options.prefs.(lead2bids_lookup{a,1})(1:end-4);
        outname = fullfile(derivatives_folder, [outname_noext suffix '.nii.gz']);
        
        % Copy, gunzip, delete gzip
        copyfile(fullfile(files(c).folder, files(c).name), outname);
        gunzip(outname, derivatives_folder)
        delete(outname)

        % Copy auxiliary files. See copy_aux_files table above.
        for aux_ix = 1:size(copy_aux_files, 1)
            if strcmpi(lead2bids_lookup{a, 1}, copy_aux_files{aux_ix, 1})
                for aux_f_ix = 1:length(copy_aux_files{aux_ix, 2})
                    aux_f_ext = ['.' copy_aux_files{aux_ix, 2}{aux_f_ix}];
                    aux_f_src = fullfile(files(c).folder, [inname_noext aux_f_ext]);
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
cd(options.uipatdirs{1}) % this should 

if exist('run','var') && run % check if lead dbs should run, not sure this is working correctly

options.dolc = 0;
options.ecog.extractsurface.do = 1;
options.ecog.extractsurface.method = 1;
options.endtolerance = 10;
options.sprungwert = 4;
options.refinesteps = 0;
options.tra_stdfactor = 0.9;
options.cor_stdfactor = 1;
options.earoot = ea_getearoot;
options.dicomimp.do = 0;
options.dicomimp.method = 1;
options.assignnii = 0;
options.normalize.do = true;
options.normalize.settings = [];
options.normalize.method = 'ea_normalize_ants';
options.normalize.methodn = 9;
options.normalize.check = false;
options.normalize.refine = 0;
options.coregmr.check = 0;
options.coregmr.method = 'SPM';
options.coregmr.do = 1;
options.overwriteapproved = 0;
options.coregct.do = true;
options.coregct.method = 'ea_coregctmri_ants';
options.coregct.methodn = 7;
options.modality = 1;
options.verbose = 3;
options.sides = [1 2];
options.doreconstruction = false;
options.maskwindow = 10;
options.automask = 1;
options.autoimprove = 0;
options.axiscontrast = 8;
options.zresolution = 10;
options.atl.genpt = false;
options.atl.normalize = 0;
options.atl.can = true;
options.atl.pt = 0;
options.atl.ptnative = false;
options.native = 0;
options.d2.col_overlay = 1;
options.d2.con_overlay = 1;
options.d2.con_color = [1 1 1];
options.d2.lab_overlay = 0;
options.d2.bbsize = 50;
options.d2.backdrop = 'MNI_ICBM_2009b_NLIN_ASYM T1';
options.d2.fid_overlay = 1;
options.d2.write = false;
options.d2.atlasopacity = 0.15;
options.d2.writeatlases = 1;
options.manualheightcorrection = false;
options.scrf.do = 1;
options.scrf.mask = 2;
options.d3.write = false;
options.d3.prolong_electrode = 2;
options.d3.verbose = 'on';
options.d3.elrendering = 1;
options.d3.exportBB = 0;
options.d3.hlactivecontacts = 0;
options.d3.showactivecontacts = 1;
options.d3.showpassivecontacts = 1;
options.d3.showisovolume = 0;
options.d3.isovscloud = 0;
options.d3.mirrorsides = 0;
options.d3.autoserver = 0;
options.d3.expdf = 0;
options.d3.writeatlases = 1;
options.numcontacts = 4;
options.entrypoint = 'STN, GPi or ViM';
options.entrypointn = 1;
options.writeoutpm = 1;
options.elmodeln = 1;
options.elmodel = 'Medtronic 3389';
options.atlasset = 'DISTAL Minimal (Ewert 2017)';
options.atlassetn = 13;
options.reconmethod = 'Refined TRAC/CORE';
options.expstatvat.do = 0;
options.fiberthresh = 10;
options.writeoutstats = 1;
options.leadprod = 'dbs';
ea_run('run',options)
end
