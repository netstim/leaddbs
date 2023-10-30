function ea_fastsurfer_seg(input, subjID, outputFolder, updateScript)
% Wrapper to run FastSurfer seg_only pipeline

if ~exist('updateScript', 'var')
    updateScript = false; % Update FastSurfer script or not
end

input = GetFullPath(input);

if isfolder(input) % Input is subj folder in LeadDBS BIDS dataset
    [~, subjID] = fileparts(input);
    bids = BIDSFetcher(regexp(input, '.*(?=derivatives.*)', 'match', 'once'));
    subj = bids.getSubj(regexprep(subjID, '^sub-', ''));
    t1 = subj.preopAnat.(subj.AnchorModality).coreg;
    outputFolder = subj.freesurferDir;
elseif isfile(input) % Input is a NIfTI file
    if isBIDSFileName(input) && contains(input, ['derivatives', filesep, 'leaddbs'])
        % NIfTI file inside LeadDBS BIDS dataset
        subjID = regexp(input, ['(?<=leaddbs\', filesep, ')sub-[^\W_]+'] , 'match', 'once');
        bids = BIDSFetcher(regexp(input, '.*(?=derivatives.*)', 'match', 'once'));
        LeadDBSDirs = bids.getLeadDBSDirs(regexprep(subjID, '^sub-', ''));
        t1 = input;
        outputFolder = LeadDBSDirs.freesurferDir;
    else
        % Any other NIfTI file
        if ~exist('subjID', 'var') || ~exist('outputFolder', 'var')
            error('Please specify ''subjID'' and ''outputFolder'' as the 2nd and 3rd parameters!');
        end
        t1 = input;
        outputFolder = GetFullPath(outputFolder);
    end
end

fastsurferFolder = fullfile(ea_prefsdir, 'fastsurfer');

% Check Conda installation
if ~ea_conda.is_installed
    ea_conda.install;
end

% Check Conda environment
condaenv = ea_conda_env('FastSurfer');
if ~condaenv.is_created
    ea_cprintf('CmdWinWarnings', 'Initializing FastSurfer reconsurf environment...\n')
    condaenv.create;
    ea_cprintf('CmdWinWarnings', 'FastSurfer conda environment initialized.\n')
elseif ~condaenv.is_up_to_date
    ea_cprintf('CmdWinWarnings', 'Updating FastSurfer reconsurf environment...\n')
    condaenv.update;
    ea_cprintf('CmdWinWarnings', 'FastSurfer conda environment initialized.\n')
end

% Check FastSurfer
runner = fullfile(fastsurferFolder, 'run_fastsurfer.sh');
if ~isfile(runner) || updateScript
    ea_cprintf('CmdWinWarnings', 'Downloading FastSurfer...\n')
    ea_delete(fastsurferFolder);
    downloadFile = fullfile(ea_prefsdir, 'FastSurfer.zip');
    websave(downloadFile, 'https://github.com/Deep-MI/FastSurfer/archive/refs/heads/stable.zip');
    unzip(downloadFile, ea_prefsdir);
    movefile(fullfile(ea_prefsdir, 'FastSurfer-stable'), fastsurferFolder);
    ea_delete(downloadFile);
    ea_cprintf('CmdWinWarnings', 'FastSurfer downloaded.\n')
end

segcmd = ['export FASTSURFER_HOME=', fastsurferFolder, ';', ...
    'bash ', runner, ' --sid ', subjID, ' --sd ', outputFolder, ' --t1 ', t1, ...
    ' --vox_size min --seg_only --viewagg_device ''cpu'''];

condaenv.system(segcmd);
