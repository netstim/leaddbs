function ea_dbsegment(input_folder, outputFolder)
% Wrapper to run DBSegment pipeline
% Note: the additional options described on
% https://github.com/luxneuroimage/DBSegment are disabled and only the
% default ones used

input = GetFullPath(input_folder);

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

dbsegmentFolder = fullfile(ea_getearoot, 'ext_libs', 'dbsegment');

% Check Conda installation
if ~ea_conda.is_installed
    ea_conda.install;
end

% Check Conda environment
condaenv = ea_conda_env('DBSegment');
if ~condaenv.is_created
    ea_cprintf('CmdWinWarnings', 'Initializing DBSegment reconsurf environment...\n')
    condaenv.create;
    ea_cprintf('CmdWinWarnings', 'DBSegment conda environment initialized.\n')
end

% % Check DBSegment
% runner = fullfile(dbsegmentFolder, 'upstream', 'run_dbsegment.sh');
% if ~isfile(runner)
%     ea_cprintf('CmdWinWarnings', 'Downloading DBSegment...\n')
%     downloadFile = fullfile(dbsegmentFolder, 'DBSegment.zip');
%     websave(downloadFile, 'https://github.com/Deep-MI/FastSurfer/archive/refs/heads/stable.zip'); % What is this? Do I need to update it?
%     unzip(downloadFile, dbsegmentFolder);
%     movefile(fullfile(dbsegmentFolder, 'DBSegment-stable'), fullfile(dbsegmentFolder, 'upstream'));
%     delete(downloadFile);
%     ea_cprintf('CmdWinWarnings', 'DBSegment downloaded.\n')
% end

segcmd = [DBSegment -i input_folder -o output_folder];


condaenv.system(segcmd);
