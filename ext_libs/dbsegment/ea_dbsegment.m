function ea_dbsegment(input)
% Wrapper to run DBSegment pipeline
% Note: the additional options described on
% https://github.com/luxneuroimage/DBSegment are disabled and only the
% default ones used

input = GetFullPath(input);

if isfolder(input) % Input is subj folder in LeadDBS BIDS dataset
    bids = BIDSFetcher(regexp(input, '.*(?=derivatives.*)', 'match', 'once'));
    [~, subjID] = fileparts(input);
    subj = bids.getSubj(regexprep(subjID, '^sub-', ''));

    t1 = subj.coreg.anat.preop.(subj.AnchorModality);
    atlasFolder = fullfile(subj.atlasDir, 'DBSegment Atlas (Baniasadi 2023)');
    inputFolder = fullfile(atlasFolder, 'input');
    outputFolder = fullfile(atlasFolder, 'output');
elseif isfile(input) % Input is a NIfTI file
    if isBIDSFileName(input) && contains(input, ['derivatives', filesep, 'leaddbs'])
        % NIfTI file inside LeadDBS BIDS dataset
        bids = BIDSFetcher(regexp(input, '.*(?=derivatives.*)', 'match', 'once'));
        subjID = regexp(input, ['(?<=leaddbs\', filesep, ')sub-[^\W_]+'] , 'match', 'once');
        subj = bids.getSubj(regexprep(subjID, '^sub-', ''));

        t1 = input;
        atlasFolder = fullfile(subj.atlasDir, 'DBSegment Atlas (Baniasadi 2023)');
        inputFolder = fullfile(atlasFolder, 'input');
        outputFolder = fullfile(atlasFolder, 'output');
    else
        % Any other NIfTI file
        t1 = input;
        randomID = char(randi([65 90],1,4));
        inputFolder = fullfile(fileparts(t1), ['DBSegmentInput_' randomID]);
        outputFolder = fullfile(fileparts(t1), ['DBSegmentOutput_' randomID]);
    end
else
    error('Input should either be a subj derivatives folder or a file!')
end

% Copy file to temporary input folder
ea_mkdir(inputFolder);
copyfile(t1, inputFolder);

% Check Conda environment
condaenv = ea_conda_env('DBSegment');
if ~condaenv.is_created
    ea_cprintf('CmdWinWarnings', 'Initializing DBSegment conda environment...\n')
    condaenv.create;
    ea_cprintf('CmdWinWarnings', 'DBSegment conda environment initialized.\n')
elseif ~condaenv.is_up_to_date
    ea_cprintf('CmdWinWarnings', 'Updating DBSegment conda environment...\n')
    condaenv.update;
    ea_cprintf('CmdWinWarnings', 'DBSegment conda environment initialized.\n')
end

% Run check DBSegment
DBSegment = fullfile(condaenv.path, 'bin', 'DBSegment');
modelPath = fullfile(condaenv.path, 'share', 'models');

segcmd = {ea_path_helper(DBSegment), ...
    '-i', ea_path_helper(inputFolder), ...
    '-o', ea_path_helper(outputFolder), ...
    '-mp', ea_path_helper(modelPath)};

status = condaenv.system(strjoin(segcmd, ' '));
if status ~= 0
    ea_error('DBSegment failed!', showdlg=false, simpleStack=true);
end

% Sort results, clean up
[~, t1FileName] = ea_niifileparts(t1);
outputFile = fullfile(outputFolder, [t1FileName '.nii.gz']);
if ~exist('atlasFolder', 'var')
    copyfile(outputFile, fullfile(fileparts(t1), [t1FileName '_seg.nii.gz']));
else
    gunzip(outputFile, atlasFolder);
end

ea_delete(inputFolder);
ea_delete(outputFolder);

% Prepare atlas
if exist('atlasFolder', 'var')
    ea_cprintf('*Comments', 'Building DBSegment atlas...\n\n');

    % Load names
    names = readcell(fullfile(ea_getearoot, 'ext_libs', 'DBSegment', 'dbsegment_label.txt'));
    names = names(:,2);

    ea_mkdir(fullfile(atlasFolder, 'lh'));
    ea_mkdir(fullfile(atlasFolder, 'rh'));

    % Export nucleus
    nii = ea_load_nii(fullfile(atlasFolder, t1FileName));
    label = setdiff(unique(nii.img), [0 1]);
    nucleus = nii;
    nucleus.dt(1) = 2;
    for i=1:length(label)
        nucleus.img(:) = 0;
        nucleus.img(nii.img==label(i)) = 1;
        if endsWith(names{label(i)}, '-L')
            nucleus.fname = fullfile(atlasFolder, 'lh', replace(names{label(i)}, '-L', '.nii'));
        elseif endsWith(names{label(i)}, '-R')
            nucleus.fname = fullfile(atlasFolder, 'rh', replace(names{label(i)}, '-R', '.nii'));
        end

        ea_write_nii(nucleus);
        ea_autocrop(nucleus.fname);
        gzip(nucleus.fname);
        ea_delete(nucleus.fname);
    end

    ea_delete(nii.fname);

    nucleusFiles = ea_regexpdir(fullfile(atlasFolder, 'lh'), '.*\.nii(\.gz)?$');

    % Get atlas_index structure
    atlases = genAtlasesStruct(nucleusFiles);

    % Create atlas
    options.root = subj.subjDir;
    options.patientname = '';
    options.native = 1;
    options.reference = subj.preopAnat.(subj.AnchorModality).coreg;
    options.atlasset = 'DBSegment Atlas (Baniasadi 2023)';
    options.atl.can=0;
    options.atl.ptnative=1;
    ea_genatlastable(atlases, atlasFolder, options);
end


function atlases = genAtlasesStruct(nucleusFiles)
% Generate template atlases struct 
if isfile(fullfile(ea_getearoot, 'ext_libs', 'DBSegment', 'atlas_index.mat'))
    load(fullfile(ea_getearoot, 'ext_libs', 'DBSegment', 'atlas_index.mat'), 'atlases');
end

[~, nucleusName, nucleusExt] = fileparts(nucleusFiles);
names = sort(strcat(nucleusName, nucleusExt))';

if ~exist('atlases', 'var') || length(names) ~= length(atlases.names) || ~all(strcmp(names, atlases.names))
    atlases = struct;
    atlases.names = names;
    atlases.types = ones(size(atlases.names)) * 3;
    atlases.threshold.type='relative_intensity';
    atlases.threshold.value=0.5;
    atlases.colormap = ea_color_wes('all', length(atlases.names));

    atlases.citation.name = 'DBSegment Atlas (Baniasadi 2023)';
    atlases.citation.short = 'Baniasadi et al. 2023';
    atlases.citation.long = {'Baniasadi, M., Petersen, M.V., Gonçalves, J., Horn, A., Vlasov, V., Hertel, F., Husch, A., 2023. DBSegment: Fast and robust segmentation of deep brain structures considering domain generalization. Hum. Brain Mapp. 44, 762–778. https://doi.org/10.1002/hbm.26097'};

    
    atlases.presets(1).label = 'Default';
    atlases.presets(1).show = find(ismember(names, {'GPI.nii.gz', 'GPE.nii.gz', 'RN.nii.gz', 'STN.nii.gz'}));
    atlases.presets(1).hide = setdiff(1:length(names), atlases.presets(1).show);
    atlases.presets(1).default = 'relative';
    atlases.defautset = 1;

    atlases.rebuild = 1;
end

