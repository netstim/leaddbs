function ea_dbsegment(input, opts)
% Wrapper to run DBSegment pipeline
% Note: the additional options described on
% https://github.com/luxneuroimage/DBSegment are disabled and only the
% default ones used

arguments
    input {mustBeTextScalar}
    opts.outputFolder {mustBeTextScalar} = ''
end

input = GetFullPath(input);

if isfolder(input) % Input is subj folder in LeadDBS BIDS dataset
    bids = BIDSFetcher(regexp(input, '.*(?=derivatives.*)', 'match', 'once'));
    [~, subjID] = fileparts(input);
    subj = bids.getSubj(regexprep(subjID, '^sub-', ''));

    t1 = subj.preopAnat.(subj.AnchorModality).coreg;
    atlasFolder = fullfile(subj.atlasDir, 'DBSegment Atlas (Baniasadi 2023)');
    ea_mkdir(atlasFolder);
    copyfile(t1, atlasFolder);
    inputFolder = atlasFolder;
    outputFolder = atlasFolder;
elseif isfile(input) % Input is a NIfTI file
    if isBIDSFileName(input) && contains(input, ['derivatives', filesep, 'leaddbs'])
        % NIfTI file inside LeadDBS BIDS dataset
        bids = BIDSFetcher(regexp(input, '.*(?=derivatives.*)', 'match', 'once'));
        subjID = regexp(input, ['(?<=leaddbs\', filesep, ')sub-[^\W_]+'] , 'match', 'once');
        subj = bids.getSubj(regexprep(subjID, '^sub-', ''));

        t1 = input;
        atlasFolder = fullfile(subj.atlasDir, 'DBSegment Atlas (Baniasadi 2023)');
        ea_mkdir(atlasFolder);
        copyfile(t1, atlasFolder);
        inputFolder = atlasFolder;
        outputFolder = atlasFolder;
    else
        % Any other NIfTI file
        if isempty(opts.outputFolder)
            error('Please specify ''outputFolder'' as the 2nd parameters!');
        end
        t1 = input;
        inputFolder = fileparts(t1);
        outputFolder = GetFullPath(opts.outputFolder);
    end
else
    error('Input should either be a subj derivatives folder or a file!')
end

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

segcmd = [DBSegment ...
    ' -i ' inputFolder ...
    ' -o ' outputFolder ...
    ' -mp ' modelPath];
condaenv.system(segcmd);

% Clean up
% TODO

% Prepare atlas
if exist('atlasFolder', 'var')
    fprintf('\Building DBSegment atlas...\n\n');
    ea_mkdir(fullfile(atlasFolder, 'lh'));
    ea_mkdir(fullfile(atlasFolder, 'rh'));

    % Sort structures
    % TODO

    leftNucleus = ea_regexpdir(fullfile(imageFolder, 'left'), '^\d+[-_].+\.nii\.gz$');
    leftNucleusNewPath = replace(leftNucleus, fullfile(imageFolder, 'left'), fullfile(imageFolder, 'atlases', 'DBSegment Atlas (Baniasadi 2023)', 'lh'));
    leftNucleusNewPath = replace(leftNucleusNewPath, regexpPattern(['\' filesep '\d+[-_]']), filesep);
    rightNucleus = ea_regexpdir(fullfile(imageFolder, 'right'), '^\d+[-_].+\.nii\.gz$');
    rightNucleusNewPath = replace(rightNucleus, fullfile(imageFolder, 'right'), fullfile(imageFolder, 'atlases', 'DBSegment Atlas (Baniasadi 2023)', 'rh'));
    rightNucleusNewPath = replace(rightNucleusNewPath, regexpPattern(['\' filesep '\d+[-_]']), filesep);
    cellfun(@(src, dst) copyfile(src, dst), leftNucleus, leftNucleusNewPath);
    cellfun(@(src, dst) copyfile(src, dst), rightNucleus, rightNucleusNewPath);
    
    % Crop nucleus
    nucleus = ea_regexpdir(fullfile(imageFolder, 'atlases', 'DBSegment Atlas (Baniasadi 2023)'), '.*\.nii\.gz$');
    cellfun(@ea_autocrop, nucleus);
    
    atlases = genAtlasesStruct(leftNucleusNewPath);
    
    options.root = subj.subjDir;
    options.patientname = '';
    options.native = 1;
    options.reference = subj.preopAnat.(subj.AnchorModality).coreg;
    options.atlasset = 'DBSegment Atlas (Baniasadi 2023)';
    options.atl.can=0;
    options.atl.ptnative=1;
    ea_genatlastable(atlases, imageFolder, options);
end


function atlases = genAtlasesStruct(nucleus)
% Generate template atlases struct 
if isfile([ea_getearoot, 'ext_libs', 'dbsegment', 'atlas_index.mat'])
    load([ea_getearoot, 'ext_libs', 'dbsegment', 'atlas_index.mat'], 'atlases');
else
    [~, nucleusName, nucleusExt] = fileparts(nucleus);
    atlases.names = sort(strcat(nucleusName, nucleusExt))';
    atlases.types = ones(size(atlases.names)) * 3;
    atlases.threshold.type='relative_intensity';
    atlases.threshold.value=0.5;
    atlases.colormap = ea_color_wes('all', length(atlases.names));
    atlases.citation.name = 'DBSegment Atlas (Baniasadi 2023)';
    atlases.citation.short = 'Baniasadi et al. 2023';
    atlases.citation.long = {'Baniasadi, M., Petersen, M.V., Gonçalves, J., Horn, A., Vlasov, V., Hertel, F., Husch, A., 2023. DBSegment: Fast and robust segmentation of deep brain structures considering domain generalization. Hum. Brain Mapp. 44, 762–778. https://doi.org/10.1002/hbm.26097'};
    atlases.presets(1).label = 'Default';
    atlases.presets(1).hide = [9, 11, 12];
    atlases.presets(1).show = setdiff(1:16, atlases.presets(1).hide);
    atlases.presets(1).default = 'relative';
    atlases.presets(2).label = 'Thalamus';
    atlases.presets(2).show = 9;
    atlases.presets(2).hide = setdiff(1:16, atlases.presets(2).show);
    atlases.presets(2).default = 'relative';
    atlases.presets(3).label = 'VL';
    atlases.presets(3).show = 11;
    atlases.presets(3).hide = setdiff(1:16, atlases.presets(3).show);
    atlases.presets(3).default = 'relative';
    atlases.defautset = 1;
    atlases.rebuild = 1;
end

