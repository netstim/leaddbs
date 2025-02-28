function ea_thomas(inputImage, imageType)
% Wrapper function to run THOMAS segmentation

arguments
    inputImage {mustBeFile}
    imageType {mustBeMember(imageType, ["t1", 'wmn', 'fgatir'])} = "t1"
end

% Check docker image
dockerImage = 'anagrammarian/thomasmerged:latest';
ea_checkDocker(dockerImage);

% Prepare parameters for docker run
inputImage = GetFullPath(inputImage);
[imageFolder, imageName, imageExt] = fileparts(inputImage);
imageName = [imageName, imageExt];

if lower(imageType) == "t1"
    typeParam = '-t1';
else
    typeParam = '';
end

%% Run segmentation via docker
fprintf('\nRunning THOMAS segmentation...\n\n');
system(['docker run ', ...
        '--volume ', ea_path_helper(imageFolder), ':/thomas '...
        '--workdir /thomas '...
        '--rm -t ', dockerImage, ' ', ...
        'bash -c "hipsthomas_csh -i ', imageName, ' ', typeParam, '"']);

% Clean up output folder
ea_delete(fullfile(imageFolder, 'temp*'));
ea_delete(fullfile(imageFolder, 'left', 'MV'));
ea_delete(fullfile(imageFolder, 'left', 'jf5-VLa.nii.gz'));
ea_delete(fullfile(imageFolder, 'left', 'san_6-VLP.nii.gz'));
ea_delete(fullfile(imageFolder, 'right', 'MV'));
ea_delete(fullfile(imageFolder, 'right', 'jf5-VLa.nii.gz'));
ea_delete(fullfile(imageFolder, 'right', 'san_6-VLP.nii.gz'));

%% Build combined parcellation
leftParc = fullfile(imageFolder, 'left', 'thomasfull.nii.gz');
rightParc = fullfile(imageFolder, 'right', 'thomasrfull.nii.gz');
if isfile(leftParc) && isfile(rightParc)
    fprintf('\nBuilding THOMAS parcellation...\n\n');
    leftParc = ea_load_nii(leftParc);
    rightParc = ea_load_nii(rightParc);

    rightParc.img(rightParc.img>0) = rightParc.img(rightParc.img>0) + 100;
    leftParc.img = leftParc.img + rightParc.img;

    ea_mkdir(fullfile(imageFolder, 'labeling'));
    leftParc.fname = fullfile(imageFolder, 'labeling', 'THOMAS (Su 2019).nii');
    leftParc.dt(1) = 2; % uint8

    ea_write_nii(leftParc);
    copyfile(fullfile(ea_getearoot, 'ext_libs', 'THOMAS', 'THOMAS (Su 2019).txt'), fullfile(imageFolder, 'labeling'));
end

%% Build atlas
fprintf('\Building THOMAS atlas...\n\n');
ea_mkdir(fullfile(imageFolder, 'atlases', 'HIPS-THOMAS Segmentation (Vidal 2024)', 'lh'));
ea_mkdir(fullfile(imageFolder, 'atlases', 'HIPS-THOMAS Segmentation (Vidal 2024)', 'rh'));

leftNucleus = ea_regexpdir(fullfile(imageFolder, 'left'), '^\d+[-_].+\.nii\.gz$');
leftNucleusNewPath = replace(leftNucleus, fullfile(imageFolder, 'left'), fullfile(imageFolder, 'atlases', 'HIPS-THOMAS Segmentation (Vidal 2024)', 'lh'));
leftNucleusNewPath = replace(leftNucleusNewPath, regexpPattern(['\' filesep '\d+[-_]']), filesep);
rightNucleus = ea_regexpdir(fullfile(imageFolder, 'right'), '^\d+[-_].+\.nii\.gz$');
rightNucleusNewPath = replace(rightNucleus, fullfile(imageFolder, 'right'), fullfile(imageFolder, 'atlases', 'HIPS-THOMAS Segmentation (Vidal 2024)', 'rh'));
rightNucleusNewPath = replace(rightNucleusNewPath, regexpPattern(['\' filesep '\d+[-_]']), filesep);
cellfun(@(src, dst) copyfile(src, dst), leftNucleus, leftNucleusNewPath);
cellfun(@(src, dst) copyfile(src, dst), rightNucleus, rightNucleusNewPath);

% Crop nucleus
nucleus = ea_regexpdir(fullfile(imageFolder, 'atlases', 'HIPS-THOMAS Segmentation (Vidal 2024)'), '.*\.nii\.gz$');
cellfun(@ea_autocrop, nucleus);

atlases = genAtlasesStruct(leftNucleusNewPath);

options.root = imageFolder;
options.patientname = '';
options.native = 1;
options.reference = inputImage;
options.atlasset = 'HIPS-THOMAS Segmentation (Vidal 2024)';
options.atl.can=0;
options.atl.ptnative=1;
ea_genatlastable(atlases, imageFolder, options);


function atlases = genAtlasesStruct(nucleus)
% Generate template atlases struct
if isfile(fullfile(ea_getearoot, 'ext_libs', 'THOMAS', 'atlas_index.mat'))
    load(fullfile(ea_getearoot, 'ext_libs', 'THOMAS', 'atlas_index.mat'), 'atlases');
end

[~, nucleusName, nucleusExt] = fileparts(nucleus);
names = sort(strcat(nucleusName, nucleusExt))';

if ~exist('atlases', 'var') || length(names) ~= length(atlases.names) || ~all(strcmp(names, atlases.names))
    atlases = struct;
    atlases.names = names;
    atlases.types = ones(size(atlases.names)) * 3;
    atlases.threshold.type='relative_intensity';
    atlases.threshold.value=0.5;
    atlases.colormap = ea_color_wes('all', length(atlases.names));

    atlases.citation.name = 'HIPS-THOMAS Segmentation (Vidal 2024)';
    atlases.citation.short = 'Vidal et al. 2024';
    atlases.citation.long = {'Vidal, J. P., Danet, L., Péran, P., Pariente, J., Bach Cuadra, M., Zahr, N. M., Barbeau, E. J., & Saranathan, M. (2024). Robust thalamic nuclei segmentation from T1-weighted MRI using polynomial intensity transformation. Brain Structure and Function. https://doi.org/10.1007/s00429-024-02777-5'};

    atlases.presets(1).label = 'Default';
    atlases.presets(1).hide = find(ismember(names, {'THALAMUS.nii.gz', 'VL.nii.gz', 'VLP.nii.gz'})); % Hide thalamus, VL and VLP, only show sub-regions
    atlases.presets(1).show = setdiff(1:length(names), atlases.presets(1).hide);
    atlases.presets(1).default = 'relative';
    atlases.defautset = 1;

    atlases.rebuild = 1;
end
