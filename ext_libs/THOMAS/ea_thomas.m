function ea_thomas(inputImage, imageType)
% Wrapper function to run THOMAS segmentation

arguments
    inputImage {mustBeFile}
    imageType {mustBeMember(imageType, ["t1", 'wmn', 'fgatir'])} = "t1"
end

% Check docker image
dockerImage = 'anagrammarian/sthomas:latest';
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
ea_cprintf('*Comments', 'Running THOMAS segmentation...\n\n');
if isunix && ~ismac
    userParam = '-u $(id -u):$(id -g) ';
else
    userParam = '';
end
system(['docker run -it --rm --name sthomas ', ...
        '-v ', ea_path_helper(imageFolder), ':/data '...
        '-w /data ', ...
        userParam, ...
        '-t ', dockerImage, ' ', ...
        'hipsthomas.sh -i ', imageName, ' -v ', typeParam, '']);

% Clean up output folder
ea_delete(fullfile(imageFolder, 'sthomas_LR_labels.nii.gz'));
ea_delete(fullfile(imageFolder, 'sthomas_LR_labels.png'));
ea_delete(fullfile(imageFolder, 'temp*'));
ea_delete(fullfile(imageFolder, 'left', 'EXTRAS'));
ea_delete(fullfile(imageFolder, 'left', '1-THALAMUS.nii.gz'));
ea_delete(fullfile(imageFolder, 'left', '4567-VL.nii.gz'));
ea_delete(fullfile(imageFolder, 'left', 'CL_L.nii.gz'));
ea_delete(fullfile(imageFolder, 'left', 'crop_*.nii.gz'));
ea_delete(fullfile(imageFolder, 'left', '33-GP.nii.gz'));
ea_delete(fullfile(imageFolder, 'left', 'nucleiVols.txt'));
ea_delete(fullfile(imageFolder, 'left', 'ocrop_*.nii.gz'));
ea_delete(fullfile(imageFolder, 'left', 'regn_L.nii.gz'));
ea_delete(fullfile(imageFolder, 'left', 'VLP.nii.gz'));
ea_delete(fullfile(imageFolder, 'right', 'EXTRAS'));
ea_delete(fullfile(imageFolder, 'right', '1-THALAMUS.nii.gz'));
ea_delete(fullfile(imageFolder, 'right', '4567-VL.nii.gz'));
ea_delete(fullfile(imageFolder, 'right', 'CL_R.nii.gz'));
ea_delete(fullfile(imageFolder, 'right', 'crop_*.nii.gz'));
ea_delete(fullfile(imageFolder, 'right', '33-GP.nii.gz'));
ea_delete(fullfile(imageFolder, 'right', 'nucleiVols.txt'));
ea_delete(fullfile(imageFolder, 'right', 'ocrop_*.nii.gz'));
ea_delete(fullfile(imageFolder, 'right', 'regn_R.nii.gz'));
ea_delete(fullfile(imageFolder, 'right', 'VLP.nii.gz'));

%% Build combined parcellation
leftParcPath = fullfile(imageFolder, 'left', 'thomasfull_L.nii.gz');
leftAmyPath = fullfile(imageFolder, 'left', '34-Amy.nii.gz');
rightParcPath = fullfile(imageFolder, 'right', 'thomasfull_R.nii.gz');
rightAmyPath = fullfile(imageFolder, 'right', '34-Amy.nii.gz');

ea_cprintf('*Comments', 'Building THOMAS parcellation...\n\n');
leftParcNii = ea_load_nii(leftParcPath);
ea_delete(leftParcPath);
rightParcNii = ea_load_nii(rightParcPath);
ea_delete(rightParcPath);
rightParcNii.img(rightParcNii.img>0) = rightParcNii.img(rightParcNii.img>0) + 100;
leftParcNii.img = leftParcNii.img + rightParcNii.img;

leftAmyNii = ea_load_nii(leftAmyPath);
rightAmyNii = ea_load_nii(rightAmyPath);
leftAmyNii.img(leftAmyNii.img==1) = 34;
rightAmyNii.img(rightAmyNii.img==1) = 134;
leftParcNii.img = leftParcNii.img + leftAmyNii.img;
leftParcNii.img = leftParcNii.img + rightAmyNii.img;

ea_mkdir(fullfile(imageFolder, 'labeling'));
leftParcNii.fname = fullfile(imageFolder, 'labeling', 'THOMAS (Su 2019).nii');
leftParcNii.dt(1) = 2; % uint8

ea_write_nii(leftParcNii);
copyfile(fullfile(ea_getearoot, 'ext_libs', 'THOMAS', 'THOMAS (Su 2019).txt'), fullfile(imageFolder, 'labeling'));

%% Build atlas
ea_cprintf('*Comments', 'Building THOMAS atlas...\n\n');
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

ea_delete(fullfile(imageFolder, 'left'));
ea_delete(fullfile(imageFolder, 'right'));

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
    % atlases.colormap = ea_color_wes('all', length(atlases.names));
    atlases.colormap = [48 18 59
                        232 213 56
                        142 10 1
                        81 249 121
                        249 186 56
                        254 152 44
                        247 113 27
                        232 75 12
                        172 250 55
                        23 217 199
                        132 254 80
                        36 236 166
                        255 237 45
                        38 189 233
                        210 48 5
                        180 27 1
                        61 55 145
                        69 91 206
                        90 165 244
                        90 205 244
                        60 157 253]/255;

    atlases.citation.name = 'HIPS-THOMAS Segmentation (Vidal 2024)';
    atlases.citation.short = 'Vidal et al. 2024';
    atlases.citation.long = {'Vidal, J. P., Danet, L., Péran, P., Pariente, J., Bach Cuadra, M., Zahr, N. M., Barbeau, E. J., & Saranathan, M. (2024). Robust thalamic nuclei segmentation from T1-weighted MRI using polynomial intensity transformation. Brain Structure and Function. https://doi.org/10.1007/s00429-024-02777-5'};

    atlases.presets(1).label = 'Default';
    atlases.presets(1).hide = [];
    atlases.presets(1).show = 1:length(names);
    atlases.presets(1).default = 'relative';
    atlases.defautset = 1;

    atlases.rebuild = 1;
end
