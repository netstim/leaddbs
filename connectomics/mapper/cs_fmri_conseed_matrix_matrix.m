function cs_fmri_conseed_matrix_matrix(dfold,cname,sfile,outputfolder)
% Create fMRI connectivity matrix based on matrix type of connectome

tic

if ~exist('dfold','var')
    dfold=''; % assume all data needed is stored here.
else
    if ~strcmp(dfold(end),filesep)
        dfold=[dfold,filesep];
    end
end

disp(['Connectome dataset: ',cname,'.']);
if contains(cname, '>')
    subset = regexprep(cname, '.*> *', '');
    cname = regexprep(cname, ' *>.*', '');
end

try
    options = evalin('caller','options');
    if strcmp(options.lcm.seeddef, 'parcellation')
        prefix = 'parc_';
        roi = options.lcm.parcSeedName;
    else
        prefix = 'seed_';
        roi = sfile;
    end
catch
    prefix = '';
    roi = sfile;
end

connLabel = ea_getConnLabel(cname, 'FullSet');

dataset = load([dfold,'fMRI',filesep,cname,filesep,'dataset_volsurf.mat']);
dataset.connLabel = connLabel;
dinfo = loadjson([dfold,'fMRI',filesep,cname,filesep,'dataset_info.json']);
dataset.type = dinfo.type;

if isempty(outputfolder)
    warning('off', 'backtrace');
    warning('Custom output folder not specified! Will save result to current folder.');
    warning('on', 'backtrace');
    outputfolder = [pwd, filesep];
elseif ~strcmp(outputfolder(end),filesep)
    outputfolder = [outputfolder,filesep];
end

if length(sfile) == 1 % Manually chosen parcellation: single NIfTI
    parc = ea_load_nii(sfile{1});
    roiID = unique(parc.img(parc.img~=0));
    roiInd = cell(size(roiID));
    for i=1:length(roiID)
        roiInd{i} = intersect(find(parc.img==roiID(i)), dataset.vol.outidx);
        roiInd{i} = sort(arrayfun(@(x) find(dataset.vol.outidx==x), roiInd{i}));
    end
else % Manually chosen seeds: multiple NIfTIs or text file with NIfTI entries
    [~, roiID] = cellfun(@ea_niifileparts, sfile, 'UniformOutput', false);
    roiInd = cell(size(roiID));
    sfile = ea_resolveseeds(sfile, dataset.vol.space);
    for i=1:length(sfile)
        seed = ea_load_nii(sfile{i});
        roiInd{i} = intersect(find(seed.img==1), dataset.vol.outidx);
        roiInd{i} = sort(arrayfun(@(x) find(dataset.vol.outidx==x), roiInd{i}));
    end
end

db = matfile([dfold,'fMRI',filesep,cname,filesep,'AllX.mat'], 'Writable', false);
connMat = nan(length(roiID));

for i=1:length(roiID)
    for j=i:length(roiID)
        block = zeros(length(roiInd{i}), length(roiInd{j}));
        for r=1:size(block,1)
            for c=1:size(block,2)
                block(r,c) = db.X(roiInd{i}(r), roiInd{j}(c));
            end
        end
        connMat(i, j) = mean(block, 'all');
    end
end

conn.roi = roi;
conn.roiID = roiID;
conn.connMat = connMat + triu(connMat', 1);
outputFile = fullfile(outputfolder, [prefix, 'conn-', connLabel, '_funcmatrix.mat']);
i = 1; 
while isfile(outputFile)
    outputFile = strrep(outputFile, '.mat', [num2str(i), '.mat']);
    i = i + 1;
end

ea_cprintf('CmdWinWarnings', 'Saved to %s\n', outputFile);
save(outputFile, '-struct', 'conn', '-v7.3');

toc


function seedfiles = ea_resolveseeds(seedfiles, spacedef)

tmpdir = ea_getleadtempdir;

spacedef.fname = fullfile(tmpdir, ['spacedef_', ea_generate_uuid, '.nii']);
ea_write_nii(spacedef);

for i=1:length(seedfiles)
    [~, ~, ext] = ea_niifileparts(seedfiles{i});
    seedfname = ea_generate_uuid;
    copyfile(seedfiles{i}, fullfile(tmpdir, [seedfname, ext]));
    if strcmp(ext, '.nii.gz')
        gunzip(fullfile(tmpdir, [seedfname, ext]));
        delete(fullfile(tmpdir, [seedfname, ext]));
    end

    ea_conformspaceto(spacedef.fname, fullfile(tmpdir,[seedfname,'.nii']),...
        0, [], fullfile(tmpdir,[seedfname,'.nii']), 0);

    seedfiles{i} = fullfile(tmpdir, [seedfname, '.nii']);
end
