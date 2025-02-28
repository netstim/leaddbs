function cs_fmri_conseed_matrix_matrix(dfold,cname,sfile,outputfolder)
% ROI based fMRI connectivity matrix based on matrix type of connectome

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
    cname = regexprep(cname, ' *>.*', '');
end

connLabel = ea_getConnLabel(cname, 'FullSet');

try
    options = evalin('caller','options');
    if strcmp(options.lcm.seeddef, 'parcellation')
        prefix = 'parc_';
        % Parcellation as roi definition
        roi = [fileparts(options.lcm.parcSeedFolder), '.nii'];
    else
        prefix = 'seed_';
        % Seeds as roi definition 
        roi = sfile;
    end
catch
    prefix = '';
    roi = sfile;
end

if isempty(outputfolder)
    ea_cprintf('CmdWinWarnings', 'Custom output folder not specified! Will save result to current folder.\n');
    outputfolder = pwd;
end

outputFile = fullfile(outputfolder, [prefix, 'conn-', connLabel, '_funcmatrix.mat']);
i = 1; 
while isfile(outputFile)
    outputFile = strrep(outputFile, '.mat', [num2str(i), '.mat']);
    i = i + 1;
end

dataset = load([dfold,'fMRI',filesep,cname,filesep,'dataset_volsurf.mat']);
dataset.connLabel = connLabel;
dinfo = loadjson([dfold,'fMRI',filesep,cname,filesep,'dataset_info.json']);
dataset.type = dinfo.type;

if length(sfile) == 1 % Manually chosen parcellation: single NIfTI
    parc = ea_load_nii(sfile{1});
    % Unique value in parcellation as roiID
    roiID = unique(parc.img(parc.img~=0));
    roiInd = cell(size(roiID));
    for i=1:length(roiID)
        % Intersect roi indices and connectome outidx
        roiInd{i} = intersect(find(parc.img==roiID(i)), dataset.vol.outidx);
        % Find corresponding indices of connectome X matrix
        roiInd{i} = sort(arrayfun(@(x) find(dataset.vol.outidx==x), roiInd{i}));
    end
else % Manually chosen seeds: multiple NIfTIs or text file with NIfTI entries
    % NIfTI file name as roiID
    [~, roiID] = cellfun(@ea_niifileparts, sfile, 'UniformOutput', false);
    roiInd = cell(size(roiID));
    seeds = ea_resolveseeds(sfile, dataset);
    for i=1:length(sfile)
        seed = seeds{i};
        % if ~ea_isbinary(seed.img)
        %     ea_error('Creating connectivity matrices using matrix connectomes is only supported for binary seeds. Please binarize input.');
        % end
        % Intersect roi indices and connectome outidx
        roiInd{i} = intersect(find(seed.img), dataset.vol.outidx);
        roiWeight{i} = seed.img(roiInd{i});
        % Find corresponding indices of connectome X matrix
        try
            roiInd{i} = arrayfun(@(x) find(dataset.vol.outidx==x), roiInd{i});
        catch % no one to one mapping between voxels and matrix entries (e.g. cifti connectomes that were converted from surface space)
            outroiInd=cell(1,length(roiInd{i}));
            ourtoiWeight=outroiInd;
            for entry=1:length(roiInd{i}) % need to go one by one
                outroiInd{entry}=find(dataset.vol.outidx==roiInd{i}(entry));
                outroiWeight{entry}=repmat(roiWeight{i}(entry),size(outroiInd{entry}));
            end
            roiInd{i}=cell2mat(outroiInd');
            roiWeight{i}=cell2mat(outroiWeight');
        end
        [roiInd{i},sorting]=sort(roiInd{i}); % do sorting of both ind and weight
        roiWeight{i}=roiWeight{i}(sorting);
        roiWeight{i}=roiWeight{i}/ea_nansum(roiWeight{i});
    end
end

disp('Memory mapping matrix connectome.');
db = matfile([dfold,'fMRI',filesep,cname,filesep,'AllX.mat'], 'Writable', false);
connMat = zeros(length(roiID));
ea_dispercent(0,'Iterating seeds');
for i=1:length(roiID)
    for j=i:length(roiID)
        % First read a larger block from connectome X matrix since ranges
        % for MatFile object must increase in equally spaced intervals
        rowBlockInd = roiInd{i}(1):roiInd{i}(end);
        colBlockInd = roiInd{j}(1):roiInd{j}(end);
        blockConnMat = db.X(rowBlockInd, colBlockInd);

        % Get the real (smaller) connectivity matrix block
        rowInd = arrayfun(@(x) find(rowBlockInd==x), roiInd{i});
        colInd = arrayfun(@(x) find(colBlockInd==x), roiInd{j});
        blockConnMat = blockConnMat(rowInd, colInd);


        % Apply weights (which sum to 1)
        blockConnMat=blockConnMat...
            .*...
            repmat(roiWeight{i},1,size(roiWeight{j},1))...
            .*...
            repmat(roiWeight{j}',size(roiWeight{i},1),1);

        % Grand mean of the connectivity matrix block
        connMat(i, j) = nansum(blockConnMat, 'all'); % nansum instead of nanmean is correct since all weighted down to a sum of 1 before
    end
    ea_dispercent(i/length(roiID));
end
ea_dispercent(1,'end');

conn.roi = roi;
conn.roiID = roiID;

% Symmetrize the roi based connectivity matrix
conn.connMat = connMat + tril(connMat', -1);

ea_cprintf('CmdWinWarnings', 'ROI based connectivity matrix saved to:\n%s\n', outputFile);
save(outputFile, '-struct', 'conn', '-v7.3');

toc


function seeds = ea_resolveseeds(seedfiles, dataset)

% harmonize with the function that cs_fmri_conseed_seed_matrix uses.
for s=1:length(seedfiles)
    seeds{s}=ea_conformseedtofmri(dataset,seedfiles{s});
end

% Old code to do the same:
% tmpdir = ea_getleadtempdir;
% 
% spacedef.fname = fullfile(tmpdir, ['spacedef_', ea_generate_uuid, '.nii']);
% ea_write_nii(spacedef);
% 
% for i=1:length(seedfiles)
%     [~, ~, ext] = ea_niifileparts(seedfiles{i});
%     seedfname = ea_generate_uuid;
%     copyfile(seedfiles{i}, fullfile(tmpdir, [seedfname, ext]));
%     if strcmp(ext, '.nii.gz')
%         gunzip(fullfile(tmpdir, [seedfname, ext]));
%         delete(fullfile(tmpdir, [seedfname, ext]));
%     end
% 
%     ea_conformspaceto(spacedef.fname, fullfile(tmpdir,[seedfname,'.nii']),...
%         0, [], fullfile(tmpdir,[seedfname,'.nii']), 0);
% 
%     seedfiles{i} = fullfile(tmpdir, [seedfname, '.nii']);
% end
