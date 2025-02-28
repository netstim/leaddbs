function ea_atlas2conn(atlasFolder, connectomeName)
% Convert atlas that contains fiber tracts into normal dMRI connectome

% Retrieve left side tracts
leftTracts = ea_regexpdir([atlasFolder, filesep, 'lh'], '^(?!\.).*\.mat$', 1);
[~, leftTractNames] = cellfun(@fileparts, leftTracts, 'Uni', 0);
leftTractNames = strcat({'Left '}, leftTractNames);

% Retrieve right side tracts
rightTracts = ea_regexpdir([atlasFolder, filesep, 'rh'], '^(?!\.).*\.mat$', 1);
[~, rightTractNames] = cellfun(@fileparts, rightTracts, 'Uni', 0);
rightTractNames = strcat({'Right '}, rightTractNames);

% Retrieve "mixed" tracts, if they exist, eg commissural tracts
mixedTracts = ea_regexpdir([atlasFolder, filesep, 'mixed'], '^(?!\.).*\.mat$', 1);
[~, mixedTractNames] = cellfun(@fileparts, mixedTracts, 'Uni', 0); % no need to add a prefix here

% Create connectome folder if necessary
connectomeFolder = [ea_getconnectomebase('dMRI'), connectomeName];
if ~isfolder(connectomeFolder)
    mkdir(connectomeFolder);
end

% Aggregate fiber tracts tp create a normal dMRI connectome
[~, ~, partition] = ea_ftr_aggregate([leftTracts; rightTracts; mixedTracts], [connectomeFolder, filesep, 'data.mat']);

% Save pathway names
pathway = repelem([leftTractNames; rightTractNames; mixedTractNames], partition);
save([connectomeFolder, filesep, 'pathway.mat'], 'pathway');
