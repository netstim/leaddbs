function ea_pathwayreport(axonActivationFolder, pathwayTable)
% Generate pathway activation report

axons = [];

% Load result from left side
if isfile([axonActivationFolder, filesep, 'axonActivation_left.mat'])
    load([axonActivationFolder, filesep, 'axonActivation_left.mat'], 'fibers');
    axons = [axons;fibers];
end

% Load result from right side
if isfile([axonActivationFolder, filesep, 'axonActivation_right.mat'])
    load([axonActivationFolder, filesep, 'axonActivation_right.mat'], 'fibers');
    axons = [axons;fibers];
end

% Load pathway lookup table
if ~exist('pathwayTable', 'var')
    pathwayTable = [ea_getconnectomebase, 'dMRI', filesep, 'McIntyre', filesep, 'pathway.mat'];
end
load(pathwayTable, 'pathway')

% Show activated pathways
axonIdx = unique(axons(axons(:,5)==1,4));
activated = unique(pathway(axonIdx));
fprintf('\nActivated pathways:\n');
fprintf('%s\n',activated{:});

% Show non-activated pathways
axonIdx = unique(axons(axons(:,5)==0,4));
nonactivated = unique(pathway(axonIdx));
fprintf('\nNon-activated pathways:\n');
fprintf('%s\n',nonactivated{:});

% Show damaged pathways
axonIdx = unique(axons(axons(:,5)==-1,4));
damaged = unique(pathway(axonIdx));
fprintf('\nDamaged pathways:\n');
fprintf('%s\n',damaged{:});

fprintf('\n');
