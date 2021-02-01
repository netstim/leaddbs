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
load(pathwayTable, 'pathway');

totalAxon = length(unique(axons(:,4)));

% Show activated pathways
axonIdx = unique(axons(axons(:,5)==1,4));
[activated, ~, ic] = unique(pathway(axonIdx));
axonCount = accumarray(ic,1);
activatedTable = table(string(activated),axonCount, 'VariableNames', ...
    {'Activated Pathways',['Axon Counts (',num2str(length(axonIdx)),'/',num2str(totalAxon),')']});
disp(activatedTable);

% Show non-activated pathways
axonIdx = unique(axons(axons(:,5)==0,4));
[nonactivated, ~, ic] = unique(pathway(axonIdx));
axonCount = accumarray(ic,1);
nonactivatedTable = table(string(nonactivated),axonCount, 'VariableNames', ...
    {'Non-activated Pathways',['Axon Counts (',num2str(length(axonIdx)),'/',num2str(totalAxon),')']});
disp(nonactivatedTable);

% Show damaged pathways
axonIdx = unique(axons(axons(:,5)==-1,4));
[damaged, ~, ic] = unique(pathway(axonIdx));
axonCount = accumarray(ic,1);
damagedTable = table(string(damaged),axonCount, 'VariableNames', ...
    {'Damaged Pathways',['Axon Counts (',num2str(length(axonIdx)),'/',num2str(totalAxon),')']});
disp(damagedTable);

fprintf('\n');
