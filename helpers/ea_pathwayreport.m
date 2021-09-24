function ea_pathwayreport(stimFolder, sortPathway, pathwayTable)
% Generate pathway activation report

fibers = [];

% Load result
results = ea_regexpdir(stimFolder, 'fiberactivation.*\.mat$', 0);
for i=1:length(results)
    result = load(results{i}, 'fibers');
    fibers = [fibers;result.fibers];
end

% Sort pathway according to the fiber counts or not
if ~exist('sortPathway', 'var')
    sortPathway = 1;
end

% Load pathway lookup table
if ~exist('pathwayTable', 'var')
    pathwayTable = [ea_getconnectomebase, 'dMRI', filesep, 'McIntyre', filesep, 'pathway.mat'];
end
load(pathwayTable, 'pathway');

totalFiber = length(unique(fibers(:,4)));

fprintf('\n\n');

% Show activated pathways
fiberIdx = unique(fibers(fibers(:,5)==1,4));
[activated, ~, ic] = unique(pathway(fiberIdx));
if sortPathway
    [fiberCount, sortInd] = sort(accumarray(ic,1), 'descend');
    activated = activated(sortInd);
else
    fiberCount = accumarray(ic,1);
end
activatedTable = table(string(activated),fiberCount, 'VariableNames', ...
    {'Activated Pathways',['Fiber Counts (',num2str(length(fiberIdx)),'/',num2str(totalFiber),')']});
disp(activatedTable);

% Show non-activated pathways
fiberIdx = unique(fibers(fibers(:,5)==0,4));
[nonactivated, ~, ic] = unique(pathway(fiberIdx));
if sortPathway
    [fiberCount, sortInd] = sort(accumarray(ic,1), 'descend');
    nonactivated = nonactivated(sortInd);
else
    fiberCount = accumarray(ic,1);
end
nonactivatedTable = table(string(nonactivated),fiberCount, 'VariableNames', ...
    {'Non-activated Pathways',['Fiber Counts (',num2str(length(fiberIdx)),'/',num2str(totalFiber),')']});
disp(nonactivatedTable);

% Show damaged pathways
fiberIdx = unique(fibers(fibers(:,5)==-1,4));
[damaged, ~, ic] = unique(pathway(fiberIdx));
if sortPathway
    [fiberCount, sortInd] = sort(accumarray(ic,1), 'descend');
    damaged = damaged(sortInd);
else
    fiberCount = accumarray(ic,1);
end
damagedTable = table(string(damaged),fiberCount, 'VariableNames', ...
    {'Damaged Pathways',['Fiber Counts (',num2str(length(fiberIdx)),'/',num2str(totalFiber),')']});
disp(damagedTable);

fprintf('\n');
