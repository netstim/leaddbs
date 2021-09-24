function ea_fiberactivation_stats(stimFolder)
% Summary of fiber activation results

results = ea_regexpdir(stimFolder, 'fiberactivation.*\.mat$', 0);

for i=1:length(results)
    parse_fiberactivation(results{i});
end


function parse_fiberactivation(resultPath)
% Load result
load(resultPath, 'fibers', 'idx');

hemiSuffix = regexp(resultPath, '(?<=fiberactivation_hemi-).+(?=\.mat$)', 'match', 'once');
hemiSuffix = strrep(hemiSuffix, '_', ' ');
hemiSuffix(1) = upper(hemiSuffix(1));

% Activated fibers
fiberIdx = unique(fibers(fibers(:,5)==1,4));
if isempty(fiberIdx)
    fiberIdx = 0;
end
fprintf('\n%s - Activated fibers: %d out of %d\n', hemiSuffix, length(fiberIdx), length(idx));

% Non-Activated fibers
fiberIdx = unique(fibers(fibers(:,5)==0,4));
if isempty(fiberIdx)
    fiberIdx = 0;
end
fprintf('\n%s - Non-activated fibers: %d out of %d\n', hemiSuffix, length(fiberIdx), length(idx));

% Damaged fibers
fiberIdx = unique(fibers(fibers(:,5)==-1,4));
if isempty(fiberIdx)
    fiberIdx = 0;
end
fprintf('\n%s - Damaged fibers: %d out of %d\n', hemiSuffix, length(fiberIdx), length(idx));
