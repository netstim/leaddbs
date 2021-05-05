function ea_fiberactivation_stats(stimFolder)
% Summary of fiber activation results

results = ea_regexpdir(stimFolder, 'fiberActivation_.*\.mat$', 0);

for i=1:length(results)
    parse_fiberactivation(results{i});
end


function parse_fiberactivation(resultPath)
% Load result
load(resultPath, 'fibers', 'idx');

[~, nameStr] = fileparts(resultPath);
nameStr = strrep(nameStr, 'fiberActivation_', '');
nameStr = regexprep(nameStr, '_', ' ', 'once');
nameStr(1) = upper(nameStr(1));

% Activated fibers
fiberIdx = unique(fibers(fibers(:,5)==1,4));
if isempty(fiberIdx)
    fiberIdx = 0;
end
fprintf('\n%s - Activated fibers: %d out of %d\n', nameStr, length(fiberIdx), length(idx));

% Non-Activated fibers
fiberIdx = unique(fibers(fibers(:,5)==0,4));
if isempty(fiberIdx)
    fiberIdx = 0;
end
fprintf('\n%s - Non-activated fibers: %d out of %d\n', nameStr, length(fiberIdx), length(idx));

% Damaged fibers
fiberIdx = unique(fibers(fibers(:,5)==-1,4));
if isempty(fiberIdx)
    fiberIdx = 0;
end
fprintf('\n%s - Damaged fibers: %d out of %d\n', nameStr, length(fiberIdx), length(idx));
