function ea_fiberactivation_stats(stimFolder)
% Summary of fiber activation results

% Load result from left side
if isfile([stimFolder, filesep, 'fiberActivation_left.mat'])
    load([stimFolder, filesep, 'fiberActivation_left.mat'], 'fibers', 'idx');
end

% Activated fibers
fiberIdx = unique(fibers(fibers(:,5)==1,4));
if isempty(fiberIdx)
    fiberIdx = 0;
end
fprintf('\nLeft side, Activated fibers: %d out of %d\n', length(fiberIdx), length(idx));

% Non-Activated fibers
fiberIdx = unique(fibers(fibers(:,5)==0,4));
if isempty(fiberIdx)
    fiberIdx = 0;
end
fprintf('\nLeft side, Non-activated fibers: %d out of %d\n', length(fiberIdx), length(idx));

% Damaged fibers
fiberIdx = unique(fibers(fibers(:,5)==-1,4));
if isempty(fiberIdx)
    fiberIdx = 0;
end
fprintf('\nLeft side, Damaged fibers: %d out of %d\n', length(fiberIdx), length(idx));

% Load result from right side
if isfile([stimFolder, filesep, 'fiberActivation_right.mat'])
    load([stimFolder, filesep, 'fiberActivation_right.mat'], 'fibers', 'idx');
end

% Activated fibers
fiberIdx = unique(fibers(fibers(:,5)==1,4));
if isempty(fiberIdx)
    fiberIdx = 0;
end
fprintf('\n\nRight side, Activated fibers: %d out of %d\n', length(fiberIdx), length(idx));

% Non-activated fibers
fiberIdx = unique(fibers(fibers(:,5)==0,4));
if isempty(fiberIdx)
    fiberIdx = 0;
end
fprintf('\nRight side, Non-activated fibers: %d out of %d\n', length(fiberIdx), length(idx));

% Damaged fibers
fiberIdx = unique(fibers(fibers(:,5)==-1,4));
if isempty(fiberIdx)
    fiberIdx = 0;
end
fprintf('\nRight side, Damaged fibers: %d out of %d\n\n', length(fiberIdx), length(idx));
