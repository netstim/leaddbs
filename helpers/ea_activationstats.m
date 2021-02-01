function ea_activationstats(axonActivationFolder)
% Summary of axon activation results

% Load result from left side
if isfile([axonActivationFolder, filesep, 'axonActivation_left.mat'])
    load([axonActivationFolder, filesep, 'axonActivation_left.mat'], 'fibers', 'idx');
end

% Activated pathways
axonIdx = unique(fibers(fibers(:,5)==1,4));
if isempty(axonIdx)
    axonIdx = 0;
end
fprintf('\nLeft side, Activated axons: %d out of %d\n', length(axonIdx), length(idx));

% Activated pathways
axonIdx = unique(fibers(fibers(:,5)==0,4));
if isempty(axonIdx)
    axonIdx = 0;
end
fprintf('\nLeft side, Non-activated axons: %d out of %d\n', length(axonIdx), length(idx));

% Activated pathways
axonIdx = unique(fibers(fibers(:,5)==-1,4));
if isempty(axonIdx)
    axonIdx = 0;
end
fprintf('\nLeft side, Damaged axons: %d out of %d\n', length(axonIdx), length(idx));

% Load result from right side
if isfile([axonActivationFolder, filesep, 'axonActivation_right.mat'])
    load([axonActivationFolder, filesep, 'axonActivation_right.mat'], 'fibers', 'idx');
end

% Activated pathways
axonIdx = unique(fibers(fibers(:,5)==1,4));
if isempty(axonIdx)
    axonIdx = 0;
end
fprintf('\n\nRight side, Activated axons: %d out of %d\n', length(axonIdx), length(idx));

% Activated pathways
axonIdx = unique(fibers(fibers(:,5)==0,4));
if isempty(axonIdx)
    axonIdx = 0;
end
fprintf('\nRight side, Non-activated axons: %d out of %d\n', length(axonIdx), length(idx));

% Activated pathways
axonIdx = unique(fibers(fibers(:,5)==-1,4));
if isempty(axonIdx)
    axonIdx = 0;
end
fprintf('\nRight side, Damaged axons: %d out of %d\n\n', length(axonIdx), length(idx));
