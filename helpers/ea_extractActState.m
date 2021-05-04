function activationStatus = ea_extractActState(resultFile, connectome)
% Extract connectome fiber activation state from OSS-DBS result
% Return a binary vector of the same size as the connectome idx,

% Load result
result = load(resultFile, 'fibers');

% Get indices of the activated fibers
fiberIdx = unique(result.fibers(result.fibers(:,5)==1,4));

% Load connectome idx to retrieve the size
connIdx = load([ea_getconnectomebase('dMRI'), connectome, filesep, 'data.mat'], 'idx');

% Construct the activation status label
activationStatus = zeros(size(connIdx.idx));
activationStatus(fiberIdx) = 1;
activationStatus = uint8(activationStatus);
