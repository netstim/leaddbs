function nWorkers = ea_permtest_parpool(maxWorkers)
% Prepare a parallel pool sized for the machine and return the number of
% workers parfor should use (0 = run serially on the client).
%
% Results never depend on the number of workers: all shuffles are fixed
% in advance by the permset, so 0, 3 or 8 workers give identical output.
% The limit is memory: every worker is a separate MATLAB process (roughly
% 1 GB of RAM before it does any work), so the default adapts to the
% machine instead of using one worker per core.
%
% maxWorkers - 'auto' (default): number of cores - 1, at most 4 (4 on a
%              12-core laptop, 3 on a 4-core one); or a number. 0 or 1
%              runs serially.
%
% Falls back to serial without the Parallel Computing Toolbox. An existing
% pool is reused if it is not larger than the cap; a larger pool is
% replaced by a capped one (same behaviour as the sweetspot permutation
% code, ea_unified_perm_cap_parpool.m).

if ~exist('maxWorkers', 'var') || isempty(maxWorkers)
    maxWorkers = 'auto';
end
if ischar(maxWorkers) || isstring(maxWorkers)
    if ~strcmpi(maxWorkers, 'auto')
        ea_error('maxWorkers must be ''auto'' or a number.');
    end
    maxWorkers = max(1, min(4, feature('numcores') - 1));
end

hasPCT = license('test', 'Distrib_Computing_Toolbox') && ~isempty(ver('parallel'));
if maxWorkers <= 1 || ~hasPCT
    nWorkers = 0;
    if maxWorkers > 1
        fprintf('Parallel Computing Toolbox not available -- running serially.\n');
    end
    return
end

p = gcp('nocreate');
if isempty(p)
    p = parpool(maxWorkers);
elseif p.NumWorkers > maxWorkers
    delete(p);
    p = parpool(maxWorkers);
end
nWorkers = p.NumWorkers;
