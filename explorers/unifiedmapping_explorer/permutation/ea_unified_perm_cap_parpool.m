function ea_unified_perm_cap_parpool(maxWorkers)
% Ensures the current parallel pool (if any) has no more than maxWorkers
% workers before a parfor loop starts. Each worker holds its own full copy
% of the ea_unifiedmapping object being permuted (including
% obj.results.sweetspotmapping.efield), so an uncapped pool (one worker per
% core) can exhaust system memory on large analyses.

p = gcp('nocreate');
if isempty(p)
    parpool(maxWorkers);
elseif p.NumWorkers > maxWorkers
    delete(p);
    parpool(maxWorkers);
end
