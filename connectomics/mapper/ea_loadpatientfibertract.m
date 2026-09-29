function [fibers, idx, voxmm, mat, fiberFile] = ea_loadpatientfibertract(patientDir, useNativeSeed, prefs)
% Load patient-specific fibers and restore their per-point fiber indices.

[fiberFile, candidates] = ea_resolvepatientfibertract(patientDir, useNativeSeed, prefs);
if isempty(fiberFile)
    ea_error(sprintf('Patient-specific fiber tract file not found. Checked:\n%s', ...
        strjoin(candidates, newline)));
end

[fibers, idx, voxmm, mat] = ea_loadfibertracts(fiberFile);

% BIDS fiber tracking stores N-by-3 coordinates plus a vector containing
% each fiber's point count. Mapper algorithms still operate on an N-by-4
% array whose fourth column identifies the fiber for every point.
if size(fibers, 2) == 3
    idx = idx(:);
    if any(~isfinite(idx)) || any(idx < 0) || any(idx ~= fix(idx)) || ...
            sum(double(idx)) ~= size(fibers, 1)
        ea_error(sprintf(['Invalid fiber index vector in:\n%s\n', ...
            'Expected nonnegative integer lengths summing to %d points.'], ...
            fiberFile, size(fibers, 1)));
    end
    fiberNumber = cast((1:numel(idx))', 'like', fibers);
    fibers(:, 4) = repelem(fiberNumber, idx);
end
