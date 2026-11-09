function nullvals = ea_unified_permutation_nulldist(permdata, obj)
% Builds the permutation null distribution for one group/side cell of the
% unified mapping explorer's 'Permutation Threshold (...)' multcompstrategy
% branch. Shuffles permdata.outcomein across patients obj.multcompNperm
% times (obj.rngseed) and recomputes permdata.statfile on permdata.valsin
% for every shuffle, in parallel when Parallel Computing Toolbox is
% available.
%
% permdata - struct with fields:
%              .valsin    V x N matrix fed into the stat test (features x patients)
%              .outcomein N x 1 outcome vector fed into the stat test
%              .statfile  name of the stat-test function (feval-able,
%                         same interface as ea_explorer_stats_*.m)
%              .H0        H0 argument passed through to the stat-test
%              .patientidx (optional) N x 1 patient of each column. With
%                         mirrored sides a patient has two columns (the
%                         original and the mirrored copy); patients are
%                         shuffled, so both copies always receive the same
%                         permuted score. Without it, every column is
%                         treated as its own patient.
% obj      - explorer object; uses obj.multcompNperm, obj.rngseed and
%            obj.multcompmaxworkers
%
% Returns:
% nullvals - V x Nperm matrix, column p is the recomputed statistic for
%            every feature under permutation p.
%
% Kept separate from ea_unified_permutation_threshold.m (which turns a null
% distribution into p-values/significance) because the 'max-statistics'
% strategy needs every cell's null distribution collected first, before it
% can build the single shared null shared across all cells -- see
% ea_unified_corrsignan.m.

Nperm = obj.multcompNperm;

[outcomePatient, patientOfColumn] = local_patientoutcome(permdata);
IpermPatient = ea_shuffle_grouped(outcomePatient, Nperm, [], [], obj.rngseed);
Iperm = IpermPatient(patientOfColumn, :);

statfile = permdata.statfile;
valsin = permdata.valsin;
H0 = permdata.H0;

nullvals = nan(size(valsin,1), Nperm);

usePar = license('test', 'Distrib_Computing_Toolbox') && ~isempty(ver('parallel'));

if usePar
    try
        % Cap the pool at obj.multcompmaxworkers, same idea as
        % ea_unified_perm_cap_parpool.m in the permutation/ subfolder, but
        % inlined here so this function has no dependency on a nested
        % folder that may not be on the MATLAB path.
        pool = gcp('nocreate');
        if isempty(pool)
            parpool(obj.multcompmaxworkers);
        elseif pool.NumWorkers > obj.multcompmaxworkers
            delete(pool);
            parpool(obj.multcompmaxworkers);
        end
        progress = 0;
        dq = parallel.pool.DataQueue;
        afterEach(dq, @(~) reportProgress());
        fprintf('Running %d permutations (parallel, up to %d workers)...\n', Nperm, obj.multcompmaxworkers);
        parfor p = 1:Nperm
            nullvals(:,p) = feval(statfile, valsin, Iperm(:,p), H0);
            send(dq, 1);
        end
    catch ME
        warning('ea_unified_permutation_nulldist:parpoolFailed', ...
            'Falling back to a serial permutation loop (parallel pool failed: %s).', ME.message);
        usePar = false;
    end
end

if ~usePar
    fprintf('Running %d permutations (serial -- Parallel Computing Toolbox not available)...\n', Nperm);
    for p = 1:Nperm
        nullvals(:,p) = feval(statfile, valsin, Iperm(:,p), H0);
    end
end

    function reportProgress()
        % Called on the client (not inside the workers) each time a
        % worker finishes one permutation, via the DataQueue above.
        progress = progress + 1;
        step = max(1, round(Nperm/20)); % ~5% increments
        if mod(progress, step) == 0 || progress == Nperm
            fprintf('Permutation progress: %d/%d (%.0f%%)\n', progress, Nperm, 100*progress/Nperm);
        end
    end
end


function [outcomePatient, patientOfColumn] = local_patientoutcome(permdata)
% One outcome value per patient, and the patient of each column of
% permdata.valsin, so that shuffling happens across patients and every
% copy of a patient (mirrored sides) receives the same permuted score.
outcome = permdata.outcomein(:);
if ~isfield(permdata, 'patientidx') || isempty(permdata.patientidx)
    outcomePatient = outcome;
    patientOfColumn = (1:numel(outcome))';
    return
end

[~, firstColumn, patientOfColumn] = unique(permdata.patientidx(:), 'stable');
outcomePatient = outcome(firstColumn);

% All copies of a patient carry the same outcome (the outcome is
% duplicated before mirroring); anything else would make a joint shuffle
% ambiguous, so refuse instead of silently picking one copy's value.
copyOutcome = outcomePatient(patientOfColumn);
copiesDiffer = ~(outcome == copyOutcome | (isnan(outcome) & isnan(copyOutcome)));
if any(copiesDiffer)
    error('ea_unified_permutation_nulldist:copiesDiffer', ...
        'The copies of %d patient(s) have different outcome values, so they cannot be shuffled together.', ...
        numel(unique(patientOfColumn(copiesDiffer))));
end
end
