function ea_fiberactivation_stats(stimFolder)
% Summary of fiber activation results

results = ea_regexpdir(stimFolder, 'fiberactivation.*\.mat$', 0);

stimParams = ea_regexpdir(stimFolder, 'stimparameters\.mat$', 0);
load(stimParams{1}, 'S');
modelLabel = ea_simModel2Label(S.model);

for i=1:length(results)
    % Load result
    load(results{i}, 'fibers', 'idx');

    hemiSuffix = regexp(results{i}, ['(?<=fiberactivation_model-',modelLabel,'_hemi-).+(?=\.mat$)'], 'match', 'once');
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
end

