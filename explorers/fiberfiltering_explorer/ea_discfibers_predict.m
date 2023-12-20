function score = ea_discfibers_predict(atlas, groupAnalysis, improvement, mode, selection, side, posneg)

if isfolder(groupAnalysis)
    load(ea_getGroupAnalysisFile(groupAnalysis), 'M');
elseif isfile(groupAnalysis)
    load(groupAnalysis, 'M');
end

if ~exist('mode', 'var') || isempty(mode)
    mode = 'bin';
end

if ~exist('selection', 'var') || isempty(selection)
    selection = 1:length(M.patient.list);
end

if ~exist('side', 'var') || isempty(side)
    side = 'both';
end

if ~exist('posneg', 'var') || isempty(posneg)
    posneg = 'both';
end

patientList = M.patient.list(selection);
stimFolder = ['stimulations', filesep, ea_nt(0), 'gs_', M.guid];
rightVATs = strcat(patientList, [filesep, stimFolder, filesep, 'vat_efield_right.nii']);
leftVATs = strcat(patientList, [filesep, stimFolder, filesep, 'vat_efield_left.nii']);

switch lower(side)
    case 'right'
        vatlist = rightVATs;
    case 'left'
        vatlist = leftVATs;
    case 'both'
        vatlist = [rightVATs;leftVATs];
end

[scoreBin, scoreSum, scoreMean, scorePeak, score5Peak] = ea_discfibers_vtascore(vatlist, atlas, 'both', posneg);

switch lower(mode)
    case 'bin'
        score = scoreBin;
    case 'sum'
        score = scoreSum;
    case 'mean'
        score = scoreMean;
    case 'peak'
        score = scorePeak;
    case '5peak'
        score = score5Peak;
end

if strcmpi(side, 'both')
    score = sum([score(1:length(score)/2), score(length(score)/2+1:end)], 2);
end

ea_corrplot(improvement, score, 'noperm', {'Discfiber Prediction', 'Empirical', 'Fiber T-scores'});
